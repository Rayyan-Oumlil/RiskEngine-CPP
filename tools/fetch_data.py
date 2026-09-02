"""Download the historical series used by the risk-measure experiments (report 7) and freeze them.

Every series is saved exactly as served by FRED (Federal Reserve Bank of St. Louis) into
data/raw/<id>.csv, and data/raw/manifest.json records, per file, the source URL, the UTC time of
extraction, the SHA-256 of the bytes, the number of observations and the date range. Experiments
only ever read the frozen files (docs/riskengine_research.md 9.1): rerun this script only to refresh
the data on purpose, and commit the result.

Series:
  NASDAQCOM  NASDAQ Composite index, daily close (1971-): the underlying. FRED's S&P 500 series
             only covers the last ten years, which would miss 1987, 2008 and 2018.
  VIXCLS     CBOE VIX (1990-): implied-volatility proxy.
  VXOCLS     CBOE VXO (1986-2021): the implied-volatility proxy before 1990 (October 1987).
  DGS3MO     3-month Treasury constant-maturity yield, percent (1981-): risk-free rate.

Usage: python3 tools/fetch_data.py
"""
import datetime
import hashlib
import json
import pathlib
import urllib.request

SERIES = ["NASDAQCOM", "VIXCLS", "VXOCLS", "DGS3MO"]
URL = "https://fred.stlouisfed.org/graph/fredgraph.csv?id={}"
OUT = pathlib.Path(__file__).resolve().parent.parent / "data" / "raw"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    manifest = {"source": "FRED, Federal Reserve Bank of St. Louis", "files": {}}
    for series in SERIES:
        url = URL.format(series)
        with urllib.request.urlopen(url, timeout=60) as response:
            data = response.read()
        lines = data.decode().strip().splitlines()
        if not lines or not lines[0].startswith("observation_date"):
            raise SystemExit(f"{series}: unexpected response, not a FRED CSV")
        observations = [line.split(",") for line in lines[1:]]
        valued = [row for row in observations if len(row) == 2 and row[1] not in ("", ".")]
        (OUT / f"{series}.csv").write_bytes(data)
        manifest["files"][f"{series}.csv"] = {
            "url": url,
            "extracted_utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "sha256": hashlib.sha256(data).hexdigest(),
            "rows": len(observations),
            "observations_with_value": len(valued),
            "first_date": valued[0][0],
            "last_date": valued[-1][0],
        }
        print(f"{series}: {len(valued)} observations, {valued[0][0]} .. {valued[-1][0]}")
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"wrote {OUT / 'manifest.json'}")


if __name__ == "__main__":
    main()
