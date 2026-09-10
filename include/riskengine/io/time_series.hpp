#pragma once

#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

// Daily time series read from the frozen FRED files in data/raw (tools/fetch_data.py).
namespace riskengine {

struct Observation {
    std::string date; // ISO 8601 (YYYY-MM-DD), so lexicographic order is chronological
    double value;
};

// True for a valid calendar date written YYYY-MM-DD (proleptic Gregorian).
inline bool is_iso_date(std::string_view d) {
    if (d.size() != 10 || d[4] != '-' || d[7] != '-') return false;
    for (std::size_t i : {0, 1, 2, 3, 5, 6, 8, 9})
        if (d[i] < '0' || d[i] > '9') return false;
    const auto num = [&](std::size_t at, std::size_t n) {
        int v = 0;
        for (std::size_t i = at; i < at + n; ++i) v = 10 * v + (d[i] - '0');
        return v;
    };
    const int year = num(0, 4), month = num(5, 2), day = num(8, 2);
    if (month < 1 || month > 12 || day < 1) return false;
    const bool leap = (year % 4 == 0 && year % 100 != 0) || year % 400 == 0;
    constexpr int kDays[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    return day <= kDays[month - 1] + (month == 2 && leap ? 1 : 0);
}

// Reads a FRED CSV ("observation_date,<ID>" then "YYYY-MM-DD,value"). Days without a value (an
// empty field or ".", e.g. market holidays) are skipped. Throws on anything else unexpected (an
// invalid calendar date, a non-numeric or non-finite value, dates out of order), so a corrupted
// file can never be read as a shorter or different series.
inline std::vector<Observation> read_fred_csv(const std::filesystem::path& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open " + path.string());
    std::string line;
    if (!std::getline(in, line) || line.rfind("observation_date,", 0) != 0)
        throw std::runtime_error(path.string() + ": not a FRED CSV");
    std::vector<Observation> series;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        const auto comma = line.find(',');
        if (comma != 10 || !is_iso_date(std::string_view(line).substr(0, 10)))
            throw std::runtime_error(path.string() + ": malformed line '" + line + "'");
        const std::string value = line.substr(comma + 1);
        if (value.empty() || value == ".") continue;
        std::size_t used = 0;
        double v = 0.0;
        try {
            v = std::stod(value, &used);
        } catch (const std::logic_error&) { // invalid_argument, out_of_range: reported below
            used = 0;
        }
        if (used != value.size() || !std::isfinite(v))
            throw std::runtime_error(path.string() + ": bad value '" + value + "'");
        if (!series.empty() && line.substr(0, 10) <= series.back().date)
            throw std::runtime_error(path.string() + ": dates not increasing at " + line.substr(0, 10));
        series.push_back({line.substr(0, 10), v});
    }
    return series;
}

// Values of `b` on the dates of `a` that `b` also has: the inner join of two date-sorted series.
struct AlignedPair {
    std::vector<std::string> dates;
    std::vector<double> a, b;
};

inline AlignedPair align(const std::vector<Observation>& a, const std::vector<Observation>& b) {
    AlignedPair out;
    std::size_t j = 0;
    for (const Observation& x : a) {
        while (j < b.size() && b[j].date < x.date) ++j;
        if (j < b.size() && b[j].date == x.date) {
            out.dates.push_back(x.date);
            out.a.push_back(x.value);
            out.b.push_back(b[j].value);
        }
    }
    return out;
}

// log(x[i + 1] / x[i]).
inline std::vector<double> log_returns(const std::vector<double>& x) {
    std::vector<double> r;
    for (std::size_t i = 1; i < x.size(); ++i) r.push_back(std::log(x[i] / x[i - 1]));
    return r;
}

} // namespace riskengine
