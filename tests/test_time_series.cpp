#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <stdexcept>
#include <string>

#include "riskengine/io/time_series.hpp"
#include "sha256.hpp"

using namespace riskengine;

namespace {

std::filesystem::path write_temp(const std::string& name, const std::string& content) {
    const auto path = std::filesystem::temp_directory_path() / name;
    std::ofstream(path, std::ios::binary) << content;
    return path;
}

} // namespace

TEST_CASE("FRED CSV: values read, holidays skipped", "[time_series]") {
    const auto path = write_temp("riskengine_fred_ok.csv",
                                 "observation_date,TEST\n2020-01-02,10.5\n2020-01-03,\n2020-01-06,.\r\n2020-01-07,11\n");
    const auto s = read_fred_csv(path);
    REQUIRE(s.size() == 2);
    CHECK(s[0].date == "2020-01-02");
    CHECK(s[0].value == 10.5);
    CHECK(s[1].date == "2020-01-07");
    CHECK(s[1].value == 11.0);
}

TEST_CASE("FRED CSV: malformed input is an error, never a shorter series", "[time_series]") {
    CHECK_THROWS_AS(read_fred_csv(write_temp("riskengine_fred_hdr.csv", "<html>\n")), std::runtime_error);
    CHECK_THROWS_AS(read_fred_csv(write_temp("riskengine_fred_val.csv", "observation_date,X\n2020-01-02,1.2x\n")),
                    std::runtime_error);
    CHECK_THROWS_AS(read_fred_csv(write_temp("riskengine_fred_ord.csv",
                                             "observation_date,X\n2020-01-03,1\n2020-01-02,2\n")),
                    std::runtime_error);
    // Not a calendar date, although the comma is in the right place.
    for (const char* date : {"2020-99-99", "2021-02-29", "2020-04-31", "2020-00-10", "20200102xx", "2020/01/02"}) {
        INFO(date);
        CHECK_THROWS_AS(read_fred_csv(write_temp("riskengine_fred_date.csv",
                                                 std::string("observation_date,X\n") + date + ",1\n")),
                        std::runtime_error);
    }
    // Non-finite or unparsable values, including those std::stod reports with its own exceptions.
    for (const char* value : {"nan", "inf", "-inf", "1e999", "abc"}) {
        INFO(value);
        CHECK_THROWS_AS(read_fred_csv(write_temp("riskengine_fred_nan.csv",
                                                 std::string("observation_date,X\n2020-01-02,") + value + "\n")),
                        std::runtime_error);
    }
    CHECK(read_fred_csv(write_temp("riskengine_fred_leap.csv", "observation_date,X\n2020-02-29,1\n")).size() == 1);
    CHECK_THROWS_AS(read_fred_csv("/nonexistent/riskengine.csv"), std::runtime_error);
}

TEST_CASE("Series are aligned on common dates; log returns", "[time_series]") {
    const std::vector<Observation> a{{"2020-01-01", 1}, {"2020-01-02", 2}, {"2020-01-03", 3}, {"2020-01-05", 5}};
    const std::vector<Observation> b{{"2020-01-02", 20}, {"2020-01-04", 40}, {"2020-01-05", 50}};
    const AlignedPair p = align(a, b);
    REQUIRE(p.dates == std::vector<std::string>{"2020-01-02", "2020-01-05"});
    CHECK(p.a == std::vector<double>{2, 5});
    CHECK(p.b == std::vector<double>{20, 50});
    const auto r = log_returns({100, 110, 99});
    CHECK(std::abs(r[0] - std::log(1.1)) < 1e-15);
    CHECK(std::abs(r[1] - std::log(0.9)) < 1e-15);
}

TEST_CASE("SHA-256 matches the FIPS 180-4 test vectors", "[time_series][data]") {
    using riskengine::test::sha256_hex;
    CHECK(sha256_hex("") == "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855");
    CHECK(sha256_hex("abc") == "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad");
    CHECK(sha256_hex("abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq") ==
          "248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c1");
}

TEST_CASE("The frozen FRED files match their manifest and parse completely", "[time_series][data]") {
    // Checksums, counts (observations_with_value) and date ranges from data/raw/manifest.json,
    // written by tools/fetch_data.py. The checksum is also looked up in the manifest itself, so the
    // test fails if either the data or the manifest changes without the other.
    const std::filesystem::path dir = RISKENGINE_DATA_DIR;
    const auto slurp = [](const std::filesystem::path& p) {
        std::ifstream in(p, std::ios::binary);
        REQUIRE(in);
        return std::string(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
    };
    const std::string manifest = slurp(dir / "manifest.json");
    struct Expected {
        const char* file;
        const char* sha256;
        std::size_t observations;
        const char* first;
        const char* last;
    };
    for (const Expected& e :
         {Expected{"NASDAQCOM.csv", "f57dabe1d01ac6bbf8eace17071443f43872c668e22d4285457533799383409f", 14025,
                   "1971-02-05", "2026-09-22"},
          Expected{"VIXCLS.csv", "675f9bcb86ee88799f8f7a8b7ee9c30f7a0be26a26df8f9549bddf9e84c4a598", 9279,
                   "1990-01-02", "2026-09-22"},
          Expected{"VXOCLS.csv", "3ab3e522630a6b404fda6cd2d7ef33834139218bf00f5a290aa089361dd334e8", 9002,
                   "1986-01-02", "2021-09-23"},
          Expected{"DGS3MO.csv", "e5b1e895403017ba74bcfcc4e592da130989784eca5f6c2197c6e7168e9d9155", 11264,
                   "1981-09-01", "2026-09-22"}}) {
        INFO(e.file);
        CHECK(manifest.find(std::string("\"sha256\": \"") + e.sha256 + "\"") != std::string::npos);
        CHECK(riskengine::test::sha256_hex(slurp(dir / e.file)) == e.sha256);
        const auto s = read_fred_csv(dir / e.file);
        CHECK(s.size() == e.observations);
        CHECK(s.front().date == e.first);
        CHECK(s.back().date == e.last);
    }
}
