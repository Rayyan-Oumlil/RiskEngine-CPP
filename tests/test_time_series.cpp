#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <stdexcept>

#include "riskengine/io/time_series.hpp"

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

TEST_CASE("The frozen FRED files parse completely and match their manifest", "[time_series][data]") {
    // Counts from data/raw/manifest.json (observations_with_value), written by tools/fetch_data.py.
    const std::filesystem::path dir = RISKENGINE_DATA_DIR;
    struct Expected {
        const char* file;
        std::size_t observations;
        const char* first;
        const char* last;
    };
    for (const Expected& e : {Expected{"NASDAQCOM.csv", 14025, "1971-02-05", "2026-09-22"},
                              Expected{"VIXCLS.csv", 9279, "1990-01-02", "2026-09-22"},
                              Expected{"VXOCLS.csv", 9002, "1986-01-02", "2021-09-23"},
                              Expected{"DGS3MO.csv", 11264, "1981-09-01", "2026-09-22"}}) {
        INFO(e.file);
        const auto s = read_fred_csv(dir / e.file);
        CHECK(s.size() == e.observations);
        CHECK(s.front().date == e.first);
        CHECK(s.back().date == e.last);
    }
}
