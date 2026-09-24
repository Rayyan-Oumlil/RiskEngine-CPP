#pragma once

#include <chrono>
#include <cstdint>
#include <cstdio>
#include <exception>
#include <filesystem>
#include <format>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <thread>
#include <type_traits>
#include <utility>
#include <vector>

#include "riskengine/git_info.hpp"

// Experiment harness: one executable = one figure or table of docs/model_risk_report.md.
//
//   int main(int argc, char** argv) {
//       return riskengine::harness::run("my_experiment", argc, argv, [](Experiment& exp) {
//           exp.param("seed", 2026);
//           auto csv = exp.csv({"x", "y"});
//           csv.row(1.0, 2.0);
//       });
//   }
//
// writes <out>/my_experiment.csv and <out>/my_experiment.meta.json (git commit, compiler, build
// type, flags, CPU, parameters, UTC time). <out> defaults to data/results and is set with --out <dir>.
// Output is byte-for-byte deterministic: %.17g round-trip formatting and '\n' line endings on
// every platform.
namespace riskengine::harness {

namespace detail {

inline std::string format_number(double x) {
    char buf[32];
    std::snprintf(buf, sizeof buf, "%.17g", x);
    return buf;
}

inline std::string_view trim(std::string_view s) {
    const auto first = s.find_first_not_of(' ');
    if (first == std::string_view::npos) return {};
    return s.substr(first, s.find_last_not_of(' ') - first + 1);
}

inline std::string json_string(std::string_view s) {
    std::string out = "\"";
    for (char c : s) {
        switch (c) {
            case '"': out += "\\\""; break;
            case '\\': out += "\\\\"; break;
            case '\n': out += "\\n"; break;
            case '\t': out += "\\t"; break;
            default: out += c;
        }
    }
    return out + "\"";
}

// CPU model, for experiments that report timings: "model name" from /proc/cpuinfo on Linux,
// "unknown" elsewhere.
inline std::string cpu_model() {
    std::ifstream cpuinfo("/proc/cpuinfo");
    for (std::string line; std::getline(cpuinfo, line);) {
        if (line.rfind("model name", 0) == 0) {
            const auto colon = line.find(':');
            if (colon != std::string::npos) return std::string(trim(std::string_view(line).substr(colon + 1)));
        }
    }
    return "unknown";
}

template <class T>
std::string to_cell(const T& value) {
    if constexpr (std::is_same_v<T, bool>) {
        return value ? "true" : "false";
    } else if constexpr (std::is_integral_v<T>) {
        return std::to_string(value);
    } else if constexpr (std::is_floating_point_v<T>) {
        return format_number(static_cast<double>(value));
    } else {
        const std::string s{std::string_view{value}};
        if (s.find_first_of(",\"\n") != std::string::npos)
            throw std::invalid_argument("CSV cell must not contain ',', '\"' or newlines: " + s);
        return s;
    }
}

template <class T>
std::string to_json(const T& value) {
    if constexpr (std::is_arithmetic_v<T>) {
        return to_cell(value);
    } else {
        return json_string(std::string_view{value});
    }
}

} // namespace detail

class CsvWriter {
public:
    CsvWriter(const std::filesystem::path& path, std::vector<std::string> columns)
        : out_(path, std::ios::binary), columns_(columns.size()) {
        if (!out_) throw std::runtime_error("cannot open " + path.string());
        write_line(columns);
    }

    template <class... Ts>
    void row(const Ts&... values) {
        if (sizeof...(Ts) != columns_)
            throw std::invalid_argument(std::format("row has {} cells, header has {}", sizeof...(Ts), columns_));
        write_line({detail::to_cell(values)...});
    }

private:
    void write_line(const std::vector<std::string>& cells) {
        for (std::size_t i = 0; i < cells.size(); ++i) out_ << (i > 0 ? "," : "") << cells[i];
        out_ << '\n';
        if (!out_) throw std::runtime_error("write failed");
    }

    std::ofstream out_;
    std::size_t columns_;
};

class Experiment {
public:
    Experiment(std::string id, std::filesystem::path out_dir) : id_(std::move(id)), out_dir_(std::move(out_dir)) {
        std::filesystem::create_directories(out_dir_);
    }

    // Records an input of the run in the metadata (seed, sample count, market parameters, ...).
    template <class T>
    void param(std::string key, const T& value) {
        params_.emplace_back(std::move(key), detail::to_json(value));
    }

    CsvWriter csv(std::vector<std::string> columns) const {
        return CsvWriter(out_dir_ / (id_ + ".csv"), std::move(columns));
    }

    void write_metadata() const {
        const auto path = out_dir_ / (id_ + ".meta.json");
        std::ofstream out(path, std::ios::binary);
        if (!out) throw std::runtime_error("cannot open " + path.string());
        const auto now = std::chrono::floor<std::chrono::seconds>(std::chrono::system_clock::now());
        out << "{\n"
            << "  \"experiment\": " << detail::json_string(id_) << ",\n"
            << "  \"git_commit\": " << detail::json_string(build_info::kGitCommit) << ",\n"
            << "  \"git_dirty\": " << (build_info::kGitDirty ? "true" : "false") << ",\n"
            << "  \"compiler\": " << detail::json_string(RISKENGINE_COMPILER) << ",\n"
            << "  \"build_type\": " << detail::json_string(RISKENGINE_BUILD_TYPE) << ",\n"
            << "  \"cxx_flags\": " << detail::json_string(detail::trim(RISKENGINE_CXX_FLAGS)) << ",\n"
            << "  \"cpu\": " << detail::json_string(detail::cpu_model()) << ",\n"
            << "  \"hardware_threads\": " << std::thread::hardware_concurrency() << ",\n"
            << "  \"created_utc\": " << detail::json_string(std::format("{:%FT%TZ}", now)) << ",\n"
            << "  \"parameters\": {";
        for (std::size_t i = 0; i < params_.size(); ++i)
            out << (i > 0 ? "," : "") << "\n    " << detail::json_string(params_[i].first) << ": " << params_[i].second;
        out << (params_.empty() ? "}" : "\n  }") << "\n}\n";
        if (!out) throw std::runtime_error("write failed: " + path.string());
    }

    const std::string& id() const { return id_; }
    // Where results go: an experiment that also writes a companion table (e.g. a summary) creates a
    // second Experiment with its own id in the same directory.
    const std::filesystem::path& out_dir() const { return out_dir_; }

private:
    std::string id_;
    std::filesystem::path out_dir_;
    std::vector<std::pair<std::string, std::string>> params_;
};

// Parses --out <dir>, runs the experiment, writes its metadata, and reports errors as exit code 1.
template <class Body>
int run(std::string id, int argc, char** argv, Body&& body) {
    std::filesystem::path out_dir = "data/results";
    for (int i = 1; i < argc; ++i) {
        const std::string_view arg = argv[i];
        if (arg == "--out" && i + 1 < argc) {
            out_dir = argv[++i];
        } else {
            std::cerr << "usage: " << argv[0] << " [--out <dir>]\n";
            return 2;
        }
    }
    try {
        Experiment exp(std::move(id), out_dir);
        body(exp);
        exp.write_metadata();
        std::cout << "wrote " << (out_dir / (exp.id() + ".csv")).string() << '\n';
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << '\n';
        return 1;
    }
}

} // namespace riskengine::harness
