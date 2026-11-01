// check.hpp -- minimal assertion counter and JSON-ish result writer for the owner tests.
#pragma once
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace chk {
inline int passed = 0, failed = 0;
inline std::vector<std::string> failures;
inline std::vector<std::pair<std::string, std::string>> results;   // key -> JSON value text
inline void check(bool ok, const std::string& what) {
    if (ok) ++passed; else { ++failed; failures.push_back(what); std::fprintf(stderr, "FAIL: %s\n", what.c_str()); }
}
inline std::string num(double v) { if (!std::isfinite(v)) return "null"; char b[64]; std::snprintf(b, sizeof b, "%.17g", v); return b; }
inline std::string str(const std::string& s) {
    std::string o = "\"";
    for (char c : s) {
        if (c == '"' || c == '\\') { o += '\\'; o += c; continue; }
        if (static_cast<unsigned char>(c) < 0x20) { char b[8]; std::snprintf(b, sizeof b, "\\u%04x", unsigned(static_cast<unsigned char>(c))); o += b; continue; }
        o += c;
    }
    return o + "\"";
}
inline void put(const std::string& k, double v) { results.push_back({k, num(v)}); }
inline void put_text(const std::string& k, const std::string& v) { results.push_back({k, str(v)}); }
inline void put_raw(const std::string& k, const std::string& json) { results.push_back({k, json}); }
inline int finish(const std::string& name, const std::string& out_dir) {
    std::ostringstream o;
    o << "{\n \"test\": " << str(name) << ",\n \"passed\": " << passed << ",\n \"failed\": " << failed << ",\n \"failures\": [";
    for (std::size_t i = 0; i < failures.size(); ++i) o << (i ? ", " : "") << str(failures[i]);
    o << "],\n \"results\": {";
    for (std::size_t i = 0; i < results.size(); ++i) o << (i ? ",\n  " : "\n  ") << str(results[i].first) << ": " << results[i].second;
    o << "\n }\n}\n";
    if (!out_dir.empty()) { std::ofstream f(out_dir + "/" + name + ".json"); f << o.str(); }
    std::printf("%s: %d passed, %d failed\n", name.c_str(), passed, failed);
    return failed ? 1 : 0;
}
}  // namespace chk
