#pragma once

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "taxonomy.hpp"

namespace homval {

struct Denoms {
    int ndom = 0;
    uint64_t n_possible_tp = 0;
    uint64_t n_possible_fp = 0;
};

inline int64_t json_int_after_key(const std::string& json, const std::string& key) {
    const std::string needle = "\"" + key + "\"";
    size_t pos = json.find(needle);
    if (pos == std::string::npos)
        die_tax("derived_info missing key " + key);
    pos = json.find(':', pos + needle.size());
    if (pos == std::string::npos)
        die_tax("derived_info malformed at " + key);
    ++pos;
    while (pos < json.size() && (json[pos] == ' ' || json[pos] == '\t'))
        ++pos;
    size_t end = pos;
    while (end < json.size() && (json[end] == '-' || std::isdigit(static_cast<unsigned char>(json[end]))))
        ++end;
    return std::stoll(json.substr(pos, end - pos));
}

inline int64_t json_truth_int(const std::string& json, const std::string& section, const std::string& truth) {
    const std::string sec = "\"" + section + "\"";
    size_t pos = json.find(sec);
    if (pos == std::string::npos)
        die_tax("derived_info missing " + section);
    pos = json.find('{', pos);
    if (pos == std::string::npos)
        die_tax("derived_info malformed " + section);
    const size_t end = json.find('}', pos);
    if (end == std::string::npos)
        die_tax("derived_info malformed " + section);
    const std::string block = json.substr(pos, end - pos + 1);
    const std::string key = "\"" + truth + "\"";
    size_t kpos = block.find(key);
    if (kpos == std::string::npos)
        die_tax("derived_info missing " + section + " for truth " + truth);
    kpos = block.find(':', kpos + key.size());
    if (kpos == std::string::npos)
        die_tax("derived_info malformed " + section + " for truth " + truth);
    ++kpos;
    while (kpos < block.size() && (block[kpos] == ' ' || block[kpos] == '\t'))
        ++kpos;
    size_t e = kpos;
    while (e < block.size() && (block[e] == '-' || std::isdigit(static_cast<unsigned char>(block[e]))))
        ++e;
    return std::stoll(block.substr(kpos, e - kpos));
}

inline Denoms load_denoms(const std::string& path, const std::string& truth) {
    const std::string json = read_file(path);
    Denoms d;
    d.ndom = static_cast<int>(json_int_after_key(json, "ndom"));
    d.n_possible_tp = static_cast<uint64_t>(json_truth_int(json, "N_possible_tp", truth));
    d.n_possible_fp = static_cast<uint64_t>(json_truth_int(json, "N_possible_fp", truth));
    return d;
}

inline std::string load_reference(const std::string& path) {
    const std::string json = read_file(path);
    return json_string_after_key(json, "reference");
}

struct QueryScore {
    bool has_tp = false;
    bool has_xf = false;
    double s_tp = 0.0;
    double s_xf = 0.0;
};

struct QueryState {
    std::unordered_map<uint32_t, QueryScore> scores;
    std::vector<uint32_t> query_order;
    int n_possible_tp = 0;
};

inline QueryState load_query_state(
    const std::string& derived_info_path,
    const std::string& truth_str,
    const Taxonomy& tax) {
    const std::string json = read_file(derived_info_path);
    const std::string key = "query_set_" + truth_str;
    const std::vector<std::string> domains = json_string_array_after_key(json, key);
    const int n_possible_tp = static_cast<int>(json_truth_int(json, "N_possible_tp", truth_str));

    QueryState state;
    state.n_possible_tp = n_possible_tp;
    state.query_order.reserve(domains.size());

    std::vector<std::string> sorted = domains;
    std::sort(sorted.begin(), sorted.end());
    for (const std::string& dom : sorted) {
        auto it = tax.dom_index.find(dom);
        if (it == tax.dom_index.end())
            die_tax("query domain not in derived_info: " + dom);
        state.query_order.push_back(it->second);
        state.scores[it->second] = QueryScore{};
    }
    return state;
}

struct LookupStrings {
    std::vector<std::string> domain_ids;
    std::vector<std::string> dom_fold;
    std::vector<std::string> dom_sf;
    std::vector<std::string> dom_fam;
};

inline LookupStrings load_lookup_strings(const std::string& path) {
    LookupStrings data;
    std::ifstream in(path);
    if (!in)
        die_tax("cannot open lookup " + path);

    std::string line;
    int line_no = 0;
    while (std::getline(in, line)) {
        ++line_no;
        trim_inplace(line);
        if (line.empty())
            continue;
        const size_t tab = line.find('\t');
        if (tab == std::string::npos)
            die_tax(path + ":" + std::to_string(line_no) + ": expected 2 tab fields");

        const std::string dom = line.substr(0, tab);
        const std::string fam = line.substr(tab + 1);
        for (const std::string& seen : data.domain_ids) {
            if (seen == dom)
                die_tax(path + ":" + std::to_string(line_no) + ": duplicate domain " + dom);
        }

        const size_t d1 = fam.find('.');
        const size_t d2 = (d1 == std::string::npos) ? std::string::npos : fam.find('.', d1 + 1);
        const size_t d3 = (d2 == std::string::npos) ? std::string::npos : fam.find('.', d2 + 1);
        if (d1 == std::string::npos || d2 == std::string::npos || d3 == std::string::npos)
            die_tax(path + ":" + std::to_string(line_no) + ": family_id must have 4 dot fields");

        data.domain_ids.push_back(dom);
        data.dom_fold.push_back(fam.substr(0, d2));
        data.dom_sf.push_back(fam.substr(0, d3));
        data.dom_fam.push_back(fam);
    }
    return data;
}

inline uint64_t ordered_pairs(uint64_t n) {
    return n * (n - 1);
}

inline std::pair<uint64_t, uint64_t> count_pair_denominators(
    const LookupStrings& data,
    const std::string& truth) {
    const uint64_t n = data.domain_ids.size();
    std::unordered_map<std::string, std::vector<size_t>> groups;

    if (truth == "fold") {
        for (size_t i = 0; i < n; ++i)
            groups[data.dom_fold[i]].push_back(i);
    } else if (truth == "superfamily" || truth == "superfamilyx") {
        for (size_t i = 0; i < n; ++i)
            groups[data.dom_sf[i]].push_back(i);
    } else if (truth == "family") {
        for (size_t i = 0; i < n; ++i)
            groups[data.dom_fam[i]].push_back(i);
    } else {
        die_tax("unknown truth standard: " + truth);
    }

    uint64_t n_tp = 0;
    for (const auto& [_, ids] : groups)
        n_tp += ordered_pairs(static_cast<uint64_t>(ids.size()));

    if (truth != "superfamilyx")
        return {n_tp, ordered_pairs(n) - n_tp};

    std::unordered_map<std::string, std::set<std::string>> fold2sfs;
    std::unordered_map<std::string, std::vector<size_t>> fold2doms;
    std::unordered_map<std::string, std::vector<size_t>> sf2doms;
    for (size_t i = 0; i < n; ++i) {
        fold2sfs[data.dom_fold[i]].insert(data.dom_sf[i]);
        fold2doms[data.dom_fold[i]].push_back(i);
        sf2doms[data.dom_sf[i]].push_back(i);
    }

    uint64_t n_ignore = 0;
    for (const auto& [fold, sfs] : fold2sfs) {
        const uint64_t n_fold = static_cast<uint64_t>(fold2doms[fold].size());
        uint64_t n_tp_in_fold = 0;
        for (const std::string& sf : sfs)
            n_tp_in_fold += ordered_pairs(static_cast<uint64_t>(sf2doms[sf].size()));
        n_ignore += ordered_pairs(n_fold) - n_tp_in_fold;
    }
    return {n_tp, ordered_pairs(n) - n_tp - n_ignore};
}

inline std::vector<std::string> query_set_topfold(const LookupStrings& data) {
    std::unordered_map<std::string, std::set<std::string>> fold2sfs;
    for (size_t i = 0; i < data.domain_ids.size(); ++i)
        fold2sfs[data.dom_fold[i]].insert(data.dom_sf[i]);

    std::vector<std::string> out;
    for (size_t i = 0; i < data.domain_ids.size(); ++i) {
        if (fold2sfs[data.dom_fold[i]].size() > 1)
            out.push_back(data.domain_ids[i]);
    }
    std::sort(out.begin(), out.end());
    return out;
}

inline std::vector<std::string> query_set_topsf(const LookupStrings& data) {
    std::unordered_map<std::string, std::set<std::string>> sf2fams;
    for (size_t i = 0; i < data.domain_ids.size(); ++i)
        sf2fams[data.dom_sf[i]].insert(data.dom_fam[i]);

    std::vector<std::string> out;
    for (size_t i = 0; i < data.domain_ids.size(); ++i) {
        if (sf2fams[data.dom_sf[i]].size() > 1)
            out.push_back(data.domain_ids[i]);
    }
    std::sort(out.begin(), out.end());
    return out;
}

inline std::string json_escape(const std::string& s) {
    std::string out;
    out.reserve(s.size() + 2);
    for (char c : s) {
        if (c == '"' || c == '\\')
            out.push_back('\\');
        out.push_back(c);
    }
    return out;
}

inline void write_json_string_array(std::ostringstream& os, const std::vector<std::string>& values) {
    os << "[\n";
    for (size_t i = 0; i < values.size(); ++i) {
        os << "    \"" << json_escape(values[i]) << "\"";
        if (i + 1 < values.size())
            os << ',';
        os << '\n';
    }
    os << "  ]";
}

inline std::string build_derived_info_json(const LookupStrings& data, const std::string& reference) {
    static const char* const TRUTH_STANDARDS[] = {
        "fold", "superfamily", "family", "superfamilyx"};
    const std::vector<std::string> qs_topfold = query_set_topfold(data);
    const std::vector<std::string> qs_topsf = query_set_topsf(data);

    std::map<std::string, int64_t> n_tp_map;
    std::map<std::string, int64_t> n_fp_map;
    for (const char* truth : TRUTH_STANDARDS) {
        const auto [n_tp, n_fp] = count_pair_denominators(data, truth);
        n_tp_map[truth] = static_cast<int64_t>(n_tp);
        n_fp_map[truth] = static_cast<int64_t>(n_fp);
    }

    std::ostringstream os;
    os << "{\n";
    os << "  \"ndom\": " << data.domain_ids.size() << ",\n";
    os << "  \"reference\": \"" << json_escape(reference) << "\",\n";
    os << "  \"domain_ids\": ";
    write_json_string_array(os, data.domain_ids);
    os << ",\n";
    os << "  \"dom_fold\": ";
    write_json_string_array(os, data.dom_fold);
    os << ",\n";
    os << "  \"dom_sf\": ";
    write_json_string_array(os, data.dom_sf);
    os << ",\n";
    os << "  \"dom_fam\": ";
    write_json_string_array(os, data.dom_fam);
    os << ",\n";
    os << "  \"query_set_topfold\": ";
    write_json_string_array(os, qs_topfold);
    os << ",\n";
    os << "  \"query_set_topsf\": ";
    write_json_string_array(os, qs_topsf);
    os << ",\n";
    os << "  \"N_possible_tp\": {\n";
    os << "    \"topfold\": " << qs_topfold.size() << ",\n";
    os << "    \"topsf\": " << qs_topsf.size() << ",\n";
    bool first = true;
    for (const char* truth : TRUTH_STANDARDS) {
        if (!first)
            os << ",\n";
        first = false;
        os << "    \"" << truth << "\": " << n_tp_map[truth];
    }
    os << "\n  },\n";
    os << "  \"N_possible_fp\": {\n";
    first = true;
    for (const char* truth : TRUTH_STANDARDS) {
        if (!first)
            os << ",\n";
        first = false;
        os << "    \"" << truth << "\": " << n_fp_map[truth];
    }
    os << "\n  }\n";
    os << "}\n";
    return os.str();
}

}  // namespace homval
