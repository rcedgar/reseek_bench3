#pragma once

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

namespace homval {

[[noreturn]] inline void die_tax(const std::string& msg) {
    std::fprintf(stderr, "error: %s\n", msg.c_str());
    std::exit(1);
}

inline void trim_inplace(std::string& s) {
    while (!s.empty() && (s.back() == '\r' || s.back() == '\n' || s.back() == ' '))
        s.pop_back();
    size_t i = 0;
    while (i < s.size() && (s[i] == ' ' || s[i] == '\t'))
        ++i;
    if (i > 0)
        s.erase(0, i);
}

inline std::string read_file(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in)
        die_tax("cannot open " + path);
    std::ostringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

struct Taxonomy {
    std::unordered_map<std::string, uint32_t> dom_index;
    std::vector<std::string> id_to_dom;
    std::vector<uint32_t> dom_fold;
    std::vector<uint32_t> dom_sf;
    std::vector<uint32_t> dom_fam;
    uint32_t ndom_lookup = 0;
};

inline uint32_t intern_string(
    const std::string& s,
    std::unordered_map<std::string, uint32_t>& pool,
    std::vector<std::string>& strings) {
    auto it = pool.find(s);
    if (it != pool.end())
        return it->second;
    const uint32_t id = static_cast<uint32_t>(strings.size());
    strings.push_back(s);
    pool.emplace(strings.back(), id);
    return id;
}

inline Taxonomy load_lookup(const std::string& path) {
    Taxonomy tax;
    std::ifstream in(path);
    if (!in)
        die_tax("cannot open lookup " + path);

    std::unordered_map<std::string, uint32_t> fold_pool, sf_pool, fam_pool;
    std::vector<std::string> fold_strings, sf_strings, fam_strings;

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
        if (tax.dom_index.count(dom))
            die_tax(path + ":" + std::to_string(line_no) + ": duplicate domain " + dom);

        const size_t d1 = fam.find('.');
        const size_t d2 = (d1 == std::string::npos) ? std::string::npos : fam.find('.', d1 + 1);
        const size_t d3 = (d2 == std::string::npos) ? std::string::npos : fam.find('.', d2 + 1);
        if (d1 == std::string::npos || d2 == std::string::npos || d3 == std::string::npos)
            die_tax(path + ":" + std::to_string(line_no) + ": family_id must have 4 dot fields");

        const std::string fold = fam.substr(0, d2);
        const std::string sf = fam.substr(0, d3);

        const uint32_t id = static_cast<uint32_t>(tax.dom_fold.size());
        tax.dom_index.emplace(dom, id);
        tax.id_to_dom.push_back(dom);
        tax.dom_fold.push_back(intern_string(fold, fold_pool, fold_strings));
        tax.dom_sf.push_back(intern_string(sf, sf_pool, sf_strings));
        tax.dom_fam.push_back(intern_string(fam, fam_pool, fam_strings));
    }
    tax.ndom_lookup = static_cast<uint32_t>(tax.dom_fold.size());
    return tax;
}

inline std::vector<std::string> json_string_array_after_key(const std::string& json, const std::string& key) {
    const std::string needle = "\"" + key + "\"";
    size_t pos = json.find(needle);
    if (pos == std::string::npos)
        die_tax("derived_info missing key " + key);
    pos = json.find('[', pos + needle.size());
    if (pos == std::string::npos)
        die_tax("derived_info malformed array for " + key);
    ++pos;

    std::vector<std::string> out;
    while (pos < json.size()) {
        while (pos < json.size() &&
               (json[pos] == ' ' || json[pos] == '\t' || json[pos] == ',' ||
                json[pos] == '\n' || json[pos] == '\r'))
            ++pos;
        if (pos >= json.size() || json[pos] == ']')
            break;
        if (json[pos] != '"')
            die_tax("derived_info expected string in " + key);
        ++pos;
        const size_t start = pos;
        while (pos < json.size() && json[pos] != '"')
            ++pos;
        out.push_back(json.substr(start, pos - start));
        if (pos < json.size())
            ++pos;
    }
    return out;
}

inline std::string json_string_after_key(const std::string& json, const std::string& key) {
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
    if (pos >= json.size() || json[pos] != '"')
        die_tax("derived_info expected string at " + key);
    ++pos;
    const size_t start = pos;
    while (pos < json.size() && json[pos] != '"')
        ++pos;
    return json.substr(start, pos - start);
}

inline std::string lookup_basename(const std::string& path) {
    std::string base = path;
    const size_t slash = base.find_last_of("/\\");
    if (slash != std::string::npos)
        base = base.substr(slash + 1);
    const size_t dot = base.rfind('.');
    if (dot != std::string::npos)
        base = base.substr(0, dot);
    return base;
}

inline Taxonomy load_taxonomy_from_derived_info(const std::string& path) {
    const std::string json = read_file(path);
    const std::vector<std::string> domain_ids = json_string_array_after_key(json, "domain_ids");
    const std::vector<std::string> folds = json_string_array_after_key(json, "dom_fold");
    const std::vector<std::string> sfs = json_string_array_after_key(json, "dom_sf");
    const std::vector<std::string> fams = json_string_array_after_key(json, "dom_fam");

    if (domain_ids.empty())
        die_tax("derived_info missing domain_ids");
    if (folds.size() != domain_ids.size() || sfs.size() != domain_ids.size() ||
        fams.size() != domain_ids.size())
        die_tax("derived_info domain annotation arrays must match domain_ids length");

    Taxonomy tax;
    std::unordered_map<std::string, uint32_t> fold_pool, sf_pool, fam_pool;
    std::vector<std::string> fold_strings, sf_strings, fam_strings;

    tax.dom_fold.reserve(domain_ids.size());
    tax.dom_sf.reserve(domain_ids.size());
    tax.dom_fam.reserve(domain_ids.size());
    tax.id_to_dom = domain_ids;
    for (size_t i = 0; i < domain_ids.size(); ++i) {
        const std::string& dom = domain_ids[i];
        if (tax.dom_index.count(dom))
            die_tax("derived_info duplicate domain " + dom);
        const uint32_t id = static_cast<uint32_t>(i);
        tax.dom_index.emplace(dom, id);
        tax.dom_fold.push_back(intern_string(folds[i], fold_pool, fold_strings));
        tax.dom_sf.push_back(intern_string(sfs[i], sf_pool, sf_strings));
        tax.dom_fam.push_back(intern_string(fams[i], fam_pool, fam_strings));
    }
    tax.ndom_lookup = static_cast<uint32_t>(domain_ids.size());
    return tax;
}

enum class Truth { Fold, Superfamily, Family, SuperfamilyX };

inline bool pair_ignore(uint32_t q, uint32_t t, Truth truth, const Taxonomy& tax) {
    if (q == t)
        return true;
    if (truth == Truth::SuperfamilyX) {
        if (tax.dom_fold[q] == tax.dom_fold[t] && tax.dom_sf[q] != tax.dom_sf[t])
            return true;
    }
    return false;
}

inline bool pair_is_tp(uint32_t q, uint32_t t, Truth truth, const Taxonomy& tax) {
    switch (truth) {
    case Truth::Fold:
        return tax.dom_fold[q] == tax.dom_fold[t];
    case Truth::Superfamily:
    case Truth::SuperfamilyX:
        return tax.dom_sf[q] == tax.dom_sf[t];
    case Truth::Family:
        return tax.dom_fam[q] == tax.dom_fam[t];
    }
    return false;
}

inline bool lookup_domain(const Taxonomy& tax, const std::string& label, uint32_t& out) {
    const auto it = tax.dom_index.find(label);
    if (it == tax.dom_index.end())
        return false;
    out = it->second;
    return true;
}

}  // namespace homval
