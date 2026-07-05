// hits_to_edf.cpp -- fast pair-list hits -> .edf summary (homval stage 1a).
//
// Two-pass streaming over the hits file:
//   pass 1: score histogram + per-query best FP score (for SFFP)
//   pass 2: count TPs strictly better than that query's best FP
// Multi-threaded chunked reads; domain strings interned to integer ids.

#include <algorithm>
#include <charconv>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <string_view>
#include <thread>
#include <unordered_map>
#include <utility>
#include <vector>

#include "hits_bin.hpp"

namespace {

enum class Truth { Fold, Superfamily, Family, SuperfamilyX };

enum class ScoreDir { HigherBetter, LowerBetter };

struct Config {
    std::string hits_path;
    std::string bin_path;
    bool use_bin = false;
    std::string lookup_path;
    std::string derived_info_path;
    std::string output_path;
    std::string truth_str;
    Truth truth = Truth::Superfamily;
    bool all_truths = false;
    std::string fields_spec = "1,2,3";
    int q_idx = 0;
    int t_idx = 1;
    int s_idx = 2;
    std::string algo;
    std::string reference;
    ScoreDir score_dir = ScoreDir::HigherBetter;
    unsigned threads = 0;
};

struct Counts { uint64_t n_tp = 0; uint64_t n_fp = 0; };

struct CurveRow {
    double score = 0.0;
    uint64_t n_tp = 0;
    uint64_t n_fp = 0;
    uint64_t cum_tp = 0;
    uint64_t cum_fp = 0;
};

struct QueryFp {
    bool has = false;
    double score = 0.0;
};

using HistMap = std::unordered_map<double, Counts>;
using FpMap = std::unordered_map<uint32_t, QueryFp>;

struct Taxonomy {
    std::unordered_map<std::string, uint32_t> dom_index;
    std::vector<uint32_t> dom_fold;
    std::vector<uint32_t> dom_sf;
    std::vector<uint32_t> dom_fam;
    uint32_t ndom_lookup = 0;
};

struct Denoms {
    int ndom = 0;
    uint64_t n_possible_tp = 0;
    uint64_t n_possible_fp = 0;
};

// ---------------------------------------------------------------------------
// Utilities
// ---------------------------------------------------------------------------

[[noreturn]] void die(const std::string& msg) {
    std::cerr << "error: " << msg << '\n';
    std::exit(1);
}

void trim_inplace(std::string& s) {
    while (!s.empty() && (s.back() == '\r' || s.back() == '\n' || s.back() == ' '))
        s.pop_back();
    size_t i = 0;
    while (i < s.size() && (s[i] == ' ' || s[i] == '\t'))
        ++i;
    if (i > 0)
        s.erase(0, i);
}

std::string read_file(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in)
        die("cannot open " + path);
    std::ostringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

bool better_score(double score, double best, ScoreDir dir) {
    if (dir == ScoreDir::LowerBetter)
        return score < best;
    return score > best;
}

bool parse_domain_label(const char* start, const char* end, std::string& out) {
    out.clear();
    while (start < end && *start != '/' && *start != '\t' && *start != ' ' && *start != '\r') {
        out.push_back(*start);
        ++start;
    }
    return !out.empty();
}

bool field_at(std::string_view line, int idx, std::string& out) {
    out.clear();
    int cur = 0;
    size_t i = 0;
    while (i < line.size()) {
        while (i < line.size() && (line[i] == ' ' || line[i] == '\t'))
            ++i;
        if (i >= line.size())
            break;
        size_t j = i;
        while (j < line.size() && line[j] != ' ' && line[j] != '\t')
            ++j;
        if (cur == idx) {
            out.assign(line.data() + i, j - i);
            return true;
        }
        ++cur;
        i = j;
    }
    return false;
}

bool parse_double_field(const std::string& s, double& out) {
    const char* begin = s.data();
    const char* end = begin + s.size();
    auto [ptr, ec] = std::from_chars(begin, end, out);
    if (ec != std::errc() || ptr != end)
        return false;
    return std::isfinite(out);
}

Truth parse_truth(const std::string& s) {
    if (s == "fold")
        return Truth::Fold;
    if (s == "superfamily")
        return Truth::Superfamily;
    if (s == "family")
        return Truth::Family;
    if (s == "superfamilyx")
        return Truth::SuperfamilyX;
    die("unknown truth standard: " + s);
}

const char* truth_name(Truth truth) {
    switch (truth) {
    case Truth::Fold:
        return "fold";
    case Truth::Superfamily:
        return "superfamily";
    case Truth::Family:
        return "family";
    case Truth::SuperfamilyX:
        return "superfamilyx";
    }
    return "unknown";
}

std::vector<Truth> all_truth_values() {
    return {Truth::Family, Truth::Superfamily, Truth::SuperfamilyX, Truth::Fold};
}

uint32_t intern_string(
    const std::string& s,
    std::unordered_map<std::string, uint32_t>& pool,
    std::vector<std::string>& strings) {
    auto it = pool.find(s);
    if (it != pool.end())
        return it->second;
    uint32_t id = static_cast<uint32_t>(strings.size());
    strings.push_back(s);
    pool.emplace(strings.back(), id);
    return id;
}

// ---------------------------------------------------------------------------
// Lookup + JSON
// ---------------------------------------------------------------------------

Taxonomy load_lookup(const std::string& path) {
    Taxonomy tax;
    std::ifstream in(path);
    if (!in)
        die("cannot open lookup " + path);

    std::unordered_map<std::string, uint32_t> fold_pool, sf_pool, fam_pool;
    std::vector<std::string> fold_strings, sf_strings, fam_strings;

    std::string line;
    int line_no = 0;
    while (std::getline(in, line)) {
        ++line_no;
        trim_inplace(line);
        if (line.empty())
            continue;
        size_t tab = line.find('\t');
        if (tab == std::string::npos)
            die(path + ":" + std::to_string(line_no) + ": expected 2 tab fields");

        std::string dom = line.substr(0, tab);
        std::string fam = line.substr(tab + 1);
        if (tax.dom_index.count(dom))
            die(path + ":" + std::to_string(line_no) + ": duplicate domain " + dom);

        size_t d1 = fam.find('.');
        size_t d2 = (d1 == std::string::npos) ? std::string::npos : fam.find('.', d1 + 1);
        size_t d3 = (d2 == std::string::npos) ? std::string::npos : fam.find('.', d2 + 1);
        if (d1 == std::string::npos || d2 == std::string::npos || d3 == std::string::npos)
            die(path + ":" + std::to_string(line_no) + ": family_id must have 4 dot fields");

        std::string fold = fam.substr(0, d2);
        std::string sf = fam.substr(0, d3);

        uint32_t id = static_cast<uint32_t>(tax.dom_fold.size());
        tax.dom_index.emplace(dom, id);
        tax.dom_fold.push_back(intern_string(fold, fold_pool, fold_strings));
        tax.dom_sf.push_back(intern_string(sf, sf_pool, sf_strings));
        tax.dom_fam.push_back(intern_string(fam, fam_pool, fam_strings));
    }
    tax.ndom_lookup = static_cast<uint32_t>(tax.dom_fold.size());
    return tax;
}

int64_t json_int_after_key(const std::string& json, const std::string& key) {
    const std::string needle = "\"" + key + "\"";
    size_t pos = json.find(needle);
    if (pos == std::string::npos)
        die("derived_info missing key " + key);
    pos = json.find(':', pos + needle.size());
    if (pos == std::string::npos)
        die("derived_info malformed at " + key);
    ++pos;
    while (pos < json.size() && (json[pos] == ' ' || json[pos] == '\t'))
        ++pos;
    size_t end = pos;
    while (end < json.size() && (json[end] == '-' || std::isdigit(static_cast<unsigned char>(json[end]))))
        ++end;
    return std::stoll(json.substr(pos, end - pos));
}

int64_t json_truth_int(const std::string& json, const std::string& section, const std::string& truth) {
    const std::string sec = "\"" + section + "\"";
    size_t pos = json.find(sec);
    if (pos == std::string::npos)
        die("derived_info missing " + section);
    pos = json.find('{', pos);
    if (pos == std::string::npos)
        die("derived_info malformed " + section);
    size_t end = json.find('}', pos);
    if (end == std::string::npos)
        die("derived_info malformed " + section);
    std::string block = json.substr(pos, end - pos + 1);
    const std::string key = "\"" + truth + "\"";
    size_t kpos = block.find(key);
    if (kpos == std::string::npos)
        die("derived_info missing " + section + " for truth " + truth);
    kpos = block.find(':', kpos + key.size());
    if (kpos == std::string::npos)
        die("derived_info malformed " + section + " for truth " + truth);
    ++kpos;
    while (kpos < block.size() && (block[kpos] == ' ' || block[kpos] == '\t'))
        ++kpos;
    size_t e = kpos;
    while (e < block.size() && (block[e] == '-' || std::isdigit(static_cast<unsigned char>(block[e]))))
        ++e;
    return std::stoll(block.substr(kpos, e - kpos));
}

Denoms load_denoms(const std::string& path, const std::string& truth) {
    const std::string json = read_file(path);
    Denoms d;
    d.ndom = static_cast<int>(json_int_after_key(json, "ndom"));
    d.n_possible_tp = static_cast<uint64_t>(json_truth_int(json, "N_possible_tp", truth));
    d.n_possible_fp = static_cast<uint64_t>(json_truth_int(json, "N_possible_fp", truth));
    return d;
}

// ---------------------------------------------------------------------------
// Pair classification
// ---------------------------------------------------------------------------

bool pair_ignore(uint32_t q, uint32_t t, Truth truth, const Taxonomy& tax) {
    if (q == t)
        return true;
    if (truth == Truth::SuperfamilyX) {
        if (tax.dom_fold[q] == tax.dom_fold[t] && tax.dom_sf[q] != tax.dom_sf[t])
            return true;
    }
    return false;
}

bool pair_is_tp(uint32_t q, uint32_t t, Truth truth, const Taxonomy& tax) {
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

bool lookup_domain(const Taxonomy& tax, const std::string& label, uint32_t& out) {
    auto it = tax.dom_index.find(label);
    if (it == tax.dom_index.end())
        return false;
    out = it->second;
    return true;
}

// ---------------------------------------------------------------------------
// Hit processing
// ---------------------------------------------------------------------------

struct Pass1Local {
    HistMap hist;
    FpMap best_fp;
    uint64_t ignored = 0;
    uint64_t considered = 0;
};

struct ChunkJob {
    size_t start = 0;
    size_t end = 0;
    size_t file_size = 0;
};

std::vector<ChunkJob> make_chunks(const std::string& path, unsigned nthreads) {
    std::ifstream in(path, std::ios::binary | std::ios::ate);
    if (!in)
        die("cannot open hits " + path);
    const size_t size = static_cast<size_t>(in.tellg());
    if (size == 0) {
        ChunkJob c;
        c.file_size = 0;
        return {c};
    }

    std::vector<size_t> bounds;
    bounds.reserve(nthreads + 1);
    bounds.push_back(0);
    for (unsigned i = 1; i < nthreads; ++i) {
        const size_t raw = (size * i) / nthreads;
        in.seekg(static_cast<std::streamoff>(raw));
        std::string skip;
        if (raw > 0)
            std::getline(in, skip);
        bounds.push_back(static_cast<size_t>(in.tellg()));
    }
    bounds.push_back(size);

    std::vector<ChunkJob> chunks;
    chunks.reserve(nthreads);
    for (unsigned i = 0; i < nthreads; ++i) {
        ChunkJob c;
        c.file_size = size;
        c.start = bounds[i];
        c.end = bounds[i + 1];
        if (c.start < c.end)
            chunks.push_back(c);
    }
    if (chunks.empty()) {
        ChunkJob c;
        c.file_size = size;
        c.start = 0;
        c.end = size;
        chunks.push_back(c);
    }
    return chunks;
}

bool read_chunk_line(std::ifstream& in, const ChunkJob& chunk, std::string& line) {
    const std::streampos pos = in.tellg();
    const bool is_last = (chunk.end == chunk.file_size);
    if (!is_last && pos >= 0 && static_cast<size_t>(pos) >= chunk.end)
        return false;
    return static_cast<bool>(std::getline(in, line));
}

void pass1_chunk(
    const Config& cfg,
    const Taxonomy& tax,
    const ChunkJob& chunk,
    Pass1Local& local) {
    std::ifstream in(cfg.hits_path, std::ios::binary);
    if (!in)
        die("cannot open hits " + cfg.hits_path);
    in.seekg(static_cast<std::streamoff>(chunk.start));

    std::string line;
    std::string fq, ft, fs;
    while (read_chunk_line(in, chunk, line)) {
        trim_inplace(line);
        if (line.empty() || line[0] == '#')
            continue;

        if (!field_at(line, cfg.q_idx, fq) || !field_at(line, cfg.t_idx, ft) || !field_at(line, cfg.s_idx, fs))
            continue;

        std::string qdom, tdom;
        if (!parse_domain_label(fq.data(), fq.data() + fq.size(), qdom) ||
            !parse_domain_label(ft.data(), ft.data() + ft.size(), tdom))
            continue;

        uint32_t qid = 0, tid = 0;
        if (!lookup_domain(tax, qdom, qid) || !lookup_domain(tax, tdom, tid)) {
            ++local.ignored;
            continue;
        }

        double score = 0.0;
        if (!parse_double_field(fs, score))
            continue;

        if (pair_ignore(qid, tid, cfg.truth, tax)) {
            ++local.ignored;
            continue;
        }

        ++local.considered;
        const bool is_tp = pair_is_tp(qid, tid, cfg.truth, tax);
        if (is_tp)
            local.hist[score].n_tp += 1;
        else {
            local.hist[score].n_fp += 1;
            QueryFp& qf = local.best_fp[qid];
            if (!qf.has || better_score(score, qf.score, cfg.score_dir)) {
                qf.has = true;
                qf.score = score;
            }
        }
    }
}

struct Pass2Local { uint64_t tp_above = 0; };

struct RecordChunk {
    size_t begin = 0;
    size_t end = 0;
};

std::vector<RecordChunk> make_record_chunks(size_t n_records, unsigned nthreads) {
    if (n_records == 0 || nthreads == 0) {
        RecordChunk c;
        c.end = n_records;
        return {c};
    }
    const size_t chunk_size = (n_records + nthreads - 1) / nthreads;
    std::vector<RecordChunk> chunks;
    chunks.reserve(nthreads);
    for (size_t begin = 0; begin < n_records; begin += chunk_size) {
        RecordChunk c;
        c.begin = begin;
        c.end = std::min(n_records, begin + chunk_size);
        if (c.begin < c.end)
            chunks.push_back(c);
    }
    if (chunks.empty()) {
        RecordChunk c;
        c.end = n_records;
        chunks.push_back(c);
    }
    return chunks;
}

void pass1_bin_range(
    const std::vector<homval::HitRecord>& hits,
    const RecordChunk& chunk,
    Truth truth,
    const Taxonomy& tax,
    ScoreDir dir,
    Pass1Local& local) {
    for (size_t i = chunk.begin; i < chunk.end; ++i) {
        const homval::HitRecord& h = hits[i];
        const uint32_t qid = h.q;
        const uint32_t tid = h.t;
        const double score = h.score;

        if (pair_ignore(qid, tid, truth, tax)) {
            ++local.ignored;
            continue;
        }

        ++local.considered;
        if (pair_is_tp(qid, tid, truth, tax))
            local.hist[score].n_tp += 1;
        else {
            local.hist[score].n_fp += 1;
            QueryFp& qf = local.best_fp[qid];
            if (!qf.has || better_score(score, qf.score, dir)) {
                qf.has = true;
                qf.score = score;
            }
        }
    }
}

void pass2_bin_range(
    const std::vector<homval::HitRecord>& hits,
    const RecordChunk& chunk,
    Truth truth,
    const Taxonomy& tax,
    const FpMap& global_fp,
    ScoreDir dir,
    Pass2Local& local) {
    for (size_t i = chunk.begin; i < chunk.end; ++i) {
        const homval::HitRecord& h = hits[i];
        const uint32_t qid = h.q;
        const uint32_t tid = h.t;
        const double score = h.score;

        if (pair_ignore(qid, tid, truth, tax))
            continue;
        if (!pair_is_tp(qid, tid, truth, tax))
            continue;

        auto it = global_fp.find(qid);
        if (it == global_fp.end() || !it->second.has)
            ++local.tp_above;
        else if (better_score(score, it->second.score, dir))
            ++local.tp_above;
    }
}

void pass2_chunk(
    const Config& cfg,
    const Taxonomy& tax,
    const FpMap& global_fp,
    const ChunkJob& chunk,
    Pass2Local& local) {
    std::ifstream in(cfg.hits_path, std::ios::binary);
    if (!in)
        die("cannot open hits " + cfg.hits_path);
    in.seekg(static_cast<std::streamoff>(chunk.start));

    std::string line;
    std::string fq, ft, fs;
    while (read_chunk_line(in, chunk, line)) {
        trim_inplace(line);
        if (line.empty() || line[0] == '#')
            continue;

        if (!field_at(line, cfg.q_idx, fq) || !field_at(line, cfg.t_idx, ft) || !field_at(line, cfg.s_idx, fs))
            continue;

        std::string qdom, tdom;
        if (!parse_domain_label(fq.data(), fq.data() + fq.size(), qdom) ||
            !parse_domain_label(ft.data(), ft.data() + ft.size(), tdom))
            continue;

        uint32_t qid = 0, tid = 0;
        if (!lookup_domain(tax, qdom, qid) || !lookup_domain(tax, tdom, tid))
            continue;

        double score = 0.0;
        if (!parse_double_field(fs, score))
            continue;

        if (pair_ignore(qid, tid, cfg.truth, tax))
            continue;

        if (!pair_is_tp(qid, tid, cfg.truth, tax))
            continue;

        auto it = global_fp.find(qid);
        if (it == global_fp.end() || !it->second.has)
            ++local.tp_above;
        else if (better_score(score, it->second.score, cfg.score_dir))
            ++local.tp_above;
    }
}

void merge_hist(HistMap& dst, const HistMap& src) {
    for (const auto& [score, c] : src) {
        dst[score].n_tp += c.n_tp;
        dst[score].n_fp += c.n_fp;
    }
}

void merge_fp(FpMap& dst, const FpMap& src, ScoreDir dir) {
    for (const auto& [qid, v] : src) {
        if (!v.has)
            continue;
        QueryFp& d = dst[qid];
        if (!d.has || better_score(v.score, d.score, dir)) {
            d = v;
        }
    }
}

// ---------------------------------------------------------------------------
// Summary statistics
// ---------------------------------------------------------------------------

void sort_scores(std::vector<double>& scores, ScoreDir dir) {
    if (dir == ScoreDir::LowerBetter)
        std::sort(scores.begin(), scores.end());
    else
        std::sort(scores.begin(), scores.end(), std::greater<double>());
}

std::vector<CurveRow> build_curve(const HistMap& hist, ScoreDir dir) {
    std::vector<double> scores;
    scores.reserve(hist.size());
    for (const auto& [s, _] : hist)
        scores.push_back(s);
    sort_scores(scores, dir);

    std::vector<CurveRow> rows;
    rows.reserve(scores.size());
    uint64_t cum_tp = 0, cum_fp = 0;
    for (double s : scores) {
        const Counts& c = hist.at(s);
        cum_tp += c.n_tp;
        cum_fp += c.n_fp;
        rows.push_back({s, c.n_tp, c.n_fp, cum_tp, cum_fp});
    }
    return rows;
}

struct Sum3Stats {
    double sepq01 = 0.0;
    double sepq1 = 0.0;
    double sepq10 = 0.0;
    double sum3 = 0.0;
};

Sum3Stats compute_sum3(const std::vector<CurveRow>& rows, int ndom, uint64_t n_possible_tp) {
    Sum3Stats out;
    if (n_possible_tp == 0 || rows.empty())
        return out;

    bool have01 = false, have1 = false, have10 = false;
    bool have_last = false;
    double last_score = 0.0;
    for (const CurveRow& r : rows) {
        if (have_last && r.score != last_score) {
            const double cov = static_cast<double>(r.cum_tp) / static_cast<double>(n_possible_tp);
            const double err = static_cast<double>(r.cum_fp) / static_cast<double>(ndom);
            if (err <= 0.1) {
                out.sepq01 = cov;
                have01 = true;
            }
            if (err <= 1.0) {
                out.sepq1 = cov;
                have1 = true;
            }
            if (err <= 10.0) {
                out.sepq10 = cov;
                have10 = true;
            }
        }
        have_last = true;
        last_score = r.score;
    }

    const double final_cov =
        static_cast<double>(rows.back().cum_tp) / static_cast<double>(n_possible_tp);
    if (!have01)
        out.sepq01 = final_cov;
    if (!have1)
        out.sepq1 = final_cov;
    if (!have10)
        out.sepq10 = final_cov;
    out.sum3 = out.sepq01 * 2.0 + out.sepq1 * 1.5 + out.sepq10;
    return out;
}

double compute_sffp(uint64_t tp_above, uint64_t n_possible_tp) {
    if (n_possible_tp == 0)
        return 0.0;
    return static_cast<double>(tp_above) / static_cast<double>(n_possible_tp);
}

double compute_pr90(const std::vector<CurveRow>& rows, uint64_t n_possible_tp) {
    if (n_possible_tp == 0 || rows.empty())
        return 0.0;
    double best = 0.0;
    for (const CurveRow& r : rows) {
        const uint64_t denom = r.cum_tp + r.cum_fp;
        if (denom == 0)
            continue;
        const double prec = static_cast<double>(r.cum_tp) / static_cast<double>(denom);
        if (prec >= 0.9) {
            const double recall = static_cast<double>(r.cum_tp) / static_cast<double>(n_possible_tp);
            if (recall > best)
                best = recall;
        }
    }
    return best;
}

void format_g(FILE* f, double v) {
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%g", v);
    std::fputs(buf, f);
}

void write_edf(
    const Config& cfg,
    const std::string& truth_str,
    const Denoms& denoms,
    const HistMap& global_hist,
    uint64_t tp_above) {
    const std::vector<CurveRow> curve = build_curve(global_hist, cfg.score_dir);
    const Sum3Stats sum3 = compute_sum3(curve, denoms.ndom, denoms.n_possible_tp);
    const double sffp = compute_sffp(tp_above, denoms.n_possible_tp);
    const double pr90 = compute_pr90(curve, denoms.n_possible_tp);

    std::cerr << "writing " << cfg.output_path << " (" << global_hist.size()
              << " distinct scores)...\n";

    FILE* out = std::fopen(cfg.output_path.c_str(), "w");
    if (!out)
        die("cannot open output " + cfg.output_path);

    std::fprintf(out, "# algo=%s\n", cfg.algo.c_str());
    if (!cfg.reference.empty())
        std::fprintf(out, "# reference=%s\n", cfg.reference.c_str());
    std::fprintf(out, "# truth=%s\n", truth_str.c_str());
    std::fprintf(out, "# fields=%s\n", cfg.fields_spec.c_str());
    std::fprintf(out, "# score_direction=%s\n",
                 cfg.score_dir == ScoreDir::LowerBetter ? "lower_better" : "higher_better");
    std::fprintf(out, "# ndom=%d\n", denoms.ndom);
    std::fprintf(out, "# N_possible_tp=%llu\n", static_cast<unsigned long long>(denoms.n_possible_tp));
    std::fprintf(out, "# N_possible_fp=%llu\n", static_cast<unsigned long long>(denoms.n_possible_fp));

    std::fprintf(out, "# SEPQ0.1=%.3f SEPQ1=%.3f SEPQ10=%.3f Sum3=%.3f SFFP=%.3f PR90=%.4g\n",
                 sum3.sepq01, sum3.sepq1, sum3.sepq10, sum3.sum3, sffp, pr90);
    std::fprintf(stderr, "# truth=%s SEPQ0.1=%.3f SEPQ1=%.3f SEPQ10=%.3f Sum3=%.3f SFFP=%.3f PR90=%.4g\n",
                 truth_str.c_str(), sum3.sepq01, sum3.sepq1, sum3.sepq10, sum3.sum3, sffp, pr90);

    std::fputs("score\tn_tp\tn_fp\tcum_tp\tcum_fp\n", out);
    for (const CurveRow& r : curve) {
        format_g(out, r.score);
        std::fprintf(out, "\t%llu\t%llu\t%llu\t%llu\n",
                     static_cast<unsigned long long>(r.n_tp),
                     static_cast<unsigned long long>(r.n_fp),
                     static_cast<unsigned long long>(r.cum_tp),
                     static_cast<unsigned long long>(r.cum_fp));
    }

    std::fclose(out);
}

unsigned effective_threads_records(size_t n_records, unsigned requested) {
    const unsigned max_by_records =
        static_cast<unsigned>(std::max<size_t>(1, n_records / 50000));
    return std::max(1u, std::min(requested, max_by_records));
}

void run_bin_passes(
    const Config& cfg,
    const Taxonomy& tax,
    const std::vector<homval::HitRecord>& hits,
    Truth truth,
    HistMap& global_hist,
    FpMap& global_fp,
    uint64_t& tp_above,
    uint64_t& ignored,
    uint64_t& considered) {
    const unsigned nworkers = effective_threads_records(hits.size(), cfg.threads);
    const std::vector<RecordChunk> chunks = make_record_chunks(hits.size(), nworkers);

    global_hist.clear();
    global_fp.clear();
    ignored = 0;
    considered = 0;
    tp_above = 0;

    {
        std::vector<Pass1Local> locals(nworkers);
        std::vector<std::thread> workers;
        workers.reserve(nworkers);
        for (unsigned t = 0; t < nworkers; ++t) {
            workers.emplace_back([&, t]() {
                pass1_bin_range(hits, chunks[t], truth, tax, cfg.score_dir, locals[t]);
            });
        }
        for (auto& w : workers)
            w.join();

        for (const Pass1Local& loc : locals) {
            merge_hist(global_hist, loc.hist);
            merge_fp(global_fp, loc.best_fp, cfg.score_dir);
            ignored += loc.ignored;
            considered += loc.considered;
        }
    }

    {
        std::vector<Pass2Local> locals(nworkers);
        std::vector<std::thread> workers;
        workers.reserve(nworkers);
        for (unsigned t = 0; t < nworkers; ++t) {
            workers.emplace_back([&, t]() {
                pass2_bin_range(hits, chunks[t], truth, tax, global_fp, cfg.score_dir, locals[t]);
            });
        }
        for (auto& w : workers)
            w.join();
        for (const Pass2Local& loc : locals)
            tp_above += loc.tp_above;
    }
}

// ---------------------------------------------------------------------------
// CLI
// ---------------------------------------------------------------------------

void parse_fields(const std::string& spec, int& q, int& t, int& s) {
    int q1 = 0, t1 = 0, s1 = 0;
    if (std::sscanf(spec.c_str(), "%d,%d,%d", &q1, &t1, &s1) != 3 || q1 < 1 || t1 < 1 || s1 < 1)
        die("fields must be three 1-based indices, e.g. 1,2,3");
    q = q1 - 1;
    t = t1 - 1;
    s = s1 - 1;
}

std::pair<std::string, std::string> parse_algo_reference(const std::string& hits_path) {
    std::string base = hits_path;
    const size_t slash = base.find_last_of("/\\");
    if (slash != std::string::npos)
        base = base.substr(slash + 1);
    for (const char* ext : {".tsv", ".txt"}) {
        const size_t n = std::strlen(ext);
        if (base.size() >= n && base.compare(base.size() - n, n, ext) == 0) {
            base.resize(base.size() - n);
            break;
        }
    }
    const size_t dot = base.find('.');
    if (dot == std::string::npos)
        return {base, ""};
    return {base.substr(0, dot), base.substr(dot + 1)};
}

Config parse_args(int argc, char** argv) {
    Config cfg;
    bool have_evalue = false;
    bool have_score = false;

    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        auto need = [&](const char* name) -> std::string {
            if (i + 1 >= argc)
                die(std::string("missing value for ") + name);
            return argv[++i];
        };
        if (arg == "--hits")
            cfg.hits_path = need("--hits");
        else if (arg == "--bin") {
            cfg.bin_path = need("--bin");
            cfg.use_bin = true;
        }
        else if (arg == "--lookup")
            cfg.lookup_path = need("--lookup");
        else if (arg == "--derived-info")
            cfg.derived_info_path = need("--derived-info");
        else if (arg == "--truth")
            cfg.truth_str = need("--truth");
        else if (arg == "--all-truths")
            cfg.all_truths = true;
        else if (arg == "--output")
            cfg.output_path = need("--output");
        else if (arg == "--fields")
            cfg.fields_spec = need("--fields");
        else if (arg == "--algo")
            cfg.algo = need("--algo");
        else if (arg == "--reference")
            cfg.reference = need("--reference");
        else if (arg == "--threads") {
            const long n = std::stol(need("--threads"));
            if (n < 1)
                die("--threads must be >= 1");
            cfg.threads = static_cast<unsigned>(n);
        }
        else if (arg == "--evalue")
            have_evalue = true;
        else if (arg == "--score")
            have_score = true;
        else if (arg == "--help" || arg == "-h") {
            std::cout <<
                "Usage: hits_to_edf (--hits PATH | --bin PATH.bin) --lookup PATH --derived-info PATH\n"
                "       (--truth {fold,superfamily,family,superfamilyx} | --all-truths) --output PATH\n"
                "       (--evalue | --score) [--fields 1,2,3] [--algo NAME] [--reference DB]\n"
                "       [--threads N]\n"
                "\n"
                "With --all-truths, --output is a prefix; writes PREFIX.fold, PREFIX.superfamily, etc.\n";
            std::exit(0);
        }
        else
            die("unknown argument: " + arg);
    }

    if (cfg.use_bin) {
        if (!cfg.hits_path.empty())
            die("use either --hits or --bin, not both");
        if (cfg.lookup_path.empty() || cfg.derived_info_path.empty() || cfg.output_path.empty())
            die("required: --bin --lookup --derived-info --output and (--evalue|--score)");
        if (cfg.all_truths) {
            if (!cfg.truth_str.empty())
                die("use either --truth or --all-truths, not both");
        } else if (cfg.truth_str.empty())
            die("required: --truth or --all-truths");
    } else {
        if (!cfg.bin_path.empty())
            die("internal: bin path set without --bin");
        if (cfg.hits_path.empty() || cfg.lookup_path.empty() || cfg.derived_info_path.empty() ||
            cfg.output_path.empty())
            die("required: --hits --lookup --derived-info --truth --output and (--evalue|--score)");
        if (cfg.all_truths)
            die("--all-truths requires --bin");
        if (cfg.truth_str.empty())
            die("required: --truth");
    }
    if (have_evalue == have_score)
        die("exactly one of --evalue or --score is required");

    if (!cfg.truth_str.empty())
        cfg.truth = parse_truth(cfg.truth_str);
    cfg.score_dir = have_evalue ? ScoreDir::LowerBetter : ScoreDir::HigherBetter;
    if (!cfg.use_bin)
        parse_fields(cfg.fields_spec, cfg.q_idx, cfg.t_idx, cfg.s_idx);

    if (!cfg.use_bin) {
        auto [def_algo, def_ref] = parse_algo_reference(cfg.hits_path);
        if (cfg.algo.empty())
            cfg.algo = def_algo;
        if (cfg.reference.empty())
            cfg.reference = def_ref;
    } else if (cfg.algo.empty())
        die("--bin mode requires --algo");
    if (cfg.threads == 0)
        cfg.threads = std::max(1u, std::thread::hardware_concurrency());
    return cfg;
}

unsigned effective_threads(unsigned requested, size_t file_size) {
    // Avoid 100+ threads on tiny inputs; ~2 MiB per thread minimum.
    const unsigned max_by_size =
        static_cast<unsigned>(std::max<size_t>(1, file_size / (2 * 1024 * 1024)));
    return std::max(1u, std::min(requested, max_by_size));
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------

int run_bin(const Config& cfg) {
    std::cerr << "\n\nloading lookup " << cfg.lookup_path << "...\n";
    const Taxonomy tax = load_lookup(cfg.lookup_path);

    std::cerr << "loading bin " << cfg.bin_path << "...\n";
    const std::vector<homval::HitRecord> hits = homval::load_hbin(cfg.bin_path, tax.ndom_lookup);
    std::cerr << "bin records: " << hits.size() << '\n';

    const std::vector<Truth> truths =
        cfg.all_truths ? all_truth_values() : std::vector<Truth>{cfg.truth};

    for (Truth truth : truths) {
        const std::string truth_str = truth_name(truth);
        const Denoms denoms = load_denoms(cfg.derived_info_path, truth_str);

        HistMap global_hist;
        FpMap global_fp;
        uint64_t ignored = 0;
        uint64_t considered = 0;
        uint64_t tp_above = 0;
        run_bin_passes(cfg, tax, hits, truth, global_hist, global_fp, tp_above, ignored, considered);

        std::cerr << "truth=" << truth_str << " ignored " << ignored
                  << " hits, considered " << considered << '\n';

        Config out_cfg = cfg;
        if (cfg.all_truths)
            out_cfg.output_path = cfg.output_path + "." + truth_str;
        write_edf(out_cfg, truth_str, denoms, global_hist, tp_above);
    }
    return 0;
}

int run_text(const Config& cfg) {
    const Denoms denoms = load_denoms(cfg.derived_info_path, cfg.truth_str);

    std::ifstream size_probe(cfg.hits_path, std::ios::binary | std::ios::ate);
    if (!size_probe)
        die("cannot open hits " + cfg.hits_path);
    const size_t hits_size = static_cast<size_t>(size_probe.tellg());
    Config run_cfg = cfg;
    run_cfg.threads = effective_threads(cfg.threads, hits_size);

    std::cerr << "\n\nloading lookup " << run_cfg.lookup_path << "...\n";
    const Taxonomy tax = load_lookup(run_cfg.lookup_path);

    const std::vector<ChunkJob> chunks = make_chunks(run_cfg.hits_path, run_cfg.threads);
    const unsigned nworkers = static_cast<unsigned>(chunks.size());

    std::cerr << "reading hits " << run_cfg.hits_path << " truth=" << run_cfg.truth_str
              << " (" << nworkers << " threads)...\n";

    HistMap global_hist;
    FpMap global_fp;
    uint64_t ignored = 0, considered = 0;

    {
        std::vector<Pass1Local> locals(nworkers);
        std::vector<std::thread> workers;
        workers.reserve(nworkers);
        for (unsigned t = 0; t < nworkers; ++t) {
            workers.emplace_back([&, t]() { pass1_chunk(run_cfg, tax, chunks[t], locals[t]); });
        }
        for (auto& w : workers)
            w.join();

        for (const Pass1Local& loc : locals) {
            merge_hist(global_hist, loc.hist);
            merge_fp(global_fp, loc.best_fp, cfg.score_dir);
            ignored += loc.ignored;
            considered += loc.considered;
        }
    }

    std::cerr << "ignored " << ignored << " hits, considered " << considered << '\n';

    uint64_t tp_above = 0;
    {
        std::vector<Pass2Local> locals(nworkers);
        std::vector<std::thread> workers;
        workers.reserve(nworkers);
        for (unsigned t = 0; t < nworkers; ++t) {
            workers.emplace_back([&, t]() {
                pass2_chunk(run_cfg, tax, global_fp, chunks[t], locals[t]);
            });
        }
        for (auto& w : workers)
            w.join();
        for (const Pass2Local& loc : locals)
            tp_above += loc.tp_above;
    }

    write_edf(run_cfg, run_cfg.truth_str, denoms, global_hist, tp_above);
    return 0;
}

int run(int argc, char** argv) {
    const Config cfg = parse_args(argc, argv);
    if (cfg.use_bin)
        return run_bin(cfg);
    return run_text(cfg);
}

}  // namespace

int main(int argc, char** argv) {
    return run(argc, argv);
}
