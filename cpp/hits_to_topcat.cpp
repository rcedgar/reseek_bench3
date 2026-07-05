// hits_to_topcat.cpp -- condense hits to leave-category-out .tcat summary (stage 1b).
//
// Single pass over text hits or binary cache; tracks per-query S_tp and S_xf.

#include <algorithm>
#include <charconv>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <string_view>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "hits_bin.hpp"

namespace {

enum class TopcatTruth { TopSf, TopFold };

enum class ScoreDir { HigherBetter, LowerBetter };

struct Config {
    std::string hits_path;
    std::string bin_path;
    bool use_bin = false;
    std::string lookup_path;
    std::string derived_info_path;
    std::string output_path;
    std::string truth_str;
    TopcatTruth truth = TopcatTruth::TopSf;
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

struct QueryScore {
    bool has_tp = false;
    bool has_xf = false;
    double s_tp = 0.0;
    double s_xf = 0.0;
};

using QueryScoreMap = std::unordered_map<uint32_t, QueryScore>;

struct Taxonomy {
    std::unordered_map<std::string, uint32_t> dom_index;
    std::vector<std::string> id_to_dom;
    std::vector<uint32_t> dom_fold;
    std::vector<uint32_t> dom_sf;
    std::vector<uint32_t> dom_fam;
    uint32_t ndom_lookup = 0;
};

struct QueryState {
    QueryScoreMap scores;
    std::vector<uint32_t> query_order;
    int n_possible_tp = 0;
};

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

TopcatTruth parse_topcat_truth(const std::string& s) {
    if (s == "topsf")
        return TopcatTruth::TopSf;
    if (s == "topfold")
        return TopcatTruth::TopFold;
    die("unknown topcat truth: " + s);
}

const char* topcat_truth_name(TopcatTruth truth) {
    switch (truth) {
    case TopcatTruth::TopSf:
        return "topsf";
    case TopcatTruth::TopFold:
        return "topfold";
    }
    return "unknown";
}

std::vector<TopcatTruth> all_topcat_truths() {
    return {TopcatTruth::TopSf, TopcatTruth::TopFold};
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
        const size_t tab = line.find('\t');
        if (tab == std::string::npos)
            die(path + ":" + std::to_string(line_no) + ": expected 2 tab fields");

        const std::string dom = line.substr(0, tab);
        const std::string fam = line.substr(tab + 1);
        if (tax.dom_index.count(dom))
            die(path + ":" + std::to_string(line_no) + ": duplicate domain " + dom);

        const size_t d1 = fam.find('.');
        const size_t d2 = (d1 == std::string::npos) ? std::string::npos : fam.find('.', d1 + 1);
        const size_t d3 = (d2 == std::string::npos) ? std::string::npos : fam.find('.', d2 + 1);
        if (d1 == std::string::npos || d2 == std::string::npos || d3 == std::string::npos)
            die(path + ":" + std::to_string(line_no) + ": family_id must have 4 dot fields");

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

int64_t json_truth_int(const std::string& json, const std::string& section, const std::string& truth) {
    const std::string sec = "\"" + section + "\"";
    size_t pos = json.find(sec);
    if (pos == std::string::npos)
        die("derived_info missing " + section);
    pos = json.find('{', pos);
    if (pos == std::string::npos)
        die("derived_info malformed " + section);
    const size_t end = json.find('}', pos);
    if (end == std::string::npos)
        die("derived_info malformed " + section);
    const std::string block = json.substr(pos, end - pos + 1);
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

std::vector<std::string> json_string_array_after_key(const std::string& json, const std::string& key) {
    const std::string needle = "\"" + key + "\"";
    size_t pos = json.find(needle);
    if (pos == std::string::npos)
        die("derived_info missing key " + key);
    pos = json.find('[', pos + needle.size());
    if (pos == std::string::npos)
        die("derived_info malformed array for " + key);
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
            die("derived_info expected string in " + key);
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

QueryState load_query_state(
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
            die("query domain not in lookup: " + dom);
        state.query_order.push_back(it->second);
        state.scores[it->second] = QueryScore{};
    }
    return state;
}

enum class Track { Ignore, Tp, Xf };

Track topcat_track(uint32_t q, uint32_t t, TopcatTruth truth, const Taxonomy& tax) {
    if (q == t)
        return Track::Ignore;
    if (truth == TopcatTruth::TopSf) {
        if (tax.dom_fam[q] == tax.dom_fam[t])
            return Track::Ignore;
        if (tax.dom_sf[q] == tax.dom_sf[t])
            return Track::Tp;
        return Track::Xf;
    }
    if (tax.dom_sf[q] == tax.dom_sf[t])
        return Track::Ignore;
    if (tax.dom_fold[q] == tax.dom_fold[t])
        return Track::Tp;
    return Track::Xf;
}

bool passes_threshold(const QueryScore& qs, bool use_tp, double threshold, ScoreDir dir) {
    if (use_tp) {
        if (!qs.has_tp)
            return false;
        if (dir == ScoreDir::LowerBetter)
            return qs.s_tp <= threshold;
        return qs.s_tp >= threshold;
    }
    if (!qs.has_xf)
        return false;
    if (dir == ScoreDir::LowerBetter)
        return qs.s_xf <= threshold;
    return qs.s_xf >= threshold;
}

bool at_least_as_good(double score, double other, ScoreDir dir) {
    if (dir == ScoreDir::LowerBetter)
        return score <= other;
    return score >= other;
}

bool is_tp_at_threshold(
    const QueryScore& qs,
    double threshold,
    ScoreDir dir) {
    if (!passes_threshold(qs, true, threshold, dir))
        return false;
    if (!passes_threshold(qs, false, threshold, dir))
        return true;
    return at_least_as_good(qs.s_tp, qs.s_xf, dir);
}

void classify_at_threshold(
    const QueryState& state,
    double threshold,
    ScoreDir dir,
    double& coverage,
    double& error) {
    if (state.n_possible_tp <= 0) {
        coverage = 0.0;
        error = 0.0;
        return;
    }

    int n_tp = 0;
    int n_fp_neg = 0;
    for (uint32_t qid : state.query_order) {
        const QueryScore& qs = state.scores.at(qid);
        if (is_tp_at_threshold(qs, threshold, dir))
            ++n_tp;
        if (passes_threshold(qs, false, threshold, dir))
            ++n_fp_neg;
    }
    const double n_q = static_cast<double>(state.n_possible_tp);
    coverage = static_cast<double>(n_tp) / n_q;
    error = static_cast<double>(n_fp_neg) / n_q;
}

std::vector<double> collect_thresholds(const QueryState& state) {
    std::unordered_set<double> seen;
    for (uint32_t qid : state.query_order) {
        const QueryScore& qs = state.scores.at(qid);
        if (qs.has_tp)
            seen.insert(qs.s_tp);
        if (qs.has_xf)
            seen.insert(qs.s_xf);
    }
    return std::vector<double>(seen.begin(), seen.end());
}

double loose_threshold(ScoreDir dir) {
    return dir == ScoreDir::LowerBetter
        ? std::numeric_limits<double>::infinity()
        : -std::numeric_limits<double>::infinity();
}

std::vector<double> sort_thresholds_strict_first(
    std::vector<double> thresholds,
    ScoreDir dir) {
    if (dir == ScoreDir::LowerBetter)
        std::sort(thresholds.begin(), thresholds.end());
    else
        std::sort(thresholds.begin(), thresholds.end(), std::greater<double>());
    return thresholds;
}

struct Top3Stats {
    double tepq0001 = 0.0;
    double tepq001 = 0.0;
    double tepq01 = 0.0;
    double top3 = 0.0;
};

Top3Stats compute_top3(const QueryState& state, ScoreDir dir) {
    Top3Stats out;
    if (state.n_possible_tp <= 0)
        return out;

    std::vector<double> taus = sort_thresholds_strict_first(collect_thresholds(state), dir);
    taus.push_back(loose_threshold(dir));

    bool have0001 = false;
    bool have001 = false;
    bool have01 = false;
    double final_cov = 0.0;

    for (double tau : taus) {
        double cov = 0.0;
        double err = 0.0;
        classify_at_threshold(state, tau, dir, cov, err);
        final_cov = cov;
        if (err >= 0.001 && !have0001) {
            out.tepq0001 = cov;
            have0001 = true;
        }
        if (err >= 0.01 && !have001) {
            out.tepq001 = cov;
            have001 = true;
        }
        if (err >= 0.1 && !have01) {
            out.tepq01 = cov;
            have01 = true;
        }
    }

    if (!have0001)
        out.tepq0001 = final_cov;
    if (!have001)
        out.tepq001 = final_cov;
    if (!have01)
        out.tepq01 = final_cov;
    out.top3 = out.tepq0001 * 2.0 + out.tepq001 * 1.5 + out.tepq01;
    return out;
}

void merge_scores(QueryScoreMap& dst, const QueryScoreMap& src, ScoreDir dir) {
    for (const auto& [qid, s] : src) {
        QueryScore& d = dst[qid];
        if (s.has_tp && (!d.has_tp || better_score(s.s_tp, d.s_tp, dir))) {
            d.has_tp = true;
            d.s_tp = s.s_tp;
        }
        if (s.has_xf && (!d.has_xf || better_score(s.s_xf, d.s_xf, dir))) {
            d.has_xf = true;
            d.s_xf = s.s_xf;
        }
    }
}

struct ScanLocal {
    QueryScoreMap scores;
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

void scan_text_chunk(
    const Config& cfg,
    const Taxonomy& tax,
    TopcatTruth truth,
    const std::unordered_set<uint32_t>& query_ids,
    const ChunkJob& chunk,
    ScanLocal& local) {
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

        auto iq = tax.dom_index.find(qdom);
        if (iq == tax.dom_index.end())
            continue;
        const uint32_t qid = iq->second;
        if (!query_ids.count(qid))
            continue;

        auto it = tax.dom_index.find(tdom);
        if (it == tax.dom_index.end())
            continue;
        const uint32_t tid = it->second;

        double score = 0.0;
        if (!parse_double_field(fs, score))
            continue;

        const Track track = topcat_track(qid, tid, truth, tax);
        if (track == Track::Ignore) {
            ++local.ignored;
            continue;
        }

        ++local.considered;
        QueryScore& qs = local.scores[qid];
        if (track == Track::Tp) {
            if (!qs.has_tp || better_score(score, qs.s_tp, cfg.score_dir)) {
                qs.has_tp = true;
                qs.s_tp = score;
            }
        } else if (!qs.has_xf || better_score(score, qs.s_xf, cfg.score_dir)) {
            qs.has_xf = true;
            qs.s_xf = score;
        }
    }
}

void scan_bin_range(
    const std::vector<homval::HitRecord>& hits,
    const RecordChunk& chunk,
    const Taxonomy& tax,
    TopcatTruth truth,
    const std::unordered_set<uint32_t>& query_ids,
    ScoreDir dir,
    ScanLocal& local) {
    for (size_t i = chunk.begin; i < chunk.end; ++i) {
        const homval::HitRecord& h = hits[i];
        const uint32_t qid = h.q;
        const uint32_t tid = h.t;
        if (!query_ids.count(qid))
            continue;

        const double score = h.score;
        const Track track = topcat_track(qid, tid, truth, tax);
        if (track == Track::Ignore) {
            ++local.ignored;
            continue;
        }

        ++local.considered;
        QueryScore& qs = local.scores[qid];
        if (track == Track::Tp) {
            if (!qs.has_tp || better_score(score, qs.s_tp, dir)) {
                qs.has_tp = true;
                qs.s_tp = score;
            }
        } else if (!qs.has_xf || better_score(score, qs.s_xf, dir)) {
            qs.has_xf = true;
            qs.s_xf = score;
        }
    }
}

unsigned effective_threads(unsigned requested, size_t work_size, size_t min_unit) {
    const unsigned max_by_work =
        static_cast<unsigned>(std::max<size_t>(1, work_size / min_unit));
    return std::max(1u, std::min(requested, max_by_work));
}

QueryState condense_text(
    const Config& cfg,
    const Taxonomy& tax,
    TopcatTruth truth,
    QueryState state) {
    std::ifstream size_probe(cfg.hits_path, std::ios::binary | std::ios::ate);
    if (!size_probe)
        die("cannot open hits " + cfg.hits_path);
    const size_t hits_size = static_cast<size_t>(size_probe.tellg());
    const unsigned nworkers = effective_threads(cfg.threads, hits_size, 2 * 1024 * 1024);
    const std::vector<ChunkJob> chunks = make_chunks(cfg.hits_path, nworkers);

    std::unordered_set<uint32_t> query_ids(state.query_order.begin(), state.query_order.end());
    std::vector<ScanLocal> locals(nworkers);
    std::vector<std::thread> workers;
    workers.reserve(nworkers);
    for (unsigned t = 0; t < nworkers; ++t) {
        workers.emplace_back([&, t]() {
            scan_text_chunk(cfg, tax, truth, query_ids, chunks[t], locals[t]);
        });
    }
    for (auto& w : workers)
        w.join();

    uint64_t ignored = 0;
    uint64_t considered = 0;
    for (ScanLocal& loc : locals) {
        merge_scores(state.scores, loc.scores, cfg.score_dir);
        ignored += loc.ignored;
        considered += loc.considered;
    }
    std::cerr << "ignored " << ignored << " hits, considered " << considered << '\n';
    return state;
}

QueryState condense_bin(
    const Config& cfg,
    const Taxonomy& tax,
    TopcatTruth truth,
    QueryState state,
    const std::vector<homval::HitRecord>& hits) {
    const unsigned nworkers = effective_threads(cfg.threads, hits.size(), 50000);
    const std::vector<RecordChunk> chunks = make_record_chunks(hits.size(), nworkers);

    std::unordered_set<uint32_t> query_ids(state.query_order.begin(), state.query_order.end());
    std::vector<ScanLocal> locals(nworkers);
    std::vector<std::thread> workers;
    workers.reserve(nworkers);
    for (unsigned t = 0; t < nworkers; ++t) {
        workers.emplace_back([&, t]() {
            scan_bin_range(hits, chunks[t], tax, truth, query_ids, cfg.score_dir, locals[t]);
        });
    }
    for (auto& w : workers)
        w.join();

    uint64_t ignored = 0;
    uint64_t considered = 0;
    for (ScanLocal& loc : locals) {
        merge_scores(state.scores, loc.scores, cfg.score_dir);
        ignored += loc.ignored;
        considered += loc.considered;
    }
    std::cerr << "ignored " << ignored << " hits, considered " << considered << '\n';
    return state;
}

void format_g(FILE* f, double v) {
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%g", v);
    std::fputs(buf, f);
}

void write_tcat(
    const Config& cfg,
    const std::string& truth_str,
    const Taxonomy& tax,
    const QueryState& state) {
    const Top3Stats top3 = compute_top3(state, cfg.score_dir);

    std::cerr << "writing " << cfg.output_path << "...\n";

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
    std::fprintf(out, "# N_possible_tp=%d\n", state.n_possible_tp);
    std::fprintf(out, "# TEPQ0.001=%.3f TEPQ0.01=%.3f TEPQ0.1=%.3f Top3=%.3f\n",
                 top3.tepq0001, top3.tepq001, top3.tepq01, top3.top3);
    std::fprintf(stderr, "# truth=%s TEPQ0.001=%.3f TEPQ0.01=%.3f TEPQ0.1=%.3f Top3=%.3f\n",
                 truth_str.c_str(), top3.tepq0001, top3.tepq001, top3.tepq01, top3.top3);

    std::fputs("domain\tS_tp\tS_xf\n", out);
    for (uint32_t qid : state.query_order) {
        const QueryScore& qs = state.scores.at(qid);
        std::fputs(tax.id_to_dom[qid].c_str(), out);
        std::fputc('\t', out);
        if (qs.has_tp)
            format_g(out, qs.s_tp);
        else
            std::fputc('.', out);
        std::fputc('\t', out);
        if (qs.has_xf)
            format_g(out, qs.s_xf);
        else
            std::fputc('.', out);
        std::fputc('\n', out);
    }
    std::fclose(out);
}

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
    for (const char* ext : {".tsv", ".txt", ".hits"}) {
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
                "Usage: hits_to_topcat (--hits PATH | --bin PATH.bin) --lookup PATH --derived-info PATH\n"
                "       (--truth {topsf,topfold} | --all-truths) --output PATH\n"
                "       (--evalue | --score) [--fields 1,2,3] [--algo NAME] [--reference DB]\n"
                "       [--threads N]\n"
                "\n"
                "With --all-truths, --output is a prefix; writes PREFIX.topsf, PREFIX.topfold.\n";
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
        cfg.truth = parse_topcat_truth(cfg.truth_str);
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

int run_text(const Config& cfg, const Taxonomy& tax) {
    QueryState state = load_query_state(cfg.derived_info_path, cfg.truth_str, tax);
    std::cerr << "reading hits " << cfg.hits_path << " truth=" << cfg.truth_str
              << " (|Q|=" << state.n_possible_tp << ")...\n";
    state = condense_text(cfg, tax, cfg.truth, state);
    write_tcat(cfg, cfg.truth_str, tax, state);
    return 0;
}

int run_bin(const Config& cfg, const Taxonomy& tax) {
    std::cerr << "loading bin " << cfg.bin_path << "...\n";
    const std::vector<homval::HitRecord> hits = homval::load_hbin(cfg.bin_path, tax.ndom_lookup);
    std::cerr << "bin records: " << hits.size() << '\n';

    const std::vector<TopcatTruth> truths =
        cfg.all_truths ? all_topcat_truths() : std::vector<TopcatTruth>{cfg.truth};

    for (TopcatTruth truth : truths) {
        const std::string truth_str = topcat_truth_name(truth);
        QueryState state = load_query_state(cfg.derived_info_path, truth_str, tax);
        std::cerr << "truth=" << truth_str << " (|Q|=" << state.n_possible_tp << ")...\n";
        state = condense_bin(cfg, tax, truth, state, hits);

        Config out_cfg = cfg;
        if (cfg.all_truths)
            out_cfg.output_path = cfg.output_path + "." + truth_str;
        write_tcat(out_cfg, truth_str, tax, state);
    }
    return 0;
}

int run(int argc, char** argv) {
    const Config cfg = parse_args(argc, argv);
    std::cerr << "\n\nloading lookup " << cfg.lookup_path << "...\n";
    const Taxonomy tax = load_lookup(cfg.lookup_path);
    if (cfg.use_bin)
        return run_bin(cfg, tax);
    return run_text(cfg, tax);
}

}  // namespace

int main(int argc, char** argv) {
    return run(argc, argv);
}
