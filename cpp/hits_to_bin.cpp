// hits_to_bin.cpp -- convert text hits TSV to compact binary for fast rescans.
//
// Record layout: uint16 q, uint16 t, float score (8 bytes each).
// Only rows where both domains appear in the lookup are stored.

#include <charconv>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <string_view>
#include <thread>
#include <unordered_map>
#include <vector>

#include "hits_bin.hpp"

namespace {

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

void parse_fields(const std::string& spec, int& q, int& t, int& s) {
    int q1 = 0, t1 = 0, s1 = 0;
    if (std::sscanf(spec.c_str(), "%d,%d,%d", &q1, &t1, &s1) != 3 || q1 < 1 || t1 < 1 || s1 < 1)
        die("fields must be three 1-based indices, e.g. 1,2,3");
    q = q1 - 1;
    t = t1 - 1;
    s = s1 - 1;
}

std::unordered_map<std::string, uint16_t> load_dom_index(const std::string& path) {
    std::unordered_map<std::string, uint16_t> dom_index;
    std::ifstream in(path);
    if (!in)
        die("cannot open lookup " + path);

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
        if (dom_index.count(dom))
            die(path + ":" + std::to_string(line_no) + ": duplicate domain " + dom);
        if (dom_index.size() >= 65535)
            die("lookup has too many domains for uint16 ids");
        dom_index.emplace(dom, static_cast<uint16_t>(dom_index.size()));
    }
    return dom_index;
}

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
    if (nthreads == 0)
        nthreads = 1;

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

struct LocalBin {
    std::vector<homval::HitRecord> hits;
    uint64_t skipped = 0;
};

void scan_chunk(
    const std::string& hits_path,
    const ChunkJob& chunk,
    const std::unordered_map<std::string, uint16_t>& dom_index,
    int q_idx,
    int t_idx,
    int s_idx,
    LocalBin& local) {
    std::ifstream in(hits_path, std::ios::binary);
    if (!in)
        die("cannot open hits " + hits_path);
    in.seekg(static_cast<std::streamoff>(chunk.start));

    std::string line;
    std::string fq, ft, fs;
    while (read_chunk_line(in, chunk, line)) {
        trim_inplace(line);
        if (line.empty() || line[0] == '#')
            continue;
        if (!field_at(line, q_idx, fq) || !field_at(line, t_idx, ft) || !field_at(line, s_idx, fs))
            continue;

        std::string qdom, tdom;
        if (!parse_domain_label(fq.data(), fq.data() + fq.size(), qdom) ||
            !parse_domain_label(ft.data(), ft.data() + ft.size(), tdom))
            continue;

        auto iq = dom_index.find(qdom);
        auto it = dom_index.find(tdom);
        if (iq == dom_index.end() || it == dom_index.end()) {
            ++local.skipped;
            continue;
        }

        double score = 0.0;
        if (!parse_double_field(fs, score))
            continue;

        local.hits.push_back({iq->second, it->second, static_cast<float>(score)});
    }
}

unsigned effective_threads(unsigned requested, size_t file_size) {
    const unsigned max_by_size =
        static_cast<unsigned>(std::max<size_t>(1, file_size / (2 * 1024 * 1024)));
    return std::max(1u, std::min(requested, max_by_size));
}

struct Args {
    std::string hits_path;
    std::string lookup_path;
    std::string output_path;
    std::string fields_spec = "1,2,3";
    int q_idx = 0;
    int t_idx = 1;
    int s_idx = 2;
    unsigned threads = 0;
};

Args parse_args(int argc, char** argv) {
    Args args;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        auto need = [&](const char* name) {
            if (i + 1 >= argc)
                die(std::string("missing value for ") + name);
            return std::string(argv[++i]);
        };
        if (arg == "--hits")
            args.hits_path = need("--hits");
        else if (arg == "--lookup")
            args.lookup_path = need("--lookup");
        else if (arg == "--output")
            args.output_path = need("--output");
        else if (arg == "--fields")
            args.fields_spec = need("--fields");
        else if (arg == "--threads") {
            const long n = std::stol(need("--threads"));
            if (n < 1)
                die("--threads must be >= 1");
            args.threads = static_cast<unsigned>(n);
        }
        else if (arg == "--help" || arg == "-h") {
            std::cout <<
                "Usage: hits_to_bin --hits PATH --lookup PATH --output PATH.bin\n"
                "       [--fields 1,2,3] [--threads N]\n";
            std::exit(0);
        }
        else
            die("unknown argument: " + arg);
    }
    if (args.hits_path.empty() || args.lookup_path.empty() || args.output_path.empty())
        die("required: --hits --lookup --output");
    parse_fields(args.fields_spec, args.q_idx, args.t_idx, args.s_idx);
    if (args.threads == 0)
        args.threads = std::max(1u, std::thread::hardware_concurrency());
    return args;
}

}  // namespace

int main(int argc, char** argv) {
    const Args args = parse_args(argc, argv);
    const auto dom_index = load_dom_index(args.lookup_path);
    const uint32_t ndom = static_cast<uint32_t>(dom_index.size());

    std::ifstream size_probe(args.hits_path, std::ios::binary | std::ios::ate);
    if (!size_probe)
        die("cannot open hits " + args.hits_path);
    const size_t hits_size = static_cast<size_t>(size_probe.tellg());
    const unsigned nworkers = effective_threads(args.threads, hits_size);
    const auto chunks = make_chunks(args.hits_path, nworkers);

    std::cerr << "hits_to_bin: " << args.hits_path << " -> " << args.output_path
              << " (" << nworkers << " threads)\n";

    std::vector<LocalBin> locals(nworkers);
    std::vector<std::thread> workers;
    workers.reserve(nworkers);
    for (unsigned t = 0; t < nworkers; ++t) {
        workers.emplace_back([&, t]() {
            scan_chunk(args.hits_path, chunks[t], dom_index,
                       args.q_idx, args.t_idx, args.s_idx, locals[t]);
        });
    }
    for (auto& w : workers)
        w.join();

    size_t total = 0;
    uint64_t skipped = 0;
    for (const LocalBin& loc : locals) {
        total += loc.hits.size();
        skipped += loc.skipped;
    }

    std::vector<homval::HitRecord> all;
    all.reserve(total);
    for (LocalBin& loc : locals) {
        all.insert(all.end(), loc.hits.begin(), loc.hits.end());
        loc.hits.clear();
    }

    homval::write_hbin(args.output_path, ndom, all);
    std::cerr << "wrote " << all.size() << " records (skipped " << skipped
              << " unknown-domain rows)\n";
    return 0;
}
