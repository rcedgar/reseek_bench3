#pragma once

#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

namespace homval {

inline constexpr char HBIN_MAGIC[4] = {'H', 'B', 'I', 'N'};
inline constexpr uint32_t HBIN_VERSION = 1;

#pragma pack(push, 1)
struct BinHeader {
    char magic[4];
    uint32_t version;
    uint32_t ndom;
    uint64_t n_hits;
};

struct HitRecord {
    uint16_t q;
    uint16_t t;
    float score;
};
#pragma pack(pop)

static_assert(sizeof(BinHeader) == 20, "BinHeader size");
static_assert(sizeof(HitRecord) == 8, "HitRecord size");

inline void die_bin(const std::string& msg) {
    std::fprintf(stderr, "error: %s\n", msg.c_str());
    std::exit(1);
}

inline std::vector<HitRecord> load_hbin(const std::string& path, uint32_t expected_ndom = 0) {
    FILE* f = std::fopen(path.c_str(), "rb");
    if (!f)
        die_bin("cannot open bin " + path);

    BinHeader hdr{};
    if (std::fread(&hdr, sizeof(hdr), 1, f) != 1)
        die_bin("truncated bin header: " + path);
    if (std::memcmp(hdr.magic, HBIN_MAGIC, 4) != 0)
        die_bin("bad bin magic: " + path);
    if (hdr.version != HBIN_VERSION)
        die_bin("unsupported bin version in " + path);
    if (expected_ndom != 0 && hdr.ndom != expected_ndom)
        die_bin("bin ndom mismatch for " + path);
    if (hdr.ndom > 65535)
        die_bin("bin ndom too large for uint16 ids: " + path);

    std::vector<HitRecord> hits(static_cast<size_t>(hdr.n_hits));
    if (hdr.n_hits > 0) {
        if (std::fread(hits.data(), sizeof(HitRecord), hits.size(), f) != hits.size())
            die_bin("truncated bin records: " + path);
    }
    std::fclose(f);
    return hits;
}

inline void write_hbin(const std::string& path, uint32_t ndom, const std::vector<HitRecord>& hits) {
    FILE* f = std::fopen(path.c_str(), "wb");
    if (!f)
        die_bin("cannot write bin " + path);

    BinHeader hdr{};
    std::memcpy(hdr.magic, HBIN_MAGIC, 4);
    hdr.version = HBIN_VERSION;
    hdr.ndom = ndom;
    hdr.n_hits = hits.size();
    if (std::fwrite(&hdr, sizeof(hdr), 1, f) != 1)
        die_bin("write failed: " + path);
    if (!hits.empty()) {
        if (std::fwrite(hits.data(), sizeof(HitRecord), hits.size(), f) != hits.size())
            die_bin("write failed: " + path);
    }
    std::fclose(f);
}

}  // namespace homval
