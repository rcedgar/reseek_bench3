// build_derived_info.cpp -- Stage 0: precompute denominators and domain annotations.
//
// Reads a lookup TSV and writes derived_info JSON used by hits_to_edf and hits_to_topcat.

#include <fstream>
#include <iostream>
#include <string>

#include "derived_info.hpp"

namespace {

struct Args {
    std::string lookup_path;
    std::string output_path;
};

[[noreturn]] void die(const std::string& msg) {
    std::cerr << "error: " << msg << '\n';
    std::exit(1);
}

Args parse_args(int argc, char** argv) {
    Args args;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        auto need = [&](const char* name) -> std::string {
            if (i + 1 >= argc)
                die(std::string("missing value for ") + name);
            return argv[++i];
        };
        if (arg == "--lookup")
            args.lookup_path = need("--lookup");
        else if (arg == "--output")
            args.output_path = need("--output");
        else if (arg == "--help" || arg == "-h") {
            std::cout <<
                "Usage: build_derived_info --lookup PATH --output PATH\n"
                "\n"
                "Precompute homval denominators and domain annotations from a lookup TSV.\n";
            std::exit(0);
        } else
            die("unknown argument: " + arg);
    }
    if (args.lookup_path.empty() || args.output_path.empty())
        die("required: --lookup --output");
    return args;
}

}  // namespace

int main(int argc, char** argv) {
    const Args args = parse_args(argc, argv);

    std::cerr << "reading " << args.lookup_path << "...\n";
    const homval::LookupStrings data = homval::load_lookup_strings(args.lookup_path);
    const std::string reference = homval::lookup_basename(args.lookup_path);
    const std::string json = homval::build_derived_info_json(data, reference);

    std::ofstream out(args.output_path);
    if (!out)
        homval::die_tax("cannot write " + args.output_path);
    out << json;
    if (!out)
        homval::die_tax("write failed: " + args.output_path);

    std::cerr << "wrote " << args.output_path << '\n';
    std::cerr << "  ndom=" << data.domain_ids.size() << '\n';
    return 0;
}
