#include "tree.hpp"
#include <set>
#include <mutex>

struct mapper_input {
    Tree* T;
    int8_t ref_nuc;
    std::string variant;
    bool print_vcf;
    std::string line;
    std::vector<std::string> words;
    std::vector<Node*>* bfs;
    std::unordered_map<std::string, size_t>* bfs_idx;
    std::vector<std::tuple<size_t, std::vector<int8_t>>> variants;
    std::vector<std::string>* variant_ids;
};

struct printer_input {
    bool print_vcf;
    std::string line;
    std::vector<std::string> words;
    int8_t ref_nuc;
    std::vector<int8_t> alts;
    std::vector<int8_t> sample_variant_id;
    std::string variant;
    size_t alt_alleles[4];
    std::vector<std::string> mutation_nodes[16];
    std::vector<size_t> mutation_clade_sizes[16];
    std::vector<char> mutation_type[16];
    std::set<std::string> flagged_leaves;
};

struct mapper_body {
    printer_input operator()(mapper_input input);
};

struct printer_body {
    size_t operator()(printer_input input);
};

