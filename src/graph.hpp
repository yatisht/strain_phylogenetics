#include "tree.hpp"
#include <set>

struct mapper_input {
    Tree* T;
    int8_t ref_nuc;
    std::string variant;
    std::vector<Node*>* bfs;
    std::unordered_map<std::string, size_t>* bfs_idx;
    std::vector<std::tuple<size_t, std::vector<int8_t>>> variants;
    std::vector<std::string>* variant_ids;
};

struct printer_input {
    size_t alt_alleles;
    std::string variant;
    std::vector<std::string> forward_mutation_nodes;
    std::vector<std::string> backward_mutation_nodes;
    std::vector<size_t> forward_mutation_clade_sizes;
    std::vector<size_t> backward_mutation_clade_sizes;
//    std::set<std::string> leaves_affected_forward;
//    std::set<std::string> leaves_affected_backward;
    std::set<std::string> flagged_leaves;
};

struct mapper_body {
    printer_input operator()(mapper_input input);
};

struct printer_body {
    size_t operator()(printer_input input);
};

