#include "graph.hpp"
#include <mutex>

std::mutex io_lock;

size_t printer_body::operator()(printer_input input) {
    if (input.variant != "") {
        size_t i;
        size_t fn = input.forward_mutation_nodes.size();
        size_t bn = input.backward_mutation_nodes.size();
        size_t fl = input.flagged_leaves.size();
        io_lock.lock();
        printf("%s\talt_alleles=%zu\tparsimony_score=%zu\tforward_mutation_nodes=", input.variant.c_str(), input.alt_alleles, \
                fn+bn); 
        i=0;
        for (auto n: input.forward_mutation_nodes) {
            if (++i != fn)
                printf("%s,", n.c_str());
            else
                printf("%s", n.c_str());
        }
        i=0;
        printf("\tbackward_mutation_nodes=");
        for (auto n: input.backward_mutation_nodes) {
            if (++i != bn)
                printf("%s,", n.c_str());
            else
                printf("%s", n.c_str());
        }
        i=0;
        printf("\tforward_mutation_clade_sizes=");
        for (auto n: input.forward_mutation_clade_sizes) {
            if (++i != fn)
                printf("%zu,", n);
            else
                printf("%zu", n);
        }
        i=0;
        printf("\tbackward_mutation_clade_sizes=");
        for (auto n: input.backward_mutation_clade_sizes) {
            if (++i != bn)
                printf("%zu,", n);
            else
                printf("%zu", n);
        }
        i=0;
        printf("\tflagged_leaves=");
        for (auto n: input.flagged_leaves) {
            if (++i != fl)
                printf("%s,", n.c_str());
            else
                printf("%s", n.c_str());
        }
        printf("\n");
        io_lock.unlock();
    }

    return 1;
}
