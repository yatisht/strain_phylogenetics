#include "graph.hpp"
#include <mutex>

std::mutex io_lock;

std::unordered_map<int8_t, char> nuc_to_char = {{0,'A'}, {1, 'C'}, {2, 'G'}, {3, 'T'}};

size_t printer_body::operator()(printer_input input) {
    if (input.variant != "") {
        size_t i;
        size_t tot, curr;
        size_t parsimony_score=0;
        int8_t n1;
        int8_t n2;

        for (i=0; i<12; i++) {
            parsimony_score += input.mutation_nodes[i].size();
        }

        io_lock.lock();
        printf("%s\t", input.variant.c_str());

        for (i=0; i<4; i++) {
            if (input.alt_alleles[i] > 0) {
                printf("%c_alt_alleles=%zu\t", nuc_to_char[i], input.alt_alleles[i]);
            }
        }
        
        printf("parsimony_score=%zu\t", parsimony_score);
        
        for (i=0; i<12; i++) {
            curr = 0;
            tot = input.mutation_nodes[i].size();
            if (tot>0) {
                n1 = i/4;
                n2 = i%4;
                printf("\t%c>%c_mutation_nodes=", nuc_to_char[n1], nuc_to_char[n2]);
                for (auto n: input.mutation_nodes[i]) {
                    if (++curr != tot)
                        printf("%s[%c],", n.c_str(), input.mutation_type[i][curr-1]);
                    else
                        printf("%s[%c]", n.c_str(), input.mutation_type[i][curr-1]);
                }
            }
        }
        
        for (i=0; i<12; i++) {
            curr = 0;
            tot = input.mutation_clade_sizes[i].size();
            if (tot>0) {
                n1 = i/4;
                n2 = i%4;
                printf("\t%c>%c_mutation_clade_sizes=", nuc_to_char[n1], nuc_to_char[n2]);
                for (auto n: input.mutation_clade_sizes[i]) {
                    if (++curr != tot)
                        printf("%zu,", n);
                    else
                        printf("%zu", n);
                }
            }
        }

        curr = 0;
        tot = input.flagged_leaves.size();
        if (tot>0) {
            printf("\tflagged_leaves=");
            for (auto n: input.flagged_leaves) {
                if (++curr != tot)
                    printf("%s,", n.c_str());
                else
                    printf("%s", n.c_str());
            }
        }

        printf("\n");
        io_lock.unlock();
    }

    return 1;
}
