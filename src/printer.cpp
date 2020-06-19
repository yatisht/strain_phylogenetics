#include "graph.hpp"

std::mutex io_lock;

std::unordered_map<int8_t, char> nuc_to_char = {{0,'A'}, {1, 'C'}, {2, 'G'}, {3, 'T'}};

size_t printer_body::operator()(printer_input input) {
    if (input.print_vcf) {
        if (input.variant == "") {
            io_lock.lock();
            printf("%s\n", input.line.c_str());
            io_lock.unlock();
        }
        else {
            size_t tot, curr;
            
            io_lock.lock();
            printf("%s\t%s\t%s\t", input.words[0].c_str(), input.words[1].c_str(), input.words[2].c_str());
            printf("%c\t", nuc_to_char[input.ref_nuc]);
            tot = input.alts.size();
            curr = 0;
            if (tot>0) {
                for (auto a: input.alts) {
                    if (++curr != tot)
                        printf("%c,", nuc_to_char[a]);
                    else
                        printf("%c\t", nuc_to_char[a]);
                }
            }
            else {
                printf(".\t");
            }
            printf("%s\t%s\t%s\t%s\t", input.words[5].c_str(), input.words[6].c_str(), input.words[7].c_str(), input.words[8].c_str());
            tot = input.sample_variant_id.size();
            curr = 0;
            if (tot>0) {
                for (auto v: input.sample_variant_id) {
                    if (++curr != tot)
                        printf("%d\t", v);
                    else
                        printf("%d",v);
                }
            }
            printf ("\n");
            io_lock.unlock();
        }
    }
    else {
        if (input.variant != "") {
            size_t i;
            size_t tot, curr;
            size_t parsimony_score=0;
            int8_t n1;
            int8_t n2;

            for (i=0; i<16; i++) {
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

            for (i=0; i<16; i++) {
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

            for (i=0; i<16; i++) {
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
    }

    return 1;
}
