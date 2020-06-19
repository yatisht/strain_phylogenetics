#include "graph.hpp"

std::mutex error_lock;

printer_input mapper_body::operator()(mapper_input input) {
    printer_input output;

    output.variant = input.variant;
    output.line = input.line;
    output.print_vcf = input.print_vcf;
    output.ref_nuc = input.ref_nuc;
    if (input.variant != "") {
        size_t num_nodes = input.bfs->size();
        std::vector<int8_t> states(num_nodes);
        std::vector<std::vector<int>> scores(num_nodes);

        for (size_t i=0; i<num_nodes; i++) {
            scores[i].resize(4);
            for (int j=0; j<4; j++) {
                scores[i][j] = 0;
            }
            states[i] = 0;
        }

        for (auto l: input.T->get_leaves()) {
            size_t node_idx = (*input.bfs_idx)[l->identifier];
            for (int j=0; j<4; j++) {
                if (j != input.ref_nuc) 
                    scores[node_idx][j] = (int) num_nodes;
            }
        }

        for (int j=0; j<4; j++) {
            output.alt_alleles[j] = 0;
        }

        for (auto v: input.variants) {
            size_t pos = std::get<0> (v);
            std::vector<int8_t> nucs = std::get<1> (v);
            std::string nid = (*input.variant_ids)[pos];
            if ((*input.bfs_idx).find(nid) != (*input.bfs_idx).end()) {
                size_t idx= (*input.bfs_idx)[nid];
                for (int j=0; j<4; j++) {
                    scores[idx][j] = (int) num_nodes;
                }
                for (auto nuc: nucs) {
                    if (nuc < 4) {
                        scores[idx][nuc] = 0;
                    }
                }
            }
            else {
                error_lock.lock();
                fprintf (stderr, "ERROR: %s in VCF not found in the tree!\n", nid.c_str());
                exit(1);
                error_lock.unlock();
            }
        }

        // Sankoff: forward pass
        for (auto it=(*input.bfs).rbegin(); it!=(*input.bfs).rend(); it++) {
            auto node = (*it);
            auto node_idx = (*input.bfs_idx)[node->identifier];

            if (!node->is_leaf()) {
                for (auto child: (*node).children) {
                    auto c_id = child->identifier;
                    auto c_idx = (*input.bfs_idx)[c_id];
                    for (int j=0; j<4; j++) {
                        int min_s = (int) num_nodes+1;
                        for (int k=0; k<4; k++) {
                            int c_s;
                            if (k==j) {
                                c_s = scores[c_idx][k];
                            }
                            else {
                                c_s = scores[c_idx][k]+1;
                            }
                            if (c_s < min_s) {
                                min_s = c_s;
                            }
                        }
                        scores[node_idx][j] += min_s;
                    }
                }
            }
        }

        // Sankoff: backward pass
        for (auto it=(*input.bfs).begin(); it!=(*input.bfs).end(); it++) {
            auto node = (*it);
            auto node_idx = (*input.bfs_idx)[node->identifier];

            int8_t par_state = 0;
            if (node->parent != NULL) {
                auto par = node->parent;
                auto par_id = par->identifier;
                auto par_idx = (*input.bfs_idx)[par_id];
                par_state = states[par_idx];
            }
            int8_t state = par_state;
            int min_s = scores[node_idx][par_state];
            for (int j=0; j<4; j++) {
                if (scores[node_idx][j] < min_s) {
                    min_s = scores[node_idx][j];
                    state = j;
                }
            }
            states[node_idx] = state;

            if (node->parent == NULL) {
                par_state = state;
            }
            auto root_state = states[(*input.bfs_idx)[input.T->root->identifier]];
            if (state != par_state) {
                auto clade_state = root_state;
                auto par_id = node->parent->identifier;
                for (auto anc: input.T->rsearch(par_id)) {
                    auto anc_state = states[(*input.bfs_idx)[anc->identifier]];
                    if (anc_state != par_state) {
                        clade_state = anc_state;
                        break;
                    }
                }
                auto mut_idx = 4*par_state + state;
                output.mutation_nodes[mut_idx].push_back(node->identifier);
                if (state != clade_state) {
                    output.mutation_type[mut_idx].push_back('F');
                }
                else {
                    output.mutation_type[mut_idx].push_back('B');
                }
            }
        }

        for (auto mut_idx=0; mut_idx<16; mut_idx++) {
            for (auto nid: output.mutation_nodes[mut_idx]) {
                auto leaves = input.T->get_leaves(nid);
                size_t clade_size = 0;
                size_t num_leaves = leaves.size();
                for (auto l: leaves) {
                    if (num_leaves < 4) {
                        output.flagged_leaves.insert(l->identifier);
                    }
                    if (states[(*input.bfs_idx)[l->identifier]] == states[(*input.bfs_idx)[nid]]) {
                        clade_size++;
                    }
                }
                output.mutation_clade_sizes[mut_idx].push_back(clade_size);
            }
        }

        auto all_leaves = input.T->get_leaves();
        for (auto l: all_leaves) {
            auto allele = states[(*input.bfs_idx)[l->identifier]];
            if (allele != input.ref_nuc) {
                output.alt_alleles[allele]++;
            }
        }

        if (input.print_vcf) {
            for (auto w: input.words) {
                output.words.push_back(w);
            }
            std::unordered_map<int8_t, int8_t> state_to_variant_id;
            int8_t curr_variant_id = 1;
            for (int8_t j=0; j<4; j++) {
                if (output.alt_alleles[j] > 0) {
                    output.alts.push_back(j);
                    state_to_variant_id[j] = curr_variant_id++;
                }
            }
            for (size_t i=0; i<(*input.variant_ids).size(); i++) {
                auto nid = (*input.variant_ids)[i];
                if (input.bfs_idx->find(nid) == input.bfs_idx->end()) {
                    error_lock.lock();
                    fprintf (stderr, "ERROR: %s in VCF not found in the tree!\n", nid.c_str());
                    exit(1);
                    error_lock.unlock();
                }
                int8_t allele = states[(*input.bfs_idx)[nid]];
                if (allele == input.ref_nuc) {
                    output.sample_variant_id.push_back(0);
                }
                else {
                    output.sample_variant_id.push_back(state_to_variant_id[allele]);
                }
            }

        }

    }

    return output;
}
