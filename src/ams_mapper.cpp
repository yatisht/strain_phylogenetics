#include "ams_graph.hpp"

std::mutex data_lock;

int mapper_body::operator()(mapper_input input) {
    if (input.variant_pos >= 0) {
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
                size_t missing_idx = 0;
                data_lock.lock();
                for (missing_idx=0; missing_idx<input.missing_samples->size(); missing_idx++) {
                    if ((*input.missing_samples)[missing_idx] == nid) {
                        break;
                    }
                }
                if (input.missing_samples->size() == missing_idx) {
                    (*input.missing_samples).push_back(nid);
                    mutation m;
                    m.position = input.variant_pos;
                    m.ref_nuc = input.ref_nuc;
                    for (auto n: nucs) {
                        m.mut_nuc.push_back(n);
                    }
                    std::vector<mutation> sample_mutations;
                    sample_mutations.push_back(m);
                    (*input.missing_sample_mutations).push_back(sample_mutations);
                }
                else {
                    mutation m;
                    m.position = input.variant_pos;
                    m.ref_nuc = input.ref_nuc;
                    for (auto n: nucs) {
                        m.mut_nuc.push_back(n);
                    }
                    (*input.missing_sample_mutations)[missing_idx].push_back(m);
                }
                data_lock.unlock();
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
            if (scores[node_idx][input.ref_nuc] == min_s) {
                state = input.ref_nuc;
            }
            states[node_idx] = state;

            if (node->parent == NULL) {
                par_state = state;
            }
            if (state != par_state) {
                if (input.node_mutations->find(node) != input.node_mutations->end()) {
                    mutation m;
                    m.position = input.variant_pos;
                    m.ref_nuc = input.ref_nuc;
                    m.mut_nuc.push_back(state);
                    (*input.node_mutations)[node].push_back(m);
                }
                else {
                    mutation m;
                    m.position = input.variant_pos;
                    m.ref_nuc = input.ref_nuc;
                    m.mut_nuc.push_back(state);
                    std::vector<mutation> node_mut;
                    node_mut.push_back(m);
                    (*input.node_mutations)[node] = node_mut;
                }
            }
        }

    }

    return 1;
}

int mapper2_body::operator()(mapper2_input input) {

    *input.set_difference = 0;

    std::vector<mutation> ancestral_mutations;
    for (auto n: input.T->rsearch(input.node->identifier)) {
        if (input.node_mutations->find(n) != input.node_mutations->end()) {
            for (auto m: (*input.node_mutations)[n]) {
                ancestral_mutations.push_back(m);
            }
        }
    }
    if (input.node_mutations->find(input.node) != input.node_mutations->end()) {
        for (auto m: (*input.node_mutations)[input.node]) {
            ancestral_mutations.push_back(m);
        }
    }

    for (auto m1: (*input.missing_sample_mutations)) {
        bool found = false;
        bool has_ref = false;
        for (auto m2: ancestral_mutations) {
            if (m1.position == m2.position) {
                assert(m2.mut_nuc.size() == 1);
                auto anc_nuc = m2.mut_nuc[0];
                for (auto nuc: m1.mut_nuc) {
                    if (nuc == anc_nuc) {
                        found = true;
                    }
                }
            }
        }
        if (!found) {
            *input.set_difference += 1;
            for (auto nuc: m1.mut_nuc) {
                if (nuc == m1.ref_nuc) {
                    has_ref = true;
                }
            }
            mutation m;
            m.position = m1.position;
            m.ref_nuc = m1.ref_nuc;
            if (has_ref) {
                m.mut_nuc.push_back(m1.ref_nuc);
            }
            else {
                m.mut_nuc.push_back(m1.mut_nuc[0]);
            }
            input.excess_mutations->push_back(m);
        }
    }

    for (auto m1: ancestral_mutations) {
        bool found = false;
        for (auto m2: (*input.missing_sample_mutations)) {
            if (m1.position == m2.position) {
                auto anc_nuc = m1.mut_nuc[0];
                for (auto nuc: m2.mut_nuc) {
                    if (nuc == anc_nuc) {
                        found = true;
                    }
                }
            }
        }
        if (!found) {
            *input.set_difference += 1;
            mutation m;
            m.position = m1.position;
            m.ref_nuc = m1.ref_nuc;
            m.mut_nuc.push_back(m1.ref_nuc);
            (*input.excess_mutations).push_back(m);
        }
    }
    
    return 1;
}

