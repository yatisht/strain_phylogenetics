#include <fstream>
#include <boost/program_options.hpp> 
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <iostream>
#include <tbb/flow_graph.h>
#include <tbb/reader_writer_lock.h>
#include <tbb/scalable_allocator.h>
#include <tbb/task_scheduler_init.h>
#include "ams_graph.hpp"
#include <omp.h>
#include "parsimony.pb.h"

namespace po = boost::program_options;

std::vector<int8_t> get_nuc_id (char c) {
    switch (c) {
        case 'a':
        case 'A': return std::vector<int8_t>{0};
                  break;
        case 'c':
        case 'C': return std::vector<int8_t>{1};
                  break;
        case 'g':
        case 'G': return std::vector<int8_t>{2};
                  break;
        case 't':
        case 'T': return std::vector<int8_t>{3};
                  break;
        case 'R': return std::vector<int8_t>{0,2};
                  break;
        case 'Y': return std::vector<int8_t>{1,3};
                  break;
        case 'S': return std::vector<int8_t>{1,2};
                  break;
        case 'W': return std::vector<int8_t>{0,3};
                  break;
        case 'K': return std::vector<int8_t>{2,3};
                  break;
        case 'M': return std::vector<int8_t>{0,1};
                  break;
        case 'B': return std::vector<int8_t>{1,2,3};
                  break;
        case 'D': return std::vector<int8_t>{0,2,3};
                  break;
        case 'H': return std::vector<int8_t>{0,1,3};
                  break;
        case 'V': return std::vector<int8_t>{0,1,2};
                  break;
        default : return std::vector<int8_t>{0,1,2,3};
                  break;
    }
}

char get_nuc_char (int8_t nuc_id) {
    switch (nuc_id) {
        case 0: return 'A';
                break;
        case 1: return 'C';
                break;
        case 2: return 'G';
                break;
        case 3: return 'T';
                break;
        default : return 'N'; 
                  break;
    }
}

int main(int argc, char** argv){

    std::string tree_filename;
    std::string din_filename;
    std::string dout_filename;
    std::string vcf_filename;
    uint32_t num_cores = tbb::task_scheduler_init::default_num_threads();
    uint32_t num_threads;
    bool collapse_tree=false;
    bool print_uncondensed_tree = false;
    size_t print_subtrees_size=0;
    po::options_description desc{"Options"};
    desc.add_options()
        ("vcf", po::value<std::string>(&vcf_filename)->required(), "Input VCF file (in uncompressed or gzip-compressed format)")
        ("tree", po::value<std::string>(&tree_filename)->default_value(""), "Input tree file")
        ("load-assignments", po::value<std::string>(&din_filename)->default_value(""), "Load existing tree and parsimonious assignments")
        ("save-assignments", po::value<std::string>(&dout_filename)->default_value(""), "Save output tree and parsimonious assignments")
        ("threads", po::value<uint32_t>(&num_threads)->default_value(num_cores), "Number of threads")
        ("collapse-final-tree", po::bool_switch(&collapse_tree), "Collapse internal nodes of the output tree with no mutations.")
        ("print_uncondensed-final-tree", po::bool_switch(&print_uncondensed_tree), "Print the final tree in uncondensed format.")
        ("print-subtrees-size", po::value<size_t>(&print_subtrees_size)->default_value(0), \
         "Print minimum set of subtrees of size equal of larger than this value covering the missing samples.")
        ("help", "Print help messages");
    
    po::options_description all_options;
    all_options.add(desc);

    //    po::positional_options_description p;
    //    p.add("tree", 1);
    //    p.add("vcf", 1);
    //    p.add("threads", 1);

    po::variables_map vm;
    try{
        po::store(po::command_line_parser(argc, argv).options(all_options).run(), vm);
        po::notify(vm);
    }
    catch(std::exception &e){
        std::cout << desc << std::endl;
        return 1;
    }

    if (print_subtrees_size == 1) {
        std::cerr << "ERROR: print-subtrees-size should be larger than 1\n";
        return 1;
    }

    if ((din_filename != "") && collapse_tree) {
        std::cerr << "ERROR: cannot load assignments and collapse tree simulaneously.\n";
        return 1;
    }

    Timer timer; 

    omp_set_num_threads(num_threads);
    omp_lock_t omplock;
    omp_init_lock(&omplock);

#if SAVE_PROFILE == 1
    Instrumentor::Get().BeginSession("test-main", "p1.json");
#endif

    Tree T;

    bool header_found = false;
    std::vector<std::string> variant_ids;
    std::vector<std::string> missing_samples;
    std::unordered_map<Node*, std::vector<mutation>> node_mutations;
    std::vector<std::vector<mutation>> missing_sample_mutations;
    size_t num_missing = 0;

    std::unordered_map<std::string, std::vector<std::string>> condensed_nodes;
    Tree condensed_T; 
    std::unordered_map<Node*, std::vector<mutation>> condensed_node_mutations;
    
    std::vector<Node*> bfs;
    std::unordered_map<std::string, size_t> bfs_idx;
    
    if (tree_filename != "") {
        T = create_tree_from_newick(tree_filename);
        bfs = T.breadth_first_expansion();

        for (size_t idx = 0; idx < bfs.size(); idx++) {
            bfs_idx[bfs[idx]->identifier] = idx;
        }

        std::ifstream infile(vcf_filename, std::ios_base::in | std::ios_base::binary);
        boost::iostreams::filtering_istream instream;
        try {
            if (vcf_filename.find(".gz\0") != std::string::npos) {
                instream.push(boost::iostreams::gzip_decompressor());
            }
            instream.push(infile);
        }
        catch(const boost::iostreams::gzip_error& e) {
            std::cout << e.what() << '\n';
        }

        fprintf(stderr, "Computing parsimonious assignments for input variants.\n"); 
        timer.Start();

        tbb::task_scheduler_init init(num_threads);

        tbb::flow::graph mapper_graph;

        tbb::flow::function_node<mapper_input, int> mapper(mapper_graph, tbb::flow::unlimited, mapper_body());
        tbb::flow::source_node <mapper_input> reader (mapper_graph,
                [&] (mapper_input &inp) -> bool {
                std::string s;
                std::getline(instream, s);
                std::vector<std::string> words;
                split(s, words);
                inp.variant_pos = -1;
                if ((not header_found) && (words.size() > 1)) {
                if (words[1] == "POS") {
                for (size_t j=9; j < words.size(); j++) {
                variant_ids.emplace_back(words[j]);
                if (bfs_idx.find(words[j]) == bfs_idx.end()) {
                missing_samples.emplace_back(words[j]);
                num_missing++;
                }
                }
                missing_sample_mutations.resize(num_missing);
                header_found = true;
                }
                }
                else if (header_found) {
                    if (words.size() != 9+variant_ids.size()) {
                        fprintf(stderr, "ERROR! Incorrect VCF format.\n");
                        exit(1);
                    }
                    std::vector<std::string> alleles;
                    alleles.clear();
                    inp.variant_pos = std::stoi(words[1]); 
                    split(words[4], ',', alleles);
                    inp.T = &T;
                    inp.bfs = &bfs;
                    inp.bfs_idx = &bfs_idx;
                    inp.variant_ids = &variant_ids;
                    inp.missing_samples = &missing_samples;
                    inp.node_mutations = &node_mutations;
                    inp.missing_sample_mutations = &missing_sample_mutations;
                    auto ref_nucs = get_nuc_id(words[3][0]);
                    assert(ref_nucs.size() == 1);
                    inp.ref_nuc = ref_nucs[0]; 
                    inp.variants.clear();
                    for (size_t j=9; j < words.size(); j++) {
                        if (isdigit(words[j][0])) {
                            int allele_id = std::stoi(words[j]);
                            if (allele_id > 0) { 
                                std::string allele = alleles[allele_id-1];
                                inp.variants.emplace_back(std::make_tuple(j-9, get_nuc_id(allele[0])));
                            }
                        }
                        else {
                            inp.variants.emplace_back(std::make_tuple(j-9, get_nuc_id('N')));
                        }
                    }
                }
        //check if reached end-of-file
        int curr_char = instream.peek();
        if(curr_char == EOF)
            return false;
        else
            return true;
                }, true );
        tbb::flow::make_edge(reader, mapper);
        mapper_graph.wait_for_all();

        fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
    }
    else if (din_filename != "") {
        fprintf(stderr, "Loading existing assignments\n");
        timer.Start();
        
        Parsimony::data data;

        std::ifstream inpfile(din_filename, std::ios::in | std::ios::binary);
        data.ParseFromIstream(&inpfile);
        inpfile.close();

        T = create_tree_from_newick_string(data.newick());

        auto dfs = T.depth_first_expansion();

        for (size_t idx = 0; idx < dfs.size(); idx++) {
            auto mutation_list = data.node_mutations(idx);
            auto node = dfs[idx];
            node_mutations.insert(std::pair<Node*, std::vector<mutation>>(node, std::vector<mutation>()));  
            for (int k = 0; k < mutation_list.mutation_size(); k++) {
                auto mut = mutation_list.mutation(k);
                mutation m;
                m.position = mut.position();
                m.ref_nuc = mut.ref_nuc();
                m.par_nuc = mut.par_nuc();
                for (int n = 0; n < mut.mut_nuc_size(); n++) {
                    m.mut_nuc.emplace_back(mut.mut_nuc(n));
                }
                node_mutations[node].emplace_back(m);
            }
        }
        
        size_t num_condensed_nodes = static_cast<size_t>(data.condensed_nodes_size());
        for (size_t idx = 0; idx < num_condensed_nodes; idx++) {
            auto cn = data.condensed_nodes(idx);
            condensed_nodes.insert(std::pair<std::string, std::vector<std::string>>(cn.node_name(), std::vector<std::string>(cn.condensed_leaves_size())));
            for (int k = 0; k < cn.condensed_leaves_size(); k++) {
                condensed_nodes[cn.node_name()][k] =cn.condensed_leaves(k);
            }
        }
        
        fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
        
        
        fprintf(stderr, "Loading VCF file\n");
        timer.Start();

        std::ifstream infile(vcf_filename, std::ios_base::in | std::ios_base::binary);
        boost::iostreams::filtering_istream instream;
        try {
            if (vcf_filename.find(".gz\0") != std::string::npos) {
                instream.push(boost::iostreams::gzip_decompressor());
            }
            instream.push(infile);
        }
        catch(const boost::iostreams::gzip_error& e) {
            std::cout << e.what() << '\n';
        }
        std::vector<size_t> missing_idx;
        std::string s;
        while (instream.peek() != EOF) {
            std::getline(instream, s);
            std::vector<std::string> words;
            split(s, words);
            if ((not header_found) && (words.size() > 1)) {
                if (words[1] == "POS") {
                    for (size_t j=9; j < words.size(); j++) {
                        variant_ids.emplace_back(words[j]);
                        if (T.get_node(words[j]) == NULL) {
                            missing_samples.emplace_back(words[j]);
                            num_missing++;
                            missing_idx.emplace_back(j);
                        }
                    }
                    missing_sample_mutations.resize(num_missing);
                    header_found = true;
                }
            }
            else if (header_found) {
                if (words.size() != 9+variant_ids.size()) {
                    fprintf(stderr, "ERROR! Incorrect VCF format. Expected %zu columns but got %zu.\n", 9+variant_ids.size(), words.size());
                    exit(1);
                }
                std::vector<std::string> alleles;
                alleles.clear();
                split(words[4], ',', alleles);
                for (auto j: missing_idx) {
                    auto iter = std::find(missing_samples.begin(), missing_samples.end(), variant_ids[j-9]);
                    if (iter != missing_samples.end()) {
                        auto mutations_iter = missing_sample_mutations.begin() + (iter - missing_samples.begin());
                        mutation m;
                        m.position = std::stoi(words[1]);
                        auto ref_nucs = get_nuc_id(words[3][0]);
                        assert(ref_nucs.size() == 1);
                        m.ref_nuc = ref_nucs[0];
                        m.par_nuc = ref_nucs[0];
                        if (isdigit(words[j][0])) {
                            int allele_id = std::stoi(words[j]);
                            if (allele_id > 0) { 
                                std::string allele = alleles[allele_id-1];
                                for (auto n: get_nuc_id(allele[0])) {
                                    m.mut_nuc.emplace_back(n);
                                }
                                (*mutations_iter).emplace_back(m);
                            }
                        }
                    }
                }
            }
        }
        fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
    }
    else {
        fprintf(stderr, "Error! No input tree or assignment file provided!\n");
        exit(1);
    }

    fprintf(stderr, "Found %zu missing samples.\n", missing_samples.size()); 

    // Timer timer;

    if (missing_samples.size() > 0) {
        fprintf(stderr, "Adding missing samples to the tree.\n");  
        
        for (size_t s=0; s<missing_samples.size(); s++) {
             timer.Start();
            auto sample = missing_samples[s];

            auto dfs = T.depth_first_expansion();
            size_t total_nodes = dfs.size();
                        
            std::vector<std::vector<mutation>> node_excess_mutations(total_nodes);
            std::vector<std::vector<mutation>> node_imputed_mutations(total_nodes);

            size_t best_level = 1e9;
            int best_set_difference = 1e9;
            size_t best_j = 0;
            Node* best_node = NULL;

#pragma omp parallel for
            for (size_t k = 0; k < total_nodes; k++) {
                mapper2_input inp;
                inp.T = &T;
                inp.node = dfs[k];
                inp.node_mutations = &node_mutations;
                inp.missing_sample_mutations = &missing_sample_mutations[s];
                inp.excess_mutations = &node_excess_mutations[k];
                inp.imputed_mutations = &node_imputed_mutations[k];
                inp.best_level = &best_level;
                inp.best_set_difference = &best_set_difference;
                inp.best_node = &best_node;
                inp.best_j =  &best_j;
                inp.j = k;

                mapper2_body(inp);
            }
            
            fprintf(stderr, "Current tree size (#nodes): %zu\tMissing sample: %s\tParsimony score: %d\n", total_nodes, sample.c_str(), \
                    best_set_difference);

            if (T.get_node(sample) == NULL) {
                if (best_node->is_leaf()) {
                    std::string nid = std::to_string(++T.curr_internal_node);
                    T.create_node(nid, best_node->parent->identifier);
                    T.create_node(sample, nid);
                    T.move_node(best_node->identifier, nid);
                    std::vector<mutation> common_mut, l1_mut, l2_mut;
                    std::vector<mutation> curr_l1_mut;

                    if (node_mutations.find(best_node) != node_mutations.end()) {
                        for (auto m1: node_mutations[best_node]) {
                            mutation m;
                            m.position = m1.position;
                            m.ref_nuc = m1.ref_nuc;
                            m.par_nuc = m1.par_nuc;
                            m.mut_nuc.emplace_back(m1.mut_nuc[0]);
                            curr_l1_mut.emplace_back(m);
                        }
                        node_mutations.erase(best_node);
                    }

                    for (auto m1: curr_l1_mut) {
                        bool found = false;
                        for (auto m2: node_excess_mutations[best_j]) {
                            if (m1.position == m2.position) {
                                if (m1.mut_nuc[0] == m2.mut_nuc[0]) {
                                    found = true;
                                    mutation m;
                                    m.position = m1.position;
                                    m.ref_nuc = m1.ref_nuc;
                                    m.par_nuc = m1.par_nuc;
                                    m.mut_nuc.emplace_back(m1.mut_nuc[0]);
                                    common_mut.emplace_back(m);
                                    break;
                                }
                            }
                        }
                        if (!found) {
                            mutation m;
                            m.position = m1.position;
                            m.ref_nuc = m1.ref_nuc;
                            m.par_nuc = m1.par_nuc;
                            m.mut_nuc.emplace_back(m1.mut_nuc[0]);
                            l1_mut.emplace_back(m);
                        }
                    }
                    for (auto m1: node_excess_mutations[best_j]) {
                        bool found = false;
                        for (auto m2: curr_l1_mut) {
                            if (m1.position == m2.position) {
                                if (m1.mut_nuc[0] == m2.mut_nuc[0]) {
                                    found = true;
                                    break;
                                }
                            }
                        }
                        if (!found) {
                            mutation m;
                            m.position = m1.position;
                            m.ref_nuc = m1.ref_nuc;
                            m.par_nuc = m1.par_nuc;
                            m.mut_nuc.emplace_back(m1.mut_nuc[0]);
                            l2_mut.emplace_back(m);
                        }
                    }

                    if (common_mut.size() > 0) {
                        node_mutations[T.get_node(nid)] = common_mut;
                    }
                    if (l1_mut.size() > 0) {
                        node_mutations[T.get_node(best_node->identifier)] = l1_mut;
                    }
                    if (l2_mut.size() > 0) {
                        node_mutations[T.get_node(sample)] = l2_mut;
                    }
                }
                else {
                    T.create_node(sample, best_node->identifier);
                    Node* node = T.get_node(sample);
                    std::vector<mutation> node_mut;
                    for (auto mut: node_excess_mutations[best_j]) {
                        mutation m;
                        m.position = mut.position;
                        m.ref_nuc = mut.ref_nuc;
                        m.par_nuc = mut.par_nuc;
                        for (auto nuc: mut.mut_nuc) {
                            m.mut_nuc.emplace_back(nuc);
                        }
                        node_mut.emplace_back(m);
                    }
                    node_mutations[node] = node_mut;
                }

                fprintf (stderr, "Imputed mutations:\t");
                size_t tot = node_imputed_mutations[best_j].size();
                for (size_t curr = 0; curr < tot; curr++) {
                    if (curr < tot-1) {
                        fprintf (stderr, "%i:%c;", node_imputed_mutations[best_j][curr].position, get_nuc_char(node_imputed_mutations[best_j][curr].mut_nuc[0]));
                    }
                    else {
                        fprintf (stderr, "%i:%c", node_imputed_mutations[best_j][curr].position, get_nuc_char(node_imputed_mutations[best_j][curr].mut_nuc[0]));
                    }
                }
                fprintf(stderr, "\n");
            }

            fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
        }

    }

    if (collapse_tree) {
        fprintf(stderr, "Collapsing final tree. \n");
        timer.Start();

        bfs.clear();
        bfs = T.breadth_first_expansion();
        
        for (size_t idx = 1; idx < bfs.size(); idx++) {
            auto mutations = node_mutations[bfs[idx]];
            if (mutations.size() == 0) {
                auto node = bfs[idx];
                auto parent = node->parent;
                auto children = node->children;
                for (auto child: children) {
                    T.move_node(child->identifier, parent->identifier);
                }
            }
        }
        
        condensed_T = create_tree_from_newick_string(get_newick_string(T, false, true));
        
        fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
        
        fprintf(stderr, "Condensing identical sequences. \n");
        timer.Start();

        bfs.clear();
        bfs = T.breadth_first_expansion();
        auto condensed_bfs = condensed_T.breadth_first_expansion();

        assert(condensed_bfs.size() == bfs.size());

        condensed_nodes.clear();

        for (size_t it = 0; it < condensed_bfs.size(); it++) {
            auto condensed_node = condensed_bfs[it];
            condensed_node_mutations.insert(std::pair<Node*, std::vector<mutation>>(condensed_node, std::vector<mutation>(node_mutations[bfs[it]].size())));
            for (size_t k = 0; k < node_mutations[bfs[it]].size(); k++) {
                condensed_node_mutations[condensed_node][k] = node_mutations[bfs[it]][k];
            }
        }

        auto tree_leaves = T.get_leaves();
        for (auto l1: tree_leaves) {
            std::vector<std::string> polytomy_nodes;

            if (std::find(missing_samples.begin(), missing_samples.end(), l1->identifier) != missing_samples.end()) {
                continue;
            }
            if (node_mutations[l1].size() > 0) {
                continue;
            }
            if (condensed_T.get_node(l1->identifier) == NULL) {
                continue;
            }

            for (auto l2: l1->parent->children) {
                if (std::find(missing_samples.begin(), missing_samples.end(), l2->identifier) != missing_samples.end()) {
                    continue;
                }
                if (l2->is_leaf() && (condensed_T.get_node(l2->identifier) != NULL) && (node_mutations[l2].size() == 0)) {
                    polytomy_nodes.push_back(l2->identifier);
                }
            }

            if (polytomy_nodes.size() > 1) {
                std::string new_node_name = "node_" + std::to_string(1+condensed_nodes.size()) + "_condensed_" + std::to_string(polytomy_nodes.size()) + "_leaves";
                auto curr_node = condensed_T.get_node(l1->identifier);
                condensed_T.create_node(new_node_name, curr_node->parent->identifier, l1->branch_length);
                auto new_node = condensed_T.get_node(new_node_name);
                condensed_node_mutations.insert(std::pair<Node*, std::vector<mutation>>(new_node, std::vector<mutation>(0)));
                condensed_nodes[new_node_name] = std::vector<std::string>(polytomy_nodes.size());

                for (size_t it = 0; it < polytomy_nodes.size(); it++) {
                    condensed_nodes[new_node_name][it] = polytomy_nodes[it];
                    condensed_T.remove_node(polytomy_nodes[it], false);
                }
            }
        }
        fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
        
        fprintf(stderr, "Printing condensed tree. \n");
        fprintf(stdout, "%s\n", get_newick_string(condensed_T, true, true).c_str());
    }

    fprintf(stderr, "Printing final tree. \n");
    fprintf(stdout, "%s\n", get_newick_string(T, true, true).c_str());

    if (print_uncondensed_tree) {
        fprintf(stderr, "Printing uncondensed final tree. \n");
        
        timer.Start();
        
        if (!collapse_tree && (condensed_nodes.size() > 0)) {
            Tree T_to_print = create_tree_from_newick_string(get_newick_string(T, false, true)); 
            for (size_t it = 0; it < condensed_nodes.size(); it++) {
                auto cn = condensed_nodes.begin();
                std::advance(cn, it);

                auto n = T_to_print.get_node(cn->first);
                auto par = (n->parent != NULL) ? n->parent : n;

                size_t num_samples = cn->second.size();

                if (num_samples > 0) {
                    T_to_print.rename_node(n->identifier, cn->second[0]);
                }
                
                for (size_t s = 1; s < num_samples; s++) {
                    T_to_print.create_node(cn->second[s], par->identifier, n->branch_length);
                }
            }
            fprintf(stdout, "%s\n", get_newick_string(T_to_print, true, true).c_str());
        }
        else {
            fprintf(stdout, "%s\n", get_newick_string(T, true, true).c_str());
        }
        
        fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
    }
    
    if (missing_samples.size() > 0) {
        fprintf(stderr, "Printing mutation paths in the new tree for missing samples.\n");  

        timer.Start();
        
        for (size_t s=0; s<missing_samples.size(); s++) {
            auto sample = missing_samples[s];
            auto sample_node = T.get_node(sample);
            std::stack<std::string> mutation_stack;
            std::string curr_node_mutation_string;

            auto curr_node_mutations = node_mutations[sample_node];
            if (curr_node_mutations.size() > 0) {
                curr_node_mutation_string = sample + ":";
                size_t num_mutations = curr_node_mutations.size();
                for (size_t k = 0; k < num_mutations; k++) {
                    curr_node_mutation_string += get_nuc_char(curr_node_mutations[k].par_nuc) + std::to_string(curr_node_mutations[k].position) + get_nuc_char(curr_node_mutations[k].mut_nuc[0]); 
                    if (k < num_mutations-1) {
                        curr_node_mutation_string += ',';
                    }
                    else {
                        curr_node_mutation_string += ';';    
                    }
                }
                mutation_stack.push(curr_node_mutation_string);
            }

            for (auto anc_node: T.rsearch(sample)) {
                curr_node_mutations = node_mutations[anc_node];
                if (curr_node_mutations.size() > 0) {
                    curr_node_mutation_string = anc_node->identifier + ":";
                    size_t num_mutations = curr_node_mutations.size();
                    for (size_t k = 0; k < num_mutations; k++) {
                        curr_node_mutation_string += get_nuc_char(curr_node_mutations[k].par_nuc) + std::to_string(curr_node_mutations[k].position) + get_nuc_char(curr_node_mutations[k].mut_nuc[0]); 
                        if (k < num_mutations-1) {
                            curr_node_mutation_string += ',';
                        }
                        else {
                            curr_node_mutation_string += ';';    
                        }
                    }
                    mutation_stack.push(curr_node_mutation_string);
                }
            }
            
            fprintf(stderr, "%s\t", sample.c_str()); 
            while (mutation_stack.size()) {
                fprintf(stderr, "%s", mutation_stack.top().c_str()); 
                mutation_stack.pop();
            }
            fprintf(stderr, "\n"); 
        }
        fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
    }

    if ((print_subtrees_size > 1) && (missing_samples.size() > 0)) {
        fprintf(stderr, "Printing subtrees for display. \n");

        timer.Start();
        
        std::vector<bool> displayed_mising_sample (missing_samples.size(), false);
        
        for (size_t i = 0; i < missing_samples.size(); i++) {
            if (displayed_mising_sample[i]) {
                continue;
            }
            
            for (auto anc: T.rsearch(missing_samples[i])) {
                size_t num_leaves = T.get_leaves(anc->identifier).size();
                if (num_leaves < print_subtrees_size) {
                    continue;
                }

                std::string newick = get_newick_string(T, anc, false, true);
                Tree new_T = create_tree_from_newick_string(newick);

                if (num_leaves > print_subtrees_size) {
                    auto last_anc = new_T.get_node(missing_samples[i]);
                    auto ancestors = new_T.rsearch(missing_samples[i]);
                    if (ancestors.size() > 1) {
                        last_anc = ancestors[ancestors.size()-2];
                    }
                    std::vector<Node*> siblings;
                    for (auto child: new_T.root->children) {
                        if (child->identifier != last_anc->identifier) {
                            siblings.push_back(child);
                        }
                    }
                    
                    for (size_t k=0; k<siblings.size(); k++) {
                        auto curr_sibling = siblings[k];
                        auto sibling_leaves = new_T.get_leaves(curr_sibling->identifier);
                        if (num_leaves-sibling_leaves.size() < print_subtrees_size) {
                            for (auto child: curr_sibling->children) {
                                siblings.push_back(child);
                            }
                        }
                        else {
                            new_T.remove_node(curr_sibling->identifier, true);
                            num_leaves -= sibling_leaves.size();
                            if (num_leaves == print_subtrees_size) {
                                break;
                            }
                        }
                    }

                    newick = get_newick_string(new_T, false, true);
                }

#pragma omp parallel for
                for (size_t j = i+1; j < missing_samples.size(); j++) {
                    if (!displayed_mising_sample[j]) {
                        if (new_T.get_node(missing_samples[j]) != NULL) {
                            displayed_mising_sample[j] = true;
                        }
                    }
                }
                
                fprintf(stdout, "%s\n", newick.c_str());
                break;
            }
        }
        fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
    }

    if (dout_filename != "") {
        fprintf(stderr, "Saving assignments. \n");

        timer.Start();

        Parsimony::data data;

        if (!collapse_tree) {
            data.set_newick(get_newick_string(T, false, true));

            auto dfs = T.depth_first_expansion();

            for (size_t idx = 0; idx < dfs.size(); idx++) {
                auto mutation_list = data.add_node_mutations();
                auto mutations = node_mutations[dfs[idx]];
                for (auto m: mutations) {
                    auto mut = mutation_list->add_mutation();
                    mut->set_position(m.position);
                    mut->set_ref_nuc(m.ref_nuc);
                    mut->set_par_nuc(m.par_nuc);
                    mut->clear_mut_nuc();
                    for (auto nuc: m.mut_nuc) {
                        mut->add_mut_nuc(nuc);
                    }
                }
            }
        }
        else {
            data.set_newick(get_newick_string(condensed_T, false, true));

            auto dfs = condensed_T.depth_first_expansion();

            for (size_t idx = 0; idx < dfs.size(); idx++) {
                auto mutation_list = data.add_node_mutations();
                auto mutations = condensed_node_mutations[dfs[idx]];
                for (auto m: mutations) {
                    auto mut = mutation_list->add_mutation();
                    mut->set_position(m.position);
                    mut->set_ref_nuc(m.ref_nuc);
                    mut->set_par_nuc(m.par_nuc);
                    mut->clear_mut_nuc();
                    for (auto nuc: m.mut_nuc) {
                        mut->add_mut_nuc(nuc);
                    }
                }
            }

            // Add condensed nodes
            for (auto cn: condensed_nodes) {
                auto cn_ptr = data.add_condensed_nodes();
                cn_ptr->set_node_name(cn.first);
                for (auto l: cn.second) {
                    cn_ptr->add_condensed_leaves(l);
                }
            }
        }

        std::ofstream outfile(dout_filename, std::ios::out | std::ios::binary);
        data.SerializeToOstream(&outfile);
        outfile.close();
        
        fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
    }
    
    google::protobuf::ShutdownProtobufLibrary();

#if SAVE_PROFILE == 1
    Instrumentor::Get().EndSession();
#endif

    return 0;
}

