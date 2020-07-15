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
//#include "../build/parsimony.pb.h"
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

int main(int argc, char** argv){

    std::string tree_filename;
    std::string din_filename;
    std::string dout_filename;
    std::string vcf_filename;
    uint32_t num_cores = tbb::task_scheduler_init::default_num_threads();
    uint32_t num_threads;
    po::options_description desc{"Options"};
    desc.add_options()
        ("vcf", po::value<std::string>(&vcf_filename)->required(), "Input VCF file (in uncompressed or gzip-compressed format)")
        ("tree", po::value<std::string>(&tree_filename)->default_value(""), "Input tree file")
        ("load-assignments", po::value<std::string>(&din_filename)->default_value(""), "Load existing tree and parsimonious assignments")
        ("save-assignments", po::value<std::string>(&dout_filename)->default_value(""), "Save output tree and parsimonious assignments")
        ("threads", po::value<uint32_t>(&num_threads)->default_value(num_cores), "Number of threads")
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


    Tree T;

    bool header_found = false;
    std::vector<std::string> variant_ids;
    std::vector<std::string> missing_samples;
    std::unordered_map<Node*, std::vector<mutation>> node_mutations;
    std::vector<std::vector<mutation>> missing_sample_mutations;
    size_t num_missing = 0;

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
    }
    else if (din_filename != "") {
        fprintf(stderr, "Loading existing assignments nad VCF file\n");
        
        Parsimony::data data;

        std::ifstream inpfile(din_filename, std::ios::in | std::ios::binary);
        data.ParseFromIstream(&inpfile);
        inpfile.close();

        T = create_tree_from_newick_string(data.newick());

        bfs.clear();
        bfs = T.breadth_first_expansion();
        
        for (size_t idx = 0; idx < bfs.size(); idx++) {
            bfs_idx[bfs[idx]->identifier] = idx;
        }

        for (size_t idx = 0; idx < bfs.size(); idx++) {
            auto mutation_list = data.node_mutations(idx);
            auto node = bfs[idx];
            node_mutations.insert(std::pair<Node*, std::vector<mutation>>(node, std::vector<mutation>()));  
            for (int k = 0; k < mutation_list.mutation_size(); k++) {
                auto mut = mutation_list.mutation(k);
                mutation m;
                m.position = mut.position();
                m.ref_nuc = mut.ref_nuc();
                for (int n = 0; n < mut.mut_nuc_size(); n++) {
                    m.mut_nuc.emplace_back(mut.mut_nuc(n));
                }
                node_mutations[node].emplace_back(m);
            }
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
        std::string s;
        while (instream.peek() != EOF) {
            std::getline(instream, s);
            std::vector<std::string> words;
            split(s, words);
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
                split(words[4], ',', alleles);
                for (size_t j=9; j < words.size(); j++) {
                    auto iter = std::find(missing_samples.begin(), missing_samples.end(), variant_ids[j-9]);
                    if (iter != missing_samples.end()) {
                        auto mutations_iter = missing_sample_mutations.begin() + (iter - missing_samples.begin());
                        mutation m;
                        m.position = std::stoi(words[1]);
                        auto ref_nucs = get_nuc_id(words[3][0]);
                        assert(ref_nucs.size() == 1);
                        m.ref_nuc = ref_nucs[0];
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
//                        else {
//                            for (auto n: get_nuc_id('N')) {
//                                m.mut_nuc.emplace_back(n);
//                            }
//                        }
                    }
                }
            }
        }
    }
    else {
        fprintf(stderr, "Error! No input tree or assignment file provided!\n");
        exit(1);
    }

    fprintf(stderr, "Computing parsimonious assignments done.\nFound %zu missing samples.\n", missing_samples.size()); 

    if (missing_samples.size() > 0) {
        fprintf(stderr, "Adding missing samples to the tree.\n");  
        
        for (size_t s=0; s<missing_samples.size(); s++) {
            auto sample = missing_samples[s];

            size_t curr = 0;
            auto dfs = T.depth_first_expansion();
            size_t total_nodes = dfs.size();
                        
            std::vector<std::vector<mutation>> node_excess_mutations(total_nodes);
            int* node_set_difference = (int*) calloc(total_nodes, sizeof(int));

            tbb::flow::graph mapper_graph2;

            tbb::flow::function_node<mapper2_input, int> mapper2(mapper_graph2, tbb::flow::unlimited, mapper2_body());
            tbb::flow::source_node <mapper2_input> splitter (mapper_graph2,
                    [&] (mapper2_input &inp) -> bool {
                    curr += 1;
                    inp.missing_sample = sample;
                    inp.T = &T;
                    inp.node = dfs[curr-1];
                    inp.node_mutations = &node_mutations;
                    inp.missing_sample_mutations = &missing_sample_mutations[s];
                    inp.excess_mutations = &node_excess_mutations[curr-1];
                    inp.set_difference = &node_set_difference[curr-1];
                    if (curr <= total_nodes) {
                        return true;
                    }
                    else {
                        return false;
                    }
                    }, true );
            tbb::flow::make_edge(splitter, mapper2);
            mapper_graph2.wait_for_all();

            Node* best_node;
            size_t best_level = 1e9;
            int best_set_difference = 1e9;
            size_t best_j = 0;
            for (size_t j=0; j< total_nodes; j++) {
                auto n = dfs[j];
                if (node_set_difference[j] < best_set_difference) {
                    best_set_difference = node_set_difference[j];
                    best_node = n;
                    best_level = n->level;
                    best_j = j;
                }
                else if (node_set_difference[j] == best_set_difference) {
                    if (n->level < best_level) {
                        best_set_difference = node_set_difference[j];
                        best_node = n;
                        best_level = n->level;
                        best_j = j;
                    }
                }
            }
            fprintf(stderr, "Current tree size (#nodes): %zu\tMissing sample: %s\tParsimony score: %d\n", total_nodes, sample.c_str(), \
                    best_set_difference);
//                    missing_sample_mutations[s].size(), best_node->identifier.c_str(), best_set_difference);
//            for (auto mut: node_excess_mutations[best_j]) {
//                fprintf(stderr, "%d:%d,", mut.position, mut.mut_nuc[0]); 
//            }
//            fprintf(stderr, "\n");

            if (T.get_node(sample) == NULL) {
                if (best_node->is_leaf()) {
                    std::string nid = std::to_string(++T.curr_internal_node);
                    T.create_node(nid, nid, best_node->parent->identifier);
                    T.create_node(sample, sample, nid);
                    T.move_node(best_node->identifier, nid);
                    std::vector<mutation> common_mut, l1_mut, l2_mut;
                    std::vector<mutation> curr_l1_mut;

                    if (node_mutations.find(best_node) != node_mutations.end()) {
                        for (auto m1: node_mutations[best_node]) {
                            mutation m;
                            m.position = m1.position;
                            m.ref_nuc = m1.ref_nuc;
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
                    T.create_node(sample, sample, best_node->identifier);
                    Node* node = T.get_node(sample);
                    std::vector<mutation> node_mut;
                    for (auto mut: node_excess_mutations[best_j]) {
                        mutation m;
                        m.position = mut.position;
                        m.ref_nuc = mut.ref_nuc;
                        for (auto nuc: mut.mut_nuc) {
                            m.mut_nuc.emplace_back(nuc);
                        }
                        node_mut.emplace_back(m);
                    }
                    node_mutations[node] = node_mut;
                }
            }

            free(node_set_difference);
        }

    }

    fprintf(stderr, "Displaying final tree: \n");
    fprintf(stdout, "%s\n", get_newick_string(T, false).c_str());

    if (dout_filename != "") {
        fprintf(stderr, "Saving assignments. \n");

        Parsimony::data data;
        data.set_newick(get_newick_string(T, false));

        bfs.clear();
        bfs = T.breadth_first_expansion();

        for (size_t idx = 0; idx < bfs.size(); idx++) {
            auto mutation_list = data.add_node_mutations();
            auto mutations = node_mutations[bfs[idx]];
            for (auto m: mutations) {
                auto mut = mutation_list->add_mutation();
                mut->set_position(m.position);
                mut->set_ref_nuc(m.ref_nuc);
                mut->clear_mut_nuc();
                for (auto nuc: m.mut_nuc) {
                    mut->add_mut_nuc(nuc);
                }
            }
        }

        std::ofstream outfile(dout_filename, std::ios::out | std::ios::binary);
        data.SerializeToOstream(&outfile);
        outfile.close();
    }
    
    google::protobuf::ShutdownProtobufLibrary();
    
    return 0;
}

