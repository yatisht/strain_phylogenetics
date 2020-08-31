#include <fstream>
#include <boost/program_options.hpp> 
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <iostream>
#include <algorithm>
#include <tbb/mutex.h>
#include <mutex>
#include <tbb/blocked_range.h>
#include <tbb/parallel_sort.h>
#include <tbb/tbb.h>
#include <unordered_map>
#include "tree.hpp"

namespace po = boost::program_options;

std::vector<Node*> get_ordered_leaves(Tree T) {
    auto dfs = T.depth_first_expansion();
    std::vector<Node*> ordered_leaves;
    for (auto n: dfs) {
        if (n->is_leaf()) {
            ordered_leaves.push_back(n);
        }
    }
    return ordered_leaves;
}

float get_average_rank (std::vector<Node*> nodes, std::unordered_map<std::string, size_t> rank_dict) {
    float num = 0;
    float den = 1e-9;

    for (auto n: nodes) {
        if (rank_dict.find(n->identifier) != rank_dict.end()) {
            num += static_cast<float>(rank_dict[n->identifier]);
            den += 1.0;
        }
    }

    return (num/den);
}

struct Node_rank {
    Node* node;
    float rank;
    Node_rank(Node* n, float r) : node(n), rank(r)
    {
    }
};

bool compare_by_rank (const Node_rank& nr1, const Node_rank nr2) {
    return nr1.rank < nr2.rank;
}

uint32_t rotate_T2 (std::vector<Node*> T1_ordered_leaves, Tree& T2) {

    std::unordered_map<std::string, size_t> rank_dict;
    
    for (size_t k=0; k<T1_ordered_leaves.size(); k++) {
        rank_dict[T1_ordered_leaves[k]->identifier] = k;
    }

    uint32_t num_rot = 0;

    auto dfs = T2.depth_first_expansion();
    size_t num_nodes = dfs.size();

    tbb::mutex data_lock;

    tbb::parallel_for( tbb::blocked_range<size_t>(0, num_nodes, 20),
            [&](tbb::blocked_range<size_t> r) {

            for (size_t k=r.begin(); k<r.end(); ++k) {
               auto curr_node = dfs[k];
               std::vector<Node_rank> children_ranks;
               for (auto c: curr_node->children) {
                   data_lock.lock();
                   auto leaves = T2.get_leaves(c->identifier);
                   data_lock.unlock();
                   float rank = get_average_rank(leaves, rank_dict);
                   children_ranks.push_back(Node_rank(c, rank));
               }
               if (!is_sorted(children_ranks.begin(), children_ranks.end(), compare_by_rank)) {
                    std::sort(children_ranks.begin(), children_ranks.end(), compare_by_rank); 
                    
                    data_lock.lock();
                    num_rot++;
                    curr_node->children.clear();
                    for (auto cr: children_ranks) {
                        curr_node->children.push_back(cr.node);
                    }
                    data_lock.unlock();
               }

            }
            });

    return num_rot;
}

int main(int argc, char** argv){

    std::string T1_filename;
    std::string T2_filename;
    std::string T1_out_filename;
    std::string T2_out_filename;
    uint32_t max_iter;

    po::options_description desc{"Options"};
    desc.add_options()
        ("T1", po::value<std::string>(&T1_filename)->required(), "Input tree1 file")
        ("T2", po::value<std::string>(&T2_filename)->required(), "Input tree2 file")
        ("T1_out", po::value<std::string>(&T1_out_filename)->required(), "Input tree1 file")
        ("T2_out", po::value<std::string>(&T2_out_filename)->required(), "Input tree2 file")
        ("max_iter", po::value<uint32_t>(&max_iter)->default_value(1e2), "Maximum number of iterations")
        ("help", "Print help messages");
    
    po::options_description all_options;
    all_options.add(desc);

    po::positional_options_description p;
    p.add("T1", 1);
    p.add("T2", 1);
    p.add("T1_out", 1);
    p.add("T2_out", 1);
    p.add("max_iter", 1);
    
    po::variables_map vm;
    try{
        po::store(po::command_line_parser(argc, argv).options(all_options).positional(p).run(), vm);
        po::notify(vm);
    }
    catch(std::exception &e){
        std::cout << desc << std::endl;
        if(vm.count("help"))
            return 0;
        else
            return 1;
    }
        
    uint32_t num_cores = tbb::task_scheduler_init::default_num_threads();
    tbb::task_scheduler_init init(num_cores);
    
    Tree T1, T2, T1_out, T2_out;
    
    T1 = TreeLib::create_tree_from_newick(T1_filename);
    T2 = TreeLib::create_tree_from_newick(T2_filename);

    uint32_t n1=1, n2=1;
    uint32_t num_iter=0;

    while ((n1+n2 > 0) && (num_iter < max_iter)) {
        auto T1_ordered_leaves = get_ordered_leaves(T1);
        n1 = rotate_T2(T1_ordered_leaves, T2);
        
        auto T2_ordered_leaves = get_ordered_leaves(T2);
        n2 = rotate_T2(T2_ordered_leaves, T1);
        
        num_iter++;
        fprintf(stderr, "===Iteration %u===\nPerformed %u rotations.\n", num_iter, n1+n2);
    }

    if ((num_iter == max_iter) && (n1+n2>0)) {
        fprintf(stderr, "\nWARNING: Reached maximum iterations before converging.\n");
    }
    
    std::ofstream ofile1, ofile2;
    ofile1.open(T1_out_filename);
    ofile2.open(T2_out_filename);

    ofile1 << TreeLib::get_newick_string(T1, false, true);
    ofile2 << TreeLib::get_newick_string(T2, false, true);

    ofile1.close();
    ofile2.close();

    return 0;
}

