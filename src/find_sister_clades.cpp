#include <fstream>
#include <boost/program_options.hpp> 
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <iostream>
#include <omp.h>
#include "tree.hpp"

namespace po = boost::program_options;

int main(int argc, char** argv){

    std::string tree_filename;
    std::string samples_filename;
    po::options_description desc{"Options"};
    desc.add_options()
        ("tree", po::value<std::string>(&tree_filename)->required(), "Input tree file")
        ("samples", po::value<std::string>(&samples_filename)->required(), "File containing missing samples.")
        ("help", "Print help messages");
    
    po::options_description all_options;
    all_options.add(desc);

    po::positional_options_description p;
    p.add("tree", 1);
    p.add("samples", 1);

    po::variables_map vm;
    try{
        po::store(po::command_line_parser(argc, argv).options(all_options).positional(p).run(), vm);
        po::notify(vm);
    }
    catch(std::exception &e){
        std::cout << desc << std::endl;
        return 1;
    }

    Tree T = create_tree_from_newick(tree_filename);

    std::ifstream infile;
    infile.open(samples_filename);

    std::vector<std::string> missing_samples;
    std::string sample;

    while (std::getline(infile, sample)) {
        assert (T.get_node(sample) != NULL);
        missing_samples.push_back(sample);
    }

    fprintf(stderr, "Printing sibling clade leaves for each sample.\n");

    for (auto s: missing_samples) {
        auto last_anc = T.get_node(s);
        std::vector<std::string> siblings;
        for (auto anc: T.rsearch(s)) {
            for (auto child: anc->children) {
                if (child != last_anc) {
                    for (auto l: T.get_leaves(child->identifier)) {
                        if (std::find(missing_samples.begin(), missing_samples.end(), l->identifier) == missing_samples.end()) {
                            siblings.push_back(l->identifier);
                        }

                    }
                }
            }
            last_anc = anc;
            if (siblings.size() > 0) {
                break;
            }
        }
        fprintf(stdout, "%s: \n", s.c_str());
        for (auto sibling: siblings) {
            fprintf(stdout, "%s\n", sibling.c_str());
        }
        fprintf(stdout, "\n");
    }

    return 0;
}

