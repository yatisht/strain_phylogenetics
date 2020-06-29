#include <fstream>
#include <boost/program_options.hpp> 
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <iostream>
#include <tbb/flow_graph.h>
#include <tbb/reader_writer_lock.h>
#include <tbb/scalable_allocator.h>
#include <tbb/task_scheduler_init.h>
//#include "zlib.h"
#include "fpa_graph.hpp"

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
    std::string vcf_filename;
    uint32_t num_cores = tbb::task_scheduler_init::default_num_threads();
    uint32_t num_threads;
    bool print_vcf = false;
    po::options_description desc{"Options"};
    desc.add_options()
        ("tree", po::value<std::string>(&tree_filename)->required(), "Input tree file")
        ("vcf", po::value<std::string>(&vcf_filename)->required(), "Input VCF file (in uncompressed or gzip-compressed format)")
        ("threads", po::value<uint32_t>(&num_threads)->default_value(num_cores), "Number of threads")
        ("print-vcf", po::bool_switch(&print_vcf), "Print VCF with variants resolved instead of printing a parsimony file.")
        ("help", "Print help messages");
    
    po::options_description all_options;
    all_options.add(desc);

    po::positional_options_description p;
    p.add("tree", 1);
    p.add("vcf", 1);
    p.add("threads", 1);

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

    auto bfs = T.breadth_first_expansion();
    std::unordered_map<std::string, size_t> bfs_idx;
    for (size_t idx = 0; idx < bfs.size(); idx++) {
        bfs_idx[bfs[idx]->identifier] = idx;
    }

    if (!print_vcf) {
        std::cout << get_newick_string(T, true) << std::endl; 
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
    
    bool header_found = false;
    std::vector<std::string> variant_ids;

    tbb::task_scheduler_init init(num_threads);

    tbb::flow::graph map_reduce_graph;
    
    tbb::flow::function_node<printer_input, size_t> printer(map_reduce_graph, tbb::flow::unlimited, printer_body());
    tbb::flow::function_node<mapper_input, printer_input> mapper(map_reduce_graph, tbb::flow::unlimited, mapper_body());
    tbb::flow::source_node <mapper_input> reader (map_reduce_graph,
            [&] (mapper_input &inp) -> bool {
            std::string s;
            std::getline(instream, s);
            std::vector<std::string> words;
            split(s, words);
            inp.variant = "";
            inp.print_vcf = print_vcf;
            if (not header_found) {
              inp.line = s;
            }
            if ((not header_found) && (words.size() > 1)) {
              if (words[1] == "POS") {
                  for (size_t j=9; j < words.size(); j++) {
                    variant_ids.push_back(words[j]);
                  }
                  header_found = true;
              }
            }
            else if (header_found) {
              if (words.size() != 9+variant_ids.size()) {
                fprintf(stderr, "ERROR! Incorrect VCF format.\n");
                exit(1);
              }
              if (print_vcf) {
                  inp.words.clear();
                  for (int j=0; j<9; j++) {
                      inp.words.push_back(words[j]);
                  }
              }
              std::vector<std::string> alleles;
              alleles.clear();
              inp.variant = words[3] + words[1] + words[4];
              split(words[4], ',', alleles);
              inp.T = &T;
              inp.bfs = &bfs;
              inp.bfs_idx = &bfs_idx;
              inp.variant_ids = &variant_ids;
              auto ref_nucs = get_nuc_id(words[3][0]);
              assert(ref_nucs.size() == 1);
              inp.ref_nuc = ref_nucs[0]; 
              inp.variants.clear();
              for (size_t j=9; j < words.size(); j++) {
                  if (isdigit(words[j][0])) {
                     int allele_id = std::stoi(words[j]);
                     if (allele_id > 0) { 
                        std::string allele = alleles[allele_id-1];
                        inp.variants.push_back(std::make_tuple(j-9, get_nuc_id(allele[0])));
                     }
                  }
                  else {
                    inp.variants.push_back(std::make_tuple(j-9, get_nuc_id('N')));
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
    tbb::flow::make_edge(mapper, printer);
    map_reduce_graph.wait_for_all();

    return 0;
}

