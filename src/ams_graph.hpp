#include "tree.hpp"
#include <set>
#include <mutex>
#include <sys/time.h>

extern std::mutex data_lock;

class Timer {
    private:
        struct timeval start_time, end_time;
    public:
        void start() {
            gettimeofday(&start_time, NULL);
        }
        long stop() {
            long useconds, seconds, mseconds;
            gettimeofday(&end_time, NULL);
            useconds = end_time.tv_usec - start_time.tv_usec;
            seconds = end_time.tv_sec - start_time.tv_sec;
            mseconds = ((seconds) * 1000 + useconds/1000.0 + 0.5);
            return mseconds;
        }
        
};

struct mutation {
    int position;
    int8_t ref_nuc;
    std::vector<int8_t> mut_nuc;
};

struct mapper_input {
    Tree* T;
    int8_t ref_nuc;
    int variant_pos;
    std::vector<Node*>* bfs;
    std::unordered_map<std::string, size_t>* bfs_idx;
    std::vector<std::tuple<size_t, std::vector<int8_t>>> variants;
    std::vector<std::string>* variant_ids;
    
    std::vector<std::string>* missing_samples;
    std::unordered_map<Node*, std::vector<mutation>>* node_mutations;
    std::vector<std::vector<mutation>>* missing_sample_mutations;
};

struct mapper_body {
    int operator()(mapper_input input);
};

struct mapper2_input {
    std::string missing_sample;
    Tree* T;
    Node* node;
    std::unordered_map<Node*, std::vector<mutation>>* node_mutations;
    std::vector<mutation>* missing_sample_mutations;
    
    int* best_set_difference;
    size_t* best_level;
    size_t j;
    size_t* best_j;
    Node** best_node;

    std::vector<mutation>* excess_mutations;
};

struct mapper2_body {
    int operator()(mapper2_input input);
};

