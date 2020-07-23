#include <fstream>
#include <unordered_map>
#include <string>
#include <sstream>
#include <vector>
#include <queue>
#include <stack>
#include <algorithm>
#include <cassert>

class Node {
    public:
        size_t level;
        std::string identifier;
        std::string tag;
        Node* parent;
        std::vector<Node*> children;

        bool is_leaf();
        bool is_root();
        Node();
        Node(std::string id, std::string t);
        Node(std::string id, std::string t, Node* p);
};

class Tree {
    private:
        void remove_node_helper (std::string nid);
        void depth_first_expansion_helper(Node* node, std::vector<Node*>& vec);
        std::unordered_map <std::string, Node*> all_nodes;
    public:
        Tree() {
            max_level = 0;
            root = NULL;
            all_nodes.clear();
        }
        
        //TODO
        Tree (Node* n);

        size_t max_level;
        Node* root;

        size_t curr_internal_node;
        size_t get_max_level ();
        std::vector<Node*> get_leaves();
        std::vector<Node*> get_leaves(std::string nid);
        void create_node (std::string identifier, std::string tag);
        void create_node (std::string identifier, std::string tag, std::string parent_id);
        Node* get_node (std::string identifier);
        bool is_ancestor (std::string anc_id, std::string nid);
        std::vector<Node*> rsearch (std::string nid);
        void remove_node (std::string nid);
        void move_node (std::string source, std::string destination);
        std::vector<Node*> breadth_first_expansion();
        std::vector<Node*> breadth_first_expansion(std::string nid);
        std::vector<Node*> depth_first_expansion();
        std::vector<Node*> depth_first_expansion(Node* node);
};

std::string get_newick_string(Tree& T, bool b);
std::string get_newick_string(Tree& T, Node* node, bool b);
Tree create_tree_from_newick (std::string filename);
Tree create_tree_from_newick_string (std::string newick_string);
void split(std::string s, char delim, std::vector<std::string>& words);
void split(std::string s, std::vector<std::string>& words);

