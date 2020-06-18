#include "tree.hpp"

bool Node::is_leaf () {
    return (children.size() == 0);
}

bool Node::is_root() {
    return (parent == NULL);
}

Node::Node() {
    level = 0;
    identifier = "";
    tag = "";
    parent = NULL;
}

Node::Node (std::string id, std::string t) {
    identifier = id;
    tag = t;
    parent = NULL;
    level = 1;
}

Node::Node (std::string id, std::string t, Node* p) {
    identifier = id;
    tag = t;
    parent = p;
    level = p->level + 1;
}

size_t Tree::get_max_level () {
    return max_level;
}

std::vector<Node*> Tree::get_leaves() {
    std::vector<Node*> leaves;
    for (auto x: all_nodes) {
        auto node = x.second;
        if (node->is_leaf()) {
            leaves.push_back(node);
        }
    }
    return leaves;
}

std::vector<Node*> Tree::get_leaves(std::string nid) {
    std::vector<Node*> leaves;
    Node* node = all_nodes[nid];

    std::queue<Node*> remaining_nodes;
    remaining_nodes.push(node);
    while (remaining_nodes.size() > 0) {
        Node* curr_node = remaining_nodes.front();
        if (curr_node->children.size() == 0)
            leaves.push_back(curr_node);
        remaining_nodes.pop();
        for (auto c: curr_node->children) {
            remaining_nodes.push(c);
        }
    }
    return leaves;
}

void Tree::create_node (std::string identifier, std::string tag) {
    all_nodes.clear();
    max_level = 1;
    Node* n = new Node(identifier, tag);
    root = n;
    all_nodes[identifier] = root;
}

void Tree::create_node (std::string identifier, std::string tag, std::string parent_id) {
    Node* par = all_nodes[parent_id];
    Node* n = new Node(identifier, tag, par);
    if (all_nodes.find(identifier) != all_nodes.end()) {
        fprintf(stderr, "Error: %s already in the tree!\n", identifier.c_str());
        exit(1);
    }
    all_nodes[identifier] = n;
    par->children.push_back(n);
    if (n->level > max_level) {
        max_level = n->level;
    }
}

Node* Tree::get_node (std::string nid) {
    if (all_nodes.find(nid) != all_nodes.end()) {
        return all_nodes[nid];
    }
    return NULL;

}

bool Tree::is_ancestor (std::string anc_id, std::string nid) {
    Node* node = all_nodes[nid];
    while (node->parent != NULL) {
        node = node->parent;
        if (node->identifier == anc_id) {
            return true;
        }
    }
    return false; 
}

std::vector<Node*> Tree::rsearch (std::string nid) {
    Node* node = all_nodes[nid];
    std::vector<Node*> ancestors;
    while (node->parent != NULL) {
        ancestors.push_back(node->parent);
        node = node->parent;
    }
    return ancestors;
}

void Tree::remove_node_helper (std::string nid) { 
    Node* source = all_nodes[nid];
    Node* curr_parent = source->parent;
    
    // Remove source from curr_parent
    size_t source_idx = 0;
    for (size_t i = 0; i < curr_parent->children.size(); i++) {
        if (curr_parent->children[i]->identifier == nid) {
            source_idx = i;
            break;
        }
    }
    curr_parent->children.erase(curr_parent->children.begin()+source_idx);
    // Remove parent if it no longer has any children
    if (curr_parent->children.size() == 0) {
        remove_node_helper (curr_parent->identifier);
    }

    //Remove source from all_nodes
    auto it = all_nodes.find(nid);
    all_nodes.erase(it);
}

void Tree::remove_node (std::string nid) { 
    remove_node_helper (nid);

    // Update max level
    size_t new_max_level = 0;
    for (auto x: all_nodes) {
        if (x.second->level > new_max_level) {
            new_max_level = x.second->level;
        }
    }
    max_level = new_max_level;
}

void Tree::move_node (std::string source_id, std::string dest_id) {
    Node* source = all_nodes[source_id];
    Node* destination = all_nodes[dest_id];
    Node* curr_parent = source->parent;

    source->parent = destination;

    // Remove source from curr_parent
    size_t source_idx = 0;
    for (size_t i = 0; i < curr_parent->children.size(); i++) {
        if (curr_parent->children[i]->identifier == source_id) {
            source_idx = i;
            break;
        }
    }
    curr_parent->children.erase(curr_parent->children.begin()+source_idx);
    if (curr_parent->children.size() == 0) {
        remove_node(curr_parent->identifier);
    }
    
    // Update levels of source descendants
    std::queue<Node*> remaining_nodes;
    remaining_nodes.push(source);
    while (remaining_nodes.size() > 0) {
        Node* curr_node = remaining_nodes.front();
        remaining_nodes.pop();
        curr_node->level = curr_node->parent->level + 1;
        for (auto c: curr_node->children) {
            remaining_nodes.push(c);
        }
    }
    
    // Update max level
    size_t new_max_level = 0;
    for (auto x: all_nodes) {
        if (x.second->level > new_max_level) {
            new_max_level = x.second->level;
        }
    }
    max_level = new_max_level;
}

std::vector<Node*> Tree::breadth_first_expansion() {
    return breadth_first_expansion(root->identifier);
}

std::vector<Node*> Tree::breadth_first_expansion(std::string nid) {
    std::vector<Node*> traversal;
    Node* node = all_nodes[nid];

    std::queue<Node*> remaining_nodes;
    remaining_nodes.push(node);
    while (remaining_nodes.size() > 0) {
        Node* curr_node = remaining_nodes.front();
        traversal.push_back(curr_node);
        remaining_nodes.pop();
        for (auto c: curr_node->children) {
            remaining_nodes.push(c);
        }
    }

    return traversal;
}

void Tree::depth_first_expansion_helper(Node* node, std::vector<Node*>& vec) {
    vec.push_back(node);
    for (auto c: node->children) {
        depth_first_expansion_helper(c, vec);
    }
}

std::vector<Node*> Tree::depth_first_expansion() {
    std::vector<Node*> traversal;
    if (root != NULL) {
        depth_first_expansion_helper(root, traversal);
    }
    return traversal;
}

std::string get_newick_string (Tree T, bool print_internal) {
    std::string newick_string = "";

    std::vector<Node*> traversal = T.depth_first_expansion();
    size_t curr_level = 0;
    bool prev_open = true;

    std::stack<std::string> node_stack;

    for (auto n: traversal) {
        size_t level = n->level;
        if (curr_level < level) {
            if (!prev_open) {
                newick_string += ",";
            }
            size_t l = level - 1;
            if (curr_level > 1) {
                l = level - curr_level;
            }
            for (size_t i=0; i < l; i++) {
                newick_string += "(";
                prev_open = true;
            }
            if (n->is_leaf()) {
                newick_string += n->identifier;
                prev_open = false;
            }
            else {
                node_stack.push(n->identifier);
            }
        }
        else if (curr_level > level) {
            prev_open = false;
            for (size_t i = level; i < curr_level; i++) {
                newick_string += ")";
                if (print_internal){
                    newick_string += node_stack.top();
                }
                node_stack.pop();
            }
            if (n->is_leaf()) {
                newick_string += "," + n->identifier;
            }
            else {
                node_stack.push(n->identifier);
            }
        }
        else {
            prev_open = false;
            if (n->is_leaf()) {
                newick_string += "," + n->identifier;
            }
            else {
                node_stack.push(n->identifier);
            }
        }
        curr_level = level;
    }
    for (size_t i = 1; i < curr_level; i++) {
        newick_string += ")";
        if (print_internal) {
            newick_string += node_stack.top();
        }
        node_stack.pop();
    }

    newick_string += ";";
    return newick_string;
}

void split (std::string s, char delim, std::vector<std::string>& words) {
    std::string curr = "";
    for (auto c: s) {
        if (c == delim) {
            words.push_back(curr);
            curr = "";
        }
        else {
            curr += c;
        }
    }
    if (curr != "") {
        words.push_back(curr);
    }
}

void split (std::string s, std::vector<std::string>& words) {
    std::string curr = "";
    std::vector<std::string> ret;
//    for (auto c: s) {
//        if ((c == ' ') || (c == '\t')) {
//            if (curr != "") {
//                ret.push_back(curr);
//            }
//            curr = "";
//        }
//        else {
//            curr += c;
//        }
//    }
//    if (curr != "") {
//        ret.push_back(curr);
//    }
    // Used to split string around spaces.
    std::istringstream ss(s);

    std::string word;
    // Traverse through all words
    while (ss >> word) {
        words.push_back(word);
    };
}

Tree create_tree_from_newick (std::string filename) {
    std::ifstream infile(filename);
    std::string newick_string;
    std::getline(infile, newick_string);
    Tree T;

    std::vector<std::string> leaves;
    std::vector<size_t> num_open;
    std::vector<size_t> num_close;

    std::vector<std::string> s1;
    split(newick_string, ',', s1);

    for (auto s: s1) {
        size_t no = 0;
        size_t nc = 0;
        bool stop = false;
        std::string leaf = "";
        for (auto c: s) {
            if (c == ':') {
                stop = true;
            }
            else if (c == '(') {
                no++;
            }
            else if (c == ')') {
                stop = true;
                nc++;
            }
            else if (!stop) {
                leaf += c;
            }
        }
        leaves.push_back(leaf);
        num_open.push_back(no);
        num_close.push_back(nc);
    }

    if (num_open.size() != num_close.size()) {
        fprintf(stderr, "ERROR: incorrect Newick format!\n");
        exit(1);
    }

    size_t curr_internal_node = 0;
    std::stack<std::string> parent_stack;

    for (size_t i=0; i<leaves.size(); i++) {
        auto leaf = leaves[i];
        auto no = num_open[i];
        auto nc = num_close[i];
        for (size_t j=0; j<no; j++) {
            std::string nid = std::to_string(++curr_internal_node);
            if (parent_stack.size() == 0) {
                T.create_node(nid, nid);
            }
            else {
                T.create_node(nid, nid, parent_stack.top());
            }
            parent_stack.push(nid);
        }
        T.create_node(leaf, leaf, parent_stack.top());
        for (size_t j=0; j<nc; j++) {
            parent_stack.pop();
        }
    }

    return T;
}