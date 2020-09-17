from treelib import Node, Tree
import sys
import argparse
from multiprocessing import Pool
import multiprocessing
import functools
from scipy.stats import entropy
import random

def create_tree(tree_filename):
    tree = Tree()
    # Read Newick file line
    f = open(tree_filename)
    line = f.readline().rstrip()
    f.close()
    # Get leaves names in s2
    s1 = line.split(',')
    s2 = [s.split(':')[0].replace('(', '').split(')')[0] for s in s1]
    stack = [(w.count('('), w.count(')')) for w in s1]
    num_open = sum([s[0] for s in stack])
    num_close = sum([s[1] for s in stack])

    if ((num_open != num_close)):
        print ('ERROR: PhyloTree in incorrect format!')
        sys.exit()

    curr_node = '0'
    curr_idx = 0
    parent_stack = []

    # Construct the tree in treelib format
    for (k, species) in enumerate(s2):
        no = stack[k][0]
        nc = stack[k][1]
        for i in range(no):
            curr_node = str(int(curr_node)+1)
            if len(parent_stack) == 0:
                tree.create_node(curr_node, curr_node)
            else:
                tree.create_node(curr_node, curr_node, parent=parent_stack[-1])
            parent_stack.append(curr_node)
        tree.create_node(species, species, parent=parent_stack[-1])
        for i in range(nc):
            nid = parent_stack[-1]
            curr_idx += 1
            parent_stack.pop()

    return tree

# Convert tree back to Newick string
def get_newick_string (tree):
    dfs = [n for n in tree.expand_tree(mode=1, sorting=False, reverse=True)]
    depths = [tree.level(n) for n in dfs]
    is_leaf = [tree.get_node(n).is_leaf() for n in dfs]
    curr_depth = -1
    newick_string = ''
    prev_open = True 
    stack = []
    for (k,s) in enumerate(dfs):
        tag = tree.get_node(s).tag
        depth = depths[k]
        if (curr_depth < depth):
            if (not prev_open):
                newick_string += ','
            for i in range(depth-max(0,curr_depth)):
                newick_string += '('
                prev_open = True
            if (is_leaf[k]):
                newick_string += tag
                prev_open = False
            else:
                stack.append(tag)
        elif (curr_depth > depth):
            prev_open = False
            for i in range(curr_depth-depth):
                newick_string += ')'
                #newick_string += stack[-1]
                stack.pop()
            if (is_leaf[k]):
                newick_string += ','+tag
            else:
                stack.append(tag)
        else:
            prev_open = False
            if (is_leaf[k]):
                newick_string += ','+tag
            else:
                stack.append(tag)
        curr_depth = depth
    for i in range(curr_depth):
        newick_string += ')'
        #newick_string += stack[-1]
        stack.pop()
    newick_string += ';'
    return newick_string

def get_all_node_ids (tree):
    return list(tree.nodes.keys())

def get_internal_node_ids (tree):
    all_node_ids = set(get_all_node_ids(tree))
    leaf_ids = set([n.identifier for n in tree.leaves()])
    return list(all_node_ids - leaf_ids)

def get_leaf_node_ids (tree):
    leaf_ids = [n.identifier for n in tree.leaves()]
    return leaf_ids

def get_leaf_node_ids_for_node (tree, nid):
    leaf_ids = [n.identifier for n in tree.leaves(nid)]
    return leaf_ids

def get_branch_split (tree, bid):
    all_leaf_ids = set([n.identifier for n in tree.leaves()])
    A = set([n.identifier for n in tree.leaves(bid)])
    B = (all_leaf_ids - A)
    return (list(A), list(B), bid)

def get_symmetric_difference (A, X, common_leaf_ids):
    d1 = len((set(A)-set(X)).intersection(common_leaf_ids)) 
    d2 = len((set(X)-set(A)).intersection(common_leaf_ids)) 
    return (d1+d2)

# returns None if clade not compatible
def where_to_insert_compatible_clade (T, C):
    curr_nodes = [T.root]
    curr_parent = '-1'
    while (len(curr_nodes) > 0):
        new_nodes = []
        chosen_n = None
        children_to_merge = []
        for n in curr_nodes:
            n_leaves = set(get_leaf_node_ids_for_node(T, n)) 
            if C.issubset(n_leaves):
                if len(C.symmetric_difference(n_leaves)) == 0:
                    # clade already present
                    return [n, None]
                chosen_n = n
                break
            elif n_leaves.issubset(C):
                children_to_merge.append(n)
            elif len(C.intersection(n_leaves)) > 0:
                return None
        if (chosen_n == None):
            return [curr_parent, children_to_merge]
        # C is contained in chosen_n
        else:
            #should we return fromm here?
            curr_parent = chosen_n
            for c in T.children(chosen_n):
                new_nodes.append(c.identifier)
        curr_nodes = new_nodes 
    return None

def get_unique_id(T):
    internal_nids = [int(i) for i in get_internal_node_ids(T)]
    return str(1+max(internal_nids))

def insert_new_compatible_clade(T, par, children_to_merge):
    if (par != None) and (len(children_to_merge) > 0):
        uid = get_unique_id(T)
        T.create_node(uid, uid, parent=par)
        for c in children_to_merge:
            T.move_node(c, uid)

def get_intersection_tree(T1, T2):
    T = Tree(tree=T1, deep=True)
    T1_bfs = [n for n in T1.expand_tree(mode=1)]
    T2_bfs = [n for n in T2.expand_tree(mode=1)]
    for nid in T1_bfs:
        X = set(get_leaf_node_ids_for_node(T,nid))
        diff = min([len(X.symmetric_difference(set( \
                    get_leaf_node_ids_for_node(T2,i)))) \
                    for i in T2_bfs])
        if diff != 0:
            par = T.parent(nid).identifier
            for c in T.children(nid):
                T.move_node(c.identifier, par)
            T.remove_subtree(nid)
    return T

def merge_trees (T1, T2):
    T = get_intersection_tree(T1,T2)
    T1_bfs = [n for n in T1.expand_tree(mode=1)]
    T2_bfs = [n for n in T2.expand_tree(mode=1)]
    for nid in T1_bfs:
        X = set(get_leaf_node_ids_for_node(T1,nid))
        ret = where_to_insert_compatible_clade(T, X)
        if ret != None:
            if ret[1] != None:
                insert_new_compatible_clade(T,ret[0],ret[1])
    for nid in T2_bfs:
        X = set(get_leaf_node_ids_for_node(T2,nid))
        ret = where_to_insert_compatible_clade(T, X)
        if ret != None:
            if ret[1] != None:
                insert_new_compatible_clade(T,ret[0],ret[1])
    return T

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Combine the clades from two '
                                     'trees to find a greedy-maximal set of '
                                     'pairwise-compatible clades that form a '
                                     '"merged" tree.')
    parser.add_argument("-T1", type=str,
                        help="tree 1 (in Newick format)")
    parser.add_argument("-T2", type=str,
                        help="tree 2 (in Newick format)")
    parser.add_argument("-T_out", type=str,
                        help="output tree filename (in Newick format)")
    parser.add_argument("-intersectOnly", type=str,
                        help="output intersection (instead of a maximal merge) of "
                        "T1 and T2 [OPTIONAL, DEFAULT=0].")
    parser.add_argument("-symmetric", type=str,
                        help="output symmetric merge of "
                        "T1 and T2 [OPTIONAL, DEFAULT=0].")

    args = vars(parser.parse_args())
    T1_filename = args.get('T1', '')
    if (T1_filename == None):
        parser.print_help()
        sys.exit(1)
    T2_filename = args.get('T2', '')
    if (T2_filename == None):
        parser.print_help()
        sys.exit(1)
    T_outfilename = args.get('T_out', '')
    if (T_outfilename == None):
        parser.print_help()
        sys.exit(1)
    intersectOnly = args.get('intersectOnly', '')
    if (intersectOnly == None):
        intersectOnly == '0'
    symmetric = args.get('symmetric', '')
    if (symmetric == None):
        symmetric == '0'

    if (symmetric == '1') and (intersectOnly == '1'):
        print('ERROR: intersectOnly and symmetric options cannot be used together.')
        sys.exit(1)
    
    T1 = create_tree(T1_filename)
    T2 = create_tree(T2_filename)

    T1_leaves = set(get_leaf_node_ids(T1))
    T2_leaves = set(get_leaf_node_ids(T2))

    if len(T1_leaves.symmetric_difference(T2_leaves)) != 0:
        print('ERROR: T1 and T2 should have the same leaf set!')
        sys.exit(1)

    if intersectOnly == '1':
        T = get_intersection_tree(T1, T2)
    elif symmetric == '1':
        merge_T1_T2 = merge_trees(T1, T2)
        merge_T2_T1 = merge_trees(T2, T1)
        T = get_intersection_tree(merge_T1_T2, merge_T2_T1)
    else:
        T = merge_trees(T1, T2)
    
    T_newick = get_newick_string(T)
    outfile = open(T_outfilename, 'w')
    print(T_newick, file=outfile)
    outfile.close()
    

