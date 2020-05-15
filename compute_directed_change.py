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
                newick_string += stack[-1]
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
        newick_string += stack[-1]
        stack.pop()
    newick_string += ';'
    return newick_string

def create_tree(tree_filename):
    tree = Tree()
    clade_roots = {}
    f = open(tree_filename)
    lines = f.readlines()
    first_line = lines[0]
    l = first_line.rstrip()
    s1 = l.split(',')
    s2 = [s.split(':')[0].replace('(', '').replace(')', '') for s in s1]
    s3 = [s.split('|')[1] if (len(s.split('|')) > 1) else s for s in s2]
    s4 = [s.replace("'", "") for s in s3]
    stack = [(w.count('('), w.count(')')) for w in s1]
    num_open = sum([s[0] for s in stack])
    num_close = sum([s[1] for s in stack])

    if ((num_open != num_close)):
        print ('ERROR: PhyloTree in incorrect format!')
        sys.exit()

    curr_node = '0'
    curr_idx = 0
    parent_stack = []

    for (k, species) in enumerate(s4):
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
            curr_idx += 1
            parent_stack.pop()
    return tree

def get_all_node_ids (tree):
    return tree.nodes.keys()

def get_internal_node_ids (tree):
    all_node_ids = set(get_all_node_ids(tree))
    leaf_ids = set([n.identifier for n in tree.leaves()])
    return list(all_node_ids - leaf_ids)

def get_leaf_node_ids (tree):
    leaf_ids = [n.identifier for n in tree.leaves()]
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

def get_split_distance(bid, tree1, tree2_splits, common_leaf_ids):
    (A, B, bid) = get_branch_split(tree1, bid)
    dist = [(get_symmetric_difference(A,X,common_leaf_ids), bid2) \
               for (X,Y,bid2) in tree2_splits] 
    min_dist = min(dist)[0]
    ret = [(0, '-1')]
    if (min_dist != len(tree1.leaves(bid))):
        ret = [x for x in dist if (x[0] == min_dist)]
    return ret

def get_entropy_weighted_split_distance(bid, tree1, tree2, tree2_splits, \
                                         common_leaf_ids):
    if (len(common_leaf_ids) == 0):
        return [(0, 0, '-1')]
    (A, B, bid) = get_branch_split(tree1, bid)
    s = len(set(A).intersection(common_leaf_ids))
    p = (1.0*s)/len(common_leaf_ids)
    q = 1.0 - p
    ent = entropy([p, q], base=2)
    dist = []
    min_dist = 1e9
    if (p == 0):
        return [(0, 0, '-1')]
    for (X,Y,bid2) in tree2_splits:
        diff = get_symmetric_difference(A,X,common_leaf_ids)
        dist.append((diff, ent*diff, bid2))
        if (diff < min_dist):
            min_dist = diff
    if min_dist > s: 
        return [(s, ent*s, '-1')]
    m = [x for x in dist if (x[0] == min_dist)]
    ret = []
    for i in range(len(m)):
        found_anc = False
        for j in range(len(m)):
            if (i != j) and (tree2.is_ancestor(m[j][2], m[i][2])):
                found_anc = True
                break
        if not found_anc:
            ret.append(m[i])
    return ret

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Compute directed change '
                                     ' for each branch of T1 wrt tree T2.')
    parser.add_argument("-T1", type=str,
                        help="tree 1 (in Newick format)")
    parser.add_argument("-T2", type=str,
                        help="tree 2 (in Newick format)")
    parser.add_argument("-CORES", type=str,
                        help="number of CPU cores (default=1) to use [OPTIONAL]")

    if len(sys.argv) <= 4:
        parser.print_help()
        sys.exit(1)

    args = vars(parser.parse_args())
    T1_filename = args.get('T1', '')
    if (T1_filename == None):
        parser.print_help()
        sys.exit(1)
    T2_filename = args.get('T2', '')
    if (T2_filename == None):
        parser.print_help()
        sys.exit(1)
    num_cores = args.get('CORES', '')
    if (num_cores == None):
        num_cores = '1'
    
    T1 = create_tree(T1_filename)
    T2 = create_tree(T2_filename)

    T1_newick = get_newick_string(T1)
    T2_newick = get_newick_string(T2)
    
    print 'T1 (in newick format) with internal nodes labelled: '
    print T1_newick 
    
    print 'T2 (in newick format) with internal nodes labelled: '
    print T2_newick 


    T1_internal_node_ids = get_internal_node_ids(T1)
    T2_internal_node_ids = get_internal_node_ids(T2)
    
    T1_leaf_node_ids = get_leaf_node_ids(T1)
    T2_leaf_node_ids = get_leaf_node_ids(T2)
    common_leaf_ids = set(T1_leaf_node_ids).intersection(set(T2_leaf_node_ids))
    
    T2_splits = [get_branch_split(T2, bid) for bid in \
                         T2_internal_node_ids] 

    pool = Pool(processes=int(num_cores))
    
    worker = functools.partial(get_entropy_weighted_split_distance, tree1=T1, \
                               tree2=T2, tree2_splits=T2_splits, \
                               common_leaf_ids = common_leaf_ids) 

    # directed change for each branch in T1
    C = pool.map(worker, T1_internal_node_ids)

    total_dist = 0.0
    print '#T1_node\tT2_best_matches\tmatching_split_distance'+\
            '\tentropy_weighted_split_distance'
    for (k, bid) in enumerate(T1_internal_node_ids):
        tag1 = T1.get_node(bid).tag
        T2_tags_arr = [T2.get_node(m[2]).tag if m[2] != '-1' else '-1' \
                       for m in C[k]] 
        T2_tags = ','.join(T2_tags_arr)
        print tag1+'\t'+T2_tags+'\t'+str(C[k][0][0])+'\t'+str(C[k][0][1])
        total_dist += C[k][0][1]

    print 'T1 size: ', len(T1_leaf_node_ids)
    print 'T2 size: ', len(T2_leaf_node_ids)
    print 'Total common leaves: ', len(list(common_leaf_ids))
    print 'Total entropy-weighted split distance: ', total_dist

    
