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

def get_minimum_subclades (S, C, tree, curr_min):
    marked = [0 for i in range(len(S))]
    index = {}
    for k,v in enumerate(S):
        index[v] = k
    subclades = 0
    for (i, s) in enumerate(S):
        if marked[i] == 0:
            uppermost = s
            marked[i] = 1
            subclades += 1
            if subclades > curr_min:
                return 1e9
            ancestors = tree.rsearch(s)
            for anc in ancestors:
                L = set([l.identifier for l in tree.leaves(anc)])
                extra = (L-S).intersection(C)
                if len(extra) == 0:
                    for m in L.intersection(S):
                        marked[index[m]] = 1
                else:
                    break
    return subclades

def get_subclade_difference(A, X, common_leaf_ids, tree1, tree2, min_dist):
    S1 = (set(A)-set(X)).intersection(common_leaf_ids)
    d1 = get_minimum_subclades (S1, common_leaf_ids, tree2, min_dist)
    S2 = (set(X)-set(A)).intersection(common_leaf_ids)
    d2 = get_minimum_subclades (S2, common_leaf_ids, tree1, min_dist)
    return d1+d2

def get_entropy_weighted_split_distance(bid, tree1, tree2, tree2_splits, \
                                         common_leaf_ids, Z):
    if (len(common_leaf_ids) == 0):
        return [(0, 0, '-1')]
    (A, B, bid) = get_branch_split(tree1, bid)
    s = len(set(A).intersection(common_leaf_ids))
    p = (1.0*s)/len(common_leaf_ids)
    q = 1.0 - p
    ent = entropy([p, q], base=2)
    #ent = 2*min([p, q])
    dist = []
    min_dist = 1e9
    min_dist_splits = []
    if (p == 0):
        return [(0, 0, '-1')]
    for (X,Y,bid2) in tree2_splits:
        diff = get_symmetric_difference(A,X,common_leaf_ids)
        dist.append((diff, ent*diff/Z, bid2))
        if (diff < min_dist):
            min_dist = diff
    if min_dist >= s: 
        return [(s, ent*s/Z, '-1')]
    ret = [x for x in dist if (x[0] == min_dist)]
    return ret

def calculate_n_divergence(all_dist, n):
    return sum(all_dist[0:n])

def get_total_entropy(tree1, tree2):
    tree1_internal_nids = get_internal_node_ids(tree1)
    tree1_leaf_node_ids = get_leaf_node_ids(tree1)
    tree2_leaf_node_ids = get_leaf_node_ids(tree2)
    common_leaf_ids = set(tree1_leaf_node_ids).intersection(set(tree2_leaf_node_ids))
    E = 0
    for bid in tree1_internal_nids:
        (A, B, bid) = get_branch_split(tree1, bid)
        s = len(set(A).intersection(common_leaf_ids))
        p = (1.0*s)/len(common_leaf_ids)
        q = 1.0 - p
        ent = entropy([p, q], base=2)
        E += ent
    return E

def get_entropy_weighted_total_distance(tree1, tree2, do_print):
    tree1_internal_node_ids = get_internal_node_ids(tree1)
    tree2_internal_node_ids = get_internal_node_ids(tree2)
    
    tree1_leaf_node_ids = get_leaf_node_ids(tree1)
    tree2_leaf_node_ids = get_leaf_node_ids(tree2)
    common_leaf_ids = set(tree1_leaf_node_ids).intersection(set(tree2_leaf_node_ids))
    
    tree2_splits = [get_branch_split(tree2, bid) for bid in \
                         tree2_internal_node_ids] 

    Z = get_total_entropy(tree1, tree2)

    pool = Pool(processes=int(num_cores))
    
    worker = functools.partial(get_entropy_weighted_split_distance, tree1=tree1, \
                               tree2=tree2, tree2_splits=tree2_splits, \
                               common_leaf_ids = common_leaf_ids, Z = Z) 

    C = pool.map(worker, tree1_internal_node_ids)

    total_dist = 0.0
    all_dist = []

    to_print = []
    for (k, bid) in enumerate(tree1_internal_node_ids):
        tag1 = tree1.get_node(bid).tag
        tree2_tags_arr = [tree2.get_node(m[2]).tag if m[2] != '-1' else '-1' \
                       for m in C[k]] 
        tree2_tags = ','.join(tree2_tags_arr)
        if (do_print):
            to_print.append(tag1+'\t'+tree2_tags+'\t'+str(C[k][0][0])+'\t'\
                            +str(C[k][0][1]))
        total_dist += C[k][0][1]
        all_dist.append(C[k][0][1])

    all_dist.sort(reverse=True)

    return total_dist, to_print

def permute_leaves(T):
    tree = Tree(tree=T, deep=True)
    leaves = get_leaf_node_ids(tree)
    shuf_leaves = leaves[:] 
    random.shuffle(shuf_leaves)
    for k,nid in enumerate(leaves):
        new_nid = shuf_leaves[k]
        tree.update_node(nid, tag=new_nid, identifier='L'+str(k))
    for k,nid in enumerate(leaves):
        new_nid = shuf_leaves[k]
        tree.update_node('L'+str(k), identifier=new_nid)
    return tree

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Compute directed change '
                                     ' for each branch of T1 wrt tree T2.')
    parser.add_argument("-T1", type=str,
                        help="tree 1 (in Newick format)")
    parser.add_argument("-T2", type=str,
                        help="tree 2 (in Newick format)")
    parser.add_argument("-permutedNorm", type=str,
                        help="use permuted trees (default=0)")
    parser.add_argument("-CORES", type=str,
                        help="number of CPU cores (default=1) to use [OPTIONAL]")

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
    permuted_norm = args.get('permutedNorm', '')
    if (permuted_norm == None):
        permuted_norm = '0'
    
    T1 = create_tree(T1_filename)
    T2 = create_tree(T2_filename)
    
    if (permuted_norm == '1'):
        T1_p = permute_leaves(T1)
        T2_p = permute_leaves(T2)

    T1_newick = get_newick_string(T1)
    T2_newick = get_newick_string(T2)
    print 'T1 (in newick format) with internal nodes labelled: '
    print T1_newick 
    print 'T2 (in newick format) with internal nodes labelled: '
    print T2_newick 

    dist_t1_t2, to_print = get_entropy_weighted_total_distance(T1, T2, True) 
    print '#T1_node\tT2_best_matches\tmatching_split_distance\tentropy_weighted_split_distance'
    print '\n'.join(to_print)
                        
    dist_t2_t1, to_print = get_entropy_weighted_total_distance(T2, T1, True) 
    print '#T2_node\tT1_best_matches\tmatching_split_distance\tentropy_weighted_split_distance'
    print '\n'.join(to_print)

    print 'D(T1,T2) = ', dist_t1_t2 
    print 'D(T2,T1) = ', dist_t2_t1 
    print 'S(T1,T2) = ', (dist_t1_t2+dist_t2_t1)/2

    if (permuted_norm == '1'):
        dist_t1_t2_p, to_print = get_entropy_weighted_total_distance(T1, T2_p, False) 
        dist_t2_t1_p, to_print = get_entropy_weighted_total_distance(T2, T1_p, False) 
        print 'D(T1,T2_p) = ', dist_t1_t2_p 
        print 'D(T2,T1_p) = ', dist_t2_t1_p 
        print 'S_p(T1,T2) = ', (dist_t1_t2+dist_t2_t1)/(dist_t1_t2_p+dist_t2_t1_p)
