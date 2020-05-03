from treelib import Node, Tree
import sys
import argparse
from multiprocessing import Pool
import multiprocessing
import functools
import matplotlib.pyplot as plt

def create_tree(tree_filename):
    tree = Tree()
    f = open(tree_filename)
    lines = f.readlines()
    first_line = lines[0]
    l = first_line.rstrip()
    s1 = l.split(',')
    s2 = [s.split(':')[0].replace('(', '').replace(')', '') for s in s1]
    s3 = l.split('#')
    s4 = [s.split(')')[1] for s in s3 if (len(s.split(')')) > 1)]
    s5 = []
    start = False
    word = ''
    for c in l:
        if c == '#':
            start = True
        elif ((start) and (c != ':')):
            word += c
        elif ((start) and (c == ':')):
            start = False
            s5.append(word)
            word = ''
    stack = [(w.count('('), w.count(')')) for w in s1]
    num_open = sum([s[0] for s in stack])
    num_close = sum([s[1] for s in stack])

    if ((num_open != num_close)):
        print ('ERROR: PhyloTree in incorrect format!')
        sys.exit()

    curr_node = '0'
    curr_idx = 0
    parent_stack = []

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
            tree.update_node(nid, tag=(nid+':'+s4[curr_idx]+':'+s5[curr_idx]))
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

def get_branch_dichotomy (tree, bid):
    all_leaf_ids = set([n.identifier for n in tree.leaves()])
    A = set([n.identifier for n in tree.leaves(bid)])
    B = (all_leaf_ids - A)
    return (list(A), list(B), bid)

def get_symmetric_difference (A, X, common_leaf_ids):
    d1 = len((set(A)-set(X)).intersection(common_leaf_ids)) 
    d2 = len((set(X)-set(A)).intersection(common_leaf_ids)) 
    return (d1+d2)

def get_directed_change(bid, tree1, tree2_dichotomies, common_leaf_ids):
    (A, B, bid) = get_branch_dichotomy(tree1, bid)
    return (min([(get_symmetric_difference(A,X,common_leaf_ids), bid2) \
                 for (X,Y,bid2) in tree2_dichotomies]))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Compute directed change '
                                     ' for each branch of T1 wrt tree T2.')
    parser.add_argument("-T1", type=str,
                        help="tree 1 (in Newick format)")
    parser.add_argument("-T2", type=str,
                        help="tree 2 (in Newick format)")

    if len(sys.argv) <= 4:
        parser.print_help()
        sys.exit(1)

    args = vars(parser.parse_args())
    T1_filename = args['T1']
    T2_filename = args['T2']
    T1 = create_tree(T1_filename)
    T2 = create_tree(T2_filename)

    print ('T1: ')
    T1.show()

    print ('T2: ')
    T2.show()

    T1_internal_node_ids = get_internal_node_ids(T1)
    T2_internal_node_ids = get_internal_node_ids(T2)
    
    T1_leaf_node_ids = get_leaf_node_ids(T1)
    T2_leaf_node_ids = get_leaf_node_ids(T2)
    common_leaf_ids = set(T1_leaf_node_ids).intersection(set(T2_leaf_node_ids))
    
    T2_dichotomies = [get_branch_dichotomy(T2, bid) for bid in \
                         T2_internal_node_ids] 

    pool = Pool(processes=multiprocessing.cpu_count())
    
    worker = functools.partial(get_directed_change, tree1=T1, \
        tree2_dichotomies=T2_dichotomies, common_leaf_ids = common_leaf_ids) 

    # directed change for each branch in T1
    C = pool.map(worker, T1_internal_node_ids)

    for (k, bid) in enumerate(T1_internal_node_ids):
        tag1 = T1.get_node(bid).tag
        tag2 = T2.get_node(C[k][1]).tag 
        T1.update_node(bid, tag=(tag1+" change="+str(C[k][0])+\
                                        " best-match="+tag2))

    print ('T1 (with directed changes): ')
    T1.show()

    print ('T1 size: '+str(len(T1_internal_node_ids)))
    print ('T2 size: '+str(len(T2_internal_node_ids)))

    clade_changes = {}

    for lid in common_leaf_ids:
        if T1.level(lid) != 0:
            parent_tag = T1.parent(lid).tag
            T1_clade = parent_tag.split()[0].split(':')[1]
            T2_clade = parent_tag.split()[2].split('=')[1].split(':')[1]
            if T1_clade != T2_clade:
                if (T1_clade != 'unassigned') and (T2_clade != 'unassigned'):
                    key = T1_clade + "<-" + T2_clade
                elif (T1_clade == 'unassigned'):
                    key =  "unassigned<-*" 
                else:
                    key =  "*<-unassigned" 
                n = clade_changes.get(key, 0)
                clade_changes[key] = n+1

    for (k,v) in clade_changes.items():
        print (k+": "+str(v))

    print (len(C))


