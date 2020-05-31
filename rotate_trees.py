from treelib import Node, Tree
import sys
from math import log
import argparse

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
    dfs = [n for n in tree.expand_tree(mode=1, sorting=False, reverse=False)]
    depths = [tree.level(n) for n in dfs]
    is_leaf = [tree.get_node(n).is_leaf() for n in dfs]
    curr_depth = -1
    newick_string = ''
    prev_open = True 
    stack = []
    leaves_ordered = []
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
                leaves_ordered.append(tag)
                prev_open = False
            else:
                stack.append(tag)
        elif (curr_depth > depth):
            prev_open = False
            for i in range(curr_depth-depth):
                newick_string += ')'
                stack.pop()
            if (is_leaf[k]):
                newick_string += ','+tag
                leaves_ordered.append(tag)
            else:
                stack.append(tag)
        else:
            prev_open = False
            if (is_leaf[k]):
                newick_string += ','+tag
                leaves_ordered.append(tag)
            else:
                stack.append(tag)
        curr_depth = depth
    for i in range(curr_depth):
        newick_string += ')'
        stack.pop()
    newick_string += ';'
    return (leaves_ordered, newick_string)

def Jaccard_similarity(A, B):
    num = len(set(A).intersection(set(B)))
    den = len(set(A).union(set(B)))
    return (1.0*num)/den

def get_avg_rank (S, rank_dict):
    #num = sum([log(1+rank_dict.get(s,0)) for s in S])
    num = sum([rank_dict.get(s,0) for s in S])
    den = 1+sum([1 if s in rank_dict.keys() else 0 for s in S])
    return (1.0*num)/den #if den==0 else 1e6

def rotate_trees(T1_leaves_order, T2):
    T1_order_dict = {}
    tree = Tree()
    for (k, v) in enumerate(T1_leaves_order):
        T1_order_dict[v] = 1.0*k
    curr_nid = '1'
    tree.create_node(curr_nid, curr_nid)
    curr_nodes = ['1']
    num_rot = 0
    while (len(curr_nodes) > 0):
        new_nodes = []
        for n in curr_nodes:
            curr_children = []
            for c in T2.children(n):
                curr_children.append(c.identifier)
            s = [(get_avg_rank([l.identifier for l in T2.leaves(c)], \
                               T1_order_dict), c) for c in \
                 curr_children]
            sorted_s = sorted(s, key=lambda x: x[0])
            if sorted_s != s:
                num_rot += 1
            for i in sorted_s:
                tree.create_node(i[1], i[1], parent=n)
                new_nodes.append(i[1])
        curr_nodes = new_nodes 
    return tree, num_rot

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Rotate trees for better '
                                     ' visualization')
    parser.add_argument("-T1", type=str,
                        help="input tree 1 (in Newick format)")
    parser.add_argument("-T2", type=str,
                        help="input tree 2 (in Newick format)")
    parser.add_argument("-T1_out", type=str,
                        help="output rotated tree 1 (in Newick format)")
    parser.add_argument("-T2_out", type=str,
                        help="output rotated tree 2 (in Newick format)")

    args = vars(parser.parse_args())

    T1_filename = args.get('T1', '')
    if (T1_filename == None):
        parser.print_help()
        sys.exit(1)
    T2_filename = args.get('T2', '')
    if (T2_filename == None):
        parser.print_help()
        sys.exit(1)
    T1_outfilename = args.get('T1_out', '')
    if (T1_outfilename == None):
        parser.print_help()
        sys.exit(1)
    T2_outfilename = args.get('T2_out', '')
    if (T2_outfilename == None):
        parser.print_help()
        sys.exit(1)


    T1 = create_tree(T1_filename)
    T2 = create_tree(T2_filename)

    n1=1
    n2=1
    while (n1+n2 > 0):
        (T1_leaves_ordered, T1_newick) = get_newick_string(T1)
        T2,n1 = rotate_trees(T1_leaves_ordered, T2)
        (T2_leaves_ordered, T2_newick) = get_newick_string(T2)
        T1,n2 = rotate_trees(T2_leaves_ordered, T1)
    
    f1 = open(T1_outfilename, 'w')
    f2 = open(T2_outfilename, 'w')
    print>>f1, get_newick_string(T1)[1]
    print>>f2, get_newick_string(T2)[1]
