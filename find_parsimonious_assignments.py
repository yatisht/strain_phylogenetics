from treelib import Node, Tree
import sys
import io
import argparse
import gzip

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

# Sankoff algorithm using a simple scoring to make most parsimonious
# assignments to tree nodes and return transitions
def get_most_parsimonious_transitions (tree, bfs, bfs_idx, var_ids):
    s = [(0, 1e6) for i in range(len(bfs_idx.keys()))]
    states = [0 for i in range(len(bfs_idx.keys()))]
    forward_mutations = []
    back_mutations = []
    leaves_affected_forward = []
    leaves_affected_backward = []
    for nid in var_ids:
        idx = bfs_idx[nid]
        s[idx] = (1e6, 0)
    #forward pass
    for nid in bfs[::-1]:
        node = tree.get_node(nid)
        node_idx = bfs_idx[nid]
        if (not node.is_leaf()):
            s_ref = 0 
            s_alt = 0 
            children = tree.children(nid)
            for child in children:
                c_id = child.identifier
                c_idx = bfs_idx[c_id]
                s_ref += min([s[c_idx][0], 1+s[c_idx][1]])
                s_alt += min([1+s[c_idx][0], s[c_idx][1]])
            s[node_idx] = (s_ref, s_alt)
    flagged_leaves = []
    #backward pass
    for nid in bfs:
        node_idx = bfs_idx[nid]
        (s_ref, s_alt) = s[node_idx]
        state = 0
        par_state = 0
        if (nid != tree.root):
            par = tree.parent(nid)
            par_id = par.identifier
            par_idx = bfs_idx[par_id]
            par_state = states[par_idx]
        if (s_ref < s_alt):
            state = 0
        elif (s_ref == s_alt):
            state = par_state
        else:
            state = 1
        states[node_idx] = state
        if (state != par_state):
            leaves = tree.leaves(nid)
            if (state == 1):
                forward_mutations.append(nid)
            else:
                back_mutations.append(nid)
    for nid in forward_mutations:
        leaves = tree.leaves(nid)
        l_f = set([l for l in leaves if states[bfs_idx[l.identifier]] == 1])
        leaves_affected_forward.append(l_f)
        if (len(leaves) < 4):
            for l in leaves:
                if (states[bfs_idx[l.identifier]] == 1):
                    flagged_leaves.append(l.identifier)
    for nid in back_mutations:
        leaves = tree.leaves(nid)
        l_b = set([l for l in leaves if states[bfs_idx[l.identifier]] == 0])
        leaves_affected_backward.append(l_b)
        if (len(leaves) < 4):
            for l in leaves:
                if (states[bfs_idx[l.identifier]] == 0):
                    flagged_leaves.append(l.identifier)
    for i in range(len(leaves_affected_forward)):
        for j in range(i):
            leaves_affected_forward[j] = leaves_affected_forward[j] - \
            leaves_affected_forward[i]
    for i in range(len(leaves_affected_backward)):
        for j in range(i):
            leaves_affected_backward[j] = leaves_affected_backward[j] - \
            leaves_affected_backward[i]
    return forward_mutations, back_mutations,\
            [str(len(l)) for l in leaves_affected_forward], \
            [str(len(l)) for l in leaves_affected_backward], \
            flagged_leaves


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Compute most parsimonious '
                                     'character assignments for each internal '
                                     'tree node given the set of variants for '
                                     'its leaves.')
    parser.add_argument("-tree", type=str,
                        help="input tree (in Newick format)")
    parser.add_argument("-vcf", type=str,
                        help="vcf file corresponding to the tree")
    parser.add_argument("-variants", type=str,
                        help="restrict to a list of variants (':' "
                        "separated)  from VCF file [OPTIONAL]")

    args = vars(parser.parse_args())

    tree_filename = args.get('tree', '')
    if (tree_filename == None):
        parser.print_help()
        sys.exit(1)
    vcf_filename = args.get('vcf', '')
    if (vcf_filename == None):
        parser.print_help()
        sys.exit(1)
    variants_args = args.get('variants', '')

    variant_list = None if (variants_args == None) else variants_args.split(':')

    tree = create_tree(tree_filename)

    labelled_tree_newick = get_newick_string(tree)
    print 'Tree (in newick format) with internal nodes labelled: '
    print labelled_tree_newick

    bfs = [n for n in tree.expand_tree(mode=2, sorting=False, reverse=True)]
    bfs_idx = {}
    for (k,v) in enumerate(bfs):
        bfs_idx[v] = k

    header_found = False
    header_line = ''
    vcf_ids = []

    total_variants = 0
    total_parsimony_score = 0

    # Check vcf file type
    if vcf_filename.endswith(".gz"):
        vcf_file = io.TextIOWrapper(io.BufferedReader(gzip.open(vcf_filename)))
    else:
        vcf_file = file(vcf_filename,'r')

    with vcf_file as f:
        for line in f:
            if (not header_found):
                if ('REF' in line):
                    header_found = True
                    words = line.split()
                    vcf_ids = words[9:]
#                    vcf_ids = [w.replace('/', '_').replace('|', '_') for w in \
#                               words[9:]]
            else:
                total_variants += 1
                print_variant = False
                words = line.split()
                variant_pos = int(words[1])
                variant = words[2]
                variant_ref = words[3]
                variant_ids = [vcf_ids[k] for (k, w) in enumerate(words[9:]) \
                               if (w.split(':')[0] == '1')]
                if ((variant_list == None) or (variant in variant_list)):
                    f_mut, b_mut, l_f, l_b, flagged_leaves = \
                            get_most_parsimonious_transitions \
                            (tree, bfs, bfs_idx, variant_ids)
                    print variant+'\talt_alleles='+str(len(variant_ids))+\
                          '\tparsimony_score='+ str(len(f_mut)+len(b_mut))+\
                          '\tforward_mutation_nodes='+','.join(f_mut)+\
                          '\tback_mutation_nodes='+','.join(b_mut)+\
                          '\tforward_mutation_clade_sizes='+','.join(l_f)+\
                          '\tback_mutation_clade_sizes='+','.join(l_b)+\
                          '\tflagged_leaves='+','.join(flagged_leaves)  

                    total_parsimony_score += len(f_mut)+len(b_mut)
                
    print 'Total leaf nodes: ', len(tree.leaves())
    print 'Total variants: ', total_variants
    print 'Total parsimony score: ', total_parsimony_score

