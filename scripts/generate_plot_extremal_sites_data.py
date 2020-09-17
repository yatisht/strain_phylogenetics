import sys
import argparse
import numpy as np
import math

def isCtoT(site):
    return any([(s[0]=='C' and s[-1] =='T') for s in site.split(',')])

def isGtoT(site):
    return any([(s[0]=='G' and s[-1] =='T') for s in site.split(',')])


def isRare(n,cutoff):
    return (n <= cutoff)

def isExtremal (s, n, cutoff, min_n_for_s, max_s_for_n):
    return ((n > cutoff) and (n==min_n_for_s[s][0]) and \
            (s==max_s_for_n[n][0]) and \
            all([(s>v[0]) for (k,v) in list(max_s_for_n.items()) if (k<n)]) and \
            all([(n<v[0]) for (k,v) in list(min_n_for_s.items()) if (k>s)]))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Identify extremal sites ')
    parser.add_argument("-in", type=str,
                        help="input parsimony file.")

    args = vars(parser.parse_args())
    parsimony_filename = args.get('in', '')
    if (parsimony_filename== None):
        parser.print_help()
        sys.exit(1)

    max_s_for_n = [dict() for x in range(3)]
    min_n_for_s = [dict() for x in range(3)]
    sites = [dict() for x in range(3)]
    saturated_s = [dict() for x in range(3)]
    all_sites = {}

    for (i,(ignoreCtoT,ignoreGtoT)) in enumerate([(0,0), (1,0), (0,1)]):
        with open(parsimony_filename,'r') as f:
            for line in f:
                if 'alt_alleles' in line:
                    words = line.split()
                    site = words[0]
                    n = 0
                    s = 0
                    for w in words:
                        if 'alt_alleles' in w:
                            n += int(w.split('=')[1])
                    for w in words:
                        if 'parsimony_score' in w:
                            s += int(w.split('=')[1])
                            break
                    s_fwd = 0
                    for w in words:
                        if 'clade_sizes' in w:
                            if (len(w.split('=')) > 1):
                                s_fwd += len([l for l in \
                                              w.split('=')[1].split(',') \
                                              if '[F]' in l])
                    s_bck = s - s_fwd
                    n = int(words[1].split('=')[1])
                    s = int(words[2].split('=')[1])
                    s_fwd = 0
                    all_sites[site] = (n, s, s_fwd, s_bck)
                    if ((i == 0) or ((i == 1) and (not isCtoT(site))) or \
                                     ((i == 2) and (not isGtoT(site)))):
                        sites[i][site] = (n, s, s_fwd, s_bck)
                        mn = min_n_for_s[i].get(s, [1e6, []])
                        ms = max_s_for_n[i].get(n, [0, []])
                        if n == s:
                            saturated_s[i][s] = 1+saturated_s[i].get(s, 0)
                        if n == mn[0]:
                            mn[1].append(site)
                            min_n_for_s[i][s] = mn
                        if s == ms[0]:
                            ms[1].append(site)
                            max_s_for_n[i][n] = ms
                        if n < mn[0]:
                            min_n_for_s[i][s] = [n, [site]]
                        if s > ms[0]:
                            max_s_for_n[i][n] = [s, [site]]
    
    cutoff = [1e9, 1e9, 1e9]
    for i in range(3):
        cutoff[i] = max([k for (k,v) in list(saturated_s[i].items()) if (v>1)])

    for site in list(all_sites.keys()):
        (n,s,s_fwd, s_bck) = all_sites[site]
        to_print = [site, str(n), str(s), str(int(isCtoT(site))), \
                    str(int(isGtoT(site))), str(int(isRare(n, cutoff[0]))),\
                    str(int(isExtremal(s,n,cutoff[0],min_n_for_s[0],max_s_for_n[0])))]
        for i in range(1,3):
            if site in list(sites[i].keys()):
                to_print.append(str(int(isRare(n, cutoff[i]))))
                to_print.append(str(int(isExtremal(s,n,cutoff[i],min_n_for_s[i],max_s_for_n[i]))))
            else:
                to_print.append('?')
                to_print.append('?')
        print('\t'.join(to_print))
        


