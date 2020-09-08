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
            all([(s>v[0]) for (k,v) in max_s_for_n.items() if (k<n)]) and \
            all([(n<v[0]) for (k,v) in min_n_for_s.items() if (k>s)]))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Identify extremal sites ')
    parser.add_argument("-in", type=str,
                        help="input parsimony file.")
    parser.add_argument("-ignoreCtoT", type=int,
                        help="set to 1 to ignore C>T sites (default=0).")
    parser.add_argument("-ignoreGtoT", type=int,
                        help="set to 1 to ignore G>T sites (default=0).")

    args = vars(parser.parse_args())
    parsimony_filename = args.get('in', '')
    if (parsimony_filename== None):
        parser.print_help()
        sys.exit(1)

    ignoreCtoT = args.get('ignoreCtoT', '0')
    if (ignoreCtoT == None):
        ignoreCtoT = '0'
    ignoreCtoT = int(ignoreCtoT)
    
    ignoreGtoT = args.get('ignoreGtoT', '0')
    if (ignoreGtoT == None):
        ignoreGtoT = '0'
    ignoreGtoT = int(ignoreGtoT)

    max_s_for_n = {}
    min_n_for_s = {}
    sites = {}
    x = []
    y = []
    saturated_s = {}
    with file(parsimony_filename,'r') as f:
        for line in f:
            if 'parsimony_score' in line:
                words = line.split()
                site = words[0]
                if ((ignoreCtoT) and (isCtoT(site)) or \
                    (ignoreGtoT) and (isGtoT(site))):
                    continue
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
                sites[site] = (n, s, s_fwd, s_bck)
                mn = min_n_for_s.get(s, [1e6, []])
                ms = max_s_for_n.get(n, [0, []])
                if n == s:
                    saturated_s[s] = 1+saturated_s.get(s, 0)
                if n == mn[0]:
                    mn[1].append(site)
                    min_n_for_s[s] = mn
                if s == ms[0]:
                    ms[1].append(site)
                    max_s_for_n[n] = ms
                if n < mn[0]:
                    min_n_for_s[s] = [n, [site]]
                if s > ms[0]:
                    max_s_for_n[n] = [s, [site]]
    
    cutoff = max([k for (k,v) in saturated_s.items() if (v>1)])
    log_base = 2
    for site in sites.keys():
        (n,s,s_fwd, s_bck) = sites[site]
        if (isExtremal(s, n, cutoff, min_n_for_s, max_s_for_n)):
            print site, 'alt_alleles='+str(n), 'parsimony_score='+str(s),\
                    'parsimony_score_forward='+str(s_fwd),\
                    'parsimony_score_backward='+str(s_bck)
            x.append(math.log(n,log_base))
            y.append(s)

    m,b = np.polyfit(np.array(x), np.array(y), 1)
    print '\nPhylogenetic instability (log-' + str(log_base)+' slope) =', m
    
