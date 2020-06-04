## Tools for analyzing and comparing SARS-CoV-2 phylogenies

This repository provides tools primarily designed for analyzing Nextstrain (http://nextstrain.org/ncov/) and other phylogenies generated for the SARS-CoV-2 genome but is also applicable for other phylogenetic applications. Detailed discussion of the tools can be found in the reference at the bottom of this page. Instructions below demonstrate the usage of these tools through examples.

### Install python prerequisites using pip
```
    $ pip install treelib numpy==1.14.6 scipy==1.0.1  
```
### Newick and Variant Call Format (VCF) input files

Example files are provided in the subdirectories tree/ and vcf/ .  Phylogenetic trees must be valid [Newick](https://en.wikipedia.org/wiki/Newick_format)-formatted files and variants must be [Variant Call Format (VCF)](https://en.wikipedia.org/wiki/Variant_Call_Format) files that include sample genotypes.  The leaf labels in the Newick files must exactly match the sample names in the VCF header line starting with "#CHROM".  Files in these formats generated from [Nextstrain](https://nextstrain.org/ncov) data are available from the [UCSC Genome Browser](https://genome.ucsc.edu/cgi-bin/hgTrackUi?db=wuhCor1&g=nextstrainSamples).

### Compute most-parsimonious state assignments in a tree for a variant list 
```
    $ python find_parsimonious_assignments.py -tree tree/pruned-sumtree-for-cog.nh -vcf vcf/pruned-sumtree-for-cog.vcf > pruned-sumtree-for-cog_PARSIMONY.txt
```
* The above command completes in under 2 sec on a Macbook Pro with Intel Core i7. It reads the tree topology of the input Newick file and assigns an internal numeric label for each internal node (ignoring the internal labels and branch lengths if already provided by the input Newick). The first two lines of the output file print the input tree with internal nodes labelled in Newick format. The output is too large to display, so we view the first 1000 characters using the command below.

```
    $ head -c 1000 pruned-sumtree-for-cog_PARSIMONY.txt
    > Tree (in newick format) with internal nodes labelled: 
    > ((Shanghai_SH0007_2020,Hangzhou_ZJU-07_2020,Wuhan_IPBCAMS-WH-03_2019,Wuhan_WH01_2019,Wuhan_WIV07_2019,Shanghai_SH0093_2020,Wuhan_HBCDC-HB-01_2019,Wuhan_IPBCAMS-WH-01_2019,Wuhan_IVDC-HB-04_2020,Wuhan_IVDC-HB-05_2019,Wuhan_WIV05_2019,Sweden_01_2020,Taiwan_2_2020,USA_CA2_2020,Shanghai_SH0040_2020,Singapore_7_2020,Australia_VIC02_2020,France_IDF0515_2020,Hangzhou_ZJU-01_2020,Shanghai_SH0037_2020,Nepal_61_2020,Wuhan_HBCDC-HB-03_2019,Wuhan_IVDC-HB-01_2019,Wuhan_WIV06_2019,Wuhan_WIV04_2019,China_WH-09_2020,Wuhan_IPBCAMS-WH-04_2019,Wuhan_IPBCAMS-WH-02_2019,Wuhan_WH03_2020,Jiangxi_IVDC-JX-002_2020,Zhejiang_WZ-01_2020,Japan_KY-V-029_2020,Hangzhou_ZJU-03_2020,Zhejiang_WZ-02_2020,Netherlands_Utrecht_12_2020,Nonthaburi_61_2020,Nonthaburi_74_2020,Hangzhou_ZJU-05_2020,Hangzhou_HZ-1_2020,Singapore_1_2020,Taiwan_NTU02_2020,England_SHEF-BFD36_2020,France_IDF0626_2020,Malaysia_MKAK-CL-2020-7554_2020,Cambodia_0012_2020,Finland_1_2020,England_200641094

```
* The VCF file specifies allele information (REF(0)/ALT(1)) for each variant in each sample of the tree. For each variant/site in the VCF file, the output file then displays the ALT allele frequency, its parsimony score, the list of nodes (comma-separated, if its length is <=4) for which the branches leading to it have acquired forward (REF\>ALT) or backward (ALT\>REF) mutation, the sizes of the clades affected by those forward and backward mutations and a list of flagged leaves which are affected by a mutation affecting 3 or fewer leaves. An example line from the file is displayed below.
```
    > G11083T	alt_alleles=61	parsimony_score=13	forward_mutation_nodes=France_IDF0515_2020,Hangzhou_ZJU-01_2020,Shanghai_SH0037_2020,21,37,48,54,160,81,Switzerland_SZ1417_2020,Brazil_SPBR-12_2020	back_mutation_nodes=Anhui_SZ005_2020,Brazil_SPBR-10_2020	forward_mutation_clade_sizes=1,1,1,26,4,6,15,2,3,1,1	back_mutation_clade_sizes=1,1	flagged_leaves=France_IDF0515_2020,Hangzhou_ZJU-01_2020,Shanghai_SH0037_2020,Anhui_SZ005_2020,Brazil_SPBR-10_2020,Switzerland_SZ1417_2020,Brazil_SPBR-12_2020

```
* The last three lines of the output file display some aggregate statistics.
```
    $ tail -n 3 pruned-sumtree-for-cog_PARSIMONY.txt
    > Total leaf nodes:  422
    > Total variants:  423
    > Total parsimony score:  498
```

### Identify extremal sites   
```
    $ python identify_extremal_sites.py -in pruned-sumtree-for-cog_PARSIMONY.txt
```
* The above command completes in under 1 sec on a Macbook Pro. It can be used for identifying and flagging extremal sites i.e. sites having exceptional parsimony scores relative to their allele frequencies and therefore also suspected to contain systematic errors. The above command identifies 6 extremal sites (C11074T, C27046T, T13402G, A3778G, G24390C, G26144T) with a phylogenetic instability value of 3.03. For the precise definition of extremal sites and phylogenetic instability, refer to our manuscript referenced at the bottom. The code also provides an ability to ignore high-frequency C\>T and G\>T mutations using optional flags.
```
    $ python identify_extremal_sites.py -in pruned-sumtree-for-cog_PARSIMONY.txt -ignoreCtoT=1 -ignoreGtoT=1
```
* The above command identifies three extremal sites (T13402G, A3778G, G24390C) with a phylogenetic instability value of 2.32.

### Plot Extremal Sites  
```
    $ python generate_plot_extremal_sites_data.py -in pruned-sumtree-for-cog_PARSIMONY.txt > plot_extremal_sites_data.txt
    $ Rscript --vanilla plot_parsimony.r plot_extremal_sites_data.txt extremal_sites_plot.pdf
```
* The above python command first creates raw input data for the extremal sites plot. Next, the R command (which should be executed after installing the *plyr* package) accepts the generated data and creates a log(allele count) by parsimony plot for all variants sites in a given vcf. It produces three plots, one of all data, one ignoring C>U mutations and one ignoring C>U and G>U mutations, as shown below. 
![Extremal](/images/extremal.png)

### Compute distances between tree pairs
```
    $ python compute_entropy_weighted_tree_distance.py -T1 tree/pruned-cog-for-sumtree.nh  -T2 tree/pruned-sumtree-for-cog.nh -CORES=2 > dist_pruned-cog-for-sumtree_pruned-sumtree-for-cog.txt
```
* The above command completes in under 2 sec on a Macbook Pro. It computes the entropy-weighted total distance (asymmetric and symmetric) between two trees pruned to contain a common set of 422 SARS-CoV-2 samples: the [COG-UK](https://www.cogconsortium.uk/data/) tree from April 24, 2020 and a consensus tree made using the [SumTrees](https://dendropy.org/programs/sumtrees.html) tool from 31 Nextstrain trees generated between March 23 and April 30, 2020. For detailed description of the distance metric, refer to our manuscript referenced at the bottom. The distance values are printed in the bottom 3 lines of the output file.
```
    $ tail -n 3 dist_pruned-cog-for-sumtree_pruned-sumtree-for-cog.txt
    > D(T1,T2) =  4.263117527652659
    > D(T2,T1) =  1.1491837172916255
    > S(T1,T2) =  2.7061506224721423
```
* The output file also provides the best-matching T2 branch(es), matching split distance and entropy-weighted matching split distance for each branch in T1. An example line in the output file displaying these values for T1 branch labelled 135 is as follows:
```
    > 135 119 6   0.09230354241562243
```

### Rotate trees for tanglegrams
* A common visualization strategy for comparing two tree topologies is to use tanglegrams i.e. plot them side-by-side with common leaves connected by straight lines. A visually appealing tanglegram is one in which corresponding clades in both trees are arranged in the same vertical order and the there is minimum crossing of connecting lines with each other. While node rotation algorithms in the context of tanglegram visualization have been implemented in the [cophylo](https://www.rdocumentation.org/packages/phytools/versions/0.7-20/topics/cophylo) and [Dendroscope3](http://dendroscope.org/) tools, these algorithms are either slow or inadequate for large SARS-CoV-2 phylogenies. We implemented a quick heuristic to produce vastly improved tanglegrams.
```
    $ python rotate_trees.py -T1 tree/pruned-sumtree-for-cog.nh -T2 tree/pruned-cog-for-sumtree.nh -T1_out rot-pruned-sumtree-for-cog.nh -T2_out rot-pruned-cog-for-sumtree.nh
```
* The above command completes in under 1 sec and produces rotated trees (rot-pruned-cog-for-sumtree.nh and rot-pruned-sumtree-for-cog.nh) with a much improved tanglegram as seen below (images generated with the help of [cophylo](https://www.rdocumentation.org/packages/phytools/versions/0.7-20/topics/cophylo) tool).

![Rot-tanglegrams](/images/tanglegrams_comparison.png)

### Pairwise merging of trees 
```
    $ python tree_merge.py -T1 tree/pruned-sumtree-for-cog.nh -T2 tree/pruned-cog-for-sumtree.nh -symmetric 1 -T_out symm-merged-sumtree-cog.nh  
```
* The above command completes in under 40 sec on a Macbook Pro and produces a merged tree (symm-merged-sumtree-cog.nh) from two input trees (pruned-sumtree-for-cog.nh and pruned-cog-for-sumtree.nh) that is maximally resolved and compatible with both input trees (refer to our manuscript referenced at the bottom for more details).  Below are the resulting tanglegrams of the resulting merged tree with the two input trees (after applying tree rotation). The above command can also be used without the symmetric flag for its asymmetric version (where first input tree is given a priority to resolve the merged tree) or using the intersectOnly flag that produces a simple consensus of the two input trees.  

![Merged-tanglegrams](/images/tanglegrams_merged.png)
### Reference
* Coming soon.
