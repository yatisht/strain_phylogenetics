## Tools for analyzing and comparing SARS-CoV2 phylogenies 

This repository provides tools primarily designed for analyzing Nextstrain (http://nextstrain.org/ncov/) and other phylogenies generated for the SARS-CoV2 genome but is also applicable for other phylogenetic applications. Detailed discussion of the tools can be found in the reference at the bottom of this page. Instructions below demonstrate the usage of these tools through examples.  

### Install prereqs
```
    $ pip install treelib
```
### Download sample newick and variant (VCF) files 
```
    $ for day in 19 20 21 22 23 24 25 26; do wget https://hgwdev.gi.ucsc.edu/~angie/nextstrainRecompute.2020-04-30/2020-04-${day}/nextstrain.nh -O tree/nextstrain-2020-04-${day}.nh; done
    $ for day in 19 20 21 22 23 24 25 26; do wget https://hgwdev.gi.ucsc.edu/~angie/nextstrainRecompute.2020-04-30/2020-04-${day}/nextstrainSamples.vcf.gz -O vcf/nextstrain-2020-04-${day}.vcf.gz; done
```
* The data wrangling for the above newick and VCF files from the Nextstrain data has been done by Angie Hinrichs (https://github.com/AngieHinrichs) at UCSC.

### Compute most-parsimonious state assignments in a tree for a variant list 
```
    $ python find_parsimonious_assignments.py -tree tree/nextstrain-2020-04-19.nh -vcf vcf/nextstrain-2020-04-19.vcf.gz > nextstrain-2020-04-19_PARSIMONY.txt 
```
* The above command completes in under *3 min* on a Macbook Pro with Intel Core i7. It reads the tree topology of the input Newick file and assigns an internal numeric label for each internal node (ignoring the internal labels and branch lengths if already provided by the input Newick). The first two lines of the output file print the input tree with internal nodes labelled in Newick format. The output is too large to display, so we view the first 1000 characters using the command below. 

```
    $ head -c 1000 nextstrain-2020-04-19_PARSIMONY.txt
    > Tree (in newick format) with internal nodes labelled: 
    > (EPI_ISL_406798|Wuhan/WH01/2019|Dec26,EPI_ISL_402120|Wuhan/IVDC-HB-04/2020|Jan1,EPI_ISL_402130|Wuhan/WIV07/2019|Dec30,EPI_ISL_403930|Wuhan/IPBCAMS-WH-03/2019|Dec30,EPI_ISL_402128|Wuhan/WIV05/2019|Dec30,EPI_ISL_402123|Wuhan/IPBCAMS-WH-01/2019|Dec24,EPI_ISL_421253|Jian/JX129/2020|Jan26,EPI_ISL_416425|Hangzhou/ZJU-07/2020|Feb3,EPI_ISL_402121|Wuhan/IVDC-HB-05/2019|Dec30,EPI_ISL_402132|Wuhan/HBCDC-HB-01/2019|Dec30,EPI_ISL_413851|Guangdong/2020XN4373-P0039/2020|Jan30,EPI_ISL_416320|Shanghai/SH0007/2020|Jan28,EPI_ISL_416389|Shanghai/SH0093/2020|Jan21,(EPI_ISL_416329|Shanghai/SH0020/2020|Feb1,EPI_ISL_413874|Guangdong/GDFS2020127-P0026/2020|Feb12)2,(EPI_ISL_402127|Wuhan/WIV02/2019|Dec30,(EPI_ISL_411955|USA/CA8/2020|Feb10,EPI_ISL_412898|Wuhan/HBCDC-HB-02/2019|Dec30)4)3,(EPI_ISL_413857|Guangdong/2020XN4448-P0002/2020|Jan31,EPI_ISL_416348|Shanghai/SH0039/2020|Feb6,(EPI_ISL_413882|Guangdong/GD2020139-P0007/2020|Feb2,EPI_ISL_414692|Guangzhou/GZM 
```
* The VCF file specifies allele information (REF(0)/ALT(1)) for each variant in each sample of the tree. For each variant/site in the VCF file, the output file then displays the ALT allele frequency, its parsimony score, the list of nodes (comma-separated, if its length is <=4) for which the branches leading to it have acquired forward (REF\>ALT) or backward (ALT\>REF) mutation and the size of the clades affected by those forward and backward mutations. An example line from the file is displayed below.
```
    > C5144T  alt_alleles=3   parsimony_score=2   forward_mutation_nodes=EPI_ISL_427354|Belgium/ULG-10103/2020|Apr6,1202  back_mutation_nodes=    forward_mutation_clade_sizes=1,2    back_mutation_clade_sizes=
```
* The last three lines of the output file display some aggregate statistics.
```
    $ tail -n 3 nextstrain-2020-04-19_PARSIMONY.txt 
    > Total leaf nodes:  8308
    > Total variants:  3246
    > Total parsimony score:  4297
```

### Identify extremal sites   
```
    $ python identify_extremal_sites.py -in nextstrain-2020-04-19_PARSIMONY.txt 
```
* The above command completes in under *1 sec* on a Macbook Pro. It can be used for identifying and flagging extremal sites i.e. sites having exceptional parsimony scores relative to their allele frequencies and therefore also suspected to contain systematic errors. The above command identifies 8 extremal sites (C11074T, G11083T, C28887T, A4050C, C6255T, T13402G, A3778G, C21575T) with a phylogenetic instability value of 3.19. For the precise definition of extremal sites and phyloegentic instability, refer to our manuscript referenced at the bottom. The code also provides an ability to ignore high-frequency C\>T and G\>T mutations using optional flags.
```
    $ python identify_extremal_sites.py -in nextstrain-2020-04-19_PARSIMONY.txt -ignoreCtoT=1 -ignoreGtoT=1
```

### Compute distances between tree pairs
```
    $ python compute_entropy_weighted_tree_distance.py -T1 tree/nextstrain-2020-04-19.nh -T2 tree/nextstrain-2020-04-20.nh -CORES=2 > dist_nextstrain-2020-04-19_nextstrain-2020-04-20.txt
```
* The above command completes in under *6 min* in a Macbook Pro. It computes the entropy-weighted total distance (asymmetric and symmetric) between two Nextstrain trees (dated 04/19/2020 and 04/20/2020). For detailed description of the distance metric, refer to our manuscript referenced at the bottom. The distance values are printed in the bottom 3 lines of the output file.
```
    $ tail -n 3 dist_nextstrain-2020-04-19_nextstrain-2020-04-20.txt
    > D(T1,T2) =  51.87327190863818
    > D(T2,T1) =  67.80089779914674
    > S(T1,T2) =  59.83708485389246
```
* The output file also provides the best-matching T2 branch(es), matching split distance and entropy-weighted matching split distance for each branch in T1. An example line in the output file displaying these values for T1 branch labelled 293 is as follows:
```
    > 293 267 3   0.0004062201871769792
```

### Rotate trees for tanglegrams
* A common visualization strategy for comparing two tree topologies is to use tanglegrams i.e. plot them side-by-side with common leaves connected by straight lines. A visually appealing tanglegram is one in which corresponding clades in both trees are arranged in the same vertical order and the there is minimum crossing of connecting lines with each other. While node rotation algorithms in the context of tanglegram visualization have been implemented in the [cophylo](https://www.rdocumentation.org/packages/phytools/versions/0.7-20/topics/cophylo) and [Dendroscope3](http://dendroscope.org/) tools, these algorithms are either slow or inadequate for large SARS-CoV2 phylogenies. We implemented a quick heuristic to produce vastly improved tanglegrams. 
```
    $ python rotate_trees.py -T1 tree/nextstrain-2020-04-19.nh -T2 tree/nextstrain-2020-04-20.nh -T1_out rot-nextstrain-2020-04-19.nh -T2_out rot-nextstrain-2020-04-20.nh
```
* The above command completes in under *1 min* and produces rotated trees (rot-nextstrain-2020-04-19.nh and rot-nextstrain-2020-04-20.nh) with a much improved tanglegram as seen below (images generated with the help of [cophylo](https://www.rdocumentation.org/packages/phytools/versions/0.7-20/topics/cophylo) tool). 

![Tanglegrams](/images/tanglegrams_comparison.png)

### Reference
* Coming soon.
