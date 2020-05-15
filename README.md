# strain_phylogenetics

* Install prereqs
```
    $ pip install treelib
```
* Download example newick files
```
    $ for day in 19 20 21 22 23 24 25 26; do wget https://hgwdev.gi.ucsc.edu/~angie/nextstrainRecompute.2020-04-30/2020-04-${day}/nextstrain.nh -O tree/nextstrain-2020-04-${day}.nh; done
```
* Download example VCF files
```
    $ for day in 19 20 21 22 23 24 25 26; do wget https://hgwdev.gi.ucsc.edu/~angie/nextstrainRecompute.2020-04-30/2020-04-${day}/nextstrainSamples.vcf.gz -O vcf/nextstrain-2020-04-${day}.vcf.gz; done
```
* Run examples
```
    $ python compute_directed_change-cog.py -T1 tree/nextstrain-2020-04-19.nh -T2 tree/nextstrain-2020-04-20.nh -CORES=2
    $ python find_parsimonious_assignments.py -tree tree/nextstrain-2020-04-19.nh -vcf vcf/nextstrain-2020-04-19.vcf.gz 
```
