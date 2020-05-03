# strain_phylogenetics

* Install prereqs
```
    $ pip install treelib
```
* Run example
```
    $ wget -rkpN -e robots=off -np -A nh,json https://hgwdev.gi.ucsc.edu/~angie/nextstrainRecompute.2020-04-30/
    $ python compute_directed_change.py -T1 T1.nh -T2 T2.nh
```
