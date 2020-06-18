
### Install prerequisites 
```
    $ sudo apt-get install libboost-all-dev
    $ sudo apt-get install cmake
```
### Build code  
```
    $ git clone https://github.com/yatisht/strain_phylogenetics
    $ cd strain_phylogenetics
    $ git checkout dev-c++
    $ git clone https://github.com/01org/tbb
    $ mkdir build
    $ cd build
    $ cmake -DTBB_ROOT=${PWD}/../tbb .. 
    $ make -j
    $ ./find_parsimonious_assignments ../tree/pruned-sumtree-for-cog.nh ../vcf/pruned-sumtree-for-cog.vcf 2 
```
