
### Install prerequisites 
```
    $ sudo apt-get install libboost-all-dev
    $ sudo apt-get install cmake
```
### Install protobuf and build code  
```
    $ git clone https://github.com/yatisht/strain_phylogenetics
    $ cd strain_phylogenetics
    $ git checkout dev-c++
    $ wget https://github.com/protocolbuffers/protobuf/releases/download/v3.12.3/protobuf-cpp-3.12.3.tar.gz
    $ tar -xvzf protobuf-cpp-3.12.3.tar.gz 
    $ cd protobuf-3.12.3
    $ ./configure --prefix=${PWD}/install
    $ make -j
    $ make install
    $ cd cmake; mkdir build; cd build;
    $ cmake ..
    $ make -j
    $ cd ../../../
    $ git clone https://github.com/01org/tbb
    $ mkdir build
    $ cd build
    $ cmake  -DTBB_ROOT=${PWD}/../tbb   -DProtobuf_INCLUDE_DIRS=${PWD}/../protobuf-3.12.3/install/include/ -DProtobuf_LIBRARIES=${PWD}/../protobuf-3.12.3/cmake/build/libprotobuf.a -DProtobuf_PATH=${PWD}/../protobuf-3.12.3/cmake/build/lib64/cmake/protobuf .. 
    $ make -j
```
### Find parsimonious assignments for a given tree and vcf file 
```
    $ ./find_parsimonious_assignments ../tree/pruned-sumtree-for-cog.nh ../vcf/pruned-sumtree-for-cog.vcf 2 
```
### Add missing samples in a vcf file from given a tree and vcf file containing both tree leaves and missing samples 
```
    $ ./add_missing_samples --tree ../tree/pruned-sumtree-for-cog.nh --vcf ../vcf/pruned-sumtree-for-cog.vcf --threads 4
```
### Add missing samples in a vcf file from given a tree and vcf file containing both tree leaves and missing samples and save the tree and node assignments to a new file 
```
    $ ./add_missing_samples --tree ../tree/pruned-sumtree-for-cog.nh --vcf ../vcf/pruned-sumtree-for-cog.vcf --threads 4 --save-assignments pruned.pb
```
### Load the tree and node assignments from a previously saved file and add new samples to the tree from the input vcf file 
```
    $ ./add_missing_samples --vcf ../vcf/pruned-sumtree-for-cog.vcf --threads 4 --load-assignments pruned.pb
```


