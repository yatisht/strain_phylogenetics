sudo apt update 
sudo apt-get --yes install build-essential
sudo apt-get --yes install libboost-all-dev
sudo apt-get --yes install libomp-dev 
sudo apt-get --yes install libprotoc-dev libprotoc-dev protobuf-compiler

git clone https://github.com/01org/tbb
mkdir build
pushd build
cmake  -DTBB_ROOT=${PWD}/../tbb  .. 
make -j
popd
    
