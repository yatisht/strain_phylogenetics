sudo apt update 
sudo apt-get --yes install build-essential
sudo apt-get --yes install cmake
sudo apt-get --yes install libboost-all-dev
sudo apt-get --yes install libomp-dev
sudo apt-get --yes install libprotoc-dev libprotoc-dev protobuf-compiler

git clone https://github.com/01org/tbb
mkdir build
cd build
cmake  -DTBB_DIR=${PWD}/../tbb  -DTBB_ROOT=${PWD}/../tbb .. 
make -j
cd ..
    
