sudo apt update 
sudo apt-get --yes install build-essential cmake libboost-all-dev libomp-dev libprotoc-dev libprotoc-dev protobuf-compiler

git clone https://github.com/01org/tbb
mkdir build
cd build
cmake  -DTBB_ROOT=${PWD}/../tbb  .. 
make -j
cd ..
    
