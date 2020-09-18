sudo yum group install -y  "Development Tools"                                                                                                                                                                  
sudo yum install -y wget
sudo yum install -y boost-devel
sudo yum install -y python3

sudo pip3 install treelib scipy numpy

# install cmake-3.18
wget https://github.com/Kitware/CMake/releases/download/v3.18.2/cmake-3.18.2.tar.gz
tar -xvzf cmake-3.18.2.tar.gz
cd cmake-3.18.2
./bootstrap --prefix=${PWD} --  -DCMAKE_USE_OPENSSL=OFF
make -j
make install
cd ..

# get TBB
wget https://github.com/oneapi-src/oneTBB/releases/download/2019_U9/tbb2019_20191006oss_lin.tgz
tar -xvzf tbb2019_20191006oss_lin.tgz

# build programs
mkdir -p build
cd build
../cmake-3.18.2/bin/cmake  -DTBB_DIR=${PWD}/../tbb2019_20191006oss -DTBB_ROOT=${PWD}/../tbb2019_20191006oss -DCMAKE_PREFIX_PATH=${PWD}/../tbb2019_20191006oss/cmake ..
make -j
cd ..
