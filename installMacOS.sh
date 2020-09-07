brew install cmake boost libomp protobuf

git clone https://github.com/01org/tbb
mkdir build
cd build
ls ${PWD}/../tbb
ls ${PWD}/../tbb/cmake
cmake -DTBB_DIR=${PWD}/../tbb/cmake -DTBB_ROOT=${PWD}/../tbb -DOpenMP_CXX_FLAGS="-Xclang -fopenmp -I/usr/local/opt/libomp/include" -DOpenMP_C_FLAGS="-Xclang -fopenmp -I/usr/local/opt/libomp/include" -DOpenMP_CXX_LIB_NAMES=libomp  -DOpenMP_C_LIB_NAMES=libomp  -DOpenMP_libomp_LIBRARY=/usr/local/opt/libomp/lib/libomp.dylib -DCMAKE_SHARED_LINKER_FLAGS="-L/usr/local/opt/libomp/lib -lomp -Wl,-rpath,/usr/local/opt/libomp/lib" ..
make -j
cd ..
