
cmake_minimum_required (VERSION 3.8)                                                                                                                                                                            

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++11 -Wall -g")
set(CMAKE_INCLUDE_CURRENT_DIR ON) 

#find_package(Protobuf REQUIRED HINTS /public/home/yatisht/strain_phylogenetics/protobuf-3.12.3/cmake/build/lib64/cmake/protobuf)
#find_package(Protobuf REQUIRED) 
find_package(Protobuf REQUIRED HINTS ${Protobuf_PATH}) 

# Print version messages
if(Protobuf_FOUND)
    message(STATUS "Using Protocol Buffers ${Protobuf_VERSION}")
    #message(STATUS "Inc dir ${Protobuf_INCLUDE_DIRS}")
    #message("Protobuf_LIBRARIES: " ${Protobuf_LIBRARIES})
    #message(STATUS "Compiler ${PROTOBUF_PROTOC_EXECUTABLE}")
endif()


include_directories(${Protobuf_INCLUDE_DIRS})
#protobuf_generate_cpp(PROTO_SRCS PROTO_HDRS addressbook.proto)

include(${TBB_ROOT}/cmake/TBBBuild.cmake)
tbb_build(TBB_ROOT ${TBB_ROOT} CONFIG_DIR TBB_DIR MAKE_ARGS tbb_cpf=1)
find_package(TBB REQUIRED tbbmalloc tbbmalloc_proxy tbb_preview)

FIND_PACKAGE(Boost COMPONENTS program_options iostreams REQUIRED)

INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})

add_executable(find_parsimonious_assignments
    src/tree.cpp
    src/fpa_mapper.cpp
    src/fpa_printer.cpp
    src/find_parsimonious_assignments.cpp
    )

add_executable(add_missing_samples
    src/tree.cpp
    src/ams_mapper.cpp
    src/add_missing_samples.cpp
    )

protobuf_generate(
    LANGUAGE cpp
    TARGET add_missing_samples 
    PROTOS parsimony.proto)

TARGET_COMPILE_OPTIONS(find_parsimonious_assignments PRIVATE -DTBB_SUPPRESS_DEPRECATED_MESSAGES)
TARGET_LINK_LIBRARIES(find_parsimonious_assignments PRIVATE stdc++  ${Boost_LIBRARIES} ${TBB_IMPORTED_TARGETS})

TARGET_COMPILE_OPTIONS(add_missing_samples PRIVATE -DTBB_SUPPRESS_DEPRECATED_MESSAGES)
TARGET_LINK_LIBRARIES(add_missing_samples PRIVATE stdc++  ${Boost_LIBRARIES} ${TBB_IMPORTED_TARGETS} ${Protobuf_LIBRARIES})
#TARGET_LINK_LIBRARIES(add_missing_samples PRIVATE stdc++  ${Boost_LIBRARIES} ${TBB_IMPORTED_TARGETS} ${PROTOBUF_LIBRARY})
#TARGET_LINK_LIBRARIES(add_missing_samples PRIVATE stdc++  ${Boost_LIBRARIES} ${TBB_IMPORTED_TARGETS})
