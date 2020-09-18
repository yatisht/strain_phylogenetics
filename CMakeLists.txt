
cmake_minimum_required (VERSION 3.8) 

project(strain_phylogenetics)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++11 -Wall -g -Wno-unused-function")
set(CMAKE_INCLUDE_CURRENT_DIR ON) 

# Print version messages
if(Protobuf_FOUND)
    message(STATUS "Using Protocol Buffers ${Protobuf_VERSION}")
endif()


include_directories(${Protobuf_INCLUDE_DIRS})

include(${TBB_DIR}/cmake/TBBBuild.cmake)
tbb_build(TBB_ROOT ${TBB_DIR} CONFIG_DIR TBB_DIR MAKE_ARGS tbb_cpf=1)
find_package(TBB REQUIRED tbbmalloc tbbmalloc_proxy tbb_preview)

FIND_PACKAGE(Boost COMPONENTS program_options iostreams REQUIRED)

find_package(OpenMP REQUIRED)

INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})

add_executable(rotate_trees
    src/tree.cpp
    src/rotate_trees.cpp
    )

add_executable(find_parsimonious_assignments
    src/tree.cpp
    src/fpa_mapper.cpp
    src/fpa_printer.cpp
    src/find_parsimonious_assignments.cpp
    )

TARGET_COMPILE_OPTIONS(rotate_trees PRIVATE -DTBB_SUPPRESS_DEPRECATED_MESSAGES)
TARGET_COMPILE_OPTIONS(find_parsimonious_assignments PRIVATE -DTBB_SUPPRESS_DEPRECATED_MESSAGES)

TARGET_LINK_LIBRARIES(rotate_trees PRIVATE stdc++  ${Boost_LIBRARIES} ${TBB_IMPORTED_TARGETS})
TARGET_LINK_LIBRARIES(find_parsimonious_assignments PRIVATE stdc++  ${Boost_LIBRARIES} ${TBB_IMPORTED_TARGETS})


