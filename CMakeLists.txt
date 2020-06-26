cmake_minimum_required (VERSION 3.8)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++11 -Wall ")

include(${TBB_ROOT}/cmake/TBBBuild.cmake)
tbb_build(TBB_ROOT ${TBB_ROOT} CONFIG_DIR TBB_DIR MAKE_ARGS tbb_cpf=1)
find_package(TBB REQUIRED tbbmalloc tbbmalloc_proxy tbb_preview)

FIND_PACKAGE(Boost COMPONENTS program_options iostreams REQUIRED)

INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})

add_executable(find_parsimonious_assignments
    src/tree.cpp
    src/mapper.cpp
    src/printer.cpp
    src/find_parsimonious_assignments.cpp
    )

TARGET_COMPILE_OPTIONS(find_parsimonious_assignments PRIVATE -DTBB_SUPPRESS_DEPRECATED_MESSAGES)
TARGET_LINK_LIBRARIES(find_parsimonious_assignments PRIVATE stdc++  ${Boost_LIBRARIES} ${TBB_IMPORTED_TARGETS})
