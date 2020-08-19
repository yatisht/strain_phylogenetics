
cmake_minimum_required (VERSION 3.8)                                                                                                                                                                            

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++11 -Wall -g")
set(CMAKE_INCLUDE_CURRENT_DIR ON) 

find_package(Protobuf REQUIRED HINTS ${Protobuf_PATH}) 

# Print version messages
if(Protobuf_FOUND)
    message(STATUS "Using Protocol Buffers ${Protobuf_VERSION}")
endif()


include_directories(${Protobuf_INCLUDE_DIRS})

include(${TBB_ROOT}/cmake/TBBBuild.cmake)
tbb_build(TBB_ROOT ${TBB_ROOT} CONFIG_DIR TBB_DIR MAKE_ARGS tbb_cpf=1)
find_package(TBB REQUIRED tbbmalloc tbbmalloc_proxy tbb_preview)

FIND_PACKAGE(Boost COMPONENTS program_options iostreams REQUIRED)

INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})

find_package(OpenMP REQUIRED)

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

add_executable(find_sister_clades
    src/tree.cpp
    src/find_sister_clades.cpp
    )

protobuf_generate(
    LANGUAGE cpp
    TARGET add_missing_samples 
    PROTOS parsimony.proto)

TARGET_COMPILE_OPTIONS(find_parsimonious_assignments PRIVATE -DTBB_SUPPRESS_DEPRECATED_MESSAGES)
TARGET_LINK_LIBRARIES(find_parsimonious_assignments PRIVATE stdc++  ${Boost_LIBRARIES} ${TBB_IMPORTED_TARGETS})

if(SAVE_PROFILE)
    if(DEBUG)
        TARGET_COMPILE_OPTIONS(add_missing_samples PRIVATE -DTBB_SUPPRESS_DEPRECATED_MESSAGES -DSAVE_PROFILE=1 -DDEBUG=1)
    else(DEBUG)
        TARGET_COMPILE_OPTIONS(add_missing_samples PRIVATE -DTBB_SUPPRESS_DEPRECATED_MESSAGES -DSAVE_PROFILE=1)
    endif(DEBUG)
else(SAVE_PROFILE)
    if(DEBUG)
        TARGET_COMPILE_OPTIONS(add_missing_samples PRIVATE -DTBB_SUPPRESS_DEPRECATED_MESSAGES -DDEBUG=1)
    else(DEBUG)
        TARGET_COMPILE_OPTIONS(add_missing_samples PRIVATE -DTBB_SUPPRESS_DEPRECATED_MESSAGES)
    endif(DEBUG)
endif(SAVE_PROFILE)

TARGET_LINK_LIBRARIES(add_missing_samples PRIVATE stdc++  ${Boost_LIBRARIES} ${TBB_IMPORTED_TARGETS} ${Protobuf_LIBRARIES} OpenMP::OpenMP_CXX)

TARGET_LINK_LIBRARIES(find_sister_clades PRIVATE stdc++  ${Boost_LIBRARIES} OpenMP::OpenMP_CXX)

