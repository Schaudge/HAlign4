cmake_minimum_required(VERSION 3.8)

project("halign4")

set(CMAKE_FIND_LIBRARY_SUFFIXES ".a" CACHE STRING "Set library suffixes to .a")
set(BUILD_SHARED_LIBS OFF)

find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    message(STATUS "OpenMP found")
    set(OPENMP ON CACHE BOOL "Enable OpenMP")
endif()

set(WFA_SRC "${CMAKE_SOURCE_DIR}/PairwiseAlignment/WFA2-lib")
set(WFA_LIB "${WFA_SRC}/lib/libwfa.a")
set(WFA_CPP_LIB "${WFA_SRC}/lib/libwfacpp.a")

add_custom_command(
    OUTPUT ${WFA_CPP_LIB}
    COMMAND ${CMAKE_COMMAND} -E make_directory ${WFA_SRC}/lib  # 创建 lib 目录
    COMMAND ${CMAKE_COMMAND} -E make_directory ${WFA_SRC}/build/cpp  # 创建 build/cpp 目录
    COMMAND make -C ${WFA_SRC} lib_wfa BUILD_EXAMPLES=0 BUILD_TOOLS=0
    WORKING_DIRECTORY ${WFA_SRC}
    COMMENT "Building libwfacpp.a in PairwiseAlignment/WFA2-lib"
)

add_custom_target(build_wfa ALL DEPENDS ${WFA_CPP_LIB})

add_library(wfa_lib STATIC IMPORTED)
add_dependencies(wfa_lib build_wfa)

set_target_properties(wfa_lib PROPERTIES IMPORTED_LOCATION "${WFA_CPP_LIB}")

include_directories(
    ${CMAKE_SOURCE_DIR}
    ${WFA_SRC}
)

set(SOURCE_FILES
    "./PairwiseAlignment/NeedlemanWunshReusable.cpp"
    "./SuffixArray/parallel_import.cpp"
    "./Utils/Arguments.cpp"
    "./Utils/Fasta.cpp"
    "./Utils/Graph.cpp"
    "./Utils/Insertion.cpp"
    "./Utils/NucleicAcidColumn.cpp"
    "./Utils/Utils.cpp"
    "./multi-thread/multi.cpp"
    "./StarAlignment/StarAligner.cpp"
    "stmsa.cpp"
)

add_executable(halign4 ${SOURCE_FILES})

add_dependencies(halign4 build_wfa)

set_target_properties(halign4 PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)

target_compile_options(halign4 PRIVATE -w -Wall -g -fPIE -O3)

if(ENABLE_DEBUG)
    target_compile_definitions(halign4 PRIVATE DEBUG)
endif()

if(ENABLE_M64)
    target_compile_definitions(halign4 PRIVATE M64)
endif()

target_link_libraries(halign4 PRIVATE wfa_lib pthread)

if(OpenMP_CXX_FOUND)
    target_link_libraries(halign4 PRIVATE OpenMP::OpenMP_CXX)
endif()
