cmake_minimum_required(VERSION 3.22)
project(alpaca)

set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
set(CMAKE_CXX_STANDARD 20)

set(CMAKE_CXX_FLAGS "-Wall -Werror=uninitialized -Werror=return-type -Werror=return-local-addr")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -fPIC -march=native")
set(CMAKE_CXX_FLAGS_DEBUG "-g -fsanitize=address,undefined -fno-sanitize-recover=all -D_GLIBCXX_ASSERTIONS")

set(external_dir ${CMAKE_SOURCE_DIR}/external)

include_directories(${external_dir}/maxpre2/src)
include_directories(${external_dir}/cadical/src)
include_directories(${external_dir}/kissat/src)
include_directories(${external_dir}/SCIPpp/include)

link_libraries(fmt)
link_libraries(scip)

link_libraries(${external_dir}/maxpre2/src/lib/libmaxpre.a)
link_libraries(${external_dir}/cadical/build/libcadical.a)
link_libraries(${external_dir}/kissat/build/libkissat.a)
link_libraries(${external_dir}/SCIPpp/libScipPP.a)

add_executable(alpaca main.cpp)