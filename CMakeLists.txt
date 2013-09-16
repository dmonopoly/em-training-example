cmake_minimum_required(VERSION 2.8)

# Options. Turn on with cmake -Dmyvarname=ON.
option(test "Build all tests." OFF)

set(PROJECT_NAME two_tag_example)
project(${PROJECT_NAME})

# Require c++11 and libc++ for more helpers like initializer_list.
set(CMAKE_CXX_FLAGS "-g -Wall -std=c++11 -stdlib=libc++")

# Libraries & Executables.
add_library(helper_lib Helper.cc)
add_library(notation_lib Notation.cc)
add_executable(${PROJECT_NAME} TwoTagExample.cc)

target_link_libraries(${PROJECT_NAME} helper_lib notation_lib)

# Testing.
if (test)
  add_executable(notation_test NotationTest.cc)
  target_link_libraries(notation_test notation_lib)
endif()
