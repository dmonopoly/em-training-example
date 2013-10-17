# CMakeLists for TwoTagExample.
cmake_minimum_required(VERSION 2.8)

# Options. Turn on with cmake -Dmyvarname=ON.
option(test "Build all tests." OFF)

set(PROJECT_NAME two_tag_example)
project(${PROJECT_NAME})

# Require c++11 and libc++ for more helpers like initializer_list.
set(CMAKE_CXX_FLAGS "-g -Wall -std=c++11 -stdlib=libc++")

# Libraries & Executables.
add_library(basic_helper_lib BasicHelper.cc)
add_library(notation_lib Notation.cc)
add_library(edge_lib Edge.cc)
add_library(nlp_helper_lib NLPHelper.cc)
add_library(trellis_aid_lib TrellisAid.cc)

add_executable(${PROJECT_NAME} TwoTagExample.cc)

target_link_libraries(${PROJECT_NAME}
                      basic_helper_lib
                      notation_lib
                      edge_lib
                      nlp_helper_lib)
#                       trellis_aid_lib)

# Testing.
if (test)
  add_executable(notation_test NotationTest.cc)
  target_link_libraries(notation_test notation_lib)
endif()
