cmake_minimum_required(VERSION 3.9)
project(MorfiTest)

set(CMAKE_CXX_STANDARD 17)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++0x")

#include_directories("${CMAKE_SOURCE_DIR}/include/googletest")
add_subdirectory(include/googletest)
#enable_testing()
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})


#include_directories("F:/Achademic/CS/emester 05/Software Engineering Project/untitled6/src")

add_executable(MorfiTest file2ObjTest.cpp)

target_link_libraries(MorfiTest morficonvert)
target_link_libraries(MorfiTest gtest)