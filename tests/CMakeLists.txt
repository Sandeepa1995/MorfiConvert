cmake_minimum_required(VERSION 3.9)
project(MorfiTest)

set(CMAKE_CXX_STANDARD 17)

if (MSVC)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /MP")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
endif()

#Pre-made
# Check if OpenMP is supported by the compiler
# If so then assign flags accordingly
find_package(OpenMP)
if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif ()

include_directories(${MorfiConvert_SOURCE_DIR}/src)

add_subdirectory(include/googletest)
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

add_executable(MorfiTest TestMorfi.cpp)

target_link_libraries(MorfiTest morficonvert)
target_link_libraries(MorfiTest gtest)
target_link_libraries(MorfiTest -fopenmp)