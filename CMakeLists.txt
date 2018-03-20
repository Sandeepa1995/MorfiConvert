cmake_minimum_required(VERSION 3.9)
project(MorfiConvert)

set(CMAKE_CXX_STANDARD 17)

#Add the src folder to include its libraries
include_directories("${CMAKE_SOURCE_DIR}/src")
link_directories(src)
add_subdirectory(src)

set(SOURCE_FILES main.cpp)
add_executable(MorfiConvert ${SOURCE_FILES})
target_link_libraries(MorfiConvert morficonvert)
