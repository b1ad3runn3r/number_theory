cmake_minimum_required(VERSION 3.20)

project("Lab 1" CXX)
set(CMAKE_CXX_STANDARD 20)

# find_package(Boost REQUIRED COMPONENTS numeric)

add_subdirectory(src/Solver)

add_executable(lab1 src/main.cpp)
target_link_libraries(lab1 solver)
