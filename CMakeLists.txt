cmake_minimum_required(VERSION 3.29)
project(SystemOfLianealEquations)

set(CMAKE_CXX_STANDARD 20)

add_executable(SystemOfLianealEquations test.cpp
       functions.cpp
        functions.h
        gauss.cpp
        squareRoots.cpp
        squareRoots.cpp
        seidel.cpp
)
