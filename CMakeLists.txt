cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(lab2 VERSION 1.0)
add_executable(lab2 lab2/lab2.cpp)

project(sanya2 VERSION 1.0)
add_executable(sanya2 lab2/sanya.cpp)

project(lab3 VERSION 1.0)
add_executable(lab3 lab3/main.cpp)

project(lab4 VERSION 1.0)
add_executable(lab4 lab4/main.cpp)

project(sanya3 VERSION 1.0)
add_executable(sanya3 lab3/sanya.cpp)