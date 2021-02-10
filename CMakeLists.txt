cmake_minimum_required(VERSION 3.17)
project(proj)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)
include_directories(.idea)

add_executable(proj
        .idea/.gitignore
        .idea/workspace.xml
        Euler.cpp
        Euler.h
        fmin.cpp
        fmin.h
        FORSYTHE.h
        lagrange.cpp
        lagrange.h
        main.cpp
        MATRIX.CPP
        MATRIX.H
        quanc8.cpp
        quanc8.h
        rkf45.cpp
        rkf45.h
        SPLINES.CPP
        SPLINES.H
        zeroin.cpp
        zeroin.h)
