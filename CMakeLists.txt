cmake_minimum_required(VERSION 3.17)
project(proj)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)
include_directories(.idea)

add_executable(proj
        .idea/.gitignore
        .idea/workspace.xml
        FORSYTHE.h
        lagrange.cpp
        lagrange.h
        main.cpp
        quanc8.cpp
        quanc8.h
        SPLINES.CPP
        SPLINES.H
        )
