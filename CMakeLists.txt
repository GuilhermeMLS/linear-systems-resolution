cmake_minimum_required(VERSION 3.14)
project(labSisLin C)

set(CMAKE_C_STANDARD 11)

include_directories(.)

add_executable(labSisLin
        labSisLin.c
        SistemasLineares.c
        SistemasLineares.h
        utils.c
        utils.h)
