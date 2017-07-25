cmake_minimum_required(VERSION 3.7)
project(yass)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES
        src/align.c
        src/align.h
        src/assemble.c
        src/assemble.h
        src/avl.c
        src/avl.h
        src/cover.pl
        src/display.c
        src/display.h
        src/global_var.c
        src/global_var.h
        src/kword.c
        src/kword.h
        src/list.c
        src/list.h
        src/main.c
        src/parse_blout.pl
        src/prdyn.c
        src/prdyn.h
        src/proba.c
        src/proba.h
        src/red_black.c
        src/red_black.h
        src/regroup.c
        src/regroup.h
        src/threads.c
        src/threads.h
        src/tuple.c
        src/tuple.h
        src/util.c
        src/util.h
        src/yass2blast.pl src/PyYASS.cpp)

add_executable(yass ${SOURCE_FILES})