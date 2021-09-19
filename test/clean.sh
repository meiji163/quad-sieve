#/!bin/bash
{
    rm CMakeCache.txt
    rm -r CMakeFiles
    rm -r build
    rm Makefile
} &> /dev/null
