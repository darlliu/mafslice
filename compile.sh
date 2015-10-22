#!/bin/bash
cd src
#clang++ -I/home/yul13/gnu/boost_1_59_0 -I/home/yul13/bin/include -L/home/yul13/bin/lib -L/home/yul13/gcc/lib indexer.hpp indexer.cpp main.cpp -lstdc++ -lm -lc -lkyotocabinet -lboost_program_options -Wl,-rpath=/home/yul13/gcc/lib -lboost_serialization -Wl,-rpath=/home/yul13/gcc/lib -lboost_filesystem -Wl,-rpath=/home/yul13/gcc/lib -O -std=c++1y -Wl, --verbose
clang++ -I/home/yul13/gnu/boost_1_59_0 -I/home/yul13/bin/include  \
    -c indexer.cpp -fpic  -std=c++1y -O
clang++  \
    -shared -o libindexer.so indexer.o -L/home/yul13/bin/lib -L/home/yul13/gcc/lib \
    -lstdc++ -lm -lc -lkyotocabinet -lboost_serialization -Wl,-rpath=/home/yul13/gcc/lib -lboost_filesystem -Wl,-rpath=/home/yul13/gcc/lib\
    -O -std=c++1y
clang++ -o ../bin/mafslice.exe main.cpp -I/home/yul13/gnu/boost_1_59_0 -L/home/yul13/gcc/lib -L/home/yul13/codebase/mycodes/mafslice/src -lindexer  -lboost_program_options -Wl,-rpath=/home/yul13/gcc/lib  -O -std=c++1y
#mv a.out ../bin/mafslice.exe
cd ..
