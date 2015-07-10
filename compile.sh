#!/bin/bash
cd src
clang++  -I/home/yul13/bin/include -L/home/yul13/bin/lib indexer.hpp indexer.cpp main.cpp -lstdc++ -lm -lc -lkyotocabinet -lboost_program_options  -lboost_serialization -lboost_filesystem -O -std=c++1y 
mv a.out ../bin/mafslice.exe
cd ..
