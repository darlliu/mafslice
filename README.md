# mafslice
multiple genome alignment database with index based access

* Index genome sequences and maf files into on disk db
* Provides fast access to underlying sequences with alignment information given a range of index
* Does not provide any sequence alignment functionality
* Python binding for the library


## Notes on compiling

The workflow is as follows:

0. Requirements are: boost 1.59+, gcc 4.8.2+ (matching boost compilation), kyotocabinet any recent version

1. Refer to compile.sh for how to compile dynamic libraries. 
    Important: Link as much as possible into the .so NOT later when compiling the python library.

2. Then edit Jamroot and boost-build.jam accordingly. Note: <token>VAL does not tolerate whitespace inbetween.

3. The structure is always cpp code -> .so (linking all necessary boost and other libraries) -> python adaptor (as simple as possible) -> python .so

