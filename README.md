# mafslice
multiple genome alignment database with index based access (C++, python library)
used internally to replace pygr

* Index genome sequences and maf files into on disk db
* Provides fast access to underlying sequences with alignment information given a range index
* Does not provide any sequence alignment functionality
* Python adapters for the library


## Notes on compiling (Manual)

The workflow is as follows:

0. Requirements are: boost 1.59+, gcc 4.8.2+ (matching boost compilation), kyotocabinet any recent version

1. Refer to compile.sh for how to compile dynamic libraries. 
    Important: Link as much as possible into the .so NOT later when compiling the python library.

2. Then edit Jamroot and boost-build.jam accordingly. Note: `<TOKEN>VAL` does not tolerate whitespace inbetween.

3. The structure is always cpp code -> .so (linking all necessary boost and other libraries) -> python adaptor (as simple as possible) -> python .so

### Alternative (cmake, recommended)

**NOTE**: This method links libraries slightly differently. Don't know how to clean up yet but it still works.

0. Set up paths in cmake (refer to `build/CMakeLists.txt`) including boost path (whole build with sources), python paths (2.7.10+), kyotocabinet path, correct c++ compiler etc

1. in `build/` do `cmake CMakeLists && make`

2. go to `bin/` to find the libraries. Note: runtime environments must bs set correctly
