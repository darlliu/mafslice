cmake_minimum_required(VERSION 2.6)
project(mafslice)

#set(CMAKE_CXX_COMPILER "clang++")
#set(CMAKE_CXX_FLAGS "-O -std=c++1y")
set(CMAKE_CXX_COMPILER "g++")
set(CMAKE_CXX_FLAGS "-O3 -std=c++1y")

set(SOURCE_DIR "../src")
#set(BOOST_ROOT "/auto/igb-libs/linux/centos/6.x/x86_64/pkgs/boost/alt-1.59.0/boost_1_59_0")
set(BOOST_ROOT "/home/yul13/gnu/boost_1_60_0")
set(BOOST_LIB_DIR "${BOOST_ROOT}/stage/lib")
set(KC_DIR "/home/yul13/bin/include")
set(KC_LIB_DIR "/home/yul13/bin/lib")
set(PYTHON_ROOT "/auto/igb-libs/linux/centos/6.x/x86_64/pkgs/python/2.7.10-gcc4.8.2")
set(PYTHON_LIBRARY "${PYTHON_ROOT}/lib")
set(PYTHON_INCLUDE_DIR "${PYTHON_ROOT}/include/python2.7")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "../bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "../bin")

#configure_file("${SOURCE_DIR}/base.hpp")

INCLUDE_DIRECTORIES(
    ${SOURCE_DIR}
    ${BOOST_ROOT}
    ${KC_DIR}
)

LINK_DIRECTORIES(
    ${BOOST_LIB_DIR}
    ${KC_LIB_DIR}
    )


message("Configuring core C++ binaries")
file (GLOB MOODS_SRC 
    "${SOURCE_DIR}/MOODS/*.h"
    "${SOURCE_DIR}/MOODS/*.cpp"
    )
add_library(moods SHARED ${MOODS_SRC})
INCLUDE_DIRECTORIES ("${SOURCE_DIR}/MOODS")
FIND_PACKAGE(Boost 1.59.0)
FIND_PACKAGE(PythonLibs 2.7.10)
add_library(indexer SHARED ${SOURCE_DIR}/indexer.cpp ${SOURCE_DIR}/mafindexer.cpp)
target_link_libraries(indexer boost_filesystem boost_serialization kyotocabinet m pthread stdc++)
add_executable(mafslice.exe ${SOURCE_DIR}/main.cpp)
target_link_libraries(mafslice.exe indexer boost_program_options)
IF(Boost_FOUND AND PYTHONLIBS_FOUND)
    message("BOOST 1.59 AND PYTHON 2.7.10 FOUND, configuring python interface")
    INCLUDE_DIRECTORIES ("${Boost_INCLUDE_DIRS}")
    INCLUDE_DIRECTORIES ( "${PYTHON_INCLUDE_DIRS}" )
    SET(Boost_USE_STATIC_LIBS OFF)
    SET(Boost_USE_MULTITHREADED ON)
    SET(Boost_USE_STATIC_RUNTIME OFF)
    FIND_PACKAGE(Boost 1.59.0 COMPONENTS python)
    add_library(motifmapdb SHARED ${SOURCE_DIR}/motifmapdb.cpp ${SOURCE_DIR}/motifmapdb.hpp)
    target_link_libraries(motifmapdb indexer moods ${Boost_LIBRARIES})
    add_library(motifmapcore SHARED ${SOURCE_DIR}/motifmapcore.cpp)
    target_link_libraries(motifmapcore motifmapdb ${Boost_LIBRARIES})
    set_target_properties(motifmapcore PROPERTIES PREFIX "" )
ELSE()
    message("Error finding boost...")
ENDIF()
