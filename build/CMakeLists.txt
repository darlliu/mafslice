cmake_minimum_required(VERSION 2.6)
project(mafslice)

#set(CMAKE_CXX_COMPILER "clang++")
#set(CMAKE_CXX_FLAGS "-O -std=c++1y")
set(CMAKE_CXX_COMPILER "g++")
set(CMAKE_CXX_FLAGS "-Os -std=c++1y")

set(SOURCE_DIR "../src")
set(BOOST_ROOT "/home/yul13/gnu/boost_1_59_0")
set(BOOST_LIB_DIR "${BOOST_ROOT}/stage/lib")
set(KC_DIR "/home/yul13/bin/include/")
set(KC_LIB_DIR "/home/yul13/bin/lib")
set(PYTHON_ROOT "/auto/igb-libs/linux/centos/6.x/x86_64/pkgs/python/Python-2.7.10")
set(PYTHON_LIBRARY "${PYTHON_ROOT}/lib")
set(PYTHON_INCLUDE_DIR "${PYTHON_ROOT}/include/python2.7")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "../bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "../bin")

#configure_file("${SOURCE_DIR}/base.hpp")

INCLUDE_DIRECTORIES(
    ${SOURCE_DIR}
    ${BOOST_DIR}
    ${KC_DIR}
)

LINK_DIRECTORIES(
    ${BOOST_LIB_DIR}
    ${KC_LIB_DIR}
    )


message("Configuring core C++ binaries")
add_library(indexer SHARED ${SOURCE_DIR}/indexer.cpp ${SOURCE_DIR}/mafindexer.cpp)
target_link_libraries(indexer boost_filesystem boost_serialization kyotocabinet m pthread stdc++)

add_executable(seqslice.exe ${SOURCE_DIR}/main.cpp)
target_link_libraries(seqslice.exe indexer boost_program_options)
add_executable(mafslice.exe ${SOURCE_DIR}/main2.cpp)
target_link_libraries(mafslice.exe indexer boost_program_options)
FIND_PACKAGE(Boost 1.59.0)
FIND_PACKAGE(PythonLibs 2.7.10)
IF(Boost_FOUND AND PYTHONLIBS_FOUND)
    message("BOOST 1.59 AND PYTHON 2.7.10 FOUND, configuring python interface")
    INCLUDE_DIRECTORIES ("${Boost_INCLUDE_DIRS}")
    INCLUDE_DIRECTORIES ( "${PYTHON_INCLUDE_DIRS}" )
    SET(Boost_USE_STATIC_LIBS OFF)
    SET(Boost_USE_MULTITHREADED ON)
    SET(Boost_USE_STATIC_RUNTIME OFF)
    FIND_PACKAGE(Boost 1.59.0 COMPONENTS python)
    add_library(mafslice SHARED ${SOURCE_DIR}/mafslice.cpp)
    target_link_libraries(mafslice indexer ${Boost_LIBRARIES} )
    set_target_properties(mafslice PROPERTIES PREFIX "" )
ELSE()
    message("Error finding boost...")
ENDIF()
