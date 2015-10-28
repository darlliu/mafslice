#ifndef BASE
#define BASE

#include <vector>
#include <algorithm>
#include <cmath>
#include <map>
#include <kcpolydb.h>
#include <iostream>
#include <fstream>
#include <string>
#include <memory>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/vector.hpp>

#define USE_BINARY_ARCHIVE 0
#define USE_DBT 0

#if USE_BINARY_ARCHIVE
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#else
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#endif

#endif
