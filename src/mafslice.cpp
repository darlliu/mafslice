// This is the mafslice.so library used for python interfacing
// done via boost.python
//
#include<boost/python.hpp>
#include"indexer.hpp"
using namespace boost::python;
//std::string (seqdb::*get)(const std::string& , const unsigned&, const unsigned&) = seqdb::get;
BOOST_PYTHON_MODULE(mafslice)
{
    class_<seqdb>("seqdb")
        .def("get", &seqdb::get3)
        .def("load", &seqdb::load_db);
};
