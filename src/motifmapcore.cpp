// This is the mafslice.so library used for python interfacing
// done via boost.python
//
#include<boost/python.hpp>
#include"motifmapdb.hpp"
using namespace boost::python;
//std::string (seqdb::*get)(const std::string& , const unsigned&, const unsigned&) = seqdb::get;
BOOST_PYTHON_MODULE(motifmapcore)
{
    class_<motifmapdb>("motifmapdb")
        .def("init", &motifmapdb::init)
        .def("compute", &motifmapdb::compute)
        .def("get_seq", &motifmapdb::get_seq)
        .def("get", &motifmapdb::get_inv);
};
