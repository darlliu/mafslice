// This is the mafslice.so library used for python interfacing
// done via boost.python
//
#include"motifmapdb.hpp"
#include <boost/python.hpp>
using namespace boost::python;
//std::string (seqdb::*get)(const std::string& , const unsigned&, const unsigned&) = seqdb::get;
BOOST_PYTHON_MODULE(motifmapcore)
{
    PyEval_InitThreads();
    class_<motifmapdb>("motifmapdb")
        .def("init", &motifmapdb::init)
        .def("compute", &motifmapdb::compute)
        .def("compute_join", &motifmapdb::compute_join)
        .def("get_seq", &motifmapdb::get_seq)
        .def("get", &motifmapdb::get_inv);
};
