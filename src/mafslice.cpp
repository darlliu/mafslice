// This is the mafslice.so library used for python interfacing
// done via boost.python
//
#include<boost/python.hpp>
#include"manager.hpp"
using namespace boost::python;
//std::string (seqdb::*get)(const std::string& , const unsigned&, const unsigned&) = seqdb::get;
BOOST_PYTHON_MODULE(mafslice)
{
    class_<manager>("manager")
        .def("init", &manager::init)
        .def("compute", &manager::compute)
        .def("get", &manager::get);
};
