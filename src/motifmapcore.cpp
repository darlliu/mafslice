// This is the mafslice.so library used for python interfacing
// done via boost.python
//
#include "motifmapdb.hpp"
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
using namespace boost::python;
// std::string (seqdb::*get)(const std::string& , const unsigned&, const
// unsigned&) = seqdb::get;
BOOST_PYTHON_MODULE(motifmapcore) {
  PyEval_InitThreads();
  class_<hit>("hit")
      .def_readwrite("pos", &hit::first)
      .def_readwrite("score", &hit::second);
  class_<std::vector<hit>>("HitVec").def(
      vector_indexing_suite<std::vector<hit>>());
  class_<std::vector<double>>("DoubleVec")
      .def(vector_indexing_suite<std::vector<double>>());
  class_<std::vector<std::vector<double>>>("DoubleMat")
      .def(vector_indexing_suite<std::vector<std::vector<double>>>());
  class_<motifmapdb>("motifmapdb")
      .def("init", &motifmapdb::init)
      .def("compute", &motifmapdb::compute)
      .def("compute_join", &motifmapdb::compute_join)
      .def("get_seq", &motifmapdb::get_seq)
      .def("add_ref", &motifmapdb::add_ref)
      .def("get", &motifmapdb::get_flanks)
      .def("get_invs", &motifmapdb::get_invs);
  class_<motifmapseq>("seqdb")
      .def("load", &motifmapseq::load)
      .def("get", &motifmapseq::get)
      .def("moods_scan", &motifmapseq::moods_scan)
      .def("moods_scan_naive", &motifmapseq::moods_naive_scan);
};
