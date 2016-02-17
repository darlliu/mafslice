#ifndef MOTIFMAPDB
#define MOTIFMAPDB
#include "MOODS/moods_scan.h"
#include "indexer.hpp"
#include "mafindexer.hpp"
#include <Python.h>
#include <mutex>
#include <queue>
#include <thread>

std::map<std::string, std::string> OUTPUTS;
std::mutex OUTPUTS_MTX;
typedef std::vector<std::vector<double>> matrix;
template <typename T>
bool CMP(std::pair<unsigned, T> p, std::pair<unsigned, T> p2) {
  return p.first < p2.first;
};

// typedef std::priority_queue<std::pair<unsigned, seqdb>,
// std::list<std::pair<unsigned, seqdb>>, decltype(CMP<seqdb>)> SEQQ;
// typedef std::priority_queue<std::pair<unsigned, mafdb>,
// std::list<std::pair<unsigned, mafdb>>, decltype(CMP<mafdb>)> MAFQ;
typedef std::map<std::string, seqdb> SEQM;
class releaseGIL {
public:
  inline releaseGIL() { save_state = PyEval_SaveThread(); }

  inline ~releaseGIL() { PyEval_RestoreThread(save_state); }

private:
  PyThreadState *save_state;
};
class motifmapcompute {
public:
  motifmapcompute(const INTERVAL_PAIR &invs, const std::string &motif,
                  const matrix &mat, const double &th, const int &cnt)
      : invs(invs), motif(motif), mat(mat), th(th), id(cnt){};
  void realign(){};
  void score(INTERVAL_PAIR &inv);
  void write();
  void routine() {
    realign();
    score(invs);
    write();
    return;
  };
  std::thread spawn() { return std::thread(&motifmapcompute::routine, this); };

private:
  INTERVAL_PAIR invs;
  std::string motif;
  matrix mat;
  std::map<std::string, double> results;
  int id;
  double th;
};

class motifmapdb {
public:
  motifmapdb(const std::string &dbname, const std::string &dbp,
             const std::string &ref)
      : dbname(dbname), dbp(dbp), ref(ref), maf(dbname, dbp, ref),
        outdir("./out/"){};
  motifmapdb() : motifmapdb("test", "./test/test.kch", "mm10"){};
  void init(const std::string &dbp_, const std::string &dbname_,
            const std::string &ref_, const std::string &out);
  std::pair<INTERVAL_PAIR, INTERVAL_PAIR>
  get_intervals(const int &l, const int &r, const int &lf, const int &rf) {
    std::pair<INTERVAL_PAIR, INTERVAL_PAIR> out;
    std::cerr << "Extracting one interval" << std::endl;
    auto pp = maf.get_intervals(l - lf, r + rf);
    bool flag = true;
    for (; pp.first != pp.second; ++pp.first) {
      auto tmp = maf.extract_intervals(*pp.first);
      auto tmp2 = maf.filter_intervals(l - lf, r + rf, tmp);
      if (flag) {
        out.first = tmp;
        out.second = tmp2;
        flag = false;
      } else {
        combine_intervals(out.first.first, tmp.first);
        combine_intervals(out.second.first, tmp2.first);
        for (auto &it : tmp.second) {
          if (out.first.second.count(it.second.ref) == 0)
            out.first.second[it.second.ref] = it.second;
          else
            combine_intervals(out.first.second[it.second.ref], it.second);
        }
        for (auto &it : tmp2.second) {
          if (out.second.second.count(it.second.ref) == 0)
            out.second.second[it.second.ref] = it.second;
          else
            combine_intervals(out.second.second[it.second.ref], it.second);
        }
      }
    }
    return out;
  };
  std::string get_flank(interval &inv, const int &lf, const int &rf);
  void flank(const int &, const int &,
             std::pair<INTERVAL_PAIR, INTERVAL_PAIR> &);
  void compute_join() {
    for (auto &c : threads) {
      c.routine();
    }
    std::cerr << "Writing outputs to json files" << std::endl;
    for (auto &it : OUTPUTS) {
      std::cerr << "Got " << it.first << std::endl;
      std::ofstream ofs(outdir + "/" + it.first + ".json");
      boost::replace_last(it.second, ",", "");
      ofs << "{" << it.second << "}";
      ofs.close();
    }
    OUTPUTS.clear();
    return;
  };
  void compute(const int &cnt, const std::string &motif, const matrix &mat,
               const double &th) {
    try {
      if (threads.size() > 100) {
        std::vector<std::thread> ths;
        std::cerr << "Releasing GIL" << std::endl;
        releaseGIL unlock = releaseGIL();
        std::cerr << "Too many threads, joining 100" << std::endl;
        for (auto &c : threads) {
          ths.push_back(c.spawn());
        }
        for (auto &th : ths)
          th.join();
        threads.clear();
      }
      std::cerr << " Creating compute object...";
      threads.push_back(motifmapcompute(inp, motif, mat, th, cnt));
    } catch (...) {
      std::cerr << "Compute failed on " << print_interval(inp.first)
                << std::endl;
    }
    return;
  };
  void add_ref(const std::string &ref) { refs.push_back(ref); }
  void get_flanks(const std::string &chr, const int &l, const int &r,
                  const int &lf, const int &rf) {
    std::cerr << "Now trying to calculate " << chr << " " << l << " " << r
              << std::endl;
    try {
      maf.set_chr(chr);
      auto invs = get_intervals(l, r, lf, rf);
      flank(lf, rf, invs);
      inp = invs.second;
    } catch (std::string s) {
      std::cerr << "Got error" << s << std::endl;
      return;
    } catch (...) {
      std::cerr << "Got unspecified exception..." << std::endl;
      return;
    }
  };
  void get_invs(const std::string &chr, const int &l, const int &r) {
    std::cerr << "Getting all intervals between " << l << " , " << r
              << std::endl;
    maf.set_chr(chr);
    auto pp = maf.get_intervals(l, r);
    for (; pp.first != pp.second; ++pp.first) {
      auto tmp = maf.extract_intervals(*pp);
      std::cerr << "Ref: " << print_interval(tmp.first) << std::endl;
      for (auto &it : tmp.second)
        std::cerr << "Maf: " << print_interval(it) << std::endl;
    }
    std::cerr << std::endl << std::endl;
  };
  std::string get_seq(const std::string &ref, const std::string &chr,
                      const int l, const int r, bool strand);

private:
  std::string dbname, dbp, ref, outdir;
  mafdb maf;
  SEQM seqm;
  INTERVAL_PAIR inp;
  std::vector<std::string> refs;
  std::vector<motifmapcompute> threads;
};

#endif
