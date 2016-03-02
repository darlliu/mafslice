#ifndef MOTIFMAPDB_CPP
#define MOTIFMAPDB_CPP
#include "motifmapdb.hpp"
#include "MOODS/moods_scan.h"
#include <algorithm>

static std::vector<double> mybg = {0.29, 0.21, 0.21, 0.29};
static std::vector<double> myth = {0.0};

void motifmapdb::init(const std::string &dbp_, const std::string &dbname_,
                      const std::string &ref_, const std::string &out) {
  dbname = dbname_;
  outdir = out;
  dbp = dbp_;
  ref = ref_;
  maf.name = dbname;
  maf.dbpath = dbp;
  maf.ref = ref;
  std::cerr << "Trying to deserialize mafdb from " << dbname << std::endl;
  maf.load_db_kch(dbp, dbname);
  std::cerr << "Trying to deserialize reference seqdb" << std::endl;
  for (auto &rf : refs) {
    std::cerr << " trying to deserialize " << rf << std::endl;
    seqm[rf] = seqdb(rf, 1e4);
    seqm[rf].load_db_kch(dbp, rf);
  }
  return;
};

std::string get_reverse_comp(const std::string &in) {
  std::string out;
  out.reserve(in.size());
  for (auto &c : in) {
    switch (c) {
    case 'A': {
      out.push_back('T');
      break;
    }
    case 'T': {
      out.push_back('A');
      break;
    }
    case 'U': {
      out.push_back('A');
      break;
    }
    case 'G': {
      out.push_back('C');
      break;
    }
    case 'C': {
      out.push_back('G');
      break;
    }
    case 'a': {
      out.push_back('t');
      break;
    }
    case 't': {
      out.push_back('a');
      break;
    }
    case 'u': {
      out.push_back('a');
      break;
    }
    case 'g': {
      out.push_back('c');
      break;
    }
    case 'c': {
      out.push_back('g');
      break;
    }
    default:
      out.push_back(c);
    }
  }
  std::reverse(out.begin(), out.end());
  return out;
};

void motifmapcompute::write() {
  std::string ss = "\"" + std::to_string(id) + "\":{";
  // ss+=inv.first.chr+":"
  //+std::to_string(inv.first.l)+","+std::to_string(inv.first.r)+"\":";
  for (auto &it : results) {
    ss += "\"" + it.first + "\":" + std::to_string(it.second) + ",";
  }
  if (results.size() > 0)
    ss.pop_back();
  ss.push_back('}');
  std::lock_guard<std::mutex> guard(OUTPUTS_MTX);
  if (OUTPUTS.count(motif) == 0) {
    OUTPUTS[motif] = "";
  }
  OUTPUTS[motif] += ss + ",\n";
  return;
}
void motifmapcompute::score(INTERVAL_PAIR &inv) {

  auto best_match = [&](const std::vector<match> &matches,
                        const interval &inv) {
    double score = matches[0].score;
    if (matches.size() == 0)
      score = th;
    for (auto &m : matches)
      if (m.score > score) {
        score = m.score;
      }
    if (results.count(inv.ref) == 0)
      results[inv.ref] = score;
    else if (results[inv.ref] < score)
      results[inv.ref] = score;
  };

  // auto matches_ = MOODS::scan::scan_dna(inv.first.seq, my_matrix,mybg,
  // std::vector<double>(){-50});
  best_match(MOODS::scan::naive_scan_dna(inv.first.seq, mat, th), inv.first);
  best_match(
      MOODS::scan::naive_scan_dna(get_reverse_comp(inv.first.seq), mat, th),
      inv.first);

  for (auto &iv : inv.second) {
    // matches_ = MOODS::scan::scan_dna(iv.seq, my_matrix,mybg, myth);
    best_match(MOODS::scan::naive_scan_dna(iv.second.seq, mat, th), iv.second);
    best_match(
        MOODS::scan::naive_scan_dna(get_reverse_comp(iv.second.seq), mat, th),
        iv.second);
  }
  return;
}
std::vector<hit>
motifmapseq::moods_scan(const std::string &seq,
                        const std::vector<std::vector<double>> &mat,
                        const std::vector<double> &bg, const double &thr) {
  auto matches =
      MOODS::scan::scan_dna(seq, mat, bg, std::vector<double>(){thr});
  auto matches2 = MOODS::scan::scan_dna(get_reserve_comp(seq), mat, bg,
                                        std::vector<double>(){thr});
  matches.insert(matches.end(), matches2.begin(), matches2.end());
  std::vector<hit> out;
  for (auto &m : matches)
    out.push_back(std::pair<int, double>{m.pos, m.score});
  return out;
};
std::vector<hit>
motifmapseq::moods_naive_scan(const std::string &seq,
                              const std::vector<std::vector<double>> &mat,
                              const double &thr) {
  auto matches = MOODS::scan::naive_scan_dna(seq, mat, thr);
  auto matches2 = MOODS::scan::naive_scan_dna(get_reserve_comp(seq), mat, thr);
  matches.insert(matches.end(), matches2.begin(), matches2.end());
  std::vector<hit> out;
  for (auto &m : matches)
    out.push_back(std::pair<int, double>{m.pos, m.score});
  return out;
};

void motifmapdb::flank(const int &lf, const int &rf,
                       std::pair<INTERVAL_PAIR, INTERVAL_PAIR> &in) {
  auto inner = [&](interval &inv, interval &inv2, const int &lf,
                   const int &rf) {
    if ((inv2.l > lf) && (inv2.r + rf + inv2.l < inv.seq.size())) {
      inv2.seq = inv.seq.substr(inv2.l - lf, inv2.r + rf + lf);
    } else {
      int start = inv.l + inv2.l;
      int stop = inv.l + inv2.l + inv2.r;
      inv2.l = start;
      inv2.r = stop;
      if (!inv2.strand) {
        auto sz = seqm[inv2.ref].sizes[inv2.chr];
        start = sz - inv.r;
        stop = sz - inv.l - 1;
        inv2.l = start;
        inv2.r = stop;
      }
      inv2.seq = get_flank(inv2, lf, rf);
    }
#if DEBUG
    std::cerr << " Flanked " << print_interval(inv2) << std::endl;
#endif
  };
  // std::cerr << "Before flank"<<print_interval(in.second.first)<<std::endl;
  inner(in.first.first, in.second.first, lf, rf);
  for (auto &it : in.first.second) {
    // std::cerr << "Before
    // flank"<<print_interval(in.second.second[i])<<std::endl;
    inner(in.first.second[it.first], in.second.second[it.first], lf, rf);
  }
  return;
}

std::string motifmapdb::get_flank(interval &inv, const int &lf, const int &rf) {
  if (seqm.count(inv.ref) == 0) {
    try {
      std::cerr << "Adding a seqdb: " << inv.ref << std::endl;
      seqm[inv.ref] = seqdb(inv.ref, 1e4);
      seqm[inv.ref].load_db_kch(dbp, inv.ref);
    } catch (std::string err) {
      std::cerr << "Got error " << err << std::endl;
      return "";
    }
  }
  // std::cerr <<" get flank .."<<inv.ref;
  if (inv.strand)
    return seqm[inv.ref].get(inv.chr, inv.l - lf, inv.r + rf);
  else
    return seqm[inv.ref].get(inv.chr, inv.l - lf, inv.r + rf);
  // std::cerr <<".. done";
}

std::string motifmapdb::get_seq(const std::string &ref, const std::string &chr,
                                const int l, const int r, bool strand) {
  auto sq = seqm[ref];
  std::string ss;
  if (strand)
    ss = sq.get(chr, l, r);
  else {
    auto sz = sq.sizes[chr];
    int start = sz - r;
    int stop = sz - l;
    ss = sq.get(chr, start, stop);
    ss = get_reverse_comp(ss);
  }
  return ss;
}

#endif
