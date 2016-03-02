#include "mafindexer.hpp"
#include <sstream>

void mafdb::import(const std::string &dirname) {
  using namespace boost::filesystem;
  auto fp = path(dirname);
  std::vector<std::string> fps, fns;
  try {
    if (exists(fp) && is_directory(fp)) {
      std::cerr << " Now reading " << dirname << std::endl;
    } else
      throw("Path does not exist : " + dirname);
  } catch (const filesystem_error &ex) {
    std::cerr << "Error opening path " << dirname << std::endl;
    throw(ex);
  }
  for (auto it = directory_iterator(fp); it != directory_iterator(); ++it) {
    if (!is_directory(it->path())) {
      auto fp = it->path().string();
      if (fp.find(".maf") == std::string::npos)
        continue;
      std::cout << "Found MAF file : " << fp << std::endl;
      fps.push_back(fp);
      fns.push_back(it->path().filename().string());

    } else {
      std::cout << "skipping directory : " << it->path().filename()
                << std::endl;
    }
  }
  std::cout << "Found files count: " << fps.size() << std::endl;
  for (std::string fn : fns) {
    boost::replace_last(fn, ".maf", "");
    chrs.push_back(fn);
  }
  if (assemble) {
    assemble_chrs();
    init_db(chrs, dbs);
    load_db_();
    return;
  }
  init_db(chrs, dbs);
  for (unsigned i = 0; i < fps.size(); ++i) {
    auto chr = chrs[i];
    set_chr(chr);
    auto chrp = fps[i];
    fapaths[chr] = chrp;
#if USE_DBT
    auto dbp = dbpath + "/" + chr + ".MSA.kct";
#else
    auto dbp = dbpath + "/" + chr + ".MSA.kch";
#endif
    dbpaths[chr] = dbp;
    import_chr();
  }
}
void mafdb::assemble_chrs() {
  std::cerr << "Assembling chromsomes from files" << std::endl;
  for (auto &chr : chrs) {
    postfixes[chr] = 1;
    for (auto &dbp : tmp_dbpaths) {
      if (dbp.find(chr + ".MSA.kch") != std::string::npos) {
        if (dbp.find(chr + ".MSA.kch.") == std::string::npos) {
          dbpaths[chr] = dbp;
          std::cerr << "Adding " << dbp << " For " << chr << std::endl;
          sizes[chr] = 2 << 20;
        } else if (dbp.find("index") == std::string::npos) {
          postfixes[chr]++;
          std::cerr << "Parts file found from " << dbp << std::endl;
        }
      }
    }
  }
  return;
};
void mafdb::import_chr(const std::string &fname) {
  using namespace boost::filesystem;
  auto fp = path(fname);
  if (exists(fp) && fname.find(".maf") != std::string::npos) {
    std::cerr << " Now reading " << fname << std::endl;
  } else
    throw("File does not exist: " + fname);
  auto fn = fp.filename().string();
  boost::replace_last(fn, ".maf", "");
  chrs.push_back(fn);
  set_chr(fn);
  fapaths[chr] = fp.string();
#if USE_DBT
  auto dbp = dbpath + "/" + chr + ".MSA.kct";
#else
  auto dbp = dbpath + "/" + chr + ".MSA.kch";
#endif
  dbpaths[chr] = dbp;
  init_db(chrs, dbs);
  import_chr();
}
void mafdb::import_chr() {
  auto chrfp = fapaths[chr];
  auto dbv = dbs[chr];
  size_t idx = 0, total = 0;
  std::ifstream ifs(chrfp);
  std::string tmp = "", line, foo, sign, seq, species;
  std::string key, val;
  std::map<std::string, std::string> data;
  std::vector<inode> inodes;
  float score;
  unsigned start, end, len, rstart, rend, rr[2];
  std::stringstream ss;
  if (!ifs.is_open())
    throw("Error Opening file!");
  ifs.seekg(0, ifs.end);
  size_t fsize = ifs.tellg();
  size_t SZLMT = 1 << 19;
  ifs.seekg(0, ifs.beg);
  std::cout << "Opened file " << chrfp << " On chr: " << chr
            << " Size: " << fsize << " Est. Recs: <" << fsize / 2000
            << std::endl;
  fsize /= 2000;
  inodes.reserve(fsize);
  if (fsize > SZLMT) {
    std::cout << "Size is too big, trimming into " << SZLMT << std::endl;
    fsize = SZLMT;
  }
  auto tune = [&](decltype(dbv[0]) &db, short pfx = 0) {
    auto dbp = dbpaths[chr];
    std::cerr << std::endl;
#if USE_DBT
    if (fsize > 65536) {
      std::cerr << " Tuning new bucket size (B+ Tree): " << fsize / 10
                << std::endl;
      db->tune_buckets(1LL * fsize / 10);
      db->tune_map(2LL << 30);
      // db->tune_defrag(8);
      std::cerr << " Tuning comparator..." << std::endl;
      db->tune_comparator(&CMPSZ);
    }
#else
    if (fsize * 4 > SZLMT) {
      std::cerr << " Tuning new bucket size (Hash table): " << fsize * 2
                << std::endl;
      db->tune_buckets(1LL * fsize * 2);
      // db->tune_options(kyotoca0binet::HashDB::TCOMPRESS);
      db->tune_map(2LL << 30);
      // db->tune_defrag(8);
    }
#endif
    auto dbp2 = dbp;
    if (pfx > 0)
      dbp2 = dbp + "." + std::to_string(pfx - 1);
    if (assemble) {
      if (!db->open(dbp2, _DB::OREADER))
        throw("Error opening DB: " + std::string(db->error().name()));

    } else {
      if (!db->open(dbp2, _DB::OWRITER | _DB::OCREATE))
        throw("Error opening DB: " + std::string(db->error().name()));
    }
    return dbp2;
  };
  auto db = dbv[0];
  tune(db);
  while (std::getline(ifs, line)) {
    foo.clear();
    if (line.size() == 0 || line[0] == '#')
      continue;
    ss.clear();
    ss.str(line);
    if (line[0] == 'a') {
      if (tmp.size() > 0) {
        // auto hs = hasher (rstart, rend);
        rr[0] = rstart;
        rr[1] = rend;
        key = std::string((char *)&rr, sizeof(rr));
        val = std::to_string(score) + " " + tmp;
        data[key] = val;
        idx++;
        if (idx % 5000 == 0) {
          save_index(chr, inodes);
          inodes.clear();
          std::cerr << "@";
          if (!assemble)
            db->set_bulk(data, false);
          std::cerr << idx << "/" << total << " .";
          data.clear();
          if (idx > SZLMT) {
            db->close();
            dbv.push_back(std::shared_ptr<_DB>(new _DB));
            db = dbv.back();
            auto dbp = tune(db, dbv.size());
            total += idx;
            std::cerr << "Creating a new db file: " << dbp << std::endl;
            idx = 0;
          }
        }
        // saving routine in here
      }
      tmp = "";
      boost::replace_last(line, "=", " ");
      ss >> foo >> foo >> score;
    } else if (line[0] == 'i') {
      continue;
    } else if (line[0] == 's') {
      ss >> foo >> species >> start >> end >> sign >> len >> seq;
      if (species.find(ref) != std::string::npos) {
        rstart = start;
        rend = end;
        inodes.push_back(inode(rstart, rend, dbv.size() - 1));
      }
      if (tmp.size() > 0)
        tmp += "\t";
      tmp += species + " " + std::to_string(start) + " " + std::to_string(end) +
             " " + sign + seq;
    } else
      continue;
  }
  {
    std::cerr << "@";
    if (!assemble)
      db->set_bulk(data, false);
    total += idx;
    std::cerr << idx << "/" << total << " .";
    data.clear();
  }
  ifs.close();
  std::cout << std::endl
            << "Loaded length " << idx << " Saving indexes..." << std::endl;
  save_index(chr, inodes);
  sizes[chr] = total;
  postfixes[chr] = dbv.size();
  inodes.clear();

  return;
}
void mafdb::save_index(const std::string &chr, const std::vector<inode> &v) {
  std::ofstream f(dbpaths[chr] + ".index",
                  std::ofstream::app | std::ofstream::binary);
  if (!f.good())
    throw("Error opening the index file!");
  // auto sz =v.size();
  // f.write((char*) &sz, sizeof(size_t));
  for (auto &it : v) {
    f.write((char *)&it.l, sizeof(unsigned));
    f.write((char *)&it.r, sizeof(unsigned));
    f.write((char *)&it.p, sizeof(short));
  }
  f.close();
  return;
}
void mafdb::clear_index(const std::string &chr) {
  auto msa = msatrees[chr];
  msa->clear();
  return;
}
void mafdb::load_index(const std::string &chr) {
  size_t cnt = 0;
  std::string key, val;
  auto msad = msadata[chr];
  auto msa = msatrees[chr];
  auto dbv = dbs[chr];
  std::cerr << "Trying to initiate MSA on " << chr << " with db counts "
            << dbv.size() << std::endl;
  class Visitor : public _DB::Visitor {
  public:
    decltype(msad) v;
    size_t cnt;
    short i;
    Visitor(decltype(msad) v_, short i_) : v(v_), cnt(0), i(i_){};
    const char *visit_full(const char *kbuf, size_t ksiz, const char *vbuf,
                           size_t vsiz, size_t *sp) {
      unsigned l = ((unsigned *)kbuf)[0];
      unsigned r = ((unsigned *)kbuf)[1];
      auto ii = inode(l, l + r, i);
      v->push_back(ii);
      cnt++;
      if (cnt % 5000 == 0)
        std::cerr << ".";
      return NOP;
    };
    const char *visit_empty(const char *kbuf, size_t ksiz, const char *vbuf,
                            size_t vsiz, size_t *sp) {
      return NOP;
    };
  };

  {
    auto dph = dbpaths[chr];
    std::ifstream f(dph + ".index", std::ifstream::in | std::ifstream::binary);
    size_t sz;
    unsigned l, r;
    short p;
    if (!f.good()) {
      std::cerr << "Error reading index file, reading from DB instead!"
                << std::endl;
      f.close();
      for (short i = 0; i < dbv.size(); ++i) {
        auto vis = Visitor(msad, i);
        dbv[i]->iterate(&vis, false);
        cnt += vis.cnt;
      }
      save_index(chr, *msad);
    } else {
      msad->reserve(sizes[chr]);
      cnt = 0;
      while (f.good()) {
        f.read((char *)&l, sizeof(unsigned));
        f.read((char *)&r, sizeof(unsigned));
        f.read((char *)&p, sizeof(short));
        msad->push_back(inode(l, l + r, p));
        if (cnt % 5000 == 0)
          std::cerr << ".";
        cnt++;
      }
      f.close();
    }
  }

  for (auto it = msad->begin(); it != msad->end(); ++it)
    msa->insert(*it);
  std::cerr << "Total inserted: " << cnt << " records" << std::endl;
#if DEBUG
  for (size_t ll = 1000; ll < msad->back().l; ll *= 2) {
    auto it = msa->lower_bound(inode(ll, ll + ll / 2));
    if (it != msa->end())
      std::cerr << "Trying to find a record >" << ll << " : " << it->l << ","
                << it->r << "@" << it->p << std::endl;
  }
#endif
  return;
}
void mafdb::init_tree() {
  if (init)
    return;
  // size_t sz_all=0;
  std::cerr << "Trying to initialize AMSet for MSA" << std::endl;
  for (auto &chr : chrs) {
    // auto sz= sizes[chr];
    // sz_all+=sz;
    msadata[chr] = std::shared_ptr<std::vector<inode>>(new std::vector<inode>);
    // msadata[chr]->reserve(sz);
    msatrees[chr] = std::shared_ptr<AMSet>(new AMSet);
    load_index(chr);
  }
  init = true;
  // std::cerr << "Total elements: " <<sz_all<< std::endl;
  return;
}
void mafdb::load_db_() {
  for (auto chr : chrs) {
    for (int i = 0; i < postfixes[chr]; ++i) {
      auto dbp = dbpaths[chr];
      if (i > 0)
        dbp = dbp + "." + std::to_string(i);
      auto db = std::shared_ptr<_DB>(new _DB);
#if USE_DBT
      db->tune_comparator(&CMPSZ);
#endif
      if (!db->open(dbp, _DB::OREADER))
        throw("open error (load db scaffold): " +
              std::string(db->error().name()) + " on " + dbp);
      dbs[chr].push_back(db);
    }
  }
  init_tree();
};
void mafdb::load_db(const std::string &fp) {
  std::cerr << "Trying to deserialize from " << fp << std::endl;
  std::ifstream ifs(fp);
  if (!ifs.is_open())
    throw("open error (load db): " + fp);
  IARCHIVE ar(ifs);
  ar >> BOOST_SERIALIZATION_NVP(*this);
  load_db_();
}
void mafdb::load_db_kch(const std::string &kdbname, const std::string &key) {
  std::cerr << "Trying to deserialize from " << kdbname << std::endl;
  std::string serial_str;

  _DB db;
  if (!db.open(kdbname, _DB::OREADER))
    throw("open error (load db kch): " + std::string(db.error().name()));
  db.get("mafdb " + key, &serial_str);

  boost::iostreams::basic_array_source<char> device(serial_str.data(),
                                                    serial_str.size());
  boost::iostreams::stream<boost::iostreams::basic_array_source<char>> s(
      device);
  BIARCHIVE ia(s);
  ia >> BOOST_SERIALIZATION_NVP(*this);
  load_db_();
};
void mafdb::export_db(const std::string &fp) {
  std::cerr << "Trying to serialize into " << fp << std::endl;
  std::ofstream ofs(fp);
  if (!ofs.is_open())
    throw("open error (export db): " + fp);
  OARCHIVE ar(ofs);
  ar << BOOST_SERIALIZATION_NVP(*this);
}
void mafdb::export_db_kch(const std::string &kdbname) {
  std::cerr << "Trying to serialize into " << kdbname << std::endl;
  std::string serial_str;
  boost::iostreams::back_insert_device<std::string> inserter(serial_str);
  boost::iostreams::stream<boost::iostreams::back_insert_device<std::string>> s(
      inserter);
  BOARCHIVE oa(s);
  oa << BOOST_SERIALIZATION_NVP(*this);
  s.flush();

  _DB db;
  if (!db.open(kdbname, _DB::OWRITER | _DB::OCREATE))
    throw("open error (export db kch): " + std::string(db.error().name()));
  db.set("mafdb " + name, serial_str);
}
std::string mafdb::get(const unsigned &l, const unsigned &r) {

#if DEBUG
  std::cerr << " Getting matches for " << l << " , " << r << std::endl;
  auto pp = get_intervals(l, r);
  for (; pp.first != pp.second; pp.first++) {
    auto hits = extract_intervals(*pp.first);
    auto filtered = filter_intervals(l, r, hits);
  }
  std::cerr << " last one at : " << pp.first->l << " , " << pp.first->r
            << std::endl;
#endif
  return "";
}
// get the content from index
INTERVAL_PAIR mafdb::extract_intervals(const inode &node) {
  INTERVAL_PAIR out;
  auto dbv = dbs[chr];
  if (node.p < 0 || node.p >= dbv.size())
    throw("Incorrect db part info!");
  auto db = dbv[node.p];
  unsigned rr[2];
  rr[0] = node.l;
  rr[1] = node.r - node.l;
  std::string key, val;
  key = std::string((char *)&rr, sizeof(rr));
  if (!db->get(key, &val)) {
    std::cerr << "No match for " << key << std::endl;
    return out;
  }
#if DEBUG
  std::cerr << " Fetched value " << val << " for " << node.p << " , " << node.l
            << " , " << node.r << std::endl;
#endif
  std::stringstream ss(val);
  std::string rf;
  float score;
  ss >> score;

  while (ss.good() && !ss.eof()) {
    interval iv;
    ss >> rf;
    auto bk = rf.find(".");
    if (bk == std::string::npos) {
      throw("Error breaking ref: " + rf);
    } else {
      iv.ref = rf.substr(0, bk);
      iv.chr = rf.substr(bk + 1);
    }
    ss >> iv.l >> iv.r >> iv.seq;
    iv.score = score;
    iv.r += iv.l;
    if (iv.seq[0] == '-') {
      iv.strand = false;
    }
    iv.seq = iv.seq.substr(1);
    if (iv.ref != ref) {
      iv.tt = matching;
      out.second[iv.ref] = iv;
    } else {
      out.first = iv;
    }
#if DEBUG
    std::cerr << "Fetched an inverval " << print_interval(iv) << std::endl;
#endif
  }
  return out;
}

// given a subinterval,
// find all intersecting actual intervals (no gaps) on the list of extracted
// intervals.
INTERVAL_PAIR mafdb::filter_intervals(const unsigned &l, const unsigned &r,
                                      INTERVAL_PAIR &hits) {
  int lshift = l > hits.first.l ? l - hits.first.l : 0;
  int dist = r > hits.first.r ? hits.first.r - l : r - l;
#if DEBUG
  std::cerr << "Lshift " << lshift << " , dist " << dist << std::endl;
#endif
  auto inner = [](const char &c, std::string &out, int &cnt, int &gap) {
    if (c != '-') {
      out.push_back(c);
      ++cnt;
      return true;
    } else
      ++gap;
    return false;
  };
  auto routine = [&](const std::string &seq) {

    std::string tmp("");
    int foo, bar;
    for (int idx = 0; idx < seq.size(); ++idx) {
      auto flag = inner(seq[idx], tmp, foo, bar);
    }
    return tmp;
  };
  auto counter = [&](const std::string &seq, int l = 0, int sz = 0) {
    int gap = 0;
    int cnt = 0;
    std::string out;
    out.reserve(sz);
    for (int i = l; i < seq.size(); ++i) {
      char c = seq[i];
      inner(c, out, cnt, gap);
      if (cnt >= sz)
        break;
    }
    return std::pair<int, std::string>(gap, out);
  };
  int gap1 = counter(hits.first.seq, 0, lshift).first;
  auto rs = counter(hits.first.seq, lshift + gap1, dist);
#if DEBUG
  std::cerr << "For bg " << l << ", " << r << ":" << hits.first.seq << " gap "
            << gap1 << " lshift " << lshift << " dist " << dist << " gap2 "
            << rs.first << " seq " << rs.second << std::endl;
#endif
  // now traverse lshift+gap for first aligned position

  INTERVAL_PAIR out;
  out.first.l = lshift;
  out.first.r = dist;
  out.first.seq = rs.second;
  out.first.chr = chr;
  out.first.ref = ref;
  out.first.strand = hits.first.strand;
  out.first.score = hits.first.score;
  hits.first.seq = routine(hits.first.seq);
  if (lshift + dist > hits.first.seq.size())
    dist = hits.first.seq.size() - lshift;
#if DEBUG
  std::cerr << "For ref: " << print_interval(out.first) << std::endl;
  ;
#endif
  for (auto &it : hits.second) {
    int t_lshift = 0, t_dist = 0, t_gap1 = 0, t_gap2 = 0;
    std::string seq1(""), seq2("");
    for (int i = 0; i < lshift + gap1; ++i) {
      char c = it.second.seq[i];
      inner(c, seq1, t_lshift, t_gap1);
    }
    for (int i = lshift + gap1; i < lshift + gap1 + dist + rs.first; ++i) {
      char c = it.second.seq[i];
      inner(c, seq2, t_dist, t_gap2);
    }
    it.second.seq = routine(it.second.seq);
    interval itt;
    itt.seq = seq2;
    itt.ref = it.second.ref;
    itt.chr = it.second.chr;
    itt.strand = it.second.strand;
    itt.score = it.second.score;
    itt.l = t_lshift;
    itt.r = t_dist; // local change
    out.second[itt.ref] = itt;
#if DEBUG
    std::cerr << " For matching seq " << print_interval(itt) << std::endl;
#endif
  }
  return out;
}
