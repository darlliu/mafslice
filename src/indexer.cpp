#include "indexer.hpp"

std::string print_interval(const interval& in)
{
    auto fm = boost::format("[%1%, %2%]: (%3%, %4%) @ %7% @ {score: %5% , strand: %6%}");
    return (fm % in.ref % in.chr % in.l % in.r % in.score % (int)in.strand % in.seq).str();
}

std::string get_reverse_comp(const std::string& in)
{
    std::string out;
    out.reserve(in.size());
    for (auto &c:in)
    {
        switch(c){
            case 'A': out.push_back('T');
            case 'T': out.push_back('A');
            case 'U': out.push_back('A');
            case 'G': out.push_back('C');
            case 'C': out.push_back('G');
            case 'a': out.push_back('t');
            case 't': out.push_back('a');
            case 'u': out.push_back('a');
            case 'g': out.push_back('c');
            case 'c': out.push_back('g');
            default: out.push_back(c);
        }
    }
};

void seqdb::import (const std::string& dirname)
{
    using namespace boost::filesystem;
    auto fp = path (dirname);
    std::vector<std::string> fps, fns, chrs;
    try
    {
        if (exists(fp) && is_directory(fp))
        {
            std::cerr << " Now reading " << dirname << std::endl;
        }
        else
            throw ("Path does not exist: "+dirname);
    }
    catch (const filesystem_error& ex)
    {
        std::cerr<< "Error opening path " << dirname <<std::endl;
        throw (ex);
    }
    for (auto it = directory_iterator(fp); it!= directory_iterator(); ++it)
    {

        if (!is_directory(it->path()))
        {
            auto fp = it->path().string();
            if (fp.find(".fa")==std::string::npos) continue;
            std::cout << "Found fasta file : " << fp << std::endl;
            fps.push_back(fp);
            fns.push_back(it->path().filename().string());

        }
        else
        {
            std::cout << "skipping directory : " << it->path().filename()<< std::endl;
        }
    }
    std::cout << "Found files count: "<<fps.size()<<std::endl;
    for (std::string fn: fns)
    {
        boost::replace_last (fn, ".fa","");
        boost::replace_last (fn, ".fasta","");
        chrs.push_back(fn);
    }
    init_db(chrs, dbs);
    for (unsigned i =0; i < fps.size(); ++i)
    {
        auto chr = chrs[i];
        set_chr(chr);
        auto chrp = fps[i];
        fapaths [chr] = chrp;
        auto dbp = dbpath+"/"+chr+".kch";
        dbpaths [chr] = dbp;
        import_chr();
    }
};
void seqdb::import_chr(const std::string& fname)
{
    using namespace boost::filesystem;
    auto fp=path(fname);
    if (exists(fp) && fname.find(".fa")!=std::string::npos)
    {
        std::cerr << " Now reading " << fname << std::endl;
    }
    else
        throw ("File does not exist: "+fname);
    auto fn = fp.filename().string();
    fapaths[chr]=fp.string();
    boost::replace_last (fn, ".fa","");
    boost::replace_last (fn, ".fasta","");
    chrs.push_back(fn);
    set_chr(fn);
    fapaths[chr]=fp.string();
#if USE_DBT
    auto dbp = dbpath+"/"+chr+".kct";
#else
    auto dbp = dbpath+"/"+chr+".kch";
#endif
    dbpaths [chr] = dbp;
    init_db(chrs, dbs);
    import_chr();
}

void seqdb::import_chr()
{
    auto chrfp = fapaths[chr];
    auto dbp = dbpaths[chr];
    auto db = dbs[chr][0];
    if (assemble){
        if (!db->open(dbp, _DB::OREADER))
            throw("open error (assemble): " + std::string(db->error().name()));

    } else {
        if (!db->open(dbp, _DB::OWRITER | _DB::OCREATE))
            throw("open error (create): " + std::string(db->error().name()));

    }
    size_t idx = 0;
    std::ifstream ifs (chrfp);
    std::string tmp="",line;
    if (!ifs.is_open()) throw("Error Opening file!");
    std::cout<< "Opened file "<<chrfp << " On chr: "<< chr <<std::endl;
    while (std::getline(ifs, line))
    {
        if (line.size()==0 || line[0] == '>') continue;
        tmp += line;
        if (tmp.size()>chunksz)
        {
            if ((!assemble)&&!db->set(std::to_string(idx), tmp.substr(0, chunksz)))
                throw( "Error adding a value to db "+ chr);
            idx+=chunksz;
            tmp = tmp.substr(chunksz);
        }
    }
    ifs.close();
    if (!assemble) db->set(std::to_string(idx), tmp);
    idx+= tmp.size();
    std::cout << "Loaded length "<<idx<<std::endl;
    sizes[chr] = idx;
    return;
};

void seqdb::import_feed()
{
    std::cerr <<"..... Now importing from stdin ....."<<std::endl;
    std::string line, tmp(""),dbp;
    _DB db, db_sz;
    auto tune_and_open = [&](_DB& db, const std::string& name)
    {
        db.tune_buckets (1LL* 100000 * 2);
        db.tune_options(kyotocabinet::HashDB::TLINEAR);
        db.tune_map(2LL << 30);
        db.tune_defrag(8);
        if (!db.open(dbpath+"/"+name+".kch", _DB::OWRITER | _DB::OCREATE))
            throw("open error (feedin): " + std::string(db.error().name()));
    };
    if (scaffold)
    {
        std::cerr <<"This is a scaffold db, tuning to 100k size automatically."<<std::endl;
        dbp= dbpath+"/"+name+".kch";
        tune_and_open(db,name);
        tune_and_open(db_sz, name+"_sizes");
    }
    size_t idx = 0;
    bool flag= false;
    auto inner = [&](const std::string& key){
        if (flag)
        {
            db.set(key, tmp);
            idx+= tmp.size();
            std::cout << "Loaded length "<<idx<<std::endl;
            if (scaffold)
            {
                db_sz.set(chr,std::to_string(idx));
            }
            else
                sizes[chr] = idx;
            idx=0;
            tmp="";
            flag=false;
        }
    };
    while (std::cin)
    {
        std::getline(std::cin, line);
        if (line.size()==0) continue;
        if (line[0]=='>')
        {
            if (scaffold)
            {
                inner(chr+std::to_string(idx));
            }
            else
            {
                dbp= dbpath+"/"+chr+".kch";
                inner(std::to_string(idx));
                if (flag) db.close();
                if (!db.open(dbp, _DB::OWRITER | _DB::OCREATE))
                {
                    throw("open error (feedin_create): " + std::string(db.error().name()));
                }
                dbpaths[chr]=dbp;
            }
            chr = line.substr(1);
            chrs.push_back(chr);
            std::cerr<<"Found breaking point "<<chr <<" . Assuming this is a chromosome-like!"<<std::endl;
            flag=true;
            continue;
        }
        tmp += line;
        if (tmp.size()>chunksz)
        {
            if (scaffold)
            {
                if (!db.set(chr+std::to_string(idx), tmp.substr(0, chunksz)))
                    throw( "Error adding a value to db "+ chr);
            } else {
                if (!db.set(std::to_string(idx), tmp.substr(0, chunksz)))
                    throw( "Error adding a value to db "+ chr);
            }
            idx+=chunksz;
            tmp = tmp.substr(chunksz);
        }

    }
    if (scaffold)
    {
        inner(chr+std::to_string(idx));
        sizes.clear();
        chrs.clear();
    }
    else
        inner(std::to_string(idx));
    return;

};

void seqdb::init_db (const std::vector<std::string>& chrs, DB& dbs)
{
    if (dbs.size()>0) close_db(dbs);
    std::cout << "Initializing DBs, length: " <<chrs.size() <<std::endl;
    for (auto chr: chrs)
    {
        dbs[chr].push_back(std::shared_ptr <_DB>(new _DB));
    }
    return;
};

std::string seqdb::get(const size_t& l, const size_t& r)
{
    //First, check the sizes and make sure that the
    if (!scaffold && (r<=l || r >= sizes[chr]))
    {
        std::cerr << "Index for slice is incorrect"<< std::endl;
        return "";
    }
    auto idx = get_index(l); //integer division on chunksz
    auto idx0 = idx;
    std::string val(""), tmp, key("");
    std::shared_ptr<_DB> db, db_sz;
    if (scaffold)
    {
        db=dbs[name][0];
        db_sz=dbs[name][1];
    }
    else
        db = dbs[chr][0];
    do
    {
        if (scaffold)
            key = chr+std::to_string(idx);
        else
            key = std::to_string(idx);
        if (!db->get(key,&tmp))
        {
            std::cerr << "Get Error : " << db->error().name() \
                << " On : "<< chr << std::to_string(idx) << std::endl;
            return "";
        }
        val += tmp;
        idx += chunksz;
    }
    while (r > idx);
    int ul =val.size();
    if (ul>(r-idx0)) return val.substr(l-idx0, r-l);;
    if (ul>(l-idx0)) return val.substr(l-idx0);
    return "";
};
void seqdb::export_db_kch(const std::string& kdbname)
{
    std::cerr << "Trying to serialize into "<< kdbname <<std::endl;
    std::string serial_str;
    boost::iostreams::back_insert_device<std::string> inserter(serial_str);
    boost::iostreams::stream<boost::iostreams::back_insert_device<std::string> > s(inserter);
    OARCHIVE oa (s);
    oa << BOOST_SERIALIZATION_NVP(*this);
    s.flush();

    _DB db;
    if (!db.open(kdbname, _DB::OWRITER|_DB::OCREATE))
        throw("open error (export db kch): " + std::string(db.error().name()));
    db.set("seqdb "+name,serial_str);
}

void seqdb::export_db(const std::string& fp )
{
    std::cerr << "Trying to serialize into "<<fp <<std::endl;
    std::ofstream ofs (fp);
    if (!ofs.is_open())
        throw("open error (export db): " + fp);
    OARCHIVE ar(ofs);
    ar << BOOST_SERIALIZATION_NVP(*this);
};

void seqdb::load_db_()
{
    if (scaffold)
    {
        auto db = std::shared_ptr <_DB>(new _DB);
        if (!db->open(dbpath+"/"+name+".kch", _DB::OREADER))
            throw("open error (load db scaffold): " + std::string(db->error().name()));
        dbs[name].push_back(db);

        auto db_sz = std::shared_ptr <_DB>(new _DB);
        if (!db_sz->open(dbpath+"/"+name+"_sizes.kch", _DB::OREADER))
            throw("open error (load db scaffold sz): " + std::string(db_sz->error().name()));
        dbs[name].push_back(db_sz);
    } else {
        for (auto it:dbpaths)
        {
            auto db = std::shared_ptr <_DB>(new _DB);
            if (!db->open(it.second, _DB::OREADER))
                throw("open error (load db increment): " + std::string(db->error().name()));
            dbs[it.first].push_back(db);
        }
    }
};

void seqdb::load_db(const std::string& fp )
{
    std::cerr << "Trying to deserialize from "<<fp <<std::endl;
    std::ifstream ifs (fp);
    if (!ifs.is_open())
        throw("open error (load db): " + fp);
    IARCHIVE ar(ifs);
    ar >> BOOST_SERIALIZATION_NVP(*this);
    load_db_();
};

void seqdb::load_db_kch(const std::string& kdbname, const std::string& key )
{
    std::cerr << "Trying to deserialize from "<< kdbname <<std::endl;
    std::string serial_str;

    _DB db;
    if (!db.open(kdbname, _DB::OREADER))
        throw("open error (load db kch): " + std::string(db.error().name()));
    db.get("seqdb "+key,&serial_str);

    boost::iostreams::basic_array_source<char> device(serial_str.data(), serial_str.size());
    boost::iostreams::stream<boost::iostreams::basic_array_source<char> > s(device);
    IARCHIVE ia(s);
    ia >> BOOST_SERIALIZATION_NVP(*this);
    load_db_();
};
