#include "indexer.hpp"
#include <bitset>
char encode_char(char* in)
{
    std::bitset<8> bs;
    //std::cerr <<"Encoding "<<bs << " from "<<in[0]<<in[1]<<std::endl;
    for (int i =0; i<2; ++i)
        switch (in[i]){
            case 'a': {bs[i*4+1]=1; }
            case 'A': {bs[i*4+2]=0; bs[i*4+3]=0; break;}
            case 't': {bs[i*4+1]=1; }
            case 'T': {bs[i*4+2]=1; bs[i*4+3]=1; break;}
            case 'c': {bs[i*4+1]=1; }
            case 'C': {bs[i*4+2]=1; bs[i*4+3]=0; break;}
            case 'g': {bs[i*4+1]=1; }
            case 'G': {bs[i*4+2]=0; bs[i*4+3]=1; break;}
            default: bs[i*4]=1;
        }
    //std::cerr <<"Got "<<bs << " and "<<bs.to_ulong()<<std::endl;
    return static_cast<unsigned char> (bs.to_ulong());
};
void decode_char(const char& in,char* out)
{
    std::bitset<8> bs ((unsigned char)in);
    //std::cerr <<"Decoding "<<bs<<"...";
    for (int i =0; i<2; ++i)
    {
        if (bs[i*4]) {out[i]='N';
        } else if (bs[i*4+2]&&!bs[i*4+3]){
            if(bs[i*4+1]) out[i]='c';
            else out[i]='C';
        } else if (!bs[i*4+2]&&bs[i*4+3]){
            if(bs[i*4+1]) out[i]='g';
            else out[i]='G';
        } else if (!bs[i*4+2]&&!bs[i*4+3]) {
            if(bs[i*4+1]) out[i]='a';
            else out[i]='A';
        } else if (bs[i*4+2]&&bs[i*4+3]){
            if(bs[i*4+1]) out[i]='t';
            else out[i]='T';
        } else throw ("Not all cases covered in decoding");
        //std::cerr <<out[i];
    }
    //std::cerr <<std::endl;
    return;
};
void encode_seq(std::string in, char* out)
{
    //out.reserve(in.size()/2);
    char tmp[2];
    for (int i =0; i< (in.size()+1)/2; ++i)
    {
        tmp[0]=in[i*2];
        if (i*2+1< in.size())
            tmp[1]=in[i*2+1];
        char c = encode_char (tmp);
        out[i]=c;
    }
    return;
};
std::string decode_seq (char* in, unsigned size)
{
    std::string out ;
    out.reserve(size*2);
    for (int i=0; i<size;++i)
    {
        char tmp[2];
        decode_char (in[i], tmp);
        out.push_back(tmp[0]);
        out.push_back(tmp[1]);
    }
    return out;
}


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
#if USE_FSQ
    for (unsigned i =0; i < fps.size(); ++i)
    {
        auto chr = chrs[i];
        set_chr(chr);
        auto chrp = fps[i];
        fapaths [chr] = chrp;
        auto dbp = dbpath+"/"+chr+".fsq";
        dbpaths [chr] = dbp;
        import_chr_fsq();
    }
    load_db_();
#else
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
#endif
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
#elif USE_DBH
    auto dbp = dbpath+"/"+chr+".kch";
#elif USE_FSQ
    auto dbp = dbpath+"/"+chr+".fsq";
#endif
    dbpaths [chr] = dbp;
#if USE_FSQ
    import_chr_fsq();
    return;
#endif
    init_db(chrs, dbs);
    import_chr();
    return;
}
void seqdb::import_chr_fsq()
{
    auto chrfp = fapaths[chr];
    auto dbp = dbpaths[chr];
    std::ofstream dbf1 (dbp, std::ofstream::binary );
    boost::replace_last(dbp, ".fsq", ".szi");
    std::ofstream dbf2 (dbp, std::ofstream::binary );
    std::ifstream ifs (chrfp);

    if (!ifs.is_open()|| !dbf1.is_open() || !dbf2.is_open()) throw("Error Opening file import sequence!");
    size_t idx = 0;
    std::string tmp="",line;
    std::cout<< "Opened file "<<chrfp << " On chr: "<< chr <<std::endl;
    while (std::getline(ifs, line))
    {
        if (line.size()==0 || line[0] == '>') continue;
        tmp += line;
    }
    ifs.close();
    auto buf = new char[(tmp.size()+1)/2];
    encode_seq(tmp, buf);
    auto pos = dbf1.tellp();
    dbf1.write(buf, (tmp.size()+1)/2);
    if ((!assemble)&&!(dbf1.good()))
        throw( "Error writing to a 4bit sequence "+ chr);
    dbf1.close();
    idx+= tmp.size();
    sizes[chr]=idx;
    shifts [chr]=pos;
    dbf2<<(chr+" "+std::to_string(pos)+" "+std::to_string(idx)+" ");
    std::cout << "Loaded length "<<idx<<std::endl;
    dbf2.close();
    delete [] buf;
    return;
};
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
#if USE_FSQ

void seqdb::import_feed()
{
    std::cerr <<"..... Now importing from stdin ....."<<std::endl;
    std::string line, tmp(""),dbp;
    std::ofstream dbf1, dbf2;
    if (scaffold)
    {
        dbp = dbpath+"/"+name+".fsq";
        dbf1.open(dbp, std::ofstream::binary);
        dbpaths[name]=dbp;
        boost::replace_last(dbp, ".fsq", ".szi");
        dbf2.open(dbp, std::ofstream::binary);
    }
    size_t idx = 0;
    bool flag= false;
    auto inner = [&]()
    {
        auto buf = new char[(tmp.size()+1)/2];
        encode_seq(tmp,buf);
        auto pos = dbf1.tellp();
        dbf1.write(buf, (tmp.size()+1)/2);
        idx+= tmp.size();
        dbf2 << (chr+" "+std::to_string(pos)+" "+std::to_string(idx)+" ");
        std::cout << "Loaded length "<<idx<<std::endl;
        idx = 0;
        dbf2.flush();
        tmp="";
        delete[] buf;
    };
    while (std::cin)
    {
        std::getline(std::cin, line);
        if (line.size()==0) continue;
        if (line[0]=='>')
        {
            if (flag)
                inner();
            chr = line.substr(1);
            std::cerr<<"Found breaking point "<<chr <<" . Assuming this is a chromosome-like!"<<std::endl;
            if (!scaffold)
            {
                dbf1.close();
                dbf2.close();
                dbp = dbpath+"/"+chr+".fsq";
                dbpaths[chr]=dbp;
                dbf1.open(dbp, std::ofstream::binary);
                boost::replace_last(dbp, ".fsq", ".szi");
                dbf2.open(dbp, std::ofstream::binary);
            }
            flag=true;
            continue;
        }
        tmp += line;
    }
    inner();
    dbf1.close();
    dbf2.close();
    load_db_();
    std::cerr << "Finished importing"<<std::endl;
    return;
}
#else
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
#endif

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
#if USE_FSQ
    if (r<l)
        throw("Interval incorrect!");
    std::shared_ptr<std::ifstream> db;
    unsigned pos=0;
    //std::cerr <<"Trying to get "<<scaffold <<" "<<l <<", "<<r<<", "<<pos <<", "<<sizes[chr]<<std::endl;
    if (scaffold)
    {
        db = dbs_fsq[name];
        pos = shifts[chr];
    }
    else
        db = dbs_fsq[chr];
    unsigned ll, rr;
    if (l%2) ll=l-1;
    else ll=l;
    if (r%2) rr=r+1;
    else rr=r;
    unsigned rr2=rr/2, ll2=ll/2;
    db->seekg(0, db->beg);//seek back
    db->seekg(pos+ll2,db->beg);
    auto buf = new char [rr2-ll2];
    //std::cerr <<"Reading file for get "<<ll2<<", "<<rr2<<" from "<<pos<<"total"<<rr2-ll2<<std::endl;
    db->read(buf, rr2-ll2);
    db->sync();
    //std::cerr <<"Decoding string for get count "<<db->gcount()<<std::endl;
    std::string decoded_seq = decode_seq(buf, rr2-ll2);
    delete[] buf;
    //std::cerr <<"returning substr from" <<decoded_seq<<std::endl;
    return decoded_seq.substr(l-ll, r-l);

#else
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
#endif
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
void seqdb::load_sizes(std::ifstream& ifs)
{
    std::string ch;
    unsigned pos;
    std::stringstream ss;
    ss<< ifs.rdbuf();
    while (ss.good())
    {
        ss >> ch;
        ss >> pos;
        shifts[ch]=pos;
        ss >> pos;
        sizes [ch]=pos;
        chrs.push_back(ch);
    }
    ifs.close();
    return;
}

void seqdb::load_db_()
{
#if USE_FSQ
    if (scaffold)
    {
        dbs_fsq[name] = std::shared_ptr<std::ifstream>
            ( new std::ifstream(dbpath+"/"+name+".fsq",
                                std::ifstream::in | std::ifstream::binary));
        dbs_fsz[name] = std::shared_ptr<std::ifstream>
            ( new std::ifstream(dbpath+"/"+name+".szi"));
        if (!dbs_fsq[name]->good() || !dbs_fsz[name]->good())
            throw ("Error opening a 4bit file at "+name);
        load_sizes(*dbs_fsz[name]);
    } else {
        for (auto it:dbpaths)
        {
            auto dbp = it.second;
            dbs_fsq[it.first] = std::shared_ptr <std::ifstream>(new std::ifstream(dbp));
            boost::replace_last(dbp, ".fsq",".szi");
            dbs_fsz[it.first] = std::shared_ptr <std::ifstream>(new std::ifstream(dbp));
            if (!dbs_fsq[it.first]->good() || !dbs_fsz[it.first]->good())
                throw ("Error opening a 4bit file at "+it.first);
            load_sizes(*dbs_fsz[it.first]);
        }
    }
    return;
#else
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
#endif
};

void seqdb::load_db(const std::string& fp )
{
    std::cerr << "Trying to deserialize from "<<fp <<" ... ";
    std::ifstream ifs (fp);
    if (!ifs.is_open())
        throw("open error (load db): " + fp);
    IARCHIVE ar(ifs);
    ar >> BOOST_SERIALIZATION_NVP(*this);
    std::cerr << "Loading SeqDB ... "<<std::endl;
    load_db_();
};

void seqdb::load_db_kch(const std::string& kdbname, const std::string& key )
{
    std::cerr << "Trying to deserialize from "<<kdbname <<" ... ";
    std::string serial_str;

    _DB db;
    if (!db.open(kdbname, _DB::OREADER))
        throw("open error (load db kch): " + std::string(db.error().name()));
    db.get("seqdb "+key,&serial_str);

    boost::iostreams::basic_array_source<char> device(serial_str.data(), serial_str.size());
    boost::iostreams::stream<boost::iostreams::basic_array_source<char> > s(device);
    IARCHIVE ia(s);
    ia >> BOOST_SERIALIZATION_NVP(*this);
    std::cerr << "Loading SeqDB ... "<<std::endl;
    load_db_();
};
