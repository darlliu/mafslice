#include"mafindexer.hpp"
#include<sstream>


bool mafdb::import (const std::string& dirname)
{
    using namespace boost::filesystem;
    auto fp = path (dirname);
    std::vector<std::string> fps, fns;
    try
    {
        if (exists(fp) && is_directory(fp))
        {
            std::cerr << " Now reading " << dirname << std::endl;
        }
        else
        {
            std::cerr << " Path does not exist: " << dirname << std::endl;
            return false ;
        }
    }
    catch (const filesystem_error& ex)
    {
        std::cerr<< "Error opening path " << dirname <<std::endl;
        return false;
    }
    for (auto it = directory_iterator(fp); it!= directory_iterator(); ++it)
    {

        if (!is_directory(it->path()))
        {
            auto fp = it->path().string();
            if (fp.find(".maf")==std::string::npos) continue;
            std::cout << "Found MAF file : " << fp << std::endl;
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
        boost::replace_last (fn, ".maf","");
        chrs.push_back(fn);
    }
    init_db(chrs, dbs);
    for (unsigned i =0; i < fps.size(); ++i)
    {
        auto chr = chrs[i];
        set_chr(chr);
        auto chrp = fps[i];
        fapaths [chr] = chrp;
#if USE_DBT
        auto dbp = dbpath+"/"+chr+".MSA.kct";
#else
        auto dbp = dbpath+"/"+chr+".MSA.kch";
#endif
        dbpaths [chr] = dbp;
        //dbpaths2 [chr] = dbp2;
        import_chr();
    }
    return true;
}
void mafdb::import_chr(const std::string& fname)
{
    using namespace boost::filesystem;
    auto fp=path(fname);
    if (exists(fp) && fname.find(".maf")!=std::string::npos)
    {
        std::cerr << " Now reading " << fname << std::endl;
    }
    else
    {
        std::cerr << "File does not exist: " <<fname<<std::endl;
        return;
    }
    auto fn = fp.filename().string();
    boost::replace_last (fn, ".maf","");
    chrs.push_back(fn);
    set_chr(fn);
    fapaths[chr]=fp.string();
#if USE_DBT
        auto dbp = dbpath+"/"+chr+".MSA.kct";
#else
        auto dbp = dbpath+"/"+chr+".MSA.kch";
#endif
    dbpaths [chr] = dbp;
    init_db(chrs, dbs);
    import_chr();
}
void mafdb::import_chr ()
{
    auto chrfp = fapaths[chr];
    auto dbv = dbs[chr];
    size_t idx = 0, total=0;
    std::ifstream ifs (chrfp);
    std::string tmp="",line, foo, sign, seq, species;
    std::string key, val;
    std::map<std::string, std::string> data;
    std::vector <inode> inodes;
    float score;
    unsigned start, end, len, rstart, rend, rr[2];
    std::stringstream ss;
    if (!ifs.is_open()) throw("Error Opening file!");
    ifs.seekg (0, ifs.end);
    size_t fsize = ifs.tellg();
    ifs.seekg (0, ifs.beg);
    std::cout<< "Opened file "<<chrfp << " On chr: "<< \
        chr <<" Size: "<<fsize<<" Est. Recs: <"<<fsize/2000<<std::endl;
    fsize/=2000;
    inodes.reserve(fsize);
    if (fsize > 1<<20 )
    {
        std::cout<< "Size is too big, trimming into "<< (1<<20)<<std::endl;
        fsize=1<<20;
    }
    auto tune = [&](decltype(dbv[0]) &db, short pfx=0){
        auto dbp = dbpaths[chr];
        std::cerr<<std::endl;
#if USE_DBT
        if (fsize>65536)
        {
            std::cerr<<" Tuning new bucket size (B+ Tree): "<<fsize/10 <<std::endl;
            db->tune_buckets (1LL* fsize/10);
            db->tune_map(2LL << 30);
            db->tune_defrag(8);
            std::cerr<<" Tuning comparator..."<<std::endl;
            db->tune_comparator(&CMPSZ);
        }
#else
        if (fsize*4>1000*1000)
        {
            std::cerr<<" Tuning new bucket size (Hash table): "<<fsize*2 <<std::endl;
            db->tune_buckets (1LL* fsize*2);
            db->tune_options(kyotocabinet::HashDB::TLINEAR);
            db->tune_map(2LL << 30);
            db->tune_defrag(8);
        }
#endif
        auto dbp2=dbp;
        if (pfx>0) dbp2 = dbp+"."+std::to_string(pfx-1);
        if (assemble){
            if (!db->open(dbp2, _DB::OREADER))
            {
                throw ("Error opening DB: "+std::string(db->error().name()));
            }

        } else {
            if (!db->open(dbp2, _DB::OWRITER | _DB::OCREATE))
            {
                throw ("Error opening DB: "+std::string(db->error().name()));
            }
        }
        return dbp2;
    };
    auto db = dbv[0];
    tune(db);
    while (std::getline(ifs, line))
    {
        if (line.size()==0 || line[0]=='#') continue;
        ss.clear();
        ss.str(line);
        if (line[0]=='a'){
            if (tmp.size()>0){
                //auto hs = hasher (rstart, rend);
                rr[0]=rstart;
                rr[1]=rend;
                key = std::string((char*)&rr, sizeof(rr));
                val = std::to_string(score)+" "+tmp;
                data[key]=val;
                idx++;
                if (idx % 5000 ==0)
                {
                    std::cerr << "@" ;
                    if(!assemble) db->set_bulk(data,false);
                    std::cerr << idx << "/"<<total<<" .";
                    data.clear();
                    if (idx > 1<<20)
                    {
                        dbv.push_back(std::shared_ptr<_DB>(new _DB));
                        db = dbv.back();
                        auto dbp = tune(db,dbv.size());
                        total+=idx;
                        std::cerr<<"Creating a new db file: "<<dbp<<std::endl;
                        idx=0;
                    }
                }
                //saving routine in here
            }
            tmp = "";
            boost::replace_last (line, "="," ");
            ss >> foo >> foo >> score;
        } else if (line[0]=='i') {
            continue;
        } else if (line[0]=='s') {
            ss >> foo >> species >> start >> end >> sign >> len >> seq;
            if (species.find(ref)!=std::string::npos){
                rstart = start;
                rend = end;
                inodes.push_back (inode (rstart, rend, dbv.size()-1));
            }
            if (tmp.size()>0) tmp+="\t";
            tmp += species+" "+std::to_string(start) + " " + std::to_string(end)\
                   +" "+sign+seq;
        } else continue;
    }
    {
        std::cerr << "@" ;
        if (!assemble) db->set_bulk(data,false);
        total+=idx;
        std::cerr << idx << "/"<<total<<" .";
        data.clear();
    }
    ifs.close();
    std::cout<<std::endl << "Loaded length "<<idx << " Saving indexes..."<<std::endl;
    save_index(chr, inodes);
    sizes[chr] = total;
    postfixes[chr] = dbv.size();
    inodes.clear();

    return;
}
void mafdb::save_index(const std::string& chr, const std::vector<inode>& v)
{
    std::ofstream f(dbpaths[chr]+".index", std::ofstream::out|std::ofstream::binary);
    if (!f.good())
        throw ("Error opening the index file!");
    auto sz =v.size();
    f.write((char*) &sz, sizeof(size_t));
    for (auto &it:v)
    {
        f.write((char*) &it.l, sizeof(unsigned));
        f.write((char*) &it.r, sizeof(unsigned));
        f.write((char*) &it.p, sizeof(short));
    }
    f.close();
    return;
}
void mafdb::clear_index(const std::string& chr)
{
    auto msa = msatrees[chr];
    msa->clear();
    return;
}
void mafdb::load_index(const std::string& chr)
{
    size_t cnt = 0;
    std::string key, val;
    auto msad = msadata[chr];
    auto msa = msatrees[chr];
    auto dbv = dbs[chr];
    std::cerr << "Trying to initiate MSA on "<<chr<<" with db counts "<<dbv.size()<<std::endl;
    class Visitor : public _DB::Visitor{
        public:
        decltype(msad) v;
        size_t cnt;
        short i;
        Visitor(decltype(msad) v_,short i_):
            v(v_), cnt(0), i(i_){};
        const char* visit_full (const char* kbuf, size_t ksiz,
                                const char* vbuf, size_t vsiz, size_t *sp)
        {
            unsigned l = ((unsigned*)kbuf)[0];
            unsigned r = ((unsigned*)kbuf)[1];
            auto ii = inode(l, l+r,i);
            v->push_back(ii) ;
            cnt++;
            if (cnt % 5000==0) std::cerr<<".";
            return NOP;
        };
        const char* visit_empty (const char* kbuf, size_t ksiz,
                                const char* vbuf, size_t vsiz, size_t *sp)
        {
            return NOP;
        };
    };

    {
        auto dph = dbpaths[chr];
        std::ifstream f(dph+".index", std::ifstream::in | std::ifstream::binary);
        size_t sz;
        unsigned l, r;
        short p;
        if (!f.good())
        {
            std::cerr << "Error reading index file, reading from DB instead!"<<std::endl;
            f.close();
            for (short i=0; i<dbv.size();++i)
            {
                auto vis = Visitor(msad, i);
                dbv[i]->iterate(&vis, false);
                cnt+=vis.cnt;
            }
            save_index(chr,*msad);
        } else {
            f.read((char*) &sz, sizeof(size_t));
            msad->reserve(sz);
            for (cnt=0; cnt<sz; ++cnt)
            {
                f.read((char*) &l, sizeof(unsigned));
                f.read((char*) &r, sizeof(unsigned));
                f.read((char*) &p, sizeof(short));
                msad->push_back(inode(l,l+r,p));
                if (cnt%5000 ==0)
                    std::cerr << ".";
            }
            f.close();
        }
    }

    for (auto it = msad->begin(); it!=msad->end();++it)
        msa->insert(*it);
    std::cerr<<"Total inserted: " <<cnt<<" records"<<std::endl;
#if DEBUG
    for (size_t ll=1000; ll<msad->back().l; ll*=2)
    {
        auto it = msa->lower_bound(inode(ll, ll+ll/2));
        if (it!=msa->end())
            std::cerr << "Trying to find a record >"<<ll<<" : "<<it->l<<","<<it->r<<"@"<<it->p<<std::endl;
    }
#endif
    return;
}
void mafdb::init_tree()
{
    if (init) return;
    //size_t sz_all=0;
    std::cerr << "Trying to initialize AMSet for MSA" << std::endl;
    for (auto &chr: chrs){
        //auto sz= sizes[chr];
        //sz_all+=sz;
        msadata[chr]=std::shared_ptr <std::vector<inode>> (new std::vector<inode>);
        //msadata[chr]->reserve(sz);
        msatrees[chr]=std::shared_ptr <AMSet> (new AMSet);
        load_index(chr);
    }
    init = true;
    //std::cerr << "Total elements: " <<sz_all<< std::endl;
    return;
}
bool mafdb::load_db (const std::string & fp)
{
    std::cerr << "Trying to deserialize from "<<fp <<std::endl;
    std::ifstream ifs (fp);
    if (!ifs.is_open())
    {
        std::cerr << "Error opening file: "<<fp <<std::endl;
        return false;
    }
    IARCHIVE ar(ifs);
    ar >> BOOST_SERIALIZATION_NVP(*this);
    for (auto chr:chrs)
    {
        for (int i =0; i<postfixes[chr]; ++i)
        {
            auto dbp=dbpaths[chr];
            if (i>0) dbp= dbp+"."+std::to_string(i);
            auto db = std::shared_ptr <_DB>(new _DB);
#if USE_DBT
            db->tune_comparator(&CMPSZ);
#endif
            if (!db->open(dbp, _DB::OREADER))
            {
                std::cerr << "open error: " << db->error().name()<< \
                    " on: "<<dbp<< std::endl;
                return false;
            }
            dbs[chr].push_back(db);
        }
    }
    init_tree();
    return true;
}
bool mafdb::export_db(const std::string& fp )
{
    std::cerr << "Trying to serialize into "<<fp <<std::endl;
    std::ofstream ofs (fp);
    if (!ofs.is_open())
    {
        std::cerr << "Error opening file! "<<std::endl;
        return false;
    }
    OARCHIVE ar(ofs);
    ar << BOOST_SERIALIZATION_NVP(*this);
    return true;
}
std::string mafdb::get(const unsigned& l , const unsigned& r)
{

#if DEBUG
    std::cerr <<" Getting matches for "<<l <<" , "<<r <<std::endl;
    auto pp = get_interval(l,r);

    std::cerr<< " last one at : "<< pp.first->l << " , " << pp.first->r <<std::endl;
#endif
    return "";
}
// get the content from index

