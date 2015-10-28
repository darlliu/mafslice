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
        //boost::replace_last (fn, ".fasta","");
        chrs.push_back(fn);
    }
    init_db(chrs, dbs);
    //init_db(chrs, dbs2);
    for (unsigned i =0; i < fps.size(); ++i)
    {
        auto chr = chrs[i];
        set_chr(chr);
        auto chrp = fps[i];
        fapaths [chr] = chrp;
#if USE_DBT
        auto dbp = dbpath+"/"+chr+".MSA.kct";
        //auto dbp2 = dbpath+"/"+chr+".MSAinfo.kct";
#else
        auto dbp = dbpath+"/"+chr+".MSA.kch";
        //auto dbp2 = dbpath+"/"+chr+".MSAinfo.kch";
#endif
        dbpaths [chr] = dbp;
        //dbpaths2 [chr] = dbp2;
        import_chr();
    }
    return true;
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
    float score;
    size_t start, end, len, rstart, rend, rr[2];
    std::stringstream ss;
    if (!ifs.is_open()) throw("Error Opening file!");
    ifs.seekg (0, ifs.end);
    size_t fsize = ifs.tellg();
    ifs.seekg (0, ifs.beg);
    std::cout<< "Opened file "<<chrfp << " On chr: "<< \
        chr <<" Size: "<<fsize<<" Est. Recs: <"<<fsize/2000<<std::endl;
    fsize/=2000;
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
        if (!db->open(dbp2, _DB::OWRITER | _DB::OCREATE))
        {
            std::cerr << "open error: " << db->error().name() << std::endl;
            throw ("Error opening DB");
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
                    db->set_bulk(data,false);
                    total+=idx;
                    std::cerr << idx << "/"<<total<<" .";
                    data.clear();
                    if (idx > 1<<20)
                    {
                        dbv.push_back(std::shared_ptr<_DB>(new _DB));
                        db = dbv.back();
                        auto dbp = tune(db,dbv.size());
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
            } else {
                if (tmp.size()>0) tmp+="\t";
                tmp += species+" "+std::to_string(start) + " " + std::to_string(end)\
                       +" "+sign+seq;
            }
        } else continue;
    }
    {
        std::cerr << "@" ;
        db->set_bulk(data,false);
        total+=idx;
        std::cerr << idx << "/"<<total<<" .";
        data.clear();
    }
    ifs.close();
    std::cout<<std::endl << "Loaded length "<<idx<<std::endl;
    sizes[chr] = total;
    postfixes[chr] = dbv.size();

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
    size_t cnt = 0, l, r;
    std::string key, val;
    auto msad = msadata[chr];
    auto msa = msatrees[chr];
    auto dbv = dbs[chr];
    std::cerr << "Trying to initiate MSA on "<<chr<<" with db counts "<<dbv.size()<<std::endl;
    for (short i=0; i<dbv.size();++i)
    {
        auto cur = dbv[i]->cursor();
        cur->jump();
        size_t lr [2] ;
        while (cur->get(&key, &val, true)){
            l = ((size_t*)key.c_str())[0];
            r = ((size_t*)key.c_str())[1];
            auto ii = inode(l, l+r,i);
            msad->push_back(ii) ;
            cnt++;
            if (cnt % 5000==0) std::cerr<<".";
        } ;
        delete cur;
    }
    for (auto it = msad->begin(); it!=msad->end();++it) msa->insert(*it);
    std::cerr << "Finished inserting, last record:"<<l<<","<<r<<std::endl;
    std::cerr <<cnt<<" records"<<std::endl;
    for (size_t ll=1000; ll<l; ll*=2)
    {
        auto it = msa->lower_bound(inode(ll, ll+ll/2));
        if (it!=msa->end())
            std::cerr << "Trying to find a record >"<<ll<<" : "<<it->l<<","<<it->r<<"@"<<it->p<<std::endl;
    }
    return;
}
void mafdb::init_tree()
{
    if (init) return;
    size_t sz_all=0;
    std::cerr << "Trying to initialize AMSet for MSA" << std::endl;
    for (auto &chr: chrs){
        auto sz= sizes[chr];
        sz_all+=sz;
        msadata[chr]=std::shared_ptr <std::vector<inode>> (new std::vector<inode>);
        msadata[chr]->reserve(sz);
        msatrees[chr]=std::shared_ptr <AMSet> (new AMSet);
        load_index(chr);
    }
    init = true;
    std::cerr << "Total elements: " <<sz_all<< std::endl;
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
};
std::string mafdb::get(const size_t& l , const size_t& r)
{
}
std::string mafdb::get(const std::string& key)
{
}
// get the content from index

