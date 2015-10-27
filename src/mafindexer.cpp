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
    init_db(chrs, dbs2);
    for (unsigned i =0; i < fps.size(); ++i)
    {
        auto chr = chrs[i];
        set_chr(chr);
        auto chrp = fps[i];
        fapaths [chr] = chrp;
        auto dbp = dbpath+"/"+chr+".MSAseq.kch";
        auto dbp2 = dbpath+"/"+chr+".MSAinfo.kch";
        dbpaths [chr] = dbp;
        dbpaths2 [chr] = dbp2;
        import_chr();
    }
    return true;
}
void mafdb::import_chr ()
{
    auto chrfp = fapaths[chr];
    auto dbp = dbpaths[chr];
    auto dbp2 = dbpaths2[chr];
    auto db = dbs[chr];
    auto db2 = dbs2[chr];
    if (!db->open(dbp, kyotocabinet::HashDB::OWRITER | kyotocabinet::HashDB::OCREATE))
    {
        std::cerr << "open error: " << db->error().name() << std::endl;
        throw ("Error opening DB");
    }
    if (!db2->open(dbp2, kyotocabinet::HashDB::OWRITER | kyotocabinet::HashDB::OCREATE))
    {
        std::cerr << "open error: " << db2->error().name() << std::endl;
        throw ("Error opening DB");
    }
    size_t idx = 0;
    std::ifstream ifs (chrfp);
    std::string tmp="",line, foo, sign, seq, species, rseq;
    float score;
    size_t start, end, len, rstart, rend;

    std::stringstream ss;

    if (!ifs.is_open()) throw("Error Opening file!");
    std::cout<< "Opened file "<<chrfp << " On chr: "<< chr <<std::endl;
    while (std::getline(ifs, line))
    {
        if (line.size()==0 || line[0]=='#') continue;
        if (line[0]=='a'){
            if (tmp.size()>0){
                db2->set(std::to_string(rstart)+","+std::to_string(rend),
                        std::to_string(score)+":"+tmp);
                //if (idx % 50 == 0){
                    //std::cerr << "Storing index info "<<rstart<<","<<rend<<":"<<tmp<<std::endl;
                //}
                idx++;
            }
            //saving routine in here
            tmp = "";
            boost::replace_last (line, "="," ");
            ss.str(std::string());
            ss.clear();
            ss << line;
            ss >> foo;
            ss >> foo;
            ss >> score;
            //if (idx % 50 == 0){
                //std::cerr <<"Line is " << line << std::endl;
                //std::cerr << "Storing score info "<<foo<<":"<<score<<std::endl;
            //}
        } else if (line[0]=='i') {
            continue;
        } else if (line[0]=='s') {
            ss.str(std::string());
            ss.clear();
            ss << line;
            ss >> foo;
            ss >> species;
            ss >> start;
            ss >> end;
            ss >> sign;
            ss >> len;
            ss >> seq;
            if (species.find(ref)!=std::string::npos){
                rstart = start;
                rend = end;
                rseq = seq;
                //auto hs = hasher (species, rstart, rend);
                auto hs = species + std::to_string(start)+std::to_string(end);
                db->set(hs, sign+seq);
            } else {
                if (tmp.size()>0) tmp+="|";
                tmp += species+","+std::to_string(start) + "," + std::to_string(end);
                //end += start;
                //auto hs = hasher (species, start, end);
                auto hs = species + std::to_string(start)+std::to_string(end);
                db->set(hs, sign+seq);
                //if (idx % 50 == 0){
                    //std::cerr <<"Line is " << line << std::endl;
                    //std::cerr << "Storing seq info "<<foo<<":"<<species<<":"<<start<<":"<<end<<","<<sign<<","<<len<<":"<<seq<<std::endl;
                //}
            }
        } else continue;
    }
    ifs.close();
    std::cout << "Loaded length "<<idx<<std::endl;
    sizes[chr] = idx;
    return;
}
void mafdb::clear_index(const std::string& chr)
{
    auto msa = msatrees[chr];
    msa->clear();
    //for (auto it=msa->begin(); it!=msa->end(); ++it)
    //{
        //msa->erase(it);
    //}
    return;
}
void mafdb::load_index(const std::string& chr)
{
    std::cerr << "Trying to initiate MSA on "<<chr<<std::endl;
    size_t cnt = 0, l, r;
    std::string key, val;
    auto db = dbs2[chr];
    auto cur = db->cursor();
    cur->jump();
    auto msad = msadata[chr];
    auto msa = msatrees[chr];
    std::stringstream ss;
    while (cur->get(&key, &val, true)){
        boost::replace_last(key,","," ");
        ss.clear();
        ss.str(key);
        ss >> l;
        ss >> r;
        auto ii = inode(l, l+r);
        msad->push_back(ii) ;
        cnt++;
    } ;
    delete cur;
    for (auto it = msad->begin(); it!=msad->end();++it) msa->insert(*it);
    std::cerr << "Finished inserting, last record:"<<l<<","<<r<<std::endl;
    std::cerr <<cnt<<" records"<<std::endl;
    auto it = msa->lower_bound(inode(1000,1050));
    if (it!=msa->end())
        std::cerr << "Trying to find a record >1000, 1050 :"<<it->l<<","<<it->r<<":"<<it->score<<std::endl;
    return;
}
void mafdb::init_tree()
{
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
    for (auto it:dbpaths)
    {
        auto db = std::shared_ptr <kyotocabinet::HashDB>(new kyotocabinet::HashDB);
        if (!db->open(it.second, kyotocabinet::HashDB::OREADER))
        {
            std::cerr << "open error: " << db->error().name() << std::endl;
            return false;
        }
        dbs[it.first]=db;
    }
    for (auto it:dbpaths2)
    {
        auto db2 = std::shared_ptr <kyotocabinet::HashDB>(new kyotocabinet::HashDB);
        if (!db2->open(it.second, kyotocabinet::HashDB::OREADER))
        {
            std::cerr << "open error: " << db2->error().name() << std::endl;
            return false;
        }
        dbs2[it.first]=db2;
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

