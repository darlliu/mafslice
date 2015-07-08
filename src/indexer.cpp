#include "indexer.hpp"

bool seqdb::import (const std::string& dirname)
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
            std::cout << "found file : " << it->path().string()<< std::endl;
            fps.push_back(it->path().string());
            fns.push_back(it->path().filename());

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
    init_db(chrs);
    for (unsigned i =0; i < fps.size(); ++i)
    {
        auto chr = chrs[i];
        set_chr(chr);
        auto chrp = fps[i];
        import_chr(chrp);
    }
    return true;
};

bool seqdb::import_chr(const std::string & chrfp)
{
    auto db = dbs[chr];
    if (!db->open(dbpath+"/"+chr+".kch", kyotocabinet::HashDB::OWRITER | kyotocabinet::HashDB::OCREATE))
    {
        std::cerr << "open error: " << db->error().name() << std::endl;
    }
    unsigned long idx = 0;
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
            if (!db->set(std::to_string(idx), tmp.substr(0, chunksz))){
                throw( "Error adding a value to db "+ chr);
            };
            idx+=chunksz;
            tmp = tmp.substr(chunksz);
        }
    }
    ifs.close();
    db->set(std::to_string(idx), tmp);
    idx+= tmp.size();
    std::cout << "Loaded length "<<idx<<std::endl;
    sizes[chr] = idx;
    return true;
};

void seqdb::init_db (const std::vector<std::string>& chrs)
{
    std::cout << "Initializing DBs, length: " <<chrs.size() <<std::endl;
    for (auto chr: chrs)
    {
        dbs[chr]=std::shared_ptr <kyotocabinet::HashDB>(new kyotocabinet::HashDB);
    }
    return;
};

std::string seqdb::get(const unsigned long& l, const unsigned long& r)
{
    //First, check the sizes and make sure that the
    if (r<=l || r >= sizes[chr])
    {
        std::cerr << "Index for slice is incorrect"<< std::endl;
        return "";
    }
    auto idx = (l/chunksz) * chunksz; //integer division on chunksz
    auto idx0 = idx;
    auto db = dbs[chr];
    std::string val(""), tmp;
    do
    {
        if (!db->get(std::to_string(idx),&tmp))
        {
            std::cerr << "Get Error : " << db->error().name() << " On : "<< chr << std::to_string(idx) << std::endl;
            return "";
        }
        val += tmp;
        idx += chunksz;
    }
    while (r > idx);
    return val.substr(l-idx0, r-l);

}
