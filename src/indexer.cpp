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
    return true;
};

void seqdb::import_chr()
{
    auto chrfp = fapaths[chr];
    auto dbp = dbpaths[chr];
    auto db = dbs[chr][0];
    if (!db->open(dbp, _DB::OWRITER | _DB::OCREATE))
    {
        std::cerr << "open error: " << db->error().name() << std::endl;
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
    if (r<=l || r >= sizes[chr])
    {
        std::cerr << "Index for slice is incorrect"<< std::endl;
        return "";
    }
    auto idx = get_index(l); //integer division on chunksz
    auto idx0 = idx;
    auto db = dbs[chr][0];
    std::string val(""), tmp;
    do
    {
        if (!db->get(std::to_string(idx),&tmp))
        {
            std::cerr << "Get Error : " << db->error().name() \
                << " On : "<< chr << std::to_string(idx) << std::endl;
            return "";
        }
        val += tmp;
        idx += chunksz;
    }
    while (r > idx);
    return val.substr(l-idx0, r-l);

};

bool seqdb::export_db(const std::string& fp )
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

bool seqdb::load_db(const std::string& fp )
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
        auto db = std::shared_ptr <_DB>(new _DB);
        if (!db->open(it.second, _DB::OREADER))
        {
            std::cerr << "open error: " << db->error().name() << std::endl;
            return false;
        }
        dbs[it.first].push_back(db);
    }
    return true;
};
