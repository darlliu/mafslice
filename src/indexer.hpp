/*
 * A base indexer for sequentially indexed contents
 * contains a map from string (by chromosome for example) to vector of indices
 * and a sequentially sorted vector with a fixed step size from indices to keys
 * and a on disk hash table from keys to actual content (chunks of sequences)
 * Functionalities are provided to query sequences from any slice of indices on major key
 */


#ifndef INDEXER
#define INDEXER
#include "base.hpp"

typedef std::vector<size_t> INDICES;
typedef std::map<std::string, INDICES> INDEXMAP;
typedef std::map<std::string, std::string> NAMES;
typedef std::map<std::string, size_t> SIZES;
typedef std::map<std::string, std::shared_ptr<kyotocabinet::HashDB>> DB;
typedef std::map<std::string, std::shared_ptr<kyotocabinet::TreeDB>> DBT;
#if USE_BINARY_ARCHIVE
typedef boost::archive::binary_oarchive OARCHIVE;
typedef boost::archive::binary_iarchive IARCHIVE;
#else
typedef boost::archive::text_oarchive OARCHIVE;
typedef boost::archive::text_iarchive IARCHIVE;
#endif

typedef enum {
    hash = 0,
    tree
} DBTYPE;
typedef enum {
    increment=0,
    custom
} INDEXTYPE; //these are optional identifiers left unused for now

class seqdb {
    friend class boost::serialization::access; //enable boost serialize and be a lazy programmer
    public:
        seqdb (const std::string & name, const size_t& sz, const std::string& dbp = "./test",
                const INDEXTYPE& idxtype=increment, const DBTYPE& dbtype=hash):
            name (name), chunksz(sz), indextype(idxtype), dbtype(dbtype),
            dbpath(dbp), chr(""){};
        seqdb(): seqdb("default", 1e4) {};
        ~seqdb()
        {
            close_db(dbs);
        };
        virtual bool import (const std::string&);
        //import from a fasta file and build a db
        virtual void import_chr();
        //import a chromosome file
        std::string get3 (const std::string& chr,
                const size_t& l, const size_t& r)
        {
            return this->get(chr, l,r);//this is an adaptor for python extension
        };
        virtual std::string get (const std::string& chr,
                const size_t& l, const size_t& r)
        {
            this->set_chr(chr);
            return this->get(l,r);
        };
        virtual std::string get (const size_t&, const size_t&);
        virtual std::string get (const std::string& key) {return "";};
        void init_db(const std::vector<std::string>&, DB&);
        void close_db(DB& dbs){
            for (auto it=dbs.begin(); it!=dbs.end();++it)
                it->second->close();
            dbs.clear();
        };
        virtual bool load_db(const std::string &);
        virtual bool export_db(const std::string &);
        void set_chr (const std::string& c) {chr = c;} ;
        void set(const std::string&, const std::string & );
        void del(const std::string&);
        virtual size_t get_index(const size_t& idx ) {return idx/chunksz * chunksz;};
        template <class archive>
            void serialize(archive & ar, const unsigned ver)
            {
                ar & name;
                ar & dbpath;
                ar & chunksz;
                ar & chrs;
                ar & indextype;
                ar & indices;
                ar & sizes;
                ar & fapaths;
                ar & dbpaths;
                //we cannot serialize the dbs -- they have to be loaded in load_db!
            };
        size_t chunk_sz(){return chunksz;};


    protected:
        std::string name, dbpath; //name of the db and the directory path for kyotocabinet
        size_t chunksz; //size of each chunk of sequence
        DBTYPE dbtype; // type of db to be intialized -> only affects db construction
        std::string chr;
        std::vector<std::string> chrs;
        INDEXTYPE indextype;
        INDEXMAP indices;
        DB dbs;
        SIZES sizes;
        NAMES fapaths, dbpaths;
};

#endif

