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

typedef enum {
    hash = 0,
    tree
} DBTYPE;
typedef enum {
    increment=0,
    custom
} INDEXTYPE;

class seqdb {
    typedef std::vector<long> INDICES;
    typedef std::map<std::string, INDICES> INDEXMAP;
    typedef std::map<std::string, std::string> NAMES;
    typedef std::map<std::string, unsigned long> SIZES;
    typedef std::map<std::string, std::shared_ptr<kyotocabinet::HashDB>> DB;
    typedef std::map<std::string, std::shared_ptr<kyotocabinet::TreeDB>> DBT;
#if USE_BINARY_ARCHIVE
    typedef boost::archive::binary_oarchive OARCHIVE;
    typedef boost::archive::binary_iarchive IARCHIVE;
#else
    typedef boost::archive::text_oarchive OARCHIVE;
    typedef boost::archive::text_iarchive IARCHIVE;
#endif
    friend class boost::serialization::access; //enable boost serialize and be a lazy programmer
    public:
        seqdb (const std::string & name, const unsigned long& sz,
                const INDEXTYPE& idxtype=increment, const DBTYPE& dbtype=hash):
            name (name), chunksz(sz), indextype(idxtype), dbtype(dbtype),
            dbpath("./test"), chr(""){};
        seqdb(): seqdb("default", 1e4) {};
        ~seqdb()
        {
            close_db();
        };

        virtual bool import (const std::string&);
        //import from a fasta file and build a db
        void import_chr(const std::string&, const std::string&);
        //import a chromosome file

        virtual std::string get (const std::string& chr,
                const int& l, const int& r)
        {
            this->set_chr(chr);
            return this->get(l,r);
        };

        virtual std::string get (const unsigned long&, const unsigned long&);
        virtual std::string get (const std::string& key) {return "";};

        void init_db(const std::vector<std::string>&);
        void close_db(){
            for (auto it=dbs.begin(); it!=dbs.end();++it)
                it->second->close();
        };
        virtual bool load_db(const std::string &);
        virtual bool export_db(const std::string &);
        void set_chr (const std::string& c) {chr = c;} ;
        void set(const std::string&, const std::string & );
        void del(const std::string&);
        virtual unsigned long get_index(const unsigned long& idx ) {return idx/chunksz * chunksz;};
        template <class archive>
            void serialize(archive & ar, const unsigned ver)
            {
                ar & name;
                ar & dbpath;
                ar & chunksz;
                ar & chr;
                ar & indextype;
                ar & indices;
                ar & sizes;
                ar & fapaths;
                ar & dbpaths;
                //we cannot serialize the dbs -- they have to be loaded in load_db!
            };
        unsigned long chunk_sz(){return chunksz;};


    private:
        std::string name, dbpath; //name of the db and the directory path for kyotocabinet
        unsigned long chunksz; //size of each chunk of sequence
        DBTYPE dbtype; // type of db to be intialized -> only affects db construction
        std::string chr;
        INDEXTYPE indextype;
        INDEXMAP indices;
        SIZES sizes;
        NAMES fapaths, dbpaths;
        bool init=false;
        DB dbs;
};

#endif

