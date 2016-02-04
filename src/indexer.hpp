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
#if USE_DBT
typedef kyotocabinet::TreeDB _DB;
#else
typedef kyotocabinet::HashDB _DB;
#endif
typedef std::map<std::string, std::vector<std::shared_ptr<_DB>>> DB;


typedef boost::archive::binary_oarchive BOARCHIVE;
typedef boost::archive::binary_iarchive BIARCHIVE;
typedef boost::archive::text_oarchive OARCHIVE;
typedef boost::archive::text_iarchive IARCHIVE;

typedef enum {
    hash = 0,
    tree
} DBTYPE;
typedef enum {
    increment=0,
    custom
} INDEXTYPE; //these are optional identifiers left unused for now
typedef enum {
    reference=0,
    matching,
    other
} INTERVALTYPE;

struct interval
{
    public:
        std::string ref, chr, seq;
        unsigned l=-1,r=-1;
        INTERVALTYPE tt = reference;
        float score = 0.0;
        bool strand=true;
};
std::string print_interval(const interval& in);

class seqdb {
    friend class boost::serialization::access; //enable boost serialize and be a lazy programmer
    public:
        seqdb (const std::string & name, const size_t& sz, const std::string& dbp = "./test",
                const INDEXTYPE& idxtype=increment):
            name (name), chunksz(sz), indextype(idxtype),
            dbpath(dbp), chr(""), assemble(false), scaffold(false){};
        seqdb(): seqdb("default", 1e4) {};
        ~seqdb()
        {
            close_db(dbs);
        };
        virtual void import (const std::string&);
        virtual void import_feed();
        //import from a fasta file and build a db
        virtual void import_chr();
        virtual void import_chr_fsq();
        //import a chromosome file
        virtual void import_chr(const std::string&);
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
        virtual std::string get (const interval& itv) {return get(itv.chr, itv.l, itv.r);};
        void init_db(const std::vector<std::string>&, DB&);
        void close_db(DB& dbs){
            for (auto it=dbs.begin(); it!=dbs.end();++it)
                for (auto &iit: it->second)
                    iit->close();
            dbs.clear();
        };
        virtual void load_db(const std::string &);
        virtual void load_db_kch(const std::string &, const std::string&);
        virtual void load_db_();
        virtual void export_db(const std::string &);
        virtual void export_db_kch(const std::string &);
        void set_chr (const std::string& c) {chr = c;} ;
        void set(const std::string&, const std::string & );
        void del(const std::string&);
        virtual size_t get_index(const size_t& idx ) {return idx/chunksz * chunksz;};
        template <class archive>
            void serialize(archive & ar, const unsigned ver)
            {
                ar & name;
                ar & scaffold;
                ar & dbpath;
                ar & chunksz;
                ar & chrs;
                ar & indextype;
                ar & indices;
                ar & sizes;
                ar & postfixes;
                ar & fapaths;
                ar & dbpaths;
                //we cannot serialize the dbs -- they have to be loaded in load_db!
            };
        size_t chunk_sz(){return chunksz;};


        std::string name, dbpath; //name of the db and the directory path for kyotocabinet
        size_t chunksz; //size of each chunk of sequence
#if USE_DBT
        DBTYPE dbtype=tree; // type of db to be intialized -> only affects db construction
#else
        DBTYPE dbtype=hash; // type of db to be intialized -> only affects db construction
#endif
        std::string chr;
        std::vector<std::string> chrs;
        INDEXTYPE indextype;
        INDEXMAP indices;
        DB dbs;
        SIZES sizes, postfixes;
        NAMES fapaths, dbpaths;
        bool assemble, scaffold;
        std::vector <std::string> tmp_dbpaths;
};

#endif

