/*
 * MAF indexer based on the normal indexer
 * 1, implement a similar parser and storer via KC
 * 2, parses maf files and store hits as records, indexed by hash(start, stop)
 * 3, implement a container based on AVL tree or B tree to retrieve the hashes
 * 4, these trees are deserialized at runtime
 *
 * Typical high level routine:
 * 1, initiate MSA indexer
 * 2, initiate Seq
 */

#ifndef MAFINDEXER
#define MAFINDEXER
#include"indexer.hpp"
#include<boost/functional/hash.hpp>
#include<boost/functional/hash/extensions.hpp>
#include<boost/intrusive/avl_set.hpp>
#include<boost/container/flat_map.hpp>
using namespace boost::intrusive;
class inode :public avl_set_base_hook < optimize_size <true> >
{
    //size_t hash;

    public:
        size_t l,r;
        float score;
        avl_set_member_hook<> member_hook_;
        inode (const size_t& ll,
               const size_t& rr,
               float sc=-1):
            l(ll), r(rr), score(sc)  {};
        friend bool operator < (const inode& l, const inode& r)
        {
            if (l.l == r.l)
                return (l.r < r.r);
            else return l.l<r.l;
        };
        friend bool operator > (const inode& l, const inode& r)
        {
            if (l.l == r.l)
                return (l.r > r.r);
            else return l.l > r.l;
        };
        friend bool operator == (const inode& l, const inode& r)
        {
            if (l.l==r.l && l.r == r.r) return true;
            return false;
            //return l.hash == r.hash;
        };

};

typedef avl_set <inode, compare<std::greater<inode>>> ASet;
typedef member_hook <inode, avl_set_member_hook<>, &inode::member_hook_> MemberOption;
typedef avl_multiset <inode, MemberOption> AMSet;

class mafdb : public seqdb
{
    typedef std::map<std::string, std::shared_ptr<AMSet>> MSAMAP;
    typedef std::map<std::string, std::shared_ptr<std::vector<inode>>> MSADATA;
    friend class boost::serialization::access; //enable boost serialize and be a lazy programmer
    public:
        mafdb (const std::string & name, const std::string& dbp, const std::string& ref)
            : seqdb(name, 1, dbp), ref(ref){};
        mafdb () : mafdb("defaultMSA","./test/maf","mm10"){ };
        ~mafdb()
        {
            close_db (dbs2);
            for (auto &chr:chrs)
                clear_index(chr);
        };
        size_t hasher(const std::string& s,const size_t& l , const size_t& r)
        {
            std::size_t seed = 0;
            boost::hash_combine(seed,l);
            boost::hash_combine(seed,r);
            boost::hash_combine(seed,s);
            return seed;
        };
        size_t hasher(const size_t& l , const size_t& r)
        {
            std::size_t seed = 0;
            boost::hash_combine(seed,l);
            boost::hash_combine(seed,r);
            return seed;
        };
        void clear_index (const std::string&);
        void load_index (const std::string&);
        void init_tree();
        // routine to save/load indices
        /*
         *  Inherited Virtuals
         */

        bool import (const std::string& dirname);
        void import_chr();
        bool load_db (const std::string & dbname);
        bool export_db (const std::string & dbname);
        //bool export_db (const std::string & dbname);
        std::string get(const size_t& l , const size_t& r);
        std::string get(const std::string& key);
        template <class archive>
            void serialize(archive & ar, const unsigned ver)
        {
            ar & boost::serialization::base_object<seqdb>(*this);
            ar & ref;
            ar & sizes;
            ar & dbpaths2;
            //ar & treesizes;
        };

    private:
        std::string ref;
        MSAMAP msatrees;
        MSADATA msadata;
        boost::hash<size_t> h_;
        boost::hash<std::string> hs_;
        //std::map<std::string, size_t> treesizes;
        DB dbs2;
        NAMES dbpaths2;
};

#endif
