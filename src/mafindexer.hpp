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
#include<kccompare.h>
#include<iterator>
#include<boost/functional/hash.hpp>
#include<boost/functional/hash/extensions.hpp>
#include<boost/intrusive/avl_set.hpp>
using namespace boost::intrusive;
class inode :public avl_set_base_hook < optimize_size <true> >
{
    //size_t hash;

    public:
        unsigned l,r;
        short p;
        avl_set_member_hook<> member_hook_;
        inode (const unsigned& ll,
               const unsigned& rr,
               short pp=0):
            l(ll), r(rr), p(pp)  {};
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


class CodedSizeTComparator :public kyotocabinet::Comparator
{
    public:
        CodedSizeTComparator(){};
        int32_t compare(const char* akbuf, size_t aksiz,
                const char* bkbuf, size_t bksiz)
        {
             unsigned l1,r1, l2, r2;
             l1= ((unsigned*)akbuf) [0];
             r1= ((unsigned*)akbuf) [1];
             l2= ((unsigned*)bkbuf) [0];
             r2= ((unsigned*)bkbuf) [1];
             if (l1<l2) return -1;
             if (l1>l2) return 1;
             if (r1>r2) return -1;
             if (r1<r2) return 1;
             return 0;
        };

};
class HashComparator :public kyotocabinet::Comparator
{
    public:
        HashComparator(){};
        int32_t compare(const char* akbuf, size_t aksiz,
                const char* bkbuf, size_t bksiz)
        {
             unsigned l1,l2;
             l1= ((unsigned*)akbuf) [0];
             l2= ((unsigned*)bkbuf) [0];
             if (l1>l2) return 1;
             else if (l1<l2) return -1;
             else if (l1>l2) return 1;
             else return 0;
        };

};
static CodedSizeTComparator CMPSZ;
static HashComparator CMPHS;
typedef avl_set <inode, compare<std::greater<inode>>> ASet;
typedef member_hook <inode, avl_set_member_hook<>, &inode::member_hook_> MemberOption;
typedef avl_multiset <inode, MemberOption> AMSet;

typedef std::map<std::string, std::shared_ptr<AMSet>> MSAMAP;
typedef std::map<std::string, std::shared_ptr<std::vector<inode>>> MSADATA;
typedef std::pair<interval, std::vector<interval>> INTERVAL_PAIR;

class mafdb : public seqdb
{
    friend class boost::serialization::access; //enable boost serialize and be a lazy programmer
    public:
        mafdb (const std::string & name, const std::string& dbp, const std::string& ref)
            : seqdb(name, 1, dbp), ref(ref), init(false){};
        mafdb () : mafdb("defaultMSA","./test/maf","mm10"){ };
        ~mafdb()
        {
            if (init) for (auto &chr:chrs)
                clear_index(chr);
        };
        unsigned hasher(const std::string& s,const unsigned& l , const unsigned& r)
        {
            std::size_t seed = 0;
            boost::hash_combine(seed,l);
            boost::hash_combine(seed,r);
            boost::hash_combine(seed,s);
            return seed;
        };
        unsigned hasher(const unsigned& l , const unsigned& r)
        {
            std::size_t seed = 0;
            boost::hash_combine(seed,l);
            boost::hash_combine(seed,r);
            return seed;
        };
        void clear_index (const std::string&);
        void load_index (const std::string&);
        void save_index(const std::string& , const std::vector<inode>&);
        void init_tree();
        // routine to save/load indices
        /*
         *  Inherited Virtuals
         */

        bool import (const std::string& dirname);
        void import_chr();
        void import_chr(const std::string&);
        bool load_db (const std::string & dbname);
        bool export_db (const std::string & dbname);
        //bool export_db (const std::string & dbname);
        INTERVAL_PAIR extract_intervals (const inode&);
        INTERVAL_PAIR filter_intervals (const unsigned&, const unsigned&,
                INTERVAL_PAIR &,  const bool& masked =true);
        void flank_intervals (const unsigned&, const unsigned&,
                std::vector<interval> &, std::vector<interval> &);

        std::string get(const unsigned& l , const unsigned& r);
        std::string get(const std::string& chr, const unsigned& l , const unsigned& r)
        {
             set_chr(chr);
             return get(l,r);
        };
        auto get_interval(const unsigned& l , const unsigned& r)
        {
            //have to define here to use auto
            if (r==0) throw("Unexpected interval r=0");
            auto msa = msatrees[chr];
            auto it = msa->upper_bound(inode(r+1, r+1));
            auto rit = it;
            int cnt = 0;
            for (; rit!=msa->begin(); --rit){
                if (rit->r > l)
                {
#if DEBUG
                    //if (cnt ==0)
                    std::cerr << " Found a match at : "<< rit->l << " , " << rit->r <<std::endl;
                    ++cnt;
#endif
                }
                else break;
            }
#if DEBUG
            std::cerr << " Found "<<cnt<< " matches ";
#endif
            //check the last one as well as begin() may be included
            //to get content just do first++ -> second
            if (rit->r > l )
                return std::pair<decltype(it) , decltype(it)> (rit, it);
            else
                return std::pair<decltype(it) , decltype(it)> (++rit, it);
        };
        template <class archive>
            void serialize(archive & ar, const unsigned ver)
        {
            ar & boost::serialization::base_object<seqdb>(*this);
            ar & ref;
            ar & sizes;
            ar & TESTS;

            //ar & treesizes;
        };

    private:
        seqdb TESTS;
        std::string ref;
        MSAMAP msatrees;
        MSADATA msadata;
        boost::hash<unsigned> h_;
        boost::hash<std::string> hs_;
        bool init;
        //std::map<std::string, size_t> treesizes;
};

#endif
