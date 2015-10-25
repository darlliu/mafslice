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
using namespace boost::intrusive;
class itree :public avl_set_base_hook < optimize_size <true> >
{
    size_t hash;
    size_t l,r;
    bool negative;
    public:
        avl_set_member_hook<> member_hook_;
        itree (const size_t& ll,
                const size_t& rr,
                const size_t & h,
                bool neg=false):
            l(ll), r(rr), hash(h), negative(neg) {};
        friend bool operator < (const itree& l, const itree& r)
        {
            if (l.l == r.l)
                return (l.r < r.r);
            else return l.l<r.l;
        };
        friend bool operator > (const itree& l, const itree& r)
        {
            if (l.l == r.l)
                return (l.r > r.r);
            else return l.l > r.l;
        };
        friend bool operator == (const itree& l, const itree& r)
        {
            //if (l.l==r.l && l.r == r.r) return true;
            //return false;
            return l.hash == r.hash;
        };

};

typedef avl_set <itree, compare<std::greater<itree>>> ITree;
typedef member_hook <itree, avl_set_member_hook<>, &itree::member_hook_> MemberOption;
typedef avl_multiset <itree, MemberOption> ATree;

class mafdb : public seqdb
{
    typedef std::map<std::string, ATree> RNGMAP;
    public:
        mafdb (const std::string & name, const std::string& dbp)
            : seqdb(name, 1, dbp){};
        mafdb () : mafdb("default","./test"){ };
        size_t hasher(const size_t& l , const size_t& r)
        {
            std::size_t seed = 0;
            boost::hash_combine(seed,l);
            boost::hash_combine(seed,r);
            return seed;
        };
        void save_indexes (const std::string&);
        void load_indexes (const std::string&);
        // routine to save/load indices
    private:
        RNGMAP indextrees;
        boost::hash< size_t > h_;
        std::vector<size_t> treesizes;
        std::map < std::string, std::vector<itree> > treeraw;

};

#endif
