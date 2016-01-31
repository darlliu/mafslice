#ifndef MANAGER
#define MANAGER
#include "indexer.hpp"
#include "mafindexer.hpp"
#include <queue>
#include <thread>
template <typename T> bool CMP (std::pair<unsigned, T> p, std::pair<unsigned, T> p2)
{return p.first<p2.first;};

//typedef std::priority_queue<std::pair<unsigned, seqdb>,
        //std::list<std::pair<unsigned, seqdb>>, decltype(CMP<seqdb>)> SEQQ;
//typedef std::priority_queue<std::pair<unsigned, mafdb>,
        //std::list<std::pair<unsigned, mafdb>>, decltype(CMP<mafdb>)> MAFQ;
typedef std::map <std::string, seqdb> SEQM;

class computer
{
    public:
        computer(const INTERVAL_PAIR& inv):
            inv (inv){};
        void realign()
        {
            std::cerr << "Reference: " << inv.first.seq<<std::endl;
            return;
        };
        void score ()
        {
            std::cerr << "Scoring!"<<inv.first.l << inv.first.r<<std::endl;
            return;
        };
        void routine()
        {
            realign();
            score();
            return;
        };
        std::thread spawn()
        {
             return std::thread(&computer::routine, this);
        };
    private:
        INTERVAL_PAIR inv;

};

class manager
{
    public:
        manager(const std::string& dbname, const std::string& dbp, const std::string& ref):
            dbname(dbname),dbp(dbp), ref(ref), maf(dbname, dbp, ref){};
        manager(): manager("test","./test/test.kch", "mm10"){};
        void init()
        {
            std::cerr << "Trying to deserialize mafdb from "<< dbname <<std::endl;
            maf.load_db_kch(dbp, dbname);
            std::cerr << "Trying to deserialize reference seqdb"<<std::endl;
            seqm[ref]=seqdb(ref, 1e4);
            seqm[ref].load_db_kch(dbp, ref);
            return;
        };
        std::pair<INTERVAL_PAIR, INTERVAL_PAIR> get(const int& l, const int & r)
        {
            std::pair<INTERVAL_PAIR, INTERVAL_PAIR> out;
            auto it = maf.get_interval(l,r);
            out.first = maf.extract_intervals(*it);
            out.second = maf.filter_intervals(l, r, out.first);
            return out;
        };
        std::string get_flank(const std::string& ref, const int& l, const int& r)
        {
            if (seqm.count(ref)==0)
            {
                std::cerr <<"Adding a seqdb: "<<ref<<std::endl;
                seqm[ref]=seqdb(ref,1e4);
                seqm[ref].load_db_kch(dbp, ref);
            }
            return seqm[ref].get(l,r);
        };
        void flank(const int& , const int& ,
                std::pair<INTERVAL_PAIR, INTERVAL_PAIR>&);
        void compute(INTERVAL_PAIR& in)
        {
            computer c(in);
            std::thread th = c.spawn();
            return;
        };


    private:
        std::string dbname, dbp, ref;
        mafdb maf;
        SEQM seqm;
};

#endif
