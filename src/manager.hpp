#ifndef MANAGER
#define MANAGER
#include "indexer.hpp"
#include "mafindexer.hpp"
#include <queue>
#include <thread>
#include "MOODS/moods_scan.h"

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
        void score ();
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
        void init(const std::string& dbp_, const std::string& dbname_, const std::string& ref_);
        std::pair<INTERVAL_PAIR, INTERVAL_PAIR> get_intervals(const int& l, const int & r)
        {
            std::pair<INTERVAL_PAIR, INTERVAL_PAIR> out;
            auto it = maf.get_interval(l,r);
            std::cerr << "Extrating ... ";
            out.first = maf.extract_intervals(*it);
            std::cerr << "... filtering... "<<std::endl;
            out.second = maf.filter_intervals(l, r, out.first);
            return out;
        };
        std::string get_flank(interval& inv, const int & lf, const int & rf);
        void flank(const int& , const int& ,
                std::pair<INTERVAL_PAIR, INTERVAL_PAIR>&);
        void compute()
        {
            try
            {
                std::cerr <<" Creating compute object...";
                computer c(inp);
                std::cerr <<" Spawning thread...";
                std::thread th = c.spawn();
                std::cerr <<" Joining thread"<<std::endl;
                th.join();
            }
            catch (...)
            {
                std::cerr <<"Compute failed on "<<print_interval(inp.first)<<std::endl;
            }
            return;
        };
        void get(const std::string& chr, const int & l,
                const int& r, const int& lf, const int& rf)
        {
            std::cerr << "Now trying to generate "<<chr<< " "<<l <<" "<<r<<std::endl;
            try
            {
                maf.set_chr(chr);
                auto invs = get_intervals(l, r);
                //flank (lf, rf, invs);
                inp = invs.second;
            }
            catch(std::string s)
            {
                std::cerr << "Got error" << s << std::endl;
                return;
            }
            catch (...)
            {
                std::cerr << "Got unspecified exception..." << std::endl;
                return;
            }
        };

    private:
        std::string dbname, dbp, ref;
        mafdb maf;
        SEQM seqm;
        INTERVAL_PAIR inp;
};

#endif
