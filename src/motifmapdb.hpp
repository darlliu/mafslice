#ifndef MOTIFMAPDB
#define MOTIFMAPDB
#include "indexer.hpp"
#include "mafindexer.hpp"
#include <queue>
#include <thread>
#include <mutex>
#include <Python.h>
#include "MOODS/moods_scan.h"

std::map <std::string, std::string > OUTPUTS;
std::mutex OUTPUTS_MTX;

template <typename T> bool CMP (std::pair<unsigned, T> p, std::pair<unsigned, T> p2)
{return p.first<p2.first;};

//typedef std::priority_queue<std::pair<unsigned, seqdb>,
        //std::list<std::pair<unsigned, seqdb>>, decltype(CMP<seqdb>)> SEQQ;
//typedef std::priority_queue<std::pair<unsigned, mafdb>,
        //std::list<std::pair<unsigned, mafdb>>, decltype(CMP<mafdb>)> MAFQ;
typedef std::map <std::string, seqdb> SEQM;
class releaseGIL{
    public:
        inline releaseGIL(){
            save_state = PyEval_SaveThread();

        }

        inline ~releaseGIL(){
            PyEval_RestoreThread(save_state);

        }
    private:
        PyThreadState *save_state;

};
class motifmapcompute
{
    public:
        motifmapcompute(const INTERVAL_PAIR& inv,
                const std::string& motif):
            inv (inv), SS(""), motif(motif){};
        void realign(){};
        void score ();
        void write ();
        void routine()
        {
            realign();
            score();
            write();
            return;
        };
        std::thread spawn()
        {
            return std::thread(&motifmapcompute::routine, this);
        };
    private:
        INTERVAL_PAIR inv;
        std::string  SS, motif;

};

class motifmapdb
{
    public:
        motifmapdb(const std::string& dbname, const std::string& dbp, const std::string& ref):
            dbname(dbname),dbp(dbp), ref(ref), maf(dbname, dbp, ref), outdir("./out/"){};
        motifmapdb(): motifmapdb("test","./test/test.kch", "mm10"){};
        void init(const std::string& dbp_, const std::string& dbname_,
                const std::string& ref_, const std::string& out);
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
        void compute_join()
        {
            for (auto &c:threads)
            {
                auto th = c.spawn();
                th.join();
            }
            std::cerr << "Writing outputs to json files"<<std::endl;
            for (auto &it:OUTPUTS)
            {
                std::cerr <<"Got "<<it.first<<std::endl;
                std::ofstream ofs (outdir+"/"+it.first+".json");
                ofs <<"{" <<it.second<<"}";
                ofs.close();
            }
        };
        void compute()
        {
            try
            {
                if (threads.size()>100)
                {
                    std::vector <std::thread> ths;
                    std::cerr << "Releasing GIL" <<std::endl;
                    releaseGIL unlock = releaseGIL();
                    std::cerr <<"Too many threads, joining 100"<<std::endl;
                    for (auto &c : threads)
                    {
                         ths.push_back(c.spawn());
                    }
                    for (auto &th:ths)
                        th.join();
                    threads.clear();
                }
                std::cerr <<" Creating compute object...";
                threads.push_back(motifmapcompute(inp,mtf));
                std::cerr <<" Spawning thread...";
            }
            catch (...)
            {
                std::cerr <<"Compute failed on "<<print_interval(inp.first)<<std::endl;
            }
            return;
        };
        void get_inv(const std::string& motif , const std::string& chr, const int & l,
            const int& r, const int& lf, const int& rf)
        {
            std::cerr << "Now trying to generate "<<chr<< " "<<l <<" "<<r<<std::endl;
            try
            {
                mtf=motif;
                maf.set_chr(chr);
                auto invs = get_intervals(l, r);
                flank (lf, rf, invs);
                inp=std::move(invs.second);
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
        std::string get_seq(const std::string& ref, const std::string& chr,
            const int l, const int r, bool strand);

    private:
        std::string dbname, dbp, ref, outdir, mtf;
        mafdb maf;
        SEQM seqm;
        INTERVAL_PAIR inp;
        std::vector <motifmapcompute> threads;
};

#endif
