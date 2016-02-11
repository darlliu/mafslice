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
typedef std::vector<std::vector<double>> matrix;
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
        motifmapcompute(const std::vector<INTERVAL_PAIR>& invs,
                const std::string& motif, const matrix& mat, const double& th, const int& cnt):
            invs (invs), SS(""), motif(motif), mat(mat), th(th), id(cnt){};
        void realign(){};
        void score (INTERVAL_PAIR& inv);
        void score_all ()
        {
            SS+="[";
            for (auto &inv: invs)
            {
                SS+="{";
                score(inv);
                SS+="},";
            }
            SS.pop_back();
            SS+="]\n";
        };
        void write ();
        void routine()
        {
            realign();
            score_all();
            write();
            return;
        };
        std::thread spawn()
        {
            return std::thread(&motifmapcompute::routine, this);
        };
    private:
        std::vector<INTERVAL_PAIR> invs;
        std::string  SS, motif;
        matrix mat;
        int id;
        double th;

};

class motifmapdb
{
    public:
        motifmapdb(const std::string& dbname, const std::string& dbp, const std::string& ref):
            dbname(dbname),dbp(dbp), ref(ref), maf(dbname, dbp, ref), outdir("./out/"){};
        motifmapdb(): motifmapdb("test","./test/test.kch", "mm10"){};
        void init(const std::string& dbp_, const std::string& dbname_,
                const std::string& ref_, const std::string& out);
        std::vector<std::pair<INTERVAL_PAIR, INTERVAL_PAIR>> get_intervals(const int& l, const int & r,
                const int& lf, const int& rf)
        {
            std::vector<std::pair<INTERVAL_PAIR, INTERVAL_PAIR>> out;
            auto pp = maf.get_intervals(l-lf,r+rf);
            std::cerr << "Extracting one interval" << std::endl;
            for (; pp.first!=pp.second; ++pp.first)
            {
                std::pair<INTERVAL_PAIR,INTERVAL_PAIR> tmp;
                tmp.first = maf.extract_intervals(*pp.first);
                tmp.second = maf.filter_intervals(l-lf, r+rf, tmp.first);
                out.push_back(tmp);
            }
            return out;
        };
        std::string get_flank(interval& inv, const int & lf, const int & rf);
        void flank(const int& , const int& ,
                std::pair<INTERVAL_PAIR, INTERVAL_PAIR>&);
        void compute_join()
        {
            for (auto &c:threads)
            {
                c.routine();
            }
            std::cerr << "Writing outputs to json files"<<std::endl;
            for (auto &it:OUTPUTS)
            {
                std::cerr <<"Got "<<it.first<<std::endl;
                std::ofstream ofs (outdir+"/"+it.first+".json");
                boost::replace_last(it.second, ",","");
                ofs <<"{" <<it.second<<"}";
                ofs.close();
            }
            OUTPUTS.clear();
            return;
        };
        void compute(const int& cnt , const std::string& motif ,
                const matrix& mat, const double& th)
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
                threads.push_back(motifmapcompute(inps,motif,mat, th, cnt));
            }
            catch (...)
            {
                std::cerr <<"Compute failed on "<<print_interval(inps[0].first)<<std::endl;
            }
            return;
        };
        void add_ref(const std::string& ref)
        {
             refs.push_back(ref);
        }
        void get_inv( const std::string& chr, const int & l,
            const int& r, const int& lf, const int& rf)
        {
            std::cerr << "Now trying to calculate "<<chr<< " "<<l <<" "<<r<<std::endl;
            try
            {
                inps.clear();
                maf.set_chr(chr);
                auto invss = get_intervals(l, r, lf, rf);
                for (auto &invs: invss)
                {
                    flank (lf, rf, invs);
                    auto inp=invs.second;
                    inps.push_back(inp);
                }
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
        std::string dbname, dbp, ref, outdir;
        mafdb maf;
        SEQM seqm;
        std::vector<INTERVAL_PAIR> inps;
        std::vector <std::string> refs;
        std::vector <motifmapcompute> threads;
};

#endif
