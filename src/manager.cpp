#ifndef MANAGER_CPP
#define MANAGER_CPP
#include "manager.hpp"
#include "MOODS/moods_scan.h"
void manager::flank(const int& lf, const int& rf,
        std::pair<INTERVAL_PAIR, INTERVAL_PAIR>& in)
{
    auto inner = [&](interval& inv, interval& inv2)
    {
        if ((inv2.l-lf>inv.l)&&(inv2.r+rf<inv.r))
        {
            inv2.seq=inv.seq.substr(inv2.l-inv.l-lf, inv2.r-inv2.l+rf);
        }
        else
        {
            inv2.seq=get_flank(inv, lf, rf);
        }

        inv2.l-=lf;
        inv2.r+=rf;
        std::cerr <<" Flanked "<<print_interval(inv2)<<std::endl;
    };
    inner(in.first.first,in.second.first);
    for (int i=0; i<in.second.second.size();++i)
    {
        inner(in.first.second[i], in.second.second[i]);
    }
    return;
}

std::string manager::get_flank(interval& inv, const int & lf, const int & rf)
{
    if (seqm.count(inv.ref)==0)
    {
        std::cerr <<"Adding a seqdb: "<<inv.ref<<std::endl;
        seqm[inv.ref]=seqdb(inv.ref,1e4);
        seqm[inv.ref].load_db_kch(dbp, inv.ref);
    }
    if (!inv.strand)
        return seqm[inv.ref].get(inv.chr, inv.l-lf,inv.r+rf);
    else
    {
        std::string val;
        seqm[inv.ref].dbs[inv.ref][1]->get(inv.chr, &val);
        unsigned sz;
        std::stringstream ss(val);
        ss>>sz;
        int start = sz+1-inv.r;
        int stop = sz+1-inv.l;
        return seqm[inv.ref].get(inv.chr,start-rf, stop+lf);
    }
}

#endif
