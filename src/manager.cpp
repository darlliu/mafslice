#ifndef MANAGER_CPP
#define MANAGER_CPP
#include "manager.hpp"
#include <algorithm>
#include "MOODS/moods_scan.h"
std::string get_reverse_comp(const std::string& in)
{
    std::string out;
    out.reserve(in.size());
    for (auto &c:in)
    {
        switch(c){
            case 'A': {out.push_back('T');break;}
            case 'T': {out.push_back('A');break;}
            case 'U': {out.push_back('A');break;}
            case 'G': {out.push_back('C');break;}
            case 'C': {out.push_back('G');break;}
            case 'a': {out.push_back('t');break;}
            case 't': {out.push_back('a');break;}
            case 'u': {out.push_back('a');break;}
            case 'g': {out.push_back('c');break;}
            case 'c': {out.push_back('g');break;}
            default: out.push_back(c);
        }
    }
    std::reverse(out.begin(), out.end());
    return out;
};
void manager::flank(const int& lf, const int& rf,
        std::pair<INTERVAL_PAIR, INTERVAL_PAIR>& in)
{
    auto inner = [&](interval& inv, interval& inv2)
    {
        if ((inv2.l-lf>inv.l)&&(inv2.r+rf<inv.r))
        {
            if (inv2.strand)
                inv2.seq=inv.seq.substr(inv2.l-inv.l-lf, inv2.r-inv2.l+rf+lf);
            else
            {
                //std::cerr <<"Getting negative strand "<<inv.r-inv2.r-rf <<" + "<<inv2.r-inv2.l+rf+lf<<std::endl;
                //std::cerr << print_interval(inv) << print_interval (inv2)<<std::endl;
                inv2.seq=inv.seq.substr(inv.r-inv2.r-rf, inv2.r-inv2.l+rf+lf);
            }
        }
        else
        {
            inv2.seq=get_flank(inv2, lf, rf);
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
        try
        {
            std::cerr <<"Adding a seqdb: "<<inv.ref<<std::endl;
            seqm[inv.ref]=seqdb(inv.ref,1e4);
            seqm[inv.ref].load_db_kch(dbp, inv.ref);
        }
        catch(std::string err)
        {
            std::cerr << "Got error "<<err <<std::endl;
            return "";
        }
    }
    if (inv.strand)
        return seqm[inv.ref].get(inv.chr, inv.l-lf,inv.r+rf);
    else
    {
        try
        {
            unsigned sz;
            if (seqm[inv.ref].scaffold)
            {
                std::string val;
                seqm[inv.ref].dbs[inv.ref][1]->get(inv.chr, &val);
                std::stringstream ss(val);
                ss>>sz;
            } else {
                sz = seqm[inv.ref].sizes[inv.chr];
            }
#if DEBUG
            std::cerr <<" ref: "<<inv.ref << " , "<<inv.chr <<" size : "<<sz<<std::endl;
#endif
            int start = sz-inv.r;
            int stop = sz-inv.l;
            return get_reverse_comp(seqm[inv.ref].get(inv.chr,start-rf, stop+lf));
        }
        catch(std::string err)
        {
            std::cerr << "Got error "<<err <<std::endl;
            return "";
        }
    }
}

#endif
