#ifndef MANAGER_CPP
#define MANAGER_CPP
#include "manager.hpp"
void manager::flank(const int& lf, const int& rf,
        std::pair<INTERVAL_PAIR, INTERVAL_PAIR>& in)
{
    auto inner = [&](interval& inv, interval& inv2)
    {
        if ((inv2.l-lf>inv.l)&&(inv2.r+rf<inv.r))
        {
            inv2.seq=inv.seq.substr(inv2.l-inv.l, inv2.r-inv2.l);
        }
        else
        {
            inv2.seq=get_flank(inv2.ref,inv2.chr,inv2.l-lf, inv2.r+rf);
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
};
#endif
