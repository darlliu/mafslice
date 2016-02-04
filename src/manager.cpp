#ifndef MANAGER_CPP
#define MANAGER_CPP
#include "manager.hpp"
#include <algorithm>
#include "MOODS/moods_scan.h"
#include "MOODS/match_types.h"
static std::vector<std::vector<std::vector<double>>> my_matrix=
{{{-1.92069227,  1.42528548,  1.53939197, -2.85922073, -2.85922073,
 -2.85922073, -1.07326227, -1.92069227, -2.85922073, -2.85922073,
 -2.85922073, -2.85922073, -2.85922073},
 { 1.89094905, -2.80559053, -2.80559053, -2.80559053, -0.6075987 ,
 -2.80559053, -2.80559053, -2.80559053, -2.80559053, -2.80559053,
 -2.80559053, -2.80559053,  2.00505554},
 {-0.6075987 , -2.80559053, -2.80559053,  2.00505554,  1.89094905,
 -0.6075987 , -1.5103015 ,  2.11079479,  2.20931053, -2.80559053,
 -2.80559053,  2.20931053, -2.80559053},
 {-2.85922073, -0.5222471 , -1.07326227, -1.07326227, -1.92069227,
 1.53939197,  1.42528548, -2.85922073, -2.85922073,  1.74364696,
 1.74364696, -2.85922073, -1.07326227}}}
;
static std::vector<std::string> refs = {"ailMel1", "anoCar2", "bosTau7", "calJac3",
"canFam3", "cavPor3", "choHof1", "chrPic1", "danRer7", "dasNov3", "dipOrd1", "echTel1", "equCab2",
"eriEur1", "felCat5", "fr3", "gadMor1", "galGal4", "gasAcu1", "gorGor3", "hetGla2", "hg19", "latCha1",
"loxAfr3", "macEug2", "melGal1", "melUnd1", "micMur1", "monDom5", "myoLuc2", "nomLeu2", "ochPri2", "oreNil2",
"ornAna1", "oryCun2", "oryLat2", "otoGar3", "oviAri1", "panTro4", "papHam1", "petMar1", "ponAbe2", "proCap1",
"pteVam1", "rheMac3", "rn5", "saiBol1", "sarHar1", "sorAra1", "speTri2", "susScr3", "taeGut1", "tarSyr1",
"tetNig2", "triMan1", "tupBel1", "turTru2", "vicPac1", "xenTro3"};
static std::vector <double> mybg= {0.29,0.21,0.21,0.29};
static std::vector <double> myth = {4.27};

void manager::init(const std::string& dbp_, const std::string& dbname_, const std::string& ref_)
{
    dbname=dbname_;
    dbp=dbp_;
    ref=ref_;
    maf.name=dbname;
    maf.dbpath=dbp;
    maf.ref=ref;
    std::cerr << "Trying to deserialize mafdb from "<< dbname <<std::endl;
    maf.load_db_kch(dbp, dbname);
    //std::cerr << "Trying to deserialize reference seqdb"<<std::endl;
    //seqm[ref]=seqdb(ref, 1e4);
    //seqm[ref].load_db_kch(dbp, ref);
    //for (auto & rf: refs)
    //{
        //std::cerr <<" trying to deserialize "<<rf<<std::endl;
        //seqm[rf]=seqdb(rf, 1e4);
        //seqm[rf].load_db_kch(dbp, rf);
    //}
    return;
};
void print_matches(std::vector<match> matches)
{
    for (auto &m: matches)
        std::cerr <<"("<<m.pos<< ": "<<m.score<<"); ";
    std::cerr<<std::endl;
    return;
};

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
void computer::score()
{
    std::cerr << "Scoring: "<<print_interval(inv.first)<<std::endl;
    auto matches = MOODS::scan::scan_dna(inv.first.seq, my_matrix,mybg, myth);
    std::cerr <<" Ref score: ";
    print_matches(matches[0]);
    for (auto & iv : inv.second)
    {
        matches = MOODS::scan::scan_dna(iv.seq, my_matrix,mybg, myth);
        if (matches[0].size()>0)
        {
            std::cerr <<" Matches for "<<print_interval(iv)<<std::endl;
            print_matches (matches[0]);
        }
    }
    return;
}

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
                inv2.seq=inv.seq.substr(inv.r-inv2.r-rf, inv2.r-inv2.l+rf+lf);
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
    std::cerr <<" get flank ..";
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
    std::cerr <<".. done";
}

#endif
