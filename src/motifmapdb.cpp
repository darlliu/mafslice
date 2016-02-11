#ifndef MOTIFMAPDB_CPP
#define MOTIFMAPDB_CPP
#include "motifmapdb.hpp"
#include <algorithm>
#include "MOODS/moods_scan.h"
#include "MOODS/match_types.h"

static std::vector <double> mybg= {0.29,0.21,0.21,0.29};
static std::vector <double> myth = {0.0};

void motifmapdb::init(const std::string& dbp_, const std::string& dbname_,
        const std::string& ref_ , const std::string& out)
{
    dbname=dbname_;
    outdir=out;
    dbp=dbp_;
    ref=ref_;
    maf.name=dbname;
    maf.dbpath=dbp;
    maf.ref=ref;
    std::cerr << "Trying to deserialize mafdb from "<< dbname <<std::endl;
    maf.load_db_kch(dbp, dbname);
    std::cerr << "Trying to deserialize reference seqdb"<<std::endl;
    for (auto & rf: refs)
    {
        std::cerr <<" trying to deserialize "<<rf<<std::endl;
        seqm[rf]=seqdb(rf, 1e4);
        seqm[rf].load_db_kch(dbp, rf);
    }
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
    std::cerr << "fw: "<<in<<std::endl;
    std::cerr << "rc: "<<out<<std::endl;
    return out;
};

void motifmapcompute::write()
{
    boost::replace_last(SS, ",","");
    std::string ss="\"";
    //ss+=inv.first.chr+":"
            //+std::to_string(inv.first.l)+","+std::to_string(inv.first.r)+"\":";
    ss+=std::to_string(id)+"\":";
    ss+="{"+SS+"}";
    std::lock_guard<std::mutex> guard (OUTPUTS_MTX);
    if (OUTPUTS.count(motif)==0)
    {
        OUTPUTS[motif]="";
    }
    OUTPUTS[motif]+=ss+",\n";
    return;
}
void motifmapcompute::score(INTERVAL_PAIR& inv)
{
    auto print_matches = [&] (std::vector<match>& matches, interval& inv)
    {
        if (matches.size()==0) return std::string("-50");
        std::stringstream ss("");
        match mx;
        mx.score = matches[0].score;
        for (auto &m: matches)
            if (m.score > mx.score)
            {
                mx.pos = m.pos;
                mx.score =m.score;
            }
#if DEBUG
        ss<<"{";
        if (inv.strand)
            ss<<"\""<<mx.pos+inv.l<< "\" : "<<mx.score<<"}";
        else
            ss<<"\""<<inv.r-mx.pos<< "\" : "<<mx.score<<"}";
#else
        ss <<mx.score;
#endif
        return ss.str();
    };
#if DEBUG
    std::cerr << "Scoring: "<<print_interval(inv.first)<<std::endl;
#endif
    //auto matches_ = MOODS::scan::scan_dna(inv.first.seq, my_matrix,mybg, std::vector<double>(){-50});
    auto matches = MOODS::scan::naive_scan_dna(inv.first.seq, mat, -50);
    auto matches2= MOODS::scan::naive_scan_dna(get_reverse_comp(inv.first.seq), mat, -50);
    matches.insert(matches.end(), matches2.begin(), matches2.end());
#if DEBUG
    std::cerr <<" Ref score: "<< print_matches(matches, inv.first) <<std::endl;
#endif
    SS+= "\"" + inv.first.ref+"\" :"+print_matches(matches, inv.first)+",";
    for (auto & iv : inv.second)
    {
        //matches_ = MOODS::scan::scan_dna(iv.seq, my_matrix,mybg, myth);
        matches = MOODS::scan::naive_scan_dna(iv.seq, mat, th);
        matches2 = MOODS::scan::naive_scan_dna(get_reverse_comp(iv.seq), mat, th);
        matches.insert(matches.end(), matches2.begin(), matches2.end());
        if (matches.size()>0)
        {
#if DEBUG
            std::cerr <<" Matches for "<<print_interval(iv)<<std::endl;
#endif
            auto ss = print_matches(matches, iv);
#if DEBUG
            std::cerr << ss <<std::endl;
#endif
            SS += "\""+iv.ref+"\" :"+ss+", ";
        }
    }
    return;
}

void motifmapdb::flank(const int& lf, const int& rf,
        std::pair<INTERVAL_PAIR, INTERVAL_PAIR>& in)
{
    auto inner = [&](interval& inv, interval& inv2,const int& lf, const int& rf)
    {
        if ((inv2.l>lf)&&(inv2.r+rf+inv2.l<inv.seq.size()))
        {
            inv2.seq=inv.seq.substr(inv2.l-lf, inv2.r+rf+lf);
        }
        else
        {
            int start = inv.l+inv2.l;
            int stop = inv.l+inv2.l+inv2.r;
            inv2.l=start;
            inv2.r=stop;
            if (!inv2.strand)
            {
                auto sz= seqm[inv2.ref].sizes[inv2.chr];
                start = sz-inv.r;
                stop = sz-inv.l-1;
                inv2.l=start;
                inv2.r=stop;
            }
            inv2.seq=get_flank(inv2, lf, rf);
        }
//#if DEBUG
        std::cerr <<" Flanked "<<print_interval(inv2)<<std::endl;
//#endif
    };
    inner(in.first.first,in.second.first,0, 0);
    for (int i=0; i<in.second.second.size();++i)
    {
        std::cerr << "Before flank"<<print_interval(in.second.second[i])<<std::endl;
        inner(in.first.second[i], in.second.second[i], lf, rf);
    }
    return;
}

std::string motifmapdb::get_flank(interval& inv, const int & lf, const int & rf)
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
    std::cerr <<" get flank .."<<inv.ref;
    if (inv.strand)
        return seqm[inv.ref].get(inv.chr, inv.l-lf,inv.r+rf);
    else
        return seqm[inv.ref].get(inv.chr, inv.l-lf,inv.r+rf);
    //std::cerr <<".. done";
}

std::string motifmapdb::get_seq(const std::string& ref, const std::string& chr,
    const int l, const int r, bool strand)
{
    auto sq = seqm[ref];
    std::string ss;
    if (strand)
        ss = sq.get(chr, l, r);
    else
    {
        auto sz = sq.sizes[chr];
        int start = sz-r;
        int stop = sz-l;
        ss = sq.get(chr, start, stop);
        ss= get_reverse_comp(ss);
    }
    return ss;
}

#endif
