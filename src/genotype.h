
#ifndef GENOTYPE_H
#define GENOTYPE_H

#include <map>
#include <string>
#include <utility>
#include <set>
#include <tuple>

#include "sequence.h"
#include "util.h"

using namespace std;

class LocVar {
public:
    LocVar();
    LocVar(std::string & read, std::pair<size_t, int> & mraPosLenReadF, std::string & ssr);
    LocVar(std::string & read, std::size_t & mraL, std::size_t & mraR, std::string & ssr);
    LocVar(std::string name, Sequence fp, Sequence rp, Sequence ff, Sequence rf, Sequence repuit, int nRep, unsigned int edCutoff, Sequence mra);
    LocVar(std::string & name, std::string & ff, std::string & rf, std::string & repuit, std::string & mra);
    void print();
    void printSnpMap();
public:
    std::string name;
    Sequence fp;
    Sequence rp;//MUST be reverse complementary seq;
    Sequence ff;
    Sequence rf;
    Sequence repuitAll;
    int repuitAllLen;//if rep1.len == rep2.len, then repAllLen = rep1.len; otherwise = rep1.len + rep2.len; if rep2.len = 0; then repalllen = 0;
    Sequence repuit;
    Sequence repuit2;
    std::string mraBase;
    int nRep;
    int edCutoff;
    Sequence mra;
    //std::vector<Sequence> mras;
    Sequence locSeq;
    Sequence effectiveSeq;
    int locLen;
    int effectiveLen;
    std::string mraName;
    std::map<int, std::string> snpsMapff;
    std::map<int, std::string> snpsMaprf;
    std::map<std::string, std::set<int>> refSnpsSetffMap;
    std::map<std::string, std::set<int>> refSnpsSetrfMap;
    int totalReads;
};

class Variance {
public:
    Variance();
    void cleanVar();
public:
    std::map<int, std::string> subMap;
    std::map<int, std::string> delMap;
    std::map<int, std::string> insMap;
};

class Genotype {
public:
    Genotype();
    void print();
public:
    LocVar baseLocVar;
    int numReads = 0;
    std::string ssrF;
    std::string ssrR;
};


class UnitedGenotype{
public:
    UnitedGenotype();
    
public:
    std::string marker;
    Sequence repuit;
    std::string mraName;
    Sequence mra;
    Sequence effectiveSeq;
    int effectiveLen;
    //std::string seqName;
    std::vector<std::string> seqNameVec;
    std::vector<int> mraSVec;
    std::vector<int> genoSVec;
    std::vector<int> numReadsVec;
    std::vector<LocVar> locVec;
    std::vector<std::string> effectiveSeqVec;
};

struct SimGeno{
    std::string geno;// true genotype, see below;
    std::string oGeno;//original genotype, including error, this is for the homo but has errors eg, 0.93 A vs 0.07 C, so the oGeno is A|C, but the genotype is A|A
    int read1;
    int read2;
    int read3;
    int read4;
    double ratio;//reads/total reads (two haplotype reads)
    bool tORf;
    bool revGeno;
    
    SimGeno(){
        geno = "";
        oGeno = "";
        read1 = 0;
        read2 = 0;
        read3 = 0;
        read4 = 0;
        ratio = 0;
        tORf = false;
        revGeno = false;
    }
};

class UnitedLocSnp {
public:
    UnitedLocSnp();
public:
    std::map<int, std::map<char, int>> preGenoMap; //pos, base ACGT, num reads;
    std::map<int, SimGeno> snpGenoMap; //pos, genotype, reads ratio;
    bool heter;//for haplotype
};

class LocSnp{
public:
    LocSnp();
        
public:
    std::string name;
    Sequence fp;
    Sequence rp;
    Sequence ref;
    std::set<int> snpPosSetHaplo;//only for haplotype snp positions. including ref snps; even if it is inconclusive but not the seq errors.
    std::set<int> snpPosSet;//actually snps including new snps, from the snpsMap;, if it is ref, it include every snps from its variants
    std::set<int> refSnpPosSet;//reference target snps
    std::map<int, std::pair<Sequence, Sequence>> snpsMap;//snps map for that sequence variants, ref and reads might be different;
    int numReads;
    int totReads;
    int totHaploReads;
    double readsRatio;//reads/(total reads - reads);
    std::string genotype;
    bool puGeno;
    void print();
    std::string getGenotype();
    UnitedLocSnp uGeno;
    std::vector<std::tuple<std::string, std::string, int>> haploVec;//seq, haplotype, num reads;
    bool isHaplotype;//inconclusive is false;
    
    std::string genoStr3;//seqerr (ratio > homo and should regarded as seq errors), inconclusive (ratio between homo and heter), homo or heter;//homo also include CC against ref AA;
    //for ref is only homo, heter and inconclusive, for each variants, it could be homo, heter, inconclusive and seqerr (if it homo, but has seq errors)
    
    double ratioHaplo; //ratio = big one / big one + small one; if 1 is homo, if > 0.9 < 1, is homo but with seqerr,  if < jetter is heter, is <homo > jetter is inconclusive.
};

struct SimSnps{
    std::set<int> snpPosSet;//including seq errors and ref snp pos;
    int numReads = 0;
    std::string snpsStr = "";//all snps 
    std::string haploStr = "";//only true haploStr;
    bool isHaplo = false;//if only genoStr8 is seqerr;
    std::string genoStr8 = "seqerr"; //homo, heter1 (high read number), heter2 (low read number), inHeter1, inHeter2, (inconclusive), indel1, indel2 and seqerr;
};

struct SimSnp{
    char snp1 = '\0';
    char snp2 = '\0';
    int reads1 = 0;
    int reads2 = 0;
    double ratio = 0.0;
    std::string color = "transparent";//green, homo with ref, red heter, orange new snps including CC (ref is AA);
    std::string genoStr3 = "inconclusive";//heter, homo, inconclusive;
};

class LocSnp2{
public:
    LocSnp2();
        
public:
    std::string name;//marker name;
    Sequence fp;
    Sequence rp;
    Sequence ref;
    std::pair<int, int> trimPos;//  trim front and back for reference and reads
    Sequence ft;//forward trimming region
    Sequence rt;//reverse trimming region;
    std::set<int> snpPosSetTrueHaplo;//for haplotype; including ref, true haplotype; not the inconclusive ones; not used
    std::set<int> snpPosSetHaplo;//only for true haplotype snp positions. including ref snps; also include the inconclusive ones; for snp table;
    std::set<int> snpPosSet;//snps poitions including seq errors; for alignment table reference
    std::set<int> refSnpPosSet;//reference target snps
    std::set<int> totPosSet;//all the snps pos including ref, combine snpPosSet and refSnpPosSet;
    int totReads;
    int maxReads;
    int totHaploReads;
    int totEffectReads;//total reads after filtering, used to calculate error rate;
    double ratioHaplo; //ratio = big one / big one + small one; if 1 is homo,  if < jetter is heter, is <homo > jetter is inconclusive.
    std::vector<std::tuple<std::string, std::string, int, double, char, char>> haploVec;//seq, haplotype, num reads, num reads/totlhaploreads, conclusive or not (if CC  against ref AA with other heter is inconclusive; eg CA -> ref CG (A|G is inconclusive)), last is indel;
    //bool isHaplotype;//inconclusive is false;
    std::string genoStr3;//seqerr (ratio > homo and should regarded as seq errors), inconclusive (ratio between homo and heter), homo or heter;//homo also include CC against ref AA;;
    bool isIndel;
    //for ref is only homo, heter and inconclusive, for each variants, it could be homo, heter, inconclusive and seqerr (if it homo, but has seq errors)
    std::map<std::string, SimSnps> genoMap;//read is the key
    std::map<int, SimSnp> snpsMap;//position, snps, only include snps in the haloptypes; also the inconclusive ones;
    void print();
    std::map<int, double> baseErrorMap;//for error rate, int: pos, percentage error rate; 
};

struct ComparatorSnp {
    bool operator()( LocSnp& a,  LocSnp& b) const {
        return a.getGenotype() < b.getGenotype();
    }
};

class Sex{
public:
    Sex();
public:
    std::string sexMarker;
    Sequence primerF;
    Sequence primerR;
    Sequence refX;
    Sequence refY;
    int readsX;
    int readsY;
    int minTotalReadsX;
    int minTotalReadsY;
    int minReadsX;//each variant of x;
    int minReadsY;//each variant of y;
    int mismatchesPF;
    int mismatchesPR;
    unsigned int mismatchesRX;
    unsigned int mismatchesRY;
    bool lengthEqual;
    double YXRatio;
    double YXRationCuttoff;
    std::string sexMF;//Male, Female, Inconclusive
    void print();
    std::string getFullRefX();
    std::string getFullRefY();
    std::vector<std::tuple<std::string, int, std::map<int, std::string>>> seqVecX;
    std::vector<std::tuple<std::string, int, std::map<int, std::string>>> seqVecY;
    std::set<int> snpsRefX;
    std::set<int> snpsRefY;
};

#endif /* GENOTYPE_H */

