
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
    Sequence rp;
    Sequence ff;
    Sequence rf;
    Sequence repuitAll;
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

class LocSnp{
public:
    LocSnp();
        
public:
    std::string name;
    Sequence fp;
    Sequence rp;
    Sequence ref;
    std::set<int> snpPosSet;
    std::map<int, std::pair<Sequence, Sequence>> snpsMap;
    int numReads;
    void print();
    std::string getGenotype();
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
    std::vector<std::tuple<std::string, int, std::map<int, std::string>>> seqVecX;
    std::vector<std::tuple<std::string, int, std::map<int, std::string>>> seqVecY;
    std::set<int> snpsRefX;
    std::set<int> snpsRefY;
};

#endif /* GENOTYPE_H */

