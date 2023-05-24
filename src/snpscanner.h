#ifndef SNPSCANNER_H
#define SNPSCANNER_H

#include <string>
#include <algorithm>
#include <vector>
#include <functional>
#include <map>
#include <queue>
#include <set>
#include <utility>
#include <sstream>
#include <algorithm>

#include "options.h"
#include "read.h"
#include "genotype.h"
#include "util.h"
#include "edlib.h"
#include "common.h"
#include "editdistance.h"
#include "sequence.h"

using namespace std;

class SnpScanner {
public:
    SnpScanner(Options * opt);
    SnpScanner(const SnpScanner& orig);
    virtual ~SnpScanner();
    
    std::string scanVar(Read* & r1);
    bool scanVar(Read* & r1, Read* & r2);
    //inline std::map<std::string, std::map<std::string, LocSnp>> getSubGenotypeMap(){return subGenotypeMap;};
    inline std::map<std::string, std::map<std::string, uint32>> getSubSeqsMap(){return subSeqsMap;};
    static std::map<std::string, std::map<std::string, LocSnp>> merge(Options * & mOptions, std::vector<std::map<std::string, std::map < std::string, uint32>>> & totalSnpSeqMapVec);
    static std::map<int, std::pair<Sequence, Sequence>> doAlignment(Options * & mOptions, std::string readName, const char* & qData, int qLength, std::string & targetName, const char* & tData, int tLength);
    static void doScanVariance(Options * & mOptions, EdlibAlignResult & result, Variance & variance, const char* & qData, const char* & tData, const int position);
    static void printVariance(Options * & mOptions, EdlibAlignResult & result, Variance & variance,
            const char* & qData, const std::string & qName, const char* & tData, const std::string & tName, const int position);
    inline Sex getSexLoc(){return *tmpSex;};
    
private:
    std::tuple<int, int, bool> doPrimerAlignment(const char* & qData, int qLength, const std::string & qName,
                     const char* & tData, int tLength, const std::string & tName, 
                     bool printAlignment);
    
private:
    Options* mOptions;
    //std::map<std::string, std::map<std::string, LocSnp>> subGenotypeMap;
    std::map<std::string, std::map<std::string, uint32>> subSeqsMap;
    LocSnp* locSnpIt;
    const char* fpData;
    int fpLength;
    const char* rpData;
    int rpLength;
    const char* target;
    int targetLength;
    const char* readSeq;
    int readLength;
    std::string readName;
    std::stringstream ss;
    Sex* tmpSex;
    std::string returnedlocus;
};

#endif /* SNPSCANNER_H */

