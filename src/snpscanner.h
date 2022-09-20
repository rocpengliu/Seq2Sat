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

using namespace std;

class SnpScanner {
public:
    SnpScanner(Options * opt);
    SnpScanner(const SnpScanner& orig);
    virtual ~SnpScanner();
    
    bool scanVar(Read* & r1);
    bool scanVar(Read* & r1, Read* & r2);
    inline std::map<std::string, std::map<std::string, LocSnp>> getSubGenotypeMap(){return subGenotypeMap;};
    static std::map<std::string, std::map<std::string, LocSnp>> merge(Options * & mOptions, std::vector<std::map<std::string, std::map < std::string, LocSnp>>> & totalGenotypeSnpMapVec);
    inline Sex getSexLoc(){return *tmpSex;};
    
private:
    std::map<int, std::pair<Sequence, Sequence>> doAlignment(const char* & qData, int qLength, const char* & tData, int tLength);
    void doScanVariance(EdlibAlignResult & result, Variance & variance, const char* & qData, const char* & tData, const int position);

    void printVariance(EdlibAlignResult & result, Variance & variance,
            const char* & qData, const std::string & qName, const char* & tData, const std::string & tName, const int position);
    
private:
    Options* mOptions;
    std::map<std::string, std::map<std::string, LocSnp>> subGenotypeMap;
    LocSnp locSnpIt;
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
};

#endif /* SNPSCANNER_H */

