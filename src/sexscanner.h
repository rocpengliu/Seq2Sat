#ifndef SEXSCANNER_H
#define SEXSCANNER_H

#include <map>
#include <utility>
#include <cmath>
#include <tuple>
#include "options.h"
#include "read.h"
#include "genotype.h"
#include "util.h"
#include "edlib.h"
#include "common.h"
#include "editdistance.h"

class SexScanner {
public:
    SexScanner(Options* opt);
    virtual ~SexScanner();
    
public:
    std::pair<bool, char> sexScan(Read* r1);
    bool sexScan(Read* r1, Read* r2);
    inline std::map<std::string, std::map<std::string, int>> getSexLoc() const {return tmpSexMap;};
    void static merge(std::vector<std::map<std::string, std::map<std::string, int>>> & totalSexLocVec, Options * & mOptions);
    void static report(Options * & mOptions);
    static std::pair<std::map<int, std::string>, bool>  doSimpleAlignment(Options * & mOptions, const char* & qData, int qLength, const char* & tData, int tLength);

private:
    std::tuple<int, int, bool> doPrimerAlignment(const char* & qData, int qLen, const std::string & qName,
                     const char* & tData, int tLen, const std::string & tName, 
                     bool printAlignment);
    
private:
    Options* mOptions;
    const char* qData;
    int qLen;
    const char* tData;
    int tLen;
    std::map<std::string, std::map<std::string, int>> tmpSexMap;
};

#endif /* SEXSCANNER_H */

