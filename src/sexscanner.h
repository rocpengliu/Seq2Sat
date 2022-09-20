#ifndef SEXSCANNER_H
#define SEXSCANNER_H

#include "options.h"
#include "read.h"

class SexScanner {
public:
    SexScanner();
    SexScanner(const SexScanner& orig);
    virtual ~SexScanner();
    
public:
    static int sexScan(Read* mergedPE, int misMatchesP = 2, int misMatchesA = 2);
    static int sexScan(Read* r1, Read* r2, int misMatchesP = 2, int misMatchesA = 2);
    
private:

};

#endif /* SEXSCANNER_H */

