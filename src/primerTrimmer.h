#ifndef PRIMERTRIMMER_H
#define PRIMERTRIMMER_H

#include <string>
#include "options.h"
#include "read.h"

class PrimerTrimmer {
public:
    PrimerTrimmer();
    ~PrimerTrimmer();
    
public:
    static bool trimByPrimerSeq(Read* r1, Read* r2, int misMatches = 4);
    static bool trimByPrimerSeq(Read* mergedPE, int misMatches = 4);
    
private:

};

#endif /* PRIMERTRIMMER_H */

