#include <set>

#include "genotype.h"

LocVar::LocVar() {
    this->name = "";
    this->fp = Sequence("");
    this->rp = Sequence("");
    this->ff = Sequence("");
    this->rf = Sequence("");
    this->repuitAll = Sequence("");
    this->repuit = Sequence("");
    this->repuit2 = Sequence("");
    this->nRep = 0;
    this->edCutoff = 0;
    this->mra = Sequence("");
    this->locSeq = Sequence("");
    this->effectiveSeq = Sequence("");
    this->locLen = 0;
    this->effectiveLen = 0;
    this->mraName = "";
    this->snpsMapff.clear();
    this->snpsMaprf.clear();
    this->refSnpsSetffMap.clear();
    this->refSnpsSetrfMap.clear();
    this->mraBase = "";
    this->totalReads = 0;
}

LocVar::LocVar(std::string & read, std::pair<size_t, int> & mraPosLenReadF, std::string & ssr){
    this->ff = Sequence(read.substr(0, mraPosLenReadF.first));
    this->rf = Sequence(read.substr(mraPosLenReadF.first + mraPosLenReadF.second));
    this->mra = Sequence(read.substr(mraPosLenReadF.first, mraPosLenReadF.second));
    this->repuit = Sequence(ssr);
    this->effectiveSeq = Sequence(read);
    this->effectiveLen = read.length();
}

LocVar::LocVar(std::string & read, std::size_t & mraL, std::size_t & mraR, std::string & ssr){
    this->ff = Sequence(read.substr(0, mraL));
    this->rf = Sequence(read.substr(mraR));
    this->mra = Sequence(read.substr(mraL, mraR - mraL));
    this->effectiveSeq = Sequence(read);
    this->effectiveLen = read.length();
    this->repuit = Sequence(ssr);
}

LocVar::LocVar(std::string & name, std::string & ff, std::string & rf, std::string & repuit, std::string & mra){
    this->name = name;
    this->ff = Sequence(ff);
    this->rf = Sequence(rf);
    this->repuit = Sequence(repuit);
    this->mra = Sequence(mra);
    this->effectiveSeq = Sequence(ff + mra + rf);
    this->effectiveLen = ff.length() + mra.length() + rf.length();
}

void LocVar::print(){
    std::string msg = "name: " + name +  " ff: " + ff.mStr + " mra: " + mra.mStr + " rf: " + rf.mStr;
    cCout(msg, 'r');
}

void LocVar::printSnpMap(){
    std::string msg = "name: " + name +  " ff: " + ff.mStr + " mra: " + mra.mStr + " rf: " + rf.mStr;
    
    msg += "\nffSnpMap:\n";
    cCout("000000000000000000000000", 'y');
    for(const auto & it : refSnpsSetffMap){
        for(const auto & it2 : it.second){
            //if(it2 > 0){
            cCout("0000000000000000000000001111111111111111111", 'y');
              msg += it.first + " : " + std::to_string(it2) + "\n";  
              cCout("0000000000000000000000002222222222222222", 'y');
            //}
        }
    }
    
    cCout("0000000000000000000000003333333333333333333", 'y');
    msg += "rfSnpMap:\n";
    for(auto & it : refSnpsSetrfMap){
        for(const auto & it2 : it.second) {
            //if (it2 > 0) {
            cCout("000000000000000000000000444444444444444444444", 'y');
                msg += it.first + " : " + std::to_string(it2) + "\n";
            //}
                cCout("000000000000000000000000555555555555555555", 'y');
        }
    }
    cCout(msg, 'g');
    cCout("0000000000000000000000006666666666666", 'y');
}

LocVar::LocVar(std::string name, Sequence fp, Sequence rp, Sequence ff, Sequence rf, Sequence repuit, int nRep, unsigned int edCutoff, Sequence mra){
    this->name = name;
    this->fp = fp;
    this->rp = rp;
    this->ff = ff;
    this->rf = rf;
    this->repuit = repuit;
    this->nRep = nRep;
    this->edCutoff = edCutoff;
    this->mra = mra;
    this->effectiveSeq = ff.mStr + mra.mStr + rf.mStr;
    this->locSeq = fp.mStr + ff.mStr + mra.mStr + rf.mStr + rp.mStr;
    this->locLen = locSeq.length();
    this->effectiveLen = effectiveSeq.length();
    
}

Variance::Variance() {
    subMap.clear();
    delMap.clear();
    insMap.clear();
}

void Variance::cleanVar() {
    subMap.clear();
    delMap.clear();
    insMap.clear();
}


Genotype::Genotype() {
    numReads = 0;
    ssrF = "";
    ssrR = "";
}

void Genotype::print() {
    cCout("Num Reads: " + std::to_string(numReads), 'r');
    baseLocVar.print();
}

UnitedGenotype::UnitedGenotype(){
    marker = "";
    repuit = Sequence("");
    mra = Sequence("");
    effectiveSeq = Sequence("");
    //seqName = "";
    effectiveLen = 0;
    mraName = "";
    mraSVec.clear();
    genoSVec.clear();
    numReadsVec.clear();
    locVec.clear();
    effectiveSeqVec.clear();
    seqNameVec.clear();
}

LocSnp::LocSnp(){
    this->name = "";
    this->fp = Sequence("");
    this->rp = Sequence("");
    this->ref = Sequence("");
    this->snpPosSet.clear();
    this->numReads = 0;
    this->snpsMap.clear();
}

void LocSnp::print(){
    std::string msg = "";
    msg = "name: " + name + " -> " + std::to_string(numReads) + "\n";
    msg += fp.mStr +  " :" + std::to_string(fp.mStr.length()) + "\n";
    msg += rp.mStr +  " :" + std::to_string(rp.mStr.length()) + "\n";
    msg += ref.mStr + " :" + std::to_string(ref.mStr.length()) + "\n";
    cCout(msg, 'r');
    msg = "";
    for(const auto & it : snpsMap){
        msg += std::to_string(it.first) + " : " + it.second.first.mStr + "|" + it.second.second.mStr + "\n";
    }
    cCout(msg, 'b');
}

std::string LocSnp::getGenotype(){
    if(snpsMap.empty()){
        return "ref";
    }
    std::string msg = "";
    for(const auto & it : snpsMap){
        msg += std::to_string(it.first) + "(" + it.second.first.mStr + "|" + it.second.second.mStr + ")";
    }
    return msg;
}

Sex::Sex(){
    this->sexMarker = "";
    this->primerF = Sequence("");
    this->primerR = Sequence("");
    this->refX = Sequence("");
    this->refY = Sequence("");
    this->readsX = 0;
    this->readsY = 0;
    this->minTotalReadsX = 10;
    this->minTotalReadsY = 10;
    this->minReadsX = 5;
    this->minReadsY = 5;
    this->mismatchesPF = 1;
    this->mismatchesPR = 1;
    this->mismatchesRX = 2;
    this->mismatchesRY = 2;
    this->lengthEqual = false;
    this->YXRatio = 0;
    this->YXRationCuttoff = 0.0001;
    this->sexMF = "Inconclusive";
    this->seqVecX.clear();
    this->seqVecY.clear();
}

void Sex::print(){
    std::cout << "name: " << sexMarker << "; primerF: " << primerF.mStr << "; primerR: " << primerR.mStr << 
            ";\n refx: " << refX.mStr << "; refy: " << refY.mStr << 
            ";\n readsx: " << readsX << "; readsy: " << readsY << "\n"; 
}

std::string Sex::getFullRefX(){
    return std::string(primerF.mStr + refX.mStr + primerR.mStr);
}

std::string Sex::getFullRefY(){
    return std::string(primerF.mStr + refY.mStr + primerR.mStr);
}