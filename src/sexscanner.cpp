#include "sexscanner.h"

SexScanner::SexScanner(Options* opt) {
    mOptions = opt;
    qData = nullptr;
    qLen = 0;
    tData = nullptr;
    tLen = 0;
    tmpSexMap.clear();
}

SexScanner::~SexScanner() {
}

std::pair<bool, char> SexScanner::sexScan(Read* r1) {
    
    std::pair<bool, char> rep{false, 'U'};
    if (mOptions->mSex.sexMarker.empty() || (r1->length() < (mOptions->mSex.primerF.length() + mOptions->mSex.primerR.length() + std::min(mOptions->mSex.refX.length(), mOptions->mSex.refY.length())))) {
        return rep;
    }
    
    int trimPosF = 0;
    bool goRP = false;
    
    int fpMismatches = (int) edit_distance(mOptions->mSex.primerF.mStr, r1->mSeq.mStr.substr(0, mOptions->mSex.primerF.length()));
    if (fpMismatches <= mOptions->mSex.mismatchesPF) {
        trimPosF = mOptions->mSex.primerF.length();
        goRP = true;
    } else {
        qData = mOptions->mSex.primerF.mStr.c_str();
        qLen = mOptions->mSex.primerF.length();
        tData = r1->mSeq.mStr.c_str();
        tLen = r1->length();
        auto endBoolF = doPrimerAlignment(qData, qLen, mOptions->mSex.sexMarker, tData, tLen, r1->mName, true);
        if (get<2>(endBoolF) && (get<1>(endBoolF) <= r1->length())) {
            fpMismatches = get<0>(endBoolF);
            if (fpMismatches <= mOptions->mSex.mismatchesPF) {
                if ((get<1>(endBoolF) + mOptions->mSex.primerR.length() + std::min(mOptions->mSex.refX.length(), mOptions->mSex.refY.length())) <= r1->length()) {
                    trimPosF = get<1>(endBoolF);
                    goRP = true;
                } else {
                    goRP = false;
                }
            } else {
                goRP = false;
            }
        } else {
            goRP = false;
        }
    }
    
    if(!goRP){
        return rep;
    }
    
    r1->trimFront(trimPosF);
    int trimLen = 0;
    int rpMismatches = (int) edit_distance(mOptions->mSex.primerR.mStr, r1->mSeq.mStr.substr(r1->mSeq.length() - mOptions->mSex.primerR.length()));
    bool goSex = false;
    if (rpMismatches <= mOptions->mSex.mismatchesPR) {
        trimLen = r1->length() - mOptions->mSex.primerR.length();
        goSex = true;
    } else {
        qData = mOptions->mSex.primerR.mStr.c_str();
        qLen = mOptions->mSex.primerR.length();
        tData = r1->mSeq.mStr.c_str();
        tLen = r1->length();
        auto endBoolR = doPrimerAlignment(qData, qLen, mOptions->mSex.sexMarker, tData, tLen, r1->mName, true);

        if (get<2>(endBoolR) && (get<1>(endBoolR) <= r1->length())) {
            rpMismatches = get<0>(endBoolR);
            if (rpMismatches <= mOptions->mSex.mismatchesPR) {
                if (get<1>(endBoolR) <= r1->length()) {
                    trimLen = get<1>(endBoolR) - mOptions->mSex.primerR.length();
                    goSex = true;
                } else {
                    goSex = false;
                }
            } else {
                goSex = false;
            }
        } else {
            goSex = false;
        }
    }

    if(!goSex || trimLen <= 0 || (trimLen < std::min(mOptions->mSex.refX.length(), mOptions->mSex.refY.length())) || 
            (trimLen > std::max(mOptions->mSex.refX.length(), mOptions->mSex.refY.length()))){
        return rep;
    }
    
    r1->resize(trimLen);
    
    if (mOptions->mSex.lengthEqual) {
        unsigned int edx = edit_distance(mOptions->mSex.refX.mStr, r1->mSeq.mStr);
        unsigned int edy = edit_distance(mOptions->mSex.refY.mStr, r1->mSeq.mStr);
        
        if(edx == edy){
            return rep;
        }

        unsigned int edmin = std::min(edx, edy);
        if (edmin == edx) {
            if (edx > mOptions->mSex.mismatchesRX) {
                return rep;
            } else {
                tmpSexMap["X"][r1->mSeq.mStr]++;
                rep = std::make_pair(true, 'X');
            }
        } else {
            if (edy > mOptions->mSex.mismatchesRY) {
                return rep;
            } else {
                tmpSexMap["Y"][r1->mSeq.mStr]++;
                rep = std::make_pair(true, 'Y');
            }
        }
    } else {
        if (r1->mSeq.length() == mOptions->mSex.refX.length()) {
            unsigned int ed = edit_distance(mOptions->mSex.refX.mStr, r1->mSeq.mStr);
            if (ed > mOptions->mSex.mismatchesRX) {
                return rep;
            } else {
                tmpSexMap["X"][r1->mSeq.mStr]++;
                rep = std::make_pair(true, 'X');
            }
        } else if (r1->mSeq.length() == mOptions->mSex.refY.length()) {
            unsigned int ed = edit_distance(mOptions->mSex.refY.mStr, r1->mSeq.mStr);
            if (ed > mOptions->mSex.mismatchesRY) {
                return rep;
            } else {
                tmpSexMap["Y"][r1->mSeq.mStr]++;
                rep = std::make_pair(true, 'Y');
            }
        } else {
            return rep;
        }
    }
    
    return rep;
}

std::tuple<int, int, bool> SexScanner::doPrimerAlignment(const char* & qData, int qLength, const std::string & qName,
        const char* & tData, int tLength, const std::string & tName, bool printAlignment) {

    EdlibAlignResult result = edlibAlign(qData, qLength, tData, tLength,
            edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));

    if (result.status == EDLIB_STATUS_OK) {

        std::set<int> snpsSet;
        std::set<int> indelSet;
        for (int i = 0; i < result.alignmentLength; i++) {
            auto cur = result.alignment[i];
            if (cur == EDLIB_EDOP_MATCH) {

            } else if (cur == EDLIB_EDOP_MISMATCH) {
                snpsSet.insert(i);
            } else if (cur == EDLIB_EDOP_INSERT) {
                indelSet.insert(i);
            } else if (cur == EDLIB_EDOP_DELETE) {
                indelSet.insert(i);
            }
        }

        int endPos = *(result.endLocations) + 1;
        edlibFreeAlignResult(result);
        if (indelSet.empty() && (snpsSet.size() <= mOptions->mLocVars.locVarOptions.maxMismatchesPSeq)) {
            return std::make_tuple(snpsSet.size(), endPos, true);
        } else {
            return std::make_tuple(0, 0, false);
        }
    } else {
        edlibFreeAlignResult(result);
        return std::make_tuple(0, 0, false);
    }
}

void SexScanner::merge(std::vector<std::map<std::string, std::map<std::string, int>>> & totalSexLocVec, Options * & mOptions) {
    std::map<std::string, int> seqMapX;
    std::map<std::string, int> seqMapY;

    for (auto & it : totalSexLocVec) {
        for (const auto & it2 : it["X"]) {
            mOptions->isPaired() ? (seqMapX[it2.first] += it2.second * 2) : (seqMapX[it2.first] += it2.second);
        }

        for (const auto & it2 : it["Y"]) {
            mOptions->isPaired() ? (seqMapY[it2.first] += it2.second * 2) : seqMapY[it2.first] += it2.second;
        }
    }

    totalSexLocVec.clear();
    totalSexLocVec.shrink_to_fit();
    
        //Calculate error rate
    std::map<int, std::map<char, int>> baseFreqMapX;
    std::map<int, std::map<char, int>> baseFreqMapY;

    const char* target;
    int targetLength;
    const char* readSeq;
    int readLength;  

    if (!seqMapX.empty()) {
        target = mOptions->mSex.refX.mStr.c_str();
        targetLength = mOptions->mSex.refX.length();
        mOptions->mSex.seqVarVecX.reserve(seqMapY.size());
        for (const auto & it : seqMapX) {
            mOptions->mSex.totReadsX += it.second;
            if (it.second > mOptions->mSex.maxReadsX) {
                mOptions->mSex.maxReadsX = it.second;
            }

            readSeq = it.first.c_str();
            readLength = it.first.length();
            auto mapPair1 = doAlignment2(mOptions, "read", readSeq, readLength, mOptions->mSex.sexMarker, target, targetLength);
            SeqVar tmpSeqVar;
            tmpSeqVar.seq = it.first;
            tmpSeqVar.numReads = it.second;

            if (mapPair1.first) {
                tmpSeqVar.snpSet = mapPair1.second;
                for (int i = 0; i < it.first.length(); i++) {
                    baseFreqMapX[i][it.first[i]] += it.second;
                }
            } else {
                tmpSeqVar.indel = true;
            }
            mOptions->mSex.seqVarVecX.push_back(tmpSeqVar);
        }
        std::sort(mOptions->mSex.seqVarVecX.begin(), mOptions->mSex.seqVarVecX.end(), [](const SeqVar &L, const SeqVar &R){ return L.numReads > R.numReads;});
        
        int keepInt = 0;
        for(int i = 0; i < mOptions->mSex.seqVarVecX.size(); i++){
            if(i < 2 || mOptions->mSex.seqVarVecX.at(i).numReads >= mOptions->mSex.minReadsSexVariant){
                mOptions->mSex.totSnpSetX.insert(mOptions->mSex.seqVarVecX.at(i).snpSet.begin(), mOptions->mSex.seqVarVecX.at(i).snpSet.end());
                keepInt++;
            } else {
                break;
            }
        }
        mOptions->mSex.seqVarVecX.erase(mOptions->mSex.seqVarVecX.begin() + keepInt, mOptions->mSex.seqVarVecX.end());
    }
    
    if (!seqMapY.empty()){
        target = mOptions->mSex.refY.mStr.c_str();
        targetLength = mOptions->mSex.refY.length();
        mOptions->mSex.seqVarVecY.reserve(seqMapY.size());
        for (auto & it : seqMapY) {
            mOptions->mSex.totReadsY += it.second;
            if(it.second > mOptions->mSex.maxReadsY){
                mOptions->mSex.maxReadsY = it.second;
            }

            readSeq = it.first.c_str();
            readLength = it.first.length();
            auto mapPairY = doAlignment2(mOptions, "read", readSeq, readLength, mOptions->mSex.sexMarker, target, targetLength);
            SeqVar tmpSeqVar;
            tmpSeqVar.seq = it.first;
            tmpSeqVar.numReads = it.second;

            if (mapPairY.first) {
                tmpSeqVar.snpSet = mapPairY.second;
                for (int i = 0; i < it.first.length(); i++) {
                    baseFreqMapY[i][it.first[i]] += it.second;
                }
            } else {
                tmpSeqVar.indel = true;
            }
            mOptions->mSex.seqVarVecY.push_back(tmpSeqVar);
        }
        mOptions->mSex.seqVarVecY.shrink_to_fit();
        std::sort(mOptions->mSex.seqVarVecY.begin(), mOptions->mSex.seqVarVecY.end(), [](const SeqVar & L, const SeqVar & R) {
            return L.numReads > R.numReads;
        });

        int keepInt = 0;
        for(int i = 0; i < mOptions->mSex.seqVarVecY.size(); i++){
            if(i < 2 || mOptions->mSex.seqVarVecY.at(i).numReads >= mOptions->mSex.minReadsSexVariant){
                mOptions->mSex.totSnpSetY.insert(mOptions->mSex.seqVarVecY.at(i).snpSet.begin(), mOptions->mSex.seqVarVecY.at(i).snpSet.end());
                keepInt++;
            } else {
                break;
            }
        }
        mOptions->mSex.seqVarVecY.erase(mOptions->mSex.seqVarVecY.begin() + keepInt, mOptions->mSex.seqVarVecY.end());
    }
    
    if(mOptions->mSex.maxReadsY >= mOptions->mSex.minReadsSexAllele){//male; or inconclusive;
        if(mOptions->mSex.maxReadsX >= mOptions->mSex.minReadsSexAllele){
            mOptions->mSex.YXRatio = getPer(mOptions->mSex.getHaploVar('y', 0).numReads, mOptions->mSex.getHaploVar('x', 0).numReads, false);
            if (mOptions->mSex.YXRatio >= mOptions->mSex.YXRationCuttoff){
                mOptions->mSex.sexMF = "Male";
                mOptions->mSex.haploStr = "homo";
            } else {
                mOptions->mSex.sexMF = "Inconclusive";
            }
        } else {
            mOptions->mSex.sexMF = "Inconclusive";
        }
    } else {//female;
        if (mOptions->mSex.maxReadsX >= mOptions->mSex.minReadsSexAllele){
            mOptions->mSex.sexMF = "Female";
            if (mOptions->mSex.getHaploVar('x', 1).numReads > 0){
                mOptions->mSex.haploRatio = getPer(mOptions->mSex.getHaploVar('x', 0).numReads, 
                (mOptions->mSex.getHaploVar('x', 0).numReads + mOptions->mSex.getHaploVar('x', 1).numReads),
                 false);

                auto mapPair2 = edit_distance2(mOptions->mSex.getHaploVar('x', 0).seq, mOptions->mSex.getHaploVar('x', 1).seq);

                if(mapPair2.first){//no indel
                    if (mapPair2.second.size() == 1){//1 snp;
                        if (mOptions->mSex.haploRatio >= mOptions->mLocSnps.mLocSnpOptions.hmPerL) {
                            mOptions->mSex.haploStr = "homo";
                        } else if (abs(mOptions->mSex.haploRatio - 0.5) <= mOptions->mLocSnps.mLocSnpOptions.htJetter) {
                            if( mOptions->mSex.getHaploVar('x', 1).numReads >= mOptions->mSex.minReadsSexVariant) {
                                mOptions->mSex.haplotype = true;
                                mOptions->mSex.haploStr = "heter";
                                mOptions->mSex.hyploSnpSetX = mapPair2.second;
                                mOptions->mSex.totSnpSetX.insert(mOptions->mSex.hyploSnpSetX.begin(), mOptions->mSex.hyploSnpSetX.end());
                            } else {
                                
                            }
                            
                        } else {
                            //mOptions->mSex.haploStr = "inconclusive";
                        }
                    } else {// >= 2 snps;
                        if (mOptions->mSex.haploRatio >= mOptions->mLocSnps.mLocSnpOptions.hmPerH){
                            mOptions->mSex.haploStr = "homo";
                        } else {
                            if( mOptions->mSex.getHaploVar('x', 1).numReads >= mOptions->mSex.minReadsSexVariant) {
                                mOptions->mSex.haplotype = true;
                                mOptions->mSex.haploStr = "heter";
                                mOptions->mSex.hyploSnpSetX = mapPair2.second;
                                mOptions->mSex.totSnpSetX.insert(mOptions->mSex.hyploSnpSetX.begin(), mOptions->mSex.hyploSnpSetX.end());
                            } else {
                                
                            }
                        }
                    }
                } else {// indel;
                    if (mOptions->mSex.haploRatio >= mOptions->mLocSnps.mLocSnpOptions.hmPerH){
                        mOptions->mSex.haploStr = "homo";
                    } else {
                        if( mOptions->mSex.getHaploVar('x', 1).numReads >= mOptions->mSex.minReadsSexVariant){
                            mOptions->mSex.haploStr = "heter";
                        } else {
                            
                        }
                    }
                    mOptions->mSex.haploIndel = true;
                }
            } else {
                mOptions->mSex.haploStr = "homo";
            }
        } else {
            mOptions->mSex.sexMF = "Inconclusive";
        }
    }
    
    if(!mOptions->mSex.getHaploVar('y', 0).indel){
        for (int i = 0; i < mOptions->mSex.getHaploVar('y', 0).seq.length(); i++) {
            mOptions->mSex.baseErrorMapY[i] = getPer((mOptions->mSex.totReadsY - baseFreqMapY[i][mOptions->mSex.getHaploVar('y', 0).seq[i]]),
                                                     mOptions->mSex.totReadsY);
        }
    }

    if (!mOptions->mSex.getHaploVar('x', 0).indel){
        if (mOptions->mSex.haploStr == "homo" && !mOptions->mSex.getHaploVar('x', 0).indel){
            for (int i = 0; i < mOptions->mSex.getHaploVar('x', 0).seq.length(); i++) {
                mOptions->mSex.baseErrorMapX[i] = getPer((mOptions->mSex.totReadsX - baseFreqMapX[i][mOptions->mSex.getHaploVar('x', 0).seq[i]]),
                                                        mOptions->mSex.totReadsX);
            }
        } else if (mOptions->mSex.haploStr == "heter") {
            if (!mOptions->mSex.getHaploVar('x', 0).indel && !mOptions->mSex.getHaploVar('x', 1).indel){
                for (int i = 0; i < mOptions->mSex.getHaploVar('x', 0).seq.length(); i++){
                    int tot = baseFreqMapX[i][mOptions->mSex.getHaploVar('x', 0).seq[i]] + baseFreqMapX[i][mOptions->mSex.getHaploVar('x', 1).seq[i]];
                    mOptions->mSex.baseErrorMapX[i] = getPer((mOptions->mSex.totReadsX - tot),
                                                             mOptions->mSex.totReadsX);
                }
            }
        } else {
            
        }
    }
}

//void SexScanner::merge2(std::vector<std::map<std::string, std::map<std::string, int>>> & totalSexLocVec, Options * & mOptions) {
//    std::map<std::string, int> seqMapX;
//    std::map<std::string, int> seqMapY;
//
//    for (auto & it : totalSexLocVec) {
//        for (const auto & it2 : it["X"]) {
//            mOptions->isPaired() ? (seqMapX[it2.first] += it2.second * 2) : (seqMapX[it2.first] += it2.second);
//        }
//
//        for (const auto & it2 : it["Y"]) {
//            mOptions->isPaired() ? (seqMapY[it2.first] += it2.second * 2) : seqMapY[it2.first] += it2.second;
//        }
//    }
//
//    totalSexLocVec.clear();
//    totalSexLocVec.shrink_to_fit();
//    
//        //Calculate error rate
//    std::map<int, std::map<char, int>> baseFreqMapX;
//    std::map<int, std::map<char, int>> baseFreqMapY;
//
//    const char* target;
//    int targetLength;
//    const char* readSeq;
//    int readLength;  
//
//    if (!seqMapX.empty()) {
//        target = mOptions->mSex.refX.mStr.c_str();
//        targetLength = mOptions->mSex.refX.length();
//        mOptions->mSex.seqVarVecX.reserve(seqMapY.size());
//        for (const auto & it : seqMapX) {
//            mOptions->mSex.totReadsX += it.second;
//            if (it.second > mOptions->mSex.maxReadsX) {
//                mOptions->mSex.maxReadsX = it.second;
//            }
//
//            readSeq = it.first.c_str();
//            readLength = it.first.length();
//            auto mapPair1 = doAlignment2(mOptions, "read", readSeq, readLength, mOptions->mSex.sexMarker, target, targetLength);
//            cCout("sexxxxxxxxxxxxxxxx111111111111111", mapPair1.second.size(), 'b');
//            SeqVar tmpSeqVar;
//            tmpSeqVar.seq = it.first;
//            tmpSeqVar.numReads = it.second;
//
//            if (mapPair1.first) {
//                tmpSeqVar.snpSet = mapPair1.second;
//                for (int i = 0; i < it.first.length(); i++) {
//                    baseFreqMapX[i][it.first[i]] += it.second;
//                }
//            } else {
//                tmpSeqVar.indel = true;
//            }
//            mOptions->mSex.seqVarVecX.push_back(tmpSeqVar);
//        }
//        std::sort(mOptions->mSex.seqVarVecX.begin(), mOptions->mSex.seqVarVecX.end(), [](const SeqVar &L, const SeqVar &R){ return L.numReads > R.numReads;});
//        
//        int keepInt = 0;
//        for(int i = 0; i < mOptions->mSex.seqVarVecX.size(); i++){
//            if(i < 2 || mOptions->mSex.seqVarVecX.at(i).numReads >= mOptions->mSex.minReadsSexVariant){
//                mOptions->mSex.totSnpSetX.insert(mOptions->mSex.seqVarVecX.at(i).snpSet.begin(), mOptions->mSex.seqVarVecX.at(i).snpSet.end());
//                keepInt++;
//            } else {
//                break;
//            }
//        }
//        mOptions->mSex.seqVarVecX.erase(mOptions->mSex.seqVarVecX.begin() + keepInt, mOptions->mSex.seqVarVecX.end());
//        
//        for(const auto & it : mOptions->mSex.seqVarVecX){
//            cCout(it.seq, it.numReads, 'g');
//        }
//    }
//    
//    if (!seqMapY.empty()){
//        target = mOptions->mSex.refY.mStr.c_str();
//        targetLength = mOptions->mSex.refY.length();
//        mOptions->mSex.seqVarVecY.reserve(seqMapY.size());
//        for (auto & it : seqMapY) {
//            mOptions->mSex.totReadsY += it.second;
//            if(it.second > mOptions->mSex.maxReadsY){
//                mOptions->mSex.maxReadsY = it.second;
//            }
//
//            readSeq = it.first.c_str();
//            readLength = it.first.length();
//            auto mapPairY = doAlignment2(mOptions, "read", readSeq, readLength, mOptions->mSex.sexMarker, target, targetLength);
//            cCout("sexYYYYYYYYYYYYYYYYYYYYY111111111111111", mapPairY.second.size(), 'b');
//            SeqVar tmpSeqVar;
//            tmpSeqVar.seq = it.first;
//            tmpSeqVar.numReads = it.second;
//
//            if (mapPairY.first) {
//                tmpSeqVar.snpSet = mapPairY.second;
//                for (int i = 0; i < it.first.length(); i++) {
//                    baseFreqMapY[i][it.first[i]] += it.second;
//                }
//            } else {
//                tmpSeqVar.indel = true;
//            }
//            mOptions->mSex.seqVarVecY.push_back(tmpSeqVar);
//        }
//        mOptions->mSex.seqVarVecY.shrink_to_fit();
//        std::sort(mOptions->mSex.seqVarVecY.begin(), mOptions->mSex.seqVarVecY.end(), [](const SeqVar & L, const SeqVar & R) {
//            return L.numReads > R.numReads;
//        });
//
//        int keepInt = 0;
//        for(int i = 0; i < mOptions->mSex.seqVarVecY.size(); i++){
//            if(i == 0 || mOptions->mSex.seqVarVecY.at(i).numReads >= mOptions->mSex.minReadsSexVariant){
//                mOptions->mSex.totSnpSetY.insert(mOptions->mSex.seqVarVecY.at(i).snpSet.begin(), mOptions->mSex.seqVarVecY.at(i).snpSet.end());
//                keepInt++;
//            } else {
//                break;
//            }
//        }
//        mOptions->mSex.seqVarVecY.erase(mOptions->mSex.seqVarVecY.begin() + keepInt, mOptions->mSex.seqVarVecY.end());
//    }
//    
//    if(mOptions->mSex.maxReadsY >= mOptions->mSex.minReadsSexAllele){//male; or inconclusive;
//        if(mOptions->mSex.maxReadsX >= mOptions->mSex.minReadsSexAllele){
//            mOptions->mSex.YXRatio = getPer(mOptions->mSex.getHaploVar('y', 0).numReads, mOptions->mSex.getHaploVar('x', 0).numReads, false);
//            if (mOptions->mSex.YXRatio >= mOptions->mSex.YXRationCuttoff){
//                mOptions->mSex.sexMF = "Male";
//                mOptions->mSex.haploStr = "homo";
//            } else {
//                mOptions->mSex.sexMF = "Inconclusive";
//            }
//        } else {
//            mOptions->mSex.sexMF = "Inconclusive";
//        }
//    } else {//female;
//        if (mOptions->mSex.maxReadsX >= mOptions->mSex.minReadsSexAllele){
//            mOptions->mSex.sexMF = "Female";
//            if (mOptions->mSex.getHaploVar('x', 1).numReads > 0){
//                mOptions->mSex.haploRatio = getPer(mOptions->mSex.getHaploVar('x', 0).numReads, 
//                (mOptions->mSex.getHaploVar('x', 0).numReads + mOptions->mSex.getHaploVar('x', 1).numReads),
//                 false);
//                const char* readSeq2 = mOptions->mSex.getHaploVar('x', 0).seq.c_str();
//                int readLength2 = mOptions->mSex.getHaploVar('x', 0).seq.length();
//                const char* target2 = mOptions->mSex.getHaploVar('x', 1).seq.c_str();
//                int targetLength2 = mOptions->mSex.getHaploVar('x', 1).seq.length();
//                
//                cCout(mOptions->mSex.getHaploVar('x', 0).seq, mOptions->mSex.getHaploVar('x', 1).seq, 'b');
//                std::pair<bool, std::set<int>> mapPair2 = doAlignment2(mOptions, "x2", readSeq2, readLength2, "x1", target2, targetLength2);
//                cCout(mOptions->mSex.haploRatio, mapPair2.second.size(), 'b');
//                if(mapPair2.first){//no indel
//                    if (mapPair2.second.size() == 1){//1 snp;
//                        if (mOptions->mSex.haploRatio >= mOptions->mLocSnps.mLocSnpOptions.hmPerL) {
//                            mOptions->mSex.haploStr = "homo";
//                        } else if (abs(mOptions->mSex.haploRatio - 0.5) <= mOptions->mLocSnps.mLocSnpOptions.htJetter) {
//                            if( mOptions->mSex.getHaploVar('x', 1).numReads >= mOptions->mSex.minReadsSexVariant) {
//                                mOptions->mSex.haplotype = true;
//                                mOptions->mSex.haploStr = "heter";
//                                mOptions->mSex.hyploSnpSetX = mapPair2.second;
//                                mOptions->mSex.totSnpSetX.insert(mOptions->mSex.hyploSnpSetX.begin(), mOptions->mSex.hyploSnpSetX.end());
//                            } else {
//                                
//                            }
//                            
//                        } else {
//                            //mOptions->mSex.haploStr = "inconclusive";
//                        }
//                    } else {// >= 2 snps;
//                        if (mOptions->mSex.haploRatio >= mOptions->mLocSnps.mLocSnpOptions.hmPerH){
//                            mOptions->mSex.haploStr = "homo";
//                        } else {
//                            if( mOptions->mSex.getHaploVar('x', 1).numReads >= mOptions->mSex.minReadsSexVariant) {
//                                mOptions->mSex.haplotype = true;
//                                mOptions->mSex.haploStr = "heter";
//                                mOptions->mSex.hyploSnpSetX = mapPair2.second;
//                                mOptions->mSex.totSnpSetX.insert(mOptions->mSex.hyploSnpSetX.begin(), mOptions->mSex.hyploSnpSetX.end());
//                            } else {
//                                
//                            }
//                        }
//                    }
//                } else {// indel;
//                    if (mOptions->mSex.haploRatio >= mOptions->mLocSnps.mLocSnpOptions.hmPerH){
//                        mOptions->mSex.haploStr = "homo";
//                    } else {
//                        if( mOptions->mSex.getHaploVar('x', 1).numReads >= mOptions->mSex.minReadsSexVariant){
//                            mOptions->mSex.haploStr = "heter";
//                        } else {
//                            
//                        }
//                    }
//                    mOptions->mSex.haploIndel = true;
//                }
//            } else {
//                mOptions->mSex.haploStr = "homo";
//            }
//        } else {
//            mOptions->mSex.sexMF = "Inconclusive";
//        }
//    }
//    
//    if(!mOptions->mSex.getHaploVar('y', 0).indel){
//        for (int i = 0; i < mOptions->mSex.getHaploVar('y', 0).seq.length(); i++) {
//            mOptions->mSex.baseErrorMapY[i] = getPer((mOptions->mSex.totReadsY - baseFreqMapY[i][mOptions->mSex.getHaploVar('y', 0).seq[i]]),
//                                                     mOptions->mSex.totReadsY);
//        }
//    }
//
//    if (!mOptions->mSex.getHaploVar('x', 0).indel){
//        if (mOptions->mSex.haploStr == "homo" && !mOptions->mSex.getHaploVar('x', 0).indel){
//            for (int i = 0; i < mOptions->mSex.getHaploVar('x', 0).seq.length(); i++) {
//                mOptions->mSex.baseErrorMapX[i] = getPer((mOptions->mSex.totReadsX - baseFreqMapX[i][mOptions->mSex.getHaploVar('x', 0).seq[i]]),
//                                                        mOptions->mSex.totReadsX);
//            }
//        } else if (mOptions->mSex.haploStr == "heter") {
//            if (!mOptions->mSex.getHaploVar('x', 0).indel && !mOptions->mSex.getHaploVar('x', 1).indel){
//                for (int i = 0; i < mOptions->mSex.getHaploVar('x', 0).seq.length(); i++){
//                    int tot = baseFreqMapX[i][mOptions->mSex.getHaploVar('x', 0).seq[i]] + baseFreqMapX[i][mOptions->mSex.getHaploVar('x', 1).seq[i]];
//                    mOptions->mSex.baseErrorMapX[i] = getPer((mOptions->mSex.totReadsX - tot),
//                                                             mOptions->mSex.totReadsX);
//                }
//            }
//        } else {
//            
//        }
//    }
//}

void SexScanner::report(Options *& mOptions) {

    std::string foutName = mOptions->prefix + "_sex_loc_id.txt";
    std::ofstream *fout = new std::ofstream();
    fout->open(foutName.c_str(), std::ofstream::out);

    if (!fout->is_open())
        error_exit("Can not open output file: " + foutName);
    if (mOptions->verbose)
        loginfo("Starting to write sex identification loc file!");

    *fout << "#SexLoc\tSexAllele\tNumReads\tRatio\tPutativeSex\tHaplotypeRatio\tPutativeHaplotype\tAlleleSeq\tSNPs\tNote\n";

    if (mOptions->mSex.sexMF == "Male") {
        *fout << mOptions->mSex.sexMarker << "\t"//loci
                << "Y" << "\t"//sex allele
                << mOptions->mSex.getHaploVar('y', 0).numReads << "\t"//reads
                << mOptions->mSex.YXRatio << "\t"//ratio
                << mOptions->mSex.sexMF << "\t" //mf
                << 1 << "\t"//MF
                << "Y" << "\t"
                << mOptions->mSex.getHaploVar('y', 0).seq << "\t"
                << (mOptions->mSex.getSnpStr('y', 0).empty() ? "NA" : mOptions->mSex.getSnpStr('y', 0)) << "\t"
                << "haplotypeY" << "\n";

        *fout << mOptions->mSex.sexMarker << "\t"
                << "X" << "\t"
                << mOptions->mSex.getHaploVar('x', 0).numReads << "\t"
                << mOptions->mSex.YXRatio << "\t"
                << mOptions->mSex.sexMF << "\t" 
                << 1 << "\t"
                << "X" << "\t"
                << mOptions->mSex.getHaploVar('x', 0).seq << "\t"
                << (mOptions->mSex.getSnpStr('x', 0).empty() ? "NA" : mOptions->mSex.getSnpStr('x', 0)) << "\t"
                << "haplotypeX" << "\n";

        for (int i = 1; i < mOptions->mSex.seqVarVecY.size(); i++) {
            *fout << mOptions->mSex.sexMarker << "\t"
                    << "Y" << "\t" 
                    << mOptions->mSex.getHaploVar('y', i).numReads << "\t"
                    << "NA" << "\t"//ratio
                    << "NA" << "\t"//MF
                    << "NA" << "\t"//
                    << "N" << "\t"//
                    << mOptions->mSex.getHaploVar('y', i).seq << "\t"
                    << (mOptions->mSex.getSnpStr('y', i).empty() ? "NA" : mOptions->mSex.getSnpStr('y', i)) << "\t"
                    << "seq_error" << "\n";
        }

        for (int i = 1; i < mOptions->mSex.seqVarVecX.size(); i++) {
            *fout << mOptions->mSex.sexMarker << "\t"
                    << "X" << "\t" 
                    << mOptions->mSex.getHaploVar('x', i).numReads << "\t"
                    << "NA" << "\t"//ratio
                    << "NA" << "\t"//MF
                    << "NA" << "\t"//
                    << "N" << "\t"//
                    << mOptions->mSex.getHaploVar('x', i).seq << "\t"
                    << (mOptions->mSex.getSnpStr('x', i).empty() ? "NA" : mOptions->mSex.getSnpStr('x', i)) << "\t"
                    << "seq_error" << "\n";
        }
        
    } else if (mOptions->mSex.sexMF == "Female") {
        
        if(mOptions->mSex.haplotype){
            *fout << mOptions->mSex.sexMarker << "\t"
                    << "X" << "\t"
                    << mOptions->mSex.getHaploVar('x', 0).numReads << "\t"
                    << mOptions->mSex.YXRatio << "\t"
                    << mOptions->mSex.sexMF << "\t"
                    << mOptions->mSex.haploRatio << "\t"
                    << (mOptions->mSex.haplotype ? "Y" : "N") << "\t"
                    << mOptions->mSex.getHaploVar('x', 0).seq << "\t"
                    << (mOptions->mSex.getSnpStr('x', 0).empty() ? "NA" : mOptions->mSex.getSnpStr('x', 0)) << "\t"
                    << "haplotypeX" << "\n";

            *fout << mOptions->mSex.sexMarker << "\t"
                    << "X" << "\t"
                    << mOptions->mSex.getHaploVar('x', 1).numReads << "\t"
                    << mOptions->mSex.YXRatio << "\t"
                    << mOptions->mSex.sexMF << "\t"
                    << mOptions->mSex.haploRatio << "\t"
                    << (mOptions->mSex.haplotype ? "Y" : "N") << "\t"
                    << mOptions->mSex.getHaploVar('x', 1).seq << "\t"
                    << (mOptions->mSex.getSnpStr('x', 1).empty() ? "NA" : mOptions->mSex.getSnpStr('x', 1)) << "\t"
                    << "haplotypeX" << "\n";
            
            for (int i = 2; i < mOptions->mSex.seqVarVecX.size(); i++) {
                *fout << mOptions->mSex.sexMarker << "\t"
                        << "X" << "\t"
                        << mOptions->mSex.getHaploVar('x', i).numReads << "\t"
                        << "NA" << "\t"//ratio
                        << "NA" << "\t"//MF
                        << "NA" << "\t"//
                        << "N" << "\t"//
                        << mOptions->mSex.getHaploVar('x', i).seq << "\t"
                        << (mOptions->mSex.getSnpStr('x', i).empty() ? "NA" : mOptions->mSex.getSnpStr('x', i)) << "\t"
                        << "seq_error" << "\n";
            }
            
        } else {
            *fout << mOptions->mSex.sexMarker << "\t"
                << "X" << "\t"
                << mOptions->mSex.getHaploVar('x', 0).numReads << "\t"
                << mOptions->mSex.YXRatio << "\t"
                << mOptions->mSex.sexMF << "\t" 
                << 1 << "\t"
                << "X" << "\t"
                << mOptions->mSex.getHaploVar('x', 0).seq << "\t"
                << (mOptions->mSex.getSnpStr('x', 0).empty() ? "NA" : mOptions->mSex.getSnpStr('x', 0)) << "\t"
                << "haplotypeX" << "\n";

            for (int i = 1; i < mOptions->mSex.seqVarVecX.size(); i++) {
                *fout << mOptions->mSex.sexMarker << "\t"
                        << "X" << "\t"
                        << mOptions->mSex.getHaploVar('x', i).numReads << "\t"
                        << "NA" << "\t"//ratio
                        << "NA" << "\t"//MF
                        << "NA" << "\t"//
                        << "N" << "\t"//
                        << mOptions->mSex.getHaploVar('x', i).seq << "\t"
                        << (mOptions->mSex.getSnpStr('x', i).empty() ? "NA" : mOptions->mSex.getSnpStr('x', i)) << "\t"
                        << "seq_error" << "\n";
            }
            
        }

    } else {
                
        if (!mOptions->mSex.seqVarVecY.empty()) {
            *fout << mOptions->mSex.sexMarker << "\t"//loci
                    << "Y" << "\t"//sex allele
                    << mOptions->mSex.getHaploVar('y', 0).numReads << "\t"//reads
                    << mOptions->mSex.YXRatio << "\t"//ratio
                    << mOptions->mSex.sexMF << "\t" //mf
                    << mOptions->mSex.haploRatio << "\t"//haplotype ratio;
                    << "N" << "\t"//tPutativeHaplotype
                    << mOptions->mSex.getHaploVar('y', 0).seq << "\t"//tAlleleSeq
                    << (mOptions->mSex.getSnpStr('y', 0).empty() ? "NA" : mOptions->mSex.getSnpStr('y', 0)) << "\t"//snps;
                    << "haplotypeY" << "\n";//note;

            for (int i = 1; i < mOptions->mSex.seqVarVecY.size(); i++) {
                *fout << mOptions->mSex.sexMarker << "\t"
                        << "Y" << "\t"
                        << mOptions->mSex.getHaploVar('y', i).numReads << "\t"
                        << "NA" << "\t"//ratio
                        << "NA" << "\t"//MF
                        << "NA" << "\t"//
                        << "N" << "\t"//
                        << mOptions->mSex.getHaploVar('y', i).seq << "\t"
                        << (mOptions->mSex.getSnpStr('y', i).empty() ? "NA" : mOptions->mSex.getSnpStr('y', i)) << "\t"
                        << "seq_error" << "\n";
            }
        }

        if (!mOptions->mSex.seqVarVecX.empty()) {
            *fout << mOptions->mSex.sexMarker << "\t"
                    << "X" << "\t"
                    << mOptions->mSex.getHaploVar('x', 0).numReads << "\t"
                    << mOptions->mSex.YXRatio << "\t"
                    << mOptions->mSex.sexMF << "\t"
                    << mOptions->mSex.haploRatio << "\t"
                    << "N" << "\t"
                    << mOptions->mSex.getHaploVar('x', 0).seq << "\t"
                    << (mOptions->mSex.getSnpStr('x', 0).empty() ? "NA" : mOptions->mSex.getSnpStr('x', 0)) << "\t"
                    << "haplotypeX" << "\n";

            for (int i = 1; i < mOptions->mSex.seqVarVecX.size(); i++) {
                *fout << mOptions->mSex.sexMarker << "\t"
                        << "X" << "\t"
                        << mOptions->mSex.getHaploVar('x', i).numReads << "\t"
                        << "NA" << "\t"//ratio
                        << "NA" << "\t"//MF
                        << "NA" << "\t"//
                        << "N" << "\t"//
                        << mOptions->mSex.getHaploVar('x', i).seq << "\t"
                        << (mOptions->mSex.getSnpStr('x', i).empty() ? "NA" : mOptions->mSex.getSnpStr('x', i)) << "\t"
                        << "seq_error" << "\n";
            }
        }
    }

    fout->flush();
    fout->clear();
    fout->close();

    foutName = mOptions->prefix + "_sex_error_rate.txt";
    fout->open(foutName.c_str(), std::ofstream::out);

    if (!fout->is_open())
        error_exit("Can not open output file: " + foutName);
    if (mOptions->verbose)
        loginfo("Starting to write sex error rate file!");
    *fout << "#Locus\tErrorRate\tTotalReads\n";
    if (!mOptions->mSex.baseErrorMapY.empty()) {
        *fout << "Y" << "\t";
        for (const auto &it : mOptions->mSex.baseErrorMapY) {
            if (it.first == mOptions->mSex.baseErrorMapY.rbegin()->first) {
                *fout << it.second;
            } else {
                *fout << it.second << ";";
            }
        }
        *fout << "\t" << mOptions->mSex.totReadsY << "\n";
    }

    if (!mOptions->mSex.baseErrorMapX.empty()) {
        *fout << "X" << "\t";
        for (const auto &it : mOptions->mSex.baseErrorMapX) {
            if (it.first == mOptions->mSex.baseErrorMapX.rbegin()->first) {
                *fout << it.second;
            } else {
                *fout << it.second << ";";
            }
        }
        *fout << "\t" << mOptions->mSex.totReadsX << "\n";
    }

    fout->flush();
    fout->clear();
    fout->close();

    if (fout) {
        delete fout;
        fout = nullptr;
    }
}

//std::pair<std::map<int, std::string>, bool> SexScanner::doSimpleAlignment(Options * & mOptions, const char* & qData, int qLength, const char* & tData, int tLength) {
//    EdlibAlignResult result = edlibAlign(qData, qLength, tData, tLength,
//            edlibNewAlignConfig(mOptions->mLocVars.locVarOptions.maxScorePrimer,
//            mOptions->mEdOptions.modeCode,
//            mOptions->mEdOptions.alignTask,
//            NULL, 0));
//
//    std::map<int, std::string> snpsMap;
//    std::set<int> indelSet;
//    if (result.status == EDLIB_STATUS_OK) {
//        for (int i = 0; i < result.alignmentLength; i++) {
//            auto cur = result.alignment[i];
//            if (cur == EDLIB_EDOP_MATCH) {
//
//            } else if (cur == EDLIB_EDOP_MISMATCH) {
//                snpsMap[i] = tData[i];
//            } else if (cur == EDLIB_EDOP_INSERT) {
//                indelSet.insert(i);
//            } else if (cur == EDLIB_EDOP_DELETE) {
//                indelSet.insert(i);
//            }
//        }
//
//    }
//    edlibFreeAlignResult(result);
//    if (indelSet.empty()) {
//        return std::make_pair(snpsMap, true);
//    } else {
//        if(qLength == tLength){
//            return std::make_pair(snpsMap, true);
//        } else {
//            return std::make_pair(snpsMap, false);
//        }
//    }
//}

std::pair<bool, std::set<int>> SexScanner::doAlignment2(Options * & mOptions, std::string readName, const char* & qData, int qLength, std::string targetName, const char* & tData, int tLength) {
    EdlibAlignResult result = edlibAlign(qData, qLength, tData, tLength,
            edlibNewAlignConfig(mOptions->mLocVars.locVarOptions.maxScorePrimer,
            EDLIB_MODE_NW,
            mOptions->mEdOptions.alignTask,
            NULL, 0));

    std::pair<bool, std::set<int>> snpsSetPair;
    
    if (result.status == EDLIB_STATUS_OK) {
        bool snps = true;
        for (int i = 0; i < result.alignmentLength; i++) {
            auto cur = result.alignment[i];
            if (cur == EDLIB_EDOP_MATCH) {

            } else if (cur == EDLIB_EDOP_MISMATCH) {
                snpsSetPair.second.insert(i);
            } else if (cur == EDLIB_EDOP_INSERT) {
                snps = snp && false;
            } else if (cur == EDLIB_EDOP_DELETE) {
                snps = snp && false;
            }
        }
        snpsSetPair.first = snps;
    }
    edlibFreeAlignResult(result);
    return snpsSetPair;
}