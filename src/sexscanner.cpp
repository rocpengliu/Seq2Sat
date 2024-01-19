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
    
    int maxReadX = 0;//determine if the reads for x is deep or shallow, which is used for filter the low abundance read variants or not default 50:
    if (!seqMapX.empty()){
        for (auto & it : seqMapX) {
            if(it.second > maxReadX){
                maxReadX = it.second;
            }
        }
    }

    int maxReadY = 0;
    if (!seqMapY.empty()) {
        for (auto & it : seqMapY) {
            if (it.second > maxReadY) {
                maxReadY = it.second;
            }
        }
    }
    
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
        for (const auto & it : seqMapX) {
              for (int i = 0; i < it.first.length(); i++) {
                  baseFreqMapX[i][it.first[i]] += it.second;
              }
            
            if(maxReadX >= mOptions->mLocSnps.mLocSnpOptions.minReads4Filter) {
                if (it.second >= mOptions->mSex.minReadsSexVariant) {
                     mOptions->mSex.readsX += it.second;                    
                    readSeq = it.first.c_str();
                    readLength = it.first.length();
                    auto snpsMapX = SexScanner::doSimpleAlignment(mOptions, readSeq, readLength, target, targetLength);
                    for (auto & it2 : snpsMapX.first) {
                        mOptions->mSex.snpsRefX.insert(it2.first);
                    }
                    mOptions->mSex.seqVecX.emplace_back(std::make_tuple(it.first, it.second, snpsMapX.first));
                }
            } else {
                mOptions->mSex.readsX += it.second;
                for (int i = 0; i < it.first.length(); i++) {
                    baseFreqMapX[i][it.first[i]] += it.second;
                }

                readSeq = it.first.c_str();
                readLength = it.first.length();
                auto snpsMapX = SexScanner::doSimpleAlignment(mOptions, readSeq, readLength, target, targetLength);
                for (auto & it2 : snpsMapX.first) {
                    mOptions->mSex.snpsRefX.insert(it2.first);
                }
                mOptions->mSex.seqVecX.emplace_back(std::make_tuple(it.first, it.second, snpsMapX.first));
            }
        }

        std::sort(mOptions->mSex.seqVecX.begin(), mOptions->mSex.seqVecX.end(),
                [](const std::tuple<std::string, int, std::map<int, std::string>> &l,
                const std::tuple<std::string, int, std::map<int, std::string>> &r) {
                    return get<1>(l) > get<1>(r);
                });

        mOptions->mSex.haploTupX = std::make_tuple(get<0>(mOptions->mSex.seqVecX.front()), get<1>(mOptions->mSex.seqVecX.front()), false);
        if(mOptions->mSex.seqVecX.size() > 1){
            mOptions->mSex.haploTupX2 = std::make_tuple(get<0>(mOptions->mSex.seqVecX[1]), get<1>(mOptions->mSex.seqVecX[1]), false);
        }
    }

    if (!seqMapY.empty()) {
        target = mOptions->mSex.refY.mStr.c_str();
        targetLength = mOptions->mSex.refY.mStr.length();
        for (auto & it : seqMapY) {
            for (int i = 0; i < it.first.length(); i++) {
                        baseFreqMapY[i][it.first[i]] += it.second;
            }
            if (maxReadY >= mOptions->mLocSnps.mLocSnpOptions.minReads4Filter) { 
                if (it.second >= mOptions->mSex.minReadsSexVariant) {
                    mOptions->mSex.readsY += it.second;    
                    readSeq = it.first.c_str();
                    readLength = it.first.length();
                    auto snpsMapY = SexScanner::doSimpleAlignment(mOptions, readSeq, readLength, target, targetLength);
                    for (auto & it2 : snpsMapY.first) {
                        mOptions->mSex.snpsRefY.insert(it2.first);
                    }
                    mOptions->mSex.seqVecY.emplace_back(std::make_tuple(it.first, it.second, snpsMapY.first));
                }

            } else {
                mOptions->mSex.readsY += it.second;
                for (int i = 0; i < it.first.length(); i++) {
                    baseFreqMapY[i][it.first[i]] += it.second;
                }
                
                readSeq = it.first.c_str();
                readLength = it.first.length();
                auto snpsMapY = SexScanner::doSimpleAlignment(mOptions, readSeq, readLength, target, targetLength);
                for (auto & it2 : snpsMapY.first) {
                    mOptions->mSex.snpsRefY.insert(it2.first);
                }
                mOptions->mSex.seqVecY.emplace_back(std::make_tuple(it.first, it.second, snpsMapY.first));
            }
        }
        std::sort(mOptions->mSex.seqVecY.begin(), mOptions->mSex.seqVecY.end(),
                [](const std::tuple<std::string, int, std::map<int, std::string>> &l,
                const std::tuple<std::string, int, std::map<int, std::string>> &r) {
                    return get<1>(l) > get<1>(r);
                });
        mOptions->mSex.haploTupY = std::make_tuple(get<0>(mOptions->mSex.seqVecY.front()), get<1>(mOptions->mSex.seqVecY.front()), false);
    }
    
    if(get<1>(mOptions->mSex.haploTupX) != 0){
        mOptions->mSex.YXRatio = std::round(((double) get<1>(mOptions->mSex.haploTupY) / (double) get<1>(mOptions->mSex.haploTupX)) * 100.00) / 100.00;

        if (mOptions->mSex.YXRatio >= mOptions->mSex.YXRationCuttoff) {
            if (get<1>(mOptions->mSex.haploTupX) >= mOptions->mSex.minReadsSexAllele &&
                    get<1>(mOptions->mSex.haploTupY) >= mOptions->mSex.minReadsSexAllele) {
                mOptions->mSex.sexMF = "Male";
                get<2>(mOptions->mSex.haploTupX) = true;
                get<2>(mOptions->mSex.haploTupX2) = false;
                get<2>(mOptions->mSex.haploTupY) = true;
            } else {
                mOptions->mSex.sexMF = "Inconclusive";
            }
        } else {
            if (get<1>(mOptions->mSex.haploTupX) >= mOptions->mSex.minReadsSexAllele) {
                mOptions->mSex.sexMF = "Female";
                if (get<1>(mOptions->mSex.haploTupX2) != 0) {
                    
                    double haploRatio = std::round(((double) get<1>(mOptions->mSex.haploTupX2) / (double) (get<1>(mOptions->mSex.haploTupX) + get<1>(mOptions->mSex.haploTupX2)))* 100.0)/100.0;
                    
                    readSeq = get<0>(mOptions->mSex.haploTupX2).c_str();
                    readLength = get<0>(mOptions->mSex.haploTupX2).length();
                    target = get<0>(mOptions->mSex.haploTupX).c_str();
                    targetLength = get<0>(mOptions->mSex.haploTupX).length();
                
                    auto snpsMapXX = SexScanner::doSimpleAlignment(mOptions, readSeq, readLength, target, targetLength);
                    
                    if (snpsMapXX.second) {//no indels
    
                        mOptions->mSex.haploRatio = haploRatio;
                    
                        if (snpsMapXX.first.size() < 2) {
                            mOptions->mLocSnps.mLocSnpOptions.hmPer = mOptions->mLocSnps.mLocSnpOptions.hmPerL;
                            
                            if(haploRatio > mOptions->mLocSnps.mLocSnpOptions.hmPer){
                                get<2>(mOptions->mSex.haploTupX) = true;
                            } else if (abs(haploRatio - 0.5) <= mOptions->mLocSnps.mLocSnpOptions.htJetter) {
                                get<2>(mOptions->mSex.haploTupX) = true;
                                get<2>(mOptions->mSex.haploTupX2) = true;
                                mOptions->mSex.haploSnpsMap = snpsMapXX.first;
                                mOptions->mSex.snpsMapX2R = get<2>(mOptions->mSex.seqVecX[1]);
                                mOptions->mSex.haplotype = true;
                            } else {
                                //inconclusive haplotype;
                                mOptions->mSex.haploSnpsMap = snpsMapXX.first;
                                mOptions->mSex.snpsMapX2R = get<2>(mOptions->mSex.seqVecX[1]);
                            }

                        } else {
                            get<2>(mOptions->mSex.haploTupX) = true;
                            mOptions->mLocSnps.mLocSnpOptions.hmPer = mOptions->mLocSnps.mLocSnpOptions.hmPerH;
                            if(haploRatio <= mOptions->mLocSnps.mLocSnpOptions.hmPer){
                                get<2>(mOptions->mSex.haploTupX2) = true;
                                mOptions->mSex.haploSnpsMap = snpsMapXX.first;
                                mOptions->mSex.snpsMapX2R = get<2>(mOptions->mSex.seqVecX[1]);
                                mOptions->mSex.haplotype = true;
                            }
                        }
                    } else {// has indels; could have a bug if indels can also be haplo
                        mOptions->mLocSnps.mLocSnpOptions.hmPer = mOptions->mLocSnps.mLocSnpOptions.hmPerH;
                        get<2>(mOptions->mSex.haploTupX) = true;
                    }

                } else {
                    get<2>(mOptions->mSex.haploTupX) = true;
                }
            } else {
                mOptions->mSex.YXRatio = 0;
                mOptions->mSex.sexMF = "Inconclusive";
            }
        }
    } else {
        mOptions->mSex.YXRatio = 0;
        mOptions->mSex.sexMF = "Inconclusive";
    }
    
    if(!baseFreqMapY.empty() && mOptions->mSex.sexMF != "Inconclusive"){
        
        std::string sexStrY = get<0>(mOptions->mSex.haploTupY);
        for(int i = 0; i < baseFreqMapY.size(); i++){
            for(const auto & ii : baseFreqMapY[i]){
                if(ii.first == sexStrY[i]){
                    mOptions->mSex.baseErrorMapY[i] = static_cast<double>((mOptions->mSex.readsY - ii.second) * 100) / static_cast<double>(mOptions->mSex.readsY);
                    break;
                }
            }
        }
    }

    if (!baseFreqMapX.empty() && mOptions->mSex.sexMF != "Inconclusive") {
        
        if(!get<2>(mOptions->mSex.haploTupX2)){//male or female;

            std::string sexStrX = get<0>(mOptions->mSex.haploTupX);
            for (int i = 0; i < baseFreqMapX.size(); i++) {
                for (const auto & ii : baseFreqMapX[i]) {
                    if (ii.first == sexStrX[i]) {
                        mOptions->mSex.baseErrorMapX[i] = static_cast<double> ( (mOptions->mSex.readsX - ii.second) * 100) / static_cast<double> (mOptions->mSex.readsX);
                        break;
                    }
                }
            }
            
        } else {//female;
            std::string sexStrX = get<0>(mOptions->mSex.haploTupX);
            std::string sexStrX2 = get<0>(mOptions->mSex.haploTupX2);
            for (int i = 0; i < baseFreqMapX.size(); i++) {
                if(sexStrX[i] == sexStrX2[i]){
                    for (const auto & ii : baseFreqMapX[i]) {
                        if (ii.first == sexStrX[i]) {
                            mOptions->mSex.baseErrorMapX[i] = static_cast<double> ( (mOptions->mSex.readsX - ii.second) * 100) / static_cast<double> (mOptions->mSex.readsX);
                            break;
                        }
                    }
                } else {
                    int sum = 0;
                    for (const auto & ii : baseFreqMapX[i]) {
                        if (ii.first == sexStrX[i] || ii.first == sexStrX2[i]) {
                            sum += ii.second;
                        }
                    }
                    mOptions->mSex.baseErrorMapX[i] = static_cast<double> ( (mOptions->mSex.readsX - sum) * 100) / static_cast<double> (mOptions->mSex.readsX);
                }
            }
        }
    }
}

void SexScanner::report(Options * & mOptions) {

    std::string foutName = mOptions->prefix + "_sex_loc_id.txt";
    std::ofstream* fout = new std::ofstream();
    fout->open(foutName.c_str(), std::ofstream::out);

    if (!fout->is_open()) error_exit("Can not open output file: " + foutName);
    if (mOptions->verbose) loginfo("Starting to write sex identification loc file!");
    
    *fout << "#SexLoc\tSexAllele\tNumReads\tRatio\tPutativeSex\tHaplotypeRatio\tPutativeHaplotype\tAlleleSeq\tSNPs\tNote\n";
    
    if(mOptions->mSex.sexMF == "Male"){
        *fout << mOptions->mSex.sexMarker << "\t" << "Y" << "\t" << get<1>(mOptions->mSex.haploTupY) << "\t" << mOptions->mSex.YXRatio << "\t" << mOptions->mSex.sexMF << "\t" << 1 << "\t" << "Y" << "\t" << get<0>(mOptions->mSex.haploTupY) << "\t";
        
        if(get<2>(mOptions->mSex.seqVecY.front()).empty()){
            *fout << "NA";
        } else {
            for (const auto & it : get<2>(mOptions->mSex.seqVecY.front())) {
                *fout << it.first << mOptions->mSex.refY.mStr[it.first] << "|" << get<0>(mOptions->mSex.seqVecY.front())[it.first] << ";";
            }
        }
        *fout << "\t" << "haplotypeY" << "\n";
        
        *fout << mOptions->mSex.sexMarker << "\t" << "X" << "\t" << get<1>(mOptions->mSex.haploTupX) << "\t" << mOptions->mSex.YXRatio << "\t" << mOptions->mSex.sexMF << "\t" << 1 << "\t" << 'Y' << "\t" << get<0>(mOptions->mSex.haploTupX) << "\t";
        if (get<2>(mOptions->mSex.seqVecX.front()).empty()) {
            *fout << "NA";
        } else {
            for (const auto & it : get<2>(mOptions->mSex.seqVecX.front())) {
                *fout << it.first << mOptions->mSex.refX.mStr[it.first] << "|" << get<0>(mOptions->mSex.seqVecX.front())[it.first] << ";";
            }
        }
        *fout << "\t" << "haplotypeX" << "\n";

        for (int i = 1; i < mOptions->mSex.seqVecY.size(); i++) {
            *fout << mOptions->mSex.sexMarker << "\t" << "Y" << "\t" << get<1>(mOptions->mSex.seqVecY[i]) << "\t" << "NA" << "\t" << "NA" << "\t" << "NA" << "\t" << "N" << "\t" << get<0>(mOptions->mSex.seqVecY[i]) << "\t";
            if (get<2>(mOptions->mSex.seqVecY[i]).empty()) {
                *fout << "NA";
            } else {
                for (const auto & it : get<2>(mOptions->mSex.seqVecY[i])) {
                    *fout << it.first << mOptions->mSex.refY.mStr[it.first] << "|" << get<0>(mOptions->mSex.seqVecY[i])[it.first] << ";";
                }
            }
            *fout << "\t" << "seq_error" << "\n";
        }
       
        for (int i = 1; i < mOptions->mSex.seqVecX.size(); i++) {
            *fout << mOptions->mSex.sexMarker << "\t" << "X" << "\t" << get<1>(mOptions->mSex.seqVecX[i]) << "\t" << "NA" << "\t" << "NA" << "\t" << "NA" << "\t" << "N" << "\t" << get<0>(mOptions->mSex.seqVecX[i]) << "\t";
            if (get<2>(mOptions->mSex.seqVecX[i]).empty()) {
                *fout << "NA";
            } else {
                for (const auto & it : get<2>(mOptions->mSex.seqVecX[i])) {
                    *fout << it.first << mOptions->mSex.refX.mStr[it.first] << "|" << get<0>(mOptions->mSex.seqVecX[i])[it.first] << ";";
                }
            }
            *fout << "\t" << "seq_error" << "\n";
        }

    } else if(mOptions->mSex.sexMF == "Female"){
        
        if(get<2>(mOptions->mSex.haploTupX2)) {
            std::set<int> hapRefSet;
            if(!get<2>(mOptions->mSex.seqVecX.front()).empty()){
                for (const auto & it : get<2>(mOptions->mSex.seqVecX.front())) {
                    hapRefSet.insert(it.first);
                }
            }
            if (!mOptions->mSex.snpsMapX2R.empty()) {
                for (const auto & it : mOptions->mSex.snpsMapX2R) {
                    hapRefSet.insert(it.first);
                }
            }
            if(!mOptions->mSex.haploSnpsMap.empty()) {
                for (const auto & it : mOptions->mSex.haploSnpsMap) {
                    hapRefSet.insert(it.first);
                }
            }
            *fout << mOptions->mSex.sexMarker << "\t" << "X" << "\t" << get<1>(mOptions->mSex.haploTupX) << "\t" << mOptions->mSex.YXRatio << "\t" << mOptions->mSex.sexMF << "\t" << mOptions->mSex.haploRatio << "\t" << (mOptions->mSex.haplotype ? "Y" : "N") << "\t" << get<0>(mOptions->mSex.haploTupX) << "\t";
            if (hapRefSet.empty()) {
                *fout << "NA";
            } else {
                for (const auto & it : hapRefSet) {
                    *fout << it << mOptions->mSex.refX.mStr[it] << "|" << get<0>(mOptions->mSex.haploTupX)[it] << get<0>(mOptions->mSex.haploTupX2)[it] << ";";
                }
            }
            *fout << "\t" << "haplotypeX1" << "\n";

            *fout << mOptions->mSex.sexMarker << "\t" << "X2" << "\t" << get<1>(mOptions->mSex.haploTupX2) << "\t" << mOptions->mSex.YXRatio << "\t" << mOptions->mSex.sexMF << "\t" << mOptions->mSex.haploRatio << "\t" << (mOptions->mSex.haplotype ? "Y" : "N") << "\t" << get<0>(mOptions->mSex.haploTupX2) << "\t";
            if (hapRefSet.empty()) {
                *fout << "NA";
            } else {
                for (const auto & it : hapRefSet) {
                    *fout << it << mOptions->mSex.refX.mStr[it] << "|" << get<0>(mOptions->mSex.haploTupX)[it] << get<0>(mOptions->mSex.haploTupX2)[it] << ";";
                }
            }
            *fout << "\t" << "haplotypeX2" << "\n";

            for (int i = 2; i < mOptions->mSex.seqVecX.size(); i++) {
                *fout << mOptions->mSex.sexMarker << "\t" << "X" << "\t" << get<1>(mOptions->mSex.seqVecX[i]) << "\t" << "NA" << "\t" << "NA" << "\t" << "NA" << "\t" << "N" << "\t" << get<0>(mOptions->mSex.seqVecX[i]) << "\t";
                if (get<2>(mOptions->mSex.seqVecX[i]).empty()) {
                    *fout << "NA";
                } else {
                    for (const auto & it : get<2>(mOptions->mSex.seqVecX[i])) {
                        *fout << it.first << mOptions->mSex.refX.mStr[it.first] << "|" << get<0>(mOptions->mSex.seqVecX[i])[it.first] << ";";
                    }
                }
                *fout << "\t" << "seq_error" << "\n";
            }

        } else {

            *fout << mOptions->mSex.sexMarker << "\t" << "X" << "\t" << get<1>(mOptions->mSex.haploTupX) << "\t" << mOptions->mSex.YXRatio << "\t" << mOptions->mSex.sexMF << "\t" << mOptions->mSex.haploRatio << "\t" << "Y" << "\t" << get<0>(mOptions->mSex.haploTupX) << "\t";
            if (get<2>(mOptions->mSex.seqVecX.front()).empty()) {
                *fout << "NA";
            } else {
                for (const auto & it : get<2>(mOptions->mSex.seqVecX.front())) {
                    *fout << it.first << mOptions->mSex.refX.mStr[it.first] << "|" << get<0>(mOptions->mSex.haploTupX)[it.first] << ";";
                }
            }
            *fout << "\t" << "haplotypeX" << "\n";

            for (int i = 1; i < mOptions->mSex.seqVecX.size(); i++) {
                *fout << mOptions->mSex.sexMarker << "\t" << "X" << "\t" << get<1>(mOptions->mSex.seqVecX[i]) << "\t" << "NA" << "\t" << "NA" << "\t" << "NA" << "\t" << "N" << "\t" << get<0>(mOptions->mSex.seqVecX[i]) << "\t";
                if (get<2>(mOptions->mSex.seqVecX[i]).empty()) {
                    *fout << "NA";
                } else {
                    for (const auto & it : get<2>(mOptions->mSex.seqVecX[i])) {
                        *fout << it.first << mOptions->mSex.refX.mStr[it.first] << "|" << get<0>(mOptions->mSex.seqVecX[i])[it.first] << ";";
                    }
                }
                *fout << "\t" << "seq_error" << "\n";
            }
        }
        
    } else {

        if (!mOptions->mSex.seqVecY.empty()) {
            *fout << mOptions->mSex.sexMarker << "\t" << "Y" << "\t" << get<1>(mOptions->mSex.haploTupY) << "\t" << mOptions->mSex.YXRatio << "\t" << mOptions->mSex.sexMF << "\t" << "NA" << "\t" << "N" << "\t" << get<0>(mOptions->mSex.haploTupY) << "\t";

            if (get<2>(mOptions->mSex.seqVecY.front()).empty()) {
                *fout << "NA";
            } else {
                for (const auto & it : get<2>(mOptions->mSex.seqVecY.front())) {
                    *fout << it.first << mOptions->mSex.refY.mStr[it.first] << "|" << get<0>(mOptions->mSex.seqVecY.front())[it.first] << ";";
                }
            }
            *fout << "\t" << "haplotypeY" << "\n";

            for (int i = 1; i < mOptions->mSex.seqVecY.size(); i++) {
                *fout << mOptions->mSex.sexMarker << "\t" << "Y" << "\t" << get<1>(mOptions->mSex.seqVecY[i]) << "\t" << "NA" << "\t" << "NA" << "\t" << "NA" << "\t" << "N" << "\t" << get<0>(mOptions->mSex.seqVecY[i]) << "\t";
                if (get<2>(mOptions->mSex.seqVecY[i]).empty()) {
                    *fout << "NA";
                } else {
                    for (const auto & it : mOptions->mSex.snpsRefY) {
                        *fout << it << mOptions->mSex.refY.mStr[it] << "|" << get<0>(mOptions->mSex.seqVecY[i])[it] << ";";
                    }
                }
                *fout << "\t" << "seq_error" << "\n";
            }
        }
        
        if(!mOptions->mSex.seqVecX.empty()){

            *fout << mOptions->mSex.sexMarker << "\t" << "X" << "\t" << get<1>(mOptions->mSex.haploTupX) << "\t" << mOptions->mSex.YXRatio << "\t" << mOptions->mSex.sexMF << "\t" << mOptions->mSex.haploRatio << "\t" << "N" << "\t" << get<0>(mOptions->mSex.haploTupX) << "\t";
            if (get<2>(mOptions->mSex.seqVecX.front()).empty()) {
                *fout << "NA";
            } else {
                for (const auto & it : get<2>(mOptions->mSex.seqVecX.front())) {
                    *fout << it.first << mOptions->mSex.refX.mStr[it.first] << "|" << get<0>(mOptions->mSex.seqVecX.front())[it.first] << ";";
                }
            }
            *fout << "\t" << "haplotypeX" << "\n";

            for (int i = 1; i < mOptions->mSex.seqVecX.size(); i++) {
                *fout << mOptions->mSex.sexMarker << "\t" << "X" << "\t" << get<1>(mOptions->mSex.seqVecX[i]) << "\t" << "NA" << "\t" << "NA" << "\t" << "NA" << "\t" << "N" << "\t" << get<0>(mOptions->mSex.seqVecX[i]) << "\t";
                if (get<2>(mOptions->mSex.seqVecX[i]).empty()) {
                    *fout << "NA";
                } else {
                    for (const auto & it : get<2>(mOptions->mSex.seqVecX[i])) {
                        *fout << it.first << mOptions->mSex.refX.mStr[it.first] << "|" << get<0>(mOptions->mSex.seqVecX[i])[it.first] << ";";
                    }
                }
                *fout << "\t" << "seq_error" << "\n";
            }
        }
    }

    fout->flush();
    fout->clear();
    fout->close();
    
    foutName = mOptions->prefix + "_sex_error_rate.txt";
    fout->open(foutName.c_str(), std::ofstream::out);

    if (!fout->is_open()) error_exit("Can not open output file: " + foutName);
    if (mOptions->verbose) loginfo("Starting to write sex error rate file!");
    *fout << "#Locus\tErrorRate\tTotalEffectiveReads\n";
    if (!mOptions->mSex.baseErrorMapY.empty()) {
        *fout << "Y" << "\t";
        for (const auto & it : mOptions->mSex.baseErrorMapY) {
            if (it.first == mOptions->mSex.baseErrorMapY.rbegin()->first) {
                *fout << it.second;
            } else {
                *fout << it.second << ";";
            }
        }
        *fout << "\t" << mOptions->mSex.readsY << "\n";
    }

    if (!mOptions->mSex.baseErrorMapX.empty()) {
        *fout << "X" << "\t";
        for (const auto & it : mOptions->mSex.baseErrorMapX) {
            if (it.first == mOptions->mSex.baseErrorMapX.rbegin()->first) {
                *fout << it.second;
            } else {
                *fout << it.second << ";";
            }
        }
        *fout << "\t" << mOptions->mSex.readsX << "\n";
    }

    fout->flush();
    fout->clear();
    fout->close();
    
    if (fout) {
        delete fout;
        fout = nullptr;
    }
}

std::pair<std::map<int, std::string>, bool> SexScanner::doSimpleAlignment(Options * & mOptions, const char* & qData, int qLength, const char* & tData, int tLength) {
    EdlibAlignResult result = edlibAlign(qData, qLength, tData, tLength,
            edlibNewAlignConfig(mOptions->mLocVars.locVarOptions.maxScorePrimer,
            mOptions->mEdOptions.modeCode,
            mOptions->mEdOptions.alignTask,
            NULL, 0));

    std::map<int, std::string> snpsMap;
    std::set<int> indelSet;
    if (result.status == EDLIB_STATUS_OK) {
        for (int i = 0; i < result.alignmentLength; i++) {
            auto cur = result.alignment[i];
            if (cur == EDLIB_EDOP_MATCH) {

            } else if (cur == EDLIB_EDOP_MISMATCH) {
                snpsMap[i] = tData[i];
            } else if (cur == EDLIB_EDOP_INSERT) {
                indelSet.insert(i);
            } else if (cur == EDLIB_EDOP_DELETE) {
                indelSet.insert(i);
            }
        }

    }
    edlibFreeAlignResult(result);
    if (indelSet.empty()) {
        return std::make_pair(snpsMap, true);
    } else {
        if(qLength == tLength){
            return std::make_pair(snpsMap, true);
        } else {
            return std::make_pair(snpsMap, false);
        }
    }
}
