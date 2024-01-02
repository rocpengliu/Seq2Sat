#include <valarray>

#include "snpscanner.h"

SnpScanner::SnpScanner(Options* opt) {
    mOptions = opt;
    //subGenotypeMap.clear();
    subSeqsMap.clear();
    target = NULL;
    targetLength = 0;
    readSeq = NULL;
    readLength = 0;
    readName = "";
    fpData = NULL;
    fpLength = 0;
    rpData = NULL;
    rpLength = 0;
    ss.str();
}

SnpScanner::SnpScanner(const SnpScanner& orig) {
}

SnpScanner::~SnpScanner() {
}

bool SnpScanner::scanVar(Read* & r1, Read* & r2) {
    return true;
}

std::string SnpScanner::scanVar(Read* & r1) {
    ss.str("");
    returnedlocus.clear();
    readSeq = r1->mSeq.mStr.c_str();
    readLength = r1->mSeq.mStr.length();
    readName = r1->mName;
    std::map<std::string, std::pair<int, int>> locMap;
    for (auto & it : mOptions->mLocSnps.refLocMap) {
        uint32 fpMismatches = edit_distance(it.second.fp.mStr, r1->mSeq.mStr.substr(0, it.second.fp.length()));
        if(fpMismatches <= mOptions->mLocSnps.mLocSnpOptions.maxMismatchesPSeq){
            uint32 rpMismatches = edit_distance(it.second.rp.mStr, r1->mSeq.mStr.substr(r1->mSeq.length() - it.second.rp.length()));
            if(rpMismatches <= mOptions->mLocSnps.mLocSnpOptions.maxMismatchesPSeq){
                locMap[it.second.name] = std::make_pair((fpMismatches + rpMismatches), readLength - it.second.fp.mStr.length() - it.second.rp.mStr.length());
            } else {
                fpData = it.second.rp.mStr.c_str();
                fpLength = it.second.rp.length();
                auto endBool = doPrimerAlignment(fpData, fpLength, it.second.name, readSeq, readLength, r1->mName, true);
                if (get<2>(endBool) && get<1>(endBool) <= readLength) {
                    locMap[it.second.name] = std::make_pair((fpMismatches + get<0>(endBool)), get<1>(endBool) - it.second.fp.mStr.length() - it.second.rp.mStr.length());
                }
            }
        }
    }
     
    if(locMap.empty()){
        return("");
    }

    std::string locName = "";
    if (locMap.size() == 1) {
        locName = locMap.begin()->first;
        //if(mOptions->debug) cCout("single value: " + locName, 'r');
    } else {
        std::vector<int> seqScoreVec;
        for (const auto & it : locMap) {
            seqScoreVec.push_back(it.second.first);
        }
        auto minValue = *std::min_element(seqScoreVec.begin(), seqScoreVec.end());
        //warning, what if there are multiple identical values
        seqScoreVec.clear();

        for (const auto & it : locMap) {
            if (it.second.first == minValue) {
                locName = it.first;
                ss.str();
                break; //warning, what if there are multiple identical values
            }
        }
    }
    
    locSnpIt = &(mOptions->mLocSnps.refLocMap[locName]);
    //locSnpIt->print();

    if (readLength < (locSnpIt->fp.length() + locSnpIt->ft.length() + locSnpIt->ref.length() + locSnpIt->rt.length() + locSnpIt->rp.length())) {
        return returnedlocus;
    } else if (readLength > (locSnpIt->fp.length() + locSnpIt->ft.length() + locSnpIt->ref.length() + locSnpIt->rt.length() + locSnpIt->rp.length())) {
        r1->resize(locSnpIt->fp.length() + locSnpIt->ft.length() + locSnpIt->ref.length() + locSnpIt->rt.length() + locSnpIt->rp.length());
    }

    //reads length is the same as ref + fp + rp now;
    r1->trimFront(locSnpIt->fp.length() + locSnpIt->ft.length());
    r1->resize(locSnpIt->ref.length());
    
    locMap.clear();

    if (mOptions->debug) cCout("detected marker: " + locSnpIt->name, 'r');

    subSeqsMap[locSnpIt->name][r1->mSeq.mStr]++;
    returnedlocus = locName;

    if (mOptions->mEdOptions.printRes) {
        cCout(ss.str(), 'g');
    }

    ss.str();
    return returnedlocus;
}

std::pair<bool,std::map<int, std::pair<Sequence, Sequence>>> SnpScanner::doAlignment(Options * & mOptions, std::string readName, const char* & qData, int qLength, std::string targetName, const char* & tData, int tLength) {
    EdlibAlignResult result = edlibAlign(qData, qLength, tData, tLength,
            edlibNewAlignConfig(mOptions->mLocVars.locVarOptions.maxScorePrimer,
            EDLIB_MODE_NW,
            mOptions->mEdOptions.alignTask,
            NULL, 0));
    
    std::pair<bool, std::map<int, std::pair<Sequence, Sequence>>> snpsMapPair;
    
    if (result.status == EDLIB_STATUS_OK) {
        if (mOptions->mEdOptions.printRes) {
            Variance variance;
            doScanVariance(mOptions, result, variance, qData, tData, *(result.endLocations));
            if (mOptions->mEdOptions.printRes) {
                //printVariance(result, variance, qData, readName, tData, targetName, *(result.endLocations));
            }

            variance.cleanVar();
        }
        bool snps = true;
        for (int i = 0; i < result.alignmentLength; i++) {
            auto cur = result.alignment[i];
            if (cur == EDLIB_EDOP_MATCH) {
                
            } else if(cur == EDLIB_EDOP_MISMATCH) {
                std::string s1(1, tData[i]);
                std::string s2(1, qData[i]);
                snpsMapPair.second[i] = std::make_pair(Sequence(s1), Sequence(s2));
            } else if(cur == EDLIB_EDOP_INSERT) {
                snps = snp && false;
            } else if(cur == EDLIB_EDOP_DELETE){
                snps = snp && false;
            }
        }
        
        snpsMapPair.first = snps;
    }
    edlibFreeAlignResult(result);
    return snpsMapPair;
}

void SnpScanner::doScanVariance(Options * & mOptions, EdlibAlignResult & result, Variance & variance,
        const char* & qData, const char* & tData, const int position) {

    int ti = -1, qi = -1;
    std::string inStr, deStr, snpStr;

    if (result.alignmentLength < 2) {
        //if (mOptions->mLocVars.locVarOptions.printRes) ss << "You must have at least 2 bp!\n";
        return;
    }

    if (mOptions->mEdOptions.modeCode == EDLIB_MODE_HW) {
        ti = position;
        for (int i = 0; i < result.alignmentLength; i++) {
            if (result.alignment[i] != EDLIB_EDOP_INSERT) {
                ti--;
            }
        }
    }

    auto pro = result.alignment[0];
    auto cur = pro;

    if (pro == EDLIB_EDOP_MATCH) {
        ti++;
        qi++;
    } else if (pro == EDLIB_EDOP_INSERT) {
        qi++;
        inStr.push_back(qData[qi]);
        if (result.alignmentLength == 1) {
            //get<1>(*vTuple).emplace_back(ti, inStr);
            variance.insMap[ti + 1] = inStr;
            inStr.clear();
        }
    } else if (pro == EDLIB_EDOP_DELETE) {
        ti++;
        deStr.push_back(tData[ti]);
        if (result.alignmentLength == 1) {
            //get<2>(*vTuple).emplace_back(ti, deStr);
            variance.delMap[ti] = deStr;
            deStr.clear();
        }
    } else if (pro == EDLIB_EDOP_MISMATCH) {
        ti++;
        qi++;
        snpStr.push_back(tData[ti]);
        snpStr.push_back('|');
        snpStr.push_back(qData[qi]);
        //get<0>(*vTuple).emplace_back(std::make_pair(ti, snpStr));
        variance.subMap[ti] = snpStr;
        snpStr.clear();
    }

    for (int i = 1; i < result.alignmentLength; i++) {
        pro = result.alignment[i - 1];
        cur = result.alignment[i];

        if (cur == EDLIB_EDOP_MATCH) {

            ti++;
            qi++;

            if (pro == EDLIB_EDOP_MATCH) {

            } else if (pro == EDLIB_EDOP_INSERT) {
                //get<1>(*vTuple).emplace_back(std::make_pair(ti, inStr));
                variance.insMap[ti - inStr.length()] = inStr;
                inStr.clear();
            } else if (pro == EDLIB_EDOP_DELETE) {
                //get<2>(*vTuple).emplace_back(std::make_pair(ti - deStr.length(), deStr));
                variance.delMap[ti - deStr.length()] = deStr;
                deStr.clear();
            } else if (pro == EDLIB_EDOP_MISMATCH) {

            }
        } else if (cur == EDLIB_EDOP_INSERT) {

            qi++;
            inStr.push_back(qData[qi]);

            if (pro == EDLIB_EDOP_DELETE) {
                //get<2>(*vTuple).emplace_back(std::make_pair(ti - deStr.length(), deStr));
                variance.delMap[ti - deStr.length()] = deStr;
                deStr.clear();
            }

            if (i == result.alignmentLength - 1) {
                //get<1>(*vTuple).emplace_back(std::make_pair(ti, inStr));
                variance.insMap[ti] = inStr;
                inStr.clear();
            }


        } else if (cur == EDLIB_EDOP_DELETE) {
            ti++;
            deStr.push_back(tData[ti]);
            if (pro == EDLIB_EDOP_INSERT) {
                //get<1>(*vTuple).emplace_back(std::make_pair(ti, inStr));
                variance.insMap[ti - inStr.length()] = inStr;
                inStr.clear();
            }

            if (i == result.alignmentLength - 1) {
                //get<2>(*vTuple).emplace_back(std::make_pair(ti - deStr.length(), deStr));
                variance.delMap[ti - deStr.length()] = deStr;
                deStr.clear();
            }


        } else if (cur == EDLIB_EDOP_MISMATCH) {

            ti++;
            qi++;

            snpStr.push_back(tData[ti]);
            snpStr.push_back('|');
            snpStr.push_back(qData[qi]);
            //get<0>(*vTuple).emplace_back(std::make_pair(ti, snpStr));
            variance.subMap[ti] = snpStr;
            snpStr.clear();

            if (pro == EDLIB_EDOP_INSERT) {
                //get<1>(*vTuple).emplace_back(std::make_pair(ti, inStr));
                variance.insMap[ti - inStr.length()] = inStr;
                inStr.clear();
            } else if (pro == EDLIB_EDOP_DELETE) {
                //get<2>(*vTuple).emplace_back(std::make_pair(ti - deStr.length(), deStr));
                variance.delMap[ti - deStr.length()] = deStr;
                deStr.clear();
            } else {

            }
        }
    }
}

//void SnpScanner::printVariance(Options * & mOptions, EdlibAlignResult & result, Variance & variance,
//        const char* & qData, const std::string & qName, const char* & tData, const std::string & tName, const int position) {
//
//    if (mOptions->mLocVars.locVarOptions.printRes) {
//        ss << "End location: " << *result.endLocations << " -> alignmentLength: " << result.alignmentLength << " with alignment distance: " << result.editDistance << "\n";
//        ss << "Ref: " << tName << " -> Query: " << qName << "\n";
//    }
//    for (int i = 0; i < result.alignmentLength; i++) {
//        if (mOptions->mLocVars.locVarOptions.printRes) ss << static_cast<unsigned> (result.alignment[i]);
//    }
//    if (mOptions->mLocVars.locVarOptions.printRes) ss << "\n";
//
//
//    for (const auto & it : variance.subMap) {
//        //std::string str = "snp: " + std::to_string(it.first) + " -> " + it.second;
//        if (mOptions->mLocVars.locVarOptions.printRes) ss << "snp: " << it.first << " -> " << it.second << "\n";
//    }
//
//    //get<0>(*vTuple).clear();
//    //get<0>(*vTuple).shrink_to_fit();
//
//    for (const auto & it : variance.insMap) {
//        //std::string str = "insertion: " + std::to_string(it.first) + " -> " + it.second;
//        if (mOptions->mLocVars.locVarOptions.printRes) ss << "insertion: " << it.first << " -> " << it.second << "\n";
//    }
//
//    //get<1>(*vTuple).clear();
//    //get<1>(*vTuple).shrink_to_fit();
//
//    for (const auto & it : variance.delMap) {
//        //std::string str = "deletion: " + std::to_string(it.first) + " -> " + it.second;
//        if (mOptions->mLocVars.locVarOptions.printRes) ss << "deletion: " << it.first << " -> " << it.second << "\n";
//    }
//
//    //get<2>(*vTuple).clear();
//    //get<2>(*vTuple).shrink_to_fit();
//
//
//    int tIdx = -1;
//    int qIdx = -1;
//    if (mOptions->mLocSnps.mLocSnpOptions.modeCode == EDLIB_MODE_HW) {
//        tIdx = *(result.endLocations);
//        for (int i = 0; i < result.alignmentLength; i++) {
//            if (result.alignment[i] != EDLIB_EDOP_INSERT) {
//                tIdx--;
//            }
//        }
//    }
//
//    for (int start = 0; start < result.alignmentLength; start += 150) {
//        // target
//        //printf("T: ");
//        if (mOptions->mLocVars.locVarOptions.printRes) ss << "T: ";
//        int startTIdx = -1;
//        for (int j = start; j < start + 150 && j < result.alignmentLength; j++) {
//            if (result.alignment[j] == EDLIB_EDOP_INSERT) {
//                //printf("-");
//                if (mOptions->mLocVars.locVarOptions.printRes) ss << "-";
//            } else {
//                //printf("%c", target[++tIdx]);
//                if (mOptions->mLocVars.locVarOptions.printRes) ss << tData[++tIdx];
//            }
//            if (j == start) {
//                startTIdx = tIdx;
//            }
//        }
//        //printf(" (%d - %d)\n", max(startTIdx, 0), tIdx);
//        if (mOptions->mLocVars.locVarOptions.printRes) {
//            ss << " (" << max(startTIdx, 0) << " - " << tIdx << ")\n";
//            ss << "   ";
//        }
//        for (int j = start; j < start + 150 && j < result.alignmentLength; j++) {
//            //printf(alignment[j] == EDLIB_EDOP_MATCH ? "|" : " ");
//            if (mOptions->mLocVars.locVarOptions.printRes) ss << (result.alignment[j] == EDLIB_EDOP_MATCH ? "|" : " ");
//        }
//        //printf("\n");
//        if (mOptions->mLocVars.locVarOptions.printRes) ss << "\n";
//
//        // query
//        //printf("Q: ");
//        if (mOptions->mLocVars.locVarOptions.printRes) ss << "Q: ";
//        int startQIdx = qIdx;
//        for (int j = start; j < start + 150 && j < result.alignmentLength; j++) {
//            if (result.alignment[j] == EDLIB_EDOP_DELETE) {
//                //printf("-");
//                if (mOptions->mLocVars.locVarOptions.printRes) ss << "-";
//            } else {
//                //printf("%c", query[++qIdx]);
//                if (mOptions->mLocVars.locVarOptions.printRes) ss << qData[++qIdx];
//            }
//            if (j == start) {
//                startQIdx = qIdx;
//            }
//        }
//        //printf(" (%d - %d)\n\n", max(startQIdx, 0), qIdx);
//        if (mOptions->mLocVars.locVarOptions.printRes) ss << " (" << max(startQIdx, 0) << " - " << qIdx << ")\n";
//    }
//}

std::tuple<int, int, bool> SnpScanner::doPrimerAlignment(const char* & qData, int qLength, const std::string & qName,
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
        if (indelSet.empty() && snpsSet.size() <= mOptions->mLocVars.locVarOptions.maxMismatchesPSeq) {
            return std::make_tuple(snpsSet.size(), endPos, true);
        } else {
            return std::make_tuple(0, 0, false);
        }
    } else {
        edlibFreeAlignResult(result);
        return std::make_tuple(0, 0, false);
    }
}

void SnpScanner::merge(Options * & mOptions, std::vector<std::map<std::string, std::map<std::string, uint32>>> & totalSnpSeqMapVec){
    //std::map<std::string, std::map<std::string, LocSnp>> allGenotypeSnpMap;
    if(totalSnpSeqMapVec.empty()){
        return;
    }

    std::map<std::string, std::map<std::string, uint32>> tmpSnpSeqsMap;
    for(const auto & it : totalSnpSeqMapVec){
        for(const auto & it2 : it){
            for(const auto & it3 : it2.second){
               tmpSnpSeqsMap[it2.first][it3.first] += (mOptions->isPaired() ?  (2 * it3.second) : it3.second);
            }
        }
    }

    std::string foutName = mOptions->prefix + "_snps_genotypes.txt";
    std::ofstream * fout = new std::ofstream();
    fout->open(foutName.c_str(), std::ofstream::out);

    if (!fout->is_open()){
        delete fout;
        fout = nullptr;
        error_exit("Can not open output file: " + foutName);
    }
    if (mOptions->verbose) loginfo("Starting to write genotype table!");
    *fout << "#Locus\tPosition\tGenotype\tNumReads\tReadsRatio\tTotalReads\tNewSnp\n";

    std::string foutName2 = mOptions->prefix + "_snps_haplotype.txt";
    std::ofstream * fout2 = new std::ofstream();
    fout2->open(foutName2.c_str(), std::ofstream::out);

    if (!fout2->is_open()) {
        delete fout2;
        fout2 = nullptr;
        error_exit("Can not open output file: " + foutName2);
    }
    if (mOptions->verbose) loginfo("Starting to write haplotype table!");
  
    *fout2 << "#Locus\tHaplotype\tNumReads\tReadsRatio\tTotalHaploReads\tTotalReads\tConclusive\tIndel\tMicroHaplotype\n";
    
    
    std::string foutName3 = mOptions->prefix + "_all_amplicon.txt";
    std::ofstream * fout3 = new std::ofstream();
    fout3->open(foutName3.c_str(), std::ofstream::out);
    *fout3 << "#Locus\tAmplicon\tNumReads\tTotalReads\n";
    
    std::string foutName4 = mOptions->prefix + "_error_rate.txt";
    std::ofstream * fout4 = new std::ofstream();
    fout4->open(foutName4.c_str(), std::ofstream::out);
    *fout4 << "#Locus\tErrorRate\tTotalEffectiveReads\n";
    
    for(const auto & it : tmpSnpSeqsMap){
        if(mOptions->mLocSnps.refLocMap.find(it.first) == mOptions->mLocSnps.refLocMap.end()){
            continue;
        }

        LocSnp2* locSnpIt = &(mOptions->mLocSnps.refLocMap[it.first]);
        std::map<int, std::map<char, int>> baseFreqMap;
        //std::set<int> posSet;
        
        for(const auto & it2 : it.second){
            locSnpIt->totReads += it2.second;
            if(it2.second > locSnpIt->maxReads){
                locSnpIt->maxReads = it2.second;
            }
        }
        
        if(locSnpIt->maxReads < mOptions->mLocSnps.mLocSnpOptions.minSeqs) continue;
        
        auto twoPeaks = getTop2MaxKeyValueVec(it.second);
        locSnpIt->totHaploReads = twoPeaks.front().second;
                  
        if(twoPeaks.size() == 1){//one homo allele
            locSnpIt->ratioHaplo = 1;
            locSnpIt->genoStr3 = "homo";
        } else if(twoPeaks.size() == 2){//two alleles
            const char* rchar1 = twoPeaks.front().first.c_str();
            const char* rchar2 = twoPeaks.back().first.c_str();
            auto mapPair = doAlignment(mOptions, "read1", rchar1, twoPeaks.front().first.length(),
                    "read2", rchar2, twoPeaks.back().first.length());

            if (mapPair.first) {
                if (mapPair.second.size() < 2) {
                    mOptions->mLocSnps.mLocSnpOptions.hmPer = mOptions->mLocSnps.mLocSnpOptions.hmPerL;
                } else {
                    mOptions->mLocSnps.mLocSnpOptions.hmPer = mOptions->mLocSnps.mLocSnpOptions.hmPerH;
                }
            } else {
                mOptions->mLocSnps.mLocSnpOptions.hmPer = mOptions->mLocSnps.mLocSnpOptions.hmPerH;
            }
            
            locSnpIt->ratioHaplo = double(twoPeaks.front().second) / (twoPeaks.front().second + twoPeaks.back().second);

            if (locSnpIt->ratioHaplo >= mOptions->mLocSnps.mLocSnpOptions.hmPer) {//homo
                twoPeaks.pop_back();
                twoPeaks.shrink_to_fit();
                locSnpIt->genoStr3 = "homo";//also include if it is heter against the ref, eg, ref: AA, target: CC;
                locSnpIt->ratioHaplo = 1;  // if the ratioHaplo for homo is not 1, it should have the seq erros
            } else if(abs(locSnpIt->ratioHaplo - 0.5) <= mOptions->mLocSnps.mLocSnpOptions.htJetter){//heter
                locSnpIt->genoStr3 = "heter";
                locSnpIt->totHaploReads += twoPeaks.back().second;
            } else  {//inconclusive;
                locSnpIt->genoStr3 = "inconclusive";
                locSnpIt->totHaploReads += twoPeaks.back().second;
            }
            
        }
        
        bool indel = false;
        std::set<int> posCorr;//for correction if there are >= 2 snps and ratio is between 0.65 - 0.90 and inconclusive one
        for (const auto & it2 : it.second) {
            bool go = false;
            if (locSnpIt->maxReads >= mOptions->mLocSnps.mLocSnpOptions.minReads4Filter) {//50
                if(it2.second >= mOptions->mLocSnps.mLocSnpOptions.minSeqs){//5
                    go = true;
                } else {
                    go = false;
                }
            } else {
                if (locSnpIt->genoStr3 == "homo") {
                    if (it2.second >= mOptions->mLocSnps.mLocSnpOptions.minSeqs) {//5
                        go = true;
                    } else {
                        go = false;
                    }
                } else {
                    if(it2.second >= mOptions->mLocSnps.mLocSnpOptions.minSeqsPer * locSnpIt->maxReads){//10%
                        go = true;
                    } else {
                        go = false;
                    }
                }
            }
            
            if(!go) continue;
            
            locSnpIt->totEffectReads += it2.second;
            *fout3 << it.first << "\t" << it2.first << "\t" << it2.second << "\t" << locSnpIt->totReads << "\n";
            
            const char* target = locSnpIt->ref.mStr.c_str();
            int targetLength = locSnpIt->ref.length();
            const char* readSeq = it2.first.c_str();
            int readLength = it2.first.length();

            SimSnps tmpSimSnps;
            tmpSimSnps.numReads = it2.second;
            tmpSimSnps.snpPosSet = locSnpIt->refSnpPosSet;
            
            auto mapPair = doAlignment(mOptions, "read", readSeq, readLength, locSnpIt->name, target, targetLength);
            
            if (mapPair.first) {//if there is no indel, that's true, including no snps
                
                //for sequence error rate;
                for(int pos = 0; pos < it2.first.length(); pos++){
                    baseFreqMap[pos][it2.first[pos]] += it2.second;
                }
                 
                if (!mapPair.second.empty()) {
                    for (auto & it3 : mapPair.second) {
                        tmpSimSnps.snpPosSet.insert(it3.first);
                        locSnpIt->snpPosSet.insert(it3.first);
                    }
                }
                
                for(const auto & it3 : tmpSimSnps.snpPosSet){
                    tmpSimSnps.snpsStr.append(std::to_string(it3 + locSnpIt->trimPos.first));
                    tmpSimSnps.snpsStr.push_back('(');
                    tmpSimSnps.snpsStr.push_back(locSnpIt->ref.mStr[it3]);
                    tmpSimSnps.snpsStr.push_back('|');
                    tmpSimSnps.snpsStr.push_back(it2.first[it3]);
                    tmpSimSnps.snpsStr.push_back(')');
                }

                if (locSnpIt->genoStr3 == "homo") {
                    if (it2.first == twoPeaks.front().first) {
                        tmpSimSnps.isHaplo = true;
                        tmpSimSnps.genoStr8 = "homo";
                        locSnpIt->snpPosSetHaplo.insert(tmpSimSnps.snpPosSet.begin(), tmpSimSnps.snpPosSet.end());
                        locSnpIt->snpPosSetTrueHaplo.insert(tmpSimSnps.snpPosSet.begin(), tmpSimSnps.snpPosSet.end());
                    } else {
                        tmpSimSnps.genoStr8 = "seqerr";
                        tmpSimSnps.isHaplo = false;
                    }

                } else if (locSnpIt->genoStr3 == "heter") {
                    if (it2.first == twoPeaks.front().first) {
                        tmpSimSnps.isHaplo = true;
                        tmpSimSnps.genoStr8 = "heter1";
                        locSnpIt->snpPosSetHaplo.insert(tmpSimSnps.snpPosSet.begin(), tmpSimSnps.snpPosSet.end());
                        locSnpIt->snpPosSetTrueHaplo.insert(tmpSimSnps.snpPosSet.begin(), tmpSimSnps.snpPosSet.end());
                    } else if (it2.first == twoPeaks.back().first) {
                        tmpSimSnps.isHaplo = true;
                        tmpSimSnps.genoStr8 = "heter2";
                        locSnpIt->snpPosSetHaplo.insert(tmpSimSnps.snpPosSet.begin(), tmpSimSnps.snpPosSet.end());
                        locSnpIt->snpPosSetTrueHaplo.insert(tmpSimSnps.snpPosSet.begin(), tmpSimSnps.snpPosSet.end());
                    } else {
                        tmpSimSnps.genoStr8 = "seqerr";
                        tmpSimSnps.isHaplo = false;
                    }
                } else {

                    if (it2.first == twoPeaks.front().first) {
                        tmpSimSnps.isHaplo = true;
                        tmpSimSnps.genoStr8 = "inHeter1";
                        locSnpIt->snpPosSetHaplo.insert(tmpSimSnps.snpPosSet.begin(), tmpSimSnps.snpPosSet.end());
                        posCorr.insert(tmpSimSnps.snpPosSet.begin(), tmpSimSnps.snpPosSet.end());
                    } else if (it2.first == twoPeaks.back().first) {
                        tmpSimSnps.isHaplo = true;
                        tmpSimSnps.genoStr8 = "inHeter2";
                        locSnpIt->snpPosSetHaplo.insert(tmpSimSnps.snpPosSet.begin(), tmpSimSnps.snpPosSet.end());
                        posCorr.insert(tmpSimSnps.snpPosSet.begin(), tmpSimSnps.snpPosSet.end());
                    } else {
                        tmpSimSnps.genoStr8 = "seqerr";
                        tmpSimSnps.isHaplo = false;
                    }
                    
                }
                locSnpIt->genoMap[it2.first] = tmpSimSnps;
            } else {//for indels
                if (locSnpIt->genoStr3 == "homo") {
                    if (it2.first == twoPeaks.front().first) {
                        tmpSimSnps.isHaplo = true;
                        tmpSimSnps.genoStr8 = "indel1";
                        indel = indel || true;
                        locSnpIt->genoMap[it2.first] = tmpSimSnps;
                    }
                } else if (locSnpIt->genoStr3 == "heter") {
                    if (it2.first == twoPeaks.front().first) {
                        tmpSimSnps.isHaplo = true;
                        tmpSimSnps.genoStr8 = "indel1";
                        indel = indel || true;
                        locSnpIt->genoMap[it2.first] = tmpSimSnps;
                    } else if (it2.first == twoPeaks.back().first) {
                        tmpSimSnps.isHaplo = true;
                        tmpSimSnps.genoStr8 = "indel2";
                        indel = indel || true;
                        locSnpIt->genoMap[it2.first] = tmpSimSnps;
                    }
                } else {
                    if (it2.first == twoPeaks.front().first) {
                        tmpSimSnps.isHaplo = true;
                        tmpSimSnps.genoStr8 = "indel1";
                        indel = indel || true;
                        locSnpIt->genoMap[it2.first] = tmpSimSnps;
                    } else if (it2.first == twoPeaks.back().first) {
                        tmpSimSnps.isHaplo = true;
                        tmpSimSnps.genoStr8 = "indel2";
                        indel = indel || true;
                        locSnpIt->genoMap[it2.first] = tmpSimSnps;
                    } 
                }
            }
        }
        
        if (locSnpIt->genoMap.empty()) continue;
        
        if (indel) {
            locSnpIt->isIndel = true;
        } else {
            if (locSnpIt->genoStr3 == "inconclusive" && (!posCorr.empty())) {
                int numSnps = 0;
                for (const auto & it2 : posCorr) {
                    if (twoPeaks.front().first[it2] != twoPeaks.back().first[it2]) {
                        numSnps++;
                    }
                }

                if (numSnps > 1) {
                    locSnpIt->snpPosSetTrueHaplo.insert(posCorr.begin(), posCorr.end());
                    locSnpIt->genoMap[twoPeaks.front().first].isHaplo = true;
                    locSnpIt->genoMap[twoPeaks.front().first].genoStr8 = "heter1";
                    locSnpIt->genoMap[twoPeaks.back().first].isHaplo = true;
                    locSnpIt->genoMap[twoPeaks.back().first].genoStr8 = "heter2";
                    locSnpIt->genoStr3 = "heter";
                }
            }
        }

        //std::set<int> totPosSet;//get total snps positions only for these true haplotypes also include the inconclusive ones
        
        locSnpIt->totPosSet.insert(locSnpIt->snpPosSet.begin(), locSnpIt->snpPosSet.end());
        
        std::string haploStr1 = "";
        std::string haploStr2 = "";

        if (!locSnpIt->isIndel) {
            for (const auto & it2 : locSnpIt->snpPosSetHaplo) {
                bool isPrint = false;
                SimSnp tmpSimSnp;
                if (twoPeaks.size() == 1) {
                    tmpSimSnp.snp1 = tmpSimSnp.snp2 = twoPeaks.front().first[it2];
                    tmpSimSnp.reads1 = tmpSimSnp.reads2 = (twoPeaks.front().second / 2);
                    tmpSimSnp.ratio = 1.0;
                    tmpSimSnp.genoStr3 = "homo";
                    if (twoPeaks.front().first[it2] == locSnpIt->ref.mStr[it2]) {
                        tmpSimSnp.color = "green";
                    } else {
                        tmpSimSnp.color = "orange";
                    }
                    isPrint = true;

                    haploStr1.push_back(twoPeaks.front().first[it2]);
                    haploStr2.push_back(twoPeaks.front().first[it2]);
                } else if (twoPeaks.size() == 2) {
                    if (twoPeaks.front().first[it2] == twoPeaks.back().first[it2]) {//AA, CC
                        tmpSimSnp.genoStr3 = "homo";
                        tmpSimSnp.reads1 = tmpSimSnp.reads2 = ((twoPeaks.front().second + twoPeaks.back().second) / 2);
                        tmpSimSnp.ratio = 1.0;
                        tmpSimSnp.snp1 = tmpSimSnp.snp2 = twoPeaks.front().first[it2];

                        if (twoPeaks.front().first[it2] == locSnpIt->ref.mStr[it2]) {//AA
                            tmpSimSnp.color = "green";
                        } else {//CC
                            tmpSimSnp.color = "orange";
                        }
                        isPrint = true;

                        haploStr1.push_back(twoPeaks.front().first[it2]);
                        haploStr2.push_back(twoPeaks.back().first[it2]);

                    } else {//AC; CA; CT; 
                        if (twoPeaks.back().first[it2] == locSnpIt->ref.mStr[it2]) {//CA;
                            tmpSimSnp.reads1 = twoPeaks.back().second;
                            tmpSimSnp.reads2 = twoPeaks.front().second;
                            tmpSimSnp.ratio = double(twoPeaks.back().second) / (twoPeaks.front().second + twoPeaks.back().second);
                            tmpSimSnp.snp1 = twoPeaks.back().first[it2];
                            tmpSimSnp.snp2 = twoPeaks.front().first[it2];

                        } else {//AC; CT;
                            tmpSimSnp.reads1 = twoPeaks.front().second;
                            tmpSimSnp.reads2 = twoPeaks.back().second;
                            tmpSimSnp.ratio = double(twoPeaks.front().second) / (twoPeaks.front().second + twoPeaks.back().second);
                            tmpSimSnp.snp1 = twoPeaks.front().first[it2];
                            tmpSimSnp.snp2 = twoPeaks.back().first[it2];
                        }

                        if (locSnpIt->genoStr3 == "heter") {
                            tmpSimSnp.color = (locSnpIt->refSnpPosSet.find(it2) == locSnpIt->refSnpPosSet.end()) ? "orange" : "red";
                            tmpSimSnp.genoStr3 = "heter";
                            isPrint = true;
                            haploStr1.push_back(twoPeaks.front().first[it2]);
                            haploStr2.push_back(twoPeaks.back().first[it2]);
                        } else if (locSnpIt->genoStr3 == "inconclusive") {
                            tmpSimSnp.color = "transparent";
                            tmpSimSnp.genoStr3 = "inconclusive";
                            isPrint = false;
                        }
                    }
                }
                locSnpIt->snpsMap[it2] = tmpSimSnp;
                if (isPrint) *fout << it.first << "\t" << (it2 + locSnpIt->trimPos.first) << "\t" << tmpSimSnp.snp1 << "|" << tmpSimSnp.snp2 << "\t" << tmpSimSnp.reads1 << "|" << tmpSimSnp.reads2 << "\t" << tmpSimSnp.ratio << "\t" << (tmpSimSnp.reads1 + tmpSimSnp.reads2) << "\t" << (locSnpIt->refSnpPosSet.find(it2) == locSnpIt->refSnpPosSet.end() ? "Y" : "N") << "\n";
            }
        }

        if (locSnpIt->isIndel) {
            if (twoPeaks.size() == 1) {
                if (locSnpIt->genoStr3 == "inconclusive") {
                    locSnpIt->haploVec.push_back(std::make_tuple(twoPeaks.front().first, haploStr1, (twoPeaks.front().second / 2), (double((twoPeaks.front().second / 2)) / locSnpIt->totHaploReads), 'N', 'Y'));
                    locSnpIt->haploVec.push_back(std::make_tuple(twoPeaks.front().first, haploStr1, (twoPeaks.front().second / 2), (double((twoPeaks.front().second / 2)) / locSnpIt->totHaploReads), 'N', 'Y'));
                } else {
                    locSnpIt->haploVec.push_back(std::make_tuple(twoPeaks.front().first, haploStr1, (twoPeaks.front().second / 2), (double((twoPeaks.front().second / 2)) / locSnpIt->totHaploReads), 'Y', 'Y'));
                    locSnpIt->haploVec.push_back(std::make_tuple(twoPeaks.front().first, haploStr1, (twoPeaks.front().second / 2), (double((twoPeaks.front().second / 2)) / locSnpIt->totHaploReads), 'Y', 'Y'));
                }
            } else if (twoPeaks.size() == 2) {
                if (locSnpIt->genoStr3 == "inconclusive") {
                    locSnpIt->haploVec.push_back(std::make_tuple(twoPeaks.front().first, haploStr1, twoPeaks.front().second, (double(twoPeaks.front().second) / locSnpIt->totHaploReads), 'N', 'Y'));
                    locSnpIt->haploVec.push_back(std::make_tuple(twoPeaks.back().first, haploStr2, twoPeaks.back().second, (double(twoPeaks.back().second) / locSnpIt->totHaploReads), 'N', 'Y'));
                } else {
                    locSnpIt->haploVec.push_back(std::make_tuple(twoPeaks.front().first, haploStr1, twoPeaks.front().second, (double(twoPeaks.front().second) / locSnpIt->totHaploReads), 'Y', 'Y'));
                    locSnpIt->haploVec.push_back(std::make_tuple(twoPeaks.back().first, haploStr2, twoPeaks.back().second, (double(twoPeaks.back().second) / locSnpIt->totHaploReads), 'Y', 'Y'));
                }
            }

        } else {
            
            if (twoPeaks.size() == 1) {
                if (locSnpIt->genoStr3 == "inconclusive") {
                    locSnpIt->haploVec.push_back(std::make_tuple(twoPeaks.front().first, haploStr1, (twoPeaks.front().second / 2), (double((twoPeaks.front().second / 2)) / locSnpIt->totHaploReads), 'N', 'N'));
                    locSnpIt->haploVec.push_back(std::make_tuple(twoPeaks.front().first, haploStr1, (twoPeaks.front().second / 2), (double((twoPeaks.front().second / 2)) / locSnpIt->totHaploReads), 'N', 'N'));
                } else {
                    locSnpIt->haploVec.push_back(std::make_tuple(twoPeaks.front().first, haploStr1, (twoPeaks.front().second / 2), (double((twoPeaks.front().second / 2)) / locSnpIt->totHaploReads), 'Y', 'N'));
                    locSnpIt->haploVec.push_back(std::make_tuple(twoPeaks.front().first, haploStr1, (twoPeaks.front().second / 2), (double((twoPeaks.front().second / 2)) / locSnpIt->totHaploReads), 'Y', 'N'));
                    
                    for(int pos = 0; pos < twoPeaks.front().first.length(); pos++){
                        int posReads = 0;
                        auto posIt = baseFreqMap[pos].find(twoPeaks.front().first[pos]);
                        if(posIt != baseFreqMap[pos].end()){
                            for(const auto & posIt2 : baseFreqMap[pos]){
                                if(posIt2.first != twoPeaks.front().first[pos]){
                                    posReads += posIt2.second;
                                }
                            }
                        }
                        locSnpIt->baseErrorMap[pos] = static_cast<double>(posReads * 100) / static_cast<double>(locSnpIt->totEffectReads);
                    }
                    
                }
            } else if (twoPeaks.size() == 2) {
                if (locSnpIt->genoStr3 == "inconclusive") {
                    locSnpIt->haploVec.push_back(std::make_tuple(twoPeaks.front().first, haploStr1, twoPeaks.front().second, (double(twoPeaks.front().second) / locSnpIt->totHaploReads), 'N', 'N'));
                    locSnpIt->haploVec.push_back(std::make_tuple(twoPeaks.back().first, haploStr2, twoPeaks.back().second, (double(twoPeaks.back().second) / locSnpIt->totHaploReads), 'N', 'N'));
                } else {
                    locSnpIt->haploVec.push_back(std::make_tuple(twoPeaks.front().first, haploStr1, twoPeaks.front().second, (double(twoPeaks.front().second) / locSnpIt->totHaploReads), 'Y', 'N'));
                    locSnpIt->haploVec.push_back(std::make_tuple(twoPeaks.back().first, haploStr2, twoPeaks.back().second, (double(twoPeaks.back().second) / locSnpIt->totHaploReads), 'Y', 'N'));
                    
                    for (int pos = 0; pos < twoPeaks.front().first.length(); pos++) {
                        int posReads = 0;

                        if (twoPeaks.front().first[pos] == twoPeaks.back().first[pos]) {
                            auto posIt = baseFreqMap[pos].find(twoPeaks.front().first[pos]);
                            if (posIt != baseFreqMap[pos].end()) {
                                for (const auto & posIt2 : baseFreqMap[pos]) {
                                    if (posIt2.first != twoPeaks.front().first[pos]) {
                                        posReads += posIt2.second;
                                    }
                                }
                            }
                        } else {
                            auto posIt = baseFreqMap[pos].find(twoPeaks.front().first[pos]);
                            auto posItb = baseFreqMap[pos].find(twoPeaks.back().first[pos]);
                            if (posIt != baseFreqMap[pos].end() && posItb != baseFreqMap[pos].end()) {
                                for (const auto & posIt2 : baseFreqMap[pos]) {
                                    if (posIt2.first != twoPeaks.front().first[pos] && posIt2.first != twoPeaks.back().first[pos]) {
                                        posReads += posIt2.second;
                                    }
                                }
                            }
                            
                        }
                        
                        locSnpIt->baseErrorMap[pos] = static_cast<double>(posReads * 100) / static_cast<double>(locSnpIt->totEffectReads);
                    }
                    
                }
            }
            
        }
        
        if(!locSnpIt->baseErrorMap.empty()){
            *fout4 << it.first << "\t";
            for(const auto & posIt : locSnpIt->baseErrorMap){
                if(posIt.first == locSnpIt->baseErrorMap.rbegin()->first){
                    *fout4 << posIt.second;
                } else {
                    *fout4 << posIt.second << ";";
                }
                
            }
            *fout4 << "\t" << locSnpIt->totEffectReads << "\n";
        }
        
        for (const auto & it2 : locSnpIt->haploVec) {
            locSnpIt->genoMap[get<0>(it2)].haploStr = get<1>(it2);
            *fout2 << it.first << "\t" << get<1>(it2) << "\t" << get<2>(it2) << "\t" << get<3>(it2) << "\t" << locSnpIt->totHaploReads << "\t" << locSnpIt->totReads << "\t" << get<4>(it2) << "\t" << get<5>(it2) << "\t" << get<0>(it2) << "\n";
        }       
    }
    
    tmpSnpSeqsMap.clear();

    fout->flush();
    fout->clear();
    fout->close();
    if (fout) {
        delete fout;
        fout = nullptr;
    }
    
    if (mOptions->verbose) loginfo("Finished writing genotype table!");

    fout2->flush();
    fout2->clear();
    fout2->close();
    if (fout2) {
        delete fout2;
        fout2 = nullptr;
    }
    
    if (mOptions->verbose) loginfo("Finished writing haplotype table!");
    
    fout3->flush();
    fout3->clear();
    fout3->close();
    if(fout3 != nullptr){
        delete fout3;
        fout3 = nullptr;
    }
    if(mOptions->verbose) loginfo("Finished writing amplicon table!");
    
    fout4->flush();
    fout4->clear();
    fout4->close();
    if(fout4 != nullptr){
        delete fout4;
        fout4 = nullptr;
    }
    if(mOptions->verbose) loginfo("Finished writing error rate table!");
    
    //return allGenotypeSnpMap;
}
