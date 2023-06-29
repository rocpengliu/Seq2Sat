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
            uint32 rpMismatches = edit_distance(it.second.rp.mStr, r1->mSeq.mStr.substr(r1->mSeq.length() - it.second.fp.length()));
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
        returnedlocus = "";
    } else {
        std::string locName = "";
        if(locMap.size() == 1){
            locName = locMap.begin()->first;
            //if(mOptions->debug) cCout("single value: " + locName, 'r');
        } else {
            std::vector<int> seqScoreVec(locMap.size());
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
        
        if(readLength != locSnpIt->ref.length()) return returnedlocus;
        
        locMap.clear();
        
        if(mOptions->debug) cCout("detected marker: " + locSnpIt->name, 'r');
        
        subSeqsMap[locSnpIt->name][r1->mSeq.mStr]++;
        returnedlocus = locName;
        
//        if (r1->mSeq.mStr.length() == locSnpIt->ref.mStr.length()) {
//
//            target = locSnpIt->ref.mStr.c_str();
//            targetLength = locSnpIt->ref.mStr.length();
//
//            std::map<std::string, std::map<std::string, LocSnp>>::iterator itGenotypeMap = subGenotypeMap.find(locSnpIt->name);//marker 
//            std::map<std::string, LocSnp> tmpGenotypeMap;
//            LocSnp* tmpGenotype = new LocSnp();
//
//            if (itGenotypeMap == subGenotypeMap.end()) {
//                if (r1->mSeq.mStr == locSnpIt->ref.mStr) {
//                    tmpGenotype->ref = locSnpIt->ref.mStr;
//                    if(mOptions->debug) cCout("identical to ref1", 'r');
//                } else {
//                    ss << "doing alignment: \n" << "read: " << r1->mSeq.mStr << "\n" << "target: " << locSnpIt->ref.mStr << "\n";
//                    tmpGenotype->snpsMap = doAlignment(readSeq, readLength, target, targetLength);
//                    tmpGenotype->ref = r1->mSeq.mStr;
//                }
//
//                tmpGenotype->name = locSnpIt->name;
//                tmpGenotype->numReads++;
//                tmpGenotypeMap[r1->mSeq.mStr] = *tmpGenotype;
//                subGenotypeMap[locSnpIt->name] = tmpGenotypeMap;
//            } else {
//                std::map<std::string, LocSnp>::iterator tmpGenotypeMapIt = itGenotypeMap->second.find(r1->mSeq.mStr);
//                if (tmpGenotypeMapIt == itGenotypeMap->second.end()) {
//                    if (r1->mSeq.mStr == locSnpIt->ref.mStr) {
//                        tmpGenotype->ref = locSnpIt->ref.mStr;
//                        if(mOptions->debug) cCout("identical to ref2", 'r');
//                    } else {
//                           ss << "doing alignment: \n" << "read: " << r1->mSeq.mStr << "\n" << "target: " << locSnpIt->ref.mStr << "\n";
//                    
//                        tmpGenotype->snpsMap = doAlignment(readSeq, readLength, target, targetLength);
//                        tmpGenotype->ref = r1->mSeq.mStr;
//                        
//                    }
//                    tmpGenotype->name = locSnpIt->name;
//                    tmpGenotype->numReads++;
//                    itGenotypeMap->second[r1->mSeq.mStr] = *tmpGenotype;
//                } else {
//                    tmpGenotypeMapIt->second.numReads++;
//                }
//            }
//            
//            if(tmpGenotype){
//                delete tmpGenotype;
//                tmpGenotype = nullptr;
//            }
//            returnedlocus = locName;
//        }
        
        if(mOptions->mEdOptions.printRes){
            cCout(ss.str(), 'g');
        }
    }
    ss.str();
    return returnedlocus;
}

std::map<int, std::pair<Sequence, Sequence>> SnpScanner::doAlignment(Options * & mOptions, std::string readName, const char* & qData, int qLength, std::string & targetName, const char* & tData, int tLength) {
    EdlibAlignResult result = edlibAlign(qData, qLength, tData, tLength,
            edlibNewAlignConfig(mOptions->mLocVars.locVarOptions.maxScorePrimer,
            EDLIB_MODE_NW,
            mOptions->mEdOptions.alignTask,
            NULL, 0));
    
    std::map<int, std::pair<Sequence, Sequence>> snpsMap;
    
    if (result.status == EDLIB_STATUS_OK) {
        if (mOptions->mEdOptions.printRes) {
            Variance variance;
            doScanVariance(mOptions, result, variance, qData, tData, *(result.endLocations));
            if (mOptions->mEdOptions.printRes) {
                //printVariance(result, variance, qData, readName, tData, targetName, *(result.endLocations));
            }

            variance.cleanVar();
        }
        for (int i = 0; i < result.alignmentLength; i++) {
            auto cur = result.alignment[i];
            if (cur == EDLIB_EDOP_MATCH) {
                
            } else if(cur == EDLIB_EDOP_MISMATCH) {
                std::string s1(1, tData[i]);
                std::string s2(1, qData[i]);
                snpsMap[i] = std::make_pair(Sequence(s1), Sequence(s2));
            } else if(cur == EDLIB_EDOP_INSERT) {
                
            } else if(cur == EDLIB_EDOP_DELETE){
                
            }
        }
        
//        for(const auto & ita : snpsMap){
//            cCout(std::to_string(ita.first) + ita.second.first.mStr + "|" + ita.second.second.mStr, 'r');
//        }
    }
    edlibFreeAlignResult(result);
    return snpsMap;
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

std::map<std::string, std::map<std::string, LocSnp>> SnpScanner::merge(Options * & mOptions, std::vector<std::map<std::string, std::map<std::string, uint32>>> & totalSnpSeqMapVec){
    
    std::map<std::string, std::map<std::string, LocSnp>> allGenotypeSnpMap;
    if(totalSnpSeqMapVec.empty()){
        return allGenotypeSnpMap;
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

    if (!fout->is_open()) error_exit("Can not open output file: " + foutName);
    if (mOptions->verbose) loginfo("Starting to write genotype table!");
    *fout << "#Locus\tPosition\tGenotype\tOriGenotype\tNumReads\tReadsRatio\tTotalReads\tNewSnp\n";

    std::string foutName2 = mOptions->prefix + "_snps_haplotype.txt";
    std::ofstream * fout2 = new std::ofstream();
    fout2->open(foutName2.c_str(), std::ofstream::out);

    if (!fout2->is_open()) error_exit("Can not open output file: " + foutName2);
    if (mOptions->verbose) loginfo("Starting to write haplotype table!");
    *fout2 << "#Locus\tHaplotype\tNumReads\tReadsRatio\tTotalReads\tSeq\n";
    
    for(const auto & it : tmpSnpSeqsMap){
        if(mOptions->mLocSnps.refLocMap.find(it.first) == mOptions->mLocSnps.refLocMap.end()){
            continue;
        }
        
        LocSnp* locSnpIt = &(mOptions->mLocSnps.refLocMap[it.first]);
        std::set<int> posSet;
        int totReads = 0, maxReads = 0;
        
        //get total and max reads
        for(const auto & it2 : it.second){
            if(it2.second >= mOptions->mLocSnps.mLocSnpOptions.minSeqs){
                totReads += it2.second;
                if(it2.second > maxReads){
                    maxReads = it2.second;
                }
            }
        }

        for(const auto & it2 : it.second){
            if(it2.second >= mOptions->mLocSnps.mLocSnpOptions.minSeqs){
                const char* target = locSnpIt->ref.mStr.c_str();
                int targetLength = locSnpIt->ref.length();
                const char* readSeq = it2.first.c_str();
                int readLength = it2.first.length();
                LocSnp tmpLocSnp;
                tmpLocSnp.ref = Sequence(it2.first);
                tmpLocSnp.numReads = it2.second;
                tmpLocSnp.readsRatio = (double) tmpLocSnp.numReads / (totReads - tmpLocSnp.numReads);
                tmpLocSnp.snpsMap = doAlignment(mOptions, "read", readSeq, readLength, locSnpIt->name, target, targetLength);
                if (!tmpLocSnp.snpsMap.empty()) {
                    for (const auto & it3 : tmpLocSnp.snpsMap) {
                        tmpLocSnp.snpPosSet.insert(it3.first);
                        posSet.insert(it3.first);
                    }
                }
                allGenotypeSnpMap[it.first][it2.first] = tmpLocSnp;
            }
        }
        
        if (!posSet.empty()) {
            locSnpIt->snpPosSet = posSet;
        }
        
        if (allGenotypeSnpMap[it.first].empty()) continue;

        std::set<int> totPosSet;
        std::set_union(locSnpIt->refSnpPosSet.begin(), locSnpIt->refSnpPosSet.end(),
                locSnpIt->snpPosSet.begin(), locSnpIt->snpPosSet.end(),
                std::inserter(totPosSet, totPosSet.begin()));

        UnitedLocSnp tmpULocSnp;
        for (const auto & it2 : totPosSet) {
            tmpULocSnp.preGenoMap[it2] = {
                {'A', 0},
                {'C', 0},
                {'G', 0},
                {'T', 0}};
            for (const auto & it3 : allGenotypeSnpMap[it.first]) {
                tmpULocSnp.preGenoMap[it2][it3.first[it2]] += it3.second.numReads;
            }
        }

        for (const auto & it2 : tmpULocSnp.preGenoMap) {
            std::vector<std::pair<char, int>> top2(2);
            std::partial_sort_copy(it2.second.begin(), it2.second.end(),
                    top2.begin(), top2.end(),
                    [](std::pair<const char, int> const & l, std::pair<const char, int> const & r) {
                        return l.second > r.second;
                    });

            SimGeno tSGeno;
            if (top2.at(1).second == 0) {
                tSGeno.geno = std::string() + top2.at(0).first + '|' + top2.at(0).first;
                tSGeno.oGeno = tSGeno.geno;
                tSGeno.read1 = top2.at(0).second;
                tSGeno.read2 = totReads - top2.at(0).second;
                tSGeno.ratio = (double) top2.at(0).second / totReads;
                if (tSGeno.ratio >= mOptions->mLocSnps.mLocSnpOptions.hmPer) {
                    tSGeno.tORf = true;
                    *fout << it.first << "\t" << it2.first << "\t" << tSGeno.geno << "\t" << tSGeno.oGeno << "\t" << top2.at(0).second << "|" << totReads << "\t" << tSGeno.ratio << "\t" << totReads << "\t" << (locSnpIt->refSnpPosSet.find(it2.first) == locSnpIt->refSnpPosSet.end() ? "Y" : "N") << "\n";
                } else {
                    tSGeno.tORf = false;
                }
            } else {
                std::string bas = locSnpIt->ref.mStr;
                if( bas[it2.first] == top2.at(0).first) {
                    tSGeno.oGeno = std::string() + top2.at(0).first + '|' + top2.at(1).first;
                    tSGeno.read1 = top2.at(0).second;
                    tSGeno.read2 = top2.at(1).second;
                } else {
                    tSGeno.oGeno = std::string() + top2.at(1).first + '|' + top2.at(0).first;
                    tSGeno.read1 = top2.at(1).second;
                    tSGeno.read2 = top2.at(0).second;
                }

                tSGeno.ratio = (double) std::max(top2.at(0).second, top2.at(1).second) / (top2.at(0).second + top2.at(1).second);
                
                if (tSGeno.ratio >= (0.5 - mOptions->mLocSnps.mLocSnpOptions.htJetter) && tSGeno.ratio <= (0.5 + mOptions->mLocSnps.mLocSnpOptions.htJetter)) {
                    tSGeno.tORf = true;
                    tmpULocSnp.heter = true;
                    tSGeno.geno = tSGeno.oGeno;
                    *fout << it.first << "\t" << it2.first << "\t" << tSGeno.geno << "\t" << top2.at(0).second << "|" << top2.at(1).second << "\t" << tSGeno.ratio << "\t" << totReads << "\t" << (locSnpIt->refSnpPosSet.find(it2.first) == locSnpIt->refSnpPosSet.end() ? "Y" : "N") << "\n";
                } else if(tSGeno.ratio >= mOptions->mLocSnps.mLocSnpOptions.hmPer){

                    if (bas[it2.first] == top2.at(0).first) {
                        tSGeno.geno = std::string() + top2.at(0).first + '|' + top2.at(0).first;
                    } else {
                        tSGeno.geno = std::string() + top2.at(1).first + '|' + top2.at(1).first;
                    }
                    
                    //tmpULocSnp.heter = false;
                    if(locSnpIt->refSnpPosSet.find(it2.first) == locSnpIt->refSnpPosSet.end()){
                        tSGeno.tORf = false;
                    } else {
                        tSGeno.tORf = true;
                        *fout << it.first << "\t" << it2.first << "\t" << tSGeno.geno << "\t" << tSGeno.oGeno << "\t" << top2.at(0).second << "|" << top2.at(1).second << "\t" << tSGeno.ratio << "\t" << totReads << "\t" << (locSnpIt->refSnpPosSet.find(it2.first) == locSnpIt->refSnpPosSet.end() ? "Y" : "N") << "\n";
                    }
                } else {
                    tSGeno.geno = tSGeno.oGeno;//for the ambiguous genotype but need to show in html report
                    tSGeno.tORf = false;
                }
                
            }
            tmpULocSnp.snpGenoMap[it2.first] = tSGeno;
        }

        locSnpIt->uGeno = tmpULocSnp;
        locSnpIt->totReads = totReads;

        std::string haploSeq1 = "";
        std::string haploSeq2 = "";
        std::string haploStr1 = "";
        std::string haploStr2 = "";
        int haploReads1 = 0;
        int haploReads2 = 0;
        if (tmpULocSnp.heter) {
            for (const auto & it3 : allGenotypeSnpMap[it.first]) {
                if (it3.second.numReads == maxReads) {
                    haploReads1 = maxReads;
                } else {
                    if (it3.second.numReads > haploReads2) {
                        haploReads2 = it3.second.numReads;
                    }
                }
            }

            for (auto & it3 : allGenotypeSnpMap[it.first]) {
                if (it3.second.numReads == haploReads1) {
                    it3.second.isHaplotype = true;
                    haploSeq1 = it3.first;
                    for (auto & it4 : locSnpIt->uGeno.snpGenoMap) {
                        if (it4.second.tORf) {
                            haploStr1.push_back(haploSeq1[it4.first]);
                        }
                    }

                } else if (it3.second.numReads == haploReads2) {
                    it3.second.isHaplotype = true;
                    haploSeq2 = it3.first;
                    for (auto & it4 : locSnpIt->uGeno.snpGenoMap) {
                        if (it4.second.tORf) {
                            haploStr2.push_back(haploSeq2[it4.first]);
                        }
                    }
                }
            }

        } else {
            for (auto & it3 : allGenotypeSnpMap[it.first]) {
                if (it3.second.numReads == maxReads) {
                    it3.second.isHaplotype = true;
                    haploSeq1 = haploSeq2 = it3.first;
                    for (auto & it4 : locSnpIt->uGeno.snpGenoMap) {
                        if (it4.second.tORf) {
                            haploStr1 += haploSeq1[it4.first];
                            haploStr2 += haploSeq2[it4.first];
                            haploReads1 = haploReads2 = (maxReads / 2);
                        }
                    }
                }
            }
        }

        locSnpIt->haploVec.push_back(std::make_tuple(haploSeq1, haploStr1, haploReads1));
        locSnpIt->haploVec.push_back(std::make_tuple(haploSeq2, haploStr2, haploReads2));

        *fout2 << it.first << "\t" << haploStr1 << "\t" << haploReads1 << "\t" << ((double) haploReads1 / totReads) << "\t" << totReads << "\t" << haploSeq1 << "\n";
        *fout2 << it.first << "\t" << haploStr2 << "\t" << haploReads2 << "\t" << ((double) haploReads2 / totReads) << "\t" << totReads << "\t" << haploSeq2 << "\n";
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
    
    return allGenotypeSnpMap;
}