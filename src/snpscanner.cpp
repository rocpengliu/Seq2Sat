#include <valarray>

#include "snpscanner.h"

SnpScanner::SnpScanner(Options* opt) {
    mOptions = opt;
    subGenotypeMap.clear();
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

bool SnpScanner::scanVar(Read* & r1) {
    readSeq = r1->mSeq.mStr.c_str();
    readLength = r1->mSeq.mStr.length();
    readName = r1->mName;
    std::map<std::string, std::pair<int, int>> locMap;
    //cCout("cccccccccccccccccccccccccccccc000000000000000000000");
    for (auto & it : mOptions->mLocSnps.refLocMap) {
       // cCout("cccccccccccccccccccccccccccccc111111111111111");
        it.second.print();
        r1->print();
        fpData = it.second.fp.mStr.c_str();
        fpLength = it.second.fp.mStr.length();
        rpData = it.second.rp.mStr.c_str();
        rpLength = it.second.rp.mStr.length();
        bool fpMatched = true;
        bool rpMatched = true;
        int fpMismatches = 0;
        int rpMismatches = 0;

        //cCout("cccccccccccccccccccccccccccccc222222222222222222222222");
        for (int i = 0; i < fpLength; i++) {
            if (fpData[i] != readSeq[i]) {
                fpMismatches++;
            }
            if (fpMismatches >= mOptions->mLocSnps.mLocSnpOptions.maxMismatchesPSeq) {
                fpMatched = false;
                break;
            }
        }

        // cCout("fpMismaches: " + std::to_string(fpMismatches), 'g');

        if (fpMatched) {
            int rlen = r1->length() - it.second.rp.mStr.length();
            for (int i = 0; i < rpLength; i++) {
                if (rpData[i] != readSeq[rlen + i]) {
                    rpMismatches++;
                }
                if (rpMismatches >= mOptions->mLocSnps.mLocSnpOptions.maxMismatchesPSeq) {
                    rpMatched = false;
                    break;
                }
            }

            if (rpMatched) {
                locMap[it.first] = std::make_pair((fpMismatches + rpMismatches), readLength - it.second.fp.mStr.length() - it.second.rp.mStr.length());
            }
        }

    }
     
    if(!locMap.empty()){
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
            //warning, what if there are multiple identity values
            seqScoreVec.clear();

            //if(mOptions->debug) cCout("minValue: " + std::to_string(minValue), 'r');
            for (const auto & it : locMap) {
                if (it.second.first == minValue) {
                    locName = it.first;
                    break; //warning, what if there are multiple identity values
                }
            }
        }
        locSnpIt = mOptions->mLocSnps.refLocMap[locName];
        locSnpIt.print();
//        cCout(locName, 'y');
//        cCout(locMap[locName].first, 'y');
//        cCout(locMap[locName].second, 'y');
//        r1->print();
//        r1->trimFront(locSnpIt.fp.length());
//        r1->resize(locMap[locName].second);
        //r1->resize(locSnpIt.ref.mStr.length());
        locMap.clear();
        
        if(mOptions->debug) cCout("detected marker: " + locSnpIt.name, 'r');
        
//        r1->print();
//        
//        cCout("0000000000000000000", 'r');
//        cCout(r1->mSeq.mStr.length(), 'g');
//        cCout(locSnpIt.ref.mStr.length(), 'r');
        if (r1->mSeq.mStr.length() == locSnpIt.ref.mStr.length()) {

            readSeq = r1->mSeq.mStr.c_str();
            readLength = r1->mSeq.mStr.length();
            target = locSnpIt.ref.mStr.c_str();
            targetLength = locSnpIt.ref.mStr.length();

            std::map<std::string, std::map < std::string, LocSnp>>::iterator itGenotypeMap = subGenotypeMap.find(locSnpIt.name);
            std::map<std::string, LocSnp> tmpGenotypeMap;
            LocSnp tmpGenotype;

            //cCout("aaaaaaaaaaaaaaaaaaaaaaaaaa");
            //bool goAlignment = false;
            if (itGenotypeMap == subGenotypeMap.end()) {
                //cCout("aaaaaaaaaaaaaaaaaaaaaaaaaa111111111111111111");
                if (r1->mSeq.mStr == locSnpIt.ref.mStr) {
                    //tmpGenotype = locSnpIt;
                    tmpGenotype.ref = locSnpIt.ref.mStr;
                    //cCout("aaaaaaaaaaaaaaaaaaaaaaaaaa2222222222222222222222");
                    if(mOptions->debug) cCout("identical to ref1", 'r');
                } else {
                    //goAlignment = true;
                   // cCout("aaaaaaaaaaaaaaaaaaaaaaaaaa33333333333333333333");
//                    target = locSnpIt.ref.mStr.c_str();
//                    targetLength = locSnpIt.ref.mStr.length();
//                    readSeq = r1->mSeq.mStr.c_str();
//                    readLength = r1->mSeq.mStr.length();
                    ss << "doing alignment: \n" << "read: " << r1->mSeq.mStr << "\n" << "target: " << locSnpIt.ref.mStr << "\n";
                    
                    tmpGenotype.snpsMap = doAlignment(readSeq, readLength, target, targetLength);
                    tmpGenotype.ref = r1->mSeq.mStr;
                    
                }
                //cCout("aaaaaaaaaaaaaaaaaaaaaaaaaa44444444444444444444");
                tmpGenotype.name = locSnpIt.name;
                tmpGenotype.numReads++;
                tmpGenotypeMap[r1->mSeq.mStr] = tmpGenotype;
                subGenotypeMap[locSnpIt.name] = tmpGenotypeMap;
            } else {
                //cCout("aaaaaaaaaaaaaaaaaaaaaaaaaa55555555555555555555555555");
                std::map<std::string, LocSnp>::iterator tmpGenotypeMapIt = itGenotypeMap->second.find(r1->mSeq.mStr);
                if (tmpGenotypeMapIt == itGenotypeMap->second.end()) {
                   // cCout("aaaaaaaaaaaaaaaaaaaaaaaaaa666666666666666666666666666");
                    if (r1->mSeq.mStr == locSnpIt.ref.mStr) {
                        //cCout("aaaaaaaaaaaaaaaaaaaaaaaaaa777777777777777777777777777777");
                        //tmpGenotype = locSnpIt;
                        tmpGenotype.ref = locSnpIt.ref.mStr;
                        if(mOptions->debug) cCout("identical to ref2", 'r');
                    } else {
                       // cCout("aaaaaaaaaaaaaaaaaaaaaaaaaa88888888888888888888888888888");
                        //goAlignment = true;
//                        target = locSnpIt.ref.mStr.c_str();
//                        targetLength = locSnpIt.ref.mStr.length();
//                        readSeq = r1->mSeq.mStr.c_str();
//                        readLength = r1->mSeq.mStr.length();
                           ss << "doing alignment: \n" << "read: " << r1->mSeq.mStr << "\n" << "target: " << locSnpIt.ref.mStr << "\n";
                    
                        tmpGenotype.snpsMap = doAlignment(readSeq, readLength, target, targetLength);
                        tmpGenotype.ref = r1->mSeq.mStr;
                        
                    }
                    //cCout("aaaaaaaaaaaaaaaaaaaaaaaaaa9999999999999999999999999");
                    tmpGenotype.name = locSnpIt.name;
                    tmpGenotype.numReads++;
                    itGenotypeMap->second[r1->mSeq.mStr] = tmpGenotype;
                } else {
                    //cCout("aaaaaaaaaaaaaaaaaaaaaaaaaa000000000000000000");
                    tmpGenotypeMapIt->second.numReads++;
                }
            }
        }
        
        if(mOptions->mLocVars.locVarOptions.printRes){
            cCout(ss.str(), 'g');
        }
    }
    ss.str();
    return true;
}

std::map<int, std::pair<Sequence, Sequence>> SnpScanner::doAlignment(const char* & qData, int qLength, const char* & tData, int tLength) {
    EdlibAlignResult result = edlibAlign(qData, qLength, tData, tLength,
            edlibNewAlignConfig(mOptions->mLocVars.locVarOptions.maxScorePrimer,
            EDLIB_MODE_NW,
            mOptions->mLocSnps.mLocSnpOptions.alignTask,
            NULL, 0));
    
    std::map<int, std::pair<Sequence, Sequence>> snpsMap;
    
    if (result.status == EDLIB_STATUS_OK) {
        if (mOptions->mLocVars.locVarOptions.printRes) {
            Variance variance;
            doScanVariance(result, variance, qData, tData, *(result.endLocations));
            if (mOptions->mLocVars.locVarOptions.printRes) {
                printVariance(result, variance, qData, readName, tData, locSnpIt.name, *(result.endLocations));
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

void SnpScanner::doScanVariance(EdlibAlignResult & result, Variance & variance,
        const char* & qData, const char* & tData, const int position) {

    int ti = -1, qi = -1;
    std::string inStr, deStr, snpStr;

    if (result.alignmentLength < 2) {
        if (mOptions->mLocVars.locVarOptions.printRes) ss << "You must have at least 2 bp!\n";
        return;
    }

    if (mOptions->mLocSnps.mLocSnpOptions.modeCode == EDLIB_MODE_HW) {
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

void SnpScanner::printVariance(EdlibAlignResult & result, Variance & variance,
        const char* & qData, const std::string & qName, const char* & tData, const std::string & tName, const int position) {

    if (mOptions->mLocVars.locVarOptions.printRes) {
        ss << "End location: " << *result.endLocations << " -> alignmentLength: " << result.alignmentLength << " with alignment distance: " << result.editDistance << "\n";
        ss << "Ref: " << tName << " -> Query: " << qName << "\n";
    }
    for (int i = 0; i < result.alignmentLength; i++) {
        if (mOptions->mLocVars.locVarOptions.printRes) ss << static_cast<unsigned> (result.alignment[i]);
    }
    if (mOptions->mLocVars.locVarOptions.printRes) ss << "\n";


    for (const auto & it : variance.subMap) {
        //std::string str = "snp: " + std::to_string(it.first) + " -> " + it.second;
        if (mOptions->mLocVars.locVarOptions.printRes) ss << "snp: " << it.first << " -> " << it.second << "\n";
    }

    //get<0>(*vTuple).clear();
    //get<0>(*vTuple).shrink_to_fit();

    for (const auto & it : variance.insMap) {
        //std::string str = "insertion: " + std::to_string(it.first) + " -> " + it.second;
        if (mOptions->mLocVars.locVarOptions.printRes) ss << "insertion: " << it.first << " -> " << it.second << "\n";
    }

    //get<1>(*vTuple).clear();
    //get<1>(*vTuple).shrink_to_fit();

    for (const auto & it : variance.delMap) {
        //std::string str = "deletion: " + std::to_string(it.first) + " -> " + it.second;
        if (mOptions->mLocVars.locVarOptions.printRes) ss << "deletion: " << it.first << " -> " << it.second << "\n";
    }

    //get<2>(*vTuple).clear();
    //get<2>(*vTuple).shrink_to_fit();


    int tIdx = -1;
    int qIdx = -1;
    if (mOptions->mLocSnps.mLocSnpOptions.modeCode == EDLIB_MODE_HW) {
        tIdx = *(result.endLocations);
        for (int i = 0; i < result.alignmentLength; i++) {
            if (result.alignment[i] != EDLIB_EDOP_INSERT) {
                tIdx--;
            }
        }
    }

    for (int start = 0; start < result.alignmentLength; start += 150) {
        // target
        //printf("T: ");
        if (mOptions->mLocVars.locVarOptions.printRes) ss << "T: ";
        int startTIdx = -1;
        for (int j = start; j < start + 150 && j < result.alignmentLength; j++) {
            if (result.alignment[j] == EDLIB_EDOP_INSERT) {
                //printf("-");
                if (mOptions->mLocVars.locVarOptions.printRes) ss << "-";
            } else {
                //printf("%c", target[++tIdx]);
                if (mOptions->mLocVars.locVarOptions.printRes) ss << tData[++tIdx];
            }
            if (j == start) {
                startTIdx = tIdx;
            }
        }
        //printf(" (%d - %d)\n", max(startTIdx, 0), tIdx);
        if (mOptions->mLocVars.locVarOptions.printRes) {
            ss << " (" << max(startTIdx, 0) << " - " << tIdx << ")\n";
            ss << "   ";
        }
        for (int j = start; j < start + 150 && j < result.alignmentLength; j++) {
            //printf(alignment[j] == EDLIB_EDOP_MATCH ? "|" : " ");
            if (mOptions->mLocVars.locVarOptions.printRes) ss << (result.alignment[j] == EDLIB_EDOP_MATCH ? "|" : " ");
        }
        //printf("\n");
        if (mOptions->mLocVars.locVarOptions.printRes) ss << "\n";

        // query
        //printf("Q: ");
        if (mOptions->mLocVars.locVarOptions.printRes) ss << "Q: ";
        int startQIdx = qIdx;
        for (int j = start; j < start + 150 && j < result.alignmentLength; j++) {
            if (result.alignment[j] == EDLIB_EDOP_DELETE) {
                //printf("-");
                if (mOptions->mLocVars.locVarOptions.printRes) ss << "-";
            } else {
                //printf("%c", query[++qIdx]);
                if (mOptions->mLocVars.locVarOptions.printRes) ss << qData[++qIdx];
            }
            if (j == start) {
                startQIdx = qIdx;
            }
        }
        //printf(" (%d - %d)\n\n", max(startQIdx, 0), qIdx);
        if (mOptions->mLocVars.locVarOptions.printRes) ss << " (" << max(startQIdx, 0) << " - " << qIdx << ")\n";
    }
}

std::map<std::string, std::map<std::string, LocSnp>> SnpScanner::merge(Options * & mOptions, std::vector<std::map<std::string, std::map < std::string, LocSnp>>> & totalGenotypeSnpMapVec){
    
    cCout("starting merging", 'r');
    std::map<std::string, std::map<std::string, LocSnp>> tmpGenotypeSnpMap;
    std::map<std::string, std::map<std::string, LocSnp>> allGenotypeSnpMap;
    if(totalGenotypeSnpMapVec.empty()){
        return allGenotypeSnpMap;
    }

    //cCout("aaaaaaaaaaaaaaaaa", 'r');
    std::set<std::string> locSet;
    for(const auto & it : totalGenotypeSnpMapVec){
        for(const auto & it2 : it){
            locSet.insert(it2.first);
        }
    }

    //cCout("aaaaaaaaaaaaaaaaa11111111111", 'r');
    
    for(const auto it : locSet){
        std::map<std::string, LocSnp> eachSnpMap;
        for(const auto & it2 : totalGenotypeSnpMapVec){
            auto it3 = it2.find(it);
            if(it3 != it2.end()){
                for(const auto & it4 : it3->second){
                    auto it5 = eachSnpMap.find(it4.first);
                    if(it5 == eachSnpMap.end()){
                        eachSnpMap[it4.first] = it4.second;
                    } else {
                        it5->second.numReads += it4.second.numReads;
                    }
                }
            }
        }
        tmpGenotypeSnpMap[it] = eachSnpMap;
    }
    
    //cCout("aaaaaaaaaaaaaaaaa2222222222222", 'r');
    std::string foutName = mOptions->prefix + "_snps_genotypes.txt";
    std::ofstream * fout = new std::ofstream();
    fout->open(foutName.c_str(), std::ofstream::out);

    if (!fout->is_open()) error_exit("Can not open output file: " + foutName);
    if (mOptions->verbose) loginfo("Starting to write genotype table!");

    *fout << "#Locus\tSnps\tNumReads\tSequence\n";
    
    for( auto & it : tmpGenotypeSnpMap){
        std::map<std::string, LocSnp> tmpSnpMap;
        for( auto & it2 : it.second){
            if(it2.second.numReads >= mOptions->mLocSnps.mLocSnpOptions.minSeqs){
                tmpSnpMap[it2.first] = it2.second;
                *fout << it.first << "\t" << it2.second.getGenotype() << "\t" << it2.second.numReads << "\t" << it2.second.ref.mStr << "\n";
            }
        }
        
        std::set<int> posSet;
        for(const auto & it2 : tmpSnpMap){
            for(const auto & it3 : it2.second.snpsMap){
                posSet.insert(it3.first);
            }
        }
        allGenotypeSnpMap[it.first] = tmpSnpMap;
        auto it2 = mOptions->mLocSnps.refLocMap.find(it.first);
        if(it2 != mOptions->mLocSnps.refLocMap.end()){
            it2->second.snpPosSet = posSet;
        }
    }

    fout->flush();
    fout->clear();
    fout->close();
    if (fout) {
        delete fout;
        fout = NULL;
    }
    if (mOptions->verbose) loginfo("Finished writing genotype table!");
    //cCout("aaaaaaaaaaaaaaaaa33333333333", 'r');
    tmpGenotypeSnpMap.clear();
    return allGenotypeSnpMap;
}