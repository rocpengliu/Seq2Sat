#include <valarray>

#include "ssrscanner.h"
#include "options.h"

SsrScanner::SsrScanner(Options* opt) {
    mOptions = opt;
    target = NULL;
    targetLength = 0;
    readSeq = NULL;
    readLength = 0;
    readName = "";
    fpData = NULL;
    fpLength = 0;
    rpData = NULL;
    rpLength = 0;
    ss.str("");
    ffEndPos = 0;
    rfStartPos = 0;
    readWithoutFF.clear();
    readWithoutFFRF.clear();
    locMra.clear();
    tmpAllGenotypeMap.clear();
    //sortedAllGenotypeMap.clear();
    enhancer = "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";
    returnedlocus.clear();
    locVarIt = nullptr;
    checkLoci = true;
    //tmpSex = &(mOptions->mSex);
    //tmpSex = mOptions->mSex;
    tmpSexMap.clear();
}

SsrScanner::~SsrScanner() {

}

bool SsrScanner::scanVar(Read* & r1, Read* & r2) {
    return true;
}

void SsrScanner::resetData() {
    readSeq = NULL;
    readLength = 0;
    target = NULL;
    targetLength = 0;
}

std::string SsrScanner::trim2ends(Read* & tmpRead, int startPos, int endPos) {
    std::string rStr = tmpRead->mSeq.mStr;
    return rStr.substr(startPos + 1, endPos - startPos - 1);

}

int SsrScanner::doSubStrAlignment(const char* & pData, int pLength, const char* & rData, int rLength,
        int & numMismaches, int & numDeletions, int & numInsertions, bool rev) {

    int trimPos = 0;

    EdlibAlignResult result = edlibAlign(pData, pLength, rData, rLength,
            edlibNewAlignConfig(mOptions->mLocVars.locVarOptions.maxScorePrimer,
            EDLIB_MODE_SHW,
            EDLIB_TASK_PATH,
            NULL, 0));

    if (result.status == EDLIB_STATUS_OK) {

        if (rev) {

            if (*result.endLocations - *result.startLocations >= pLength - mOptions->mLocVars.locVarOptions.maxDeletionRPrimer) {

                for (int i = 0; i < result.alignmentLength; i++) {
                    auto j = static_cast<unsigned> (result.alignment[i]);
                    if (j == 3) {
                        numMismaches++;
                    } else if (j == 2) {
                        numDeletions++;
                    } else if (j == 1) {
                        numInsertions++;
                    }
                }

                if (numMismaches <= mOptions->mLocVars.locVarOptions.maxMismatchesRPrimer &&
                        numDeletions <= mOptions->mLocVars.locVarOptions.maxDeletionRPrimer &&
                        numInsertions <= mOptions->mLocVars.locVarOptions.maxInsertionRPrimer) {

                    trimPos = rLength - (*result.endLocations + 1);
                    if (mOptions->mLocVars.locVarOptions.printRes) {
                        ss << "rev_RPrimer: " << pData << " (" << pLength << " residues) against rev_read " << readName << ": (" << rLength << " residues)" << ": score = " << result.editDistance << "\n";
                        //cCout();
                        printAlignment(pData, rData, result.alignment, result.alignmentLength, *result.endLocations, mOptions->mLocVars.locVarOptions.modeCode);
                    }
                }

            }

        } else {
            if (*result.endLocations - *result.startLocations <= pLength + mOptions->mLocVars.locVarOptions.maxDeletionRPrimer &&
                    rLength - mOptions->mLocVars.locVarOptions.maxDeletionRPrimer <= *result.endLocations &&
                    *result.endLocations < rLength) {

                for (int i = 0; i < result.alignmentLength; i++) {
                    auto j = static_cast<unsigned> (result.alignment[i]);
                    if (j == 3) {
                        numMismaches++;
                    } else if (j == 2) {
                        numDeletions++;
                    } else if (j == 1) {
                        numInsertions++;
                    }
                }

                if (numMismaches <= mOptions->mLocVars.locVarOptions.maxMismatchesRPrimer &&
                        numDeletions <= mOptions->mLocVars.locVarOptions.maxDeletionRPrimer &&
                        numInsertions <= mOptions->mLocVars.locVarOptions.maxInsertionRPrimer) {
                    trimPos = *result.startLocations;
                }

                if (mOptions->mLocVars.locVarOptions.printRes) {
                    ss << "RPrimer: " << pData << " (" << pLength << " residues) against read " << readName << ": (" << rLength << " residues)" << ": score = " << result.editDistance << "\n";
                    //cCout();
                    printAlignment(pData, rData, result.alignment, result.alignmentLength, *result.endLocations, mOptions->mLocVars.locVarOptions.modeCode);
                }
            }
        }
    }
    edlibFreeAlignResult(result);
    return trimPos;
}

std::map<int, std::string> SsrScanner::doSimpleAlignment(Options * & mOptions, const char* & qData, int qLength, const char* & tData, int tLength) {
    EdlibAlignResult result = edlibAlign(qData, qLength, tData, tLength,
            edlibNewAlignConfig(mOptions->mLocVars.locVarOptions.maxScorePrimer,
            mOptions->mLocVars.locVarOptions.modeCode,
            mOptions->mLocVars.locVarOptions.alignTask,
            NULL, 0));

    std::map<int, std::string> snpsMap;

    if (result.status == EDLIB_STATUS_OK) {
        for (int i = 0; i < result.alignmentLength; i++) {
            auto cur = result.alignment[i];
            if (cur == EDLIB_EDOP_MATCH) {

            } else if (cur == EDLIB_EDOP_MISMATCH) {
                snpsMap[i] = tData[i];
            } else if (cur == EDLIB_EDOP_INSERT) {

            } else if (cur == EDLIB_EDOP_DELETE) {

            }
        }

    }
    edlibFreeAlignResult(result);
    return snpsMap;
}

void SsrScanner::printAlignment(const char* & query, const char* & target,
        const unsigned char* alignment, const int alignmentLength,
        const int position, const EdlibAlignMode modeCode) {
    int tIdx = -1;
    int qIdx = -1;
    if (modeCode == EDLIB_MODE_HW) {
        tIdx = position;
        for (int i = 0; i < alignmentLength; i++) {
            if (alignment[i] != EDLIB_EDOP_INSERT)
                tIdx--;
        }
    }
    for (int start = 0; start < alignmentLength; start += 50) {
        // target
        ss << "T: ";
        int startTIdx = -1;
        for (int j = start; j < start + 50 && j < alignmentLength; j++) {
            if (alignment[j] == EDLIB_EDOP_INSERT)
                ss << "-";
            else
                ss << target[++tIdx];
            if (j == start)
                startTIdx = tIdx;
        }
        ss << " (" << max(startTIdx, 0) << " - " << tIdx << ")\n";

        // match / mismatch
        ss << "   ";
        for (int j = start; j < start + 50 && j < alignmentLength; j++) {
            ss << (alignment[j] == EDLIB_EDOP_MATCH ? "|" : " ");
        }
        ss << "\n";

        // query
        ss << "Q: ";
        int startQIdx = qIdx;
        for (int j = start; j < start + 50 && j < alignmentLength; j++) {
            if (alignment[j] == EDLIB_EDOP_DELETE)
                ss << "-";
            else
                ss << query[++qIdx];
            if (j == start)
                startQIdx = qIdx;
        }
        ss << " (" << max(startQIdx, 0) << " - " << qIdx << ")\n\n";
    }
}

int SsrScanner::doAlignment(const char* & qData, int qLength, const std::string & qName,
        const char* & tData, int tLength, const std::string & tName, Variance & variance, bool printAlignment) {
    int endPos = 0;
    EdlibAlignResult result = edlibAlign(qData, qLength, tData, tLength,
            edlibNewAlignConfig(mOptions->mLocVars.locVarOptions.maxScorePrimer,
            mOptions->mLocVars.locVarOptions.modeCode,
            mOptions->mLocVars.locVarOptions.alignTask,
            NULL, 0));

    if (result.status == EDLIB_STATUS_OK) {
        endPos = *result.endLocations;
        doScanVariance(result, variance, qData, tData, *(result.endLocations));
        if (mOptions->mLocVars.locVarOptions.printRes && printAlignment) {
            printVariance(result, variance, qData, qName, tData, tName, *(result.endLocations));
        }
    }
    edlibFreeAlignResult(result);
    return endPos;
}

void SsrScanner::doScanVariance(EdlibAlignResult & result, Variance & variance,
        const char* & qData, const char* & tData, const int position) {

    int ti = -1, qi = -1;
    std::string inStr, deStr, snpStr;

    if (result.alignmentLength < 2) {
        if (mOptions->mLocVars.locVarOptions.printRes) ss << "You must have at least 2 bp!\n";
        return;
    }

    if (mOptions->mLocVars.locVarOptions.modeCode == EDLIB_MODE_HW) {
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

void SsrScanner::printVariance(EdlibAlignResult & result, Variance & variance,
        const char* & qData, const std::string & qName, const char* & tData, const std::string & tName, const int position) {

    if (mOptions->mLocVars.locVarOptions.printRes) {
        ss << "End location: " << *result.endLocations << " -> alignmentLength: " << result.alignmentLength << " with alignment distance: " << result.editDistance << "\n";
        ss << "Ref: " << tName << " -> Query: " << qName << "\n";
        ss << "Repeat seq: " << locVarIt->mra.mStr << " | " << locVarIt->mra.reverseComplement().mStr << "\n";
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
    if (mOptions->mLocVars.locVarOptions.modeCode == EDLIB_MODE_HW) {
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

std::pair<size_t, int> SsrScanner::analyzeMRA(std::string rawStr, //rawStr must be passed by value not reference
        const std::string & ssr, std::size_t & ffTrimPos, std::size_t & rfTrimPos) {

    if (ffTrimPos + mOptions->mLocVars.locVarOptions.minNSSRUnit * (locVarIt->repuit.mStr.length() + locVarIt->repuit2.mStr.length()) > rfTrimPos) {
        ffTrimPos = 0;
        rfTrimPos = rawStr.length() - 1;
    }

    if (rfTrimPos != 0) {
        rawStr = rawStr.substr(0, rfTrimPos);
    }
    rawStr = rawStr.substr(ffTrimPos);

    std::pair<size_t, int> posSSRlenPair(0, 0);
    std::map<size_t, int> posMap;
    size_t startPos = rawStr.find(ssr);
    size_t endPos = 0;
    if (startPos == std::string::npos) {
        return posSSRlenPair;
    }
    while (startPos != std::string::npos) {
        endPos = rawStr.find(ssr, startPos + ssr.size());
        if (endPos == std::string::npos) {
            if (posMap.empty()) {
                posMap[startPos] = 1;
            } else {
                auto it = posMap.rbegin();
                if (it->first + it->second * ssr.length() == startPos) {
                    it->second++;
                } else {
                    posMap[startPos] = 1;
                }
            }
            break;
        } else {
            if (endPos - startPos == ssr.size()) {
                if (posMap.empty()) {
                    posMap[startPos]++;
                } else {
                    posMap.rbegin()->second++;
                }

            } else {
                if (posMap.empty()) {
                    posMap[startPos]++;
                } else {
                    posMap.rbegin()->second++;
                }
                if (endPos != std::string::npos) {
                    posMap[endPos] = 0;
                }
            }
        }
        startPos = endPos;
    }

    size_t coreKey = 0;
    int coreCount = 0;
    std::vector<std::pair<size_t, int>> preVec;
    for (const auto & it : posMap) {
        //cCout(it.first, 'b');
        //cCout(it.second, 'm');
        preVec.emplace_back(std::make_pair(it.first, it.second));
        if (it.second > coreCount) {
            coreCount = it.second;
            coreKey = it.first;
        }
    }

    std::map<size_t, int> posMapFil = {
        {coreKey, coreCount}
    };
    int coreIndex = distance(posMap.begin(), posMap.find(coreKey));
    posMap.clear();
    int index = coreIndex;
    while (index > 0) {
        auto cur = preVec.at(index);
        index--;
        auto pre = preVec.at(index);
        if (cur.first - (pre.first + pre.second * ssr.length()) == ssr.length()) {
            posMapFil[pre.first] = pre.second;
        } else {
            break;
        }
    }
    index = coreIndex;

    while (index < preVec.size() - 1) {
        auto cur = preVec.at(index);
        index++;
        auto pos = preVec.at(index);
        if (pos.first - (cur.first + cur.second * ssr.length()) == ssr.length()) {
            posMapFil[pos.first] = pos.second;
        } else {
            break;
        }
    }
    preVec.clear();

    std::vector<std::pair<size_t, int>> posVec;
    for (const auto & it : posMapFil) {
        posVec.emplace_back(std::make_pair(it.first, it.second));
    }

    coreIndex = distance(posMapFil.begin(), posMapFil.find(coreKey));
    std::pair<size_t, int> corePair = posVec.at(coreIndex);
    posSSRlenPair = std::make_pair(corePair.first, corePair.second * ssr.length());
    int posIndex = coreIndex;
    while (posIndex > 0) {
        auto cur = posVec.at(posIndex);
        posIndex--;
        auto pre = posVec.at(posIndex);
        if (cur.first - (pre.first + pre.second * ssr.length()) == ssr.length()) {
            auto testSeq = rawStr.substr(pre.first, cur.first - pre.first);
            auto testSeqPos = locVarIt->ff.mStr.find(testSeq);
            if (testSeqPos == std::string::npos) {
                if (locVarIt->ff.mStr.length() >= testSeq.length()) {
                    unsigned int ed = edit_distance(testSeq, locVarIt->ff.mStr.substr(locVarIt->ff.mStr.length() - testSeq.length(), testSeq.length()));
                    if (ed <= locVarIt->edCutoff) {
                        break;
                    } else {
                        posSSRlenPair = std::make_pair(pre.first, posSSRlenPair.first - pre.first + posSSRlenPair.second);
                    }
                } else {
                    posSSRlenPair = std::make_pair(pre.first, posSSRlenPair.first - pre.first + posSSRlenPair.second);
                }
            } else {
                if (testSeqPos == locVarIt->ff.mStr.length() - testSeq.length()) {
                    break;
                } else {
                    posSSRlenPair = std::make_pair(pre.first, posSSRlenPair.first - pre.first + posSSRlenPair.second);
                }
            }

        } else {
            break;
        }
    }

    posIndex = coreIndex;
    while (posIndex < posVec.size() - 1) {
        auto cur = posVec.at(posIndex);
        posIndex++;
        auto suc = posVec.at(posIndex);
        if (suc.first - (cur.first + cur.second * ssr.length()) == ssr.length()) {
            auto testSeq = rawStr.substr(cur.first + cur.second * ssr.length(),
                    suc.first + suc.second * ssr.length() - (cur.first + cur.second * ssr.length()));
            auto testSeqPos = locVarIt->rf.mStr.find(testSeq);

            if (testSeqPos == std::string::npos) {

                if (locVarIt->rf.mStr.length() >= testSeq.length()) {
                    unsigned int ed = edit_distance(testSeq, locVarIt->rf.mStr.substr(0, testSeq.length()));
                    if (ed <= locVarIt->edCutoff) {
                        break;
                    } else {
                        posSSRlenPair = std::make_pair(posSSRlenPair.first, suc.first + suc.second * ssr.length() - posSSRlenPair.first);
                    }
                } else {
                    posSSRlenPair = std::make_pair(posSSRlenPair.first, suc.first + suc.second * ssr.length() - posSSRlenPair.first);
                }

            } else {
                if (testSeqPos == 0) {
                    break;
                } else {
                    posSSRlenPair = std::make_pair(posSSRlenPair.first, suc.first + suc.second * ssr.length() - posSSRlenPair.first);
                }
            }
        } else {
            break;
        }
    }

    posMapFil.clear();
    posVec.clear();

    if (mOptions->mLocVars.locVarOptions.printRes) {
        std::string msg = "rawStr: " + rawStr;
        ss << msg << "\n";
        //cCout(msg, 'r');
        msg = "mra: " + rawStr.substr(posSSRlenPair.first, posSSRlenPair.second);
        ss << msg << "\n";
        //cCout(msg, 'b');
    }
    //posSSRlenPair = std::make_pair(posSSRlenPair.first );
    if (rfTrimPos == 0) {
        posSSRlenPair = std::make_pair(posSSRlenPair.first + ffTrimPos, posSSRlenPair.second);
    } else {
        posSSRlenPair = std::make_pair(posSSRlenPair.first + ffTrimPos, rfTrimPos - (posSSRlenPair.first + ffTrimPos));
    }

    return posSSRlenPair;
}

void SsrScanner::analyzeMRASub(Genotype & genotype, std::map<int, std::string> & subMap) {
    std::vector<size_t> subMapVec;
    genotype.baseLocVar.mra.mStr = locVarIt->mra.mStr;
    size_t si = 0;

    for (const auto & it : subMap) {
        if (it.first >= locVarIt->ff.length() && it.first <= locVarIt->ff.length() + locVarIt->mra.length() - 1) {
            si = (size_t) (it.first - locVarIt->ff.length());
            genotype.baseLocVar.mra.mStr.replace(si, 1, it.second.substr(2, 1));
            subMapVec.emplace_back(si);
        }
    }

    genotype.baseLocVar.name = locVarIt->name + "_";
    if (!subMapVec.empty()) {
        int count = 0;
        std::string tmpStr = "";
        std::string tmpGenotypeStr = "";
        size_t startPosSi = 0;
        for (int i = 0; i < locVarIt->nRep; i++) {
            startPosSi = (size_t) (i * locVarIt->repuit.length());
            tmpStr = genotype.baseLocVar.mra.mStr.substr(startPosSi, locVarIt->repuit.length());
            if (tmpStr == locVarIt->repuit.mStr) {
                count++;
                tmpGenotypeStr = "(" + locVarIt->repuit.mStr + ")" + std::to_string(count);
            } else {
                genotype.baseLocVar.name.append(tmpGenotypeStr);
                tmpGenotypeStr.clear();
                count = 0;
                genotype.baseLocVar.name.append(tmpStr);
            }
        }
        if (count != 0) {
            genotype.baseLocVar.name.append(tmpGenotypeStr);
        }
    } else {
        genotype.baseLocVar.name.append("(" + locVarIt->repuit.mStr + ")" + std::to_string(locVarIt->nRep));
    }
}

bool SsrScanner::checkMRASub(Variance & fVar) {
    bool shouldAnalysis = false;

    for (const auto & it : fVar.subMap) {
        if (it.first >= locVarIt->ff.length() && it.first <= locVarIt->ff.length() + locVarIt->mra.length() - 1) {
            ss << "sub111: " << it.first << " : " << it.second << "\n";
            shouldAnalysis = true;
        }
    }

    for (const auto & it : fVar.insMap) {
        if (it.first >= locVarIt->ff.length() && it.first <= locVarIt->ff.length() + locVarIt->mra.length() - 1) {
            ss << "insert111: " << it.first << " : " << it.second << "\n";
            shouldAnalysis = false;
            //return shouldAnalysis;
        }
    }

    for (const auto & it : fVar.delMap) {
        if (it.first >= locVarIt->ff.length() && it.first <= locVarIt->ff.length() + locVarIt->mra.length() - 1) {
            shouldAnalysis = false;
            ss << "del111: " << it.first << " : " << it.second << "\n";
            //return shouldAnalysis;
        }
    }

    return shouldAnalysis;
}

std::map<std::string, std::map<std::string, Genotype>> SsrScanner::merge(std::vector<std::map<std::string, std::map<std::string, Genotype>>> & totalGenotypeSsrMapVec) {

    std::map<std::string, std::map < std::string, Genotype>> allGenotypeSsrMap;
    if (totalGenotypeSsrMapVec.empty()) {
        return allGenotypeSsrMap;
    }

    std::set<std::string> locSet;
    for (const auto & it : totalGenotypeSsrMapVec) {//std::vector<std::map<std::string, std::map< std::string, Genotype>>>
        for (const auto & it2 : it) {
            locSet.insert(it2.first);
        }
    }

    for (const auto & it : locSet) {
        std::map<std::string, Genotype> eachGenotypeMap;
        for (const auto & it2 : totalGenotypeSsrMapVec) {
            auto it3 = it2.find(it);
            if (it3 != it2.end()) {
                for (const auto & it4 : it3->second) {
                    auto it5 = eachGenotypeMap.find(it4.first);
                    if (it5 == eachGenotypeMap.end()) {
                        eachGenotypeMap[it4.first] = it4.second;
                    } else {
                        it5->second.numReads += it4.second.numReads;
                    }
                }
            }
        }
        allGenotypeSsrMap[it] = eachGenotypeMap;
    }
    return allGenotypeSsrMap;
}

void SsrScanner::merge(std::vector<std::map<std::string, std::map<std::string, int>>> & totalSexLocVec, Options * & mOptions){
    std::map<std::string, int> seqMapX;
    std::map<std::string, int> seqMapY;

    for (auto & it : totalSexLocVec) {
        for (const auto & it2 : it["X"]) {
            seqMapX[it2.first] += it2.second;
        }

        for (const auto & it2 : it["Y"]) {
            seqMapY[it2.first] += it2.second;
        }
    }
    
    totalSexLocVec.clear();
    totalSexLocVec.shrink_to_fit();

    const char* target;
    int targetLength;
    const char* readSeq;
    int readLength;
    std::string tmpStr;
    
    if(!seqMapX.empty()){
        tmpStr = mOptions->mSex.getFullRefX();
        target = tmpStr.c_str();
        targetLength = tmpStr.length();
        for(auto & it : seqMapX){
            if(it.second >= mOptions->mSex.minReadsX){
                mOptions->mSex.readsX += it.second;
                readSeq = it.first.c_str();
                readLength = it.first.length();
                auto snpsMapX = SsrScanner::doSimpleAlignment(mOptions, readSeq, readLength, target, targetLength);
                for (auto & it2 : snpsMapX) {
                    mOptions->mSex.snpsRefX.insert(it2.first);
                }   
                mOptions->mSex.seqVecX.emplace_back(std::make_tuple(it.first, it.second, snpsMapX));
            }
        }

        std::sort(mOptions->mSex.seqVecX.begin(), mOptions->mSex.seqVecX.end(),
                [](const std::tuple<std::string, int, std::map<int, std::string>> &l,
                const std::tuple<std::string, int, std::map<int, std::string>> &r) {
                    return get<1>(l) > get<1>(r);
                });
    }
    
    if(!seqMapY.empty()) {
        tmpStr = mOptions->mSex.getFullRefY();
        target = tmpStr.c_str();
        targetLength = tmpStr.length();
        for(auto & it : seqMapY){
            if(it.second >= mOptions->mSex.minReadsY){
                mOptions->mSex.readsY += it.second;
                readSeq = it.first.c_str();
                readLength = it.first.length();
                auto snpsMapY = SsrScanner::doSimpleAlignment(mOptions, readSeq, readLength, target, targetLength);
                for (auto & it2 : snpsMapY) {
                    mOptions->mSex.snpsRefY.insert(it2.first);
                }   
                mOptions->mSex.seqVecY.emplace_back(std::make_tuple(it.first, it.second, snpsMapY));
            }
        }
        std::sort(mOptions->mSex.seqVecY.begin(), mOptions->mSex.seqVecY.end(),
                [](const std::tuple<std::string, int, std::map<int, std::string>> &l,
                const std::tuple<std::string, int, std::map<int, std::string>> &r) {
                    return get<1>(l) > get<1>(r);
                });
    }
}

//std::vector<std::map<std::string, std::vector<std::pair<std::string, Genotype>>>> SsrScanner::merge(std::vector<std::map<std::string, std::map<std::string, Genotype>>> & totalGenotypeSsrMapVec) {
//
//    std::vector<std::map<std::string, std::vector<std::pair<std::string, Genotype>>>> sortedGenotypeSsrMapVec;
//    std::map<std::string, std::map<std::string, Genotype>> finalGenotypeMap;
//    if (totalGenotypeSsrMapVec.empty()) {
//        return sortedGenotypeSsrMapVec;
//    } 
//
//    std::set<std::string> locSet;
//    for (const auto & it : totalGenotypeSsrMapVec) {//std::vector<std::map<std::string, std::map< std::string, Genotype>>>
//        for (const auto & itt : it) {//std::map<std::string, std::map< std::string, Genotype>>
//            locSet.insert(itt.first);
//        }
//    }
//
//    for (const auto & it : locSet) {
//        std::set<std::string> seqSet;
//        for (const auto & its : totalGenotypeSsrMapVec) {//std::vector<std::map<std::string, std::map< std::string, Genotype>>>
//            auto itt = its.find(it); //std::map<std::string, std::map< std::string, Genotype>>
//            if (itt != its.end()) {
//                for (const auto & itf : itt->second) {
//                    seqSet.insert(itf.first);
//                }
//            }
//        }
//
//        std::map<std::string, Genotype> tmpFinalGenotypeMap;
//        for (const auto & it2 : seqSet) {
//            //tmpFinalGenotypeMap[its] = tmpGenotype();
//            Genotype tmpFinalGenotype;
//            for (const auto & it3 : list) {////std::vector<std::map<std::string, std::map< std::string, Genotype>>>
//                auto it4 = it3.find(it); //std::map<std::string, std::map< std::string, Genotype>>
//                if (it4 != it3.end()) {
//                    auto it5 = it4->second.find(it2); //// std::map< std::string, Genotype>>
//                    if (it5 != it4->second.end()) {
//                        if (tmpFinalGenotype.baseLocVar.mra.mStr.empty()) {
//                            tmpFinalGenotype.baseLocVar = it5->second.baseLocVar;
//                        } else {
//                            tmpFinalGenotype.numReads += it5->second.numReads;
//                        }
//                    }
//                }
//            }
//            tmpFinalGenotypeMap[it2] = tmpFinalGenotype;
//        }
//        finalGenotypeMap[it] = tmpFinalGenotypeMap;
//    }
//
//    locSet.clear();
//    list.clear();
//    list.shrink_to_fit();
//    return finalGenotypeMap;
//}

std::vector<std::map<std::string, std::vector<std::pair<std::string, Genotype>>>> SsrScanner::report(Options * & mOptions, std::map<std::string, std::map<std::string, Genotype>> &allGenotypeMap) {

    std::map<std::string, std::vector<std::pair < std::string, Genotype>>> sortedAllGenotypeMap;
    std::map<std::string, std::vector<std::pair < std::string, Genotype>>> sortedAllGenotypeMraMap;
    std::vector<std::map < std::string, std::vector<std::pair < std::string, Genotype>>>> sortedAllGenotypeMapVec{sortedAllGenotypeMap, sortedAllGenotypeMraMap};
    const char* target;
    int targetLength;
    const char* readSeq;
    int readLength;
    for (auto & it : allGenotypeMap) {//marker, seq, genotypes;
        auto locVarIt = &(mOptions->mLocVars.refLocMap.find(it.first)->second);
        int maxReads = 0;
        for (const auto & it2 : it.second) {
            if (it2.second.numReads > maxReads) {
                maxReads = it2.second.numReads;
            }
        }

        std::set<int> ffSet;
        std::set<int> rfSet;
        std::vector<std::pair < std::string, Genotype>> tmpVec;
        tmpVec.reserve(it.second.size());
        std::vector<std::pair<std::string, Genotype>> tmpMraVec;
        tmpMraVec.reserve(it.second.size());
        for (auto & it2 : it.second) {
            if (it2.second.numReads >= mOptions->mLocVars.locVarOptions.minSeqs &&
                    (it2.second.numReads * 100 / maxReads) >= mOptions->mLocVars.locVarOptions.minSeqsPer && 
                    it2.second.baseLocVar.mra.mStr.length() >= (it2.second.baseLocVar.repuit.mStr.length() + it2.second.baseLocVar.repuit2.mStr.length()) * mOptions->mLocVars.locVarOptions.minNSSRUnit) {
                
                auto tmpGeno = getGenotype(it2.second.baseLocVar.mra.mStr, locVarIt->repuit.mStr, locVarIt->repuit2.mStr);
                auto maxNRep = getMaxNumRep(tmpGeno);
                if (maxNRep >= mOptions->mLocVars.locVarOptions.minNSSRUnit) {
                    target = locVarIt->ff.mStr.c_str();
                    targetLength = locVarIt->ff.mStr.length();
                    readSeq = it2.second.baseLocVar.ff.mStr.c_str();
                    readLength = it2.second.baseLocVar.ff.mStr.length();
                    it2.second.baseLocVar.snpsMapff = SsrScanner::doSimpleAlignment(mOptions, readSeq, readLength, target, targetLength);
                    for (auto & it3 : it2.second.baseLocVar.snpsMapff) {
                        ffSet.insert(it3.first);
                    }           
                    target = locVarIt->rf.mStr.c_str();
                    targetLength = locVarIt->rf.mStr.length();
                    readSeq = it2.second.baseLocVar.rf.mStr.c_str();
                    readLength = it2.second.baseLocVar.rf.mStr.length();
                    it2.second.baseLocVar.snpsMaprf = SsrScanner::doSimpleAlignment(mOptions, readSeq, readLength, target, targetLength);

                    for (auto & it3 : it2.second.baseLocVar.snpsMaprf) {
                        rfSet.insert(it3.first);
                    }
                    mOptions->mLocVars.refLocMap[it.first] = *locVarIt;
                    it2.second.baseLocVar.mraName = tmpGeno;
                    it2.second.baseLocVar.repuitAll = locVarIt->repuitAll;
                    it2.second.baseLocVar.mraBase = getMraBase(it2.second.baseLocVar.mraName);
                    tmpVec.push_back(std::make_pair(it2.first, it2.second));
                    tmpMraVec.push_back(std::make_pair(it2.first, it2.second));
                }
            }
        }
        mOptions->mLocVars.refLocMap[it.first].refSnpsSetffMap[basename(mOptions->prefix)] = ffSet;
        mOptions->mLocVars.refLocMap[it.first].refSnpsSetrfMap[basename(mOptions->prefix)] = rfSet;
        ffSet.clear();
        rfSet.clear();
        std::sort(tmpVec.begin(), tmpVec.end(),
                [](const std::pair<std::string, Genotype> & l, const std::pair<std::string, Genotype> & r) {
                    if (l.second.baseLocVar.effectiveLen != r.second.baseLocVar.effectiveLen) {
                        return l.second.baseLocVar.effectiveLen < r.second.baseLocVar.effectiveLen;
                    }
                    return l.first.length() < r.first.length();
                });
        sortedAllGenotypeMapVec.at(0)[it.first] = tmpVec;
        tmpVec.clear();
        tmpVec.shrink_to_fit();

        std::sort(tmpMraVec.begin(), tmpMraVec.end(),
                [](const std::pair<std::string, Genotype> & l, const std::pair<std::string, Genotype> & r) {
                    if (l.second.baseLocVar.mra.mStr.length() != r.second.baseLocVar.mra.mStr.length()) {
                        return l.second.baseLocVar.mra.mStr.length() < r.second.baseLocVar.mra.mStr.length();
                    }
                    return l.first.length() < r.first.length();
                });

        sortedAllGenotypeMapVec.at(1)[it.first] = tmpMraVec;
        tmpMraVec.clear();
        tmpMraVec.shrink_to_fit();
    }
    allGenotypeMap.clear();
    
//    std::string foutName = mOptions->prefix + "_ssr_genotypes.txt";
//    std::ofstream * fout = new std::ofstream();
//    fout->open(foutName.c_str(), std::ofstream::out);
//
//    if (!fout->is_open()) error_exit("Can not open output file: " + foutName);
//    if (mOptions->verbose) loginfo("Starting to write genotype table!");
//
//    *fout << "#Locus\tMicrosatellite\tMRABase\tMRAName\tMRASize\tGenotype\tNumReads\tFF\tMRA\tRF\tSnpsFF\tSnpsRF\n";
//
//    for ( auto & it : sortedAllGenotypeMapVec.at(0)) {
//        auto locVarIt = mOptions->mLocVars.refLocMap[it.first]; 
//        for (auto & it2 : it.second) {            
//            *fout << it.first << "\t" << it2.second.baseLocVar.repuitAll.mStr << "\t" <<
//                    it2.second.baseLocVar.mraBase << "\t" <<
//                    it2.second.baseLocVar.mraName << "\t" << it2.second.baseLocVar.mra.mStr.length() << "\t" <<
//                    it2.second.baseLocVar.effectiveLen << "\t" << it2.second.numReads << "\t" <<
//                    //(it2.second.baseLocVar.totalReads > 0 ? "Y" : "N") << "\t" <<
//                    it2.second.baseLocVar.ff.mStr << "\t" << it2.second.baseLocVar.mra.mStr << "\t"
//                    << it2.second.baseLocVar.rf.mStr << "\t";
//            if (!locVarIt.refSnpsSetffMap[basename(mOptions->prefix)].empty()) {
//                for (const auto & snps : locVarIt.refSnpsSetffMap[basename(mOptions->prefix)]) {
//                    *fout << snps << locVarIt.effectiveSeq.mStr[snps] << "|" << it2.second.baseLocVar.effectiveSeq.mStr[snps] << ";";
//                }
//            } else {
//                *fout << "NA";
//            }
//            *fout << "\t";
//
//            if (!locVarIt.refSnpsSetrfMap[basename(mOptions->prefix)].empty()) {
//                for (const auto & snps : locVarIt.refSnpsSetrfMap[basename(mOptions->prefix)]) {
//                    *fout << snps << locVarIt.effectiveSeq.mStr[snps] << "|" << it2.second.baseLocVar.effectiveSeq.mStr[snps] << ";";
//                }
//            } else {
//                *fout << "NA";
//            }
//            *fout << "\n";
//        }
//    }
//
//    fout->flush();
//    fout->clear();
//    if (fout) {
//        delete fout;
//        fout = NULL;
//    }

    std::string foutName = mOptions->prefix + "_genotypes_mra.txt";
    std::ofstream* fout = new std::ofstream();
    fout->open(foutName.c_str(), std::ofstream::out);

    if (!fout->is_open()) error_exit("Can not open output file: " + foutName);
    if (mOptions->verbose) loginfo("Starting to write genotype table!");

    *fout << "#Sample\tLocus\tMicrosatellite\tMRABase\tMRAName\tMRASize\tAllele\tNumReads\tPutativeAllele\tFF\tMRA\tRF\tSnpsFF\tSnpsRF\n";
    for (auto & it : sortedAllGenotypeMapVec.at(1)) {//marker, seq, genotype
        auto locVarIt = &(mOptions->mLocVars.refLocMap[it.first]);
        std::map<int, int> genoReadsMap;
        std::map<int, int> genoReadsMap2;
        std::map<int, int> genoReadsMapT;
    
        for(const auto & it2 : it.second){
            genoReadsMap[it2.second.baseLocVar.effectiveLen] += it2.second.numReads;
        }
        
        if(genoReadsMap.size() == 1){//only one peak;
            genoReadsMapT[genoReadsMap.begin()->first] = genoReadsMap.begin()->second;
        } else if (genoReadsMap.size() > 1){// > 1 peak;
            auto maxGenoReads = getMaxKeyValue(genoReadsMap);//key max value pair;
            genoReadsMapT[maxGenoReads.first] = maxGenoReads.second;

            int occ = 0;
            for (const auto & it2 : genoReadsMap) {
                if (it2.second == maxGenoReads.second) {
                    occ++;
                }
            }

            if (occ == 1) {
                for (const auto & it2 : genoReadsMap) {
                    if (it2.first > maxGenoReads.first) {
                        genoReadsMap2[it2.first] = it2.second;
                    }
                }

                if (!genoReadsMap2.empty()) {
                    auto maxGenoReads2 = getMaxKeyValue(genoReadsMap2, true);
                    int repul = locVarIt->repuit.mStr.length();
                    if (maxGenoReads2.first - maxGenoReads.first == repul) {
                        if ((double) maxGenoReads2.second / maxGenoReads.second >= mOptions->mLocVars.locVarOptions.hlRatio1) {
                            genoReadsMapT[maxGenoReads2.first] = maxGenoReads2.second;
                        }
                    } else if (maxGenoReads2.first - maxGenoReads.first == repul * 2) {
                        if ((double) maxGenoReads2.second / maxGenoReads.second >= mOptions->mLocVars.locVarOptions.hlRatio2) {
                            genoReadsMapT[maxGenoReads2.first] = maxGenoReads2.second;
                        }
                    } else if (maxGenoReads2.first - maxGenoReads.first == repul * 3) {
                        if ((double) maxGenoReads2.second / maxGenoReads.second >= mOptions->mLocVars.locVarOptions.hlRatio2 / 2) {
                            genoReadsMapT[maxGenoReads2.first] = maxGenoReads2.second;
                        }
                    } else {
                        if ((double) maxGenoReads2.second / maxGenoReads.second >= mOptions->mLocVars.locVarOptions.hlRatio2 / 4) {
                            genoReadsMapT[maxGenoReads2.first] = maxGenoReads2.second;
                        }
                    }
                }
                
                genoReadsMap2.clear();
                
            } else if (occ > 1) {
                for (const auto & it2 : genoReadsMap) {
                    if (it2.second == maxGenoReads.second && it2.first != maxGenoReads.first) {
                        genoReadsMapT[it2.first] = maxGenoReads.second;
                        break;
                    }
                }
            }
        }
        
        genoReadsMap.clear();
        
        for ( auto & it2 : it.second) {
            auto tmpVar = genoReadsMapT.find(it2.second.baseLocVar.effectiveLen);
            it2.second.baseLocVar.totalReads = (tmpVar == genoReadsMapT.end() ? 0 : tmpVar->second);
            
//            std::cout << "it2.second.baseLocVar.totalReads: " << it2.second.baseLocVar.totalReads << "\n";
            *fout << mOptions->prefix << "\t" << it.first << "\t" << it2.second.baseLocVar.repuitAll.mStr << "\t" <<
                    it2.second.baseLocVar.mraBase << "\t" <<
                    it2.second.baseLocVar.mraName << "\t" << it2.second.baseLocVar.mra.mStr.length() << "\t" <<
                    it2.second.baseLocVar.effectiveLen << "\t" << it2.second.numReads << "\t" <<
                    (it2.second.baseLocVar.totalReads > 0 ? "Y" : "N") << "\t" <<
                    it2.second.baseLocVar.ff.mStr << "\t" << it2.second.baseLocVar.mra.mStr << "\t"
                    << it2.second.baseLocVar.rf.mStr << "\t";
            if (!locVarIt->refSnpsSetffMap[basename(mOptions->prefix)].empty()) {
                for (const auto & snps : locVarIt->refSnpsSetffMap[basename(mOptions->prefix)]) {
                    *fout << snps << locVarIt->effectiveSeq.mStr[snps] << "|" << it2.second.baseLocVar.effectiveSeq.mStr[snps] << ";";
                }
            } else {
                *fout << "NA";
            }
            *fout << "\t";

            if (!locVarIt->refSnpsSetrfMap[basename(mOptions->prefix)].empty()) {
                for (const auto & snps : locVarIt->refSnpsSetrfMap[basename(mOptions->prefix)]) {
                    *fout << snps << locVarIt->effectiveSeq.mStr[snps] << "|" << it2.second.baseLocVar.effectiveSeq.mStr[snps] << ";";
                }
            } else {
                *fout << "NA";
            }
            *fout << "\n";
        }
        genoReadsMapT.clear();
    }

    fout->flush();
    fout->clear();
    fout->close();
    if (fout) {
        delete fout;
        fout = NULL;
    }
    if (mOptions->verbose) loginfo("End writing genotype table!");

    if (!mOptions->mSex.sexMarker.empty()) {
        
        if(mOptions->mSex.readsX != 0){
            mOptions->mSex.YXRatio = std::round(((double) mOptions->mSex.readsY / (double) mOptions->mSex.readsX) * 100.0) / 100.0;
            if(mOptions->mSex.YXRationCuttoff < mOptions->mSex.YXRatio){
                if(mOptions->mSex.minTotalReadsX < mOptions->mSex.readsX){
                    if(mOptions->mSex.minTotalReadsY < mOptions->mSex.readsY){
                        mOptions->mSex.sexMF = "Male"; 
                    } else {
                        mOptions->mSex.sexMF = "Female"; 
                    }
                } else {
                    mOptions->mSex.sexMF = "Inconclusive"; 
                }
                
            } else {
                if (mOptions->mSex.minTotalReadsX < mOptions->mSex.readsX) {
                    mOptions->mSex.sexMF = "Female";
                } else {
                    mOptions->mSex.YXRatio = 0;
                    mOptions->mSex.sexMF = "Inconclusive";
                }
            }
        } else {
            mOptions->mSex.YXRatio = 0;
            mOptions->mSex.sexMF = "Inconclusive";
        }
        
        //std::cout << "sexLoc: " << mOptions->mSex.readsY << " : " << mOptions->mSex.readsX << " -> " << mOptions->mSex.YXRatio << "\n";
        
        std::string foutName = mOptions->prefix + "_sex_loc_id.txt";
        std::ofstream* fout = new std::ofstream();
        fout->open(foutName.c_str(), std::ofstream::out);

        if (!fout->is_open()) error_exit("Can not open output file: " + foutName);
        if (mOptions->verbose) loginfo("Starting to write sex identification loc file!");
        
        *fout << "#SexLoc\tNumReadsX\tNumReadsY\tRatio\tPutativeSex\tAlleleX\tSnpsX\tAlleleY\tSnpsY\tNote\n";
        
        *fout << mOptions->mSex.sexMarker << "\t" << mOptions->mSex.readsX << "\t" << mOptions->mSex.readsY << "\t" << mOptions->mSex.YXRatio << "\t" << 
                mOptions->mSex.sexMF << "\t" << mOptions->mSex.getFullRefX() << "\t";
        if(mOptions->mSex.snpsRefX.empty()){
            *fout << "NA\t";
        } else {
            for(const auto & its : mOptions->mSex.snpsRefX){
                *fout << its << ";";
            }
            *fout << "\t";
        }
        
        *fout << mOptions->mSex.getFullRefY() << "\t";

        if (mOptions->mSex.snpsRefY.empty()) {
            *fout << "NA\t";
        } else {
            for (const auto & its : mOptions->mSex.snpsRefY) {
                *fout << its << ";";
            }
            *fout << "\t";
        }
        
        *fout << "total\n";

        std::map<int, std::string> tmpSnpMap;
        if (!mOptions->mSex.seqVecX.empty()){
            for (const auto & its : mOptions->mSex.seqVecX) {
                *fout << "X\t" << get<1>(its) << "\t0\t0\tNA\t" << get<0>(its) << "\t";
                tmpSnpMap = get<2>(its);
                if (tmpSnpMap.empty()) {
                    *fout << "NA";
                } else {
                    for (const auto & itm : tmpSnpMap) {
                        *fout << itm.first << mOptions->mSex.getFullRefX()[itm.first] << "|" << get<0>(its)[itm.first] << ";";
                    }
                }
                *fout << "\tNA\tNA\teach\n";
            }
        }

        if (!mOptions->mSex.seqVecY.empty()) {
            for (const auto & its : mOptions->mSex.seqVecY) {
                *fout << "Y\t0\t" << get<1>(its) << "\t0\tNA\tNA\tNA\t" << get<0>(its) << "\t";
                tmpSnpMap = get<2>(its);
                if (tmpSnpMap.empty()) {
                    *fout << "NA";
                } else {
                    for (const auto & itm : tmpSnpMap) {
                        *fout << itm.first << mOptions->mSex.getFullRefY()[itm.first] << "|" << get<0>(its)[itm.first] << ";";
                    }
                }
                *fout << "\teach\n";
            }
        }
        
        tmpSnpMap.clear();
        
        fout->flush();
        fout->clear();
        fout->close();
        if(fout){
            delete fout;
            fout = NULL;
        }
    }

    return sortedAllGenotypeMapVec;
}

std::vector<std::pair<std::string, Genotype>> SsrScanner::sortGenotypeMap(std::map<std::string, Genotype> & genoMap) {
    std::vector<std::pair < std::string, Genotype>> sVec;

    std::copy(genoMap.begin(), genoMap.end(),
            std::back_inserter<std::vector<std::pair < std::string, Genotype>>>(sVec));

    std::sort(sVec.begin(), sVec.end(),
            [](const std::pair<std::string, Genotype> & l, const std::pair<std::string, Genotype> & r) {
                if (l.second.baseLocVar.effectiveLen != r.second.baseLocVar.effectiveLen) {
                    return l.second.baseLocVar.effectiveLen < r.second.baseLocVar.effectiveLen;
                }
                return l.first.length() < r.first.length();
            });

    return sVec;
}

std::size_t SsrScanner::mutationMatch(std::string & target, std::string query, bool rev) {
    std::string orig = query;
    std::size_t pos = 0;
    for (int i = 0; i < query.length(); i++) {
        for (int j = 0; j < 3; j++) {
            if (query[i] != ATCG_BASES[j]) {
                query[i] = ATCG_BASES[j];
                pos = rev ? target.rfind(query) : target.find(query);
                if (pos != std::string::npos) {
                    return pos;
                }
                query = orig;
            }
        }
    }
    return 0;
}

void SsrScanner::preAnalyze(Read* & r1, std::size_t & ffpos, std::size_t & rfpos, bool & mraAnalyze) {
    if (ffpos != std::string::npos && rfpos != std::string::npos) {
        if(ffpos + locVarIt->ff.mStr.length() + mOptions->mLocVars.locVarOptions.minNSSRUnit * (locVarIt->repuit.mStr.length() + locVarIt->repuit2.mStr.length()) <= rfpos) {
            std::pair<size_t, int> tmp = analyzeMRA(r1->mSeq.mStr, locVarIt->repuit.mStr, ffpos, rfpos);
            if (tmp.second >= mOptions->mLocVars.locVarOptions.minNSSRUnit * locVarIt->repuit.mStr.length()) {
                ffpos += locVarIt->ff.mStr.length();
            } else {
                mraAnalyze = true;
            }
        } else {
            ffpos = 0;
            rfpos = r1->mSeq.length() - 1;
            mraAnalyze = true;
        }
    } else if (ffpos == std::string::npos && rfpos != std::string::npos) {
        if (locVarIt->ff.mStr.empty()) {
            ffpos = 0;
            if (ffpos + mOptions->mLocVars.locVarOptions.minNSSRUnit * (locVarIt->repuit.mStr.length() + locVarIt->repuit2.mStr.length()) > rfpos) {
                rfpos = r1->mSeq.mStr.length() - 1;
                mraAnalyze = true;
            }
        } else {
            std::string minff = locVarIt->ff.mStr.substr(locVarIt->ff.mStr.length() - std::min(mOptions->mLocVars.locVarOptions.minMatchesFR, static_cast<int>(locVarIt->ff.mStr.length())));
            ffpos = r1->mSeq.mStr.find(minff);
            if (ffpos != std::string::npos) {
                
                if(ffpos + minff.length() + mOptions->mLocVars.locVarOptions.minNSSRUnit * (locVarIt->repuit.mStr.length() + locVarIt->repuit2.mStr.length()) <= rfpos){
                    ffpos += minff.length();
                } else {
                    ffpos = 0;
                    mraAnalyze = true;
                }
            } else {
                ffpos = mutationMatch(r1->mSeq.mStr, minff);
                if (ffpos == 0) {
                    mraAnalyze = true;
                } else {
                    if (ffpos + minff.length() + mOptions->mLocVars.locVarOptions.minNSSRUnit * (locVarIt->repuit.mStr.length() + locVarIt->repuit2.mStr.length()) <= rfpos) {
                        ffpos += minff.length();
                    } else {
                        ffpos = 0;
                        mraAnalyze = true;
                    }
                }
            }
        }
    } else if (ffpos != std::string::npos && rfpos == std::string::npos) {
        if(ffpos + locVarIt->ff.mStr.length() + mOptions->mLocVars.locVarOptions.minNSSRUnit * (locVarIt->repuit.mStr.length() + locVarIt->repuit2.mStr.length()) <= r1->mSeq.mStr.length()) {
            ffpos += locVarIt->ff.mStr.length();
            if (locVarIt->rf.mStr.empty()) {
                rfpos = r1->mSeq.mStr.length() - 1;
            } else {
                std::string minrf = locVarIt->rf.mStr.substr(0, mOptions->mLocVars.locVarOptions.minMatchesFR);
                rfpos = r1->mSeq.mStr.rfind(minrf);
                if (rfpos != std::string::npos) {
                    if(ffpos + mOptions->mLocVars.locVarOptions.minNSSRUnit * (locVarIt->repuit.mStr.length() + locVarIt->repuit2.mStr.length()) > rfpos){
                        rfpos = r1->mSeq.mStr.length() - 1;
                        mraAnalyze = true;
                    }
                } else {
                    rfpos = mutationMatch(r1->mSeq.mStr, minrf, true);
                    if (rfpos == 0) {
                        mraAnalyze = true;
                    }
                }
            }
        } else {
            ffpos = 0;
            mraAnalyze = true;
        }
    } else if (ffpos == std::string::npos && rfpos == std::string::npos) {

        if (locVarIt->ff.mStr.empty() && locVarIt->rf.mStr.empty()) {
            ffpos = 0;
            rfpos = r1->mSeq.mStr.length() - 1;
        } else if (!locVarIt->ff.mStr.empty() && locVarIt->rf.mStr.empty()) {
            rfpos = r1->mSeq.mStr.length() - 1;
            std::string minff = locVarIt->ff.mStr.substr(locVarIt->ff.mStr.length() - std::min(mOptions->mLocVars.locVarOptions.minMatchesFR, static_cast<int>(locVarIt->ff.mStr.length())));
            ffpos = r1->mSeq.mStr.find(minff);
            if (ffpos != std::string::npos) {
               if(ffpos + minff.length() + mOptions->mLocVars.locVarOptions.minNSSRUnit * (locVarIt->repuit.mStr.length() + locVarIt->repuit2.mStr.length()) <= rfpos){
                   ffpos += minff.length();
               } else {
                   ffpos = 0;
                   mraAnalyze = true;
               }
            } else {
                ffpos = mutationMatch(r1->mSeq.mStr, minff);
                if (ffpos == 0) {
                    mraAnalyze = true;
                } else {
                    if (ffpos + minff.length() + mOptions->mLocVars.locVarOptions.minNSSRUnit * (locVarIt->repuit.mStr.length() + locVarIt->repuit2.mStr.length()) <= rfpos) {
                        ffpos += minff.length();
                    } else {
                        ffpos = 0;
                        mraAnalyze = true;
                    }
                }
            }
        } else if (locVarIt->ff.mStr.empty() && !locVarIt->rf.mStr.empty()) {
            ffpos = 0;
            std::string minrf = locVarIt->rf.mStr.substr(0, mOptions->mLocVars.locVarOptions.minMatchesFR);
            rfpos = r1->mSeq.mStr.rfind(minrf);
            if (rfpos != std::string::npos) {
                
                if(ffpos + mOptions->mLocVars.locVarOptions.minNSSRUnit * (locVarIt->repuit.mStr.length() + locVarIt->repuit2.mStr.length()) > rfpos){
                    rfpos = r1->mSeq.mStr.length() - 1;
                    mraAnalyze = true;
                } 
                
            } else {
                rfpos = mutationMatch(r1->mSeq.mStr, minrf, true);
                if (rfpos == 0) {
                    mraAnalyze = true;
                } else {
                    if (ffpos + mOptions->mLocVars.locVarOptions.minNSSRUnit * (locVarIt->repuit.mStr.length() + locVarIt->repuit2.mStr.length()) > rfpos) {
                        rfpos = r1->mSeq.mStr.length() - 1;
                        mraAnalyze = true;
                    }
                }
            }

        } else if (!locVarIt->ff.mStr.empty() && !locVarIt->rf.mStr.empty()) {
            
            std::string minff = locVarIt->ff.mStr.substr(locVarIt->ff.mStr.length() - std::min(mOptions->mLocVars.locVarOptions.minMatchesFR, static_cast<int>(locVarIt->ff.mStr.length())));
            ffpos = r1->mSeq.mStr.find(minff);

            std::string minrf = locVarIt->rf.mStr.substr(0, mOptions->mLocVars.locVarOptions.minMatchesFR);
            rfpos = r1->mSeq.mStr.rfind(minrf);

            if (ffpos != std::string::npos && rfpos != std::string::npos) {
                if(ffpos + minff.length() + mOptions->mLocVars.locVarOptions.minNSSRUnit * (locVarIt->repuit.mStr.length() + locVarIt->repuit2.mStr.length()) <= rfpos){
                    ffpos += minff.length();
                } else {
                    ffpos = 0;
                    mraAnalyze = true;
                }
            } else if (ffpos != std::string::npos && rfpos == std::string::npos) {
                if(ffpos + minff.length() + mOptions->mLocVars.locVarOptions.minNSSRUnit * (locVarIt->repuit.mStr.length() + locVarIt->repuit2.mStr.length()) <= rfpos) {
                    ffpos += minff.length();
                    rfpos = mutationMatch(r1->mSeq.mStr, minrf, true);
                    if (rfpos == 0) {
                        mraAnalyze = true;
                    } else {
                        if(ffpos + mOptions->mLocVars.locVarOptions.minNSSRUnit * (locVarIt->repuit.mStr.length() + locVarIt->repuit2.mStr.length()) > rfpos){
                            rfpos = r1->mSeq.mStr.length() - 1;
                            mraAnalyze = true;
                        }
                    }
                } else {
                    ffpos = 0;
                    mraAnalyze = true;
                }
            } else if (ffpos == std::string::npos && rfpos != std::string::npos) {
                ffpos = mutationMatch(r1->mSeq.mStr, minff);
                if (ffpos == 0) {
                    mraAnalyze = true;
                } else {
                    if(ffpos + minff.length() + mOptions->mLocVars.locVarOptions.minNSSRUnit * (locVarIt->repuit.mStr.length() + locVarIt->repuit2.mStr.length()) <= rfpos){
                        ffpos += minff.length();
                    } else {
                        ffpos = 0;
                        mraAnalyze = true;
                    }
                }
            } else if (ffpos == std::string::npos && rfpos == std::string::npos) {

                ffpos = mutationMatch(r1->mSeq.mStr, minff);
                rfpos = mutationMatch(r1->mSeq.mStr, minrf, true);
                if (ffpos == 0 && rfpos == 0) {
                    rfpos = r1->mSeq.mStr.length();
                    mraAnalyze = true;
                } else if(ffpos == 0 && rfpos != 0) {
                    if (ffpos + minff.length() + mOptions->mLocVars.locVarOptions.minNSSRUnit * (locVarIt->repuit.mStr.length() + locVarIt->repuit2.mStr.length()) > rfpos) {
                        rfpos = r1->mSeq.mStr.length();
                        mraAnalyze = true;
                    }
                } else if(ffpos != 0 && rfpos == 0){
                    if (ffpos + minff.length() + mOptions->mLocVars.locVarOptions.minNSSRUnit * (locVarIt->repuit.mStr.length() + locVarIt->repuit2.mStr.length()) > rfpos) {
                        ffpos = 0;
                        rfpos = r1->mSeq.mStr.length();
                        mraAnalyze = true;
                    }
                } else {
                    if (ffpos + minff.length() + mOptions->mLocVars.locVarOptions.minNSSRUnit * (locVarIt->repuit.mStr.length() + locVarIt->repuit2.mStr.length()) > rfpos) {
                        ffpos = 0;
                        rfpos = r1->mSeq.mStr.length();
                        mraAnalyze = true;
                    }
                }
            }
        }
    }
}

std::string SsrScanner::scanVar(Read* & r1) {
    ss.str("");
    readSeq = r1->mSeq.mStr.c_str();
    readLength = r1->mSeq.length();
    readName = r1->mName;
    returnedlocus.clear();
    
    if(mOptions->mSex.sexMarker.empty()){
        checkLoci = true;
    } else {
        fpData = mOptions->mSex.primerF.mStr.c_str();
        rpData = mOptions->mSex.primerR.mStr.c_str();
        bool fpMatched = true;
        bool rpMatched = true;
        int fpMismatches = 0;
        int rpMismatches = 0;
        
        for(int i = 0; i < mOptions->mSex.primerF.length(); i++){
            if (fpData[i] != readSeq[i]) {
                fpMismatches++;
            }
            
            if(fpMismatches > mOptions->mSex.mismatchesPF){
                fpMatched = false;
                break;
            }
        }

        if(fpMatched){
            int rlen = r1->length() - mOptions->mSex.primerR.length();
            for(int i = 0; i < mOptions->mSex.primerR.length(); i++){
                if (rpData[i] != readSeq[rlen + i]){
                    rpMismatches++;
                }
                if(rpMismatches > mOptions->mSex.mismatchesPR){
                    rpMatched = false;
                    break;
                }
            }
            if(rpMatched){
                
                if(mOptions->mSex.lengthEqual){
                    unsigned int edx = edit_distance(mOptions->mSex.refX.mStr,
                                r1->mSeq.mStr.substr(mOptions->mSex.primerF.length(),
                                (r1->mSeq.length() - mOptions->mSex.primerF.length() - mOptions->mSex.primerR.length())));
                    
                    unsigned int edy = edit_distance(mOptions->mSex.refY.mStr,
                                r1->mSeq.mStr.substr(mOptions->mSex.primerF.length(),
                                (r1->mSeq.length() - mOptions->mSex.primerF.length() - mOptions->mSex.primerR.length())));
                    
                    if(edx == edy){
                        //checkLoci = true;
                        returnedlocus = mOptions->mSex.sexMarker + "_failed";
                    } else {
                        unsigned int edmin = std::min(edx, edy);
                        if (edmin == edx) {
                            if (edx > mOptions->mSex.mismatchesRX) {
                                //checkLoci = true;
                                returnedlocus = mOptions->mSex.sexMarker + "_failed";
                            } else {
                                checkLoci = false;
                                //tmpSex.readsX++;
                                //tmpSex.seqMapX[r1->mSeq.mStr]++;
                                tmpSexMap["X"][r1->mSeq.mStr]++;
                                //returnedlocus = mOptions->mSex.sexMarker;
                                returnedlocus = mOptions->mSex.sexMarker + "_true";
                                
                            }
                        } else {
                            if (edy > mOptions->mSex.mismatchesRY) {
                                //checkLoci = true;
                                returnedlocus = mOptions->mSex.sexMarker + "_failed";
                            } else {
                                checkLoci = false;
                                //tmpSex.readsY++;
                                //tmpSex.seqMapY[r1->mSeq.mStr]++;
                                tmpSexMap["Y"][r1->mSeq.mStr]++;
                                returnedlocus = mOptions->mSex.sexMarker + "_true";
                                //returnedlocus = mOptions->mSex.sexMarker;
                            }
                        }
                    }
                } else {
                    if (r1->mSeq.length() - mOptions->mSex.primerF.length() - mOptions->mSex.primerR.length() == mOptions->mSex.refX.length()) {
                        unsigned int ed = edit_distance(mOptions->mSex.refX.mStr,
                                r1->mSeq.mStr.substr(mOptions->mSex.primerF.length(),
                                (r1->mSeq.length() - mOptions->mSex.primerF.length() - mOptions->mSex.primerR.length())));
                        if (ed > mOptions->mSex.mismatchesRX) {
                            //checkLoci = true;
                            returnedlocus = mOptions->mSex.sexMarker + "_failed";
                        } else {
                            checkLoci = false;
                            //tmpSex.readsX++;
                            //tmpSex.seqMapX[r1->mSeq.mStr]++;
                            tmpSexMap["X"][r1->mSeq.mStr]++;
                            returnedlocus = mOptions->mSex.sexMarker + "_true";
                            //returnedlocus = mOptions->mSex.sexMarker;
                        }
                    } else if (r1->mSeq.length() - mOptions->mSex.primerF.length() - mOptions->mSex.primerR.length() == mOptions->mSex.refY.length()) {
                        unsigned int ed = edit_distance(mOptions->mSex.refY.mStr,
                                r1->mSeq.mStr.substr(mOptions->mSex.primerF.length(),
                                (r1->mSeq.length() - mOptions->mSex.primerF.length() - mOptions->mSex.primerR.length())));
                        if (ed > mOptions->mSex.mismatchesRY) {
                            //checkLoci = true;
                            returnedlocus = mOptions->mSex.sexMarker + "_failed";
                        } else {
                            checkLoci = false;
                            //tmpSex.readsY++;
                            //tmpSex.seqMapY[r1->mSeq.mStr]++;
                            tmpSexMap["Y"][r1->mSeq.mStr]++;
                            returnedlocus = mOptions->mSex.sexMarker + "_true";
                            //returnedlocus = mOptions->mSex.sexMarker;
                        }
                    } else {
                        //checkLoci = true;
                        returnedlocus = mOptions->mSex.sexMarker + "_failed";
                    }
                }
                
                
            } else {
                checkLoci = true;
            }
        } else {
            checkLoci = true;
        }
    }

    if (checkLoci) {
        std::map<std::string, std::pair<int, int>> locMap; //loci, seq score, trimpos;
        for (auto & it : mOptions->mLocVars.refLocMap) {
            fpData = it.second.fp.mStr.c_str();
            fpLength = it.second.fp.mStr.length();
            rpData = it.second.rp.mStr.c_str();
            rpLength = it.second.rp.mStr.length();
            bool fpMatched = true;
            bool rpMatched = true;
            int fpMismatches = 0;
            int rpMismatches = 0;
            int numMismaches = 0;
            int numDeletions = 0;
            int numInsertions = 0;
            for (int i = 0; i < it.second.fp.mStr.length(); i++) {
                if (fpData[i] != readSeq[i]) {
                    fpMismatches++;

                    if (fpMismatches > mOptions->mLocVars.locVarOptions.maxMismatchesPSeq) {
                        fpMatched = false;
                        break;
                    }
                }
            }
            if (fpMatched) {
                int rlen = r1->length() - it.second.rp.mStr.length();
                for (int i = 0; i < it.second.rp.mStr.length(); i++) {
                    if (rpData[i] != readSeq[rlen + i]) {
                        rpMismatches++;
                    }
                    if (rpMismatches > mOptions->mLocVars.locVarOptions.maxMismatchesPSeq) {
                        rpMatched = false;
                        break;
                    }
                }

                if (rpMatched) {
                    locMap[it.second.name] = std::make_pair((fpMismatches + rpMismatches), readLength - it.second.fp.mStr.length() - it.second.rp.mStr.length());
                    //should punish if there are mismatches in fp; 
                    //and could not be the best if there are indel in the head of the rp
                }
            }
        }

        if (!locMap.empty()) {
            std::string locName = "";
            if (locMap.size() == 1) {
                locName = locMap.begin()->first;
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
                        break; //warning, what if there are multiple identical values
                    }
                }
            }

            locVarIt = &(mOptions->mLocVars.refLocMap[locName]);
            r1->trimFront(locVarIt->fp.length());
            r1->resize(locMap[locName].second);
            locMap.clear();
            std::map<std::string, std::map<std::string, Genotype>>::iterator itGenotypeMap = tmpAllGenotypeMap.find(locVarIt->name);
            std::map<std::string, Genotype> tmpGenotypeMap;
            Genotype* tmpGenotype = new Genotype();
            if (itGenotypeMap == tmpAllGenotypeMap.end()) {
                if (r1->mSeq.mStr == locVarIt->effectiveSeq.mStr) {//identical seq with ref
                    tmpGenotype->baseLocVar = *locVarIt;
                    tmpGenotype->numReads++;
                    tmpGenotypeMap[r1->mSeq.mStr] = *tmpGenotype;
                    tmpAllGenotypeMap[locVarIt->name] = tmpGenotypeMap;
                    returnedlocus = locVarIt->name + "_true";
                    if (mOptions->mLocVars.locVarOptions.printRes) {
                        //ss << "Genotype: " << tmpGenotype.numReads << " : " << tmpGenotype.baseLocVar.mra.mStr << "\n";
                    }
                } else {
                    bool mraAnalyze = false;
                    std::size_t ffpos = r1->mSeq.mStr.find(locVarIt->ff.mStr);
                    std::size_t rfpos = r1->mSeq.mStr.rfind(locVarIt->rf.mStr);
                    preAnalyze(r1, ffpos, rfpos, mraAnalyze);
                    if (mraAnalyze) {
                        //return the pos and the length of mra after the pos
                        std::pair<size_t, int> mraPosLenReadF = analyzeMRA(r1->mSeq.mStr, locVarIt->repuit.mStr, ffpos, rfpos);

                        ss << "analyzeMRA: " << mraPosLenReadF.first << " : " << mraPosLenReadF.second << " -> " << r1->mSeq.mStr.substr(mraPosLenReadF.first, mraPosLenReadF.second) << "\n";

                        if (mraPosLenReadF.second > 0) {
                            tmpGenotype->baseLocVar = LocVar(r1->mSeq.mStr, mraPosLenReadF, locVarIt->repuit.mStr);
                            tmpGenotype->numReads++;
                            tmpGenotypeMap[r1->mSeq.mStr] = *tmpGenotype;
                            tmpAllGenotypeMap[locVarIt->name] = tmpGenotypeMap;
                            if (mOptions->mLocVars.locVarOptions.printRes) {
                                //ss << "Genotype: " << tmpGenotype.numReads << " : " << tmpGenotype.baseLocVar.mra.mStr << "\n";
                            }
                            returnedlocus = locVarIt->name + "_true";
                        } else {
                            returnedlocus = locVarIt->name + "_failed";
                        }
                    } else {
                        tmpGenotype->baseLocVar = LocVar(r1->mSeq.mStr, ffpos, rfpos, locVarIt->repuit.mStr);
                        tmpGenotype->numReads++;
                        tmpGenotypeMap[r1->mSeq.mStr] = *tmpGenotype;
                        tmpAllGenotypeMap[locVarIt->name] = tmpGenotypeMap;
                        returnedlocus = locVarIt->name + "_true";
                    }

                    if (mOptions->mLocVars.locVarOptions.printRes) {
                        // cCout(ss.str(), 'y');
                    }
                }

            } else {

                std::map<std::string, Genotype>::iterator tmpGenotypeMapIT = itGenotypeMap->second.find(r1->mSeq.mStr);

                if (tmpGenotypeMapIT == itGenotypeMap->second.end()) {

                    if (r1->mSeq.mStr == locVarIt->effectiveSeq.mStr) {//identical seq with ref
                        tmpGenotype->baseLocVar = *locVarIt;
                        tmpGenotype->numReads++;
                        itGenotypeMap->second[r1->mSeq.mStr] = *tmpGenotype;
                        if (mOptions->mLocVars.locVarOptions.printRes) {
                            //ss << "Genotype: " << tmpGenotype.numReads << " : " << tmpGenotype.baseLocVar.mra.mStr << "\n";
                        }
                        returnedlocus = locVarIt->name + "_true";
                    } else {
                        bool mraAnalyze = false;
                        std::size_t ffpos = r1->mSeq.mStr.find(locVarIt->ff.mStr);
                        std::size_t rfpos = r1->mSeq.mStr.rfind(locVarIt->rf.mStr);
                        preAnalyze(r1, ffpos, rfpos, mraAnalyze);
                        if (mraAnalyze) {
                            //return the pos and the length of mra after the pos
                            std::pair<size_t, int> mraPosLenReadF = analyzeMRA(r1->mSeq.mStr, locVarIt->repuit.mStr, ffpos, rfpos);
                            if (mraPosLenReadF.second > 0) {
                                tmpGenotype->baseLocVar = LocVar(r1->mSeq.mStr, mraPosLenReadF, locVarIt->repuit.mStr);
                                tmpGenotype->numReads++;
                                itGenotypeMap->second[r1->mSeq.mStr] = *tmpGenotype;
                                if (mOptions->mLocVars.locVarOptions.printRes) {
                                    //ss << "Genotype: " << tmpGenotype.numReads << " : " << tmpGenotype.baseLocVar.mra.mStr << "\n";
                                }
                                returnedlocus = locVarIt->name + "_true";
                            } else {
                                returnedlocus = locVarIt->name + "_failed";
                            } 
                        } else {
                            tmpGenotype->baseLocVar = LocVar(r1->mSeq.mStr, ffpos, rfpos, locVarIt->repuit.mStr);
                            tmpGenotype->numReads++;
                            itGenotypeMap->second[r1->mSeq.mStr] = *tmpGenotype;
                            returnedlocus = locVarIt->name + "_true";
                        }
                        if (mOptions->mLocVars.locVarOptions.printRes) {
                            //cCout(ss.str(), 'y');
                        }
                    }
                } else {
                    tmpGenotypeMapIT->second.numReads++;
                    returnedlocus = locVarIt->name + "_true";
                }
            }

            if (tmpGenotype) {
                delete tmpGenotype;
                tmpGenotype = NULL;
            }
        }
    }
    ss.str();
    return returnedlocus;
}
