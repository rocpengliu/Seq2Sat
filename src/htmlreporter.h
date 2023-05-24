#ifndef HTML_REPORTER_H
#define HTML_REPORTER_H

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <map>
#include <fstream>
#include <vector>
#include <utility>
#include <set>
#include <algorithm>
#include "options.h"
#include "stats.h"
#include "filterresult.h"
#include "ssrscanner.h"
#include "genotype.h"

using namespace std;

struct Comparator {
    bool operator()(const LocSnp& a, const LocSnp& b) const {
        return a.snpPosSet < b.snpPosSet;
    }
};

class HtmlReporter{
public:
    HtmlReporter(Options* opt);
    ~HtmlReporter();
    void setDupHist(int* dupHist, double* dupMeanGC, double dupRate);
    void setInsertHist(long* insertHist, int insertSizePeak);
    void report(std::vector<std::map<std::string, std::vector<std::pair<std::string, Genotype>>>> &sortedAllGenotypeMapVec,
            std::map<std::string, std::map<std::string, LocSnp>> & allSnpsMap,
            FilterResult* result, Stats* preStats1, Stats* postStats1, Stats* preStats2 = NULL, Stats* postStats2 = NULL);

    static void outputRow(ofstream& ofs, string key, long value);
    static void outputRow(ofstream& ofs, string key, string value);
    static void outputRow(ofstream& ofs, std::string & marker, std::vector<std::pair<std::string, Genotype>> & outGenotype, Options*& mOptions);
    static void outputRow(ofstream& ofs, std::string & marker, std::map<std::string, LocSnp> & snpsMap, int & totReads, std::set<int> & refSet);
    static void outputRow(ofstream& ofs, std::map<int, SimGeno> & snpGenoMap, std::set<int> & refSet, int & totReads);
    static string formatNumber(long number);
    static string getPercents(long numerator, long denominator);
    static std::string highligher(std::string & str, std::map<int, std::string> & snpsMap);
    static std::string highligher(std::string & str, std::set<int> & snpsSet);
    static std::string highligher(LocSnp & locSnp, bool ref, std::set<int> & refSet);
    
private:
    const string getCurrentSystemTime();
    void printHeader(ofstream& ofs);
    void printCSS(ofstream& ofs);
    void printJS(ofstream& ofs);
    void printFooter(ofstream& ofs);
    void reportDuplication(ofstream& ofs);
    void reportInsertSize(ofstream& ofs, int isizeLimit);
    void printSummary(ofstream& ofs, FilterResult* result, Stats* preStats1, Stats* postStats1, Stats* preStats2, Stats* postStats2);
    void reportAllGenotype(ofstream& ofs, std::vector<std::map<std::string, std::vector<std::pair<std::string, Genotype>>>> & sortedAllGenotypeMapVec);
    void reportAllSnps(ofstream& ofs, std::map<std::string, std::map<std::string, LocSnp>> & allSnpsMap);
    void reportEachGenotype(ofstream& ofs, std::string marker, 
                            std::vector<int> & x_vec, std::map< std::string, std::vector<int>> & stackMap, 
                            std::map< std::string, std::vector<std::string>> & stackYlabMap, std::vector<double> & bar_width_vec, 
                            std::vector<std::pair<std::string, Genotype>> & outGenotype,
                            std::vector<int> & xmra_vec, std::map< std::string, std::vector<int>> & stackMraMap, 
                            std::map< std::string, std::vector<std::string>> & stackYLabMMap, 
                            std::vector<double> & barmra_width_vec, 
                            std::vector<std::pair<std::string, Genotype>> & outGenotypeMra);
    void reportEachSnpGenotype(ofstream& ofs, std::string marker, std::map<std::string, LocSnp> & snpsMap,
            std::vector<std::string> & x_vec, std::vector<int> & y_vec, std::vector<double> & bar_width_vec, 
            std::vector<std::string> & x_vec_2, std::vector<int> & y_vec_2, std::vector<double> & bar_width_vec_2,
            int & totReads, bool & twoTrace);
    void reportSnpAlignmentTable(ofstream& ofs, std::string marker, std::map<std::string, LocSnp> & snpsMap, int & totReads);
    void reportSnpTablePlot(ofstream& ofs, std::string marker, std::string & divName, int & totReads);
    void reportSex(ofstream & ofs);
    
private:
    Options* mOptions;
    int* mDupHist;
    double* mDupMeanGC;
    double mDupRate;
    long* mInsertHist;
    int mInsertSizePeak;
};


#endif