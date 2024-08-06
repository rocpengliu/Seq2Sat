#ifndef HTMLREPORTERALL_H
#define HTMLREPORTERALL_H

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <map>
#include <fstream>
#include <vector>
#include <utility>
#include <set>
#include "options.h"
#include "stats.h"
#include "filterresult.h"
#include "ssrscanner.h"
#include "genotype.h"

using namespace std;

class HtmlReporterAll {
public:
    HtmlReporterAll(Options* opt);
    HtmlReporterAll(const HtmlReporterAll& orig);
    virtual ~HtmlReporterAll();
    void report();
    static void outputRow(ofstream& ofs, string key, long value);
    static void outputRow(ofstream& ofs, string key, string value);
    static void outputRow(ofstream& ofs, std::string & marker, std::vector<std::pair<std::string, Genotype>> & outGenotype);
    static void outputRow(ofstream& ofs, std::string & marker, std::map<std::string, LocSnp> & snpsMap);
    static string formatNumber(long number);
    static string getPercents(long numerator, long denominator);
    static std::string highligher(std::string & str, std::map<int, std::string> & snpsMap);
    static std::string highligherSet(std::string & str, std::set<int> & snpsSet);
    static std::string highligher(LocSnp & locSnp, bool ref = false);
    
private:
    void reportEachSample(ofstream& ofs, Sample & sample);
    void reportEachGenotype(ofstream& ofs, std::string sample, std::string marker,
            std::vector<int> & x_vec, std::map< std::string, std::vector<int>> &stackMap,
            std::map< std::string, std::vector<std::string>> &stackYlabMap, std::vector<double> & bar_width_vec,
            std::vector<std::pair<std::string, Genotype>> &outGenotype,
            std::vector<int> & xmra_vec, std::map< std::string, std::vector<int>> &stackMraMap,
            std::map< std::string, std::vector<std::string>> &stackYLabMMap,
            std::vector<double> & barmra_width_vec,
            std::vector<std::pair<std::string, Genotype>> &outGenotypeMra);
    const string getCurrentSystemTime();
    void printHeader(ofstream& ofs);
    void printCSS(ofstream& ofs);
    void printJS(ofstream& ofs);
    void printFooter(ofstream& ofs);
    
private:
    Options* mOptions;
};

#endif /* HTMLREPORTERALL_H */

