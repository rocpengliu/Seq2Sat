#ifndef OPTIONS_H
#define OPTIONS_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <map>

#include "edlib.h"
#include "genotype.h"

using namespace std;

#define UMI_LOC_NONE 0
#define UMI_LOC_INDEX1 1
#define UMI_LOC_INDEX2 2
#define UMI_LOC_READ1 3
#define UMI_LOC_READ2 4
#define UMI_LOC_PER_INDEX 5
#define UMI_LOC_PER_READ 6

class DuplicationOptions {
public:
    DuplicationOptions() {
        enabled = true;
        keylen = 12;
        histSize = 32;
    }
public:
    bool enabled;
    int keylen;
    int histSize;
};


class LowComplexityFilterOptions {
public:
    LowComplexityFilterOptions() {
        enabled = false;
        threshold = 0.3;
    }
public:
    bool enabled;
    double threshold;
};

class PolyGTrimmerOptions {
public:
    PolyGTrimmerOptions() {
        enabled = false;
        minLen = 10;
    }
public:
    bool enabled;
    int minLen;
};

class PolyXTrimmerOptions {
public:
    PolyXTrimmerOptions() {
        enabled = false;
        minLen = 10;
    }
public:
    bool enabled;
    int minLen;
};

class UMIOptions {
public:
    UMIOptions() {
        enabled = false;
        location = UMI_LOC_NONE;
        length = 0;
        skip = 0;
    }
public:
    bool enabled;
    int location;
    int length;
    int skip;
    string prefix;
    string separator;
};

class CorrectionOptions {
public:
    CorrectionOptions() {
        enabled = true;
    }
public:
    bool enabled;
};

class QualityCutOptions {
public:
    QualityCutOptions() {
        enabledFront = false;
        enabledTail = false;
        enabledRight = false;
        windowSizeShared = 4;
        qualityShared = 20;
        windowSizeFront = windowSizeShared;
        qualityFront = qualityShared;
        windowSizeTail = windowSizeShared;
        qualityTail = qualityShared;
        windowSizeRight = windowSizeShared;
        qualityRight = qualityShared;
    }
public:
    // enable 5' cutting by quality
    bool enabledFront;
    // enable 3' cutting by quality
    bool enabledTail;
    // enable agressive cutting mode
    bool enabledRight;
    // the sliding window size
    int windowSizeShared;
    // the mean quality requirement
    int qualityShared;
    // the sliding window size for cutting by quality in 5'
    int windowSizeFront;
    // the mean quality requirement for cutting by quality in 5'
    int qualityFront;
    // the sliding window size for cutting by quality in 3'
    int windowSizeTail;
    // the mean quality requirement for cutting by quality in 3'
    int qualityTail;
    // the sliding window size for cutting by quality in aggressive mode
    int windowSizeRight;
    // the mean quality requirement for cutting by quality in aggressive mode
    int qualityRight;
};

class AdapterOptions {
public:
    AdapterOptions() {
        enabled = true;
        hasSeqR1 = false;
        hasSeqR2 = false;
        detectAdapterForPE = false;
        hasFasta = false;
    }

public:
    bool enabled;
    string sequence;
    string sequenceR2;
    string detectedAdapter1;
    string detectedAdapter2;
    vector<string> seqsInFasta;
    string fastaFile;
    bool hasSeqR1;
    bool hasSeqR2;
    bool hasFasta;
    bool detectAdapterForPE;
};

class TrimmingOptions {
public:
    TrimmingOptions() {
        front1 = 0;
        tail1 = 0;
        front2 = 0;
        tail2 = 0;
        maxLen1 = 0;
        maxLen2 = 0;
    }
public:
    // trimming first cycles for read1
    int front1;
    // trimming last cycles for read1
    int tail1;
    // trimming first cycles for read2
    int front2;
    // trimming last cycles for read2
    int tail2;
    // max length of read1
    int maxLen1;
    // max length of read2
    int maxLen2;
};

class QualityFilteringOptions {
public:
    QualityFilteringOptions() {
        enabled = true;
        // '0' = Q15
        qualifiedQual = '0';
        unqualifiedPercentLimit = 40;
        nBaseLimit = 5;
    }
public:
    // quality filter enabled
    bool enabled;
    // if a base's quality phred score < qualifiedPhred, then it's considered as a low_qual_base
    char qualifiedQual;
    // if low_qual_base_num > lowQualLimit, then discard this read
    int unqualifiedPercentLimit;
    // if n_base_number > nBaseLimit, then discard this read
    int nBaseLimit;
    // if average qual score < avgQualReq, then discard this read
    int avgQualReq;
};

class ReadLengthFilteringOptions {
public:
    ReadLengthFilteringOptions() {
        enabled = false;
        requiredLength = 50;
        maxLength = 0;
    }
public:
    // length filter enabled
    bool enabled;
    // if read_length < requiredLength, then this read is discard
    int requiredLength;
    // length limit, 0 for no limitation
    int maxLength;
};

class EdOptions{
public:
    EdOptions(){
    mode = "HW";
    silent = false;
    printRes = true;
    format = "NICE";
    findAlignment = false;
    findStartLocation = false;
    modeCode = EDLIB_MODE_HW;
    alignTask = EDLIB_TASK_DISTANCE;
    maxDeletionRPrimer = 4;
    maxInsertionRPrimer = 4;
    maxMismatchesRPrimer = 4;
    };
    
public:
    std::string mode;
    bool findAlignment;
    bool findStartLocation;
    bool silent;
    bool printRes;
    std::string format;
    EdlibAlignMode modeCode;
    EdlibAlignTask alignTask;
    int maxDeletionRPrimer;
    //T: ATCTCTTCATATATATCCTC-T--CT (156 - 178)//read
    //   ||||||||| ||||| |||| |  ||
    //Q: ATCTCTTCAAATATAGCCTCATTGCT (0 - 25)//primer
    int maxInsertionRPrimer;
    //T: ACTGGTAAATACATC-TTCCCCACT (143 - 166)read
    //   ||||||||||||||| ||||| |||
    //Q: ACTGGTAAATACATCATTCCC-ACT (0 - 23)primer
    int maxMismatchesRPrimer;
    //T: TGACTACAGGCTTGCACCG (115 - 133)read
    //   |||| ||||||||||| |
    //Q: GGACTCCAGGCTTGCACTG (0 - 18)primer
};

class LocVarOptions{
public:
    LocVarOptions(){
        maxMismatchesPSeq = 2;
        minSeqs = 5;
        minSeqsPer = 5;
        minWarningSeqs = 50;
        minMatchesFR = 6;
        maxMismatchesPer4FR = 0.3;
        maxScore = -1;
        maxScorePrimer = -1;
        numBestSeqs = 0;
        coreRep = 1;
        minNSSRUnit = 3;
        hlRatio1 = 0.4;
        hlRatio2 = 0.2;
        varRatio = 1.5;
        heterRatio = 0.3;
    };
public:
    int maxMismatchesPSeq;
    int minSeqs;
    int minWarningSeqs;
    int minSeqsPer;//against largest peak; 5/%
    int minMatchesFR;//minimum matches bps for both flanking regions
    double maxMismatchesPer4FR;//percentage for minimum matches bps for both flanking regions;
    int maxScore;
    int maxScorePrimer;
    int numBestSeqs;
    int coreRep;
    int minNSSRUnit;
    double hlRatio1; //if two locus size is 1 ssr unit;
    double hlRatio2; //if two locus size is 2 ssr unit;
    double varRatio; //if there are variations in mra, flanking regions, the n of one allele against the other, the ideal should be 1, but default one is 1.5
    double heterRatio;//how to determine variations in mra, flanking regions are true, and if two max peak size are same, if it > heterratio, the variation based genotype is true
};

class LocSnpOptions{
public:
    LocSnpOptions(){
        maxMismatchesPSeq = 2;
        maxScorePrimer = -1;
        minSeqs = 5;
        minSeqsPer = 0.1;
        hmPerH = 0.9;
        hmPerL = 0.85;
        hmPer = 0.0;
        htJetter = 0.25;
        minReads4Filter = 50;
        maxRows4Align = 6;
    };
public:
    uint32 maxMismatchesPSeq;
    uint32 maxScorePrimer;
    uint32 minSeqs;
    double minSeqsPer;//against largest peak; 10/%
    int minReads4Filter;
    double hmPerL, hmPerH, htJetter;
    double hmPer; //the actual one, it is either hmPerL or hmPerH, depending on how many true snps there.
    int maxRows4Align;
};

class MultiLocVars{
public:
    MultiLocVars(){
        refLocMap.clear();
    };
    
public:
    std::map<std::string, LocVar> refLocMap;//change to pointer to reduce memory usage// sample -> locvar
    LocVarOptions locVarOptions;
};

class MultiLocSnps{
public:
    MultiLocSnps(){
        refLocMap.clear();
    };
    
public:
    std::map<std::string, LocSnp2> refLocMap;//change to pointer to reduce memory usage
    LocSnpOptions mLocSnpOptions;
};

enum varType {
    ssr,
    snp
};

class Sample{
public:
    Sample(){
        prefix = "";
        path = "";
        in1 = "";
        in2 = "";
        feature = "";
        cmd = "";
        numRawReads = 0;
        numCleanReads = 0;
        sortedAllGenotypeMapVec.clear();
    }
public:
    std::string prefix;
    std::string path;
    std::string in1;
    std::string in2;
    std::string feature;
    time_t startTime;
    time_t endTime;
    std::string cmd;
    int numRawReads;
    int numCleanReads;
    std::vector<std::map <std::string, std::vector<std::pair<std::string, Genotype>>>> sortedAllGenotypeMapVec;
    //void reset(){sortedAllGenotypeMapVec.clear();};
};

class Options{
public:
    Options();
    void init();
    bool isPaired();
    bool validate();
    bool adapterCuttingEnabled();
    bool polyXTrimmingEnabled();
    string getAdapter1();
    string getAdapter2();
    bool shallDetectAdapter(bool isR2 = false);
    void loadFastaAdapters();
    void readLocFile();
    void readSexLoc();
    void parseSampleTable();
public:
    // file name of read1 input
    string in1;
    // file name of read2 input
    string in2;
    // file name of read1 output
    string out1;
    // file name of read2 output
    string out2;
    //output failed reads
    bool outFR;
    //file name of failed reads;
    string outFRFile;
    // genome FASTA file
    string genomeFile;
    // kmer FASTA file
    string kmerFile;
    // kmer FASTA file
    string kmerCollectionFile;
    // json file
    string jsonFile;
    // html file
    string htmlFile;
    // html report title
    string reportTitle;
    // compression level
    int compression;
    // the input file is using phred64 quality scoring
    bool phred64;
    // do not rewrite existing files
    bool dontOverwrite;
    // read STDIN
    bool inputFromSTDIN;
    // write STDOUT
    bool outputToSTDOUT;
    // the input R1 file is interleaved
    bool interleavedInput;
    // only process first N reads
    int readsToProcess;
    // worker thread number
    int thread;
    // trimming options
    TrimmingOptions trim;
    // quality filtering options
    QualityFilteringOptions qualfilter;
    // length filtering options
    ReadLengthFilteringOptions lengthFilter;
    // adapter options
    AdapterOptions adapter;
    // options for quality cutting
    QualityCutOptions qualityCut;
    // options for base correction
    CorrectionOptions correction;
    // options for UMI
    UMIOptions umi;
    // 3' end polyG trimming, default for Illumina NextSeq/NovaSeq
    PolyGTrimmerOptions polyGTrim;
    // 3' end polyX trimming
    PolyXTrimmerOptions polyXTrim;
    int seqLen1;
    int seqLen2;
    // low complexity filtering
    LowComplexityFilterOptions complexityFilter;
    // options for duplication profiling
    DuplicationOptions duplicate;
    // options for duplication profiling
    int insertSizeMax;
    // overlap analysis threshold
    int overlapRequire;
    int overlapDiffLimit;
    int overlapDiffPercentLimit;
    // output debug information
    bool verbose;
    // the length of KMER, default is 25
    bool debug;
    bool revCom;
    varType mVarType;
    std::string var;
    std::string sampleTable;
    
    int minAmpliconEffectiveLen;//length without primers for both sex and nonsex loci
    MultiLocVars mLocVars;
    MultiLocSnps mLocSnps;
    EdOptions mEdOptions;
    string prefix;
    std::string locFile;//loci name, forward primer, reverse primer (RC), forward flanking region, reverse flanking region, ssr repeat unit, number of repeat, mra; must be separated by tab;
    std::string sexFile;//sex id name, forward primer, reverse primer (reverse complement), x ref, y ref; must separated by tab;
    
    //merge overlapped PE read;
    bool mergerOverlappedPE;
    std::vector<Sample> samples;
    Sex mSex;
    std::map<std::string, Sex> sexMap;

    bool noErrorPlot;// not to show error plot in html report for sex and snp id;
    bool noSnpPlot; // not to show snp plot in html report;
    bool noPlot;// not to plot any figures;
    bool nanopore_default;
};

#endif
