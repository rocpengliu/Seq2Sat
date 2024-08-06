#include "threadconfig.h"

ThreadConfig::ThreadConfig(Options* opt, int threadId, bool paired){
    mOptions = opt;
    mThreadId = threadId;
    mPreStats1 = new Stats(mOptions, false);
    mPostStats1 = new Stats(mOptions, false);
    if(paired){
        mPreStats2 = new Stats(mOptions, true);
        mPostStats2 = new Stats(mOptions, true);
    } else {
        mPreStats2 = NULL;
        mPostStats2 = NULL;
    }
    mWriter1 = NULL;
    mWriter2 = NULL;

    mFilterResult = new FilterResult(opt, paired);
    mCanBeStopped = false;
    mSsrScanner = new SsrScanner(mOptions);
    mSnpScanner = new SnpScanner(mOptions);
    mSexScanner = new SexScanner(mOptions);
}

ThreadConfig::~ThreadConfig() {
    cleanup();
    if(mSnpScanner){
        delete mSnpScanner;
        mSnpScanner = NULL;
    }
    if(mSsrScanner){
        delete mSsrScanner;
        mSsrScanner = NULL;
    }
    if(mSexScanner){
        delete mSexScanner;
        mSexScanner = nullptr;
    }
    if(mPreStats1){
        delete mPreStats1;
        mPreStats1 = NULL;
    }
    if(mPostStats1){
        delete mPostStats1;
        mPostStats1 = NULL;
    }
    if(mPreStats2){
        delete mPreStats2;
        mPreStats2 = NULL;
    }
    if(mPostStats2){
        delete mPostStats2;
        mPostStats2 = NULL;
    }
    if(mFilterResult){
        delete mFilterResult;
        mFilterResult = NULL;
    }
}

void ThreadConfig::cleanup() {
    deleteWriter();
}

void ThreadConfig::deleteWriter() {
    if(mWriter1 != NULL) {
        delete mWriter1;
        mWriter1 = NULL;
    }
    if(mWriter2 != NULL) {
        delete mWriter2;
        mWriter2 = NULL;
    }
}

void ThreadConfig::initWriter(string filename1) {
    deleteWriter();
    mWriter1 = new Writer(filename1, mOptions->compression);
}

void ThreadConfig::initWriter(string filename1, string filename2) {
    deleteWriter();
    mWriter1 = new Writer(filename1, mOptions->compression);
    mWriter2 = new Writer(filename2, mOptions->compression);
}

void ThreadConfig::initWriter(ofstream* stream) {
    deleteWriter();
    mWriter1 = new Writer(stream);
}

void ThreadConfig::initWriter(ofstream* stream1, ofstream* stream2) {
    deleteWriter();
    mWriter1 = new Writer(stream1);
    mWriter2 = new Writer(stream2);
}

void ThreadConfig::initWriter(gzFile gzfile) {
    deleteWriter();
    mWriter1 = new Writer(gzfile);
}

void ThreadConfig::initWriter(gzFile gzfile1, gzFile gzfile2) {
    deleteWriter();
    mWriter1 = new Writer(gzfile1);
    mWriter2 = new Writer(gzfile2);
}

void ThreadConfig::addFilterResult(int result, int readNum) {
    mFilterResult->addFilterResult(result, readNum);
}

void ThreadConfig::addMergedPairs(int pairs) {
    mFilterResult->addMergedPairs(pairs);
}

void ThreadConfig::markProcessed(long readNum) {
}

bool ThreadConfig::canBeStopped() {
    return mCanBeStopped;
}