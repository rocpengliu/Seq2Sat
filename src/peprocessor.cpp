#include "peprocessor.h"

PairEndProcessor::PairEndProcessor(Options* opt){
    mOptions = opt;
    mProduceFinished = false;
    mFinishedThreads = 0;
    mFilter = new Filter(opt);
    mZipFile1 = NULL;
    mZipFile2 = NULL;
    mUmiProcessor = new UmiProcessor(opt);

    int isizeBufLen = mOptions->insertSizeMax + 1;
    mInsertSizeHist = new long[isizeBufLen];
    memset(mInsertSizeHist, 0, sizeof(long)*isizeBufLen);
    mLeftWriter =  NULL;
    mRightWriter = NULL;
    mFailedWriter = NULL;

    mDuplicate = NULL;
    if(mOptions->duplicate.enabled) {
        mDuplicate = new Duplicate(mOptions);
    }
}

PairEndProcessor::~PairEndProcessor() {
    delete[] mInsertSizeHist;
    if(mFilter){
        delete mFilter;
        mFilter = NULL;
    }
    if(mDuplicate) {
        delete mDuplicate;
        mDuplicate = NULL;
    }

    if (mUmiProcessor) {
        delete mUmiProcessor;
        mUmiProcessor = NULL;
    }
}

void PairEndProcessor::initOutput() {
//    if(mOptions->out1.empty())
//        return;
    
    if (!mOptions->out1.empty()) {
        mLeftWriter = new WriterThread(mOptions, mOptions->out1);
    }
    
    if(!mOptions->out2.empty()){
        mRightWriter = new WriterThread(mOptions, mOptions->out2);
    }
    
    if(!mOptions->outFRFile.empty()){
        mFailedWriter = new WriterThread(mOptions, mOptions->outFRFile);
    }
}

void PairEndProcessor::closeOutput() {
    if(mLeftWriter) {
        delete mLeftWriter;
        mLeftWriter = NULL;
    }
    if(mRightWriter) {
        delete mRightWriter;
        mRightWriter = NULL;
    }
    if(mFailedWriter){
        delete mFailedWriter;
        mFailedWriter = NULL;
    }
}

void PairEndProcessor::initConfig(ThreadConfig* config) {
    if(mOptions->out1.empty())
        return;
}

bool PairEndProcessor::process(){
    initOutput();
    initPackRepository();
    std::thread producer(std::bind(&PairEndProcessor::producerTask, this));
    //TODO: get the correct cycles
    int cycle = 151;
    ThreadConfig** configs = new ThreadConfig*[mOptions->thread];
    for(int t=0; t<mOptions->thread; t++){
        configs[t] = new ThreadConfig(mOptions, t, true);
        initConfig(configs[t]);
    }

    std::thread** threads = new thread*[mOptions->thread];
    for(int t=0; t<mOptions->thread; t++){
        threads[t] = new std::thread(std::bind(&PairEndProcessor::consumerTask, this, configs[t]));
    }

    std::thread* leftWriterThread = NULL;
    std::thread* rightWriterThread = NULL;
    std::thread* failedWriterThread = NULL;
    if(mLeftWriter)
        leftWriterThread = new std::thread(std::bind(&PairEndProcessor::writeTask, this, mLeftWriter));
    if(mRightWriter)
        rightWriterThread = new std::thread(std::bind(&PairEndProcessor::writeTask, this, mRightWriter));
    if(mFailedWriter){
        failedWriterThread = new std::thread(std::bind(&PairEndProcessor::writeTask, this, mFailedWriter));
    }
    
    producer.join();
    for(int t=0; t<mOptions->thread; t++){
        threads[t]->join();
    }

    if(leftWriterThread)
        leftWriterThread->join();
    if(rightWriterThread)
        rightWriterThread->join();
    if(failedWriterThread)
        failedWriterThread->join();

    destroyPackRepository();
    
    if(mOptions->verbose)
        loginfo("start to generate reports\n");

    // merge stats and filter results
    vector<Stats*> preStats1;
    vector<Stats*> postStats1;
    vector<Stats*> preStats2;
    vector<Stats*> postStats2;
    vector<FilterResult*> filterResults;
    std::vector<std::map<std::string, std::map<std::string, Genotype>>> totalGenotypeSsrMapVec;
    std::vector<std::map<std::string, std::map<std::string, uint32>>> totalSnpSeqMapVec;
    std::vector<std::map<std::string, std::map<std::string, int>>> totalSexLocVec;
    totalSexLocVec.reserve(mOptions->thread);
    if (mOptions->mVarType == ssr) {
        totalGenotypeSsrMapVec.reserve(mOptions->thread);
    } else if (mOptions->mVarType == snp) {
        totalSnpSeqMapVec.reserve(mOptions->thread);
    }
    for(int t=0; t<mOptions->thread; t++){
        preStats1.push_back(configs[t]->getPreStats1());
        postStats1.push_back(configs[t]->getPostStats1());
        preStats2.push_back(configs[t]->getPreStats2());
        postStats2.push_back(configs[t]->getPostStats2());
        filterResults.push_back(configs[t]->getFilterResult());
        if(!mOptions->mSex.sexMarker.empty()){
            totalSexLocVec.emplace_back(configs[t]->getSexScanner()->getSexLoc());
        }
        if (mOptions->mVarType == ssr) {
            totalGenotypeSsrMapVec.push_back(configs[t]->getSsrScanner()->getGenotypeMap());
        } else if(mOptions->mVarType == snp){
            totalSnpSeqMapVec.push_back(configs[t]->getSnpScanner()->getSubSeqsMap());
        }
    }
    Stats* finalPreStats1 = Stats::merge(preStats1);
    Stats* finalPostStats1 = Stats::merge(postStats1);
    Stats* finalPreStats2 = Stats::merge(preStats2);
    Stats* finalPostStats2 = Stats::merge(postStats2);
    FilterResult* finalFilterResult = FilterResult::merge(filterResults);

    std::map<std::string, std::map< std::string, Genotype>> allGenotypeMap;//marker, seq, geno
    std::vector<std::map <std::string, std::vector<std::pair<std::string, Genotype>>>> sortedAllGenotypeMapVec;

    if (!mOptions->mSex.sexMarker.empty()) {
        SexScanner::merge(totalSexLocVec, mOptions);
        SexScanner::report(mOptions);
    }
    
    if (mOptions->mVarType == ssr) {
         allGenotypeMap = SsrScanner::merge(totalGenotypeSsrMapVec);
         sortedAllGenotypeMapVec = SsrScanner::report(mOptions, allGenotypeMap);
         if(!mOptions->samples.empty()){
             for(auto & it : mOptions->samples){
                 if(it.prefix == mOptions->prefix){
                     //it.sortedAllGenotypeMapVec = sortedAllGenotypeMapVec;
                     break;
                 }
             }
         }
    } else if (mOptions->mVarType == snp) {
        SnpScanner::merge2(mOptions, totalSnpSeqMapVec);
    }

    int* dupHist = NULL;
    double* dupMeanTlen = NULL;
    double* dupMeanGC = NULL;
    double dupRate = 0.0;
    if(mOptions->duplicate.enabled) {
        dupHist = new int[mOptions->duplicate.histSize];
        memset(dupHist, 0, sizeof(int) * mOptions->duplicate.histSize);
        dupMeanGC = new double[mOptions->duplicate.histSize];
        memset(dupMeanGC, 0, sizeof(double) * mOptions->duplicate.histSize);
        dupRate = mDuplicate->statAll(dupHist, dupMeanGC, mOptions->duplicate.histSize);
        std::cerr << endl;
        std::cerr << "Duplication rate: " << dupRate * 100.0 << "%" << endl;
    }

    // insert size distribution
    int peakInsertSize = getPeakInsertSize();
    std::cerr << endl;
    std::cerr << "Insert size peak (evaluated by paired-end reads): " << peakInsertSize << endl;

    // make JSON report
    
    JsonReporter jr(mOptions);
    jr.setDupHist(dupHist, dupMeanGC, dupRate);
    jr.setInsertHist(mInsertSizeHist, peakInsertSize);
    jr.report(sortedAllGenotypeMapVec, finalFilterResult, finalPreStats1, finalPostStats1, finalPreStats2, finalPostStats2);
    std::cerr << "Finished Json report" << endl;
    // make HTML report
    HtmlReporter hr(mOptions);
    hr.setDupHist(dupHist, dupMeanGC, dupRate);
    hr.setInsertHist(mInsertSizeHist, peakInsertSize);
    hr.report(sortedAllGenotypeMapVec, finalFilterResult, finalPreStats1, finalPostStats1, finalPreStats2, finalPostStats2);
    std::cerr << "Finished Html report" << endl;

    // clean up
    for(int t=0; t<mOptions->thread; t++){
        delete threads[t];
        threads[t] = NULL;
        delete configs[t];
        configs[t] = NULL;
    }
    delete finalPreStats1;
    delete finalPostStats1;
    delete finalPreStats2;
    delete finalPostStats2;
    delete finalFilterResult;
    if(mOptions->duplicate.enabled) {
        delete[] dupHist;
        delete[] dupMeanGC;
    }
    delete[] threads;
    delete[] configs;
    if(leftWriterThread)
        delete leftWriterThread;
    if(rightWriterThread)
        delete rightWriterThread;
    if(failedWriterThread)
        delete failedWriterThread;
    closeOutput();
    return true;
}

int PairEndProcessor::getPeakInsertSize() {
    int peak = 0;
    long maxCount = -1;
    for(int i=0; i<mOptions->insertSizeMax; i++) {
        if(mInsertSizeHist[i] > maxCount) {
            peak = i;
            maxCount = mInsertSizeHist[i];
        }
    }
    return peak;
}

bool PairEndProcessor::processPairEnd(ReadPairPack* pack, ThreadConfig* config){
    string outstr1;
    string outstr2;
    string singleOutput;
    string locus = "";
    string failedOutput;
    int readPassed = 0;
    for(int p=0;p<pack->count;p++){
        ReadPair* pair = pack->data[p];
        Read* or1 = pair->mLeft;
        Read* or2 = pair->mRight;

        int lowQualNum1 = 0;
        int nBaseNum1 = 0;
        int lowQualNum2 = 0;
        int nBaseNum2 = 0;

        // stats the original read before trimming
        config->getPreStats1()->statRead(or1);
        config->getPreStats2()->statRead(or2);

        // handling the duplication profiling
        if(mDuplicate)
            mDuplicate->statPair(or1, or2);

        // umi processing
        if(mOptions->umi.enabled)
            mUmiProcessor->process(or1, or2);

        // trim in head and tail, and apply quality cut in sliding window
        int frontTrimmed1 = 0;
        int frontTrimmed2 = 0;
//        Read* r1 = mFilter->trimAndCut(or1, mOptions->trim.front1, mOptions->trim.tail1, frontTrimmed1);
//        Read* r2 = mFilter->trimAndCut(or2, mOptions->trim.front2, mOptions->trim.tail2, frontTrimmed2);
        Read* r1 = mFilter->trimAndCut(or1, 0, mOptions->trim.tail1, frontTrimmed1);
        //if(r1 = NULL) delete or1; or1 = NULL;
        Read* r2 = mFilter->trimAndCut(or2, 0, mOptions->trim.tail2, frontTrimmed2);
        //if(r2 = NULL) delete or2; or2 = NULL;
        if(r1 != NULL && r2!=NULL) {
            if(mOptions->polyGTrim.enabled)
                PolyX::trimPolyG(r1, r2, config->getFilterResult(), mOptions->polyGTrim.minLen);
        }
        bool isizeEvaluated = false;
        if(r1 != NULL && r2!=NULL && (mOptions->adapter.enabled || mOptions->correction.enabled)){
            OverlapResult ov = OverlapAnalysis::analyze(r1, r2, mOptions->overlapDiffLimit, mOptions->overlapRequire, mOptions->overlapDiffPercentLimit/100.0);
            // we only use thread 0 to evaluae ISIZE
            if(config->getThreadId() == 0) {
                statInsertSize(r1, r2, ov, frontTrimmed1, frontTrimmed2);
                isizeEvaluated = true;
            }
            if(mOptions->correction.enabled) {
                BaseCorrector::correctByOverlapAnalysis(r1, r2, config->getFilterResult(), ov);
            }
            if(mOptions->adapter.enabled) {
                bool trimmed = AdapterTrimmer::trimByOverlapAnalysis(r1, r2, config->getFilterResult(), ov, frontTrimmed1, frontTrimmed2);
                bool trimmed1 = trimmed;
                bool trimmed2 = trimmed;
                if(!trimmed){
                    if(mOptions->adapter.hasSeqR1)
                        trimmed1 = AdapterTrimmer::trimBySequence(r1, config->getFilterResult(), mOptions->adapter.sequence, false);
                    if(mOptions->adapter.hasSeqR2)
                        trimmed2 = AdapterTrimmer::trimBySequence(r2, config->getFilterResult(), mOptions->adapter.sequenceR2, true);
                }
                if(mOptions->adapter.hasFasta) {
                    AdapterTrimmer::trimByMultiSequences(r1, config->getFilterResult(), mOptions->adapter.seqsInFasta, false, !trimmed1);
                    AdapterTrimmer::trimByMultiSequences(r2, config->getFilterResult(), mOptions->adapter.seqsInFasta, true, !trimmed2);
                }
            }
        }

        if(config->getThreadId() == 0 && !isizeEvaluated && r1 != NULL && r2!=NULL) {
            OverlapResult ov = OverlapAnalysis::analyze(r1, r2, mOptions->overlapDiffLimit, mOptions->overlapRequire, mOptions->overlapDiffPercentLimit/100.0);
            statInsertSize(r1, r2, ov, frontTrimmed1, frontTrimmed2);
            isizeEvaluated = true;
        }

        if(r1 != NULL && r2!=NULL) {
            if(mOptions->polyXTrim.enabled)
                PolyX::trimPolyX(r1, r2, config->getFilterResult(), mOptions->polyXTrim.minLen);
        }

        if(r1 != NULL && r2!=NULL) {
            if( mOptions->trim.maxLen1 > 0 && mOptions->trim.maxLen1 < r1->length())
                r1->resize(mOptions->trim.maxLen1);
            if( mOptions->trim.maxLen2 > 0 && mOptions->trim.maxLen2 < r2->length())
                r2->resize(mOptions->trim.maxLen2);
        }


        int result1 = mFilter->passFilter(r1);
        int result2 = mFilter->passFilter(r2);

        config->addFilterResult(max(result1, result2), 2);

        if( r1 != NULL &&  result1 == PASS_FILTER && r2 != NULL && result2 == PASS_FILTER ) {
            if(mOptions->mergerOverlappedPE){
                OverlapResult ov = OverlapAnalysis::analyze(r1, r2, mOptions->overlapDiffLimit, mOptions->overlapRequire, mOptions->overlapDiffPercentLimit / 100.0);
                if(ov.overlapped){
                    Read* merged = OverlapAnalysis::merge(r1, r2, ov);
                    int result = mFilter->passFilter(merged);
                    if (result == PASS_FILTER) {
                        locus.clear();
                        bool goGeno = false;
                        if(!mOptions->mSex.sexMarker.empty()){
                            auto rep = config->getSexScanner()->sexScan(merged);
                            if (rep.first) {
                                locus = rep.second;
                                goGeno = false;
                            } else {
                                Read* revReads = merged->reverseComplement();
                                rep = config->getSexScanner()->sexScan(revReads);
                                if(rep.first){
                                    locus = rep.second;
                                    goGeno = false;
                                } else {
                                    goGeno = true;
                                }
                                delete revReads;
                            }
                        } else {
                            goGeno = true;
                        }

                        if (goGeno) {
                            if (mOptions->mVarType == ssr) {
                                locus = config->getSsrScanner()->scanVar(merged);
                                size_t found = locus.find("_failed");
                                if (found != std::string::npos && found == (locus.length() - 7)) {
                                    Read* revReads = merged->reverseComplement();
                                    locus = config->getSsrScanner()->scanVar(revReads);
                                    delete revReads;
                                }
                            } else {
                                locus = config->getSnpScanner()->scanVar(merged);
                                if(locus.empty()){
                                    Read* revReads = merged->reverseComplement();
                                    locus = config->getSnpScanner()->scanVar(revReads);
                                    delete revReads;
                                }
                            }
                        }

                        if (!locus.empty()) {
                            failedOutput += merged->toStringWithTag(locus);//bug here
                        } else {
                            outstr1 += r1->toString();
                            outstr2 += r2->toString();
                        }
                    }
                    
                    delete merged;
                }
            }

            // stats the read after filtering
            config->getPostStats1()->statRead(r1);
            config->getPostStats2()->statRead(r2);

            readPassed++;
        }

        delete pair; pair = NULL;
        // if no trimming applied, r1 should be identical to or1
        if(r1 != or1 && r1 != NULL){
            delete r1;
            r1 = NULL;
        } else if(r1 == NULL){
            delete or1;
            or1 = NULL;
        }
        // if no trimming applied, r2 should be identical to or2
        if(r2 != or2 && r2 != NULL){
            delete r2;
            r2 = NULL;
        } else if(r2 == NULL){
            delete or2;
            or2 = NULL;
        }
    }

    //if(mOptions->debug) cCout("consuming reads pack end");
    mOutputMtx.lock();
    // normal output by left/right writer thread
    if(mRightWriter && mLeftWriter && (!outstr1.empty() || !outstr2.empty())) {
        // write PE
        char* ldata = new char[outstr1.size()];
        memcpy(ldata, outstr1.c_str(), outstr1.size());
        mLeftWriter->input(ldata, outstr1.size());

        char* rdata = new char[outstr2.size()];
        memcpy(rdata, outstr2.c_str(), outstr2.size());
        mRightWriter->input(rdata, outstr2.size());
    } else if(mLeftWriter && !singleOutput.empty()) {
        // write singleOutput
        char* ldata = new char[singleOutput.size()];
        memcpy(ldata, singleOutput.c_str(), singleOutput.size());
        mLeftWriter->input(ldata, singleOutput.size());
    } 
    
    if(mFailedWriter && !failedOutput.empty()){
        char* fdata = new char[failedOutput.size()];
        memcpy(fdata, failedOutput.c_str(), failedOutput.size());
        mFailedWriter->input(fdata, failedOutput.size());
    }
    mOutputMtx.unlock();
    config->markProcessed(pack->count);
    delete[] pack->data;
    delete pack;
    return true;
}
    
void PairEndProcessor::statInsertSize(Read* r1, Read* r2, OverlapResult& ov, int frontTrimmed1, int frontTrimmed2) {
    int isize = mOptions->insertSizeMax;
    if(ov.overlapped) {
        if(ov.offset > 0)
            isize = r1->length() + r2->length() - ov.overlap_len + frontTrimmed1 + frontTrimmed2;
        else
            isize = ov.overlap_len + frontTrimmed1 + frontTrimmed2;
    }

    if(isize > mOptions->insertSizeMax)
        isize = mOptions->insertSizeMax;

    mInsertSizeHist[isize]++;
}

bool PairEndProcessor::processRead(Read* r, ReadPair* originalPair, bool reversed) {
    // do something here
    return true;
}

void PairEndProcessor::initPackRepository() {
    mRepo.packBuffer = new ReadPairPack*[PACK_NUM_LIMIT];
    memset(mRepo.packBuffer, 0, sizeof(ReadPairPack*)*PACK_NUM_LIMIT);
    mRepo.writePos = 0;
    mRepo.readPos = 0;
    
}

void PairEndProcessor::destroyPackRepository() {
    delete[] mRepo.packBuffer;
    mRepo.packBuffer = NULL;
}

void PairEndProcessor::producePack(ReadPairPack* pack){
    mRepo.packBuffer[mRepo.writePos] = pack;
    mRepo.writePos++;
}

void PairEndProcessor::consumePack(ThreadConfig* config){
    ReadPairPack* data;
    mInputMtx.lock();
    while(mRepo.writePos <= mRepo.readPos) {
        usleep(1000);
        if(mProduceFinished) {
            mInputMtx.unlock();
            return;
        }
    }
    data = mRepo.packBuffer[mRepo.readPos];
    mRepo.readPos++;
    mInputMtx.unlock();
    processPairEnd(data, config);
}

void PairEndProcessor::producerTask(){
    if(mOptions->verbose)
        loginfo("start to load data");
    long lastReported = 0;
    int slept = 0;
    long readNum = 0;
    bool splitSizeReEvaluated = false;
    ReadPair** data = new ReadPair*[PACK_SIZE];
    memset(data, 0, sizeof(ReadPair*)*PACK_SIZE);
    FastqReaderPair reader(mOptions->in1, mOptions->in2, true, mOptions->phred64, mOptions->interleavedInput);
    int count=0;
    bool needToBreak = false;
    if(mOptions->debug) cCout("start to read data: " + mOptions->in1 + " and " + mOptions->in2);
    while(true){
        ReadPair* read = reader.read();
        // TODO: put needToBreak here is just a WAR for resolve some unidentified dead lock issue 
        if(!read || needToBreak){
            // the last pack
            ReadPairPack* pack = new ReadPairPack;
            pack->data = data;
            pack->count = count;
            producePack(pack);
            //delete[] data;
            data = NULL;
            if(read) {
                delete read;
                read = NULL;
            }
            break;
        }
        data[count] = read;
        count++;
        // configured to process only first N reads
        if(mOptions->readsToProcess >0 && count + readNum >= mOptions->readsToProcess) {
            needToBreak = true;
        }
        if(mOptions->verbose && count + readNum >= lastReported + 1000000) {
            lastReported = count + readNum;
            string msg = "loaded " + to_string((lastReported/1000000)) + "M read pairs";
            loginfo(msg);
        }
        // a full pack
        if(count == PACK_SIZE || needToBreak){
            ReadPairPack* pack = new ReadPairPack;
            pack->data = data;
            pack->count = count;
            producePack(pack);
            //re-initialize data for next pack
            data = new ReadPair*[PACK_SIZE];
            memset(data, 0, sizeof(ReadPair*)*PACK_SIZE);
            // if the consumer is far behind this producer, sleep and wait to limit memory usage
            while(mRepo.writePos - mRepo.readPos > PACK_IN_MEM_LIMIT){
                slept++;
                usleep(1000);
            }
            readNum += count;
            // if the writer threads are far behind this producer, sleep and wait
            // check this only when necessary
            if(readNum % (PACK_SIZE * PACK_IN_MEM_LIMIT) == 0 && mLeftWriter) {
                while( (mLeftWriter && mLeftWriter->bufferLength() > PACK_IN_MEM_LIMIT) || (mRightWriter && mRightWriter->bufferLength() > PACK_IN_MEM_LIMIT) || (mFailedWriter && mFailedWriter->bufferLength() > PACK_IN_MEM_LIMIT)){
                    slept++;
                    usleep(1000);
                }
            }
            // reset count to 0
            count = 0;
        }
    }

    //std::unique_lock<std::mutex> lock(mRepo.readCounterMtx);
    mProduceFinished = true;
    if(mOptions->verbose)
        loginfo("all reads loaded, start to monitor thread status");
    //lock.unlock();

    // if the last data initialized is not used, free it
    if(data != NULL){
        for(int i = 0; i < PACK_SIZE; ++i){
            if(data[i]!= NULL){
                delete data[i];
                data[i] = NULL;
            }
        }
        delete[] data;
        data = NULL;
    }
}

void PairEndProcessor::consumerTask(ThreadConfig* config){
    while(true) {
        if(config->canBeStopped()){
            mFinishedThreads++;
            break;
        }
        while(mRepo.writePos <= mRepo.readPos) {
            if(mProduceFinished)
                break;
            usleep(1000);
        }
        //std::unique_lock<std::mutex> lock(mRepo.readCounterMtx);
        if(mProduceFinished && mRepo.writePos == mRepo.readPos){
            mFinishedThreads++;
            if(mOptions->verbose) {
                string msg = "thread " + to_string(config->getThreadId() + 1) + " data processing completed";
                loginfo(msg);
            }
            //lock.unlock();
            break;
        }
        if(mProduceFinished){
            if(mOptions->verbose) {
                string msg = "thread " + to_string(config->getThreadId() + 1) + " is processing the " + to_string(mRepo.readPos) + " / " + to_string(mRepo.writePos) + " pack";
                loginfo(msg);
            }
            consumePack(config);
            //lock.unlock();
        } else {
            //lock.unlock();
            consumePack(config);
        }
    }

    if(mFinishedThreads == mOptions->thread) {
        if(mLeftWriter)
            mLeftWriter->setInputCompleted();
        if(mRightWriter)
            mRightWriter->setInputCompleted();
        if(mFailedWriter)
            mFailedWriter->setInputCompleted();
    }
    
    if(mOptions->verbose) {
        string msg = "thread " + to_string(config->getThreadId() + 1) + " finished";
        loginfo(msg);
    }
}

void PairEndProcessor::writeTask(WriterThread* config){
    while(true) {
        if(config->isCompleted()){
            // last check for possible threading related issue
            config->output();
            break;
        }
        config->output();
    }

    if(mOptions->verbose) {
        string msg = config->getFilename() + " writer finished";
        loginfo(msg);
    }
}