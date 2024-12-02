#include "evaluator.h"
#include "fastqreader.h"
#include <map>
#include <memory.h>
#include "nucleotidetree.h"
#include "knownadapters.h"

Evaluator::Evaluator(Options* opt){
    mOptions = opt;
}


Evaluator::~Evaluator(){
}

bool Evaluator::isTwoColorSystem() {
    FastqReader reader(mOptions->in1);

    Read* r = reader.read();

    if(!r)
        return false;

    // NEXTSEQ500, NEXTSEQ 550, NOVASEQ
    if(starts_with(r->mName, "@NS") || starts_with(r->mName, "@NB") || starts_with(r->mName, "@A0")) {
        delete r;
        return true;
    }

    delete r;
    return false;
}

void Evaluator::evaluateSeqLen() {
    if(!mOptions->in1.empty())
        mOptions->seqLen1 = computeSeqLen(mOptions->in1);
    if(!mOptions->in2.empty())
        mOptions->seqLen2 = computeSeqLen(mOptions->in2);
}

int Evaluator::computeSeqLen(string filename) {
    FastqReader reader(filename);

    long records = 0;
    bool reachedEOF = false;

    // get seqlen
    int seqlen=0;
    while(records < 1000) {
        Read* r = reader.read();
        if(!r) {
            reachedEOF = true;
            break;
        }
        int rlen = r->length();
        if(rlen > seqlen)
            seqlen = rlen;
        records ++;
        delete r;
    }

    return seqlen;
}

void Evaluator::evaluateReadNum(long& readNum) {
    FastqReader reader(mOptions->in1);

    const long READ_LIMIT = 512*1024;
    const long BASE_LIMIT = 151 * 512*1024;
    long records = 0;
    long bases = 0;
    size_t firstReadPos = 0;

    size_t bytesRead;
    size_t bytesTotal;

    bool reachedEOF = false;
    bool first = true;
    while(records < READ_LIMIT && bases < BASE_LIMIT) {
        Read* r = reader.read();
        if(!r) {
            reachedEOF = true;
            break;
        }
        if(first) {
            reader.getBytes(bytesRead, bytesTotal);
            firstReadPos = bytesRead;
            first = false;
        }
        records++;
        bases += r->length();
        delete r;
    }

    readNum = 0;
    if(reachedEOF){
        readNum = records;
    } else if(records>0) {
        // by the way, update readNum so we don't need to evaluate it if splitting output is enabled
        reader.getBytes(bytesRead, bytesTotal);
        double bytesPerRead = (double)(bytesRead - firstReadPos) / (double) records;
        // increase it by 1% since the evaluation is usually a bit lower due to bad quality causes lower compression rate
        readNum = (long) (bytesTotal*1.01 / bytesPerRead);
    }
}

string Evaluator::evalAdapterAndReadNum(long& readNum, bool isR2) {
    string filename = mOptions->in1;
    if(isR2)
        filename = mOptions->in2;
    FastqReader reader(filename);
    // stat up to 256K reads
    const long READ_LIMIT = 256*1024;
    const long BASE_LIMIT = 151 * READ_LIMIT;
    long records = 0;
    long bases = 0;
    size_t firstReadPos = 0;

    size_t bytesRead;
    size_t bytesTotal;

    Read** loadedReads = new Read*[READ_LIMIT];
    memset(loadedReads, 0, sizeof(Read*)*READ_LIMIT);
    bool reachedEOF = false;
    bool first = true;
    while(records < READ_LIMIT && bases < BASE_LIMIT) {
        Read* r = reader.read();
        if(!r) {
            reachedEOF = true;
            break;
        }
        if (r->length() >= 100000) continue;
        if (first) {
            reader.getBytes(bytesRead, bytesTotal);
            firstReadPos = bytesRead;
            first = false;
        }
        int rlen = r->length();
        bases += rlen;
        loadedReads[records] = r;
        records++;
    }

    readNum = 0;
    if(reachedEOF){
        readNum = records;
    } else if(records>0) {
        // by the way, update readNum so we don't need to evaluate it if splitting output is enabled
        reader.getBytes(bytesRead, bytesTotal);
        double bytesPerRead = (double)(bytesRead - firstReadPos) / (double) records;
        // increase it by 1% since the evaluation is usually a bit lower due to bad quality causes lower compression rate
        readNum = (long) (bytesTotal*1.01 / bytesPerRead);
    }
    // we need at least 10000 valid records to evaluate
    if(records < 10000) {
        for(int r=0; r<records; r++) {
            delete loadedReads[r];
            loadedReads[r] = NULL;
        }
        delete[] loadedReads;
        return "";
    }
    // we have to shift last cycle for evaluation since it is so noisy, especially for Illumina data
    const int shiftTail = max(1, mOptions->trim.tail1);

    // why we add trim_tail here? since the last cycle are usually with low quality and should be trimmed
    const int keylen = 10;
    int size = 1 << (keylen*2 );
    unsigned int* counts = new unsigned int[size];
    memset(counts, 0, sizeof(unsigned int)*size);
    for(int i=0; i<records; i++) {
        Read* r = loadedReads[i];
        const char* data = r->mSeq.mStr.c_str();
        int key = -1;
        for(int pos = 20; pos <= r->length()-keylen-shiftTail; pos++) {
            key = seq2int(r->mSeq.mStr, pos, keylen, key);
            if(key >= 0) {
                counts[key]++;
            }
        }
    }
    // set AAAAAAAAAA = 0;
    counts[0] = 0;

    // get the top N
    const int topnum = 10;
    int topkeys[topnum] = {0};
    long total = 0;
    for(int k=0; k<size; k++) {
        int atcg[4] = {0};
        for(int i=0; i<keylen; i++) {
            int baseOfBit = (k >> (i*2)) & 0x03;
            atcg[baseOfBit]++;
        }
        bool lowComplexity = false;
        for(int b=0; b<4; b++) {
            if(atcg[b] >= keylen-4)
                lowComplexity=true;
        }
        if(lowComplexity)
            continue;
        // too many GC
        if(atcg[2] + atcg[3] >= keylen-2)
            continue;

        // starts with GGGG
        if( k>>12 == 0xff)
            continue;

        unsigned int val = counts[k];
        total += val;
        for(int t=topnum-1; t>=0; t--) {
            // reach the middle
            if(val < counts[topkeys[t]]){
                if(t<topnum-1) {
                    for(int m=topnum-1; m>t+1; m--) {
                        topkeys[m] = topkeys[m-1];
                    }
                    topkeys[t+1] = k;
                }
                break;
            } else if(t == 0) { // reach the top
                for(int m=topnum-1; m>t; m--) {
                    topkeys[m] = topkeys[m-1];
                }
                topkeys[t] = k;
            }
        }
    }
    const int FOLD_THRESHOLD = 20;
    for(int t=0; t<topnum; t++) {
        int key = topkeys[t];
        string seq = int2seq(key, keylen);
        if(key == 0)
            continue;
        long count = counts[key];
        if(count<10 || count*size < total * FOLD_THRESHOLD)
            break;
        // skip low complexity seq
        int diff = 0;
        for(int s=0; s<seq.length() - 1; s++) {
            if(seq[s] != seq[s+1])
                diff++;
        }
        if(diff <3){
            continue;
        }
        string adapter = getAdapterWithSeed(key, loadedReads, records, keylen);
        if(!adapter.empty()){
            delete[] counts;
            for(int r=0; r<records; r++) {
                delete loadedReads[r];
                loadedReads[r] = NULL;
            }
            delete[] loadedReads;
            return adapter;
        }
    }
    delete[] counts;
    for(int r=0; r<records; r++) {
        delete loadedReads[r];
        loadedReads[r] = NULL;
    }
    delete[] loadedReads;
    return "";

}

string Evaluator::getAdapterWithSeed(int seed, Read** loadedReads, long records, int keylen) {
    // we have to shift last cycle for evaluation since it is so noisy, especially for Illumina data
    const int shiftTail = max(1, mOptions->trim.tail1);
    NucleotideTree forwardTree(mOptions);
    // forward search
    for(int i=0; i<records; i++) {
        Read* r = loadedReads[i];
        const char* data = r->mSeq.mStr.c_str();
        int key = -1;
        for(int pos = 20; pos <= r->length()-keylen-shiftTail; pos++) {
            key = seq2int(r->mSeq.mStr, pos, keylen, key);
            if(key == seed) {
                forwardTree.addSeq(r->mSeq.mStr.substr(pos+keylen, r->length()-keylen-shiftTail-pos));
            }
        }
    }
    bool reachedLeaf = true;
    string forwardPath = forwardTree.getDominantPath(reachedLeaf);
    NucleotideTree backwardTree(mOptions);
    // backward search
    for(int i=0; i<records; i++) {
        Read* r = loadedReads[i];
        const char* data = r->mSeq.mStr.c_str();
        int key = -1;
        for(int pos = 20; pos <= r->length()-keylen-shiftTail; pos++) {
            key = seq2int(r->mSeq.mStr, pos, keylen, key);
            if(key == seed) {
                string seq =  r->mSeq.mStr.substr(0, pos);
                string rcseq = reverse(seq);
                backwardTree.addSeq(rcseq);
            }
        }
    }
    string backwardPath = backwardTree.getDominantPath(reachedLeaf);
    string adapter = reverse(backwardPath) + int2seq(seed, keylen) + forwardPath;
    if(adapter.length()>60)
        adapter.resize(60);
    string matchedAdapter = matchKnownAdapter(adapter);
    if(!matchedAdapter.empty()) {
        map<string, string> knownAdapters = getKnownAdapter();
        //std::cerr << knownAdapters[matchedAdapter] << endl << matchedAdapter << endl;
        return matchedAdapter;
    } else {
        if(reachedLeaf) {
            //std::cerr << adapter << endl;
            return adapter;
        } else {
            return "";
        }
    }
}

string Evaluator::matchKnownAdapter(string seq) {
    map<string, string> knownAdapters = getKnownAdapter();
    map<string, string>::iterator iter;
    for(iter = knownAdapters.begin(); iter != knownAdapters.end(); iter++) {
        string adapter = iter->first;
        string desc = iter->second;
        if(seq.length()<adapter.length()) {
            continue;
        }
        int diff = 0;
        for(int i=0; i<adapter.length() && i<seq.length(); i++) {
            if(adapter[i] != seq[i])
                diff++;
        }
        if(diff == 0)
            return adapter;
    }
    return "";
}

string Evaluator::int2seq(unsigned int val, int seqlen) {
    char bases[4] = {'A', 'T', 'C', 'G'};
    string ret(seqlen, 'N');
    int done = 0;
    while(done < seqlen) {
        ret[seqlen - done - 1] = bases[val & 0x03];
        val = (val >> 2);
        done++;
    }
    return ret;
}

int Evaluator::seq2int(string& seq, int pos, int keylen, int lastVal) {
    int rlen = seq.length();
    if(lastVal >= 0) {
        const int mask = (1 << (keylen*2 )) - 1;
        int key = (lastVal<<2) & mask;
        char base = seq[pos + keylen - 1];
        switch (base) {
            case 'A':
                key += 0;
                break;
            case 'T':
                key += 1;
                break;
            case 'C':
                key += 2;
                break;
            case 'G':
                key += 3;
                break;
            default:
                // N or anything else
                return -1;
        }
        return key;
    } else {
        int key = 0;
        for(int i=pos; i<keylen+pos; i++) {
            key = (key << 2);
            char base = seq[i];
            switch (base) {
                case 'A':
                    key += 0;
                    break;
                case 'T':
                    key += 1;
                    break;
                case 'C':
                    key += 2;
                    break;
                case 'G':
                    key += 3;
                    break;
                default:
                    // N or anything else
                    return -1;
            }
        }
        return key;
    }
}

std::string Evaluator::getMarkerDirection(Options*& opt) {
    std::cerr << "starting to assess the marker primer direction" << std::endl;
    unsigned int records = 0;
    std::map<std::string, std::map<std::string, int>> mapped_marker_freq;// marker, forward/reverse/both, number of records
    if (opt->isPaired()) {
        FastqReaderPair reader(opt->in1, opt->in2);
        while (records < 10000) {
             ReadPair* r = reader.read();
             if (!r) {
                 break;
             }
             std::pair<std::string, std::string> mar_direct; // marker, direction
             int num_mismatches = 100;
             if(opt->var == "ssr"){
                for(auto & it : opt->mLocVars.refLocMap){
                    if(r->mLeft->mSeq.length() < (it.second.fp.mStr.length() + it.second.rp.mStr.length())){
                        continue;
                    }
                    auto pf = doSimpleAlignment(opt, it.second.fp.mStr, r->mLeft);
                    if (pf > opt->mLocVars.locVarOptions.maxMismatchesPSeq){
                        continue;
                    }
                    if(r->mRight->mSeq.length() < (it.second.fp.mStr.length() + it.second.rp.mStr.length())){
                        continue;
                    }
                    auto prf = doSimpleAlignment(opt, it.second.rp.mStr, r->mRight);
                    auto prr = doSimpleAlignment(opt, it.second.rp.reverseComplement().mStr, r->mRight);
                    if(std::max(prf, prr) > opt->mLocVars.locVarOptions.maxMismatchesPSeq){
                        continue;
                    } else {
                        if (prf < prr){
                            if (prf < num_mismatches){
                                num_mismatches = prf;
                                mar_direct = std::make_pair(it.second.name, "forward");
                            }
                        } else {
                            if(prr < num_mismatches){
                                num_mismatches = prr;
                                mar_direct = std::make_pair(it.second.name, "reverse");
                            }
                        }
                    }
                }
             } else {
                for(auto & it : opt->mLocSnps.refLocMap){
                    if(r->mLeft->mSeq.length() < (it.second.fp.mStr.length() + it.second.rp.mStr.length())){
                        continue;
                    }
                    auto pf = doSimpleAlignment(opt, it.second.fp.mStr, r->mLeft);
                    if (pf > opt->mLocVars.locVarOptions.maxMismatchesPSeq){
                        continue;
                    }
                    if(r->mRight->mSeq.length() < (it.second.fp.mStr.length() + it.second.rp.mStr.length())){
                        continue;
                    }
                    auto prf = doSimpleAlignment(opt, it.second.rp.mStr, r->mRight);
                    auto prr = doSimpleAlignment(opt, it.second.rp.reverseComplement().mStr, r->mRight);
                    if(std::max(prf, prr) > opt->mLocVars.locVarOptions.maxMismatchesPSeq){
                        continue;
                    } else {
                        if (prf < prr){
                            if (prf < num_mismatches){
                                num_mismatches = prf;
                                mar_direct = std::make_pair(it.second.name, "forward");
                            }
                        } else {
                            if(prr < num_mismatches){
                                num_mismatches = prr;
                                mar_direct = std::make_pair(it.second.name, "reverse");
                            }
                        }
                    }
                }
             }
             records++;
             delete r;
         }
    } else {
        FastqReader reader(opt->in1);
        while(records < 10000){
            Read* r = reader.read();
            if(!r){
                break;
            }
            std::pair<std::string, std::string> mar_direct; // marker, direction
            int num_mismatches = 100;
            if (opt->var == "ssr") {
                for(auto & it : opt->mLocVars.refLocMap){
                    std::cout << "marker " << it.first << std::endl;
                    if(r->mSeq.length() < (it.second.fp.mStr.length() + it.second.rp.mStr.length())){
                        continue;
                    }
                    auto pf = doSimpleAlignment(opt, it.second.fp.mStr, r);
                    if (pf > opt->mLocVars.locVarOptions.maxMismatchesPSeq){
                        continue;
                    }
                    auto prf = doSimpleAlignment(opt, it.second.rp.mStr, r);
                    auto prr = doSimpleAlignment(opt, it.second.rp.reverseComplement().mStr, r);
                    if(std::max(prf, prr) > opt->mLocVars.locVarOptions.maxMismatchesPSeq){
                        continue;
                    } else {
                        if (prf < prr){
                            if (prf < num_mismatches){
                                num_mismatches = prf;
                                mar_direct = std::make_pair(it.first, "forward");
                            }
                        } else {
                            if(prr < num_mismatches){
                                num_mismatches = prr;
                                mar_direct = std::make_pair(it.first, "reverse");
                            }
                        }
                    }
                }
            } else {
                for(auto & it : opt->mLocSnps.refLocMap){
                    if(r->mSeq.length() < (it.second.fp.mStr.length() + it.second.rp.mStr.length())){
                        continue;
                    }
                    auto pf = doSimpleAlignment(opt, it.second.fp.mStr, r);
                    if (pf > opt->mLocVars.locVarOptions.maxMismatchesPSeq){
                        continue;
                    }
                    auto prf = doSimpleAlignment(opt, it.second.rp.mStr, r);
                    auto prr = doSimpleAlignment(opt, it.second.rp.reverseComplement().mStr, r);

                    if(std::max(prf, prr) > opt->mLocVars.locVarOptions.maxMismatchesPSeq){
                        continue;
                    } else {
                        if (prf < prr){
                            if (prf < num_mismatches){
                                num_mismatches = prf;
                                mar_direct = std::make_pair(it.second.name, "forward");
                            }
                        } else {
                            if(prr < num_mismatches){
                                num_mismatches = prr;
                                mar_direct = std::make_pair(it.second.name, "reverse");
                            }
                        }
                    }
                }
            }
            delete r;
            records++;
            
            if(!mar_direct.first.empty()){
                mapped_marker_freq[mar_direct.first][mar_direct.second]++;
            }
        }
    }
    if(mapped_marker_freq.empty()){
        return "both";
    }
    std::map<std::string, int> mar_direction_map;//forward resverse int
    for (const auto& it : mapped_marker_freq) {
        for(const auto & itt : it.second){
            std::cout << it.first << " " << itt.first << " " << itt.second << "\n";
        }
        auto forward_it = it.second.find("forward");
        int forward_n = 0;
        if (forward_it != it.second.end()) {
            forward_n = forward_it->second;
        }
        auto reverse_it = it.second.find("reverse");
        int reverse_n = 0;
        if (reverse_it != it.second.end()) {
            reverse_n = reverse_it->second;
        }
        auto forward_per = getPer(forward_n, reverse_n + forward_n);
        if (forward_per >= 0.8){
            mar_direction_map["forward"]++;
        } else if(forward_per <= 0.2){
            mar_direction_map["reverse"]++;
        } else {
            mar_direction_map["both"]++;
        }
    }
    if(mar_direction_map.empty()){
        return "both";
    }
    auto max_dict = getMaxKeyValue(mar_direction_map);
    std::cerr << "the marker direction has been assessed: " << max_dict.first << std::endl;
    return max_dict.first;
}

 std::string Evaluator::getSexMarker(Options*& opt) {
     std::string marker = "";
     long records = 0;
     std::map<std::string, double> sexEvaMap;
     if (opt->isPaired()) {
         FastqReaderPair reader(opt->in1, opt->in2);
         while (records < 100000) {
             ReadPair* r = reader.read();
             if (!r) {
                 break;
             }
             std::map<std::string, double> tsexMap;
             for (auto& it : opt->sexMap) {
                 std::string goLR = "";  // left, right, both;
                 if (r->mLeft->length() > (it.second.primerF.length() + it.second.primerR.length())) {
                     if (r->mRight->length() > (it.second.primerF.length() + it.second.primerR.length())) {
                         goLR = "both";
                     } else {
                         goLR = "left";
                     }
                 } else {
                     if (r->mRight->length() > (it.second.primerF.length() + it.second.primerR.length())) {
                         goLR = "right";
                     }
                 }
                 if (goLR.empty()) {
                     continue;
                 } else if (goLR == "both") {
                     auto pf = doSimpleAlignment(opt, it.second.primerF.mStr, r->mLeft);
                     auto pr = doSimpleAlignment(opt, it.second.primerR.mStr, r->mLeft);

                     auto pfr = doSimpleAlignment(opt, it.second.primerF.reverseComplement().mStr, r->mLeft);
                     auto prr = doSimpleAlignment(opt, it.second.primerR.reverseComplement().mStr, r->mLeft);

                     if ((pf + pr) < (pfr + prr)) {
                        if(pf <= opt->mLocVars.locVarOptions.maxMismatchesPSeq && pr <= opt->mLocVars.locVarOptions.maxMismatchesPSeq){
                            tsexMap[it.first] = (pf + pr) / (it.second.primerF.length() + it.second.primerR.length());
                        }
                     } else {
                        if(pfr <= opt->mLocVars.locVarOptions.maxMismatchesPSeq && prr <= opt->mLocVars.locVarOptions.maxMismatchesPSeq){
                            tsexMap[it.first] = (pfr + prr) / (it.second.primerF.length() + it.second.primerR.length());
                        }
                     }

                     auto rpf = doSimpleAlignment(opt, it.second.primerF.mStr, r->mRight);
                     auto rpr = doSimpleAlignment(opt, it.second.primerR.mStr, r->mRight);

                     auto rpfr = doSimpleAlignment(opt, it.second.primerF.reverseComplement().mStr, r->mRight);
                     auto rprr = doSimpleAlignment(opt, it.second.primerR.reverseComplement().mStr, r->mRight);

                     if ((rpf + rpr) < (rpfr + rprr)) {
                        if(rpf <= opt->mLocVars.locVarOptions.maxMismatchesPSeq && rpr <= opt->mLocVars.locVarOptions.maxMismatchesPSeq){
                            tsexMap[it.first] = (tsexMap[it.first] + (rpf + rpr) / (it.second.primerF.length() + it.second.primerR.length())) / 2;
                        }
                     } else {
                        if(rpfr <= opt->mLocVars.locVarOptions.maxMismatchesPSeq && rprr <= opt->mLocVars.locVarOptions.maxMismatchesPSeq){
                            tsexMap[it.first] = (tsexMap[it.first] + (rpfr + rprr) / (it.second.primerF.length() + it.second.primerR.length())) / 2;
                        }
                     }
                 } else {
                     Read* tmpr = nullptr;
                     if (goLR == "left") {
                         tmpr = r->mLeft;
                     } else if (goLR == "right") {
                         tmpr = r->mRight;
                     }

                     auto pf = doSimpleAlignment(opt, it.second.primerF.mStr, tmpr);
                     auto pr = doSimpleAlignment(opt, it.second.primerR.mStr, tmpr);

                     auto pfr = doSimpleAlignment(opt, it.second.primerF.reverseComplement().mStr, tmpr);
                     auto prr = doSimpleAlignment(opt, it.second.primerR.reverseComplement().mStr, tmpr);

                      if ((pf + pr) < (pfr + prr)) {
                        if(pf <= opt->mLocVars.locVarOptions.maxMismatchesPSeq && pr <= opt->mLocVars.locVarOptions.maxMismatchesPSeq){
                            tsexMap[it.first] = (pf + pr) / (it.second.primerF.length() + it.second.primerR.length());
                        }
                     } else {
                        if(pfr <= opt->mLocVars.locVarOptions.maxMismatchesPSeq && prr <= opt->mLocVars.locVarOptions.maxMismatchesPSeq){
                            tsexMap[it.first] = (pfr + prr) / (it.second.primerF.length() + it.second.primerR.length());
                        }
                     }
                     //delete tmpr;
                     tmpr = nullptr;
                 }
             }
             records++;
             delete r;
              if (tsexMap.empty()) continue;
              auto minVK = getMinKeyValue(tsexMap);
              sexEvaMap[minVK.first] = ((sexEvaMap[minVK.first] + minVK.second) / 2);
         }
    } else {
        FastqReader reader(opt->in1);
        while (records < 100000) {
            Read* r = reader.read();
            if (!r) {
                break;
            }

            std::map<std::string, double> tsexMap;
            for (auto& it : opt->sexMap) {
                if (r->length() < (it.second.primerF.length() + it.second.primerR.length())) {
                    continue;
                }
                auto pf = doSimpleAlignment(opt, it.second.primerF.mStr, r);
                auto pr = doSimpleAlignment(opt, it.second.primerR.mStr, r);

                auto pfr = doSimpleAlignment(opt, it.second.primerF.reverseComplement().mStr, r);
                auto prr = doSimpleAlignment(opt, it.second.primerR.reverseComplement().mStr, r);

                if ((pf + pr) < (pfr + prr)) {
                        if(pf <= opt->mLocVars.locVarOptions.maxMismatchesPSeq && pr <= opt->mLocVars.locVarOptions.maxMismatchesPSeq){
                            tsexMap[it.first] = (pf + pr) / (it.second.primerF.length() + it.second.primerR.length());
                        }
                     } else {
                        if(pfr <= opt->mLocVars.locVarOptions.maxMismatchesPSeq && prr <= opt->mLocVars.locVarOptions.maxMismatchesPSeq){
                            tsexMap[it.first] = (pfr + prr) / (it.second.primerF.length() + it.second.primerR.length());
                        }
                }
            }
            delete r;
            records++;
            if (tsexMap.empty()) continue;
            auto minVK = getMinKeyValue(tsexMap);
            sexEvaMap[minVK.first] = ((sexEvaMap[minVK.first] + minVK.second) / 2);
        }
    }

    if (!sexEvaMap.empty()) {
        marker = getMinKeyValue(sexEvaMap).first;
    }
    return marker;
}
 int Evaluator::doSimpleAlignment(Options*& opt, const string& query, Read*& r) {
     const char* tData = r->mSeq.mStr.c_str();
     const char* qData = query.c_str();
     EdlibAlignResult result =
         edlibAlign(qData, query.length(), tData, r->mSeq.length(), 
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
         edlibFreeAlignResult(result);
         if (indelSet.empty()) {
             return snpsSet.size();
         } else {
             return query.length();
         }
    } else {
        edlibFreeAlignResult(result);
        return query.length();
    }
}

bool Evaluator::test() {
    Evaluator eval(NULL);
    string s = "ATCGATCGAT";
    std::cerr << eval.int2seq(eval.seq2int(s, 0, 10, -1), 10) << endl;
    return eval.int2seq(eval.seq2int(s, 0, 10, -1), 10) == s;
}