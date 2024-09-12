#include <stdio.h>
#include <time.h>
#include <sstream>
#include "util.h"
#include "options.h"
#include "processor.h"
#include "evaluator.h"
#include "fastqreader.h"
#include "unittest.h"
#include "cmdline.h"
#include "htmlreporterall.h"

string command;
mutex logmtx;

int main(int argc, char* argv[]){
    // display version info if no argument is given
    if(argc == 1) {
        std::cerr << "Seq2Sat: an ultra-fast and comprehensive genetic variation identification tool for raw sequencing data." << endl << "version " << SEQ2SAT_VER << endl;
    }
    if (argc == 2 && strcmp(argv[1], "test")==0){
        UnitTest tester;
        tester.run();
        return 0;
    }
    if (argc == 2 && (strcmp(argv[1], "-v")==0 || strcmp(argv[1], "--version")==0)){
        std::cerr << "seq2sat " << SEQ2SAT_VER << endl;
        return 0;
    }
    cmdline::parser cmd;
    cmd.add<string>("var", 0, "genetic variance, must be either microsatellite/ssr or snp", false, "");
    cmd.add<string>("in1", 'i', "read1 input file name", false, "");
    cmd.add<string>("in2", 'I', "read2 input file name", false, "");
    cmd.add<string>("prefix", 'X', "prefix name for output files, eg: sample01", false, "");
    cmd.add("outFReads", 0, "If specified, off-target reads will be outputed in a file");
    cmd.add<string>("out1", 'o', "file name to store read1 with on-target sequences", false, "");
    cmd.add<string>("out2", 'O', "file name to store read2 with on-target sequences", false, "");
    cmd.add("nanopore_default", 0, "If specified, use the default settings (mainly about Q score) of nanopore sequencing");

    cmd.add<string>("loc", 0, "loci file containing loci names, 5'primer sequence, reverse complement of 3'primer sequence, 5'flank region, 3'flank region, repeat unit and reference microsatellite repeat array, separated by '\t", false, "");
    cmd.add("revCom", 0, "if your reverse primer sequence in the loc file is not reverse complentary, please specify it");
    cmd.add<int>("minSeqs", 0, "minimum number of reads for a genotype, default: 5", false, 5);
    cmd.add<int>("maxMismatchesPSeq", 0, "maximum mismatches for primer sequences 2", false, 2);
    cmd.add("noPlot", 0, "If specified, do not plot");
    
    cmd.add<string>("mode", 0, "specify the sequence alignment mode: NW (default) | HW | SHW", false, "NW");
    cmd.add<int>("maxScore", 0, "specify the maximum score of sequence alignment with sore > maxScore will be discarded, default value is -1, and no sequence will be discarded.",false, -1);
    cmd.add<int>("numBestSeqs", 0, "Score will be calculated only for N best sequences (best = with smallest score). If N = 0 then all sequences will be calculated.", 0);
    cmd.add("notFindAlignment", 0, "If specified, alignment path will be not found and printed. This may significantly speed up the calculation");
    cmd.add("findStartLocation", 0, "If specified, start locations will be found and printed. Each start location corresponds to one end location. This may somewhat slow down the calculation, but is still faster then finding alignment path and does not consume any extra memory.");
    cmd.add<string>("format", 0, "NICE|CIG_STD|CIG_EXT  Format that will be used to print alignment path, can be used only with -p. NICE will give visually attractive format, CIG_STD will give standard cigar format and CIG_EXT will give extended cigar format. [default: NICE]", false, "NICE");
    cmd.add<int>("core", 0, "Core part of calculation will be repeated N times. This is useful only for performance measurement, when single execution is too short to measure.", false, 1);
    cmd.add("silentAlignment", 0, "If specified, there will be no score or alignment output");
    cmd.add("printResults", 0, "If specified, alignment results will be printed but with only 1 thread");
    
    //for ssr;
    cmd.add<double>("maxMismatchesPer4FR", 0, "maximum percentage mismatches for the forward and reverse flanking regions, default 0.3 (30%) ", false, 0.3);
    cmd.add<int>("minSeqsPercentage", 0, "minimum percentage (%) reads against largest peak for a genotype, default: 5 (5%)", false, 5);
    cmd.add<int>("minWarningSeqs", 0, "minimum number of reads for warning a genotype, default: 50", false, 50);
    cmd.add<double>("hlRatio1", 0, "ratio of loci sizes of largest and second largest numbers of reads when the length difference = 1 ssr unit, default: 0.4", false, 0.4);
    cmd.add<double>("hlRatio2", 0, "ratio of loci sizes of largest and second largest numbers of reads when the length difference = 2 ssr unit, default: 0.2", false, 0.2);
    cmd.add<double>("maxVarRatio", 0, "ratio of two heter alleles based on variations either in flanking regions or MRA, the ideal is 1, default: 1.5", false, 1.5);
    
    //for snp;
    cmd.add<int>("htJetter", 0, "jetter rate for heter loci, eg 25 means the percentage of reads with SNPs to total reads are from 25 - 75 %, must be coupled with hmPer, default: 25", false, 25);
    cmd.add<int>("hmPerH", 0, "allele is considered as homo when its reads against total reads is > 90 % and there are at least 2 true SNPs, must be coupled with htJetter and hmPerL, must be > hmPerL. default: 90", false, 90);
    cmd.add<int>("hmPerL", 0, "allele is considered as homo when its reads against total reads is > 80 % and there are at most 1 true SNP, must be coupled with htJetter and hmPerH, must be < hmPerH and > htJetter + 50. default: 85", false, 85);
    cmd.add<int>("minSeqsPerSnp", 0, "minimum percentage (%) reads against largest peak for a genotype, default: 10 (10%)", false, 10);
    cmd.add<int>("minReads4Filter", 0, "minimum reads for filtering read variant. if the maximum reads of haplotype is more than this, the low abundance read variants will be filtered, otherwise will be kept. This is used for the shallow sequencing. default: 50", false, 50);
    cmd.add<int>("maxRows4Align", 0, "maximum rows for alignment table, must be > 2 rows. default: 6", false, 6);
    cmd.add("noSnpPlot", 0, "If specified, do not plot SNPs");
    cmd.add("noErrorPlot", 0, "If specified, do not plot error rate");
    
    //for sex
    cmd.add<string>("sex", 0, "sex loci file containing sex locus names, 5'primer sequence, reverse complement of 3'primer sequence, X/Z reference sequence, Y/W reference sequence, separated by '\t", false, "");
    cmd.add<unsigned int>("maxMismatchesSexPSeq", 0, "maximum number of mismatches for sex primers, default: 2", false, 2);
    cmd.add<unsigned int>("maxMismatchesSexRefSeq", 0, "maximum number of mismatches for sex reference sequences, default: 2", false, 2);
    cmd.add<double>("yxRatio", 0, "minimum ratio of numbers of reads Y/X to W/Z, default: 0.001", false, 0.001);
    cmd.add<int>("minReadsSexAllele", 0, "minimum number of reads assigned to each sex allele of X or Y; default: 10", false, 10);
    cmd.add<int>("minReadsSexVariant", 0, "minimum number of reads assigned to each sex variant of X or Y; default: 5", false, 5);
    
    cmd.add("debug", 0, "If specified, print debug");
    
    //cmd.add("dont_merge_overlapped_PE", 0, "don't merge the overlapped PE reads; this is off by default");
    
    // reporting
    cmd.add<string>("json", 'j', "the json format report file name", false, "seq2sat.json");
    cmd.add<string>("html", 'h', "the html format report file name", false, "seq2sat.html");
    cmd.add<string>("report_title", 'R', "should be quoted with \' or \", default is \"seq2sat report\"", false, "seq2sat report");

    // threading
    cmd.add<int>("thread", 'w', "worker thread number, default is 4", false, 4);

    // other I/O
    cmd.add("phred64", '6', "indicate the input is using phred64 scoring (it'll be converted to phred33, so the output will still be phred33)");
    cmd.add<int>("compression", 'z', "compression level for gzip output (1 ~ 9). 1 is fastest, 9 is smallest, default is 4.", false, 4);
    cmd.add("stdin", 0, "input from STDIN. If the STDIN is interleaved paired-end FASTQ, please also add --interleaved_in.");
    cmd.add("stdout", 0, "stream passing-filters reads to STDOUT. This option will result in interleaved FASTQ output for paired-end output. Disabled by default.");
    cmd.add("interleaved_in", 0, "indicate that <in1> is an interleaved FASTQ which contains both read1 and read2. Disabled by default.");
    cmd.add<int>("reads_to_process", 0, "specify how many reads/pairs to be processed. Default 0 means process all reads.", false, 0);
    cmd.add("dont_overwrite", 0, "don't overwrite existing files. Overwritting is allowed by default.");
    cmd.add("verbose", 'V', "output verbose log information (i.e. when every 1M reads are processed).");

    // adapter
    cmd.add("disable_adapter_trimming", 'A', "adapter trimming is enabled by default. If this option is specified, adapter trimming is disabled");
    cmd.add<string>("adapter_sequence", 'a', "the adapter for read1. For SE data, if not specified, the adapter will be auto-detected. For PE data, this is used if R1/R2 are found not overlapped.", false, "auto");
    cmd.add<string>("adapter_sequence_r2", 0, "the adapter for read2 (PE data only). This is used if R1/R2 are found not overlapped. If not specified, it will be the same as <adapter_sequence>", false, "auto");
    cmd.add<string>("adapter_fasta", 0, "specify a FASTA file to trim both read1 and read2 (if PE) by all the sequences in this FASTA file", false, "");
    cmd.add("detect_adapter_for_pe", 0, "by default, the auto-detection for adapter is for SE data input only, turn on this option to enable it for PE data.");

    // trimming
    cmd.add<int>("trim_front1", 'f', "trimming how many bases in front for read1, default is 0", false, 0);
    cmd.add<int>("trim_tail1", 't', "trimming how many bases in tail for read1, default is 0", false, 0);
    cmd.add<int>("max_len1", 'b', "if read1 is longer than max_len1, then trim read1 at its tail to make it as long as max_len1. Default 0 means no limitation", false, 0);
    cmd.add<int>("trim_front2", 'F', "trimming how many bases in front for read2. If it's not specified, it will follow read1's settings", false, 0);
    cmd.add<int>("trim_tail2", 'T', "trimming how many bases in tail for read2. If it's not specified, it will follow read1's settings", false, 0);
    cmd.add<int>("max_len2", 'B', "if read2 is longer than max_len2, then trim read2 at its tail to make it as long as max_len2. Default 0 means no limitation. If it's not specified, it will follow read1's settings", false, 0);

    // polyG tail trimming
    cmd.add<int>("poly_g_min_len", 0, "the minimum length to detect polyG in the read tail. 10 by default.", false, 10);
    cmd.add("disable_trim_poly_g", 'G', "disable polyG tail trimming, by default trimming is automatically enabled for Illumina NextSeq/NovaSeq data");
    
    // polyX tail trimming
    cmd.add("trim_poly_x", 'x', "enable polyX trimming in 3' ends.");
    cmd.add<int>("poly_x_min_len", 0, "the minimum length to detect polyX in the read tail. 10 by default.", false, 10);

    // cutting by quality
    cmd.add("cut_front", '5', "move a sliding window from front (5') to tail, drop the bases in the window if its mean quality < threshold, stop otherwise.");
    cmd.add("cut_tail", '3', "move a sliding window from tail (3') to front, drop the bases in the window if its mean quality < threshold, stop otherwise.");
    cmd.add("cut_right", 'r', "move a sliding window from front to tail, if meet one window with mean quality < threshold, drop the bases in the window and the right part, and then stop.");
    cmd.add<int>("cut_window_size", 'W', "the window size option shared by cut_front, cut_tail or cut_sliding. Range: 1~1000, default: 4", false, 4);
    cmd.add<int>("cut_mean_quality", 'M', "the mean quality requirement option shared by cut_front, cut_tail or cut_sliding. Range: 1~36 default: 20 (Q20)", false, 20);
    cmd.add<int>("cut_front_window_size", 0, "the window size option of cut_front, default to cut_window_size if not specified", false, 4);
    cmd.add<int>("cut_front_mean_quality", 0, "the mean quality requirement option for cut_front, default to cut_mean_quality if not specified", false, 20);
    cmd.add<int>("cut_tail_window_size", 0, "the window size option of cut_tail, default to cut_window_size if not specified", false, 4);
    cmd.add<int>("cut_tail_mean_quality", 0, "the mean quality requirement option for cut_tail, default to cut_mean_quality if not specified", false, 20);
    cmd.add<int>("cut_right_window_size", 0, "the window size option of cut_right, default to cut_window_size if not specified", false, 4);
    cmd.add<int>("cut_right_mean_quality", 0, "the mean quality requirement option for cut_right, default to cut_mean_quality if not specified", false, 20);

    // quality filtering
    cmd.add("disable_quality_filtering", 'Q', "quality filtering is enabled by default. If this option is specified, quality filtering is disabled");
    cmd.add<int>("qualified_quality_phred", 'q', "the quality value that a base is qualified. Default 20 (for Illumina and 10 for nanopore) means phred quality >=Q20 is qualified.", false, 20);
    cmd.add<int>("unqualified_percent_limit", 'u', "how many percents of bases are allowed to be unqualified (0~100). Default 40 means 40%", false, 40);
    cmd.add<int>("n_base_limit", 'n', "if one read's number of N base is >n_base_limit, then this read/pair is discarded. Default is 5", false, 5);
    cmd.add<int>("average_qual", 'e', "if one read's average quality score <avg_qual, then this read/pair is discarded. Default 0 means no requirement", false, 0);

    // length filtering
    cmd.add("disable_length_filtering", 'L', "length filtering is enabled by default. If this option is specified, length filtering is disabled");
    cmd.add<int>("length_required", 'l', "reads shorter than length_required will be discarded, default is 30.", false, 30);
    cmd.add<int>("length_limit", 0, "reads longer than length_limit will be discarded, default 0 means no limitation.", false, 0);

    // low complexity filtering
    cmd.add("low_complexity_filter", 'y', "enable low complexity filter. The complexity is defined as the percentage of base that is different from its next base (base[i] != base[i+1]).");
    cmd.add<int>("complexity_threshold", 'Y', "the threshold for low complexity filter (0~100). Default is 30, which means 30% complexity is required.", false, 30);

    // filter by indexes
    cmd.add<string>("filter_by_index1", 0, "specify a file contains a list of barcodes of index1 to be filtered out, one barcode per line", false, "");
    cmd.add<string>("filter_by_index2", 0, "specify a file contains a list of barcodes of index2 to be filtered out, one barcode per line", false, "");
    cmd.add<int>("filter_by_index_threshold", 0, "the allowed difference of index barcode for index filtering, default 0 means completely identical.", false, 0);
    
    // base correction in overlapped regions of paired end data
    cmd.add("no_correction", 'C', "disable base correction in overlapped regions (only for PE data), default is enabled");
    cmd.add<int>("overlap_len_require", 0, "the minimum length to detect overlapped region of PE reads. This will affect overlap analysis based PE merge, adapter trimming and correction. 30 by default.", false, 30);
    cmd.add<int>("overlap_diff_limit", 0, "the maximum number of mismatched bases to detect overlapped region of PE reads. This will affect overlap analysis based PE merge, adapter trimming and correction. 5 by default.", false, 5);
    cmd.add<int>("overlap_diff_percent_limit", 0, "the maximum percentage of mismatched bases to detect overlapped region of PE reads. This will affect overlap analysis based PE merge, adapter trimming and correction. Default 20 means 20%.", false, 20);

    // umi
    cmd.add("umi", 'U', "enable unique molecular identifier (UMI) preprocessing");
    cmd.add<string>("umi_loc", 0, "specify the location of UMI, can be (index1/index2/read1/read2/per_index/per_read, default is none", false, "");
    cmd.add<int>("umi_len", 0, "if the UMI is in read1/read2, its length should be provided", false, 0);
    cmd.add<string>("umi_prefix", 0, "if specified, an underline will be used to connect prefix and UMI (i.e. prefix=UMI, UMI=AATTCG, final=UMI_AATTCG). No prefix by default", false, "");
    cmd.add<int>("umi_skip", 0, "if the UMI is in read1/read2, seq2sat can skip several bases following UMI, default is 0", false, 0);

    cmd.parse_check(argc, argv);

    if(argc == 1) {
        std::cerr << cmd.usage() <<endl;
        return 0;
    }

    Options * opt = new Options();

    string seq2satProgPath = string(argv[0]);
    string seq2satDir = dirname(seq2satProgPath);
    
    // I/O
    opt->var = cmd.get<string>("var");
    opt->mEdOptions.mode = cmd.get<string>("mode");
    opt->mEdOptions.findAlignment = !cmd.exist("notFindAlignment");
    opt->mEdOptions.findStartLocation = cmd.exist("findStartLocation");
    opt->mEdOptions.silent = cmd.exist("silentAlignment");
    opt->mEdOptions.printRes = cmd.exist("printResults");
    opt->mEdOptions.format = cmd.get<string>("format");

    opt->sexFile = cmd.get<string>("sex");
    opt->mSex.mismatchesPF = opt->mSex.mismatchesPR = cmd.get<unsigned int>("maxMismatchesSexPSeq");
    opt->mSex.mismatchesRX = opt->mSex.mismatchesRY = cmd.get<unsigned int>("maxMismatchesSexRefSeq");
    opt->mSex.YXRationCuttoff = cmd.get<double>("yxRatio");
    opt->mSex.minReadsSexAllele = cmd.get<int>("minReadsSexAllele");
    opt->mSex.minReadsSexVariant = cmd.get<int>("minReadsSexVariant");

    opt->mLocVars.locVarOptions.maxScore = cmd.get<int>("maxScore");
    opt->mLocVars.locVarOptions.numBestSeqs = cmd.get<int>("numBestSeqs");
    opt->mLocVars.locVarOptions.coreRep = cmd.get<int>("core");
    opt->noPlot = cmd.exist("noPlot");
        
    if (opt->var == "snp" || opt->sexFile != "") {
        opt->mLocSnps.mLocSnpOptions.minSeqs = cmd.get<int>("minSeqs");
        opt->mLocSnps.mLocSnpOptions.minSeqsPer = (double) cmd.get<int>("minSeqsPerSnp") / 100.00;
        opt->mLocSnps.mLocSnpOptions.htJetter = (double) cmd.get<int>("htJetter") / 100.00;
        opt->mLocSnps.mLocSnpOptions.hmPerH = (double) cmd.get<int>("hmPerH") / 100.00;
        opt->mLocSnps.mLocSnpOptions.hmPerL = (double) cmd.get<int>("hmPerL") / 100.00;
        opt->mLocSnps.mLocSnpOptions.minReads4Filter = cmd.get<int>("minReads4Filter");
        opt->mLocSnps.mLocSnpOptions.maxRows4Align = cmd.get<int>("maxRows4Align");
        opt->noSnpPlot = cmd.exist("noSnpPlot");
        opt->noErrorPlot = cmd.exist("noErrorPlot");
    } else {
        opt->mLocVars.locVarOptions.maxMismatchesPSeq = cmd.get<int>("maxMismatchesPSeq");
        opt->mLocVars.locVarOptions.maxMismatchesPer4FR = cmd.get<double>("maxMismatchesPer4FR");
        opt->mLocVars.locVarOptions.minSeqs = cmd.get<int>("minSeqs");
        opt->mLocVars.locVarOptions.minWarningSeqs = cmd.get<int>("minWarningSeqs");
        opt->mLocVars.locVarOptions.minSeqsPer = cmd.get<int>("minSeqsPercentage");
        opt->mLocVars.locVarOptions.hlRatio1 = cmd.get<double>("hlRatio1");
        opt->mLocVars.locVarOptions.hlRatio2 = cmd.get<double>("hlRatio2");
        opt->mLocVars.locVarOptions.varRatio = cmd.get<double>("maxVarRatio");
    }
    opt->locFile = cmd.get<string>("loc");
    opt->revCom = cmd.exist("revCom");
    opt->nanopore_default = cmd.exist("nanopore_default");
    
    //opt->mergerOverlappedPE = cmd.exist("dont_merge_overlapped_PE") ? false : true;
    
    opt->compression = cmd.get<int>("compression");
    opt->readsToProcess = cmd.get<int>("reads_to_process");
    opt->phred64 = cmd.exist("phred64");
    opt->dontOverwrite = cmd.exist("dont_overwrite");
    opt->inputFromSTDIN = cmd.exist("stdin");
    opt->outputToSTDOUT = cmd.exist("stdout");
    opt->interleavedInput = cmd.exist("interleaved_in");
    opt->verbose = cmd.exist("verbose");
    opt->debug = cmd.exist("debug");

    // adapter cutting
    opt->adapter.enabled = !cmd.exist("disable_adapter_trimming");
    opt->adapter.detectAdapterForPE = cmd.exist("detect_adapter_for_pe");
    opt->adapter.sequence = cmd.get<string>("adapter_sequence");
    opt->adapter.sequenceR2 = cmd.get<string>("adapter_sequence_r2");
    opt->adapter.fastaFile = cmd.get<string>("adapter_fasta");
    if(opt->adapter.sequenceR2=="auto" && !opt->adapter.detectAdapterForPE && opt->adapter.sequence != "auto") {
        opt->adapter.sequenceR2 = opt->adapter.sequence;
    }
    if(!opt->adapter.fastaFile.empty()) {
        opt->loadFastaAdapters();
    }

    // trimming
    opt->trim.front1 = cmd.get<int>("trim_front1");
    opt->trim.tail1 = cmd.get<int>("trim_tail1");
    opt->trim.maxLen1 = cmd.get<int>("max_len1");
    // read2 settings follows read1 if it's not specified
    if(cmd.exist("trim_front2"))
        opt->trim.front2 = cmd.get<int>("trim_front2");
    else
        opt->trim.front2 = opt->trim.front1;
    
    if(cmd.exist("trim_tail2"))
        opt->trim.tail2 = cmd.get<int>("trim_tail2");
    else
        opt->trim.tail2 = opt->trim.tail1;
    
    if(cmd.exist("max_len2"))
        opt->trim.maxLen2 = cmd.get<int>("max_len2");
    else
        opt->trim.maxLen2 = opt->trim.maxLen1;

    // polyG tail trimming
    if(cmd.exist("disable_trim_poly_g")) {
        opt->polyGTrim.enabled = false;
    }
    opt->polyGTrim.minLen = cmd.get<int>("poly_g_min_len");

    // polyX tail trimming
    if(cmd.exist("trim_poly_x")) {
        opt->polyXTrim.enabled = true;
    }
    opt->polyXTrim.minLen = cmd.get<int>("poly_x_min_len");

    // sliding window cutting by quality
    opt->qualityCut.enabledFront = cmd.exist("cut_front");
    opt->qualityCut.enabledTail = cmd.exist("cut_tail");
    opt->qualityCut.enabledRight = cmd.exist("cut_right");

    opt->qualityCut.windowSizeShared = cmd.get<int>("cut_window_size");
    if(opt->nanopore_default){
        opt->qualityCut.qualityShared = 10;
    } else {
        opt->qualityCut.qualityShared = cmd.get<int>("cut_mean_quality");
    }

    if(cmd.exist("cut_front_window_size"))
        opt->qualityCut.windowSizeFront = cmd.get<int>("cut_front_window_size");
    else
        opt->qualityCut.windowSizeFront = opt->qualityCut.windowSizeShared;
    if(cmd.exist("cut_front_mean_quality"))
        if(opt->nanopore_default){
            opt->qualityCut.qualityFront = 10;
        } else {
            opt->qualityCut.qualityFront = cmd.get<int>("cut_front_mean_quality");
        }
    else
        opt->qualityCut.qualityFront = opt->qualityCut.qualityShared;

    if(cmd.exist("cut_tail_window_size"))
        opt->qualityCut.windowSizeTail = cmd.get<int>("cut_tail_window_size");
    else
        opt->qualityCut.windowSizeTail = opt->qualityCut.windowSizeShared;
    if(cmd.exist("cut_tail_mean_quality"))
        if(opt->nanopore_default){
            opt->qualityCut.qualityTail = 10;
        } else {
            opt->qualityCut.qualityTail = cmd.get<int>("cut_tail_mean_quality");
        }
    else
        opt->qualityCut.qualityTail = opt->qualityCut.qualityShared;

    if(cmd.exist("cut_right_window_size"))
        opt->qualityCut.windowSizeRight = cmd.get<int>("cut_right_window_size");
    else
        opt->qualityCut.windowSizeRight = opt->qualityCut.windowSizeShared;
    if(cmd.exist("cut_right_mean_quality"))
        if(opt->nanopore_default){
            opt->qualityCut.qualityRight = 10;
        } else {
            opt->qualityCut.qualityRight = cmd.get<int>("cut_right_mean_quality");
        }
    else
        opt->qualityCut.qualityRight = opt->qualityCut.qualityShared;

    // raise a warning if cutting option is not enabled but -W/-M is enabled
    if(!opt->qualityCut.enabledFront && !opt->qualityCut.enabledTail && !opt->qualityCut.enabledRight) {
        if(cmd.exist("cut_window_size") || cmd.exist("cut_mean_quality") 
            || cmd.exist("cut_front_window_size") || cmd.exist("cut_front_mean_quality") 
            || cmd.exist("cut_tail_window_size") || cmd.exist("cut_tail_mean_quality") 
            || cmd.exist("cut_right_window_size") || cmd.exist("cut_right_mean_quality"))
            std::cerr << "WARNING: you specified the options for cutting by quality, but forogt to enable any of cut_front/cut_tail/cut_right. This will have no effect." << endl;
    }

    // quality filtering
    opt->qualfilter.enabled = !cmd.exist("disable_quality_filtering");
    if(opt->nanopore_default){
        opt->qualfilter.qualifiedQual = num2qual(10);
    } else {
         opt->qualfilter.qualifiedQual = num2qual(cmd.get<int>("qualified_quality_phred"));
    }
   
    opt->qualfilter.unqualifiedPercentLimit = cmd.get<int>("unqualified_percent_limit");
    opt->qualfilter.avgQualReq = cmd.get<int>("average_qual");
    opt->qualfilter.nBaseLimit = cmd.get<int>("n_base_limit");

    // length filtering
    opt->lengthFilter.enabled = !cmd.exist("disable_length_filtering");
    opt->lengthFilter.requiredLength = cmd.get<int>("length_required");
    opt->lengthFilter.maxLength = cmd.get<int>("length_limit");

    // low complexity filter
    opt->complexityFilter.enabled = cmd.exist("low_complexity_filter");
    opt->complexityFilter.threshold = (min(100, max(0, cmd.get<int>("complexity_threshold")))) / 100.0;

    // overlap correction
    opt->correction.enabled = !cmd.exist("no_correction");
    opt->overlapRequire = cmd.get<int>("overlap_len_require");
    opt->overlapDiffLimit = cmd.get<int>("overlap_diff_limit");
    opt->overlapDiffPercentLimit = cmd.get<int>("overlap_diff_percent_limit");

    // threading
    opt->thread = cmd.get<int>("thread");

    // umi
    opt->umi.enabled = cmd.exist("umi");
    opt->umi.length = cmd.get<int>("umi_len");
    opt->umi.prefix = cmd.get<string>("umi_prefix");
    opt->umi.skip = cmd.get<int>("umi_skip");
    if(opt->umi.enabled) {
        string umiLoc = cmd.get<string>("umi_loc");
        str2lower(umiLoc);
        if(umiLoc.empty())
            error_exit("You've enabled UMI by (--umi), you should specify the UMI location by (--umi_loc)");
        if(umiLoc != "index1" && umiLoc != "index2" && umiLoc != "read1" && umiLoc != "read2" && umiLoc != "per_index" && umiLoc != "per_read") {
            error_exit("UMI location can only be index1/index2/read1/read2/per_index/per_read");
        }
        if(!opt->isPaired() && (umiLoc == "index2" || umiLoc == "read2"))
            error_exit("You specified the UMI location as " + umiLoc + ", but the input data is not paired end.");
        if(opt->umi.length == 0 && (umiLoc == "read1" || umiLoc == "read2" ||  umiLoc == "per_read"))
            error_exit("You specified the UMI location as " + umiLoc + ", but the length is not specified (--umi_len).");
        if(umiLoc == "index1") {
            opt->umi.location = UMI_LOC_INDEX1;
        } else if(umiLoc == "index2") {
            opt->umi.location = UMI_LOC_INDEX2;
        } else if(umiLoc == "read1") {
            opt->umi.location = UMI_LOC_READ1;
        } else if(umiLoc == "read2") {
            opt->umi.location = UMI_LOC_READ2;
        } else if(umiLoc == "per_index") {
            opt->umi.location = UMI_LOC_PER_INDEX;
        } else if(umiLoc == "per_read") {
            opt->umi.location = UMI_LOC_PER_READ;
        }
    }

    stringstream ss;
    for(int i=0;i<argc;i++){
        ss << argv[i] << " ";
    }
    command = ss.str();
    bool supportEvaluation = !opt->inputFromSTDIN && opt->in1!="/dev/stdin";
    time_t t1 = time(NULL);
    //opt->sampleTable = cmd.get<string>("sampleTable");
    if (opt->sampleTable.empty()) {
        opt->in1 = cmd.get<string>("in1");
        opt->in2 = cmd.get<string>("in2");
        opt->out1 = cmd.get<string>("out1");
        opt->out2 = cmd.get<string>("out2");
        opt->prefix = cmd.get<string>("prefix");
        opt->outFRFile = opt->prefix + "_output_reads.fastq.gz";
        opt->jsonFile = cmd.get<string>("json");
        opt->htmlFile = cmd.get<string>("html");
        opt->reportTitle = cmd.get<string>("report_title");
        
        Evaluator eva(opt);
        if (supportEvaluation) {
            eva.evaluateSeqLen();
        }

        long readNum = 0;

        // using evaluator to guess how many reads in total
        if (opt->shallDetectAdapter(false)) {
            if (!supportEvaluation) {
                std::cerr << "Adapter auto-detection is disabled for STDIN mode" << endl;
            } else {
                std::cerr << "Detecting adapter sequence for read1..." << endl;
                string adapt = eva.evalAdapterAndReadNum(readNum, false);
                if (adapt.length() > 60)
                    adapt.resize(0, 60);
                if (adapt.length() > 0) {
                    opt->adapter.sequence = adapt;
                    opt->adapter.detectedAdapter1 = adapt;
                } else {
                    std::cerr << "No adapter detected for read1" << endl;
                    opt->adapter.sequence = "";
                }
                std::cerr << endl;
            }
        }
        if (opt->shallDetectAdapter(true)) {
            if (!supportEvaluation) {
                std::cerr << "Adapter auto-detection is disabled for STDIN mode" << endl;
            } else {
                std::cerr << "Detecting adapter sequence for read2..." << endl;
                string adapt = eva.evalAdapterAndReadNum(readNum, true);
                if (adapt.length() > 60)
                    adapt.resize(0, 60);
                if (adapt.length() > 0) {
                    opt->adapter.sequenceR2 = adapt;
                    opt->adapter.detectedAdapter2 = adapt;
                } else {
                    std::cerr << "No adapter detected for read2" << endl;
                    opt->adapter.sequenceR2 = "";
                }
                std::cerr << endl;
            }
        }
        opt->validate();
        if(opt->debug) cCout(opt->prefix + " " + opt->in1 + " " + opt->in2);
        // using evaluator to check if it's two color system
        if (!cmd.exist("disable_trim_poly_g") && supportEvaluation) {
            bool twoColorSystem = eva.isTwoColorSystem();
            if (twoColorSystem) {
                opt->polyGTrim.enabled = true;
            }
        }
        opt->readLocFile();
        if(!opt->sexFile.empty()){
             opt->readSexLoc();
             if(opt->sexMap.size() > 1){
                 if (opt->verbose) loginfo("reading sex markers:" + std::to_string(opt->sexMap.size()));
                 std::string sexMarker = Evaluator::getSexMarker(opt);
                 if(sexMarker.empty()){
                     opt->sexFile.clear();
                     std::cerr << "No valid sex marker is detected!" << endl;
                 } else {
                     opt->mSex = opt->sexMap[sexMarker];
                     if(opt->verbose) {
                         std::cerr << "Sex marker " << sexMarker << " is detected!" << endl;
                     }
                 }
             }
        }

        Processor* p = new Processor(opt);
        p->process();

        if (p) {
            delete p;
            p = NULL;
        }

        time_t t2 = time(NULL);

        std::cerr << endl << "JSON report: " << opt->jsonFile << endl;
        std::cerr << "HTML report: " << opt->htmlFile << endl;
        std::cerr << endl << command << endl;
        std::cerr << "seq2sat v" << SEQ2SAT_VER << ", time used: " << (t2) - t1 << " seconds" << endl;

    } else {
        opt->parseSampleTable();
        for (auto & it : opt->samples) {
            time_t t2 = time(NULL);
            opt->prefix = it.prefix;
            opt->in1 = it.in1;
            opt->in2 = it.in2;
            opt->outFRFile = opt->prefix + "_output_reads.fastq.gz";
            //opt->out1 = it.prefix + "_R1.fastq.gz";
           // opt->out2 = it.prefix + "_R2.fastq.gz";
            opt->jsonFile = it.prefix + ".json";
            opt->htmlFile = it.prefix + ".html";
            opt->reportTitle = it.prefix;
            if(opt->verbose) cCout("Processing sample: " + basename(opt->prefix));
            
            Evaluator eva(opt);
            if (supportEvaluation) {
                eva.evaluateSeqLen();
            }

            long readNum = 0;

            // using evaluator to guess how many reads in total
            if (opt->shallDetectAdapter(false)) {
                if (!supportEvaluation) {
                    //std::cerr << "Adapter auto-detection is disabled for STDIN mode" << endl;
                } else {
                    //std::cerr << "Detecting adapter sequence for read1..." << endl;
                    string adapt = eva.evalAdapterAndReadNum(readNum, false);
                    if (adapt.length() > 60)
                        adapt.resize(0, 60);
                    if (adapt.length() > 0) {
                        opt->adapter.sequence = adapt;
                        opt->adapter.detectedAdapter1 = adapt;
                    } else {
                        //std::cerr << "No adapter detected for read1" << endl;
                        opt->adapter.sequence = "";
                    }
                    std::cerr << endl;
                }
            }
            if (opt->shallDetectAdapter(true)) {
                if (!supportEvaluation) {
                    //std::cerr << "Adapter auto-detection is disabled for STDIN mode" << endl;
                } else {
                    //std::cerr << "Detecting adapter sequence for read2..." << endl;
                    string adapt = eva.evalAdapterAndReadNum(readNum, true);
                    if (adapt.length() > 60)
                        adapt.resize(0, 60);
                    if (adapt.length() > 0) {
                        opt->adapter.sequenceR2 = adapt;
                        opt->adapter.detectedAdapter2 = adapt;
                    } else {
                        //std::cerr << "No adapter detected for read2" << endl;
                        opt->adapter.sequenceR2 = "";
                    }
                    std::cerr << endl;
                }
            }

            opt->validate();

            if(opt->debug) cCout(opt->prefix + " " + opt->in1 + " " + opt->in2);
            // using evaluator to check if it's two color system
            if (!cmd.exist("disable_trim_poly_g") && supportEvaluation) {
                bool twoColorSystem = eva.isTwoColorSystem();
                if (twoColorSystem) {
                    opt->polyGTrim.enabled = true;
                }
            }

            opt->readLocFile();
            opt->readSexLoc();
            
            Processor p(opt);
            p.process();
            
            time_t t3 = time(NULL);

            std::cerr << endl << "JSON report: " << opt->jsonFile << endl;
            std::cerr << "HTML report: " << opt->htmlFile << endl;
            std::cerr << endl << command << endl;
            std::cerr << "seq2sat v" << SEQ2SAT_VER << ", time used: " << (t3) - t2 << " seconds" << endl;
        }

        //HtmlReporterAll hra(opt);
        //hra.report();
    }

    if(opt){
        delete opt;
        opt = NULL;
    }
    time_t t4 = time(NULL);
    std::cerr << "seq2sat v" << SEQ2SAT_VER << ", time used: " << (t4) - t1 << " seconds" << endl;
    return 0;
}
