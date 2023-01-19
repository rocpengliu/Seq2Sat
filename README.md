# Seq2Sat


## Description 

**Seq2Sat** is a novel, multifunctional, ultra-fast and all-in-one software to score microsatellite/SSR from target sequencing.

## Key features of Seq2Sat

* **Multifunctional**: Seq2Sat can conduct genotyping for both microsatellite/SSR and SNPs from both target sequencing and GBS/RAD-seq.


* **Ultra-fast**: ~ 5 second/sample using single thread.

* **Extremely low memory cost**: negligible ~ 0.1 MB RAM.

* **All-in-one**: Seq2Sat directly takes raw target sequencing reads as input and output genotype table without any intermediate file writing and loading, making I/O very efficient.

* **Easy to use**: requires minimal programing skills.


## Getting started
### Step 1. Install the software
Seq2Sat  is written in C/C++11 and can be installed on Linux or Mac OS X (with Xcode and Xcode Command Line Tools installed). 
We have tested Seq2Sat on Ubuntu (16.04 LTS and above) and macOS Catalina.

```
git clone 
cd Seq2Sat
make
```

### Step 2. Run a small test 





### Running Seq2Sat
#### 1.Preparing sample.txt file

This file consists of 3 columns and separated by '\t'. The first column is the prefix name of each sample, and the second is the forward reads file, the thrid column is the reverse reads file. If you have a single-end (SE) reads, remove the reverse reads (third) column. It looks like this:
```
A1.CE2-S1	A1.CE2-S1_R1.fastq.gz	A1.CE2-S1_R1.fastq.gz
A3.CE1-M2	A3.CE1-M2_R1.fastq.gz	A3.CE1-M2_R1.fastq.gz
A4.CE1-H5	A4.CE1-H5_R1.fastq.gz	A4.CE1-H5_R1.fastq.gz
B1.CE2-S2	B1.CE2-S2_R1.fastq.gz	B1.CE2-S2_R1.fastq.gz
B3.CE1-M3	B3.CE1-M3_R1.fastq.gz	B3.CE1-M3_R1.fastq.gz
C1.CE2-S3	C1.CE2-S3_R1.fastq.gz	C1.CE2-S3_R1.fastq.gz
C3.CE1-M4	C3.CE1-M4_R1.fastq.gz	C3.CE1-M4_R1.fastq.gz
```

## Seq2Sat full usage options
```
Use seq2fun or seq2fun --help to show the full usage options

  options:

 var, genetic variance, must be either microsatellite/ssr or snp
 sampleTable, 's', "Sample table consisting of sample prefix/name, forward sequence file, reverse sequence file if paired-end, and must be separated by '\t', it is highly recommended if you have multiple samples
 in1, 'i', "read1 input file name
 in2, 'I', "read2 input file name
 prefix, 'X', "prefix name for output files, eg: sample01
 out1, 'o', "file name to store read1 with on-target sequences
 out2, 'O', "file name to store read2 with on-target sequences
    
  loc, 0, "loci file containing loci names, 5'primer sequence, reverse complement of 3'primer sequence, 5'flank region, 3'flank region, repeat unit and reference microsatellite repeat array, separated by '\t"
    maxMismatchesPSeq", 0, "maximum mismatches for primer sequences 2
    minSeqs", 0, "minimum number of reads for a genotype, default: 10
    minSeqsPercentage", 0, "minimum percentage (%) reads against largest peak for a genotype, default: 10 (10%)
    mode", 0, "specify the sequence alignment mode: NW (default) | HW | SHW
    maxScore", 0, "specify the maximum score of sequence alignment with sore > maxScore will be discarded, default value is -1, and no sequence will be discarded
    numBestSeqs", 0, "Score will be calculated only for N best sequences (best = with smallest score). If N = 0 then all sequences will be calculated
    notFindAlignment", 0, "If specified, alignment path will be not found and printed. This may significantly speed up the calculation
    findStartLocation", 0, "If specified, start locations will be found and printed. Each start location corresponds to one end location. This may somewhat slow down the calculation, but is still faster then finding alignment path and does not consume any extra memory

    format", 0, "NICE|CIG_STD|CIG_EXT  Format that will be used to print alignment path, can be used only with -p. NICE will give visually attractive format, CIG_STD will give standard cigar format and CIG_EXT will give extended cigar format. [default: NICE]
    core", 0, "Core part of calculation will be repeated N times. This is useful only for performance measurement, when single execution is too short to measure
    silentAlignment", 0, "If specified, there will be no score or alignment output
    printResults", 0, "If specified, alignment results will be printed but with only 1 thread
    debug", 0, "If specified, print debug
    
    dont_merge_overlapped_PE", 0, "don't merge the overlapped PE reads; this is off by default
    
    // reporting
    json", 'j', "the json format report file name
    html", 'h', "the html format report file name
    report_title", 'R', "should be quoted with \' or \", default is \"seq2sat report\"

    // threading
    thread", 'w', "worker thread number, default is 4

    // qother I/O
    phred64", '6', "indicate the input is using phred64 scoring (it'll be converted to phred33, so the output will still be phred33
    compression", 'z', "compression level for gzip output (1 ~ 9). 1 is fastest, 9 is smallest, default is 4
    stdin", 0, "input from STDIN. If the STDIN is interleaved paired-end FASTQ, please also add --interleaved_in
    stdout", 0, "stream passing-filters reads to STDOUT. This option will result in interleaved FASTQ output for paired-end output. Disabled by default
    interleaved_in", 0, "indicate that <in1> is an interleaved FASTQ which contains both read1 and read2. Disabled by default
    reads_to_process", 0, "specify how many reads/pairs to be processed. Default 0 means process all reads.
    dont_overwrite", 0, "don't overwrite existing files. Overwritting is allowed by default
    verbose", 'V', "output verbose log information (i.e. when every 1M reads are processed).

    // adapter
    disable_adapter_trimming", 'A', "adapter trimming is enabled by default. If this option is specified, adapter trimming is disabled
    adapter_sequence", 'a', "the adapter for read1. For SE data, if not specified, the adapter will be auto-detected. For PE data, this is used if R1/R2 are found not overlapped.
    adapter_sequence_r2", 0, "the adapter for read2 (PE data only). This is used if R1/R2 are found not overlapped. If not specified, it will be the same as <adapter_sequence>
    adapter_fasta", 0, "specify a FASTA file to trim both read1 and read2 (if PE) by all the sequences in this FASTA file
    detect_adapter_for_pe", 0, "by default, the auto-detection for adapter is for SE data input only, turn on this option to enable it for PE data

    // trimming
    trim_front1", 'f', "trimming how many bases in front for read1, default is 0
    trim_tail1", 't', "trimming how many bases in tail for read1, default is 0
    max_len1", 'b', "if read1 is longer than max_len1, then trim read1 at its tail to make it as long as max_len1. Default 0 means no limitation
    trim_front2", 'F', "trimming how many bases in front for read2. If it's not specified, it will follow read1's settings
    trim_tail2", 'T', "trimming how many bases in tail for read2. If it's not specified, it will follow read1's settings
    max_len2", 'B', "if read2 is longer than max_len2, then trim read2 at its tail to make it as long as max_len2. Default 0 means no limitation. If it's not specified, it will follow read1's settings

    // polyG tail trimming
    poly_g_min_len", 0, "the minimum length to detect polyG in the read tail. 10 by default
    disable_trim_poly_g", 'G', "disable polyG tail trimming, by default trimming is automatically enabled for Illumina NextSeq/NovaSeq data
    
    // polyX tail trimming
    trim_poly_x", 'x', "enable polyX trimming in 3' ends.
    poly_x_min_len", 0, "the minimum length to detect polyX in the read tail. 10 by default.

    // cutting by quality
    cut_front", '5', "move a sliding window from front (5') to tail, drop the bases in the window if its mean quality < threshold, stop otherwise.
    cut_tail", '3', "move a sliding window from tail (3') to front, drop the bases in the window if its mean quality < threshold, stop otherwise.
    cut_right", 'r', "move a sliding window from front to tail, if meet one window with mean quality < threshold, drop the bases in the window and the right part, and then stop.
    cut_window_size", 'W', "the window size option shared by cut_front, cut_tail or cut_sliding. Range: 1~1000, default: 4
    cut_mean_quality", 'M', "the mean quality requirement option shared by cut_front, cut_tail or cut_sliding. Range: 1~36 default: 20 (Q20)
    
    cut_front_window_size", 0, "the window size option of cut_front, default to cut_window_size if not specified
    cut_front_mean_quality", 0, "the mean quality requirement option for cut_front, default to cut_mean_quality if not specified
    cut_tail_window_size", 0, "the window size option of cut_tail, default to cut_window_size if not specified"
    cut_tail_mean_quality", 0, "the mean quality requirement option for cut_tail, default to cut_mean_quality if not specified
    cut_right_window_size", 0, "the window size option of cut_right, default to cut_window_size if not specified
    cut_right_mean_quality", 0, "the mean quality requirement option for cut_right, default to cut_mean_quality if not specified

    // quality filtering
    disable_quality_filtering", 'Q', "quality filtering is enabled by default. If this option is specified, quality filtering is disabled");
    qualified_quality_phred", 'q', "the quality value that a base is qualified. Default 20 means phred quality >=Q15 is qualified.
    unqualified_percent_limit", 'u', "how many percents of bases are allowed to be unqualified (0~100). Default 40 means 40%
    n_base_limit", 'n', "if one read's number of N base is >n_base_limit, then this read/pair is discarded. Default is 5
    average_qual", 'e', "if one read's average quality score <avg_qual, then this read/pair is discarded. Default 0 means no requirement

    // length filtering
    cmd.add("disable_length_filtering", 'L', "length filtering is enabled by default. If this option is specified, length filtering is disabled");
    length_required", 'l', "reads shorter than length_required will be discarded, default is 50.
    length_limit", 0, "reads longer than length_limit will be discarded, default 0 means no limitation.

    // low complexity filtering
    low_complexity_filter", 'y', "enable low complexity filter. The complexity is defined as the percentage of base that is different from its next base (base[i] != base[i+1]).
    complexity_threshold", 'Y', "the threshold for low complexity filter (0~100). Default is 30, which means 30% complexity is required.

    // filter by indexes
    filter_by_index1", 0, "specify a file contains a list of barcodes of index1 to be filtered out, one barcode per line
    filter_by_index2", 0, "specify a file contains a list of barcodes of index2 to be filtered out, one barcode per line
    filter_by_index_threshold", 0, "the allowed difference of index barcode for index filtering, default 0 means completely identical.
    
    // base correction in overlapped regions of paired end data
    no_correction", 'C', "disable base correction in overlapped regions (only for PE data), default is enabled
    overlap_len_require", 0, "the minimum length to detect overlapped region of PE reads. This will affect overlap analysis based PE merge, adapter trimming and correction. 30 by default.
    overlap_diff_limit", 0, "the maximum number of mismatched bases to detect overlapped region of PE reads. This will affect overlap analysis based PE merge, adapter trimming and correction. 5 by default.
    overlap_diff_percent_limit", 0, "the maximum percentage of mismatched bases to detect overlapped region of PE reads. This will affect overlap analysis based PE merge, adapter trimming and correction. Default 20 means 20%.

    // umi
    umi", 'U', "enable unique molecular identifier (UMI) preprocessing
    umi_loc", 0, "specify the location of UMI, can be (index1/index2/read1/read2/per_index/per_read, default is none
    umi_len", 0, "if the UMI is in read1/read2, its length should be provided
    umi_prefix", 0, "if specified, an underline will be used to connect prefix and UMI (i.e. prefix=UMI, UMI=AATTCG, final=UMI_AATTCG). No prefix by default
    umi_skip", 0, "if the UMI is in read1/read2, seq2sat can skip several bases following UMI, default is 0

```


## Bugs or feature requests

To inform us of any bugs or requests, please open a new issue or send an email to peng.liu@ec.gc.ca

## Seq2Sate History & Updates
07/27/2022 - added functions to generate html report
05/03/2022 - first realse
