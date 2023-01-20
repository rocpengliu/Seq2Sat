# Seq2Sat


## Description 

**Seq2Sat** is a novel, multifunctional, ultra-fast and all-in-one software to score microsatellite/SSR from target sequencing.

## Key features of Seq2Sat

* **Multifunctional**: Seq2Sat can conduct genotyping for both microsatellite/SSR and SNPs from both target sequencing and GBS/RAD-seq.


* **Ultra-fast**: ~ 1 second/sample using single thread.

* **Extremely low memory cost**: negligible ~ 0.1 MB RAM.

* **All-in-one**: Seq2Sat directly takes raw target sequencing reads as input and output genotype table without any intermediate file writing and loading, making I/O very efficient.

* **Easy to use**: requires minimal programing skills.


## Getting started
### Step 1. Install the software
Seq2Sat  is written in C/C++11 and can be installed on Linux or Mac OS X (with Xcode and Xcode Command Line Tools installed). 
We have tested Seq2Sat on Ubuntu (16.04 LTS and above) and macOS Catalina.

```
git clone https://github.com/ecogenomicscanada/Seq2Sat
cd Seq2Sat
make
```

### Step 2. Run a small test 

```
cd testdata;
cat sample.txt| while read i j k l; do ../seq2sat --prefix ${i} -i ${j} -I ${k} --loc loc.txt --sex sexLoc.txt --var ssr -w 8 -V; done;
```

### View result
```
Please view the sex identification, genotypes and reads quality in a [html report](https://github.com/ecogenomicscanada/Seq2Sat/tree/master/res/) for one sample.
```

## Seq2Sat full usage options
```
Use seq2sat or seq2sat --help to show the full usage options
options:
      --var                           genetic variance, must be either microsatellite/ssr or snp (string [=])
  -i, --in1                           read1 input file name (string [=])
  -I, --in2                           read2 input file name (string [=])
  -X, --prefix                        prefix name for output files, eg: sample01 (string [=])
      --outFReads                     If specified, off-target reads will be outputed in a file
  -o, --out1                          file name to store read1 with on-target sequences (string [=])
  -O, --out2                          file name to store read2 with on-target sequences (string [=])
      --loc                           loci file containing loci names, 5'primer sequence, reverse complement of 3'primer sequence, 5'flank region, 3'flank region, repeat unit and reference microsatellite repeat array, separated by '	 (string [=])
      --maxMismatchesPSeq             maximum mismatches for primer sequences 2 (int [=2])
      --maxMismatchesPer4FR           maximum percentage mismatches for the forward and reverse flanking regions, default 0.3 (30%)  (double [=0.1])
      --minSeqs                       minimum number of reads for a genotype, default: 5 (int [=5])
      --minWarningSeqs                minimum number of reads for warning a genotype, default: 50 (int [=50])
      --minSeqsPercentage             minimum percentage (%) reads against largest peak for a genotype, default: 5 (5%) (int [=5])
      --hlRatio1                      ratio of loci sizes of largest and second largest numbers of reads when the length difference = 1 ssr unit, default: 0.4 (double [=0.4])
      --hlRatio2                      ratio of loci sizes of largest and second largest numbers of reads when the length difference = 2 ssr unit, default: 0.2 (double [=0.2])
      --maxVarRatio                   ratio of two heter alleles based on variations either in flanking regions or MRA, the ideal is 1, default: 1.5 (double [=1.5])
      --mode                          specify the sequence alignment mode: NW (default) | HW | SHW (string [=NW])
      --maxScore                      specify the maximum score of sequence alignment with sore > maxScore will be discarded, default value is -1, and no sequence will be discarded. (int [=-1])
      --numBestSeqs                   Score will be calculated only for N best sequences (best = with smallest score). If N = 0 then all sequences will be calculated. (int [=0])
      --notFindAlignment              If specified, alignment path will be not found and printed. This may significantly speed up the calculation
      --findStartLocation             If specified, start locations will be found and printed. Each start location corresponds to one end location. This may somewhat slow down the calculation, but is still faster then finding alignment path and does not consume any extra memory.
      --format                        NICE|CIG_STD|CIG_EXT  Format that will be used to print alignment path, can be used only with -p. NICE will give visually attractive format, CIG_STD will give standard cigar format and CIG_EXT will give extended cigar format. [default: NICE] (string [=NICE])
      --core                          Core part of calculation will be repeated N times. This is useful only for performance measurement, when single execution is too short to measure. (int [=1])
      --silentAlignment               If specified, there will be no score or alignment output
      --printResults                  If specified, alignment results will be printed but with only 1 thread
      --sex                           sex loci file containing sex locus names, 5'primer sequence, reverse complement of 3'primer sequence, X/Z reference sequence, Y/W reference sequence, separated by ' (string [=])
      --maxMismatchesSexPSeq          maximum number of mismatches for sex primers, default: 2 (unsigned int [=2])
      --maxMismatchesSexRefSeq        maximum number of mismatches for sex reference sequences, default: 2 (unsigned int [=2])
      --yxRatio                       minimum ratio of numbers of reads Y/X to W/Z, default: 0.001 (double [=0.001])
      --minTotalReadsX                minimum number of reads assigned to X; default: 10 (int [=10])
      --minTotalReadsY                minimum number of reads assigned to Y; default: 10 (int [=10])
      --minReadsX                     minimum number of reads assigned to each variant of X; default: 5 (int [=5])
      --minReadsY                     minimum number of reads assigned to each variant of Y; default: 5 (int [=5])
      --debug                         If specified, print debug
      --dont_merge_overlapped_PE      don't merge the overlapped PE reads; this is off by default
  -j, --json                          the json format report file name (string [=seq2sat.json])
  -h, --html                          the html format report file name (string [=seq2sat.html])
  -R, --report_title                  should be quoted with ' or ", default is "seq2sat report" (string [=seq2sat report])
  -w, --thread                        worker thread number, default is 4 (int [=4])
  -6, --phred64                       indicate the input is using phred64 scoring (it'll be converted to phred33, so the output will still be phred33)
  -z, --compression                   compression level for gzip output (1 ~ 9). 1 is fastest, 9 is smallest, default is 4. (int [=4])
      --stdin                         input from STDIN. If the STDIN is interleaved paired-end FASTQ, please also add --interleaved_in.
      --stdout                        stream passing-filters reads to STDOUT. This option will result in interleaved FASTQ output for paired-end output. Disabled by default.
      --interleaved_in                indicate that <in1> is an interleaved FASTQ which contains both read1 and read2. Disabled by default.
      --reads_to_process              specify how many reads/pairs to be processed. Default 0 means process all reads. (int [=0])
      --dont_overwrite                don't overwrite existing files. Overwritting is allowed by default.
  -V, --verbose                       output verbose log information (i.e. when every 1M reads are processed).
  -A, --disable_adapter_trimming      adapter trimming is enabled by default. If this option is specified, adapter trimming is disabled
  -a, --adapter_sequence              the adapter for read1. For SE data, if not specified, the adapter will be auto-detected. For PE data, this is used if R1/R2 are found not overlapped. (string [=auto])
      --adapter_sequence_r2           the adapter for read2 (PE data only). This is used if R1/R2 are found not overlapped. If not specified, it will be the same as <adapter_sequence> (string [=auto])
      --adapter_fasta                 specify a FASTA file to trim both read1 and read2 (if PE) by all the sequences in this FASTA file (string [=])
      --detect_adapter_for_pe         by default, the auto-detection for adapter is for SE data input only, turn on this option to enable it for PE data.
  -f, --trim_front1                   trimming how many bases in front for read1, default is 0 (int [=0])
  -t, --trim_tail1                    trimming how many bases in tail for read1, default is 0 (int [=0])
  -b, --max_len1                      if read1 is longer than max_len1, then trim read1 at its tail to make it as long as max_len1. Default 0 means no limitation (int [=0])
  -F, --trim_front2                   trimming how many bases in front for read2. If it's not specified, it will follow read1's settings (int [=0])
  -T, --trim_tail2                    trimming how many bases in tail for read2. If it's not specified, it will follow read1's settings (int [=0])
  -B, --max_len2                      if read2 is longer than max_len2, then trim read2 at its tail to make it as long as max_len2. Default 0 means no limitation. If it's not specified, it will follow read1's settings (int [=0])
      --poly_g_min_len                the minimum length to detect polyG in the read tail. 10 by default. (int [=10])
  -G, --disable_trim_poly_g           disable polyG tail trimming, by default trimming is automatically enabled for Illumina NextSeq/NovaSeq data
  -x, --trim_poly_x                   enable polyX trimming in 3' ends.
      --poly_x_min_len                the minimum length to detect polyX in the read tail. 10 by default. (int [=10])
  -5, --cut_front                     move a sliding window from front (5') to tail, drop the bases in the window if its mean quality < threshold, stop otherwise.
  -3, --cut_tail                      move a sliding window from tail (3') to front, drop the bases in the window if its mean quality < threshold, stop otherwise.
  -r, --cut_right                     move a sliding window from front to tail, if meet one window with mean quality < threshold, drop the bases in the window and the right part, and then stop.
  -W, --cut_window_size               the window size option shared by cut_front, cut_tail or cut_sliding. Range: 1~1000, default: 4 (int [=4])
  -M, --cut_mean_quality              the mean quality requirement option shared by cut_front, cut_tail or cut_sliding. Range: 1~36 default: 20 (Q20) (int [=20])
      --cut_front_window_size         the window size option of cut_front, default to cut_window_size if not specified (int [=4])
      --cut_front_mean_quality        the mean quality requirement option for cut_front, default to cut_mean_quality if not specified (int [=20])
      --cut_tail_window_size          the window size option of cut_tail, default to cut_window_size if not specified (int [=4])
      --cut_tail_mean_quality         the mean quality requirement option for cut_tail, default to cut_mean_quality if not specified (int [=20])
      --cut_right_window_size         the window size option of cut_right, default to cut_window_size if not specified (int [=4])
      --cut_right_mean_quality        the mean quality requirement option for cut_right, default to cut_mean_quality if not specified (int [=20])
  -Q, --disable_quality_filtering     quality filtering is enabled by default. If this option is specified, quality filtering is disabled
  -q, --qualified_quality_phred       the quality value that a base is qualified. Default 20 means phred quality >=Q15 is qualified. (int [=20])
  -u, --unqualified_percent_limit     how many percents of bases are allowed to be unqualified (0~100). Default 40 means 40% (int [=40])
  -n, --n_base_limit                  if one read's number of N base is >n_base_limit, then this read/pair is discarded. Default is 5 (int [=5])
  -e, --average_qual                  if one read's average quality score <avg_qual, then this read/pair is discarded. Default 0 means no requirement (int [=0])
  -L, --disable_length_filtering      length filtering is enabled by default. If this option is specified, length filtering is disabled
  -l, --length_required               reads shorter than length_required will be discarded, default is 30. (int [=30])
      --length_limit                  reads longer than length_limit will be discarded, default 0 means no limitation. (int [=0])
  -y, --low_complexity_filter         enable low complexity filter. The complexity is defined as the percentage of base that is different from its next base (base[i] != base[i+1]).
  -Y, --complexity_threshold          the threshold for low complexity filter (0~100). Default is 30, which means 30% complexity is required. (int [=30])
      --filter_by_index1              specify a file contains a list of barcodes of index1 to be filtered out, one barcode per line (string [=])
      --filter_by_index2              specify a file contains a list of barcodes of index2 to be filtered out, one barcode per line (string [=])
      --filter_by_index_threshold     the allowed difference of index barcode for index filtering, default 0 means completely identical. (int [=0])
  -C, --no_correction                 disable base correction in overlapped regions (only for PE data), default is enabled
      --overlap_len_require           the minimum length to detect overlapped region of PE reads. This will affect overlap analysis based PE merge, adapter trimming and correction. 30 by default. (int [=30])
      --overlap_diff_limit            the maximum number of mismatched bases to detect overlapped region of PE reads. This will affect overlap analysis based PE merge, adapter trimming and correction. 5 by default. (int [=5])
      --overlap_diff_percent_limit    the maximum percentage of mismatched bases to detect overlapped region of PE reads. This will affect overlap analysis based PE merge, adapter trimming and correction. Default 20 means 20%. (int [=20])
  -U, --umi                           enable unique molecular identifier (UMI) preprocessing
      --umi_loc                       specify the location of UMI, can be (index1/index2/read1/read2/per_index/per_read, default is none (string [=])
      --umi_len                       if the UMI is in read1/read2, its length should be provided (int [=0])
      --umi_prefix                    if specified, an underline will be used to connect prefix and UMI (i.e. prefix=UMI, UMI=AATTCG, final=UMI_AATTCG). No prefix by default (string [=])
      --umi_skip                      if the UMI is in read1/read2, seq2sat can skip several bases following UMI, default is 0 (int [=0])
  -?, --help                          print this message


```


## Bugs or feature requests

To inform us of any bugs or requests, please open a new issue or send an email to peng.liu@ec.gc.ca

## Seq2Sate History & Updates
12/20/2022 - first version released  
10/11/2022 - added functions to generate 3d plot  
08/10/2022 - added functions for sex identification   
07/27/2022 - added functions to generate html report  
05/03/2022 - first realse
