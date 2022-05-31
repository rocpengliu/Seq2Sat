#include "htmlreporter.h"
#include <chrono>
#include <memory.h>
#include <valarray>

extern string command;

HtmlReporter::HtmlReporter(Options* opt){
    mOptions = opt;
    mDupHist = NULL;
    mDupRate = 0.0;
}

HtmlReporter::~HtmlReporter(){
}

void HtmlReporter::setDupHist(int* dupHist, double* dupMeanGC, double dupRate) {
    mDupHist = dupHist;
    mDupMeanGC = dupMeanGC;
    mDupRate = dupRate;
}

void HtmlReporter::setInsertHist(long* insertHist, int insertSizePeak) {
    mInsertHist = insertHist;
    mInsertSizePeak = insertSizePeak;
}

void HtmlReporter::outputRow(ofstream& ofs, string key, long v) {
    ofs << "<tr><td class='col1'>" + key + "</td><td class='col2'>" + to_string(v) + "</td></tr>\n";
}

void HtmlReporter::outputRow(ofstream& ofs, string key, string v) {
    ofs << "<tr><td class='col1'>" + key + "</td><td class='col2'>" + v + "</td></tr>\n";
}

void HtmlReporter::outputRow(ofstream& ofs, std::string & marker, std::vector<std::pair<std::string, Genotype>> & outGenotype) {
    for (auto & it : outGenotype) {
        ofs << "<tr>";
        ofs << "<td>" + marker + "</td>" +
                "<td>" + it.second.baseLocVar.repuitAll.mStr + "</td>" +
                "<td>" + it.second.baseLocVar.mraBase + "</td>" +
                "<td>" + std::to_string(it.second.baseLocVar.mra.mStr.length()) + "</td>" +
                "<td>" + std::to_string(it.second.baseLocVar.effectiveLen) + "</td>" +
                "<td>" + std::to_string(it.second.numReads) + "</td>" +
                "<td align='right'>" + highligher(it.second.baseLocVar.ff.mStr, it.second.baseLocVar.snpsMapff) + "</td>" +
                "<td align='center'>" + it.second.baseLocVar.mraName + "</td>" +
                "<td align='left'>" + highligher(it.second.baseLocVar.rf.mStr, it.second.baseLocVar.snpsMaprf) + "</td>";
        ofs << "</tr>\n";
    }
}

void HtmlReporter::outputRow(ofstream& ofs, std::string & marker, std::map<std::string, LocSnp> & snpsMap) {
    for (auto & it : snpsMap) {
        ofs << "<tr>";
        ofs << "<td>" + marker + "</td>" +
                "<td>" + it.second.getGenotype() + "</td>" +
                "<td>" + std::to_string(it.second.numReads) + "</td>" +
                "<td>" + std::to_string(it.second.ref.mStr.length()) + "</td>" +
                "<td align='center'>" + highligher(it.second) + "</td>";
        ofs << "</tr>\n";
    }
}

string HtmlReporter::formatNumber(long number) {
    double num = (double)number;
    string unit[6] = {"", "K", "M", "G", "T", "P"};
    int order = 0;
    while (num > 1000.0) {
        order += 1;
        num /= 1000.0;
    }

    if (order == 0)
        return to_string(number);
    else
        return to_string(num) + " " + unit[order];
}

string HtmlReporter::getPercents(long numerator, long denominator) {
    if(denominator == 0)
        return "0.0";
    else
        return to_string((double)numerator * 100.0 / (double)denominator);
}

void HtmlReporter::printSummary(ofstream& ofs, FilterResult* result, Stats* preStats1, Stats* postStats1, Stats* preStats2, Stats* postStats2) {
    long pre_total_reads = preStats1->getReads();
    if(preStats2)
        pre_total_reads += preStats2->getReads();

    long pre_total_bases = preStats1->getBases();
    if(preStats2)
        pre_total_bases += preStats2->getBases();

    long pre_q20_bases = preStats1->getQ20();
    if(preStats2)
        pre_q20_bases += preStats2->getQ20();

    long pre_q30_bases = preStats1->getQ30();
    if(preStats2)
        pre_q30_bases += preStats2->getQ30();

    long pre_total_gc = preStats1->getGCNumber();
    if(preStats2)
        pre_total_gc += preStats2->getGCNumber();

    long post_total_reads = postStats1->getReads();
    if(postStats2)
        post_total_reads += postStats2->getReads();

    long post_total_bases = postStats1->getBases();
    if(postStats2)
        post_total_bases += postStats2->getBases();

    long post_q20_bases = postStats1->getQ20();
    if(postStats2)
        post_q20_bases += postStats2->getQ20();

    long post_q30_bases = postStats1->getQ30();
    if(postStats2)
        post_q30_bases += postStats2->getQ30();

    long post_total_gc = postStats1->getGCNumber();
    if(postStats2)
        post_total_gc += postStats2->getGCNumber();

    string sequencingInfo  = mOptions->isPaired()?"paired end":"single end";
    if(mOptions->isPaired()) {
        sequencingInfo += " (" + to_string(preStats1->getCycles()) + " cycles + " + to_string(preStats2->getCycles()) + " cycles)";
    } else {
        sequencingInfo += " (" + to_string(preStats1->getCycles()) + " cycles)";
    }

    ofs << endl;
    ofs << "<div class='section_div'>\n";
    ofs << "<div class='section_title' onclick=showOrHide('summary')><a name='summary'>Data QC summary <font color='#88CCFF' > (click to show/hide) </font></a></div>\n";
    ofs << "<div id='summary'>\n";

    ofs << "<div class='subsection_title' onclick=showOrHide('general')>General</div>\n";
    ofs << "<div id='general'>\n";
    ofs << "<table class='summary_table'>\n";
    outputRow(ofs, "sequencing:", sequencingInfo);

    // report read length change
    if(mOptions->isPaired()) {
        outputRow(ofs, "mean length before filtering:", to_string(preStats1->getMeanLength()) + "bp, " + to_string(preStats2->getMeanLength()) + "bp");
    } else  {
        outputRow(ofs, "mean length before filtering:", to_string(preStats1->getMeanLength()) + "bp");
        outputRow(ofs, "mean length after filtering:", to_string(postStats1->getMeanLength()) + "bp");
    }

    if(mOptions->duplicate.enabled) {
        string dupStr = to_string(mDupRate*100) + "%";
        if(!mOptions->isPaired())
            dupStr += " (may be overestimated since this is SE data)";
        outputRow(ofs, "duplication rate:", dupStr);
    }
    if(mOptions->isPaired()) {
        outputRow(ofs, "Insert size peak:", mInsertSizePeak);
    }
    if(mOptions->adapterCuttingEnabled()) {
        if(!mOptions->adapter.detectedAdapter1.empty())
            outputRow(ofs, "Detected read1 adapter:", mOptions->adapter.detectedAdapter1);
        if(!mOptions->adapter.detectedAdapter2.empty())
            outputRow(ofs, "Detected read2 adapter:", mOptions->adapter.detectedAdapter2);
    }
    ofs << "</table>\n";
    ofs << "</div>\n";

    ofs << "<div class='subsection_title' onclick=showOrHide('before_filtering_summary')>Original data</div>\n";
    ofs << "<div id='before_filtering_summary'>\n";
    ofs << "<table class='summary_table'>\n";
    outputRow(ofs, "total reads:", formatNumber(pre_total_reads));
    outputRow(ofs, "total bases:", formatNumber(pre_total_bases));
    outputRow(ofs, "Q20 bases:", formatNumber(pre_q20_bases) + " (" + getPercents(pre_q20_bases,pre_total_bases) + "%)");
    outputRow(ofs, "Q30 bases:", formatNumber(pre_q30_bases) + " (" + getPercents(pre_q30_bases, pre_total_bases) + "%)");
    outputRow(ofs, "GC content:", getPercents(pre_total_gc,pre_total_bases) + "%");
    ofs << "</table>\n";
    ofs << "</div>\n";

    ofs << "<div class='subsection_title' onclick=showOrHide('after_filtering_summary')>Clean data used for detection</div>\n";
    ofs << "<div id='after_filtering_summary'>\n";
    ofs << "<table class='summary_table'>\n";
    outputRow(ofs, "total reads:", formatNumber(post_total_reads));
    outputRow(ofs, "total bases:", formatNumber(post_total_bases));
    outputRow(ofs, "Q20 bases:", formatNumber(post_q20_bases) + " (" + getPercents(post_q20_bases, post_total_bases) + "%)");
    outputRow(ofs, "Q30 bases:", formatNumber(post_q30_bases) + " (" + getPercents(post_q30_bases, post_total_bases) + "%)");
    outputRow(ofs, "GC content:", getPercents(post_total_gc,post_total_bases) + "%");
    ofs << "</table>\n";
    ofs << "</div>\n";

    if(result) {
        ofs << "<div class='subsection_title' onclick=showOrHide('filtering_result')>Filtering result</div>\n";
        ofs << "<div id='filtering_result'>\n";
        result -> reportHtml(ofs, pre_total_reads, pre_total_bases);
        ofs << "</div>\n";
    }

    ofs << "</div>\n";
    ofs << "</div>\n";

    if(result && mOptions->adapterCuttingEnabled()) {
        ofs << "<div class='section_div'>\n";
        ofs << "<div class='section_title' onclick=showOrHide('adapters')><a name='summary'>Adapters <font color='#88CCFF' > (click to show/hide) </font></a></div>\n";
        ofs << "<div id='adapters' style='display:none'>\n";

        result->reportAdapterHtml(ofs, pre_total_bases);

        ofs << "</div>\n";
        ofs << "</div>\n";
    }

    if(mOptions->duplicate.enabled) {
        ofs << "<div class='section_div'>\n";
        ofs << "<div class='section_title' onclick=showOrHide('duplication')><a name='summary'>Duplication <font color='#88CCFF' > (click to show/hide) </font></a></div>\n";
        ofs << "<div id='duplication' style='display:none'>\n";

        reportDuplication(ofs);

        ofs << "</div>\n";
        ofs << "</div>\n";
    }

    if(mOptions->isPaired()) {
        ofs << "<div class='section_div'>\n";
        ofs << "<div class='section_title' onclick=showOrHide('insert_size')><a name='summary'>Insert size estimation <font color='#88CCFF' > (click to show/hide) </font></a></div>\n";
        ofs << "<div id='insert_size' style='display:none'>\n";

        reportInsertSize(ofs, preStats1->getCycles() + preStats2->getCycles() - mOptions->overlapRequire);

        ofs << "</div>\n";
        ofs << "</div>\n";
    }
}

void HtmlReporter::reportInsertSize(ofstream& ofs, int isizeLimit) {
    if(isizeLimit<1)
        isizeLimit = 1;
    int total = min(mOptions->insertSizeMax, isizeLimit);
    long *x = new long[total];
    double allCount = 0;
    for(int i=0; i<total; i++) {
        x[i] = i;
        allCount += mInsertHist[i];
    }
    allCount += mInsertHist[mOptions->insertSizeMax];
    double* percents = new double[total];
    memset(percents, 0, sizeof(double)*total);
    if(allCount > 0) {
        for(int i=0; i<total; i++) {
            percents[i] = (double)mInsertHist[i] * 100.0 / (double)allCount;
        }
    }

    double unknownPercents = (double)mInsertHist[mOptions->insertSizeMax] * 100.0 / (double)allCount;

    ofs << "<div id='insert_size_figure'>\n";
    ofs << "<div class='figure' id='plot_insert_size' style='height:400px;'></div>\n";
    ofs << "</div>\n";

    ofs << "<div class='sub_section_tips'>This estimation is based on paired-end overlap analysis, and there are ";
    ofs << to_string(unknownPercents);
    ofs << "% reads found not overlapped. <br /> The nonoverlapped read pairs may have insert size &lt;" << mOptions->overlapRequire;
    ofs << " or &gt;" << isizeLimit;
    ofs << ", or contain too much sequencing errors to be detected as overlapped.";
    ofs <<"</div>\n";
    
    ofs << "\n<script type=\"text/javascript\">" << endl;
    string json_str = "var data=[";

    json_str += "{";
    json_str += "x:[" + Stats::list2string(x, total) + "],";
    json_str += "y:[" + Stats::list2string(percents, total) + "],";
    json_str += "name: 'Percent (%)  ',";
    json_str += "type:'bar',";
    json_str += "line:{color:'rgba(128,0,128,1.0)', width:1}\n";
    json_str += "}";

    json_str += "];\n";

    json_str += "var layout={title:'Insert size distribution (" + to_string(unknownPercents) + "% reads are with unknown length)', xaxis:{title:'Insert size'}, yaxis:{title:'Read percent (%)'}};\n";
    json_str += "Plotly.newPlot('plot_insert_size', data, layout);\n";

    ofs << json_str;
    ofs << "</script>" << endl;

    delete[] x;
    delete[] percents;
}

//void HtmlReporter::printGenomeCoverage(ofstream& ofs, Genomes* g) {
//    ofs << "<div class='section_div'>\n";
//    ofs << "<div class='section_title' onclick=showOrHide('genome_coverage')><a name='result'>Genome coverages for file: <I>" << mOptions->genomeFile << "</I><font color='#88CCFF' > (click to show/hide) </font></a></div>\n";
//    ofs << "<div id='sort_by_div' style='text-align:center;padding:10px;font-size:12px;'>Order by: <font color='#FF6600'>Coverage rate</font> | <a href='javascript:switch_sort();'>Bases on target</a></div>\n";
//    ofs << "<div id='genome_coverage'>\n";
//
//    g->reportHtml(ofs);
//
//    ofs << "</div>\n"; //result
//
//    ofs << "</div>\n"; // section_div
//}
//
//void HtmlReporter::reportKmerCollection(ofstream& ofs, KmerCollection* kc) {
//    ofs << "<div class='section_div'>\n";
//    ofs << "<div class='section_title' onclick=showOrHide('kcr')><a name='result'>Detection result for k-mer collection file: <I>" << mOptions->kmerCollectionFile << "</I><font color='#88CCFF' > (click to show/hide) </font></a></div>\n";
//    ofs << "<div id='kcr'>\n";
//
//    ofs << "<div id='kcr_result'>\n";
//    
//    kc->reportHTML(ofs);
//
//    ofs << "</div>\n"; //kcr_result
//    ofs << "</div>\n"; //kcr
//
//    ofs << "</div>\n"; // section_div
//}
//
//void HtmlReporter::reportKmerHits(ofstream& ofs, Kmer* kmer) {
//    ofs << "<div id='kmer_hits_figure'>\n";
//    ofs << "<div class='figure' id='plot_kmer_hits' style='height:300px;width:98%'></div>\n";
//    ofs << "</div>\n";
//    
//    ofs << "\n<script type=\"text/javascript\">" << endl;
//    string json_str = "var data=[";
//
//    json_str += "{";
//    json_str += "x:[" + kmer->getPlotX() + "],";
//    json_str += "y:[" + kmer->getPlotY() + "],";
//    json_str += "name: 'Hit',";
//    json_str += "type:'bar',";
//    json_str += "line:{color:'rgba(128,0,128,1.0)', width:1}\n";
//    json_str += "}";
//
//    json_str += "];\n";
//
//    json_str += "var layout={title:'Unique k-mer hits (" + to_string(kmer->getKmerCount()) + " k-mer keys)', xaxis:{tickangle:60, tickfont:{size: 8,color: '#bc6f98'}},yaxis:{title:'Hit'}};\n";
//    json_str += "Plotly.newPlot('plot_kmer_hits', data, layout);\n";
//
//    ofs << json_str;
//    ofs << "</script>" << endl;
//}

void HtmlReporter::reportDuplication(ofstream& ofs) {

    ofs << "<div id='duplication_figure'>\n";
    ofs << "<div class='figure' id='plot_duplication' style='height:400px;'></div>\n";
    ofs << "</div>\n";
    
    ofs << "\n<script type=\"text/javascript\">" << endl;
    string json_str = "var data=[";

    int total = mOptions->duplicate.histSize - 2;
    long *x = new long[total];
    double allCount = 0;
    for(int i=0; i<total; i++) {
        x[i] = i+1;
        allCount += mDupHist[i+1];
    }
    double* percents = new double[total];
    memset(percents, 0, sizeof(double)*total);
    if(allCount > 0) {
        for(int i=0; i<total; i++) {
            percents[i] = (double)mDupHist[i+1] * 100.0 / (double)allCount;
        }
    }
    int maxGC = total;
    double* gc = new double[total];
    for(int i=0; i<total; i++) {
        gc[i] = (double)mDupMeanGC[i+1] * 100.0;
        // GC ratio will be not accurate if no enough reads to average
        if(percents[i] <= 0.05 && maxGC == total)
            maxGC = i;
    }

    json_str += "{";
    json_str += "x:[" + Stats::list2string(x, total) + "],";
    json_str += "y:[" + Stats::list2string(percents, total) + "],";
    json_str += "name: 'Read percent (%)  ',";
    json_str += "type:'bar',";
    json_str += "line:{color:'rgba(128,0,128,1.0)', width:1}\n";
    json_str += "},";

    json_str += "{";
    json_str += "x:[" + Stats::list2string(x, maxGC) + "],";
    json_str += "y:[" + Stats::list2string(gc, maxGC) + "],";
    json_str += "name: 'Mean GC ratio (%)  ',";
    json_str += "mode:'lines',";
    json_str += "line:{color:'rgba(255,0,128,1.0)', width:2}\n";
    json_str += "}";

    json_str += "];\n";

    json_str += "var layout={title:'duplication rate (" + to_string(mDupRate*100.0) + "%)', xaxis:{title:'duplication level'}, yaxis:{title:'Read percent (%) & GC ratio'}};\n";
    json_str += "Plotly.newPlot('plot_duplication', data, layout);\n";

    ofs << json_str;
    ofs << "</script>" << endl;

    delete[] x;
    delete[] percents;
    delete[] gc;
}

//void HtmlReporter::printDetectionResult(ofstream& ofs, Kmer* kmer) {
//    ofs << "<div class='section_div'>\n";
//    ofs << "<div class='section_title' onclick=showOrHide('result')><a name='result'>Detection result for target unique k-mer file: <I>" << mOptions->kmerFile << "</I><font color='#88CCFF' > (click to show/hide) </font></a></div>\n";
//    ofs << "<div id='result'>\n";
//
//    ofs << "<div id='detection_result'>\n";
//    ofs << "<table class='summary_table' style='width:800px'>\n";
//    string result;
//    if(kmer->getMeanHit() >= mOptions->positiveThreshold)
//        result = "<font color='red'><B>POSITIVE<B></font>";
//    else
//        result = "NEGATIVE";
//    outputRow(ofs, "Detection result for target k-mer file:", result);
//    outputRow(ofs, "Mean depth of k-mer coverage:", to_string(kmer->getMeanHit()));
//    outputRow(ofs, "Threshold to be positive:", to_string(mOptions->positiveThreshold));
//    ofs << "</table>\n";
//
//    ofs << "</div>\n"; //detection_result
//
//    reportKmerHits(ofs, kmer);
//
//    ofs << "</div>\n"; //result
//
//    ofs << "</div>\n"; // section_div
//}

void HtmlReporter::report(std::vector<std::map<std::string, std::vector<std::pair<std::string, Genotype>>>> & sortedAllGenotypeMapVec,
        std::map<std::string, std::map<std::string, LocSnp>> & allSnpsMap,
        FilterResult* result, Stats* preStats1, Stats* postStats1, Stats* preStats2, Stats* postStats2) {
    ofstream ofs;
    ofs.open(mOptions->htmlFile, ifstream::out);

    printHeader(ofs);

    ofs << "<h1 style='text-align:left;'><a href='https://github.com/seq2sat' target='_blank' style='color:#009900;text-decoration:none;'>Seq2Sat Report</a </h1>"<<endl;
    ofs << "<div style='font-size:12px;font-weight:normal;text-align:left;color:#666666;padding:5px;'>" << "Sample: " << basename(mOptions->prefix) << "</div>" << endl;

    ofs << "<div class='section_div'>\n";
    ofs << "<div class='section_title' onclick=showOrHide('genotype')><a name='genotype'>All genotypes <font color='#88CCFF' > (click to show/hide) </font></a></div>\n";
    ofs << "<div id='genotype'  style='display:none'>\n";
    
    if(mOptions->mVarType == ssr){
        reportAllGenotype(ofs, sortedAllGenotypeMapVec);
    } else {
        reportAllSnps(ofs, allSnpsMap);
    }

    ofs << "</div>\n";
    ofs << "</div>\n";

    printSummary(ofs, result, preStats1, postStats1, preStats2, postStats2);

    ofs << "<div class='section_div'>\n";
    ofs << "<div class='section_title' onclick=showOrHide('before_filtering')><a name='summary'>Original data <font color='#88CCFF' > (click to show/hide) </font></a></div>\n";
    ofs << "<div id='before_filtering'  style='display:none'>\n";

    if(preStats1) {
        preStats1 -> reportHtml(ofs, "Original data", "read1");
    }

    if(preStats2) {
        preStats2 -> reportHtml(ofs, "Original data", "read2");
    }

    ofs << "</div>\n";
    ofs << "</div>\n";

    ofs << "<div class='section_div'>\n";
    ofs << "<div class='section_title' onclick=showOrHide('after_filtering')><a name='summary'>Clean data used for detection <font color='#88CCFF' > (click to show/hide) </font></a></div>\n";
    ofs << "<div id='after_filtering'  style='display:none'>\n";

    if(postStats1) {  
        string name = "read1";
        postStats1 -> reportHtml(ofs, "Clean data used for detection", name);
    }

    if(postStats2) {
        postStats2 -> reportHtml(ofs, "Clean data used for detection", "read2");
    }

    ofs << "</div>\n";
    ofs << "</div>\n";

    printFooter(ofs);

}

void HtmlReporter::reportAllSnps(ofstream& ofs, std::map<std::string, std::map<std::string, LocSnp>> & allSnpsMap){
    for(auto & it : allSnpsMap){
        std::vector<std::string> x_vec;
        x_vec.reserve(it.second.size());
        std::vector<int> y_vec;
        y_vec.reserve(it.second.size());
        std::vector<double> bar_width_vec;
        bar_width_vec.reserve(it.second.size());
        for(auto & it2 : it.second){
            x_vec.push_back(it2.second.getGenotype());
            y_vec.push_back(it2.second.numReads);
            bar_width_vec.push_back(0.5);
            it2.second.print();
        }
        reportEachSnpGenotype(ofs, it.first, it.second, x_vec, y_vec, bar_width_vec);
        x_vec.clear();
        x_vec.shrink_to_fit();
        y_vec.clear();
        y_vec.shrink_to_fit();
        bar_width_vec.clear();
        bar_width_vec.shrink_to_fit();
    }
}

void HtmlReporter::reportEachSnpGenotype(ofstream& ofs, std::string marker, std::map<std::string, LocSnp> & snpsMap,
            std::vector<std::string> & x_vec, std::vector<int> & y_vec, std::vector<double> & bar_width_vec){
    
    std::string subsection = "Marker: " + marker;
    std::string divName = replace(subsection, " ", "_");
    divName = replace(divName, ":", "_");
    std::string title = marker;

    ofs << "<div class='subsection_title' onclick=showOrHide('" + divName + "')><a name='" + subsection + "'>" + subsection + "<font color='#88CCFF' > (click to show/hide) </font></a></div>\n";
    //ofs << "<div class='subsection_title'><a title='click to hide/show' onclick=showOrHide('" + divName + "')>" + subsection + "</a></div>\n";
    ofs << "<div id='" + divName + "'>\n";
    ofs << "<div class='sub_section_tips'>Value of each genotype will be shown on mouse over.</div>\n";
    
    ofs << "<div class='figure' id='plot_" + divName + "'></div>\n";

    ofs << "<div class='sub_section_tips'>SNPs are highlighted with red background while SNPs in reference are highlighted with orange background</div>\n";
    ofs << "<pre overflow: scroll>\n";
    ofs << "<table class='summary_table' style='width:100%'>\n";
    ofs <<  "<tr style='background:#cccccc'> <td>Marker</td><td>Genotype</td><td>N. of Reads</td><td>Length</td><td align='center'>Sequence</td></tr>\n";
    auto it = mOptions->mLocSnps.refLocMap.find(marker);
    if(it != mOptions->mLocSnps.refLocMap.end()) {
        ofs << "<tr style='color:blue'>";
        ofs << "<td>Reference</td>" <<
                "<td>" + it->second.getGenotype() + "</td>" <<
                "<td>N.A.</td>" <<
                "<td>" + std::to_string(it->second.ref.mStr.length()) + "</td>" <<
                "<td align='center'>" + highligher(it->second, true) + "</td>";
        ofs << "</tr>";
        outputRow(ofs, marker, snpsMap);
    }
    
    ofs << "</table>\n";
    ofs << "</pre>\n";
    ofs << "</div>\n";
    
    ofs << "\n<script type=\"text/javascript\">" << endl;
    
    string json_str = "var data=[";
    json_str += "{";
    json_str += "x:[" + Stats::list2string(x_vec, x_vec.size()) + "],";
    json_str += "y:[" + Stats::list2string(y_vec, y_vec.size()) + "],";
    json_str += "text: [" + Stats::list2string(x_vec, x_vec.size()) + "],";
    json_str += "width: [" + Stats::list2string(bar_width_vec, bar_width_vec.size()) + "],";
    json_str += "type:'bar', textposition: 'auto'";
    json_str += "}";
    json_str += "];\n";
    json_str += "var layout={title:'" + title + "',";
    json_str += "xaxis:{tickmode: 'array', tickvals:[" + Stats::list2string(x_vec, x_vec.size()) + "],  title:'" + "SNP" + "', automargin: true},";
    json_str += "yaxis:{title:'Number of reads', automargin: true}, ";
    json_str += "barmode: 'stack'};\n";
    json_str += "Plotly.newPlot('plot_" + divName + "', data, layout);\n";
    ofs << json_str;
    ofs << "</script>" << endl;
    
}

void HtmlReporter::reportAllGenotype(ofstream& ofs, std::vector<std::map<std::string, std::vector<std::pair<std::string, Genotype>>>> & sortedAllGenotypeMapVec){
    std::map<std::string, std::vector<UnitedGenotype>> uGenoMap;
    for (const auto & it : sortedAllGenotypeMapVec.at(0)) {
            std::multimap<int, Genotype> tmpGenoMMap;
            for(const auto & it2 : it.second){
                tmpGenoMMap.insert(std::make_pair(it2.second.baseLocVar.effectiveLen, it2.second));
            }
            std::vector<UnitedGenotype> tmpUGenoVec;
            for(auto it2 = tmpGenoMMap.begin(); it2 != tmpGenoMMap.end(); it2 = tmpGenoMMap.upper_bound(it2->first)){
                auto genoSize = it2->first;
                auto geno = tmpGenoMMap.equal_range(genoSize);
                UnitedGenotype uGeno;
                uGeno.marker = it.first;
                uGeno.effectiveLen = genoSize;
              
                for(std::multimap<int, Genotype>::iterator it3 = geno.first; it3 != geno.second; it3++){
                    uGeno.locVec.push_back(it3->second.baseLocVar);
                    uGeno.mraSVec.push_back(it3->second.baseLocVar.mra.mStr.length());
                    uGeno.numReadsVec.push_back(it3->second.numReads);
                    uGeno.effectiveSeqVec.push_back(it3->second.baseLocVar.effectiveSeq.mStr);
                    uGeno.seqNameVec.push_back(it3->second.baseLocVar.ff.mStr + it3->second.baseLocVar.mraName + it3->second.baseLocVar.rf.mStr);
               
                }
                tmpUGenoVec.push_back(uGeno);
            }
        uGenoMap[it.first] = tmpUGenoVec;
    }

    std::map<std::string, std::vector < UnitedGenotype>> uGenoMraMap;
    for (const auto & it : sortedAllGenotypeMapVec.at(1)) {
        std::multimap<int, Genotype> tmpGenoMMap;
        for (const auto & it2 : it.second) {
            tmpGenoMMap.insert(std::make_pair(it2.second.baseLocVar.mra.mStr.length(), it2.second));
        }
        
        std::vector<UnitedGenotype> tmpUGenoVec;
        for (auto it2 = tmpGenoMMap.begin(); it2 != tmpGenoMMap.end(); it2 = tmpGenoMMap.upper_bound(it2->first)) {
            auto mraSize = it2->first;
            auto geno = tmpGenoMMap.equal_range(mraSize);

            UnitedGenotype uGeno;
            uGeno.marker = it.first;
            uGeno.repuit = Sequence(it2->second.baseLocVar.repuit.mStr);
            uGeno.mraName = it2->second.baseLocVar.mraName;
            uGeno.mra = Sequence(it2->second.baseLocVar.mra.mStr);
            for (std::multimap<int, Genotype>::iterator it3 = geno.first; it3 != geno.second; it3++) {
                uGeno.locVec.push_back(it3->second.baseLocVar);
                uGeno.mraSVec.push_back(it3->second.baseLocVar.mra.mStr.length());
                uGeno.numReadsVec.push_back(it3->second.numReads);
            }
            tmpUGenoVec.push_back(uGeno);
        }
        uGenoMraMap[it.first] = tmpUGenoVec;
    }
    
    for (const auto & it : uGenoMap) {
        auto itt = sortedAllGenotypeMapVec.at(0).find(it.first);
        std::vector<int> x_vec;
        x_vec.reserve(it.second.size());
        std::vector<double> bar_width_vec;
        bar_width_vec.reserve(it.second.size());
        int nStacks = 0;
        for (const auto & it2 : it.second) {
            x_vec.push_back(it2.effectiveLen);
            if (it2.effectiveSeqVec.size() > nStacks) {
                nStacks = it2.effectiveSeqVec.size();
            }
            bar_width_vec.push_back(0.5);
        }
        std::map< std::string, std::vector<int>> stackMap;
        std::map< std::string, std::vector<std::string>> stackYlabMap;
        for(int i = 0; i < nStacks; i++){
            std::vector<int> numReadsVec;
            numReadsVec.reserve(nStacks);
            std::vector<std::string> yLabVec;
            yLabVec.reserve(nStacks);
            for(const auto & it2 : it.second){
                int numReads = 0;
                std::string read = "";
                if(i < it2.numReadsVec.size()){
                    numReads = it2.numReadsVec.at(i);
                    //read = it2.seqName;
                    read = it2.seqNameVec.at(i);
                }
                numReadsVec.push_back(numReads);
                yLabVec.push_back(read);
            }
            stackMap["Genotype" + std::to_string(i)] = numReadsVec;
            stackYlabMap["Genotype" + std::to_string(i)] = yLabVec;
        }
        
        auto itm = sortedAllGenotypeMapVec.at(1).find(it.first);
        auto itk = uGenoMraMap.find(it.first);
        
        std::vector<int> xmra_vec;
        xmra_vec.reserve(itk->second.size());
//        std::vector<std::string> ymra_label_vec;
//        ymra_label_vec.reserve(itk->second.size());
        std::vector<double> barmra_width_vec;
        barmra_width_vec.reserve(itk->second.size());
        int nmStacks = 0;
        for (const auto & it2 : itk->second) {
            xmra_vec.push_back(it2.mra.mStr.length());
            if (it2.locVec.size() > nmStacks) {
                nmStacks = it2.locVec.size();
            }
            //ymra_label_vec.push_back(it2.mraName);
            barmra_width_vec.push_back(0.5);
        }
        std::map< std::string, std::vector<int>> stackMraMap;
        std::map< std::string, std::vector<std::string>> stackYLabMMap;
        for (int i = 0; i < nmStacks; i++) {
            std::vector<int> numReadsVec;
            numReadsVec.reserve(nStacks);
            std::vector<std::string> yLabVec;
            yLabVec.reserve(nStacks);
            for (const auto & it2 : itk->second) {
                int numReads = i < it2.numReadsVec.size() ? it2.numReadsVec.at(i) : 0;
                numReadsVec.push_back(numReads);
                std::string mra = i < it2.numReadsVec.size() ? it2.locVec.at(i).mraName : "";
                yLabVec.push_back(mra);
            }
            stackMraMap["MRA" + std::to_string(i)] = numReadsVec;
            stackYLabMMap["MRA" + std::to_string(i)] = yLabVec;
        }
        reportEachGenotype(ofs, it.first, x_vec, stackMap, stackYlabMap, bar_width_vec, itt->second,
                           xmra_vec, stackMraMap, stackYLabMMap, barmra_width_vec, itm->second);
        x_vec.clear();
        x_vec.shrink_to_fit();
        bar_width_vec.clear();
        bar_width_vec.shrink_to_fit();
        stackMap.clear();
        stackYlabMap.clear();

        xmra_vec.clear();
        xmra_vec.shrink_to_fit();
        barmra_width_vec.clear();
        barmra_width_vec.shrink_to_fit();
        stackMraMap.clear();
        stackYLabMMap.clear();
    }
    uGenoMap.clear();
    uGenoMraMap.clear();
    sortedAllGenotypeMapVec.clear();
}

void HtmlReporter::reportEachGenotype(ofstream& ofs, std::string marker, 
                            std::vector<int> & x_vec, std::map< std::string, std::vector<int>> & stackMap, 
                            std::map< std::string, std::vector<std::string>> & stackYlabMap, std::vector<double> & bar_width_vec, 
                            std::vector<std::pair<std::string, Genotype>> & outGenotype,
                            std::vector<int> & xmra_vec, std::map< std::string, std::vector<int>> & stackMraMap, 
                            std::map< std::string, std::vector<std::string>> & stackYLabMMap, 
                            std::vector<double> & barmra_width_vec, 
                            std::vector<std::pair<std::string, Genotype>> & outGenotypeMra){
    std::string subsection = "Marker: " + marker;
    std::string divName = replace(subsection, " ", "_");
    divName = replace(divName, ":", "_");
    std::string title = marker;

    ofs << "<div class='subsection_title' onclick=showOrHide('" + divName + "')><a name='" + subsection + "'>" + subsection + "<font color='#88CCFF' > (click to show/hide) </font></a></div>\n";
    //ofs << "<div class='subsection_title'><a title='click to hide/show' onclick=showOrHide('" + divName + "')>" + subsection + "</a></div>\n";
    ofs << "<div id='" + divName + "'>\n";
    ofs << "<div class='sub_section_tips'>Value of each genotype will be shown on mouse over.</div>\n";
    
    ofs << "<div class='left'>\n";
    ofs << "<div class='figure' id='plot_" + divName + "'></div>\n";
    ofs << "</div>\n";
    ofs << "<div class='right'>\n";
    ofs << "<div class='figure' id='plot_mra" + divName + "'></div>\n";
    ofs << "</div>";
    
    
    //ofs << "<div class='figure' id='plot_" + divName + "'></div>\n";

    ofs << "<div class='sub_section_tips'>SNPs are highlighted with red background.</div>\n";
    ofs << "<table class='summary_table'>\n";
    ofs <<  "<tr style='background:#cccccc'> <td>Marker</td><td>Repeat unit</td><td>MRA base</td><td>MRA size</td><td>Genotype</td><td>N. of Reads</td><td align = 'right'>Forward flanking region</td><td align='center'>MRA</td><td align='left'>Reverse flanking region</td></tr>\n";
    auto it = mOptions->mLocVars.refLocMap.find(marker);
    if(it != mOptions->mLocVars.refLocMap.end()) {
        ofs << "<tr style='color:blue'>";
        ofs << "<td>Reference</td>" <<
                "<td>" + it->second.repuitAll.mStr + "</td>" <<
                "<td>" + it->second.mraBase+ "</td>" <<
                "<td>" + std::to_string(it->second.mra.length()) + "</td>" <<
                "<td>" + std::to_string(it->second.effectiveLen) + "</td>" <<
                "<td>N.A.</td>" <<
                "<td align='right'>" + highligher(it->second.ff.mStr, it->second.refSnpsSetffMap[basename(mOptions->prefix)]) + "</td>" << //could use <xmp>
                "<td align='center'>" + it->second.mraName + "</td>" <<
                "<td align='left'>" + highligher(it->second.rf.mStr, it->second.refSnpsSetrfMap[basename(mOptions->prefix)]) + "</td>";
        ofs << "</tr>";
        outputRow(ofs, marker, outGenotypeMra);
    }
    
    ofs << "</table>\n";
    
    ofs << "</div>\n";
    
    ofs << "\n<script type=\"text/javascript\">" << endl;
    
    string json_str = "";
    for(auto & it : stackMap) {
        auto itt = stackYlabMap.find(it.first);
        json_str += "var " + it.first + " = {";
        json_str += "x:[" + Stats::list2string(x_vec, x_vec.size()) + "],";
        json_str += "y:[" + Stats::list2string(it.second, it.second.size()) + "],";
        json_str += "text: [" + Stats::list2string(itt->second, itt->second.size()) + "],";
        json_str += "width: [" + Stats::list2string(barmra_width_vec, barmra_width_vec.size()) + "],";
        json_str += "type:'bar', textposition: 'auto',";
        json_str += "};\n";
    }
    
    json_str += "var data = [" + Stats::list2string(stackMap, stackMap.size()) + "];\n";
    json_str += "var layout={title:'" + title + "',";
    json_str += "xaxis:{tickmode: 'array', tickvals:[" + Stats::list2string(x_vec, x_vec.size()) + "],  title:'" + "Allele size (bp)" + "', automargin: true},";
    json_str += "yaxis:{title:'Number of reads', automargin: true}, ";
    json_str += "barmode: 'stack'};\n";
    json_str += "Plotly.newPlot('plot_" + divName + "', data, layout);\n";
    
    for(auto & it : stackMraMap) {
        auto itt = stackYLabMMap.find(it.first);
        json_str += "var " + it.first + " = {";
        json_str += "x:[" + Stats::list2string(xmra_vec, xmra_vec.size()) + "],";
        json_str += "y:[" + Stats::list2string(it.second, it.second.size()) + "],";
        json_str += "text: [" + Stats::list2string(itt->second, itt->second.size()) + "],";
        json_str += "width: [" + Stats::list2string(barmra_width_vec, barmra_width_vec.size()) + "],";
        json_str += "type:'bar', textposition: 'auto',";
        json_str += "};\n";
    }
    
    json_str += "var data = [" + Stats::list2string(stackMraMap, stackMraMap.size()) + "];\n";
    json_str += "var layout={title:'" + title + "',";
    json_str += "xaxis:{tickmode: 'array', tickvals:[" + Stats::list2string(xmra_vec, xmra_vec.size()) + "],  title:'" + "MRA size (bp)" + "', automargin: true},";
    json_str += "yaxis:{title:'Number of reads', automargin: true}, ";
    json_str += "barmode: 'stack'};\n";
    json_str += "Plotly.newPlot('plot_mra" + divName + "', data, layout);\n";
    ofs << json_str;
    ofs << "</script>" << endl;
}

std::string HtmlReporter::highligher(std::string & str, std::map<int, std::string> & snpsMap){
    if(snpsMap.empty()){
        return str;
    }
    std::string hstr = str;
    for(auto it = snpsMap.rbegin(); it != snpsMap.rend(); it++){
        hstr.insert(it->first + 1, "</mark>");
        hstr.insert(it->first, "<mark>");
    }
    return hstr;
}

std::string HtmlReporter::highligher(std::string & str, std::set<int> & snpsSet){
    if(snpsSet.empty()){
        return str;
    }
    std::string hstr = str;
    for(auto it = snpsSet.rbegin(); it != snpsSet.rend(); it++){
        hstr.insert(*it + 1, "</mark>");
        hstr.insert(*it, "<mark>");
    }
    return hstr;
}

std::string HtmlReporter::highligher(LocSnp & locSnp, bool ref){
    if(ref) {
        
        std::map<int, bool> hiliMap;
        for(const auto & it : locSnp.snpPosSet){
            hiliMap[it] = false;
        }
        
        for(const auto & it : locSnp.snpsMap){
            hiliMap[it.first] = true;
        }
        
        std::string hstr = locSnp.ref.mStr;
        
        for(auto it = hiliMap.rbegin(); it != hiliMap.rend(); it++) {
            if (it->second) {
                hstr.insert(it->first + 1, "</mark2>");
                hstr.insert(it->first, "<mark2>");
            } else {
                hstr.insert(it->first + 1, "</mark>");
                hstr.insert(it->first, "<mark>");
            }
        }
        return hstr;
    } else {
        if (locSnp.snpsMap.empty()) {
            return locSnp.ref.mStr;
        } else {
            std::string hstr = locSnp.ref.mStr;
            for (auto it = locSnp.snpsMap.rbegin(); it != locSnp.snpsMap.rend(); it++) {
                hstr.insert(it->first + 1, "</mark>");
                hstr.insert(it->first, "<mark>");
            }
            return hstr;
        }
    }
}

void HtmlReporter::printHeader(ofstream& ofs){
    ofs << "<html><head><meta http-equiv=\"content-type\" content=\"text/html;charset=utf-8\" />";
    ofs << "<title>Seq2Sat report at " + getCurrentSystemTime() + " </title>";
    printJS(ofs);
    printCSS(ofs);
    ofs << "</head>";
    ofs << "<body><div id='container'>";
}

void HtmlReporter::printCSS(ofstream& ofs){
    ofs << "<style type=\"text/css\">" << endl;
    //ofs << "td {border:1px solid #dddddd;padding:5px;font-size:12px;}" << endl;
    ofs << "td {border:1px; solid #dddddd;padding:5px;font-size:12px; width:1px; white-space:nowrap; border:1px solid gray;}" << endl;
    //ofs << "table {border:1px solid #999999;padding:2x;border-collapse:collapse; width:800px;}" << endl;
    ofs << "table {border:1px solid #999999;padding:2px;border-collapse:collapse; table-layout:auto; border:1px solid gray;}" << endl;
    ofs << ".col1 {width:320px; font-weight:bold;}" << endl;
    ofs << ".adapter_col {width:500px; font-size:10px;}" << endl;
    ofs << "img {padding:30px;}" << endl;
    ofs << "#menu {font-family:Consolas, 'Liberation Mono', Menlo, Courier, monospace;}" << endl;
    ofs << "#menu a {color:#0366d6; font-size:18px;font-weight:600;line-height:28px;text-decoration:none;font-family:-apple-system, BlinkMacSystemFont, 'Segoe UI', Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol'}" << endl;
    ofs << "a:visited {color: #999999}" << endl;
    ofs << ".alignleft {text-align:left;}" << endl;
    ofs << ".alignright {text-align:right;}" << endl;
    ofs << ".figure {width:auto; height:auto;}" << endl;
    //ofs << ".figure {width:100%; height:100%;}" << endl;
    //ofs << ".figure {width:800px;height:600px;}" << endl;
    ofs << ".header {color:#ffffff;padding:1px;height:20px;background:#000000;}" << endl;
    ofs << ".sub_section_div {font-size:13px;padding-left:10px;text-align:left; margin-top:10px;}" << endl;
    ofs << ".sub_section_title {color:#ffffff;font-size:13px;padding-left:10px;text-align:left;background:#009900; margin-top:10px; width:20%;}" << endl;
    ofs << ".section_title {color:#ffffff;font-size:14px;padding:7px;text-align:left;background:#009900; margin-top:10px;}" << endl;
    ofs << ".subsection_title {font-size:12px;padding:12px;margin-top:10px;text-align:left;color:blue}" << endl;
    ofs << "#container {text-align:center;padding:3px 3px 3px 10px;font-family:Arail,'Liberation Mono', Menlo, Courier, monospace;}" << endl;
    ofs << ".menu_item {text-align:left;padding-top:5px;font-size:18px;}" << endl;
    ofs << ".highlight {text-align:left;padding-top:30px;padding-bottom:30px;font-size:20px;line-height:35px;}" << endl;
    ofs << "#helper {text-align:left;border:1px dotted #fafafa;color:#777777;font-size:12px;}" << endl;
    ofs << "#footer {text-align:left;padding:15px;color:#ffffff;font-size:10px;background:#009900;font-family:Arail,'Liberation Mono', Menlo, Courier, monospace;}" << endl;
    ofs << ".kmer_table {text-align:center;font-size:8px;padding:2px;}" << endl;
    ofs << ".kmer_table td{text-align:center;font-size:8px;padding:0px;color:#ffffff}" << endl;
    ofs << ".sub_section_tips {color:#999999;font-size:10px;padding-left:12px;padding-bottom:3px;text-align:left;}" << endl;
    ofs << ".left, .right{display: inline-block}" << endl;
    ofs << "mark{background-color: red; color: black;}" << endl;
    ofs << "mark2{background-color: orange; color: black;}" << endl;
    ofs << "pre{overflow: auto; width:0; min-width:100%;}" << endl;
    ofs << "</style>" << endl;
}

void HtmlReporter::printJS(ofstream& ofs){
    ofs << "<script src='https://www.seq2fun.ca/resources/javascript/plotly-1.2.0.min.js'></script>" << endl;
    ofs << "\n<script type=\"text/javascript\">" << endl;
    ofs << "    function showOrHide(divname) {" << endl;
    ofs << "        div = document.getElementById(divname);" << endl;
    ofs << "        if(div.style.display == 'none')" << endl;
    ofs << "            div.style.display = 'block';" << endl;
    ofs << "        else" << endl;
    ofs << "            div.style.display = 'none';" << endl;
    ofs << "    }" << endl;
    ofs << "</script>" << endl;
}

const string HtmlReporter::getCurrentSystemTime(){
  auto tt = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
  struct tm* ptm = localtime(&tt);
  char date[60] = {0};
  sprintf(date, "%d-%02d-%02d      %02d:%02d:%02d",
    (int)ptm->tm_year + 1900,(int)ptm->tm_mon + 1,(int)ptm->tm_mday,
    (int)ptm->tm_hour,(int)ptm->tm_min,(int)ptm->tm_sec);
  return std::string(date);
}

void HtmlReporter::printFooter(ofstream& ofs){
    ofs << "\n</div>" << endl;
    ofs << "<div id='footer'> ";
    ofs << "<p>"<<command<<"</p>";
    ofs << "Seq2Sat " << SEQ2SAT_VER << ", at " << getCurrentSystemTime() << " </div>";
    ofs << "</body></html>";
}
