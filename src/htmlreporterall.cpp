#include "htmlreporterall.h"

extern string command;

HtmlReporterAll::HtmlReporterAll(Options* opt) {
    mOptions = opt;
}

HtmlReporterAll::HtmlReporterAll(const HtmlReporterAll& orig) {
}

HtmlReporterAll::~HtmlReporterAll() {
}

void HtmlReporterAll::report() {
    
    ofstream ofs;
    ofs.open(dirname(mOptions->prefix) + "Seq2Sat_all_samples.html");
    printHeader(ofs);

    ofs << "<h1 style='text-align:left;'><a href='https://github.com/seq2sat' target='_blank' style='color:#009900;text-decoration:none;'>Seq2Sat report</a </h1>" << endl;
    //string intro = "Created by <a href='https://github.com/OpenGene/fastv' style='color:#1F77B4'>fastv</a> v" + string(FASTV_VER)+ ", " + " an ultra-fast tool for fast identification of SARS-CoV-2 and other microbes from sequencing data";
    //ofs << "<div style='font-size:10px;font-weight:normal;text-align:left;color:#666666;padding:5px;'>" << basename(mOptions->prefix) << "</div>" << endl;
    ofs << "<div class='section_div'>\n";
    ofs << "<div class='section_title' onclick=showOrHide('AllSamples')><a name='genotype'>All Samples <font color='#88CCFF' > (click to show/hide) </font></a></div>\n";
    ofs << "<div id='AllSamples'  style='display:none'>\n";
    for (Sample & smp : mOptions->samples) {
        reportEachSample(ofs, smp);
    }
    ofs << "</div>\n";
    ofs << "</div>\n";
    mOptions->samples.clear();
    mOptions->samples.shrink_to_fit();
    printFooter(ofs);
}

void HtmlReporterAll::reportEachSample(ofstream& ofs, Sample & sample){
    std::string smpNm = basename(sample.prefix);
    ofs << "<div class='sub_section_div'>\n";
    ofs << "<div class='sub_section_title' onclick=showOrHide('" + smpNm + "')><a name='" + smpNm + "'>" + smpNm + "<font color='#88CCFF' > (click to show/hide) </font></a></div>\n";
    ofs << "<div id='" + smpNm + "'  style='display:none'>\n";

    std::map<std::string, std::vector < UnitedGenotype>> uGenoMap;
    
    for (const auto & it : sample.sortedAllGenotypeMapVec.at(0)) {
        std::multimap<int, Genotype> tmpGenoMMap;
        for (const auto & it2 : it.second) {
            tmpGenoMMap.insert(std::make_pair(it2.second.baseLocVar.effectiveLen, it2.second));
        }
        std::vector<UnitedGenotype> tmpUGenoVec;
        for (auto it2 = tmpGenoMMap.begin(); it2 != tmpGenoMMap.end(); it2 = tmpGenoMMap.upper_bound(it2->first)) {
            auto genoSize = it2->first;
            auto geno = tmpGenoMMap.equal_range(genoSize);
            UnitedGenotype uGeno;
            uGeno.marker = it.first;
            uGeno.effectiveLen = genoSize;
            for (std::multimap<int, Genotype>::iterator it3 = geno.first; it3 != geno.second; it3++) {
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
    for (const auto & it : sample.sortedAllGenotypeMapVec.at(1)) {
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
        auto itt = sample.sortedAllGenotypeMapVec.at(0).find(it.first);
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
        std::map< std::string, std::vector < std::string>> stackYlabMap;
        for (int i = 0; i < nStacks; i++) {
            std::vector<int> numReadsVec;
            numReadsVec.reserve(nStacks);
            std::vector<std::string> yLabVec;
            yLabVec.reserve(nStacks);
            for (const auto & it2 : it.second) {
                int numReads = 0;
                std::string read = "";
                if (i < it2.numReadsVec.size()) {
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

        auto itm = sample.sortedAllGenotypeMapVec.at(1).find(it.first);
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
        std::map< std::string, std::vector < std::string>> stackYLabMMap;
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
        reportEachGenotype(ofs, smpNm, it.first, x_vec, stackMap, stackYlabMap, bar_width_vec, itt->second,
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
    
    ofs << "</div>\n";
    ofs << "</div>\n";
}

void HtmlReporterAll::reportEachGenotype(ofstream& ofs, std::string sample, std::string marker, 
                            std::vector<int> & x_vec, std::map< std::string, std::vector<int>> & stackMap, 
                            std::map< std::string, std::vector<std::string>> & stackYlabMap, std::vector<double> & bar_width_vec, 
                            std::vector<std::pair<std::string, Genotype>> & outGenotype,
                            std::vector<int> & xmra_vec, std::map< std::string, std::vector<int>> & stackMraMap, 
                            std::map< std::string, std::vector<std::string>> & stackYLabMMap, 
                            std::vector<double> & barmra_width_vec, 
                            std::vector<std::pair<std::string, Genotype>> & outGenotypeMra){
    std::string subsection = "Sample: " + sample +  " Marker: " + marker;
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

    ofs << "<div class='sub_section_tips'>The lengths of both Forward and reverse primer are not included in allele size. SNPs are highlighted with red background.</div>\n";
    ofs << "<pre overflow: scroll>\n";
    //ofs << "<table class='summary_table' style='width:100%'>\n";
    ofs << "<table class='summary_table'>\n";
    ofs <<  "<tr style='background:#cccccc'> <td>Marker</td><td>Repeat unit</td><td>MRA base</td><td>MRA size</td><td>Allele size</td><td>N. of Reads</td><td align = 'right'>FF region</td><td align='center'>MRA</td><td align='left'>RF region</td></tr>\n";
    auto it = mOptions->mLocVars.refLocMap.find(marker);
    if(it != mOptions->mLocVars.refLocMap.end()) {
        ofs << "<tr style='color:blue'>";
        ofs << "<td>Reference</td>" <<
                "<td>" + it->second.repuitAll.mStr + "</td>" <<
                "<td>" + it->second.mraBase+ "</td>" <<
                "<td>" + std::to_string(it->second.mra.length()) + "</td>" <<
                "<td>" + std::to_string(it->second.effectiveLen) + "</td>" <<
                "<td>N.A.</td>" <<
                "<td align='right'>" + highligherSet(it->second.ff.mStr, it->second.refSnpsSetffMap[sample]) + "</td>" << //could use <xmp>
                "<td align='center'>" + it->second.mraName + "</td>" <<
                "<td align='left'>" + highligherSet(it->second.rf.mStr, it->second.refSnpsSetrfMap[sample]) + "</td>";
        ofs << "</tr>";
        outputRow(ofs, marker, outGenotypeMra);
    }
    
    ofs << "</table>\n";
    ofs << "</pre>\n";
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

void HtmlReporterAll::outputRow(ofstream& ofs, string key, long v) {
    ofs << "<tr><td class='col1'>" + key + "</td><td class='col2'>" + to_string(v) + "</td></tr>\n";
}

void HtmlReporterAll::outputRow(ofstream& ofs, string key, string v) {
    ofs << "<tr><td class='col1'>" + key + "</td><td class='col2'>" + v + "</td></tr>\n";
}

void HtmlReporterAll::outputRow(ofstream& ofs, std::string & marker, std::vector<std::pair<std::string, Genotype>> & outGenotype) {
    for (auto & it : outGenotype) {
        ofs << "<tr>";
        ofs << "<td>" + marker + "</td>" +
                "<td>" + it.second.baseLocVar.repuitAll.mStr + "</td>" +
                "<td>" + it.second.baseLocVar.mraBase + "</td>" +
                "<td>" + std::to_string(it.second.baseLocVar.mra.mStr.length()) + "</td>" +
                "<td bgcolor=" + (it.second.baseLocVar.totalReads > 0 ? "'green'>" : "'transparent'>") +
                std::to_string(it.second.baseLocVar.effectiveLen) + "</td>" +
                "<td>" + std::to_string(it.second.numReads) + "</td>" +
                "<td align='right'>" + highligher(it.second.baseLocVar.ff.mStr, it.second.baseLocVar.snpsMapff) + "</td>" +
                "<td align='center'>" + it.second.baseLocVar.mraName + "</td>" +
                "<td align='left'>" + highligher(it.second.baseLocVar.rf.mStr, it.second.baseLocVar.snpsMaprf) + "</td>";
//        ofs << "<td width=4%>" + marker + "</td>" +
//                "<td width=6%>" + it.second.baseLocVar.repuitAll.mStr + "</td>" +
//                "<td width=8%>" + it.second.baseLocVar.mraBase + "</td>" +
//                "<td width=4%>" + std::to_string(it.second.baseLocVar.mra.mStr.length()) + "</td>" +
//                "<td width=4%>" + std::to_string(it.second.baseLocVar.effectiveLen) + "</td>" +
//                "<td width=6%>" + std::to_string(it.second.numReads) + "</td>" +
//                "<td width=30% align='right'>" + highligher(it.second.baseLocVar.ff.mStr, it.second.baseLocVar.snpsMapff) + "</td>" +
//                "<td width=8% align='center'>" + it.second.baseLocVar.mraName + "</td>" +
//                "<td width=30% align='left'>" + highligher(it.second.baseLocVar.rf.mStr, it.second.baseLocVar.snpsMaprf) + "</td>";
        ofs << "</tr>\n";
    }
}

void HtmlReporterAll::outputRow(ofstream& ofs, std::string & marker, std::map<std::string, LocSnp> & snpsMap) {
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

std::string HtmlReporterAll::highligher(std::string & str, std::map<int, std::string> & snpsMap){
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

std::string HtmlReporterAll::highligherSet(std::string & str, std::set<int> & snpsSet){
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

std::string HtmlReporterAll::highligher(LocSnp & locSnp, bool ref){
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

string HtmlReporterAll::formatNumber(long number) {
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

string HtmlReporterAll::getPercents(long numerator, long denominator) {
    if(denominator == 0)
        return "0.0";
    else
        return to_string((double)numerator * 100.0 / (double)denominator);
}

void HtmlReporterAll::printHeader(ofstream& ofs){
    ofs << "<html><head><meta http-equiv=\"content-type\" content=\"text/html;charset=utf-8\" />";
    ofs << "<title>Seq2Sat report at " + getCurrentSystemTime() + " </title>";
    printJS(ofs);
    printCSS(ofs);
    ofs << "</head>";
    ofs << "<body><div id='container'>";
}

void HtmlReporterAll::printCSS(ofstream& ofs){
    ofs << "<style type=\"text/css\">" << endl;
    //ofs << "td {border:1px solid #dddddd;padding:5px;font-size:12px;}" << endl;
    ofs << "td {border:1px; solid #dddddd;padding:5px;font-size:12px; width:1px; white-space:nowrap; border:1px solid gray;}" << endl;
    //ofs << "table {border:1px solid #999999;padding:2x;border-collapse:collapse; width:800px;}" << endl;
    ofs << "table {border:1px solid #999999;padding:2x;border-collapse:collapse; table-layout:auto; border:1px solid gray;}" << endl;
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

void HtmlReporterAll::printJS(ofstream& ofs){
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

const string HtmlReporterAll::getCurrentSystemTime(){
  auto tt = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
  struct tm* ptm = localtime(&tt);
  char date[60] = {0};
  sprintf(date, "%d-%02d-%02d      %02d:%02d:%02d",
    (int)ptm->tm_year + 1900,(int)ptm->tm_mon + 1,(int)ptm->tm_mday,
    (int)ptm->tm_hour,(int)ptm->tm_min,(int)ptm->tm_sec);
  return std::string(date);
}

void HtmlReporterAll::printFooter(ofstream& ofs){
    ofs << "\n</div>" << endl;
    ofs << "<div id='footer'> ";
    ofs << "<p>"<<command<<"</p>";
    ofs << "Seq2Sat " << SEQ2SAT_VER << ", at " << getCurrentSystemTime() << " </div>";
    ofs << "</body></html>";
}
