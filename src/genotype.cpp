#include "genotype.hpp"
#include <regex>
#include <cmath>
#include <numeric>
#include <string>

using namespace std;
using namespace tree;
using namespace genotype;

enum variationType {
    SNP = 1,
    INS = 2,
    DEL = 4
};

double phred_complement(double q) {
    double p = pow(10, (-q / 10));
    return -10 * log10(1 - p);
}

void to_upper(string& str) {
    for (char& c : str) {
        c = toupper(c);
    }
}

static int getIndexFromNucleotide(char nuc) {
    switch(nuc) {
        case 'A':
            return 0;
        case 'C':
            return 1;
        case 'G':
            return 2;
        case 'T':
            return 3;
        case '*':
            return 4;
        default:
            return 5;
    }
    return 5;
}

static char getNucleotideFromIndex(int index) {
    switch(index) {
        case 0:
            return 'A';
        case 1:
            return 'C';
        case 2:
            return 'G';
        case 3:
            return 'T';
        case 4:
            return '*';
        default:
            return 'N';
    }
}

double likelihood(
    int genotype_idx,
    const vector< vector<double> >& read_errs,
    const map<string, vector<double> >& deletions,
    const map<string, vector<double> >& insertions,
    int variation_type
) {
    vector<double> genotype_probs;
    vector<double> variants_probs;

    for (int i = 0; i < read_errs.size(); i++) {
        const auto& row = read_errs[i];
        if ((variation_type & variationType::SNP) && (genotype_idx == i)) {
            for (const auto& prob : row) {
                genotype_probs.push_back(phred_complement(prob));
            }
        } else {
            variants_probs.insert(variants_probs.end(), row.begin(), row.end());
        }
    }

    int ins_i = 0;
    for (const auto& insertion : insertions) {
        if ((variation_type & variationType::INS) && (genotype_idx == ins_i)) {
            for (const auto& prob : insertion.second) {
                genotype_probs.push_back(phred_complement(prob));
            }
        } else {
            variants_probs.insert(variants_probs.end(), insertion.second.begin(), insertion.second.end());
        }
        ins_i++;
    }
    
    int del_i = 0;
    for (const auto& deletion : deletions) {
        if ((variation_type & variationType::DEL) && (genotype_idx == del_i)) {
            for (const auto& prob : deletion.second) {
                genotype_probs.push_back(phred_complement(prob));
            }
        } else {
            variants_probs.insert(variants_probs.end(), deletion.second.begin(), deletion.second.end());
        }
        del_i++;
    }
    
    double genotype_prob = accumulate(genotype_probs.begin(), genotype_probs.end(), 0.0);
    double variant_prob = accumulate(variants_probs.begin(), variants_probs.end(), 0.0);
    
    return genotype_prob + variant_prob;

}

vector<double> genotype_likelihoods(
    const vector< vector<double> >& read_errs, const map<string, vector<double> >& deletions,
    const map<string, vector<double> >& insertions, const int8_t& site_info
) {
    vector<double> likelihoods;
    likelihoods.resize(5);
    for (auto i = 0; i < 5; ++i) {
        likelihoods[i] = numeric_limits<double>::max();
    }

    auto variation_types = site_info & 7;
    auto ref_nuc = site_info >> 3;

    for (int i = 0; i < read_errs.size(); ++i) {
        if (read_errs[i].empty() && i != ref_nuc) { continue; }
        likelihoods[i] = likelihood(i, read_errs, deletions, insertions, variationType::SNP);
    }

    if (variation_types & variationType::INS) {
        for (auto i = 0; i < insertions.size(); i++) {
            likelihoods.push_back(likelihood(i, read_errs, deletions, insertions, variationType::INS));
        }
    }

    if (variation_types & variationType::DEL) {
        for (auto i = 0; i < deletions.size(); i++) {
            likelihoods.push_back(likelihood(i, read_errs, deletions, insertions, variationType::DEL));
        }
    }

    return likelihoods;
}

vector<double> genotype_posteriors(
    const vector<double>& likelihoods, const map<string, vector<double> >& deletions,
    const map<string, vector<double> >& insertions, const int8_t& site_info, const mutationMatrices& mutmat
) {
    vector<double> posteriors;
    posteriors.resize(4);
    for (auto i = 0; i < 4; i++) {
        posteriors[i] = numeric_limits<double>::max();
    }

    auto ref_nuc = site_info >> 3;
    for (auto i = 0; i < 4; ++i) {
        if (likelihoods[i] != numeric_limits<double>::max()) {
            if (ref_nuc == 5) {
                posteriors[i] = likelihoods[i];
            } else {
                posteriors[i] = likelihoods[i] + mutmat.submat[ref_nuc][i];
            }
        }
    }

    size_t insertion_idx = 0;
    for (const auto& insertion : insertions) {
        if (insertion.first.size() > mutmat.insmat.size() - 1) {
            posteriors.push_back(likelihoods[5 + insertion_idx] + mutmat.insmat.back());
        } else {
            posteriors.push_back(likelihoods[5 + insertion_idx] + mutmat.insmat[insertion.first.size()]);
        }
        
        insertion_idx++; 
    }

    size_t deletion_idx = 0;
    for (const auto& deletion : deletions) {
        if (deletion.first.size() > mutmat.delmat.size() - 1) {
            posteriors.push_back(likelihoods[5 + insertion_idx + deletion_idx] + mutmat.delmat.back());
        } else {
            posteriors.push_back(likelihoods[5 + insertion_idx + deletion_idx] + mutmat.delmat[deletion.first.size()]);
        }

        insertion_idx++;
    }

    double min_score = *min_element(posteriors.begin(), posteriors.end());
    for (int i = 0; i < posteriors.size(); ++i) {
        posteriors[i] -= min_score;
    }
    return posteriors;
}

genotype::VariationSite::VariationSite(
    size_t sid, char ref, size_t position, int variation_types, const string& nucs,
    const vector<string>& insertion_seqs, const vector<string>& deletion_seqs, const string& errors,
    const mutationMatrices& mutMat
) {
    this->site_id = sid;
    this->ref_position = position;
    size_t offset = 0;
    this->site_info = (getIndexFromNucleotide(ref) << 3) + variation_types;
    this->read_errs.resize(5);
    assert(errors.size() == nucs.size() + insertion_seqs.size() + deletion_seqs.size());


    for (auto i = 0; i < nucs.size(); ++i) {
        if (getIndexFromNucleotide(nucs[i]) == 5) {
            this->read_errs[4].push_back(double(errors[i]) - 33.0);
        } else {
            this->read_errs[getIndexFromNucleotide(nucs[i])].push_back(double(errors[i]) - 33.0);
        }
    }
    offset += nucs.size();


    if (variation_types & variationType::INS) {
        for (auto i = 0; i < insertion_seqs.size(); ++i) {
            this->insertions[insertion_seqs[i]].push_back(double(errors[i + offset] - 33.0));
        }
        offset += insertion_seqs.size();
    }

    if (variation_types & variationType::DEL) {
        for (auto i = 0; i < deletion_seqs.size(); i++) {
            this->deletions[deletion_seqs[i]].push_back(double(errors[i + offset] - 33.0));
        }
    }

    this->likelihoods = genotype_likelihoods(this->read_errs, this->deletions, this->insertions, this->site_info);
    this->posteriors = genotype_posteriors(this->likelihoods, this->deletions, this->insertions, this->site_info, mutMat);
    
    for (size_t i = 0; i < 4; i++) {
        read_depth.push_back(this->read_errs[i].size());
    }
    for (const auto& ins : this->insertions) {
        read_depth.push_back(ins.second.size());
    }
    for (const auto& del : this->deletions) {
        read_depth.push_back(del.second.size());
    }

    for (size_t i = 0; i < this->posteriors.size(); i++) {
        if (this->posteriors[i] == 0.0) {
            this->most_probable_idx = i;
            break;
        }
    }
}

int parse_readbases(
    string readbase_string, const string& readbase_errors,
    char ref_nuc, string& nucs, vector<string>& insertion_seqs,
    vector<string>& deletion_seqs, string& errs 
) {
    int variation_types = 0;

    regex extra_regex("\\^.{1}|\\$");
    readbase_string = regex_replace(readbase_string, extra_regex, "");

    string snp_errs, ins_errs, del_errs;
    size_t cur_start = 0, cur_idx = 0;

    string basePairs = string("ATCGatcg*");

    while (cur_start < readbase_string.size()){

        bool is_last = (cur_start == readbase_string.size() - 1);


        if(basePairs.find(readbase_string[cur_start]) != std::string::npos){
            //SNP
            nucs += toupper(readbase_string[cur_start]);
            variation_types |= variationType::SNP;
            snp_errs += readbase_errors[cur_idx];
            cur_start += 1;
            cur_idx++;
        }else if(readbase_string[cur_start] == '.' || readbase_string[cur_start] == ',') {
            if (is_last){
                //SNP
                nucs += ref_nuc;
                snp_errs += readbase_errors[cur_idx];
                cur_start += 1;
                cur_idx++;
            }else{
                if(readbase_string[cur_start+1] == '-'){
                    //DEL
                    variation_types |= variationType::DEL;
                    del_errs += readbase_errors[cur_idx];
                    int indel_size = 0;
                    
                    cur_start+=2;
                    while(std::isdigit((int)readbase_string[cur_start])) {
                        
                        indel_size *= 10;
                        indel_size += (int)readbase_string[cur_start] - 48;
                        cur_start++;
                    }

                    string seq = readbase_string.substr(cur_start, indel_size);
                    to_upper(seq);
                    deletion_seqs.push_back(seq);
                    cur_start += indel_size;
                    
                    cur_idx++;
                }else if(readbase_string[cur_start+1] == '+'){

                    //INS
                    variation_types |= variationType::INS;
                    ins_errs += readbase_errors[cur_idx];
                    int indel_size = 0;

                    cur_start+=2;
                    while(std::isdigit((int)readbase_string[cur_start])) {
                        indel_size *= 10;
                        indel_size += (int)readbase_string[cur_start] - 48;
                        cur_start++;
                    }

                    string seq = readbase_string.substr(cur_start, indel_size);
                    to_upper(seq);
                    insertion_seqs.push_back(seq);

                    cur_start += indel_size;
                    
                    cur_idx++;
                }else{
                    //SNP
                    nucs += ref_nuc;
                    snp_errs += readbase_errors[cur_idx];
                    cur_start += 1;
                    cur_idx++;
                }
            }
            
        }else if(readbase_string[cur_start] == '-'){ //This means the same bp has an insertion and a deletion, for now we ignore the deletion TODO
            //SUBSTITION

            //variation_types |= variationType::DEL;
            //del_errs += readbase_errors[cur_idx];
            int indel_size = 0;
            
            cur_start++;
            while(std::isdigit((int)readbase_string[cur_start])) {
            
                indel_size *= 10;
                indel_size += (int)readbase_string[cur_start] - 48;
                cur_start++;
            }

            //string seq = readbase_string.substr(cur_start, indel_size);
            //to_upper(seq);
            //deletion_seqs.push_back(seq);
            cur_start += indel_size;
                    
            //cur_idx++;
        }else{
            cur_start++;
        }
    }

    errs = snp_errs + ins_errs + del_errs;
    return variation_types;
}


pair< vector<VariationSite>, pair<size_t, size_t> > genotype::getVariantSites(std::istream& fin, const mutationMatrices& mutMat) {
    regex variant_pattern("[ACGTacgt\\*]+");
    vector<VariationSite> candidateVariants;
    pair<size_t, size_t> maskRange(numeric_limits<size_t>::max(), 0);
    size_t site_id = 0;
    string line;

    while(getline(fin, line)) {
        

        vector<string> fields;
        stringSplit(line, '\t', fields);
        string readbases_string = fields[4];
        string readbases_errors = fields[5];
        size_t coverage = stoul(fields[3]);
        size_t position = stoul(fields[1]) - 1;
        if (position < maskRange.first)  { maskRange.first  = position; }
        if (position > maskRange.second) { maskRange.second = position; }
        if ((coverage > 0) && regex_search(readbases_string, variant_pattern)) {
            char ref_nuc = fields[2][0];
            string errs;
            string nucs;
            vector<string> insertion_seqs;
            vector<string> deletion_seqs;

            int variation_types = parse_readbases(readbases_string, readbases_errors, ref_nuc, nucs, insertion_seqs, deletion_seqs, errs);
            candidateVariants.emplace_back(VariationSite(
                site_id, ref_nuc, position, variation_types, nucs,
                insertion_seqs, deletion_seqs, errs, mutMat
            ));
            site_id++;

        }
        
    }
    
    return make_pair(candidateVariants, maskRange);
}

static void printVCFLine(const VariationSite& site, std::ofstream& fout) {
    size_t position  = site.ref_position + 1;
    int ref_nuc_idx  = site.site_info >> 3;
    size_t readDepth = 0;
    vector<string> altAlleles;
    vector<size_t> ad;
    vector<int> pl;
    string refAllele;
    int gt;
    

    string quality;
    if (ref_nuc_idx == 5) {
        quality = ".";
    } else {
        quality = to_string(int(site.posteriors[ref_nuc_idx]));
    }

    for (const auto& depth : site.read_depth) {
        readDepth += depth;
    }
   
    refAllele += getNucleotideFromIndex(ref_nuc_idx);

    // find longest deletion
    size_t ldl = 0; // longest deletion length
    string lds;     // longest deletion string
    for (const auto& del : site.deletions) {
        if (del.first.size() > ldl) {
            ldl = del.first.size();
            lds = del.first;
        }
    }

    if (!lds.empty()) {
        refAllele += lds;
    }

    // depth and likelihood for reference
    if (ref_nuc_idx == 5) {
        ad.push_back(0);
        pl.push_back(-1);
    } else {
        ad.push_back(site.read_depth[ref_nuc_idx]);
        pl.push_back(site.posteriors[ref_nuc_idx]);
    }
    if (pl.back() == 0.0) {
        gt = 0;
    }

    // substitutions
    for (int i = 0; i < 4; i++) {
        if (i == ref_nuc_idx) {
            continue;
        }
        if (site.posteriors[i] != numeric_limits<double>::max()) {
            altAlleles.push_back(getNucleotideFromIndex(i) + lds);
            ad.push_back(site.read_depth[i]);
            pl.push_back(site.posteriors[i]);
            if (pl.back() == 0.0) {
                gt = altAlleles.size();
            }
        }
    }

    // insertions
    size_t indelIdx = 0;
    for (const auto& ins : site.insertions) {
        altAlleles.push_back(getNucleotideFromIndex(ref_nuc_idx) + ins.first + lds);
        ad.push_back(site.read_depth[4 + indelIdx]);
        pl.push_back(site.posteriors[4 + indelIdx]);
        if (pl.back() == 0.0) {
            gt = altAlleles.size();
        }
        indelIdx += 1;
    }

    // deletions
    for (const auto& del : site.deletions) {
        size_t delSize = del.first.size();
        if (delSize == ldl) {
            altAlleles.push_back(refAllele.substr(0, 1));
        } else {
            altAlleles.push_back(getNucleotideFromIndex(ref_nuc_idx) + lds.substr(delSize, ldl - delSize));
        }
        ad.push_back(site.read_depth[4 + indelIdx]);
        pl.push_back(site.posteriors[4 + indelIdx]);
        if (pl.back() == 0.0) {
            gt = altAlleles.size();
        }
        indelIdx += 1;
    }

    fout << "ref" << "\t"                              // #CHROM
         << position << "\t"                           // POS
         << "." << "\t"                                // ID
         << refAllele << "\t";                         // REF
    for (size_t i = 0; i < altAlleles.size() - 1; i++) {
        fout << altAlleles[i] << ",";
    }                                    
    fout << altAlleles[altAlleles.size() - 1] << "\t"; // ALT
    fout << quality << "\t"                            // QUAL
         << "." << "\t"                                // FILTER
         << "DP=" << readDepth << "\t"                 // INFO
         << "GT:AD:GP" << "\t"                         // FORMAT
         << gt << ":";                                 // SAMPLE
    for (size_t i = 0; i < ad.size() - 1; i++) {
        fout << ad[i] << ",";
    }
    fout << ad[ad.size() - 1] << ":";
    for (size_t i = 0; i < pl.size() - 1; i++) {
        if (pl[i] == -1) {
            fout << ".,";
            continue;
        }
        fout << pl[i] << ",";
    }
    fout << pl[pl.size() - 1] << endl;
}

void genotype::printSamplePlacementVCF(std::istream& fin, const mutationMatrices& mutMat, bool variantOnly, size_t maskSize, std::ofstream& fout) {
    pair< vector<VariationSite>, pair<size_t, size_t> > variantSites = getVariantSites(fin, mutMat);
    const vector<VariationSite> candidateVariants = variantSites.first;
    if (candidateVariants.empty()) {
        return;
    }
    
    fout << "##fileformat=VCFv4.3\n"
         << "##contig=<ID=ref>\n"
         << "##INFO=<ID=DP,Number=1,Type=Integer,Description=Read Depth>\n"
         << "##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>\n"
         << "##FORMAT=<ID=AD,Number=1,Type=String,Description=Read depth for each allele>\n"
         << "##FORMAT=<ID=GP,Number=1,Type=String,Description=Genotype posterior probabilities in phred scale>\n";

    fout << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n";
    
    vector<VariationSite> groupSites;

    for (const auto& curSite : candidateVariants) {
        bool skip = true;
        for (int i = 0; i < curSite.posteriors.size(); i++) {
            if ((curSite.site_info >> 3) == 5 || (i != (curSite.site_info >> 3) && curSite.posteriors[i] != numeric_limits<double>::max())) {
                skip = false;
                break;
            }
        }
        if (skip) {
            continue;
        }

        if (variantOnly && curSite.most_probable_idx == (curSite.site_info >> 3)) {
            continue;
        }

        if (curSite.ref_position >= variantSites.second.first + maskSize && curSite.ref_position <= variantSites.second.second - maskSize) {
            printVCFLine(curSite, fout);
        }
    }
}