#include "genotype.hpp"
#include <regex>
#include <cmath>
#include <numeric>
#include <string>
#include <iomanip>

using namespace std;
using namespace   seed_annotated_tree;
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

double third_phred(double q) {
    double p = pow(10, (-q / 10)) / 3.0;
    return -10 * log10(p);
}

void to_upper(string& str) {
    for (char& c : str) {
        c = toupper(c);
    }
}

static int getIndexFromNucleotide(char nuc) {
    switch(nuc) {
        case 'A':
        case 'a':
            return 0;
        case 'C':
        case 'c':
            return 1;
        case 'G':
        case 'g':
            return 2;
        case 'T':
        case 't':
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

inline std::vector<std::vector<double>> phredMatrix2ProbMatrix(const std::vector<std::vector<double>>& phredMatrix) {
  std::vector<std::vector<double>> prob_matrix = phredMatrix;
  for(int i = 0; i < prob_matrix.size(); i++) {
    for(int j = 0; j < prob_matrix[i].size(); j++) {
      prob_matrix[i][j] = pow(10, -prob_matrix[i][j] / 10);
    }
  }
  return prob_matrix;
}

inline std::vector<std::vector<double>> probMatrix2PhredMatrix(const std::vector<std::vector<double>>& probMatrix) {
  std::vector<std::vector<double>> phred_matrix = probMatrix;
  for(int i = 0; i < phred_matrix.size(); i++) {
    for(int j = 0; j < phred_matrix[i].size(); j++) {
      phred_matrix[i][j] = -10 * log10(phred_matrix[i][j]);
    }
  }
  return phred_matrix;
}

inline double getAverageMutationRate(const std::vector<std::vector<double>>& matrix) {
  if (matrix.empty() || matrix.size() != matrix[0].size()) throw std::invalid_argument("Matrix must be square and non-empty.");

  double sum = 0.0;
  size_t count = 0;
  for (size_t i = 0; i < matrix.size(); ++i) {
    for (size_t j = 0; j < matrix.size(); ++j) {
      if (i != j) {
        sum += matrix[i][j];
        ++count;
      }
    }
  }

  if (count == 0) throw std::runtime_error("No off-diagonal elements to calculate average.");
  return sum / count;
}

vector<std::vector<double>> genotype::scaleMutationSpectrum(const mutationMatrices& mutMat, double mutationRate) {
  vector<std::vector<double>> scaled_submat_phred = mutMat.submat;
  vector<std::vector<double>> scaled_submat_prob = phredMatrix2ProbMatrix(scaled_submat_phred);
  double avg_mutation_rate = getAverageMutationRate(scaled_submat_prob);
  double scale_factor = mutationRate / avg_mutation_rate;

  std::vector<double> scaled_submat_prob_row_sums(scaled_submat_prob.size());
  for(int i = 0; i < scaled_submat_prob.size(); i++) {
    for(int j = 0; j < scaled_submat_prob[i].size(); j++) {
      if (i == j) continue;
      scaled_submat_prob[i][j] *= scale_factor;
      scaled_submat_prob_row_sums[i] += scaled_submat_prob[i][j];
    }
  }

  for (int i = 0; i < scaled_submat_prob.size(); i++) {
    scaled_submat_prob[i][i] = 1 - scaled_submat_prob_row_sums[i];
  }

  return probMatrix2PhredMatrix(scaled_submat_prob);
}

std::vector<char> parse_alts(const std::string& alts_str) {
  std::vector<char> alts;
  std::istringstream alts_stream(alts_str);
  std::string alt;
  while (std::getline(alts_stream, alt, ',')) {
    if (alt.size() > 1) {
      throw std::runtime_error("Error: alt allel parsing error..");
    }
    alts.push_back(alt[0]);
  }
  return alts;
}

std::tuple<int, std::vector<double>, std::vector<int>, std::string, std::string> parse_sample_formats(const std::string& sample_formats_str) {
  std::vector<std::string> sample_formats;
  std::istringstream sample_formats_stream(sample_formats_str);
  std::string sample_format;
  while (std::getline(sample_formats_stream, sample_format, ':')) {
    sample_formats.push_back(sample_format);
  }

  int gt = std::stoi(sample_formats[0]);
  std::vector<double> pls;
  std::istringstream pls_stream(sample_formats[1]);
  std::string pl;
  while (std::getline(pls_stream, pl, ',')) {
    pls.push_back(std::stod(pl));
  }

  std::vector<int> ads;
  std::istringstream ads_stream(sample_formats[2]);
  std::string ad;
  while (std::getline(ads_stream, ad, ',')) {
    ads.push_back(std::stoi(ad));
  }

  return std::make_tuple(gt, pls, ads, sample_formats[1], sample_formats[2]);
}

std::string genotype::applyMutationSpectrum(const std::string& line, const std::vector<std::vector<double>>& scaled_submat) {
  std::vector<std::string> fields;
  std::istringstream line_stream(line);
  std::string field;
  while (std::getline(line_stream, field, '\t')) {
    fields.push_back(field);
  }

  if (fields.size() < 10 || fields[0] == "#CHROM") {
    return line;
  }
  
  if (fields.size() != 10) {
    throw std::runtime_error("Couldn't parse VCF. Unrecognized number of fields.");
  }

  if (fields[4] == ".") {
    return "";
  } else if (fields[7].substr(0, 2) != "DP") {
    if (fields[9][0] == '0') return "";
    else                     return line;
  } else if (getIndexFromNucleotide(fields[3][0]) > 3) {
    if (fields[9][0] == '0') return "";
    else                     return line;
  }

  if (fields[3].size() > 1) {
    throw std::runtime_error("Error: reference allele parsing error.");
  }

  int ref_nuc_idx = getIndexFromNucleotide(fields[3][0]);
  std::vector<char> alts = parse_alts(fields[4]);


  auto [gt, pls, ads, pls_string, ads_string] = parse_sample_formats(fields[9]);

  std::vector<double> gls;
  if (alts.size() + 1 == pls.size()) {
    gls = pls;
  } else {
    for (int i = 0; i < pls.size(); i += 2) {
      gls.push_back(pls[i]);
      if (gls.size() == alts.size() + 1) {
        break;
      }
    }
  }

  gls[0] += scaled_submat[ref_nuc_idx][ref_nuc_idx];
  for (int i = 1; i < gls.size(); i++) {
    gls[i] = gls[i] + scaled_submat[ref_nuc_idx][getIndexFromNucleotide(alts[i - 1])];
  }

  double min_gl = *std::min_element(gls.begin(), gls.end());
  int min_gl_index;
  for (int i = 0; i < gls.size(); i++) {
    gls[i] -= min_gl;
    if (gls[i] == 0) min_gl_index = i;
  }

  if (min_gl_index == 0) return "";

  gt = min_gl_index;
  double qual = gls[0];

  std::stringstream ss;
  ss << fields[0] << "\t" << fields[1] << "\t" << fields[2] << "\t" << fields[3] << "\t" << fields[4] << "\t" << std::fixed << std::setprecision(4) << qual << "\t" << fields[6] << "\t" << fields[7] << "\t" << fields[8] << "\t" << gt << ":" << pls_string << ":" << ads_string;
  return ss.str();
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
          // for (const auto& prob : row) {
          //   variants_probs.push_back(third_phred(prob));
          // }
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
    const map<string, vector<double> >& insertions, const int8_t& site_info, const char& ref_nuc
) {
    vector<double> likelihoods;
    likelihoods.resize(5);
    for (auto i = 0; i < 5; ++i) {
        likelihoods[i] = numeric_limits<double>::max();
    }

    auto variation_types = site_info & 7;
    auto ref_nuc_idx = site_info >> 3;

    for (int i = 0; i < read_errs.size(); ++i) {
        if (read_errs[i].empty() && i != ref_nuc_idx) { continue; }
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
      int64_t insSize = insertion.first.size();

      if (mutmat.insmat.find(insSize) != mutmat.insmat.end()) {
        if (mutmat.insmat.at(insSize) > mutmat.maxInsLogProb) {
          posteriors.push_back(likelihoods[5 + insertion_idx] + mutmat.maxInsLogProb);
        } else {
          posteriors.push_back(likelihoods[5 + insertion_idx] + mutmat.insmat.at(insSize));
        }
      } else {
        posteriors.push_back(likelihoods[5 + insertion_idx] + mutmat.maxInsLogProb);
      }

      insertion_idx++; 
    }

    size_t deletion_idx = 0;
    for (const auto& deletion : deletions) {
      int64_t delSize = deletion.first.size();

      if (mutmat.delmat.find(delSize) != mutmat.delmat.end()) {
        if (mutmat.delmat.at(delSize) > mutmat.maxDelLogProb) {
          posteriors.push_back(likelihoods[5 + insertion_idx + deletion_idx] + mutmat.maxDelLogProb);
        } else {
            posteriors.push_back(likelihoods[5 + insertion_idx + deletion_idx] + mutmat.delmat.at(delSize));
        }
      } else {
        posteriors.push_back(likelihoods[5 + insertion_idx + deletion_idx] + mutmat.maxDelLogProb);
      }

      deletion_idx++;
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
    this->ref_nuc = ref;
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

    this->likelihoods = genotype_likelihoods(this->read_errs, this->deletions, this->insertions, this->site_info, this->ref_nuc);
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

static double get_qual_for_ambiguous_ref(const VariationSite& site) {
  double qual = 0.0;
  for (int i = 0; i < 4; i++) {
      const auto& row = site.read_errs[i];
      qual += accumulate(row.begin(), row.end(), 0.0);
  }

  for (const auto& insertion : site.insertions) {
    qual += accumulate(insertion.second.begin(), insertion.second.end(), 0.0);
  }

  for (const auto& deletion : site.deletions) {
    qual += accumulate(deletion.second.begin(), deletion.second.end(), 0.0);
  }

  vector<double> filtered_likelihoods;
  for (size_t i = 0; i < site.likelihoods.size(); ++i) {
      if (i != 3 || site.likelihoods[i] != std::numeric_limits<double>::max()) {
          filtered_likelihoods.push_back(site.likelihoods[i]);
      }
  }


  return qual - *min_element(filtered_likelihoods.begin(), filtered_likelihoods.end());
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
        quality = to_string(int(get_qual_for_ambiguous_ref(site)));
    } else {
        quality = to_string(int(site.posteriors[ref_nuc_idx]));
    }
   
    refAllele += site.ref_nuc;

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
        pl.push_back(stoi(quality));
    } else {
        ad.push_back(site.read_depth[ref_nuc_idx]);
        pl.push_back(site.posteriors[ref_nuc_idx]);
        readDepth += site.read_depth[ref_nuc_idx];
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
            readDepth += site.read_depth[i];
            if (pl.back() == 0.0) {
                gt = altAlleles.size();
            }
        }
    }

    // insertions
    size_t indelIdx = 0;
    for (const auto& ins : site.insertions) {
        altAlleles.push_back(site.ref_nuc + ins.first + lds);
        ad.push_back(site.read_depth[4 + indelIdx]);
        pl.push_back(site.posteriors[4 + indelIdx]);
        readDepth += site.read_depth[4 + indelIdx];
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
            altAlleles.push_back(site.ref_nuc + lds.substr(delSize, ldl - delSize));
        }
        ad.push_back(site.read_depth[4 + indelIdx]);
        pl.push_back(site.posteriors[4 + indelIdx]);
        readDepth += site.read_depth[4 + indelIdx];
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

void genotype::printSamplePlacementVCF(std::istream& fin, const mutationMatrices& mutMat, bool keep_alts, size_t maskSize, std::ofstream& fout) {
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

        if (!keep_alts && curSite.most_probable_idx == (curSite.site_info >> 3)) {
            continue;
        }

        if (curSite.ref_position >= variantSites.second.first + maskSize && curSite.ref_position <= variantSites.second.second - maskSize) {
            printVCFLine(curSite, fout);
        }
    }
}