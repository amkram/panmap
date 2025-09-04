#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/lzma.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <algorithm>
#include <iostream>
#include <vector>
#include <random>
#include <functional>
#include <set>
#include "panmanUtils.hpp"
#include "genotyping.hpp"


namespace po = boost::program_options;
namespace fs = boost::filesystem;

void makeFasta(const std::string& name, const std::string& seq, const std::string& path);
void makeDir(const std::string& path);
std::vector<int> genMutNum(const std::vector<double>& mutNum_double, std::mt19937& gen);
std::string readFastaSequence(const std::string& fastaPath);
std::string readFastaHeader(const std::string& fastaPath);
void sim(panmanUtils::Tree* T, const std::string& refNode, const std::string& out_dir, const std::string& prefix,
    const std::string& mut_spec_type, const std::vector<double>& num, const std::pair<int, int>& indel_len, const std::string& model,
    int n_reads, int rep, const genotyping::mutationMatrices& mutMat, unsigned seed, int cpus, bool no_reads, const std::string& in_fasta, bool include_internal);
void scaleMutationMatrices(genotyping::mutationMatrices& mutMat, double mutation_rate);

int main(int argc, char *argv[]) {
    std::cout << "What is my purpose?\nYou pass butter" << std::endl;
    try {
        po::options_description desc("Options");
        desc.add_options()
            ("help,h", "Produce help message")
            ("in_panman",    po::value<std::string>()->required(), "Path to tree.pman file (required).")
            ("in_mm",  po::value<std::string>()->default_value(""), "Use input mutation matrix file to model mutations")            
            ("in_fasta",     po::value<std::string>()->default_value(""), "Path to input FASTA genome (optional).")
            ("out_dir",   po::value<std::string>()->required(), "Output directory (required)")
            ("mutnum",    po::value<std::vector<double>>()->multitoken(), "Number of mutations for snp, insertion, and deletion [10 0 0].")
            ("mm_type", po::value<std::string>()->default_value("snp"), "Type of mutation matrix to use. Options: snp, indel, both. Default: snp.")
            ("n_replicates",       po::value<int>()->default_value(1), "Number of replicates to simulate [1].")
            ("n_reads",   po::value<int>()->default_value(100), "Number of reads to simulate [100].")
            ("model",     po::value<std::string>()->default_value("NovaSeq"), "InSilicoSeq error model [HiSeq]. Options: HiSeq, NextSeq, NovaSeq, MiSeq. For detail, visit InSilicoSeq github (https://github.com/HadrienG/InSilicoSeq).")
            ("perfect",   "Generate perfect reads with no sequencing errors (overrides --model)")
            ("cpus",      po::value<int>()->default_value(8), "Number of CPUs to use [8].")
            ("seed",      po::value<std::string>()->default_value("RANDOM"), "Random seed for simulation [default: random]")
            ("ref",       po::value<std::string>()->default_value("RANDOM"), "Reference node ID to use [default: random]")
            ("prefix",    po::value<std::string>()->default_value("sim"), "Output prefix [default: sim]")
            ("mutation_rate", po::value<double>()->default_value(-1), "Mutation rate to scale matrices [-1 = no scaling]")
            ("indel_len", po::value<std::vector<int>>()->multitoken(), "Min and max indel length [1 9]")
            ("include_internal", "Include internal nodes in random selection (default: leaves only)")
            ("no-reads",  "Skip read simulation")
        ;

        po::variables_map vm;
        po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
        if (vm.count("help")) {
            std::cout << desc << "\n";
            return 0;
        }
        po::notify(vm);
        
        // input variables
        std::string panmatPath    = vm["in_panman"].as<std::string>();
        std::string mut_spec      = vm["in_mm"].as<std::string>();
        std::string mut_spec_type = vm["mm_type"].as<std::string>();
        std::string out_dir       = vm["out_dir"].as<std::string>();
        std::string refNode       = vm["ref"].as<std::string>();
        std::string prefix        = vm["prefix"].as<std::string>();
        std::string model         = vm["model"].as<std::string>();
        std::string seedstr       = vm["seed"].as<std::string>();
        std::string in_fasta      = vm["in_fasta"].as<std::string>();
        double mutation_rate      = vm["mutation_rate"].as<double>();
        int n_reads               = vm["n_reads"].as<int>();
        int rep                   = vm["n_replicates"].as<int>();
        int cpus                  = vm["cpus"].as<int>();
        bool no_reads             = vm.count("no-reads") > 0;
        bool perfect              = vm.count("perfect") > 0;
        bool include_internal     = vm.count("include_internal") > 0;

        // Override model if perfect flag is set
        if (perfect) {
            model = "perfect";
        }

        // Check mut_spec
    genotyping::mutationMatrices mutMat = genotyping::mutationMatrices();
        if (!mut_spec.empty()) {
            if (fs::exists(mut_spec)) {
                std::ifstream mminf(mut_spec);
                genotyping::fillMutationMatricesFromFile(mutMat, mminf);
                mminf.close();
            } else {
                throw std::invalid_argument("--in_mm input file doesn't exist");
            }
        }

        if (mutation_rate != -1) {
          scaleMutationMatrices(mutMat, mutation_rate);
        }

        // Check out_dir input
        if (!fs::is_directory(fs::path(out_dir))) {
            throw std::invalid_argument("--out_dir input is not a valid directory");
        }

        std::ofstream logFile(out_dir + "/log.txt");
        logFile << "Running simulation with following parameters:\n"
                << "PanMAT: " << panmatPath << "\n"
                << "Reference node: " << refNode << "\n"
                << "Output directory: " << out_dir << "\n"
                << "Output prefix: " << prefix << "\n"
                << "Mutation specification: " << mut_spec << "\n"
                << "Number of reads: " << n_reads << "\n"
                << "Number of replicates: " << rep << "\n"
                << "Error model: " << model << "\n"
                << "Number of CPUs: " << cpus << "\n"
                << "Random seed: " << seedstr << "\n"
                << "No reads: " << no_reads << "\n";

        // Check multitoken inputs
        std::pair<int, int> indel_len;
        std::vector<double> mutnum_double;
        if (!vm.count("indel_len")) {
            indel_len = {1, 9};
        } else {
            const auto& indel_len_vec = vm["indel_len"].as<std::vector<int>>();
            if (indel_len_vec.size() != 2) {
                throw std::invalid_argument("--indel_len must have 2 inputs");
            }
            indel_len.first  = indel_len_vec[0];
            indel_len.second = indel_len_vec[1];
            if (indel_len.first < 1) {
                throw std::invalid_argument("--indel_len must be greater than or equal to 1");
            }
        }
        if (!vm.count("mutnum")) {
            mutnum_double = {10.0, 0.0, 0.0};
        } else {
            mutnum_double = vm["mutnum"].as<std::vector<double>>();
            if (mutnum_double.size() != 3) {
                throw std::invalid_argument("--mutnum must have 3 inputs");
            }
        }
        if (mut_spec_type != "snp" && mut_spec_type != "indel" && mut_spec_type != "both") {
            throw std::invalid_argument("--mm_type only supports snp, indel, or both");
        }

        // check model input
        const std::vector<std::string> acceptedModels = {"NextSeq", "NovaSeq", "HiSeq", "MiSeq", "perfect"};
        if (std::find(acceptedModels.begin(), acceptedModels.end(), model) == acceptedModels.end()) {
            throw std::invalid_argument("Unknown --model input, please double check");
        }


        // read treeeee
        // If the data structure loaded into memory is a PanMAT, it is pointed to by T
        panmanUtils::Tree *T = nullptr;
        // If the data structure loaded into memory is a PanMAN, it is pointed to by TG
        panmanUtils::TreeGroup *TG = nullptr;

        try {
          // Try to load PanMAN file directly into memory
          std::ifstream inputFile(panmatPath);
          boost::iostreams::filtering_streambuf< boost::iostreams::input> inPMATBuffer;
          inPMATBuffer.push(boost::iostreams::lzma_decompressor());
          inPMATBuffer.push(inputFile);
          std::istream inputStream(&inPMATBuffer);


          TG = new panmanUtils::TreeGroup(inputStream);

          
          inputFile.close();
        } catch (const std::exception &e) {
          std::cerr << "Attempting to load as PanMAT...\n";
          try {
            std::ifstream inputFile(panmatPath);
            boost::iostreams::filtering_streambuf< boost::iostreams::input> inPMATBuffer;
            inPMATBuffer.push(boost::iostreams::lzma_decompressor());
            inPMATBuffer.push(inputFile);
            std::istream inputStream(&inPMATBuffer);


            T = new panmanUtils::Tree(inputStream);


          } catch (const std::exception &e) {
            std::cerr << "Error: failed to load PanMAT file" << std::endl;
            return 1;
          }
        }
        if (TG != nullptr) {
          T = &(TG->trees[0]);
        }



        // check node
        if (refNode != "RANDOM" && T->allNodes.find(refNode) == T->allNodes.end()) {
            throw std::invalid_argument("Couldn't find --ref node on tree");
        }
        // GO time
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        if (seedstr != "RANDOM") {
            try {
                std::hash<std::string> hasher;
                seed = (unsigned)hasher(seedstr);
            } catch (const std::invalid_argument& e) {
                throw std::invalid_argument("Invalid seed value: cannot convert to unsigned");
            } catch (const std::out_of_range& e) {
                throw std::out_of_range("Invalid seed value: out of range for unsigned");
            }
        }
        logFile << "Using seed: " << seed << "\n";
        logFile.close();
        sim(T, refNode, out_dir, prefix, mut_spec_type, mutnum_double, indel_len, model, n_reads, rep, mutMat, seed, cpus, no_reads, in_fasta, include_internal);

    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}

void makeFasta(const std::string& name, const std::string& seq, const std::string& path) {
    if (!fs::exists(path)) {
        std::ofstream faos(path);
        faos << '>' << name << '\n';
        size_t linesize = 80;
        for (size_t i = 0; i < seq.size(); i += linesize) {
            faos << seq.substr(i, std::min(linesize, seq.size() - i)) << '\n';
        }
        faos.close();
    }
}

void makeDir(const std::string& path) {
    if (!fs::exists(path)) {
        fs::create_directory(path);
    }
}

std::string readFastaSequence(const std::string& fastaPath) {
    std::ifstream file(fastaPath);
    std::string line, sequence;
    
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open FASTA file: " + fastaPath);
    }
    
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') continue; // Skip header lines
        sequence += line;
    }
    
    file.close();
    return sequence;
}

std::string readFastaHeader(const std::string& fastaPath) {
    std::ifstream file(fastaPath);
    std::string line;
    
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open FASTA file: " + fastaPath);
    }
    
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            file.close();
            // Remove '>' and any trailing whitespace/comments
            std::string header = line.substr(1);
            size_t space_pos = header.find(' ');
            if (space_pos != std::string::npos) {
                header = header.substr(0, space_pos);
            }
            return header;
        }
    }
    
    file.close();
    throw std::runtime_error("No FASTA header found in file: " + fastaPath);
}

char getRandomChar(const std::vector<char>& charList, std::mt19937& gen) {
    std::uniform_int_distribution<> distr(0, charList.size() - 1);
    int index = distr(gen);

    return charList[index];
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

double getMinDouble(const std::vector<double>& doubles) {
    if (doubles.empty()) {
        return -1.0;
    }

    double minDouble = doubles[0];
    for (size_t i = 1; i < doubles.size(); i++) {
        minDouble = std::min(minDouble, doubles[i]);
    }

    return minDouble;

}


char subNuc(char ref, const genotyping::mutationMatrices& mutMat, std::mt19937& gen) {
    std::vector<char> bases = {'A', 'C', 'G', 'T'};
    int refIdx = getIndexFromNucleotide(ref);
    if (refIdx > 3) {
        return getRandomChar(bases, gen);
    }

    if (!mutMat.filled) {
        bases.erase(bases.begin() + refIdx);
        return getRandomChar(bases, gen);
    }

    std::vector<double> probs = mutMat.submat[refIdx];
    bases.erase(bases.begin() + refIdx);
    probs.erase(probs.begin() + refIdx);
    
    double sum = 0;
    for (size_t i = 0; i < probs.size(); i++) {
      probs[i] = pow(10, -probs[i] / 10);
      sum += probs[i];
    }
    for (size_t i = 0; i < probs.size(); i++) {
      probs[i] /= sum;
    }

    std::discrete_distribution<> distr(probs.begin(), probs.end());
    int index = distr(gen);
    return bases[index];
}

std::vector<double> convertMap(const std::unordered_map<long, double> &in) {
    int maxKey = std::max_element(in.begin(), in.end(), 
                                  [](const std::pair<long, double>& a, const std::pair<long, double>& b) {
                                      return a.first < b.first;
                                  })->first;
    std::vector<double> outVector(maxKey + 1);
    for (const auto& pair : in) {
        outVector[pair.first] = pair.second;
    }
    return outVector;
}
size_t genLen(const std::pair<int, int>& indel_len, const genotyping::mutationMatrices& mutMat, int type, std::mt19937& gen) {
    // Always use uniform distribution for indel lengths for now to test
    std::uniform_int_distribution<> distribLen(indel_len.first, indel_len.second);
    return distribLen(gen);
}
    // Original mutation matrix code (commented out for testing)
    /*
    if (!mutMat.filled) {
        std::uniform_int_distribution<> distribLen(indel_len.first, indel_len.second);
        return distribLen(gen);
    }

    std::vector<double> probs;
    switch (type) {
        case 2:
            probs = convertMap(mutMat.insmat);
            break;
        case 4:
            probs = convertMap(mutMat.delmat);
            break;
        default:
            break;
    }
    probs.erase(probs.begin());
    std::vector<int> wgts;
    std::vector<size_t> lens;
    double minProb = getMinDouble(probs);
    for (size_t i = indel_len.first; i < indel_len.second+1; i++) {
        lens.push_back(i);
        wgts.push_back(pow(10, (minProb - probs[i-1]) / 10) * 1000);
    }
    std::discrete_distribution<> distr(wgts.begin(), wgts.end());
    int index = distr(gen);
    return lens[index];
    */
    

void genMut(const std::string& curNode, const std::string& seq, const fs::path& fastaOut, const fs::path& vcfOut,
    const genotyping::mutationMatrices& mutMat, const std::vector<double> mutnum_double, const std::pair<int, int> indel_len,
    size_t beg, size_t end, std::mt19937& gen, unsigned seed)
{
    if (fs::exists(fastaOut) && fs::exists(vcfOut)) {
        return;
    }

    std::vector<char> bases = {'A', 'C', 'G', 'T'};
    std::ofstream faos(fastaOut.string(), std::ofstream::trunc);
    std::ofstream vros(vcfOut.string(), std::ofstream::trunc);

    std::string nseq = seq;
    std::vector< std::tuple<int, int, int> > muts;
    std::vector<std::string> vref;
    std::vector<int> vtp;
    std::vector<int> num = genMutNum(mutnum_double, gen);
    for (int i = 0; i < num.size(); i++) {
        for (int j = 0; j < num[i]; j++) {
            switch (i) {
                case 0:
                    vtp.push_back(1);
                    break;
                case 1:
                    vtp.push_back(2);
                    break;
                case 2:
                    vtp.push_back(4);
                    break;
            }
        }
    }
    std::uniform_int_distribution<> distribPos(beg, seq.size() - end);

    int c = 0;
    int sumnum = num[0] + num[1] + num[2];
    while (muts.size() < sumnum) {
        int curPos  = distribPos(gen);
        int varType = vtp[c % vtp.size()];
        int length;
        if (varType == 1) {
            length = 1;
        } else if (varType == 2) {
            length = genLen(indel_len, mutMat, 2, gen);
        } else {
            length = genLen(indel_len, mutMat, 4, gen);
        }
        
        bool posConflict = false;
        for (const auto& mut : muts) {
            if (get<1>(mut) == 4 && abs(curPos - get<0>(mut)) < (2 * get<2>(mut))) {
                posConflict = true;
                break;
            } else if (abs(curPos - get<0>(mut)) == 0) {
                posConflict = true;
                break;
            }

        }
        if (posConflict) {
            continue;
        }
        muts.emplace_back(std::make_tuple(curPos, varType, length));
        c++;
    }
    
    sort(muts.begin(), muts.end(), 
        [](const std::tuple<int, int, int>& a, const std::tuple<int, int, int>& b) {
            return get<0>(a) < std::get<0>(b);
        }
    );

    vros << "##fileformat=VCFv4.3\n"
        //  << "##contig=<ID=" + curNode + ">\n"
         << "##contig=<ID=ref>\n"
         << "##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>\n";
    vros << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + fastaOut.stem().string() + "\n";
    int offset = 0;
    for (int i = 0; i < sumnum; i++) {
        int pos = get<0>(muts[i]);
        int var = get<1>(muts[i]);
        
        switch (var) {
            case 1:
                {
                    vros << "ref\t" << std::to_string(pos + 1) + "\t.\t";
                    // vros << curNode << "\t" << std::to_string(pos + 1) + "\t.\t";
                    char ref  = seq[pos];
                    char mut  = subNuc(ref, mutMat, gen);
                    nseq[pos + offset] = mut;
                    vros << ref << "\t" << mut << "\t.\t.\t.\tGT\t1\n";
                    break;
                }
            case 2:
                {
                    vros << "ref\t" << std::to_string(pos + 1) + "\t.\t";
                    // vros << curNode << "\t" << std::to_string(pos + 1) + "\t.\t";
                    char ref = seq[pos];
                    std::string inss;
                    for (size_t j = 0; j < get<2>(muts[i]); j++) {
                        inss += getRandomChar(bases, gen);
                    }
                    nseq = nseq.substr(0, pos + offset + 1) + inss + nseq.substr(pos + offset + 1, nseq.size() - (pos + offset + 1));
                    vros << ref << "\t" << ref << inss << "\t.\t.\t.\tGT\t1\n";
                    offset += get<2>(muts[i]);
                    break;
                }
            case 4:
                {
                    vros << "ref\t" << std::to_string(pos + 1) + "\t.\t";
                    // vros << curNode << "\t" << std::to_string(pos + 1) + "\t.\t";
                    std::string ref = seq.substr(pos, get<2>(muts[i]) + 1);
                    char   del = seq[pos];
                    nseq = nseq.substr(0, pos + offset + 1) + nseq.substr(pos + offset + get<2>(muts[i]) + 1, nseq.size() - (pos + offset + get<2>(muts[i]) + 1));
                    vros << ref << "\t" << del << "\t.\t.\t.\tGT\t1\n";
                    offset -= get<2>(muts[i]);
                    break;
                }
        }
    }

    // Count different types of mutations applied and track indel lengths
    int snp_count = 0, insertion_count = 0, deletion_count = 0;
    std::vector<int> insertion_lengths, deletion_lengths;
    
    for (const auto& mut : muts) {
        switch (get<1>(mut)) {
            case 1: snp_count++; break;
            case 2: 
                insertion_count++; 
                insertion_lengths.push_back(get<2>(mut));
                break; 
            case 4: 
                deletion_count++; 
                deletion_lengths.push_back(get<2>(mut));
                break;
        }
    }
    
    // Log the mutations applied with indel lengths
    std::cout << "Applied " << snp_count << " SNPs, " << insertion_count << " insertions";
    if (!insertion_lengths.empty()) {
        std::cout << " (lengths: ";
        for (size_t i = 0; i < insertion_lengths.size(); i++) {
            if (i > 0) std::cout << ",";
            std::cout << insertion_lengths[i];
        }
        std::cout << ")";
    }
    std::cout << ", " << deletion_count << " deletions";
    if (!deletion_lengths.empty()) {
        std::cout << " (lengths: ";
        for (size_t i = 0; i < deletion_lengths.size(); i++) {
            if (i > 0) std::cout << ",";
            std::cout << deletion_lengths[i];
        }
        std::cout << ")";
    }
    std::cout << " to node " << curNode << std::endl;

    faos << '>' << curNode << '\n';
    size_t linesize = 80;
    for (size_t i = 0; i < nseq.size(); i += linesize) {
        faos << nseq.substr(i, std::min(linesize, nseq.size() - i)) << '\n';
    }
    
    vros.close();
    faos.close();
}

struct snpInfo {
  size_t pos;
  char ref;
  char mut;
};

void genMutSNP(
  const std::string& curNode, const std::string& seq, const fs::path& fastaOut, const fs::path& vcfOut,
    const genotyping::mutationMatrices& mutMat, const std::vector<double>& mutnum_double, std::mt19937& gen, std::vector<std::discrete_distribution<>>& distributions, std::vector<char>& bases
) {
  std::string nseq = seq;
  std::vector<snpInfo> snps;
  
  // Check if we should apply any mutations at all
  std::vector<int> num = genMutNum(mutnum_double, gen);
  int total_mutations = num[0]; // Only SNPs for genMutSNP
  
  if (total_mutations > 0) {
    // Apply mutations to random positions
    std::uniform_int_distribution<> pos_dist(0, seq.size() - 1);
    std::set<size_t> mutated_positions;
    
    for (int mut_count = 0; mut_count < total_mutations && mutated_positions.size() < seq.size(); mut_count++) {
      size_t pos;
      do {
        pos = pos_dist(gen);
      } while (mutated_positions.count(pos) > 0); // Avoid mutating the same position twice
      
      mutated_positions.insert(pos);
      
      const char& ref = seq[pos];
      int refidx = getIndexFromNucleotide(ref);
      if (refidx > 3) continue; // Skip invalid nucleotides
      
      int mutidx = distributions[refidx](gen);
      if (mutidx == refidx) continue; // Skip if mutation is same as reference
      
      char mut = bases[mutidx];
      nseq[pos] = mut;
      snps.push_back({pos, ref, mut});
    }
  }

  std::ofstream faos(fastaOut.string(), std::ofstream::trunc);
  std::ofstream vros(vcfOut.string(), std::ofstream::trunc);

  // Log the mutations applied
  std::cout << "Applied " << snps.size() << " SNPs to node " << curNode << std::endl;

  vros << "##fileformat=VCFv4.3\n"
      //  << "##contig=<ID=" + curNode + ">\n"
        << "##contig=<ID=ref>\n"
        << "##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>\n";
  vros << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + fastaOut.stem().string() + "\n";
  for (const auto& snp : snps) {
    vros << "ref\t" << std::to_string(snp.pos + 1) + "\t.\t" << snp.ref << "\t" << snp.mut << "\t.\t.\t.\tGT\t1\n";
  }

  faos << '>' << fastaOut.stem().string() << '\n';
  size_t linesize = 80;
  for (size_t i = 0; i < nseq.size(); i += linesize) {
      faos << nseq.substr(i, std::min(linesize, nseq.size() - i)) << '\n';
  }
  
  vros.close();
  faos.close();
}

void simReads(const fs::path& fastaOut, const fs::path& outReadsObj, const std::string& model, int n_reads, int cpus, unsigned seed) {
    std::string cmd;
    if (model == "perfect") {
        cmd = "iss generate --mode perfect --genomes \'" + fastaOut.string() + "\'"
            + " -n " + std::to_string(n_reads) + " --output \'" + (outReadsObj / fastaOut.stem()).string() + "\'" + " --cpus " + std::to_string(cpus) + " --seed " + std::to_string(seed);
    } else {
        cmd = "iss generate --model " + model + " --genomes \'" + fastaOut.string() + "\'"
            + " -n " + std::to_string(n_reads) + " --output \'" + (outReadsObj / fastaOut.stem()).string() + "\'" + " --cpus " + std::to_string(cpus) + " --seed " + std::to_string(seed);
    }
    std::cout << "iss cmd: " << cmd << std::endl;
    system(cmd.c_str());

    std::string clean = "rm " + (outReadsObj / fs::path("*.iss.tmp*vcf")).string();
    system(clean.c_str());
}

void sim(panmanUtils::Tree* T, const std::string& refNode, const std::string& outDir, const std::string& prefix,
    const std::string& mut_spec_type, const std::vector<double>& num, const std::pair<int, int>& indel_len, const std::string& model,
    int n_reads, int rep, const genotyping::mutationMatrices& mutMat, unsigned seed, int cpus, bool no_reads, const std::string& in_fasta, bool include_internal)
{
    fs::path outDirObj = outDir;
    fs::path outRefFastaObj = outDir / fs::path(prefix + "_refFasta");
    fs::path outVarFastaObj = outDir / fs::path(prefix + "_varFasta");
    fs::path outVCFTrueObj  = outDir / fs::path(prefix + "_vcfTrue");
    fs::path outReadsObj    = outDir / fs::path(prefix + "_reads");

    makeDir(outRefFastaObj.string());
    makeDir(outVarFastaObj.string());
    makeDir(outVCFTrueObj.string());
    makeDir(outReadsObj.string());

    std::mt19937 gen(seed);
    std::vector<char> bases = {'A', 'C', 'G', 'T'};
    std::vector<std::discrete_distribution<>> distributions(4);

    for (size_t i = 0; i < 4; i++) {
      std::vector<double> probs = mutMat.submat[i];
      double sum = 0;
      for (size_t j = 0; j < 4; j++) {
        probs[j] = pow(10, -probs[j] / 10);
        sum += probs[j];
      }
      for (size_t j = 0; j < 4; j++) {
        probs[j] /= sum;
      }
      distributions[i] = std::discrete_distribution<>(probs.begin(), probs.end());
    }
    
    // Log requested mutations for this simulation
    std::cout << "=== SIMULATION PARAMETERS ===" << std::endl;
    std::cout << "Requested mutations: SNPs=" << num[0] << ", Insertions=" << num[1] << ", Deletions=" << num[2] << std::endl;
    std::cout << "Mutation matrix type: " << mut_spec_type << std::endl;
    std::cout << "Number of replicates: " << rep << std::endl;
    std::cout << "=============================" << std::endl;
    
    // Only enumerate nodes if we will need them for random selection
    std::vector<std::string> nodeNames;
    bool will_use_random_selection = (in_fasta.empty() || !fs::exists(in_fasta)) && (refNode == "RANDOM");
    if (will_use_random_selection) {
        for (const auto& pair : T->allNodes) {
            // Check if we should include internal nodes
            if (include_internal || T->allNodes[pair.first]->children.empty()) {
                // Include all nodes if include_internal=true, or only leaves if include_internal=false
                nodeNames.push_back(pair.first);
            }
        }
        std::cout << "Available nodes for random selection: " << nodeNames.size() 
                  << " (include_internal=" << include_internal << ")" << std::endl;
        if (nodeNames.empty()) {
            throw std::runtime_error("No suitable nodes found for random selection");
        }
    }
    
    for (int i = 0; i < rep; i++) {
                std::string curNode;
                std::string curNodeID;
                if (!in_fasta.empty() && fs::exists(in_fasta)) {
                        // Use node name from FASTA header
                        curNodeID = readFastaHeader(in_fasta);
                        curNode = curNodeID;
                        // Sanitize filename by replacing problematic characters
                        std::replace(curNode.begin(), curNode.end(), '/', '_');
                        std::replace(curNode.begin(), curNode.end(), ' ', '_');
                        std::replace(curNode.begin(), curNode.end(), '|', '_');
                        std::replace(curNode.begin(), curNode.end(), ':', '_');
                        std::replace(curNode.begin(), curNode.end(), '?', '_');
                        std::replace(curNode.begin(), curNode.end(), '*', '_');
                        std::replace(curNode.begin(), curNode.end(), '<', '_');
                        std::replace(curNode.begin(), curNode.end(), '>', '_');
                        std::replace(curNode.begin(), curNode.end(), '"', '_');
                        std::cout << "curNodeID (from FASTA): " << curNodeID << std::endl;
                        std::cout << "Making reference fasta for " << curNode << std::endl;
                        std::string refSequence = readFastaSequence(in_fasta);
                        // Make reference fasta - use original curNodeID in header, sanitized curNode for filename
                        makeFasta(curNodeID, refSequence, (outRefFastaObj / fs::path(curNode + ".fa")).string());
                        // Make variant fasta and variant vcf
                        fs::path fastaOut = outVarFastaObj / fs::path(curNode + ".var.fa");
                        fs::path vcfOut  = outVCFTrueObj  / fs::path(curNode + ".var.vcf");
                        if (mut_spec_type == "snp") {
                            genMutSNP(curNodeID, refSequence, fastaOut, vcfOut, mutMat, num, gen, distributions, bases);
                        } else {
                            genMut(curNodeID, refSequence, fastaOut, vcfOut, mutMat, num, indel_len, 500, 500, gen, seed);
                        }
                        // Make reads using InSilicoSeq
                        if (!no_reads) {
                                simReads(fastaOut, outReadsObj, model, n_reads, cpus, seed);
                        }
                } else {
                        // Use random or specified node from panman
                        if (refNode == "RANDOM") {
                                // Generate a different random selection for each replicate
                                std::mt19937 rng(seed + i); // Different seed for each replicate
                                std::uniform_int_distribution<> dis(0, nodeNames.size() - 1);
                                curNode = nodeNames[dis(rng)];
                        } else {
                                curNode = refNode;
                        }
                        curNodeID = curNode;
                        // Sanitize filename by replacing problematic characters
                        std::replace(curNode.begin(), curNode.end(), '/', '_');
                        std::replace(curNode.begin(), curNode.end(), ' ', '_');
                        std::replace(curNode.begin(), curNode.end(), '|', '_');
                        std::replace(curNode.begin(), curNode.end(), ':', '_');
                        std::replace(curNode.begin(), curNode.end(), '?', '_');
                        std::replace(curNode.begin(), curNode.end(), '*', '_');
                        std::replace(curNode.begin(), curNode.end(), '<', '_');
                        std::replace(curNode.begin(), curNode.end(), '>', '_');
                        std::replace(curNode.begin(), curNode.end(), '"', '_');
                        std::cout << "curNodeID: " << curNodeID << std::endl;
                        std::cout << "Making reference fasta for " << curNode << std::endl;
                        std::string refSequence = T->getStringFromReference(curNodeID, false);
                        // Make reference fasta - use original curNodeID in header, sanitized curNode for filename
                        makeFasta(curNodeID, refSequence, (outRefFastaObj / fs::path(curNode + ".fa")).string());
                        // Make variant fasta and variant vcf
                        fs::path fastaOut;
                        fs::path vcfOut;
                        if (refNode == "RANDOM") {
                                fastaOut = outVarFastaObj / fs::path(curNode + ".var.fa");
                                vcfOut  = outVCFTrueObj  / fs::path(curNode + ".var.vcf");
                        } else {
                                fastaOut = outVarFastaObj / fs::path(curNode + ".var." + std::to_string(i) + ".fa");
                                vcfOut  = outVCFTrueObj  / fs::path(curNode + ".var." + std::to_string(i) + ".vcf");
                        }
                        if (mut_spec_type == "snp") {
                            genMutSNP(curNodeID, refSequence, fastaOut, vcfOut, mutMat, num, gen, distributions, bases);
                        } else {
                            genMut(curNodeID, refSequence, fastaOut, vcfOut, mutMat, num, indel_len, 500, 500, gen, seed);
                        }
                        // Make reads using InSilicoSeq
                        if (!no_reads) {
                                simReads(fastaOut, outReadsObj, model, n_reads, cpus, seed);
                        }
                }
    }
}

bool close_to_int(double num_double) {
    int num_int = int(num_double + 0.5);
    double num_int_double = double(num_int);
    if (abs(num_int_double - num_double) < 0.0001) {
        return true;
    }
    return false;
}

int genInt(double num_double, std::mt19937& gen) {
    std::uniform_real_distribution<> dist(0, 1);
    int i = int(num_double);
    double prob_to_add = num_double - double(i);
    if (dist(gen) < prob_to_add) {
        ++i;
    }
    return i;
}

std::vector<int> genMutNum(const std::vector<double>& mutNum_double, std::mt19937& gen) {
    std::vector<int> mutNum_int;
    for (double num : mutNum_double) {
        if (close_to_int(num)) {
            mutNum_int.push_back(num + 0.5);
        } else {
            mutNum_int.push_back(genInt(num, gen));
        }
    }
    return mutNum_int;
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

void scaleMutationMatrices(genotyping::mutationMatrices& mutMat, double mutation_rate) {
  std::vector<std::vector<double>> scaled_submat_phred = mutMat.submat;
  std::vector<std::vector<double>> scaled_submat_prob = phredMatrix2ProbMatrix(scaled_submat_phred);
  double avg_mutation_rate = getAverageMutationRate(scaled_submat_prob);
  double scale_factor = mutation_rate / avg_mutation_rate;

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

  mutMat.submat = probMatrix2PhredMatrix(scaled_submat_prob);
}