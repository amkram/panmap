#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/lzma.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <algorithm>
#include <iostream>
#include <vector>
#include <random>
#include <functional>
#include "panmanUtils.hpp"
#include "../seed_annotated_tree.hpp"


namespace po = boost::program_options;
namespace fs = boost::filesystem;

void makeFasta(const std::string& name, const std::string& seq, const std::string& path);
void makeDir(const std::string& path);
std::vector<int> genMutNum(const std::vector<double>& mutNum_double, std::mt19937& gen);
void sim(panmanUtils::Tree* T, const std::string& refNode, const std::string& out_dir, const std::string& prefix,
  const std::string& mut_spec_type, const std::vector<double>& num, const std::pair<int, int>& indel_len, const std::string& model,
  int n_reads, int rep, const seed_annotated_tree::mutationMatrices& mutMat, unsigned seed, int cpus, bool no_reads, bool leaf_node_only, bool sim_ref_reads);
void scaleMutationMatrices(seed_annotated_tree::mutationMatrices& mutMat, double mutation_rate);

int main(int argc, char *argv[]) {
    std::cout << "What is my purpose?\nYou pass butter" << std::endl;
    try {
        po::options_description desc("Options");
        desc.add_options()
            ("help,h", "Produce help message")
            ("panmat",    po::value<std::string>()->required(), "Path to tree.pmat file (required).")
            ("ref",       po::value<std::string>()->required(), "Reference node name (required). Enter RANDOM to use a random node on the tree, and a new random node is selected without replacement for each replicate.")
            ("out_dir",   po::value<std::string>()->required(), "Output directory (required)")
            ("prefix",    po::value<std::string>()->default_value("prefix"), "Output prefix")
            ("mutnum",    po::value<std::vector<double>>()->multitoken(), "Number of mutations for snp, insertion, and deletion [10 0 0].")
            ("indel_len", po::value<std::vector<int>>()->multitoken(), "Min and max indel length [1 9]. Uniform distribution.")
            ("mut_spec",  po::value<std::string>()->default_value(""), "Use input mutation matrix file to model mutations")
            ("mut_spec_type", po::value<std::string>()->default_value(""), "Type of mutation matrix to use. Options: snp, indel, both. Default: none. Currently only supports snp. When using snp, the number of mutations and reference nucleotide to mutate are also modeled by the mutation matrix.")
            ("mutation_rate", po::value<double>()->default_value(-1), "Mutation rate for snp. Default: no scaling")
            ("rep",       po::value<int>()->default_value(1), "Number of replicates to simulate [1].")
            ("n_reads",   po::value<int>()->default_value(2000), "Number of reads to simulate [2000].")
            ("model",     po::value<std::string>()->default_value("NovaSeq"), "InSilicoSeq error model [HiSeq]. Options: HiSeq, NextSeq, NovaSeq, MiSeq. For detail, visit InSilicoSeq github (https://github.com/HadrienG/InSilicoSeq).")
            ("leaf-node-only", po::bool_switch()->default_value(false), "Only simulate reads on leaf nodes")
            ("sim-ref-reads", po::bool_switch()->default_value(false), "Simulate reads on the reference node")
            ("no-reads",  po::bool_switch()->default_value(false), "Do not simulate reads")
            ("cpus",      po::value<int>()->default_value(1), "Number of CPUs to use [1].")
            ("seed",      po::value<std::string>()->default_value("RANDOM"), "Random seed for simulation [default: random]")
        ;

        po::positional_options_description p;
        p.add("panmat", 1);
        p.add("ref", 1);

        po::variables_map vm;
        po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
        if (vm.count("help")) {
            std::cout << desc << "\n";
            return 0;
        }
        po::notify(vm);
        
        // input variables
        std::string panmatPath    = vm["panmat"].as<std::string>();
        std::string mut_spec      = vm["mut_spec"].as<std::string>();
        std::string mut_spec_type = vm["mut_spec_type"].as<std::string>();
        std::string out_dir       = vm["out_dir"].as<std::string>();
        std::string refNode       = vm["ref"].as<std::string>();
        std::string prefix        = vm["prefix"].as<std::string>();
        std::string model         = vm["model"].as<std::string>();
        std::string seedstr       = vm["seed"].as<std::string>();
        double mutation_rate      = vm["mutation_rate"].as<double>();
        int n_reads               = vm["n_reads"].as<int>();
        int rep                   = vm["rep"].as<int>();
        int cpus                  = vm["cpus"].as<int>();
        bool no_reads             = vm["no-reads"].as<bool>();
        bool leaf_node_only       = vm["leaf-node-only"].as<bool>();
        bool sim_ref_reads        = vm["sim-ref-reads"].as<bool>();
        // Check mut_spec
       seed_annotated_tree::mutationMatrices mutMat = seed_annotated_tree::mutationMatrices();
        if (out_dir != "") {
            if (fs::exists(mut_spec)) {
                std::ifstream mminf(mut_spec);
               seed_annotated_tree::fillMutationMatricesFromFile(mutMat, mminf);
                mminf.close();
            } else {
                throw std::invalid_argument("--mut_sepc input file doesn't exist");
            }
        }

        if (mutation_rate != -1) {
          scaleMutationMatrices(mutMat, mutation_rate);
        }

        // Check out_dir input
        if (!fs::is_directory(fs::path(out_dir))) {
            throw std::invalid_argument("--out_dir input s not a valid directory");
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
                << "No reads: " << no_reads << "\n"
                << "Leaf node only: " << leaf_node_only << "\n";

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
                throw std::invalid_argument("--num must have 3 inputs");
            }
        }
        if (mut_spec_type != "" && mut_spec_type != "snp") {
            throw std::invalid_argument("--mut_spec_type only supports snp");
        }

        // check model input
        const std::vector<std::string> acceptedModels = {"NextSeq", "NovaSeq", "HiSeq", "MiSeq"};
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
        sim(T, refNode, out_dir, prefix, mut_spec_type, mutnum_double, indel_len, model, n_reads, rep, mutMat, seed, cpus, no_reads, leaf_node_only, sim_ref_reads);

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


char subNuc(char ref, const seed_annotated_tree::mutationMatrices& mutMat, std::mt19937& gen) {
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
size_t genLen(const std::pair<int, int>& indel_len, const seed_annotated_tree::mutationMatrices& mutMat, int type, std::mt19937& gen) {
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

}

void genMut(const std::string& curNode, const std::string& seq, const fs::path& fastaOut, const fs::path& vcfOut,
    const seed_annotated_tree::mutationMatrices& mutMat, const std::vector<double> mutnum_double, const std::pair<int, int> indel_len,
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

    faos << '>' << fastaOut.stem().string() << '\n';
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
  const seed_annotated_tree::mutationMatrices& mutMat, std::mt19937& gen, std::vector<std::discrete_distribution<>>& distributions, std::vector<char>& bases
) {
  std::string nseq = seq;
  std::vector<snpInfo> snps;
  for (size_t i = 0; i < seq.size(); i++) {
    const char& ref = seq[i];
    int refidx = getIndexFromNucleotide(ref);
    if (refidx > 3) continue;
    int mutidx = distributions[refidx](gen);
    if (mutidx == refidx) continue;
    char mut = bases[mutidx];
    nseq[i] = mut;
    snps.push_back({i, ref, mut});
  }

  std::ofstream faos(fastaOut.string(), std::ofstream::trunc);
  std::ofstream vros(vcfOut.string(), std::ofstream::trunc);


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
    std::string cmd = "iss generate --model " + model + " --genomes \'" + fastaOut.string() + "\'"
        + " -n " + std::to_string(n_reads) + " --output \'" + (outReadsObj / fastaOut.stem()).string() + "\'" + " --cpus " + std::to_string(cpus) + " --seed " + std::to_string(seed);
    std::cout << "iss cmd: " << cmd << std::endl;
    system(cmd.c_str());

    std::string clean = "rm " + (outReadsObj / fs::path("*.iss.tmp*vcf")).string();
    system(clean.c_str());
}

void sim(panmanUtils::Tree* T, const std::string& refNode, const std::string& outDir, const std::string& prefix,
    const std::string& mut_spec_type, const std::vector<double>& num, const std::pair<int, int>& indel_len, const std::string& model,
    int n_reads, int rep, const seed_annotated_tree::mutationMatrices& mutMat, unsigned seed, int cpus, bool no_reads, bool leaf_node_only, bool sim_ref_reads)
{
    fs::path outDirObj = outDir;
    fs::path outRefFastaObj = outDir / fs::path(prefix + "_refFasta");
    fs::path outVarFastaObj = outDir / fs::path(prefix + "_varFasta");
    fs::path outVCFTrueObj  = outDir / fs::path(prefix + "_vcfTrue");
    fs::path outReadsObj    = outDir / fs::path(prefix + "_reads");
    fs::path outRefReadsObj = outDir / fs::path(prefix + "_refReads");

    makeDir(outRefFastaObj.string());
    makeDir(outVarFastaObj.string());
    makeDir(outVCFTrueObj.string());
    makeDir(outReadsObj.string());

    if (sim_ref_reads) {
      makeDir(outRefReadsObj.string());
    }

    std::vector<std::string> nodeNames;
    if (refNode == "RANDOM") {
        for (const auto& pair : T->allNodes) {
          if (leaf_node_only) {
            if (pair.second->children.empty()) {
              nodeNames.push_back(pair.first);
            }
          } else {
            nodeNames.push_back(pair.first);
          }
        }
        std::default_random_engine rng(seed);
        std::shuffle(begin(nodeNames), end(nodeNames), rng);
    }

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
    
    for (int i = 0; i < rep; i++) {
      std::string curNode;
      if (refNode == "RANDOM") {
          curNode = nodeNames[i];
      } else {
          curNode = refNode;
      }

      std::string curNodeID = curNode;
      std::replace(curNode.begin(), curNode.end(), '/', '_');
      std::replace(curNode.begin(), curNode.end(), ' ', '_');
      std::cout << "curNodeID: " << curNodeID << std::endl;
      std::cout << "Making reference fasta for " << curNode << std::endl;
      // Make reference fasta
      makeFasta(curNode, T->getStringFromReference(curNodeID, false), (outRefFastaObj / fs::path(curNode + ".fa")).string());

      // Make variant fasta and vairant vcf
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
        genMutSNP(curNode, T->getStringFromReference(curNodeID, false), fastaOut, vcfOut, mutMat, gen, distributions, bases);
      } else {
        genMut(curNode, T->getStringFromReference(curNodeID, false), fastaOut, vcfOut, mutMat, num, indel_len, 500, 500, gen, seed);
      }

      // Make reads using InSilicoSeq
      if (!no_reads) {
          simReads(fastaOut, outReadsObj, model, n_reads, cpus, seed);
          if (sim_ref_reads) {
            fs::path refFastaOut = outRefFastaObj / fs::path(curNode + ".fa");
            simReads(refFastaOut, outRefReadsObj, model, n_reads, cpus, seed);
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

void scaleMutationMatrices(seed_annotated_tree::mutationMatrices& mutMat, double mutation_rate) {
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