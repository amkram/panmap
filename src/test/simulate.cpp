#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <algorithm>
#include <iostream>
#include <vector>
#include <random>
#include "panmanUtils.hpp"
#include "../tree.hpp"


namespace po = boost::program_options;
namespace fs = boost::filesystem;

void makeFasta(const std::string& name, const std::string& seq, const std::string& path);
void makeDir(const std::string& path);
std::vector<int> genMutNum(const std::vector<double>& mutNum_double);
void sim(panmanUtils::Tree* T, const std::string& refNode, const std::string& out_dir, const std::string& prefix,
    const std::vector<double>& num, const std::pair<int, int>& indel_len, const std::string& model,
    int n_reads, int rep, const seed_annotated_tree::mutationMatrices& mutMat);

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
            ("rep",       po::value<int>()->default_value(1), "Number of replicates to simulate [1].")
            ("n_reads",   po::value<int>()->default_value(2000), "Number of reads to simulate [2000].")
            ("model",     po::value<std::string>()->default_value("NovaSeq"), "InSilicoSeq error model [HiSeq]. Options: HiSeq, NextSeq, NovaSeq, MiSeq. For detail, visit InSilicoSeq github (https://github.com/HadrienG/InSilicoSeq).")
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
        std::string panmatPath = vm["panmat"].as<std::string>();
        std::string mut_spec   = vm["mut_spec"].as<std::string>();
        std::string out_dir    = vm["out_dir"].as<std::string>();
        std::string refNode    = vm["ref"].as<std::string>();
        std::string prefix     = vm["prefix"].as<std::string>();
        std::string model      = vm["model"].as<std::string>();
        int n_reads            = vm["n_reads"].as<int>();
        int rep                = vm["rep"].as<int>();
        
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


        // Check out_dir input
        if (!fs::is_directory(fs::path(out_dir))) {
            throw std::invalid_argument("--out_dir input s not a valid directory");
        }

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

        // check model input
        const std::vector<std::string> acceptedModels = {"NextSeq", "NovaSeq", "HiSeq", "MiSeq"};
        if (std::find(acceptedModels.begin(), acceptedModels.end(), model) == acceptedModels.end()) {
            throw std::invalid_argument("Unknown --model input, please double check");
        }


        // read treeeee
        std::ifstream ifs(panmatPath);
        boost::iostreams::filtering_streambuf<boost::iostreams::input> b;
        b.push(boost::iostreams::gzip_decompressor());
        b.push(ifs);
        std::istream is(&b);
        panmanUtils::Tree* T = new panmanUtils::Tree(is);

        // check node
        if (refNode != "RANDOM" && T->allNodes.find(refNode) == T->allNodes.end()) {
            throw std::invalid_argument("Couldn't find --ref node on tree");
        }
        // GO time
        sim(T, refNode, out_dir, prefix, mutnum_double, indel_len, model, n_reads, rep, mutMat);

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

char getRandomChar(const std::vector<char>& charList) {
    std::random_device rd;
    std::mt19937 gen(rd());
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

char getRandomCharWithWeights(const std::vector<char>& chars, const std::vector<int>& weights) {
    // Create a random device and generator
    std::random_device rd;
    std::mt19937 gen(rd());

    // Create a discrete distribution with the given weights
    std::discrete_distribution<> distr(weights.begin(), weights.end());

    // Get a random index based on the weights
    int index = distr(gen);

    // Return the character at the random index
    return chars[index];
}

char subNuc(char ref, const seed_annotated_tree::mutationMatrices& mutMat) {
    std::vector<char> bases = {'A', 'C', 'G', 'T'};
    int refIdx = getIndexFromNucleotide(ref);
    if (refIdx > 3) {
        return getRandomChar(bases);
    }

    if (!mutMat.filled) {
        bases.erase(bases.begin() + refIdx);
        return getRandomChar(bases);
    }

    std::vector<double> probs = mutMat.submat[refIdx];
    bases.erase(bases.begin() + refIdx);
    probs.erase(probs.begin() + refIdx);
    
    std::vector<int> wgts;
    double minProb = getMinDouble(probs);
    for (const auto& prob : probs) {
        double deci = pow(10, (minProb - prob) / 10);
        wgts.push_back(deci * 1000);
    }

    return getRandomCharWithWeights(bases, wgts);
}

size_t genLen(const std::pair<int, int>& indel_len, const seed_annotated_tree::mutationMatrices& mutMat, int type) {
    std::random_device rd;
    std::mt19937 gen(rd());

    if (!mutMat.filled) {
        std::uniform_int_distribution<> distribLen(indel_len.first, indel_len.second);
        return distribLen(gen);
    }

    std::vector<double> probs;
    switch (type) {
        case 2:
            probs = mutMat.insmat;
            break;
        case 4:
            probs = mutMat.delmat;
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
    size_t beg, size_t end)
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
    std::vector<int> num = genMutNum(mutnum_double);
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
    std::random_device rd;
    std::mt19937 gen(rd());
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
            length = genLen(indel_len, mutMat, 2);
        } else {
            length = genLen(indel_len, mutMat, 4);
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
                    char mut  = subNuc(ref, mutMat);
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
                        inss += getRandomChar(bases);
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

void simReads(const fs::path& fastaOut, const fs::path& outReadsObj, const std::string& model, int n_reads) {
    std::string cmd = "iss generate --model " + model + " --genomes " + fastaOut.string()
        + " -n " + std::to_string(n_reads) + " --output " + (outReadsObj / fastaOut.stem()).string() + " --cpus 2";
    system(cmd.c_str());

    std::string clean = "rm " + (outReadsObj / fs::path("*.iss.tmp*vcf")).string();
    system(clean.c_str());
}

void sim(panmanUtils::Tree* T, const std::string& refNode, const std::string& outDir, const std::string& prefix,
    const std::vector<double>& num, const std::pair<int, int>& indel_len, const std::string& model,
    int n_reads, int rep, const seed_annotated_tree::mutationMatrices& mutMat)
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

    std::vector<std::string> nodeNames;
    if (refNode == "RANDOM") {
        for (const auto& pair : T->allNodes) {
            nodeNames.push_back(pair.first);
        }
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        std::default_random_engine rng(seed);
        std::shuffle(begin(nodeNames), end(nodeNames), rng);
    }

    for (int i = 0; i < rep; i++) {
        std::string curNode;
        if (refNode == "RANDOM") {
            curNode = nodeNames[i];
        } else {
            curNode = refNode;
        }
        // Make reference fasta
        makeFasta(curNode, T->getStringFromReference(curNode, false), (outRefFastaObj / fs::path(curNode + ".fa")).string());

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
        genMut(curNode, T->getStringFromReference(curNode, false), fastaOut, vcfOut, mutMat, num, indel_len, 500, 500);

        // Make reads using InSilicoSeq
        simReads(fastaOut, outReadsObj, model, n_reads);
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

int genInt(double num_double) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dist(0, 1);
    int i = int(num_double);
    double prob_to_add = num_double - double(i);
    if (dist(gen) < prob_to_add) {
        ++i;
    }
    return i;
}

std::vector<int> genMutNum(const std::vector<double>& mutNum_double) {
    std::vector<int> mutNum_int;
    for (double num : mutNum_double) {
        if (close_to_int(num)) {
            mutNum_int.push_back(num + 0.5);
        } else {
            mutNum_int.push_back(genInt(num));
        }
    }
    return mutNum_int;
}