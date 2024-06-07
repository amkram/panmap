#include "kseq.h"
#include "PangenomeMAT.hpp"
#include "pmi.hpp"
#include "place.hpp"
#include "mgsr.hpp"
#include "tree.hpp"
#include "seeding.hpp"
#include <iostream>
#include <string>
#include <cstdio>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <tbb/global_control.h>

using namespace pmi;
namespace po = boost::program_options;
namespace fs = boost::filesystem;

/* Helpers */
void promptAndPlace(Tree *T, const tree::mutationMatrices& mutMat, const int32_t k, const int32_t s, const std::string &indexFile, const std::string &pmatFile, std::string &reads1File, std::string &reads2File, const std::string& prefix, const bool& makeSam, const bool& makeBam, const bool& makeMPileup, const bool& makeVCF, const bool& makeRef, const bool prompt, bool use_root) {

    std::cin.clear();
    std::fflush(stdin);
    using namespace std;
    // Check if reads were supplied by user       
    if (reads1File == "" && reads2File == "" && !prompt) {
        cout << "Can't place sample because no reads were provided and -f was specified." << std::endl;
        exit(1);
    }
    // Prompt for fastq files
    if (reads1File == "") {
        while (!reads1File.size()) {
            cout << "First FASTQ path: ";
            getline(cin, reads1File);
        }
        reads2File = "";
        cout << "[Second FASTQ path]: ";
        getline(cin, reads2File);

    }
    if (reads1File == "q" || reads2File == "q") {
        exit(0);
    }
    std::ifstream ifs(indexFile);

    place::placeIsolate(ifs, mutMat, reads1File, reads2File, prefix, makeSam, makeBam, makeMPileup, makeVCF, makeRef, T, use_root);
}

void promptAndIndex(Tree *T, const bool prompt, const std::string &indexFile) {
    using namespace std;
    int32_t k, s, j;
    cout << "Syncmer parameter k: ";
    cin >> k;
    cout << "Syncmer parameter s: ";
    cin >> s;
    cout << "Seed-mer parameter j: ";
    cin >> j;
    cout << "Building (" << k << ", " << s << ", " << j << ") index ..." << std::endl;
    seedIndex index;
    pmi::build(index, T, j, k, s);
    std::cout << "Writing to " << indexFile << "..." << std::endl;
    std::ofstream fout(indexFile);
    fout << index.outStream.str();
    std::cout << "Done." << std::endl;
}

/* Parse command line input and run logic */
int main(int argc, char *argv[]) {
    try {
        po::options_description desc("Options");
        desc.add_options()
            ("panmat", po::value<std::string>()->required(), "Path to tree.pmat file (required)")
            ("reads1", po::value<std::string>(), "Path to first (or single-end) fastq file (optional)")
            ("reads2", po::value<std::string>(), "Path to second paired-end fastq file (optional)")
            ("index,i", po::value<std::string>(), "Path to index file")
            ("mutmat,u", po::value<std::string>(), "Path to mutation matrix file")
            ("params,p", po::value<std::vector<int32_t>>()->multitoken(), "Specify syncmer parameters k,s for indexing (e.g. for k=13, s=8 use -p 13,8)")
            ("f", "Proceed without prompting for confirmation. Applies to index construction if -p is specified or no index is found at ( /path/to/pmat.pmi ). Applies to sample placement if reads are provided.")
            ("window,w", po::value<size_t>()->default_value(20), "Mutation Matrix Build Specific: window size for selecting starting and end sequence positions for building mutation spectrum")
            ("percent_identity", po::value<double>()->default_value(0.80), "Mutation Matrix Build Specific: percent identity threshold for selecting starting and end sequence positions for building mutation spectrum")
            ("prefix", po::value<std::string>(), "output prefix")
            ("s", "produce sam output")
            ("b", "produce bam output")
            ("m", "produce mpileup output")
            ("v", "produce vcf output")
            ("r", "produce reference fasta output")
            /* accio */
            ("accio", "Accio!")
            ("accioParams", po::value<std::vector<int32_t>>()->multitoken(), "Accio parameters --accioParams k s l")
            ("cpus", po::value<size_t>()->default_value(1), "cpus to use")
            /* accio */
            ("g", "Use root as reference, do not place") //FIXME for debugging purposes, remove for release
        ;

        po::positional_options_description p;
        p.add("panmat", 1);
        p.add("reads1", 1);
        p.add("reads2", 1);
        po::variables_map vm;
        po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
        po::notify(vm);

        fs::path pmatFilePathObject = fs::path(vm["panmat"].as<std::string>());
        std::string pmatFile = fs::canonical(pmatFilePathObject).string();
        
        std::cout << "   ┏━┳━●\033[36;1m\033[1m pan \033[0m" << std::endl;
        std::cout << "  ━┫ ┗━━━●\033[36;1m\033[1m map\033[0m\033[32;1m\033[0m" << std::endl;
        std::cout << "   ┗━\033[36;1m v0.0\033[0m ━●\033[32;1m\033[0m" << std::endl;

        std::ifstream ifs(pmatFile);
        boost::iostreams::filtering_streambuf< boost::iostreams::input> b;
        b.push(boost::iostreams::gzip_decompressor());
        b.push(ifs);
        std::istream is(&b);

        auto T = new PangenomeMAT::Tree(is);

        std::cout << "\nUsing tree: \e[3;1m" << pmatFile << "\e[0m  (" << T->allNodes.size() << " nodes)" << std::endl;
        
        std::string prefix     = "";
        std::string indexFile  = "";
        std::string mutmatFile = "";
        std::string reads1File = "";
        std::string reads2File = "";
        if (vm.count("reads1")) {
            reads1File = vm["reads1"].as<std::string>();
        }
        if (vm.count("reads2")) {
            reads2File = vm["reads2"].as<std::string>();
        }
        if (vm.count("prefix")) {
            prefix = vm["prefix"].as<std::string>();
        }

        // TODO, add some sort of error if all of these are empty, cus then why are we placing
        bool makeSam = false;
        bool makeBam = false;
        bool makeMPileup = false;
        bool makeVCF = false;
        bool makeRef = false;
        if (vm.count("s")) {
            makeSam = true;
        }
        if (vm.count("b")) {
            makeBam = true;
        }
        if (vm.count("m")) {
            makeMPileup = true;
        }
        if (vm.count("v")) {
            makeVCF = true;
        }
        if (vm.count("r")) {
            makeRef = true;
        }

        // Mutation matrix (mm) file
        tree::mutationMatrices mutMat = tree::mutationMatrices();
        std::string defaultMutmatPath = pmatFile + ".mm";
        if (vm.count("mutmat")) {
            mutmatFile = fs::canonical(fs::path(vm["mutmat"].as<std::string>())).string();
            std::ifstream mminf(mutmatFile);
            tree::fillMutationMatricesFromFile(mutMat, mminf);
            mminf.close();
            std::cout << "Using mutation matrix file: " << mutmatFile << std::endl;
        } else if (fs::exists(defaultMutmatPath)) {
            std::ifstream mminf(defaultMutmatPath);
            tree::fillMutationMatricesFromFile(mutMat, mminf);
            mminf.close();
            std::cout << "Mutation matrix file detected, using mutation matrix file: " << defaultMutmatPath << std::endl;
        } else {
            std::cout << "No mutation matrix file detected, building mutation matrix ..." << std::endl; 
            // build mutation matrix
            tree::fillMutationMatricesFromTree(mutMat, T, vm["window"].as<size_t>(), vm["percent_identity"].as<double>());
            // write to file
            std::cout << "Writing to " << defaultMutmatPath << " ..." << std::endl;
            std::ofstream mmfout(defaultMutmatPath);
            tree::writeMutationMatrices(mutMat, mmfout);
            mmfout.close();
            mutmatFile = defaultMutmatPath;
        }

        if (vm.count("accio")) {
            /* accio */
            tbb::global_control c(tbb::global_control::max_allowed_parallelism, vm["cpus"].as<size_t>());            
            int32_t accioK = 10;
            int32_t accioS = 5;
            int32_t accioL = 3;
            if (vm.count("accioParams")) {
                auto accioParamValues = vm["accioParams"].as<std::vector<int32_t>>();
                int32_t accioK = accioParamValues[0];
                int32_t accioS = accioParamValues[1];
                int32_t accioL = accioParamValues[2];
            }
            std::string defaultKmiPath = pmatFile + "." + std::to_string(accioK) + "_" + std::to_string(accioS) + "_" + std::to_string(accioL) + ".kmi";
            std::string defaultSmiPath = pmatFile + "." + std::to_string(accioK) + "_" + std::to_string(accioS) + ".smi";
            if (!fs::exists(defaultKmiPath)) {
                std::cerr << "kmi file not detected... building kmi file using input parameters" << std::endl;
                pmi::seedIndex smiIndex;
                std::stringstream seedmersOutStream;
                mgsr::buildSeedmer(smiIndex, T, accioL, accioK, accioS, seedmersOutStream);
                std::cerr << "Writing to " << defaultKmiPath << "..." << std::endl;
                std::ofstream smfout(defaultKmiPath);
                smfout << seedmersOutStream.str();
                std::cerr << "Writing to " << defaultSmiPath << "..." << std::endl;
                std::ofstream fout(defaultSmiPath);
                fout << smiIndex.outStream.str();
            }
            int maximumGap = 10;
            int minimumCount = 0;
            int minimumScore = 0;
            double errorRate = 0.005;
            bool confidence = false;
            int32_t roundsRemove = 2;
            double removeThreshold = 0.005;
            place::placeAccio(T, mutMat, accioK, accioS, accioL, defaultKmiPath, reads1File, reads2File, prefix, makeSam, makeBam, makeMPileup, makeVCF, makeRef,
                maximumGap, minimumCount, minimumScore, errorRate, confidence, roundsRemove, removeThreshold);

            /* accio */
        } else {
            int32_t k = -1;
            int32_t s = -1;
            bool prompt = !vm.count("f");
            bool use_root = vm.count("g");

            if (vm.count("params")) {
                auto k_s_values = vm["params"].as<std::vector<int32_t>>();
                if (k_s_values.size() != 2) {
                    std::cerr << " Error: -r/--reindex requires two integer values (k and s), e.g. -r 13,8." << std::endl;
                    return 1;
                } else {
                    k = k_s_values[0];
                    s = k_s_values[1];
                }
            } else if (vm.count("index")) {
                indexFile = vm["index"].as<std::string>();
                std::cout << "Using seed index: \e[3;1m" << indexFile << "\e[0m" << std::endl;          
            } else {
                bool keep = false;
                std::string defaultIndexPath = pmatFile + ".pmi";
                std::string inp;
                if (boost::filesystem::exists(defaultIndexPath)) {
                    if (prompt){
                        std::cout << "Use existing index? [Y/n]";
                        getline(std::cin, inp);
                        if (inp == "Y" || inp == "y" || inp == "") {
                            keep = true;
                        }
                    }else{
                        keep = true;
                    }
                } else {
                    std::cout << "\nNo index found." << std::endl;
                }
                if (!keep) {
                    promptAndIndex(T, prompt, defaultIndexPath);
                }
                indexFile = defaultIndexPath;
            }

            promptAndPlace(T, mutMat, k, s, indexFile, pmatFile, reads1File, reads2File, prefix, makeSam, makeBam, makeMPileup, makeVCF, makeRef, prompt, use_root);

        }

    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}

