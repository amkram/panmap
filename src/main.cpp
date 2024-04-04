#include "kseq.h"
#include "pmi.hpp"
#include "genotype.hpp"
#include "place.hpp"
#include <iostream>
#include <string>
#include <cstdio>
#include <panmat/PangenomeMAT.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

using namespace pmi;
namespace po = boost::program_options;
namespace fs = boost::filesystem;

/* Helpers */
void promptAndPlace(Tree *T, const int32_t k, const int32_t s, const std::string &indexFile, const std::string &pmatFile, std::string &reads1File, std::string &reads2File, std::string &samFileName, std::string &bamFileName, std::string &mpileupFileName, std::string &vcfFileName, std::string &refFileName, const bool prompt) {
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
        string reads2File = "";
        cout << "[Second FASTQ path]: ";
        getline(cin, reads2File);
    }
    if (reads1File == "q" || reads2File == "q") {
        exit(0);
    }
    std::ifstream ifs(indexFile);

    place::placeIsolate(ifs, reads1File, reads2File, samFileName, bamFileName, mpileupFileName, vcfFileName, refFileName, T);

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
            ("params,p", po::value<std::vector<int32_t>>()->multitoken(), "Specify syncmer parameters k,s for indexing (e.g. for k=13, s=8 use -p 13,8)")
            ("f", "Proceed without prompting for confirmation. Applies to index construction if -p is specified or no index is found at ( /path/to/pmat.pmi ). Applies to sample placement if reads are provided.")
            ("s", po::value<std::string>(), "Path to sam output")
            ("b", po::value<std::string>(), "Path to bam output")
            ("m", po::value<std::string>(), "Path to mpileup output")
            ("v", po::value<std::string>(), "Path to vcf output")
            ("r", po::value<std::string>(), "Path to reference fasta output")
        ;

        po::positional_options_description p;
        p.add("panmat", 1);
        //p.add("reads1", 1);
        //p.add("reads2", 1);
        po::variables_map vm;
        po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
        po::notify(vm);
        std::string pmatFile = vm["panmat"].as<std::string>();
        
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

        std::string indexFile = "";
        std::string reads1File = "";
        std::string reads2File = "";
        if (vm.count("reads1")) {
            reads1File = vm["reads1"].as<std::string>();
        }
        if (vm.count("reads2")) {
            reads2File = vm["reads2"].as<std::string>();
        }

        int32_t k = -1;
        int32_t s = -1;
        bool prompt = !vm.count("f");

        if (vm.count("params")) {
            auto k_s_values = vm["params"].as<std::vector<int32_t>>();
            if (k_s_values.size() != 2) {
                std::cerr << "× Error: -r/--reindex requires two integer values (k and s), e.g. -r 13,8." << std::endl;
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
                std::cout << "Use existing index? [Y/n]";
                getline(std::cin, inp);
                if (inp == "Y" || inp == "y" || inp == "") {
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


        //TODO, add some sort of error if all of these are empty, cus then why are we placing
        std::string samFileName = "";
        std::string bamFileName = "";
        std::string mpileupFileName = "";
        std::string vcfFileName = "";
        std::string refFileName = "";
        if (vm.count("s")) {
            samFileName = vm["s"].as<std::string>();
        }
        if (vm.count("b")) {
            bamFileName = vm["b"].as<std::string>();
        }
        if (vm.count("m")) {
            mpileupFileName = vm["m"].as<std::string>();
        }
        if (vm.count("v")) {
            vcfFileName = vm["v"].as<std::string>();
        }
        if (vm.count("r")) {
            refFileName = vm["r"].as<std::string>();
        }

        promptAndPlace(T, k, s, indexFile, pmatFile, reads1File, reads2File, samFileName, bamFileName, mpileupFileName, vcfFileName, refFileName, prompt);

    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}

