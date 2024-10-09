#include "genotype.hpp"
#include "place.hpp"
#include "pmi.hpp"
#include "seed_annotated_tree.hpp"
#include <panmanUtils.hpp>
#include <boost/filesystem.hpp>
#include <boost/iostreams/filter/lzma.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/program_options.hpp>
#include <csignal>
#include <cstdio>
#include <iostream>
#include <minimap2/kseq.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <tbb/global_control.h>
#include "docopt.h"
#include "index.capnp.h"
#include <capnp/message.h>
#include <capnp/serialize-packed.h>
#include <fcntl.h>
#include <chrono>
#include <thread>
// #include <nlohmann/json.hpp>
#include <cstdlib>
#include <thread>
#include <future>
#include <iostream>


using namespace pmi;
namespace po = boost::program_options;
namespace fs = boost::filesystem;
namespace sat = seed_annotated_tree;

static const char USAGE[] = 
R"(panmap -- v0.0 ⟗ ⟗

Assemble haplotypes from reads and a pangenome guide tree.

Default outputs (see -o for more options):
    - BAM alignment to the reference selected by panmap (panmap.bam)
    - VCF of variant calls (panmap.vcf)
    - FASTA of the assembled consensus haplotype (panmap.fasta)
    - Log of the panmap run (panmap.log)

Usage:
  panmap <guide> [<reads1.fastq>] [<reads2.fastq>] [options]

<guide>         Path to pangenome reference tree (pan-MAT or pan-MAN file).
<reads1.fastq>  [Optional] Path to first FASTQ file.
<reads2.fastq>  [Optional] Path to second FASTQ file (for paired-end reads).

Preset options:
  --fast                    Equivalent to --subsample-reads 5.0 --subsample-seeds 0.8
  --superfast               Equivalent to --subsample-reads 0.5 --subsample-seeds 0.5
  --accurate                Equivalent to --subsample-reads 100.0 --subsample-seeds 1.0

Input/output options:
  -p --prefix <prefix>      Prefix for output files. [default: panmap]
  -o <outputs>              List of outputs to generate.
                              Accepted values:
                                  placement / p:  Save phylogenetic placement results to <prefix>.placement
                                  assembly / a:   Save consensus assembly to <prefix>.assembly.fa
                                  reference / r:  Save panmap's selected reference to <prefix>.reference.fa
                                  spectrum / c:   Save mutation spectrum matrix to <prefix>.mm
                                  mpileup / m:    Save read pileup to <prefix>.mpileup
                                  sam / s:        Save aligned reads to <prefix>.sam
                                  bam / b:        Save aligned reads to <prefix>.bam
                                  vcf / v:        Save variant calls and likelihoods to <prefix>.vcf
                                  all / A:        Save all possible outputs, each to <prefix>.<ext>
                              [default: bam,vcf,assembly]

  -i --index <path>         Provide a precomputed panmap index. If not specified, index
                                is loaded from <guide.pmat>.pmi, if it exists, otherwise
                                it is built with parameters <k> and <s> (see below). [default: ]
  
  -m --mutmat <path>        Provide a precomputed mutation spectrum matrix instead of computing
                                one from the tree. if not specified, one is computed from the tree.
                                Overrides --prior. [default: ]

Seeding/alignment options:
  -k <k>                             Length of k-mer seeds. [default: 19]
  -s <s>                             Length of s-mers for seed (syncmer) selection. [default: 8]
  -R --subsample-reads <coverage>    Subsample reads for placement to approximately <coverage> depth. [default: 20.0]
  -S --subsample-seeds <proportion>  Use only <proportion> of total seeds for placement. [default: 1.0]
  -P --prior                         Compute and use a mutation spectrum prior for genotyping.
  -f --reindex                       Don't load index from disk, build it from scratch.

Other options:
<<<<<<< HEAD
  -r --place-per-read        placement per read (panmama).
=======
>>>>>>> 1c2118f82a824df6281ca487b72e9982fb396f47
  -c --cpus <num>            Number of CPUs to use. [default: 1]
  -x --stop-after <stage>    Stop after the specified stage. Accepted values:
                                    indexing / i:   Stop after seed indexing
                                    placement / p:  Stop after placement
                                    mapping / m:    Stop after read mapping
                                    genotyping / g: Stop after computing genotype likelihoods
                                    assembly / a:   Stop after consensus assembly
                                 [default: ]
  -Q <seed>                 Integer seed for random number generation. [default: 42]
  -V --version              Show version.
  -h --help                 Show this screen.
  -D --dump                 Dump all seeds to file.
  -X --dump-real            Dump true seeds to file.

Placement per read options:
  --place-per-read       Place reads per read (panmama)

  --maximum-gap <gap>              Maximum gap between matches. [default: 10]
  --minimum-count <count>          Minimum count of seeds in a match. [default: 0]
  --minimum-score <score>          Minimum score of seeds in a match. [default: 0]
  --error-rate <rate>              Error rate of kminmer [default: 0.005]

  --redo-read-threshold           Re-chain reads when the number of kminmers need to be updated exceeds this threshold [default: 5]
  --recalculate-score             Recalculate chain scores for every read at each node (Sometimes, although rare, a
                                    read's kminmer matches are not affected but the coordinate of the kminmer changes,
                                    which can affect the chaining score. This typically affects the scores very slightly.)
  --rescue-duplicates             Rescue duplicate seeds.
  --rescue-duplicates-threshold   Rescue reads with number of duplicates not greater than <threshold>. [default: 5]

  --filter-round <round>          Maximum number of rounds to filter low probability haplotypes before EM. [default: 5]
  --check-frequency <frequency>   Number of iterations between each filter-round. [default: 20]
  --remove-iteration <count>      Remove haplotypes that has probability less than insig-prob for more than remove-count consecutive iterations. [default: 20]
  --insig-prob <prob>             As described in --remove-iteration. (default is calculated from total number of nodes, where default = 1 / (total number of nodes * 10))

  --rounds-remove                 Number of rounds to clean up and remove low probability haplotypes after EM. [default: 3]
  --remove-threshold              Remove haplotypes with probability less than this threshold during rounds-remove. [default: 0.005]

  --leaf-nodes-only               Only consider leaf nodes when placing reads.
)";


using namespace std;


void writeCapnp(::capnp::MallocMessageBuilder &message, std::string &filename) {
  int fd = open(filename.c_str(), O_WRONLY | O_CREAT, 0644);

  if (fd < 0) {
    perror("Failed to open proto file for writing");
    return;
  }

  try {
    capnp::writePackedMessageToFd(fd, message);
  } catch (const std::exception &e) {
    std::cerr << "Failed to write Cap'n Proto message: " << e.what() << std::endl;
  }

  close(fd);
}

std::unique_ptr<::capnp::PackedFdMessageReader> readCapnp(std::string &filename) {
  int fd = open(filename.c_str(), O_RDONLY);
  ::capnp::ReaderOptions options = {(uint64_t) -1, 64};
  auto message = std::make_unique<::capnp::PackedFdMessageReader>(fd, options);
  return message;
}

panmanUtils::Tree* loadPanmanOrPanmat(const std::string &pmatFile) {
    // If the data structure loaded into memory is a PanMAT, it is pointed to by T
    panmanUtils::Tree *T = nullptr;
    // If the data structure loaded into memory is a PanMAN, it is pointed to by TG
    panmanUtils::TreeGroup *TG = nullptr;

    try {
        // Try to load PanMAN file directly into memory
        std::ifstream inputFile(pmatFile);
        boost::iostreams::filtering_streambuf< boost::iostreams::input> inPMATBuffer;
        inPMATBuffer.push(boost::iostreams::lzma_decompressor());
        inPMATBuffer.push(inputFile);
        std::istream inputStream(&inPMATBuffer);


        TG = new panmanUtils::TreeGroup(inputStream);

        
        inputFile.close();
    } catch (const std::exception &e) {
        std::cerr << "Attempting to load as PanMAT...\n";
        try {
            std::ifstream inputFile(pmatFile);
            boost::iostreams::filtering_streambuf< boost::iostreams::input> inPMATBuffer;
            inPMATBuffer.push(boost::iostreams::lzma_decompressor());
            inPMATBuffer.push(inputFile);
            std::istream inputStream(&inPMATBuffer);


            T = new panmanUtils::Tree(inputStream);


        } catch (const std::exception &e) {
            return nullptr;
        }
    }
    if (TG != nullptr) {
      return &(TG->trees[0]);
    }
    return T;
}

void log(const std::string& message) {
    std::ofstream logFile("panmap.log", std::ios_base::app);
    logFile << message << std::endl;
    std::cout << message << std::endl;
}

int main(int argc, const char** argv) {

    std::map<std::string, docopt::value> args = docopt::docopt(USAGE, { argv + 1, argv + argc }, true, "panmap 0.0");
    tbb::global_control c(tbb::global_control::max_allowed_parallelism, std::stoi(args["--cpus"].asString()));
    std::string guide = args["<guide>"].asString();
    std::string reads1 = args["<reads1.fastq>"] ? args["<reads1.fastq>"].asString() : "";
    std::string reads2 = args["<reads2.fastq>"] ? args["<reads2.fastq>"].asString() : "";
    std::string prefix = args["--prefix"] ? args["--prefix"].asString() : "panmap";
    std::string outputs = args["-o"] && args["-o"].isString() ? args["-o"].asString() : "bam,vcf,assembly";



    std::vector<std::string> outputs_seperated;

    std::string token;
    std::stringstream ss(outputs);
    while (std::getline(ss, token, ',')) {
        outputs_seperated.push_back(token);
    }

    std::string refFileName = "";
    std::string samFileName = "";
    std::string bamFileName = "";
    std::string mpileupFileName = ""; 
    std::string vcfFileName = "";

    
    for (const auto& output : outputs_seperated) {
        if(output.size() == 1){
          switch(output[0]){
            case 'r':
              refFileName = prefix + ".reference.fa";
              break;
            case 's':
              samFileName = prefix + ".sam";
              break;
            case 'b':
              bamFileName = prefix + ".bam";
              break;
            case 'm':
              mpileupFileName = prefix + ".mpileup";
              break;
            case 'v':
              vcfFileName = prefix + ".vcf";
              break;
            case 'A':
              refFileName = prefix + ".reference.fa";
              samFileName = prefix + ".sam";
              bamFileName = prefix + ".bam";
              mpileupFileName = prefix + ".mpileup";
              vcfFileName = prefix + ".vcf";
              break;
          }
        }else{
          
            if(output == "reference"){
              refFileName = prefix + ".reference.fa";
            }else if(output == "sam"){
              samFileName = prefix + ".sam";
            }else if(output == "bam"){
              bamFileName = prefix + ".bam";
            }else if(output == "mpileup"){
              mpileupFileName = prefix + ".mpileup";
            }else if(output == "vcf"){
              vcfFileName = prefix + ".vcf";
            }else if(output == "all"){
              refFileName = prefix + ".reference.fa";
              samFileName = prefix + ".sam";
              bamFileName = prefix + ".bam";
              mpileupFileName = prefix + ".mpileup";
              vcfFileName = prefix + ".vcf";
            }
        }
    }


    double subsample_reads = std::stod(args["--subsample-reads"].asString());
    double subsample_seeds = std::stod(args["--subsample-seeds"].asString());
    bool reindex = args["--reindex"] && args["--reindex"].isBool() ? args["--reindex"].asBool() : false;
    bool prior = args["--prior"] && args["--prior"].isBool() ? args["--prior"].asBool() : false;
    bool placement_per_read = args["--place-per-read"] && args["--place-per-read"].isBool() ? args["--place-per-read"].asBool() : false;
  
    int k = std::stoi(args["-k"].asString());
    int s = std::stoi(args["-s"].asString());
    std::string index_path = args["--index"] ? args["--index"].asString() : "";
    std::cout << std::endl;
    std::cout << "   ┏━┳━● pan" << std::endl;
    std::cout << "  ━┫ ┗━━━● map" << std::endl;
    std::cout << "   ┗━ v0.0 ━●" << std::endl << std::endl;

    panmanUtils::Tree *T = loadPanmanOrPanmat(guide);
    if (T == nullptr) {
      std::cerr << "Failed to load guide panMAN/panMAT.\n";
      return 1;
    }


    log("--- Settings ---");
    log("Pangenome: " + guide + " (" + std::to_string(T->allNodes.size()) + " nodes)");
    if (!reads1.empty()) {
      log("Reads: " + reads1 + (reads2.empty() ? "" : " + " + reads2));
    } else {
      log("Reads: <none>");
    }
    log("Output prefix: " + prefix);
    log("Outputs: " + outputs);
    log("Subsample reads: " + std::to_string(subsample_reads));
    log("Subsample seeds: " + std::to_string(subsample_seeds));
    log("Reindex: " + std::to_string(reindex));
    log("k-mer length: " + std::to_string(k));
    log("s-mer length: " + std::to_string(s));


    bool build = true;
    ::capnp::MallocMessageBuilder outMessage;
    std::unique_ptr<::capnp::PackedFdMessageReader> inMessage;
    std::string default_index_path = guide + ".pmi";

    if (!index_path.empty() && !reindex) {
      log("Index path: " + index_path);
    } else {
      if (!reindex && fs::exists(default_index_path)) {
        log("Index loaded from: " + default_index_path);
        inMessage = readCapnp(default_index_path);
        build = false;
      } else if (reindex) {
        log("Reindexing.");
      } else {
        log("Index not found at: " + default_index_path + ", will build.");
      }
    }
    if (args["--stop-after"]) {
        std::string stop_after = args["--stop-after"].asString();
        log("Will stop after stage: " + stop_after);
    }

    if (build) {
      // build
      log("--- Run ---");
      log("Indexing...");

      // capnp index object
      Index::Builder index = outMessage.initRoot<Index>();
      
      index.setK(k);
      index.setS(s);
      index.setT(1);
      index.setOpen(false);
      index.setL(3);

      auto start = std::chrono::high_resolution_clock::now();
      pmi::build(T, index);
      auto end = std::chrono::high_resolution_clock::now();
      auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
      log("Build time: " + std::to_string(duration.count()) + " milliseconds");

      writeCapnp(outMessage, default_index_path);
    }




    sat::mutationMatrices mutMat;
    std::string mutmat_path = args["--mutmat"] ? args["--mutmat"].asString() : "";
    std::string default_mutmat_path = guide + ".mm";
    if (!mutmat_path.empty()) {
      log("Loading mutation matrices from: " + mutmat_path);
      std::ifstream mutmat_file(mutmat_path);
      sat::fillMutationMatricesFromFile(mutMat, mutmat_file);
    } else if (fs::exists(default_mutmat_path)) {
      log("Loading default mutation matrices from: " + default_mutmat_path);
      std::ifstream mutmat_file(default_mutmat_path);
      sat::fillMutationMatricesFromFile(mutMat, mutmat_file);
    } else {
      log("No mutation matrices found, building...");
      sat::fillMutationMatricesFromTree_test(mutMat, T, default_mutmat_path);
    }


    // Placement
    log("Reading...");
    inMessage = readCapnp(default_index_path);
    Index::Reader index_input = inMessage->getRoot<Index>();

    log("Placing...");
    auto start = std::chrono::high_resolution_clock::now();
    if (placement_per_read) {
      int maximumGap                = std::stoi(args["--maximum-gap"].asString());
      int minimumCount              = std::stoi(args["--minimum-count"].asString());
      int minimumScore              = std::stoi(args["--minimum-score"].asString());
      double errorRate              = std::stod(args["--error-rate"].asString());
      int redoReadThreshold         = std::stoi(args["--redo-read-threshold"].asString());
      bool recalculateScore         = args["--recalculate-score"] && args["--recalculate-score"].isBool() ? args["--recalculate-score"].asBool() : false;
      bool rescueDuplicates         = args["--rescue-duplicates"] && args["--rescue-duplicates"].isBool() ? args["--rescue-duplicates"].asBool() : false;
      int rescueDuplicatesThreshold = std::stoi(args["--rescue-duplicates-threshold"].asString());
      int filterRound               = std::stoi(args["--filter-round"].asString());
      int checkFrequency            = std::stoi(args["--check-frequency"].asString());
      int removeIteration           = std::stoi(args["--remove-iteration"].asString());
      double insigProb              = std::stod(args["--insig-prob"].asString());
      int roundsRemove              = std::stoi(args["--rounds-remove"].asString());
      double removeThreshold        = std::stod(args["--remove-threshold"].asString());
      bool leafNodesOnly            = args["--leaf-nodes-only"] && args["--leaf-nodes-only"].isBool() ? args["--leaf-nodes-only"].asBool() : false;

      pmi::place_per_read(
        T, index_input, reads1, reads2, maximumGap, minimumCount, minimumScore,
        errorRate, redoReadThreshold, recalculateScore, rescueDuplicates,
        rescueDuplicatesThreshold, filterRound, checkFrequency, removeIteration,
        insigProb, roundsRemove, removeThreshold, leafNodesOnly);
    } else {
      pmi::place(T, index_input, reads1, reads2, mutMat, refFileName, samFileName, bamFileName, mpileupFileName, vcfFileName);
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    log("Placement time: " + std::to_string(duration.count()) + " milliseconds");



    // Mapping
    //log("Mapping...");
    // Mapping logic here

    // Genotyping
    //log("Genotyping...");
    // Genotyping logic here

    




    // Assembly
    //log("Assembly...");
    // Assembly logic here

    log("panmap run completed.");

    // // Keep the application running
    // std::cout << "Press Ctrl+C to exit." << std::endl;
    // std::signal(SIGINT, [](int signum) {
    //     std::cout << "\nExiting panmap..." << std::endl;
    //     std::exit(signum);
    // });

    // // Infinite loop to keep the application running
    // while (true) {
    //     std::this_thread::sleep_for(std::chrono::seconds(1));
    // }

    return 0;
}
