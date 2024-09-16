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


using namespace pmi;
namespace po = boost::program_options;
namespace fs = boost::filesystem;

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

Seeding/alignment options:
  -k <k>                             Length of k-mer seeds. [default: 19]
  -s <s>                             Length of s-mers for seed (syncmer) selection. [default: 8]
  -R --subsample-reads <coverage>    Subsample reads for placement to approximately <coverage> depth. [default: 20.0]
  -S --subsample-seeds <proportion>  Use only <proportion> of total seeds for placement. [default: 1.0]
  -P --prior                         Compute and use a mutation spectrum prior for genotyping.
  -f --reindex                       Don't load index from disk, build it from scratch.

  -M --mutmat <path>                 Provide a mutation spectrum matrix instead of computing.
                                       one from the tree. Overrides --prior. [default: ]
  -I --identity-threshold <cutoff>   Identity cutoff for mutation spectrum, effective with --prior. [default: 0.80]
Other options:
  -c --cpus <num>            Number of CPUs to use. [default: 1]
  -x --stop-after <stage>        Stop after the specified stage. Accepted values:
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

    double subsample_reads = std::stod(args["--subsample-reads"].asString());
    double subsample_seeds = std::stod(args["--subsample-seeds"].asString());
    bool reindex = args["--reindex"] && args["--reindex"].isBool() ? args["--reindex"].asBool() : false;
    bool prior = args["--prior"] && args["--prior"].isBool() ? args["--prior"].asBool() : false;
  
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
    if (!index_path.empty() && !reindex) {
      log("Index path: " + index_path);
    } else {
      std::string default_index_path = guide + ".pmi";
      if (!reindex && fs::exists(default_index_path)) {
        log("Index loaded from: " + default_index_path);
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
    log("");
    log("--- Run ---");
    log("Indexing...");

    // capnp index object

    // ::capnp::MallocMessageBuilder message;
    // Index::Builder index = message.initRoot<Index>();
    // index.setK(k);
    // index.setS(s);
    // index.setL(3);
    // auto start = std::chrono::high_resolution_clock::now();
    // pmi::build(T, index);
    // auto end = std::chrono::high_resolution_clock::now();
    // auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    // log("Build time: " + std::to_string(duration.count()) + " milliseconds");

    // std::string tst = "atest.pmi"; 
    // writeCapnp(message, tst);

    // Placement
    log("Reading...");
    std::string tst = "rsv_4k.pmi";
    auto message = readCapnp(tst);
    Index::Reader index_input = message->getRoot<Index>();

    log("Placing...");
    pmi::place(T, index_input, reads1, reads2);


    // Mapping
    log("Mapping...");
    // Mapping logic here

    // Genotyping
    log("Genotyping...");
    // Genotyping logic here

    // Assembly
    log("Assembly...");
    // Assembly logic here

    log("panmap run completed.");
  return 0;
}