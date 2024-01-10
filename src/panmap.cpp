#include "kseq.h"
#include "PangenomeMAT.hpp"
#include "pmi.hpp"
#include "place.hpp"
#include <iostream>
#include <string>
#include <cstdio>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

using namespace pmi;
namespace po = boost::program_options;
namespace fs = boost::filesystem;

void loadReads(const std::string& reads1, const std::string& reads2) {
    // Implement your loadReads function here
    std::cout << "Loading reads..." << std::endl;
}

void promptAndPlace(Tree *T, const std::string &indexFile, const std::string &pmatFile, std::string &reads1File, std::string &reads2File, const bool prompt) {
    int32_t k = 15;
    int32_t s = 6;
    std::cin.clear();
    std::fflush(stdin);
    using namespace std;
    string userInput = "";
    // Confirm unless -f is specified
    if (prompt) {
        cout << "Place sample now? (y)es / (q)uit: ";
        getline(cin, userInput);
        if (!(userInput == "Y" || userInput == "y")) {
            exit(0);
        }
    }
    // Check if reads were supplied by user       
    if (reads1File == "" && reads2File == "" && !prompt) {
        cout << "Can't place sample because no reads were provided and -f was specified." << std::endl;
        exit(1);
    }
    // Prompt for fastq files
    if (reads1File == "") {
        while (!reads1File.size()) {
            cout << "Enter path to first paired-end reads FASTQ (or single-end): ";
            getline(cin, reads1File);
        }
        string reads2File = "";
        cout << "Enter path to second paired-end reads FASTQ (leave empty for single-end): ";
        getline(cin, reads2File);
    }
    cout << "\n" << reads1File << std::endl;
    if (reads2File.size()) {
        cout << ", " << reads2File << endl;
    }
    std::ifstream ifs(indexFile);
    place::placeIsolate(ifs, reads1File, reads2File, T, k, s);
}
void promptAndIndex(Tree *T, int32_t k, int32_t s, const bool prompt, const std::string &indexFile) {
    using namespace std;
    string userInput = "";
    char nl;
    if (prompt && (k == -1 || s == -1)) {
        cout << "Create index? Y(es) / q(uit) ➜ ";
        getline(cin, userInput);
        if (!(userInput == "Y" || userInput == "y")) {
            exit(0);
        }
        cout << "Enter integer value for k ➜ ";
        cin >> k;
        cout << "Enter integer value for s ➜ ";
        cin >> s;
        cout << "Building (" << k << ", " << s << ") index ..." << std::endl;
    } else if (prompt) {
        cout << "\nIndex now with " << "(k,s) = (" << k << "," << s << ")?\n(Y)es / (q)uit ➜ ";
        getline(cin, userInput);
        if (!(userInput == "Y" || userInput == "y" || userInput == "")) {
            exit(0);
        }
    } else if (k == -1 || s == -1) {
        cout << "Can't build index because -f was specified and no parameters provided with -p." << std::endl;
        exit(1);
    }
    seedIndex index;
    pmi::build(index, T, k, s);
    std::cout << "Writing to " << indexFile << "..." << std::endl;
    std::ofstream fout(indexFile);
    pmi::write(fout, T, index);
    std::cout << "Done." << std::endl;
}

int main(int argc, char *argv[]) {
    try {
        po::options_description desc("Options");
        desc.add_options()
            ("panmat", po::value<std::string>()->required(), "Path to tree.pmat file (required)")
            ("reads1", po::value<std::string>(), "Path to first (or single-end) fastq file (optional)")
            ("reads2", po::value<std::string>(), "Path to second paired-end fastq file (optional)")
            ("index,i", po::value<std::string>(), "Path to index file")
            ("params,p", po::value<std::vector<int32_t>>()->multitoken(), "Specify syncmer parameters k,s for indexing (e.g. for k=13, s=8 use -p 13,8)")
            ("f", "Proceed without prompting for confirmation. Applies to index construction if -p is specified or no index is found at ( /path/to/pmat.pmi ). Applies to sample placement if reads are provided.");

        po::positional_options_description p;
        p.add("panmat", 1);
        p.add("reads1", 1);
        p.add("reads2", 1);
        po::variables_map vm;
        po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
        po::notify(vm);
        std::string pmatFile = vm["panmat"].as<std::string>();

        std::cout << "  ╭───────────────────╮" << std::endl;
        std::cout << "  │   ┏━┳━● pan       │" << std::endl;
        std::cout << "  │   ┃ ┣━━━● map  ◠◡ │" << std::endl;
        std::cout << "  │ ━━┫ ┗━━━━━━━○ ⤶   │" << std::endl;
        std::cout << "  │   ┗━━━● v0.0      │" << std::endl;
        std::cout << "  ╰───────────────────╯" << std::endl;
/* tidbits
 * ⋌⋋
 *  ━̲━̲̲━̲̲━̲━̲̲━̲
 * 
 * 〰〰〰〰  
 * ╭──╮   ⎂⎂⎂⎂ a̲̅  𝐚̢̲̲̅𝗰̲̅𝐠̲̅𝘁̲̅>< 𝐚̢̲̲̅𝗰̲̅𝐠̲̅𝘁̲̅≽
 * 
 * 
 *   |̲̅̇ ╏͡ |̲̅̇|̲̅̇╳  
 */ 
//  * 
//  * 𝐚̢̲̲̅𝗰̲̅𝐠̲̅𝘁̲̅≽  𝚊𝚌𝚐𝚝 ᵍᴳᴬᴬªᵃᴳᵍᵀᵗ ⍜⍜⍜ ṯṱṰ 
//         ⎺⎺⎻⎼⎽ ⏋⏌⬭⬬⬣𛱱𛱰𛱲ꜙ︙ ov╥╦╨╩╧╤═┯┷〜〰﹋﹌﹏̳̿̿𐀸𐀳𐀭𐀙𐀝𐀋𐀇𐀅⫰⫯‾s℄ˈˌa̲̅c̲̅ ⎁
//  */


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
                std::cout << "Index found at " << defaultIndexPath << std::endl;
                std::cout << "Use it? (Y)es / (r)ebuild / (q)uit: ";
                getline(std::cin, inp);
                if (inp == "Y" || inp == "y" || inp == "") {
                    keep = true;
                } else if (inp == "Q" || inp == "q") {
                    std::exit(0);
                }
            } else {
                std::cout << "No index at \e[3;1m" << pmatFile << ".pmi\e[0m" << std::endl;
            }
            if (!keep) {
                promptAndIndex(T, k, s, prompt, defaultIndexPath);
            }
        }

        promptAndPlace(T, indexFile, pmatFile, reads1File, reads2File, prompt);

    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}




// struct greater
// {
//     template<class T>
//     bool operator()(T const &a, T const &b) const { return a > b; }
// };

// // Index construction
// void PangenomeMAT::discardSeeds(std::vector<kmer_t> &inSeeds, const tbb::concurrent_vector<std::pair<int32_t, int32_t>>& B, std::string &gappedSequence, tbb::concurrent_unordered_map<std::string, kmer_t> &to_insert, seedIndex &index, std::string nid, size_t k) {
//     // B should be sorted by start
//     std::mutex mtx;
//     auto cmp = [](kmer_t a, kmer_t b) { return a.idx > b.idx; };
//     std::set<kmer_t, decltype(cmp)> dels(cmp);

//     tbb::parallel_for((size_t)0, inSeeds.size(), [&](size_t i) {
//         kmer_t &s = inSeeds[i];
//         for (int32_t j = 0; j < B.size(); j++) {
//             const auto &b = B[j];
//             if (b.first > s.pos) {
//                 break;
//             } else if (s.gappedEnd <= b.second) {
//                 mtx.lock();
//                 auto it = to_insert.find(s.seq);
//                 if (it == to_insert.end()) {
//                     s.idx = i;
//                     dels.insert(s);
//                     mtx.unlock();
//                 } else {
//                     to_insert.unsafe_erase(it);
//                     mtx.unlock();
//                     s.pos = it->second.pos;
//                     s.gappedEnd = it->second.gappedEnd;
//                 }
//             }
//         }
//     });
    
//     std::stack<int32_t> rmDel;
//     for (const kmer_t &d : dels) {
//         rmDel.push(d.idx);
//     }

//     index.deletions[nid] = std::vector<kmer_t>(dels.begin(), dels.end());
    
//     removeIndices(inSeeds, rmDel);

// }

// bool compareTuples(const std::pair<int32_t, int32_t>& a, const std::pair<int32_t, int32_t>& b) {
//     return a.first < b.first;
// }

// void PangenomeMAT::applyMutations(Tree *T, sequence_t &sequence, Node *node, mutInfo &mutInfo, int32_t k) {
//     blockMutationInfo_t &blockInfo = mutInfo.blockInfo;
//     nucMutationInfo_t &nucInfo = mutInfo.nucInfo;
//     blockExists_t &blockExists = mutInfo.blockExists;
//     blockStrand_t &blockStrand = mutInfo.blockStrand;

//     for(const auto &mutation: node->blockMutation) {
//         int32_t primaryBlockId = mutation.primaryBlockId;
//         int32_t secondaryBlockId = mutation.secondaryBlockId;
//         bool type = mutation.blockMutInfo;
//         bool inversion = mutation.inversion;

//         int32_t blockStart = T->getGlobalCoordinate(primaryBlockId, secondaryBlockId, 0, -1);
//         //char c = getNucleotideFromCode(((T->blocks[primaryBlockId].consensusSeq[0]) >> (4*(7 - k))) & 15);
//         // nucInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, 0, 0, c, c, blockStart, 0));
//         if(type == 1) {
//             // insertion
//             bool oldStrand;
//             bool oldMut;
//             if(secondaryBlockId != -1) {
//                 oldStrand = blockStrand[primaryBlockId].second[secondaryBlockId];
//                 oldMut = blockExists[primaryBlockId].second[secondaryBlockId];
//                 blockExists[primaryBlockId].second[secondaryBlockId] = true;
//                 // if insertion of inverted block takes place, the strand is backwards
//                 blockStrand[primaryBlockId].second[secondaryBlockId] = !inversion;
//             } else {
//                 oldStrand = blockStrand[primaryBlockId].first;
//                 oldMut = blockExists[primaryBlockId].first;
//                 blockExists[primaryBlockId].first = true;

//                 // if insertion of inverted block takes place, the strand is backwards
//                 blockStrand[primaryBlockId].first = !inversion;
//             }
//             blockInfo.push_back( std::make_tuple(mutation.primaryBlockId, mutation.secondaryBlockId, oldMut, oldStrand, true, !inversion) );
//         } else {
//             bool oldMut;
//             bool oldStrand;
//             if(inversion) {
//                 // This means that this is not a deletion, but instead an inversion
//                 if(secondaryBlockId != -1) {
//                     oldStrand = blockStrand[primaryBlockId].second[secondaryBlockId];
//                     oldMut = blockExists[primaryBlockId].second[secondaryBlockId];
//                     blockStrand[primaryBlockId].second[secondaryBlockId] = !oldStrand;
//                 } else {
//                     oldStrand = blockStrand[primaryBlockId].first;
//                     oldMut = blockExists[primaryBlockId].first;
//                     blockStrand[primaryBlockId].first = !oldStrand;
//                 }
//                 if(oldMut != true) {
//                     std::cout << "Problem in PanMAT generation" << std::endl;
//                 }
//             blockInfo.push_back( std::make_tuple(mutation.primaryBlockId, mutation.secondaryBlockId, oldMut, oldStrand, oldMut, !oldStrand) );
//             } else {
//                 // Actually a deletion

//                 if(secondaryBlockId != -1) {
//                     oldStrand = blockStrand[primaryBlockId].second[secondaryBlockId];
//                     oldMut = blockExists[primaryBlockId].second[secondaryBlockId];
//                     blockExists[primaryBlockId].second[secondaryBlockId] = false;

//                     // resetting strand to true during deletion
//                     blockStrand[primaryBlockId].second[secondaryBlockId] = true;
//                 } else {
//                     oldStrand = blockStrand[primaryBlockId].first;
//                     oldMut = blockExists[primaryBlockId].first;
//                     blockExists[primaryBlockId].first = false;

//                     // resetting strand to true during deletion
//                     blockStrand[primaryBlockId].first = true;
//                 }
//             }
//             blockInfo.push_back( std::make_tuple(mutation.primaryBlockId, mutation.secondaryBlockId, oldMut, oldStrand, false, true) );
//         }
//     }

//     // Nuc mutations
//     for(size_t i = 0; i < node->nucMutation.size(); i++) {
//         int32_t primaryBlockId = node->nucMutation[i].primaryBlockId;
//         int32_t secondaryBlockId = node->nucMutation[i].secondaryBlockId;

//         int32_t nucPosition = node->nucMutation[i].nucPosition;
//         int32_t nucGapPosition = node->nucMutation[i].nucGapPosition;
//         uint32_t type = (node->nucMutation[i].mutInfo & 0x7);
//         char newVal = '-';
//         size_t globalCoord = T->getGlobalCoordinate(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition);
//         // switch(type) {
//         //     case PangenomeMAT::NucMutationType::NS:
//         //         std::cout << "NS " << " block:" << primaryBlockId << " ntPos:" << nucPosition << " ntGap:" << nucGapPosition << " globalPos:" << globalCoord << std::endl;
//         //         break;
//         //     case PangenomeMAT::NucMutationType::NI:
//         //         std::cout << "NI " << " block:" << primaryBlockId << " ntPos:" << nucPosition << " ntGap:" << nucGapPosition << " globalPos:" << globalCoord << std::endl;
//         //         break;
//         //     case PangenomeMAT::NucMutationType::ND:
//         //         std::cout << "ND " << " block:" << primaryBlockId << " ntPos:" << nucPosition << " ntGap:" << nucGapPosition << " globalPos:" << globalCoord << std::endl;
//         //         break;
//         //     case PangenomeMAT::NucMutationType::NSNPS:
//         //         std::cout << "NSNPS " << " block:" << primaryBlockId << " ntPos:" << nucPosition << " ntGap:" << nucGapPosition << " globalPos:" << globalCoord << std::endl;
//         //         break;
//         //     case PangenomeMAT::NucMutationType::NSNPI:
//         //         std::cout << "NSNPI " << " block:" << primaryBlockId << " ntPos:" << nucPosition << " ntGap:" << nucGapPosition << " globalPos:" << globalCoord << std::endl;
//         //         break;
//         //     case PangenomeMAT::NucMutationType::NSNPD:
//         //         std::cout << "NSNPD " << " block:" << primaryBlockId << " ntPos:" << nucPosition << " ntGap:" << nucGapPosition << " globalPos:" << globalCoord << std::endl;
//         //         break;
//         // }
//         // std::cout << "" << mutInfo.nucInfo.size() << " " << node->nucMutation.size() << "\n";
        

//         if(type < 3) {
//             // Either S, I or D

//             int len = ((node->nucMutation[i].mutInfo) >> 4);

//             if(type == PangenomeMAT::NucMutationType::NS) {
//                 // Substitution
//                 if(secondaryBlockId != -1) {
//                     if(nucGapPosition != -1) {
//                         for(int j = 0; j < len; j++) {
//                             char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition+j];
//                             newVal = PangenomeMAT::getNucleotideFromCode(((node->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
//                             sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition+j] = newVal;
//                             nucInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition+j, oldVal, newVal, globalCoord, len));
//                         }
//                     } else {
//                         for(int j = 0; j < len; j++) {
//                             char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition + j].first;
//                             newVal = PangenomeMAT::getNucleotideFromCode(((node->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
//                             sequence[primaryBlockId].second[secondaryBlockId][nucPosition + j].first = newVal;
//                             nucInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition + j, nucGapPosition, oldVal, newVal, globalCoord, len));
//                         }

//                     }
//                 } else {
//                     if(nucGapPosition != -1) {
//                         for(int j = 0; j < len; j++) {
//                             char oldVal = sequence[primaryBlockId].first[nucPosition].second[nucGapPosition+j];
//                             newVal = PangenomeMAT::getNucleotideFromCode(((node->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
//                             sequence[primaryBlockId].first[nucPosition].second[nucGapPosition+j] = newVal;
//                             nucInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition+j, oldVal, newVal, globalCoord, len));   
//                         }
//                     } else {
//                         for(int j = 0; j < len; j++) {
//                             char oldVal = sequence[primaryBlockId].first[nucPosition+j].first;
//                             newVal = PangenomeMAT::getNucleotideFromCode(((node->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
//                             sequence[primaryBlockId].first[nucPosition+j].first = newVal;
//                             nucInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition + j, nucGapPosition, oldVal, newVal, globalCoord, len));   
//                         }
//                     }
//                 }
//             }
//             else if(type == PangenomeMAT::NucMutationType::NI) {
//                 // Insertion
                
//                 if(secondaryBlockId != -1) {
//                     if(nucGapPosition != -1) {
//                         for(int j = 0; j < len; j++) {
//                             char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition + j];
//                             newVal = PangenomeMAT::getNucleotideFromCode(((node->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
//                             sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition+j] = newVal;
//                             nucInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition+j, oldVal, newVal, globalCoord, len));
//                         }
//                     } else {
//                         for(int j = 0; j < len; j++) {
//                             char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition + j].first;
//                             newVal = PangenomeMAT::getNucleotideFromCode(((node->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
//                             sequence[primaryBlockId].second[secondaryBlockId][nucPosition + j].first = newVal;
//                             nucInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition + j, nucGapPosition, oldVal, newVal, globalCoord, len));
//                         }

//                     }
//                 } else {
//                     if(nucGapPosition != -1) {
//                         for(int j = 0; j < len; j++) {
//                             char oldVal = sequence[primaryBlockId].first[nucPosition].second[nucGapPosition+j];
//                             newVal = PangenomeMAT::getNucleotideFromCode(((node->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
//                             sequence[primaryBlockId].first[nucPosition].second[nucGapPosition+j] = newVal;
//                             nucInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition+j, oldVal, newVal, globalCoord, len));
//                         }
//                     } else {
//                         for(int j = 0; j < len; j++) {
//                             char oldVal = sequence[primaryBlockId].first[nucPosition+j].first;
//                             const int nucCode = ((node->nucMutation[i].nucs) >> (4*(5-j))) & 0xF;
//                             //std::cout << "nucCode: " << nucCode << std::endl;
//                             newVal = PangenomeMAT::getNucleotideFromCode(((node->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
//                             sequence[primaryBlockId].first[nucPosition+j].first = newVal;
//                             nucInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition + j, nucGapPosition, oldVal, newVal, globalCoord, len));   
//                         }
//                     }
//                 }
//             }
//             else if(type == PangenomeMAT::NucMutationType::ND) {
//                 // Deletion
                
//                 if(secondaryBlockId != -1) {
//                     if(nucGapPosition != -1) {
//                         for(int j = 0; j < len; j++) {
//                             char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition+j];
//                             sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition+j] = '-';
//                             nucInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition+j, oldVal, '-', globalCoord, len));
//                         }
//                     } else {
//                         for(int j = 0; j < len; j++) {
//                             char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition + j].first;
//                             sequence[primaryBlockId].second[secondaryBlockId][nucPosition + j].first = '-';
//                             nucInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition + j, nucGapPosition, oldVal, '-', globalCoord, len));
//                         }

//                     }
//                 } else {
//                     if(nucGapPosition != -1) {
//                         for(int j = 0; j < len; j++) {
//                             char oldVal = sequence[primaryBlockId].first[nucPosition].second[nucGapPosition+j];
//                             sequence[primaryBlockId].first[nucPosition].second[nucGapPosition+j] = '-';
//                             nucInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition+j, oldVal, '-', globalCoord, len));
//                         }
//                     } else {
//                         for(int j = 0; j < len; j++) {
//                             char oldVal = sequence[primaryBlockId].first[nucPosition+j].first;
//                             sequence[primaryBlockId].first[nucPosition+j].first = '-';
//                             nucInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition + j, nucGapPosition, oldVal, '-', globalCoord, 0));
//                         }
//                     }
//                 }
//             }
//         } 
//         else {
//             int len = 0;
//             if(type == PangenomeMAT::NucMutationType::NSNPS) {
//                 // SNP Substitution
                
//                 newVal = PangenomeMAT::getNucleotideFromCode(((node->nucMutation[i].nucs) >> 20) & 0xF);
//                 if(secondaryBlockId != -1) {
//                     if(nucGapPosition != -1) {
//                         char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition];
//                         sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition] = newVal;
//                         nucInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, newVal, globalCoord, len));
//                     } else {
//                         char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition].first;
//                         sequence[primaryBlockId].second[secondaryBlockId][nucPosition].first = newVal;
//                         nucInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, newVal, globalCoord, len));
//                     }
//                 } else {
//                     if(nucGapPosition != -1) {
//                         char oldVal = sequence[primaryBlockId].first[nucPosition].second[nucGapPosition];
//                         sequence[primaryBlockId].first[nucPosition].second[nucGapPosition] = newVal;
//                         nucInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, newVal, globalCoord, len));
//                     } else {
//                         char oldVal = sequence[primaryBlockId].first[nucPosition].first;
//                         sequence[primaryBlockId].first[nucPosition].first = newVal;
//                         nucInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, newVal, globalCoord, len));   
//                     }
//                 }
//             }
//             else if(type == PangenomeMAT::NucMutationType::NSNPI) {
//                 // SNP Insertion
//                 len = 1;
//                 newVal = PangenomeMAT::getNucleotideFromCode(((node->nucMutation[i].nucs) >> 20) & 0xF);
//                 if(secondaryBlockId != -1) {
//                     if(nucGapPosition != -1) {
//                         char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition];
//                         sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition] = newVal;
//                         nucInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, newVal, globalCoord, len));
//                     } else {
//                         char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition].first;
//                         sequence[primaryBlockId].second[secondaryBlockId][nucPosition].first = newVal;
//                         nucInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, newVal, globalCoord, len));
//                     }
//                 } else {
//                     if(nucGapPosition != -1) {
//                         char oldVal = sequence[primaryBlockId].first[nucPosition].second[nucGapPosition];
//                         sequence[primaryBlockId].first[nucPosition].second[nucGapPosition] = newVal;
//                         nucInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, newVal, globalCoord, len));   
//                     } else {
//                         char oldVal = sequence[primaryBlockId].first[nucPosition].first;
//                         sequence[primaryBlockId].first[nucPosition].first = newVal;
//                         nucInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, newVal, globalCoord, len));   
//                     }
//                 }
//             }
//             else if(type == PangenomeMAT::NucMutationType::NSNPD) {
//                 // SNP Deletion
//                 if(secondaryBlockId != -1) {
//                     if(nucGapPosition != -1) {
//                         char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition];
//                         sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition] = '-';
//                         nucInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, '-', globalCoord, len));
//                     } else {
//                         char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition].first;
//                         sequence[primaryBlockId].second[secondaryBlockId][nucPosition].first = '-';
//                         nucInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, '-', globalCoord, len));
//                     }
//                 } else {
//                     if(nucGapPosition != -1) {
//                         char oldVal = sequence[primaryBlockId].first[nucPosition].second[nucGapPosition];
//                         sequence[primaryBlockId].first[nucPosition].second[nucGapPosition] = '-';
//                         nucInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, '-', globalCoord, len));
//                     } else {
//                         char oldVal = sequence[primaryBlockId].first[nucPosition].first;
//                         sequence[primaryBlockId].first[nucPosition].first = '-';
//                         nucInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, '-', globalCoord, len));
//                     }
//                 }
//             }
//         }
//     }
// }

// void PangenomeMAT::revertMutations(Tree *T, sequence_t &sequence, Node *node, mutInfo &mutInfo) {
//     blockMutationInfo_t &blockInfo = mutInfo.blockInfo;
//     nucMutationInfo_t &nucInfo = mutInfo.nucInfo;
//     blockExists_t &blockExists = mutInfo.blockExists;
//     blockStrand_t &blockStrand = mutInfo.blockStrand;

//     for(auto it = blockInfo.rbegin(); it != blockInfo.rend(); it++) {
//         auto mutation = *it;
//         if(std::get<1>(mutation) != -1) {
//             blockExists[std::get<0>(mutation)].second[std::get<1>(mutation)] = std::get<2>(mutation);
//             blockStrand[std::get<0>(mutation)].second[std::get<1>(mutation)] = std::get<3>(mutation);
//         } else {
//             blockExists[std::get<0>(mutation)].first = std::get<2>(mutation);
//             blockStrand[std::get<0>(mutation)].first = std::get<3>(mutation);
//         }
//     }

//     // Undo nuc mutations when current node and its subtree have been processed
//     for(auto it = nucInfo.rbegin(); it != nucInfo.rend(); it++) {
//         auto mutation = *it;
//         if(std::get<1>(mutation) != -1) {
//             if(std::get<3>(mutation) != -1) {
//                 sequence[std::get<0>(mutation)].second[std::get<1>(mutation)][std::get<2>(mutation)].second[std::get<3>(mutation)] = std::get<4>(mutation);
//             } else {
//                 sequence[std::get<0>(mutation)].second[std::get<1>(mutation)][std::get<2>(mutation)].first = std::get<4>(mutation);
//             }
//         } else {
//             if(std::get<3>(mutation) != -1) {
//                 sequence[std::get<0>(mutation)].first[std::get<2>(mutation)].second[std::get<3>(mutation)] = std::get<4>(mutation);
//             } else {
//                 sequence[std::get<0>(mutation)].first[std::get<2>(mutation)].first = std::get<4>(mutation);
//             }
//         }
//     }
// }


// std::vector<PangenomeMAT::Node *> PangenomeMAT::getParallelStartNodes(Node *node, int32_t numThreads) {

//     int32_t perLevelCtr = 0;

//     std::queue<Node *> bfsQueue;
//     std::vector<Node *> levelNodes;
//     std::vector<Node *> prevLevelNodes;


//     size_t prevLev = 0;
    
//     bfsQueue.push(node);
//     levelNodes.push_back(node);

//     while(!bfsQueue.empty()) {
//         Node* current = bfsQueue.front();
//         bfsQueue.pop();
//         levelNodes.push_back(current);

//         if(current->level != prevLev) {
//             if (perLevelCtr > numThreads) {
//                 return prevLevelNodes;
//             }
//             perLevelCtr = 0;
//             prevLev = current->level;
//             prevLevelNodes = levelNodes;
//             levelNodes.clear();
//         }
//         perLevelCtr++;

//         for(auto child: current->children) {
//             bfsQueue.push(child);
//         }
//     }
//     return levelNodes;
// }


// void PangenomeMAT::applyBlockMutations(Tree *T, Node *node, int32_t k, blockMutationInfo_t &blockInfo, blockExists_t &blockExists, blockStrand_t &blockStrand) {

//     for(const auto &mutation: node->blockMutation) {
//         int32_t primaryBlockId = mutation.primaryBlockId;
//         int32_t secondaryBlockId = mutation.secondaryBlockId;
//         bool type = mutation.blockMutInfo;
//         bool inversion = mutation.inversion;

//         int32_t blockStart = T->getGlobalCoordinate(primaryBlockId, secondaryBlockId, 0, -1);
//         //char c = getNucleotideFromCode(((T->blocks[primaryBlockId].consensusSeq[0]) >> (4*(7 - k))) & 15);
//         // nucInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, 0, 0, c, c, blockStart, 0));
//         if(type == 1) {
//             // insertion
//             bool oldStrand;
//             bool oldMut;
//             if(secondaryBlockId != -1) {
//                 oldStrand = blockStrand[primaryBlockId].second[secondaryBlockId];
//                 oldMut = blockExists[primaryBlockId].second[secondaryBlockId];
//                 blockExists[primaryBlockId].second[secondaryBlockId] = true;
//                 // if insertion of inverted block takes place, the strand is backwards
//                 blockStrand[primaryBlockId].second[secondaryBlockId] = !inversion;
//             } else {
//                 oldStrand = blockStrand[primaryBlockId].first;
//                 oldMut = blockExists[primaryBlockId].first;
//                 blockExists[primaryBlockId].first = true;

//                 // if insertion of inverted block takes place, the strand is backwards
//                 blockStrand[primaryBlockId].first = !inversion;
//             }
//             blockInfo.push_back( std::make_tuple(mutation.primaryBlockId, mutation.secondaryBlockId, oldMut, oldStrand, true, !inversion) );
//         } else {
//             bool oldMut;
//             bool oldStrand;
//             if(inversion) {
//                 // This means that this is not a deletion, but instead an inversion
//                 if(secondaryBlockId != -1) {
//                     oldStrand = blockStrand[primaryBlockId].second[secondaryBlockId];
//                     oldMut = blockExists[primaryBlockId].second[secondaryBlockId];
//                     blockStrand[primaryBlockId].second[secondaryBlockId] = !oldStrand;
//                 } else {
//                     oldStrand = blockStrand[primaryBlockId].first;
//                     oldMut = blockExists[primaryBlockId].first;
//                     blockStrand[primaryBlockId].first = !oldStrand;
//                 }
//                 if(oldMut != true) {
//                     std::cout << "Problem in PanMAT generation" << std::endl;
//                 }
//             blockInfo.push_back( std::make_tuple(mutation.primaryBlockId, mutation.secondaryBlockId, oldMut, oldStrand, oldMut, !oldStrand) );
//             } else {
//                 // Actually a deletion

//                 if(secondaryBlockId != -1) {
//                     oldStrand = blockStrand[primaryBlockId].second[secondaryBlockId];
//                     oldMut = blockExists[primaryBlockId].second[secondaryBlockId];
//                     blockExists[primaryBlockId].second[secondaryBlockId] = false;

//                     // resetting strand to true during deletion
//                     blockStrand[primaryBlockId].second[secondaryBlockId] = true;
//                 } else {
//                     oldStrand = blockStrand[primaryBlockId].first;
//                     oldMut = blockExists[primaryBlockId].first;
//                     blockExists[primaryBlockId].first = false;

//                     // resetting strand to true during deletion
//                     blockStrand[primaryBlockId].first = true;
//                 }
//             }
//             blockInfo.push_back( std::make_tuple(mutation.primaryBlockId, mutation.secondaryBlockId, oldMut, oldStrand, false, true) );
//         }
//     }
// }

// void PangenomeMAT::revertBlockMutations(Tree *T, Node *node, int32_t k, blockMutationInfo_t &blockInfo, blockExists_t &blockExists, blockStrand_t &blockStrand) {
//     for(auto it = blockInfo.rbegin(); it != blockInfo.rend(); it++) {
//         auto mutation = *it;
//         blockExists[std::get<0>(mutation)].first = std::get<2>(mutation);
//         blockStrand[std::get<0>(mutation)].first = std::get<3>(mutation);
//     }
// }

