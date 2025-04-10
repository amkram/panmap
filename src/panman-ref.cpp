
void panmanUtils::Tree::printMutationsNew(std::ostream& fout) {

    // Get reference sequence
    sequence_t rootSequence;
    blockExists_t rootBlockExists;
    blockStrand_t rootBlockStrand;
    getSequenceFromReference(rootSequence, rootBlockExists, rootBlockStrand, root->identifier);

    // sequence_t st;
    // blockExists_t bt;
    // blockStrand_t bst;
    // getSequenceFromReference(st, bt, bst, "USA/CA-CDC-QDX21497008/2021|MW666944.1|2021-01-27");

    // std::cout << st[55].first[132].second[40] << std::endl;


    tbb::concurrent_map< std::tuple< int, int, int >, size_t > panMATCoordinateToGlobal;
    tbb::concurrent_map< std::tuple< int, int, int >, char > rootCurrentCharacter;

    tbb::concurrent_unordered_map< std::string,
        std::vector< std::tuple< char, size_t, char, char, bool > > > nodeMutations;

    tbb::concurrent_unordered_map< size_t, bool > rootPresentBlocks;

    tbb::concurrent_map< std::tuple< int, int, int >, bool > isGapCoordinate;

    // convert PanMAT coordinate to global reference coordinate
    size_t rootCtr = 0;
    for(size_t i = 0; i < rootSequence.size(); i++) {
        if(rootBlockExists[i].first) {
            rootPresentBlocks[i] = true;
            if(rootBlockStrand[i].first) {
                for(size_t j = 0; j < rootSequence[i].first.size(); j++) {
                    for(size_t k = 0; k < rootSequence[i].first[j].second.size(); k++) {
                        panMATCoordinateToGlobal[std::make_tuple(i,j,k)] = rootCtr;
                        if(rootSequence[i].first[j].second[k] != '-' && rootSequence[i].first[j].second[k] != 'x') {
                            rootCurrentCharacter[std::make_tuple(i,j,k)] = rootSequence[i].first[j].second[k];
                            isGapCoordinate[std::make_tuple(i,j,k)] = false;
                            rootCtr++;
                        } else {
                            rootCurrentCharacter[std::make_tuple(i,j,k)] = '-';
                            isGapCoordinate[std::make_tuple(i,j,k)] = true;
                        }
                    }
                    panMATCoordinateToGlobal[std::make_tuple(i,j,-1)] = rootCtr;
                    if(rootSequence[i].first[j].first != '-' && rootSequence[i].first[j].first != 'x') {
                        rootCurrentCharacter[std::make_tuple(i,j,-1)] = rootSequence[i].first[j].first;
                        isGapCoordinate[std::make_tuple(i,j,-1)] = false;
                        rootCtr++;
                    } else {
                        rootCurrentCharacter[std::make_tuple(i,j,-1)] = '-';
                        isGapCoordinate[std::make_tuple(i,j,-1)] = true;
                    }
                }
            } else {
                for(size_t j = rootSequence[i].first.size() - 1; j + 1 > 0; j--) {
                    panMATCoordinateToGlobal[std::make_tuple(i,j,-1)] = rootCtr;
                    if(rootSequence[i].first[j].first != '-' && rootSequence[i].first[j].first != 'x') {
                        rootCurrentCharacter[std::make_tuple(i,j,-1)] = rootSequence[i].first[j].first;
                        isGapCoordinate[std::make_tuple(i,j,-1)] = false;
                        rootCtr++;
                    } else {
                        rootCurrentCharacter[std::make_tuple(i,j,-1)] = '-';
                        isGapCoordinate[std::make_tuple(i,j,-1)] = true;
                    }
                    for(size_t k = rootSequence[i].first[j].second.size() - 1; k + 1 > 0; k--) {
                        panMATCoordinateToGlobal[std::make_tuple(i,j,k)] = rootCtr;
                        if(rootSequence[i].first[j].second[k] != '-' && rootSequence[i].first[j].second[k] != 'x') {
                            rootCurrentCharacter[std::make_tuple(i,j,k)] = rootSequence[i].first[j].second[k];
                            // if(rootCtr == 240) {
                            //     std::cout << rootSequence[i].first[j].first << " " << i << " " << j << " " << k << std::endl;
                            // }
                            isGapCoordinate[std::make_tuple(i,j,k)] = false;
                            rootCtr++;
                        } else {
                            rootCurrentCharacter[std::make_tuple(i,j,k)] = '-';
                            isGapCoordinate[std::make_tuple(i,j,k)] = true;
                        }
                    }
                }
            }
        } else {
            rootPresentBlocks[i] = false;
            if(rootBlockStrand[i].first) {
                for(size_t j = 0; j < rootSequence[i].first.size(); j++) {
                    for(size_t k = 0; k < rootSequence[i].first[j].second.size(); k++) {
                        panMATCoordinateToGlobal[std::make_tuple(i,j,k)] = rootCtr;
                        rootCurrentCharacter[std::make_tuple(i,j,k)] = '-';
                        isGapCoordinate[std::make_tuple(i,j,k)] = false;
                    }
                    panMATCoordinateToGlobal[std::make_tuple(i,j,-1)] = rootCtr;
                    rootCurrentCharacter[std::make_tuple(i,j,-1)] = '-';
                    isGapCoordinate[std::make_tuple(i,j,-1)] = false;
                }
            } else {
                for(size_t j = rootSequence[i].first.size() - 1; j + 1 > 0; j--) {
                    panMATCoordinateToGlobal[std::make_tuple(i,j,-1)] = rootCtr;
                    rootCurrentCharacter[std::make_tuple(i,j,-1)] = '-';
                    isGapCoordinate[std::make_tuple(i,j,-1)] = false;
                    for(size_t k = rootSequence[i].first[j].second.size() - 1; k + 1 > 0; k--) {
                        panMATCoordinateToGlobal[std::make_tuple(i,j,k)] = rootCtr;
                        rootCurrentCharacter[std::make_tuple(i,j,k)] = '-';
                        isGapCoordinate[std::make_tuple(i,j,k)] = false;
                    }
                }
            }
        }
    }

    tbb::concurrent_map< std::tuple< std::string, int, int, int >, char > seqChar;

    tbb::parallel_for_each(allNodes, [&](auto u) {
        sequence_t st;
        blockExists_t bt;
        blockStrand_t bst;
        getSequenceFromReference(st, bt, bst, u.first);

        for(size_t i = 0; i < st.size(); i++) {
            if(bst[i].first) {
                for(size_t j = 0; j < st[i].first.size(); j++) {
                    for(size_t k = 0; k < st[i].first[j].second.size(); k++) {
                        if(st[i].first[j].second[k] != '-' && st[i].first[j].second[k] != 'x') {
                            seqChar[std::make_tuple(u.first,i,j,k)] = st[i].first[j].second[k];
                        } else {
                            seqChar[std::make_tuple(u.first,i,j,k)] = '-';
                        }
                    }
                    if(st[i].first[j].first != '-' && st[i].first[j].first != 'x' && bt[i].first) {
                        seqChar[std::make_tuple(u.first,i,j,-1)] = st[i].first[j].first;
                    } else {
                        seqChar[std::make_tuple(u.first,i,j,-1)] = '-';
                    }
                }
            } else {
                for(size_t j = st[i].first.size() - 1; j + 1 > 0; j--) {
                    if(st[i].first[j].first != '-' && st[i].first[j].first != 'x') {
                        seqChar[std::make_tuple(u.first,i,j,-1)] = rootSequence[i].first[j].first;
                    } else {
                        seqChar[std::make_tuple(u.first,i,j,-1)] = '-';
                    }
                    for(size_t k = st[i].first[j].second.size() - 1; k + 1 > 0; k--) {
                        if(st[i].first[j].second[k] != '-' && st[i].first[j].second[k] != 'x') {
                            seqChar[std::make_tuple(u.first,i,j,k)] = st[i].first[j].second[k];
                        } else {
                            seqChar[std::make_tuple(u.first,i,j,k)] = '-';
                        }
                    }
                }
            }
        }
    });

    nodeMutations[root->identifier];

    // Compute mutations for each of the other sequences
    tbb::parallel_for_each(allNodes, [&](auto u) {
        if(u.first == root->identifier) {
            return;
        }

        tbb::concurrent_map< std::tuple< int, int, int >, char > currentCharacter = rootCurrentCharacter;
        tbb::concurrent_unordered_map< size_t, bool > presentBlocks = rootPresentBlocks;

        nodeMutations[u.first];

        Node* it = u.second;
        std::vector< panmanUtils::Node* > path;

        while(it != root) {
            path.push_back(it);
            it = it->parent;
        }

        std::vector< std::pair< size_t, std::tuple< char, size_t, char, char, bool > > > currentNodeMutations;

        for(auto node = path.rbegin(); node != path.rend(); node++) {
            for(auto mutation: (*node)->blockMutation) {
                int32_t primaryBlockId = mutation.primaryBlockId;
                if(rootPresentBlocks.find(primaryBlockId) == rootPresentBlocks.end()) {
                    continue;
                }

                bool type = mutation.blockMutInfo;
                bool inversion = mutation.inversion;
                if(type == 1) {
                    // insertion
                    presentBlocks[primaryBlockId] = true;
                    if(inversion) {
                        std::cout << "INVERTED BLOCK FOUND" << std::endl;
                    }
                } else {
                    if(inversion) {
                        // This means that this is not a deletion, but instead an inversion
                        std::cout << "INVERSION FOUND" << std::endl;
                    } else {
                        // Actually a deletion
                        presentBlocks[primaryBlockId] = false;
                    }
                }
            }
        }

        for(auto node = path.rend()-1; node != path.rend(); node++) {
            for(size_t i = 0; i < (*node)->nucMutation.size(); i++) {
                int32_t primaryBlockId = (*node)->nucMutation[i].primaryBlockId;
                if(rootPresentBlocks.find(primaryBlockId) == rootPresentBlocks.end()) {
                    continue;
                }

                int32_t nucPosition = (*node)->nucMutation[i].nucPosition;
                int32_t nucGapPosition = (*node)->nucMutation[i].nucGapPosition;
                uint32_t type = ((*node)->nucMutation[i].mutInfo & 0x7);
                char newVal = '-';

                if(type < 3) {
                    int len = (((*node)->nucMutation[i].mutInfo) >> 4);

                    if(type == panmanUtils::NucMutationType::NS) {
                        if(nucGapPosition != -1) {
                            for(int j = 0; j < len; j++) {
                                newVal = panmanUtils::getNucleotideFromCode((((*node)->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                                if(presentBlocks[primaryBlockId]) {
                                    // char oldVal = currentCharacter[std::make_tuple(primaryBlockId, nucPosition, nucGapPosition+j)];
                                    char oldVal = seqChar[std::make_tuple((*node)->parent->identifier,primaryBlockId, nucPosition, nucGapPosition+j)];
                                    if(oldVal == '-' || oldVal == 'x') {
                                        continue;
                                    }
                                    // if(node == path.rend()-1)
                                    currentNodeMutations.push_back(std::make_pair(primaryBlockId, std::make_tuple('S', panMATCoordinateToGlobal[std::make_tuple(primaryBlockId, nucPosition, nucGapPosition+j)], oldVal, newVal, isGapCoordinate[std::make_tuple(primaryBlockId, nucPosition, nucGapPosition+j)])));
                                    // nodeMutations[u.first].push_back(std::make_tuple('S', panMATCoordinateToGlobal[std::make_tuple(primaryBlockId, nucPosition, nucGapPosition+j)], oldVal, newVal));
                                }
                                currentCharacter[std::make_tuple(primaryBlockId, nucPosition, nucGapPosition+j)] = newVal;
                            }
                        } else {
                            for(int j = 0; j < len; j++) {
                                newVal = panmanUtils::getNucleotideFromCode((((*node)->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                                if(presentBlocks[primaryBlockId]) {
                                    char oldVal = seqChar[std::make_tuple((*node)->parent->identifier,primaryBlockId, nucPosition+j, -1)];
                                    // char oldVal = currentCharacter[std::make_tuple(primaryBlockId, nucPosition + j, -1)];
                                    // if(u.first == "Denmark/DCGC-504971/2022|OX187739.1|2022-04-28" && panMATCoordinateToGlobal[std::make_tuple(primaryBlockId, nucPosition + j, -1)] == 185) {
                                    //     std::cout << oldVal << " " << newVal << " " << primaryBlockId << " " << nucPosition+j << " " << -1 << std::endl;
                                    // }
                                    if(oldVal == '-' || oldVal == 'x') {
                                        // std::cout << "NOT ACTUALLY A SUBSTITUTION" << std::endl;
                                        continue;
                                    }
                                    // if(oldVal == newVal) {
                                    //     std::cout << primaryBlockId << " " << nucPosition << " " << nucGapPosition+j << " " << oldVal << " " << newVal << std::endl;
                                    // }
                                    // if(node == path.rend()-1)
                                    currentNodeMutations.push_back(std::make_pair(primaryBlockId, std::make_tuple('S', panMATCoordinateToGlobal[std::make_tuple(primaryBlockId, nucPosition + j, -1)], oldVal, newVal, isGapCoordinate[std::make_tuple(primaryBlockId, nucPosition+j, -1)])));
                                    // nodeMutations[u.first].push_back(std::make_tuple('S', panMATCoordinateToGlobal[std::make_tuple(primaryBlockId, nucPosition + j, -1)], oldVal, newVal));
                                }
                                currentCharacter[std::make_tuple(primaryBlockId, nucPosition + j, -1)] = newVal;
                            }
                        }
                    } else if(type == panmanUtils::NucMutationType::NI) {
                        if(nucGapPosition != -1) {
                            for(int j = 0; j < len; j++) {
                                newVal = panmanUtils::getNucleotideFromCode((((*node)->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                                // if(node == path.rend()-1)
                                currentNodeMutations.push_back(std::make_pair(primaryBlockId, std::make_tuple('I', panMATCoordinateToGlobal[std::make_tuple(primaryBlockId, nucPosition, nucGapPosition + j)], '-', newVal, isGapCoordinate[std::make_tuple(primaryBlockId, nucPosition, nucGapPosition + j)])));
                            }
                        } else {
                            for(int j = 0; j < len; j++) {
                                newVal = panmanUtils::getNucleotideFromCode((((*node)->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                                // if(node == path.rend()-1)
                                currentNodeMutations.push_back(std::make_pair(primaryBlockId, std::make_tuple('I', panMATCoordinateToGlobal[std::make_tuple(primaryBlockId, nucPosition + j, -1)], '-', newVal, isGapCoordinate[std::make_tuple(primaryBlockId, nucPosition + j, -1)])));

                            }
                        }
                    } else if(type == panmanUtils::NucMutationType::ND) {
                        if(nucGapPosition != -1) {
                            for(int j = 0; j < len; j++) {
                                char oldVal = seqChar[std::make_tuple((*node)->parent->identifier,primaryBlockId, nucPosition, nucGapPosition + j)];
                                // if(node == path.rend()-1){
                                currentNodeMutations.push_back(std::make_pair(primaryBlockId, std::make_tuple('D', panMATCoordinateToGlobal[std::make_tuple(primaryBlockId, nucPosition, nucGapPosition + j)], oldVal, '-', isGapCoordinate[std::make_tuple(primaryBlockId, nucPosition, nucGapPosition + j)])));
                                // }
                            }
                        } else {
                            for(int j = 0; j < len; j++) {
                                char oldVal = seqChar[std::make_tuple((*node)->parent->identifier,primaryBlockId, nucPosition + j, nucGapPosition)];
                                // if(node == path.rend()-1){
                                currentNodeMutations.push_back(std::make_pair(primaryBlockId, std::make_tuple('D', panMATCoordinateToGlobal[std::make_tuple(primaryBlockId, nucPosition + j, -1)], oldVal, '-', isGapCoordinate[std::make_tuple(primaryBlockId, nucPosition + j, -1)])));
                                // }
                            }
                        }
                    }
                } else if(type == panmanUtils::NucMutationType::NSNPS) {
                    newVal = panmanUtils::getNucleotideFromCode((((*node)->nucMutation[i].nucs) >> 20) & 0xF);
                    if(nucGapPosition != -1) {
                        if(presentBlocks[primaryBlockId]) {
                            char oldVal = seqChar[std::make_tuple((*node)->parent->identifier,primaryBlockId, nucPosition, nucGapPosition)];
                            // char oldVal = currentCharacter[std::make_tuple(primaryBlockId, nucPosition, nucGapPosition)];
                            if(oldVal == '-' || oldVal == 'x') {
                                std::cout << "NOT ACTUALLY A SUBSTITUTION" << std::endl;
                            }
                            // if(node == path.rend()-1)
                            currentNodeMutations.push_back(std::make_pair(primaryBlockId, std::make_tuple('S', panMATCoordinateToGlobal[std::make_tuple(primaryBlockId, nucPosition, nucGapPosition)], oldVal, newVal, isGapCoordinate[std::make_tuple(primaryBlockId, nucPosition, nucGapPosition)])));
                            // nodeMutations[u.first].push_back(std::make_tuple('S', panMATCoordinateToGlobal[std::make_tuple(primaryBlockId, nucPosition, nucGapPosition)], oldVal, newVal));
                        }
                        currentCharacter[std::make_tuple(primaryBlockId, nucPosition, nucGapPosition)] = newVal;
                    } else {
                        if(presentBlocks[primaryBlockId]) {
                            char oldVal = seqChar[std::make_tuple((*node)->parent->identifier,primaryBlockId, nucPosition, -1)];
                            // char oldVal = currentCharacter[std::make_tuple(primaryBlockId, nucPosition, -1)];
                            if(oldVal == '-' || oldVal == 'x') {
                                std::cout << "NOT ACTUALLY A SUBSTITUTION" << std::endl;
                            }
                            // if(node == path.rend()-1)
                            currentNodeMutations.push_back(std::make_pair(primaryBlockId, std::make_tuple('S', panMATCoordinateToGlobal[std::make_tuple(primaryBlockId, nucPosition, -1)], oldVal, newVal, isGapCoordinate[std::make_tuple(primaryBlockId, nucPosition, -1)])));
                            // nodeMutations[u.first].push_back(std::make_tuple('S', panMATCoordinateToGlobal[std::make_tuple(primaryBlockId, nucPosition, -1)], oldVal, newVal));
                        }
                        currentCharacter[std::make_tuple(primaryBlockId, nucPosition, -1)] = newVal;
                    }
                }
            }
        }

        for(auto mut: currentNodeMutations) {
            if(presentBlocks.find(mut.first) != presentBlocks.end()) {
                nodeMutations[u.first].push_back(mut.second);
            }
        }
    });

    for(auto& u: nodeMutations) {
        // print all substitutions first
        fout << "Substitutions:\t";
        fout << u.first << '\t';
        for(auto v: u.second) {
            if(std::get<0>(v) == 'S') {
                fout << " > " << (std::get<4>(v) ? "g" : "") << std::get<2>(v) << std::get<1>(v)+1 << std::get<3>(v);
            }
        }
        fout << '\n';

        fout << "Insertions:\t";
        fout << u.first << '\t';
        // print insertions
        for(auto v: u.second) {
            if(std::get<0>(v) == 'I') {
                fout << " > " << (std::get<4>(v) ? "g" : "") << std::get<1>(v)+1 << std::get<3>(v);
            }
        }
        fout << '\n';

        fout << "Deletions:\t";
        fout << u.first << '\t';
        // print deletions
        for(auto v: u.second) {
            if(std::get<0>(v) == 'D') {
                fout << " > " << (std::get<4>(v) ? "g" : "") << std::get<1>(v)+1 << std::get<2>(v);
            }
        }
        fout << '\n';
    }

}

struct tuple_hash {
    template <class T1, class T2, class T3>
    std::size_t operator() (const std::tuple<T1, T2, T3>& tuple) const {
        auto hash1 = std::hash<T1>{}(std::get<0>(tuple));
        auto hash2 = std::hash<T2>{}(std::get<1>(tuple));
        auto hash3 = std::hash<T3>{}(std::get<2>(tuple));
        return hash1 ^ hash2 ^ hash3;
    }
};

struct tuple_equal {
    template <class T1, class T2, class T3>
    bool operator() (const std::tuple<T1, T2, T3>& lhs, const std::tuple<T1, T2, T3>& rhs) const {
        return lhs == rhs;
    }
};

void printMutationsNewHelper(panmanUtils::Node* node, std::unordered_map<std::tuple<int, int, int>, size_t, tuple_hash, tuple_equal>refPanMATToGlobalCoord, std::string& foutHelp) {

    if (node == nullptr) {
        fprintf(stderr, "Node is null\n");
        return;
    }
    if (node->nucMutation.size() == 0) {
        return;
    }
    foutHelp += node->identifier + ":\t";
    auto mutation = node->nucMutation;
    for (int i=0; i<mutation.size(); i++) {
        int32_t primaryBlockId = mutation[i].primaryBlockId;
        int32_t nucPosition = mutation[i].nucPosition;
        int32_t nucGapPosition = mutation[i].nucGapPosition;
        uint32_t type = (mutation[i].mutInfo & 0x7);
        // if (type != panmanUtils::NucMutationType::NSNPS && type != panmanUtils::NucMutationType::NS) {
        //     continue;
        // }
        char nucType;
        switch (type) {
            case panmanUtils::NucMutationType::NS:
                nucType = 'S';
                break;
            case panmanUtils::NucMutationType::NI:
                nucType = 'I';
                break;
            case panmanUtils::NucMutationType::ND:
                nucType = 'D';
                break;
            case panmanUtils::NucMutationType::NSNPS:
                nucType = 'S';
                break;
            case panmanUtils::NucMutationType::NSNPD:
                nucType = 'D';
                break;
            case panmanUtils::NucMutationType::NSNPI:
                nucType = 'I';
                break;
            default:
                nucType = 'E';
                break;
        }

        int len = ((mutation[i].mutInfo) >> 4);

        
        for (auto j=0;j<len;j++){
            size_t globalCoord;
            // if (nucGapPosition == -1){
            //     globalCoord = refPanMATToGlobalCoord[std::make_tuple(primaryBlockId, nucPosition+j, nucGapPosition)];
            //     foutHelp += nucType;
            // } else {
            //     globalCoord = refPanMATToGlobalCoord[std::make_tuple(primaryBlockId, nucPosition, nucGapPosition+j)];
            //     foutHelp += "g" + nucType;
            // }
            globalCoord = nucPosition+j;
            char newVal = panmanUtils::getNucleotideFromCode(((node->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
            if (newVal != 'N') {
                foutHelp += nucType;
                foutHelp += std::to_string(globalCoord) + ",";
            }
        }

        
    }

    foutHelp += "\n";
}

void panmanUtils::Tree::printMutationsNew(std::ostream& fout, std::string &referenceString) {
    // Get root sequence
    sequence_t rootSequence;
    blockExists_t rootBlockExists;
    blockStrand_t rootBlockStrand;
    getSequenceFromReference(rootSequence, rootBlockExists, rootBlockStrand, root->identifier);

    // Get reference coordinate
    std::unordered_map<std::tuple<int, int, int>, size_t, tuple_hash, tuple_equal> refPanMATToGlobalCoord;
    size_t refCtr = 0;
    size_t MSACtr = 0;
    for(size_t i = 0; i < rootSequence.size(); i++) {
        if(rootBlockExists[i].first) {
            if(rootBlockStrand[i].first) {
                for(size_t j = 0; j < rootSequence[i].first.size(); j++) {
                    for(size_t k = 0; k < rootSequence[i].first[j].second.size(); k++) {
                        refPanMATToGlobalCoord[std::make_tuple(i,j,k)] = refCtr;
                        if(referenceString[MSACtr] != '-' && referenceString[MSACtr] != 'x') {
                            refCtr++;
                        }
                        MSACtr++;
                    }
                    refPanMATToGlobalCoord[std::make_tuple(i,j,-1)] = refCtr;
                    if(referenceString[MSACtr] != '-' && referenceString[MSACtr] != 'x') {
                        refCtr++;
                    }
                    MSACtr++;
                }
            } else {
                for(size_t j = rootSequence[i].first.size() - 1; j + 1 > 0; j--) {
                    refPanMATToGlobalCoord[std::make_tuple(i,j,-1)] = refCtr;
                    if(referenceString[MSACtr] != '-' && referenceString[MSACtr] != 'x') {
                        refCtr++;
                    }
                    MSACtr++;
                    for(size_t k = rootSequence[i].first[j].second.size() - 1; k + 1 > 0; k--) {
                        refPanMATToGlobalCoord[std::make_tuple(i,j,k)] = refCtr;
                        if(referenceString[MSACtr] != '-' && referenceString[MSACtr] != 'x') {
                            refCtr++;
                        }
                        MSACtr++;
                    }
                }
            }
        } else {
            if(rootBlockStrand[i].first) {
                for(size_t j = 0; j < rootSequence[i].first.size(); j++) {
                    for(size_t k = 0; k < rootSequence[i].first[j].second.size(); k++) {
                        refPanMATToGlobalCoord[std::make_tuple(i,j,k)] = refCtr;
                        MSACtr++;
                    }
                    refPanMATToGlobalCoord[std::make_tuple(i,j,-1)] = refCtr;
                    MSACtr++;
                }
            } else {
                for(size_t j = rootSequence[i].first.size() - 1; j + 1 > 0; j--) {
                    refPanMATToGlobalCoord[std::make_tuple(i,j,-1)] = refCtr;
                    MSACtr++;
                    for(size_t k = rootSequence[i].first[j].second.size() - 1; k + 1 > 0; k--) {
                        refPanMATToGlobalCoord[std::make_tuple(i,j,k)] = refCtr;
                        MSACtr++;
                    }
                }
            }
        }
    }
    std::unordered_map< std::string, std::mutex > nodeMutexes;
    for(auto u: allNodes) {
        nodeMutexes[u.first];
    }    
    // for (auto node: allNodes) {
    tbb::parallel_for_each(allNodes, [&](auto node) {
        std::string foutHelp = "";
        printMutationsNewHelper(node.second, refPanMATToGlobalCoord, foutHelp);
        nodeMutexes[node.first].lock();
        fout << foutHelp;
        nodeMutexes[node.first].unlock();
    });
    
}


void panmanUtils::Tree::printMutationsNew(std::ostream& fout, std::vector<std::string>& nodesReq, std::string& referenceString) {

    // Get root sequence
    sequence_t rootSequence;
    blockExists_t rootBlockExists;
    blockStrand_t rootBlockStrand;
    getSequenceFromReference(rootSequence, rootBlockExists, rootBlockStrand, root->identifier);

    // Get reference coordinate
    tbb::concurrent_map< std::tuple< int, int, int >, size_t > refPanMATToGlobalCoord;
    size_t refCtr = 0;
    size_t MSACtr = 0;
    for(size_t i = 0; i < rootSequence.size(); i++) {
        if(rootBlockExists[i].first) {
            if(rootBlockStrand[i].first) {
                for(size_t j = 0; j < rootSequence[i].first.size(); j++) {
                    for(size_t k = 0; k < rootSequence[i].first[j].second.size(); k++) {
                        refPanMATToGlobalCoord[std::make_tuple(i,j,k)] = refCtr;
                        if(referenceString[MSACtr] != '-' && referenceString[MSACtr] != 'x') {
                            refCtr++;
                        }
                        MSACtr++;
                    }
                    refPanMATToGlobalCoord[std::make_tuple(i,j,-1)] = refCtr;
                    if(referenceString[MSACtr] != '-' && referenceString[MSACtr] != 'x') {
                        refCtr++;
                    }
                    MSACtr++;
                }
            } else {
                for(size_t j = rootSequence[i].first.size() - 1; j + 1 > 0; j--) {
                    refPanMATToGlobalCoord[std::make_tuple(i,j,-1)] = refCtr;
                    if(referenceString[MSACtr] != '-' && referenceString[MSACtr] != 'x') {
                        refCtr++;
                    }
                    MSACtr++;
                    for(size_t k = rootSequence[i].first[j].second.size() - 1; k + 1 > 0; k--) {
                        refPanMATToGlobalCoord[std::make_tuple(i,j,k)] = refCtr;
                        if(referenceString[MSACtr] != '-' && referenceString[MSACtr] != 'x') {
                            refCtr++;
                        }
                        MSACtr++;
                    }
                }
            }
        } else {
            if(rootBlockStrand[i].first) {
                for(size_t j = 0; j < rootSequence[i].first.size(); j++) {
                    for(size_t k = 0; k < rootSequence[i].first[j].second.size(); k++) {
                        refPanMATToGlobalCoord[std::make_tuple(i,j,k)] = refCtr;
                        MSACtr++;
                    }
                    refPanMATToGlobalCoord[std::make_tuple(i,j,-1)] = refCtr;
                    MSACtr++;
                }
            } else {
                for(size_t j = rootSequence[i].first.size() - 1; j + 1 > 0; j--) {
                    refPanMATToGlobalCoord[std::make_tuple(i,j,-1)] = refCtr;
                    MSACtr++;
                    for(size_t k = rootSequence[i].first[j].second.size() - 1; k + 1 > 0; k--) {
                        refPanMATToGlobalCoord[std::make_tuple(i,j,k)] = refCtr;
                        MSACtr++;
                    }
                }
            }
        }
    }


    std::cout << refPanMATToGlobalCoord.size() << std::endl;    

    tbb::concurrent_map< std::tuple< int, int, int >, size_t > panMATCoordinateToGlobal;
    tbb::concurrent_map< std::tuple< int, int, int >, char > rootCurrentCharacter;

    tbb::concurrent_unordered_map< size_t, bool > rootPresentBlocks;

    tbb::concurrent_map< std::tuple< int, int, int >, bool > isGapCoordinate;

    // convert PanMAT coordinate to global reference coordinate
    size_t rootCtr = 0;
    for(size_t i = 0; i < rootSequence.size(); i++) {
        if(rootBlockExists[i].first) {
            rootPresentBlocks[i] = true;
            if(rootBlockStrand[i].first) {
                for(size_t j = 0; j < rootSequence[i].first.size(); j++) {
                    for(size_t k = 0; k < rootSequence[i].first[j].second.size(); k++) {
                        panMATCoordinateToGlobal[std::make_tuple(i,j,k)] = rootCtr;
                        if(rootSequence[i].first[j].second[k] != '-' && rootSequence[i].first[j].second[k] != 'x') {
                            rootCurrentCharacter[std::make_tuple(i,j,k)] = rootSequence[i].first[j].second[k];
                            isGapCoordinate[std::make_tuple(i,j,k)] = false;
                            rootCtr++;
                        } else {
                            rootCurrentCharacter[std::make_tuple(i,j,k)] = '-';
                            isGapCoordinate[std::make_tuple(i,j,k)] = true;
                        }
                    }
                    panMATCoordinateToGlobal[std::make_tuple(i,j,-1)] = rootCtr;
                    if(rootSequence[i].first[j].first != '-' && rootSequence[i].first[j].first != 'x') {
                        rootCurrentCharacter[std::make_tuple(i,j,-1)] = rootSequence[i].first[j].first;
                        isGapCoordinate[std::make_tuple(i,j,-1)] = false;
                        rootCtr++;
                    } else {
                        rootCurrentCharacter[std::make_tuple(i,j,-1)] = '-';
                        isGapCoordinate[std::make_tuple(i,j,-1)] = true;
                    }
                }
            } else {
                for(size_t j = rootSequence[i].first.size() - 1; j + 1 > 0; j--) {
                    panMATCoordinateToGlobal[std::make_tuple(i,j,-1)] = rootCtr;
                    if(rootSequence[i].first[j].first != '-' && rootSequence[i].first[j].first != 'x') {
                        rootCurrentCharacter[std::make_tuple(i,j,-1)] = rootSequence[i].first[j].first;
                        isGapCoordinate[std::make_tuple(i,j,-1)] = false;
                        rootCtr++;
                    } else {
                        rootCurrentCharacter[std::make_tuple(i,j,-1)] = '-';
                        isGapCoordinate[std::make_tuple(i,j,-1)] = true;
                    }
                    for(size_t k = rootSequence[i].first[j].second.size() - 1; k + 1 > 0; k--) {
                        panMATCoordinateToGlobal[std::make_tuple(i,j,k)] = rootCtr;
                        if(rootSequence[i].first[j].second[k] != '-' && rootSequence[i].first[j].second[k] != 'x') {
                            rootCurrentCharacter[std::make_tuple(i,j,k)] = rootSequence[i].first[j].second[k];
                            // if(rootCtr == 240) {
                            //     std::cout << rootSequence[i].first[j].first << " " << i << " " << j << " " << k << std::endl;
                            // }
                            isGapCoordinate[std::make_tuple(i,j,k)] = false;
                            rootCtr++;
                        } else {
                            rootCurrentCharacter[std::make_tuple(i,j,k)] = '-';
                            isGapCoordinate[std::make_tuple(i,j,k)] = true;
                        }
                    }
                }
            }
        } else {
            rootPresentBlocks[i] = false;
            if(rootBlockStrand[i].first) {
                for(size_t j = 0; j < rootSequence[i].first.size(); j++) {
                    for(size_t k = 0; k < rootSequence[i].first[j].second.size(); k++) {
                        panMATCoordinateToGlobal[std::make_tuple(i,j,k)] = rootCtr;
                        rootCurrentCharacter[std::make_tuple(i,j,k)] = '-';
                        isGapCoordinate[std::make_tuple(i,j,k)] = false;
                    }
                    panMATCoordinateToGlobal[std::make_tuple(i,j,-1)] = rootCtr;
                    rootCurrentCharacter[std::make_tuple(i,j,-1)] = '-';
                    isGapCoordinate[std::make_tuple(i,j,-1)] = false;
                }
            } else {
                for(size_t j = rootSequence[i].first.size() - 1; j + 1 > 0; j--) {
                    panMATCoordinateToGlobal[std::make_tuple(i,j,-1)] = rootCtr;
                    rootCurrentCharacter[std::make_tuple(i,j,-1)] = '-';
                    isGapCoordinate[std::make_tuple(i,j,-1)] = false;
                    for(size_t k = rootSequence[i].first[j].second.size() - 1; k + 1 > 0; k--) {
                        panMATCoordinateToGlobal[std::make_tuple(i,j,k)] = rootCtr;
                        rootCurrentCharacter[std::make_tuple(i,j,k)] = '-';
                        isGapCoordinate[std::make_tuple(i,j,k)] = false;
                    }
                }
            }
        }
    }

    // tbb::parallel_for_each(nodesReq, [&](auto w) {
    int c=0;
    for (auto w: nodesReq) {
        std::cout << c++ << std::endl;
        if (allNodes.find(w) == allNodes.end()) {
            std::cerr << "Could not find node " << w << std::endl;
            continue;
        }
        std::pair<std::string, Node*> u = std::make_pair(w, allNodes[w]);

        tbb::concurrent_map< std::tuple< std::string, int, int, int >, char > seqChar;

        Node* it = u.second;
        std::vector< panmanUtils::Node* > path;

        while(it != root) {
            path.push_back(it);
            it = it->parent;
        }
        path.push_back(root);
        tbb::concurrent_map< std::tuple< int, int, int >, char > currentCharacter = rootCurrentCharacter;
        tbb::concurrent_unordered_map< size_t, bool > presentBlocks = rootPresentBlocks;
        for(auto node = path.rbegin(); node != path.rend(); node++) {
            for(auto mutation: (*node)->blockMutation) {
                int32_t primaryBlockId = mutation.primaryBlockId;
                if(rootPresentBlocks.find(primaryBlockId) == rootPresentBlocks.end()) {
                    continue;
                }

                bool type = mutation.blockMutInfo;
                bool inversion = mutation.inversion;
                if(type == 1) {
                    // insertion
                    presentBlocks[primaryBlockId] = true;
                    if(inversion) {
                        std::cout << "INVERTED BLOCK FOUND" << std::endl;
                    }
                } else {
                    if(inversion) {
                        // This means that this is not a deletion, but instead an inversion
                        std::cout << "INVERSION FOUND" << std::endl;
                    } else {
                        // Actually a deletion
                        presentBlocks[primaryBlockId] = false;
                    }
                }
            }
        }

        // std::cout << "Seq char" << std::endl;
        for(auto node = path.rbegin(); node != path.rend(); node++) {
            sequence_t st;
            blockExists_t bt;
            blockStrand_t bst;
            getSequenceFromReference(st, bt, bst, (*node)->identifier);

            for(size_t i = 0; i < st.size(); i++) {
                if(bst[i].first) {
                    for(size_t j = 0; j < st[i].first.size(); j++) {
                        for(size_t k = 0; k < st[i].first[j].second.size(); k++) {
                            if(st[i].first[j].second[k] != '-' && st[i].first[j].second[k] != 'x') {
                                seqChar[std::make_tuple(u.first,i,j,k)] = st[i].first[j].second[k];
                            } else {
                                seqChar[std::make_tuple(u.first,i,j,k)] = '-';
                            }
                        }
                        if(st[i].first[j].first != '-' && st[i].first[j].first != 'x' && bt[i].first) {
                            seqChar[std::make_tuple(u.first,i,j,-1)] = st[i].first[j].first;
                        } else {
                            seqChar[std::make_tuple(u.first,i,j,-1)] = '-';
                        }
                    }
                } else {
                    for(size_t j = st[i].first.size() - 1; j + 1 > 0; j--) {
                        if(st[i].first[j].first != '-' && st[i].first[j].first != 'x') {
                            seqChar[std::make_tuple(u.first,i,j,-1)] = rootSequence[i].first[j].first;
                        } else {
                            seqChar[std::make_tuple(u.first,i,j,-1)] = '-';
                        }
                        for(size_t k = st[i].first[j].second.size() - 1; k + 1 > 0; k--) {
                            if(st[i].first[j].second[k] != '-' && st[i].first[j].second[k] != 'x') {
                                seqChar[std::make_tuple(u.first,i,j,k)] = st[i].first[j].second[k];
                            } else {
                                seqChar[std::make_tuple(u.first,i,j,k)] = '-';
                            }
                        }
                    }
                }
            }
            // std::cout << seqChar.size() << std::endl;
        }

        std::vector< std::pair< std::string, std::tuple< std::string, std::string > > > currentNodeMutations;
        for(auto node = path.rbegin(); node != path.rend(); node++) {
        // for (auto omega = 0; omega < 1; omega++) {
            // std::cout << (*node)->identifier << std::endl;
            // auto node = &path[omega];
            if ((*node)->identifier == root->identifier) continue;
            for(size_t i = 0; i < (*node)->nucMutation.size(); i++) {
                int32_t primaryBlockId = (*node)->nucMutation[i].primaryBlockId;
                if(rootPresentBlocks.find(primaryBlockId) == rootPresentBlocks.end()) {
                    continue;
                }

                int32_t nucPosition = (*node)->nucMutation[i].nucPosition;
                int32_t nucGapPosition = (*node)->nucMutation[i].nucGapPosition;
                uint32_t type = ((*node)->nucMutation[i].mutInfo & 0x7);
                char newVal = '-';
                // std::cout << "mutation count: " << i << " " <<
                //                 type << " " << nucPosition << " " << nucGapPosition << " " << (((*node)->nucMutation[i].mutInfo) >> 4) <<std::endl;

                if(type < 3) {
                    int len = (((*node)->nucMutation[i].mutInfo) >> 4);

                    if(type == panmanUtils::NucMutationType::NS) {
                        if(nucGapPosition != -1) {
                            for(int j = 0; j < len; j++) {
                                newVal = seqChar[std::make_tuple((*node)->identifier,primaryBlockId, nucPosition, nucGapPosition+j)];
                                if(presentBlocks[primaryBlockId]) {
                                    char oldVal = seqChar[std::make_tuple((*node)->parent->identifier,primaryBlockId, nucPosition, nucGapPosition+j)];
                                    if(oldVal == '-' || oldVal == 'x') {
                                        continue;
                                    }
                                    std::string currMutType = "";
                                    if (isGapCoordinate[std::make_tuple(primaryBlockId, nucPosition, nucGapPosition+j)]) currMutType += "g";
                                    currMutType += "S";
                                    std::string currMut = oldVal + std::to_string(refPanMATToGlobalCoord[std::make_tuple(primaryBlockId, nucPosition, nucGapPosition+j)]) + newVal;
                                    currentNodeMutations.push_back(std::make_pair((*node)->identifier, std::make_tuple(currMutType, currMut)));
                                    // nodeMutations[u.first].push_back(std::make_tuple('S', panMATCoordinateToGlobal[std::make_tuple(primaryBlockId, nucPosition, nucGapPosition+j)], oldVal, newVal));
                                }
                            }
                        } else {
                            for(int j = 0; j < len; j++) {
                                newVal = seqChar[std::make_tuple((*node)->identifier,primaryBlockId, nucPosition+j, -1)];
                                if(presentBlocks[primaryBlockId]) {
                                    char oldVal = seqChar[std::make_tuple((*node)->parent->identifier,primaryBlockId, nucPosition+j, -1)];
                                    // }
                                    if(oldVal == '-' || oldVal == 'x') {
                                        continue;
                                    }
                                    std::string currMutType = "";
                                    if (isGapCoordinate[std::make_tuple(primaryBlockId, nucPosition+j, -1)]) currMutType += "g";
                                    currMutType += "S";
                                    std::string currMut = oldVal + std::to_string(refPanMATToGlobalCoord[std::make_tuple(primaryBlockId, nucPosition + j, -1)]) + newVal;
                                    currentNodeMutations.push_back(std::make_pair((*node)->identifier, std::make_tuple(currMutType, currMut)));
                                    // nodeMutations[u.first].push_back(std::make_tuple('S', panMATCoordinateToGlobal[std::make_tuple(primaryBlockId, nucPosition + j, -1)], oldVal, newVal));
                                }
                            }
                        }
                    } else if(type == panmanUtils::NucMutationType::NI) {
                        if(nucGapPosition != -1) {
                            for(int j = 0; j < len; j++) {
                                newVal = seqChar[std::make_tuple((*node)->identifier,primaryBlockId, nucPosition, nucGapPosition + j)];
                                std::string currMutType = "";
                                if (isGapCoordinate[std::make_tuple(primaryBlockId, nucPosition, nucGapPosition + j)]) currMutType += "g";
                                currMutType += "I";
                                std::string currMut = "-" + std::to_string(refPanMATToGlobalCoord[std::make_tuple(primaryBlockId, nucPosition, nucGapPosition + j)]) + newVal;
                                currentNodeMutations.push_back(std::make_pair((*node)->identifier, std::make_tuple(currMutType, currMut)));
                            }
                        } else {
                            for(int j = 0; j < len; j++) {
                                // std::cout << (*node)->identifier << " " << primaryBlockId << " " << (nucPosition+j) << " " << (seqChar.find(std::make_tuple((*node)->identifier,primaryBlockId, nucPosition+j, -1)) == seqChar.end()) <<
                                // " " << (isGapCoordinate.find(std::make_tuple(primaryBlockId, nucPosition + j, -1)) == isGapCoordinate.end()) << std::endl;
                                newVal = seqChar[std::make_tuple((*node)->identifier,primaryBlockId, nucPosition+j, -1)];
                                std::string currMutType = "";
                                if (isGapCoordinate[std::make_tuple(primaryBlockId, nucPosition + j, -1)]) currMutType += "g";
                                currMutType += "I";
                                std::string currMut = "-" + std::to_string(refPanMATToGlobalCoord[std::make_tuple(primaryBlockId, nucPosition + j, -1)]) + newVal;
                                // std::cout << currMutType << " " << currMut << std::endl;
                                currentNodeMutations.push_back(std::make_pair((*node)->identifier, std::make_tuple(currMutType, currMut)));

                            }
                        }
                    } else if(type == panmanUtils::NucMutationType::ND) {
                        if(nucGapPosition != -1) {
                            for(int j = 0; j < len; j++) {
                                char oldVal = seqChar[std::make_tuple((*node)->parent->identifier,primaryBlockId, nucPosition, nucGapPosition + j)];
                                std::string currMutType = "";
                                if (isGapCoordinate[std::make_tuple(primaryBlockId, nucPosition, nucGapPosition + j)]) currMutType += "g";
                                currMutType += "D";
                                std::string currMut = oldVal + std::to_string(refPanMATToGlobalCoord[std::make_tuple(primaryBlockId, nucPosition, nucGapPosition + j)]) + "-";
                                currentNodeMutations.push_back(std::make_pair((*node)->identifier, std::make_tuple(currMutType, currMut)));
                                // }
                            }
                        } else {
                            for(int j = 0; j < len; j++) {
                                char oldVal = seqChar[std::make_tuple((*node)->parent->identifier,primaryBlockId, nucPosition + j, nucGapPosition)];
                                std::string currMutType = "";
                                if (isGapCoordinate[std::make_tuple(primaryBlockId, nucPosition + j, -1)]) currMutType += "g";
                                currMutType += "D";
                                std::string currMut = oldVal + std::to_string(refPanMATToGlobalCoord[std::make_tuple(primaryBlockId, nucPosition + j, -1)]) + "-";
                                currentNodeMutations.push_back(std::make_pair((*node)->identifier, std::make_tuple(currMutType, currMut)));
                                // }
                            }
                        }
                    }
                } else if(type == panmanUtils::NucMutationType::NSNPS) {
                    newVal = seqChar[std::make_tuple((*node)->identifier,primaryBlockId, nucPosition, nucGapPosition)];
                    if(nucGapPosition != -1) {
                        if(presentBlocks[primaryBlockId]) {
                            char oldVal = seqChar[std::make_tuple((*node)->parent->identifier,primaryBlockId, nucPosition, nucGapPosition)];
                            if(oldVal == '-' || oldVal == 'x') {
                                std::cout << "NOT ACTUALLY A SUBSTITUTION" << std::endl;
                            }
                            std::string currMutType = "";
                            if (isGapCoordinate[std::make_tuple(primaryBlockId, nucPosition, nucGapPosition)]) currMutType += "g";
                            currMutType += "S";
                            std::string currMut = oldVal + std::to_string(refPanMATToGlobalCoord[std::make_tuple(primaryBlockId, nucPosition, nucGapPosition)]) + newVal;
                            currentNodeMutations.push_back(std::make_pair((*node)->identifier, std::make_tuple(currMutType, currMut)));
                        }
                    } else {
                        if(presentBlocks[primaryBlockId]) {
                            char oldVal = seqChar[std::make_tuple((*node)->parent->identifier,primaryBlockId, nucPosition, -1)];
                            if(oldVal == '-' || oldVal == 'x') {
                                std::cout << "NOT ACTUALLY A SUBSTITUTION" << std::endl;
                            }
                            std::string currMutType = "";
                            if (isGapCoordinate[std::make_tuple(primaryBlockId, nucPosition, -1)]) currMutType += "g";
                            currMutType += "S";
                            std::string currMut = oldVal + std::to_string(refPanMATToGlobalCoord[std::make_tuple(primaryBlockId, nucPosition, -1)]) + newVal;
                            currentNodeMutations.push_back(std::make_pair((*node)->identifier, std::make_tuple(currMutType, currMut)));
                        }
                    }
                }
            }
            
        }
        
        // std::cout << "Writing mutations" << std::endl;
        auto nodeName = root->identifier;
        fout << u.first << "\tD:\n";
        for(auto mut: currentNodeMutations) {
            // if(presentBlocks.find(mut.first) != presentBlocks.end()) {
                if (nodeName != mut.first) { 
                    fout << "\n\t>" + mut.first + "\t";
                    nodeName = mut.first;
                }
                if (std::get<0>(mut.second) == "D" || std::get<0>(mut.second) == "gD") {
                    fout << std::get<1>(mut.second) << "\t";
                }
            // }
        }
        fout << '\n';

    // });
    }

    

    
        
    

    // for(auto& u: nodeMutations) {
    //     // print all substitutions first
    //     // fout << "Substitutions:\t";
    //     fout << u.first << '\t';
    //     for(auto v: u.second) {
    //         if(std::get<0>(v) == 'S') {
    //             fout << " > " << (std::get<4>(v) ? "g" : "") << std::get<2>(v) << std::get<1>(v)+1 << std::get<3>(v);
    //         }
    //     }
    //     fout << '\n';

        // fout << "Insertions:\t";
        // fout << u.first << '\t';
        // // print insertions
        // for(auto v: u.second) {
        //     if(std::get<0>(v) == 'I') {
        //         fout << " > " << (std::get<4>(v) ? "g" : "") << std::get<1>(v)+1 << std::get<3>(v);
        //     }
        // }
        // fout << '\n';

        // fout << "Deletions:\t";
        // fout << u.first << '\t';
        // // print deletions
        // for(auto v: u.second) {
        //     if(std::get<0>(v) == 'D') {
        //         fout << " > " << (std::get<4>(v) ? "g" : "") << std::get<1>(v)+1 << std::get<2>(v);
        //     }
        // }
        // fout << '\n';
    // }

}
