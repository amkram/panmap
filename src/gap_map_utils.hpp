#pragma once

#include <map>
#include <vector>
#include <utility>
#include <algorithm>

namespace gap_map {

// Templated gap map operations - works with both int64_t and uint64_t
// T = coordinate type (int64_t or uint64_t)

template<typename T>
void updateGapMapStep(
    std::map<T, T>& gapMap,
    T start,
    T end,
    bool toGap,
    std::vector<std::pair<bool, std::pair<T, T>>>& backtrack,
    std::vector<std::pair<bool, std::pair<T, T>>>& gapMapUpdates,
    bool recordGapMapUpdates = true
) {
    auto rightIt = gapMap.upper_bound(start);
    auto leftIt = (rightIt == gapMap.begin()) ? gapMap.end() : std::prev(rightIt);

    bool rightItExists = rightIt != gapMap.end();
    bool leftItExists = leftIt != gapMap.end();

    if (toGap) {
        // add gap range
        if (gapMap.empty()) {
            gapMap[start] = end;
            backtrack.emplace_back(true, std::make_pair(start, end));
            if (recordGapMapUpdates) {
                gapMapUpdates.emplace_back(false, std::make_pair(start, end));
            }
            return;
        }
        
        decltype(rightIt) curIt;

        // curIt starts outside of any range
        if (!leftItExists || (!rightItExists && start > leftIt->second) || 
            (leftItExists && start > leftIt->second && rightItExists && start < rightIt->first)) {
            if (leftItExists && start == leftIt->second + 1) {
                // 1 base after left range and merge with left
                curIt = leftIt;
                backtrack.emplace_back(false, std::make_pair(curIt->first, curIt->second));
                curIt->second = end;
                if (recordGapMapUpdates) {
                    gapMapUpdates.emplace_back(false, std::make_pair(curIt->first, curIt->second));
                }
            } else {
                // insert new range
                auto tmpIt = gapMap.emplace(start, end);
                curIt = tmpIt.first;
                backtrack.emplace_back(true, std::make_pair(curIt->first, curIt->second));
                if (recordGapMapUpdates) {
                    gapMapUpdates.emplace_back(false, std::make_pair(curIt->first, curIt->second));
                }
            }
        } else {
            curIt = leftIt;
            if (end <= curIt->second) {
                return;
            }
            backtrack.emplace_back(false, std::make_pair(curIt->first, curIt->second));
            curIt->second = end;
            if (recordGapMapUpdates) {
                gapMapUpdates.emplace_back(false, std::make_pair(curIt->first, curIt->second));
            }
        }

        auto nextIt = std::next(curIt);
        while (true) {
            if (nextIt == gapMap.end()) {
                break;
            }

            if (nextIt->second <= curIt->second) {
                auto tmpIt = nextIt;
                nextIt = std::next(nextIt);
                backtrack.emplace_back(false, std::make_pair(tmpIt->first, tmpIt->second));
                if (recordGapMapUpdates) {
                    gapMapUpdates.emplace_back(true, std::make_pair(tmpIt->first, tmpIt->second));
                }
                gapMap.erase(tmpIt);
            } else if (nextIt->first <= end + 1) {
                backtrack.emplace_back(false, std::make_pair(curIt->first, curIt->second));
                curIt->second = nextIt->second;
                backtrack.emplace_back(false, std::make_pair(nextIt->first, nextIt->second));
                if (recordGapMapUpdates) {
                    gapMapUpdates.emplace_back(false, std::make_pair(curIt->first, curIt->second));
                    gapMapUpdates.emplace_back(true, std::make_pair(nextIt->first, nextIt->second));
                }
                gapMap.erase(nextIt);
                break;
            } else {
                break;
            }
        }
    } else {
        // remove gap range
        if (gapMap.empty() || (!leftItExists && end < rightIt->first) || 
            (!rightItExists && start > leftIt->second)) {
            return;
        }

        decltype(rightIt) curIt;
        decltype(rightIt) nextIt;
        if (!leftItExists || (leftItExists && start > leftIt->second && rightItExists && start < rightIt->first)) {
            // curIt starts outside of any range
            curIt = rightIt;

            if (end < curIt->first) {
                return;
            }

            // ends within the curIt range
            if (end <= curIt->second) {
                if (end == curIt->second) {
                    backtrack.emplace_back(false, std::make_pair(curIt->first, curIt->second));
                    if (recordGapMapUpdates) {
                        gapMapUpdates.emplace_back(true, std::make_pair(curIt->first, curIt->second));
                    }
                    gapMap.erase(curIt);
                } else {
                    gapMap[end+1] = curIt->second;
                    backtrack.emplace_back(true, std::make_pair(end+1, curIt->second));
                    backtrack.emplace_back(false, std::make_pair(curIt->first, curIt->second));
                    if (recordGapMapUpdates) {
                        gapMapUpdates.emplace_back(false, std::make_pair(end+1, curIt->second));
                        gapMapUpdates.emplace_back(true, std::make_pair(curIt->first, curIt->second));
                    }
                    gapMap.erase(curIt);
                }
                return;
            } else {
                nextIt = std::next(curIt);
                backtrack.emplace_back(false, std::make_pair(curIt->first, curIt->second));
                if (recordGapMapUpdates) {
                    gapMapUpdates.emplace_back(true, std::make_pair(curIt->first, curIt->second));
                }
                gapMap.erase(curIt);
            }
            
        } else {
            // curIt starts inside of a range
            curIt = leftIt;
            
            if (end <= curIt->second) {
                // contained in the curIt range
                if (start == curIt->first && end == curIt->second) {
                    backtrack.emplace_back(false, std::make_pair(curIt->first, curIt->second));
                    if (recordGapMapUpdates) {
                        gapMapUpdates.emplace_back(true, std::make_pair(curIt->first, curIt->second));
                    }
                    gapMap.erase(curIt);
                } else if (start == curIt->first) {
                    gapMap[end + 1] = curIt->second;
                    backtrack.emplace_back(true, std::make_pair(end+1, curIt->second));
                    backtrack.emplace_back(false, std::make_pair(curIt->first, curIt->second));
                    if (recordGapMapUpdates) {
                        gapMapUpdates.emplace_back(false, std::make_pair(end+1, curIt->second));
                        gapMapUpdates.emplace_back(true, std::make_pair(curIt->first, curIt->second));
                    }
                    gapMap.erase(curIt);
                } else if (end == curIt->second) {
                    backtrack.emplace_back(false, std::make_pair(curIt->first, curIt->second));
                    curIt->second = start - 1;
                    if (recordGapMapUpdates) {
                        gapMapUpdates.emplace_back(false, std::make_pair(curIt->first, curIt->second));
                    }
                } else {
                    gapMap[end + 1] = curIt->second;
                    backtrack.emplace_back(true, std::make_pair(end+1, curIt->second));
                    backtrack.emplace_back(false, std::make_pair(curIt->first, curIt->second));
                    if (recordGapMapUpdates) {
                        gapMapUpdates.emplace_back(false, std::make_pair(end+1, curIt->second));
                        gapMapUpdates.emplace_back(false, std::make_pair(curIt->first, start-1));
                    }
                    curIt->second = start - 1;
                }
                return;
            } else {
                if (start == curIt->first) {
                    nextIt = std::next(curIt);
                    backtrack.emplace_back(false, std::make_pair(curIt->first, curIt->second));
                    if (recordGapMapUpdates) {
                        gapMapUpdates.emplace_back(true, std::make_pair(curIt->first, curIt->second));
                    }
                    gapMap.erase(curIt);
                } else {
                    backtrack.emplace_back(false, std::make_pair(curIt->first, curIt->second));
                    curIt->second = start - 1;
                    if (recordGapMapUpdates) {
                        gapMapUpdates.emplace_back(false, std::make_pair(curIt->first, curIt->second));
                    }
                    nextIt = std::next(curIt);
                }
            }
        }

        while (true) {
            if (nextIt == gapMap.end()) {
                break;
            }

            if (nextIt->first > end) {
                break;
            } else if (nextIt->second <= end) {
                auto tmpIt = nextIt;
                nextIt = std::next(nextIt);
                backtrack.emplace_back(false, std::make_pair(tmpIt->first, tmpIt->second));
                if (recordGapMapUpdates) {
                    gapMapUpdates.emplace_back(true, std::make_pair(tmpIt->first, tmpIt->second));
                }
                gapMap.erase(tmpIt);
            } else {
                gapMap[end + 1] = nextIt->second;
                backtrack.emplace_back(true, std::make_pair(end+1, nextIt->second));
                backtrack.emplace_back(false, std::make_pair(nextIt->first, nextIt->second));
                if (recordGapMapUpdates) {
                    gapMapUpdates.emplace_back(false, std::make_pair(end+1, nextIt->second));
                    gapMapUpdates.emplace_back(true, std::make_pair(nextIt->first, nextIt->second));
                }
                gapMap.erase(nextIt);
                break;
            }
        }
    }
}

// Convenience overload that takes update as a pair (for genotyping.cpp compatibility)
template<typename T>
void updateGapMapStep(
    std::map<T, T>& gapMap,
    const std::pair<bool, std::pair<T, T>>& update,
    std::vector<std::pair<bool, std::pair<T, T>>>& backtrack,
    std::vector<std::pair<bool, std::pair<T, T>>>& gapMapUpdates,
    bool recordGapMapUpdates = true
) {
    updateGapMapStep(gapMap, update.second.first, update.second.second, update.first, 
                     backtrack, gapMapUpdates, recordGapMapUpdates);
}

template<typename T>
void invertGapMap(
    std::map<T, T>& gapMap,
    const std::pair<T, T>& invertRange,
    std::vector<std::pair<bool, std::pair<T, T>>>& backtrack,
    std::vector<std::pair<bool, std::pair<T, T>>>& gapMapUpdates
) {
    const auto& [start, end] = invertRange;

    auto rightIt = gapMap.upper_bound(start);
    auto leftIt = (rightIt == gapMap.begin()) ? gapMap.end() : std::prev(rightIt);

    bool rightItExists = rightIt != gapMap.end();
    bool leftItExists = leftIt != gapMap.end();

    // completely inside or outside a gap range -> do nothing
    if (gapMap.empty() ||
        (!leftItExists && end < rightIt->first) ||
        (!rightItExists && start > leftIt->second)) {
        return;
    }

    // Early exit for range completely between two gap ranges (index_single_mode optimization)
    if (!leftItExists || (leftItExists && start > leftIt->second && rightItExists && start < rightIt->first)) {
        if (rightItExists && end < rightIt->first) {
            return;
        }
    }
    
    //                    gaps            beg      end
    std::vector<std::pair<bool, std::pair<int64_t, int64_t>>> blockRuns;
    if (!leftItExists || (leftItExists && start > leftIt->second && rightItExists && start < rightIt->first)) {
        // start outside of a range
        auto curIt = rightIt;
        blockRuns.emplace_back(false, std::make_pair(static_cast<int64_t>(start), 
                                                      static_cast<int64_t>(curIt->first) - 1));
        if (end <= curIt->second) {
            blockRuns.emplace_back(true, std::make_pair(blockRuns.back().second.second + 1, 
                                                        static_cast<int64_t>(end)));
        } else {
            blockRuns.emplace_back(true, std::make_pair(blockRuns.back().second.second + 1, 
                                                        static_cast<int64_t>(curIt->second)));
            curIt = std::next(curIt);
            while (true) {
                if (curIt == gapMap.end()) {
                    blockRuns.emplace_back(false, std::make_pair(blockRuns.back().second.second + 1, 
                                                                  static_cast<int64_t>(end)));
                    break;
                } else if (end < curIt->first) {
                    blockRuns.emplace_back(false, std::make_pair(blockRuns.back().second.second + 1, 
                                                                  static_cast<int64_t>(end)));
                    break;
                } else if (end > curIt->second) {
                    blockRuns.emplace_back(false, std::make_pair(blockRuns.back().second.second + 1, 
                                                                  static_cast<int64_t>(curIt->first) - 1));
                    blockRuns.emplace_back(true, std::make_pair(static_cast<int64_t>(curIt->first), 
                                                                static_cast<int64_t>(curIt->second)));
                    curIt = std::next(curIt);
                } else {
                    blockRuns.emplace_back(false, std::make_pair(blockRuns.back().second.second + 1, 
                                                                  static_cast<int64_t>(curIt->first) - 1));
                    blockRuns.emplace_back(true, std::make_pair(static_cast<int64_t>(curIt->first), 
                                                                static_cast<int64_t>(end)));
                    break;
                }
            }
        }
    } else {
        // start inside of a range
        auto curIt = leftIt;
        blockRuns.emplace_back(true, std::make_pair(static_cast<int64_t>(start), 
                                                    static_cast<int64_t>(curIt->second)));
        curIt = std::next(curIt);
        while (true) {
            if (curIt == gapMap.end()) {
                blockRuns.emplace_back(false, std::make_pair(blockRuns.back().second.second + 1, 
                                                              static_cast<int64_t>(end)));
                break;
            } else if (end < curIt->first) {
                blockRuns.emplace_back(false, std::make_pair(blockRuns.back().second.second + 1, 
                                                              static_cast<int64_t>(end)));
                break;
            } else if (end > curIt->second) {
                blockRuns.emplace_back(false, std::make_pair(blockRuns.back().second.second + 1, 
                                                              static_cast<int64_t>(curIt->first) - 1));
                blockRuns.emplace_back(true, std::make_pair(static_cast<int64_t>(curIt->first), 
                                                            static_cast<int64_t>(curIt->second)));
                curIt = std::next(curIt);
            } else {
                blockRuns.emplace_back(false, std::make_pair(blockRuns.back().second.second + 1, 
                                                              static_cast<int64_t>(curIt->first) - 1));
                blockRuns.emplace_back(true, std::make_pair(static_cast<int64_t>(curIt->first), 
                                                            static_cast<int64_t>(end)));
                break;
            }
        }
    }

    int64_t curBeg = blockRuns.front().second.first;
    for (auto it = blockRuns.rbegin(); it != blockRuns.rend(); ++it) {
        int64_t curEnd = curBeg + (it->second.second - it->second.first);
        updateGapMapStep(gapMap, static_cast<T>(curBeg), static_cast<T>(curEnd), it->first, 
                         backtrack, gapMapUpdates, false);
        curBeg = curEnd + 1;
    }
}

template<typename T>
void revertGapMapChanges(
    std::vector<std::pair<bool, std::pair<T, T>>>& backtrack,
    std::map<T, T>& gapMap
) {
    for (auto it = backtrack.rbegin(); it != backtrack.rend(); ++it) {
        const auto& [del, range] = *it;
        if (del) {
            gapMap.erase(range.first);
        } else {
            gapMap[range.first] = range.second;
        }
    }
}

} // namespace gap_map
