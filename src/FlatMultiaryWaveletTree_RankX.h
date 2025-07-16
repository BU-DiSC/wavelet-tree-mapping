#ifndef FLAT_MULTIARY_WAVELET_TREE_H
#define FLAT_MULTIARY_WAVELET_TREE_H

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <algorithm> // for std::min

using namespace std;

class FlatMultiaryWaveletTree {
public:
    const int arity = 256;  // 256-way tree
    int originalSize;
    int treeSize;
    int numLevels;

    uint8_t **dataArr;
    void **rankMatrixArr;

private:
    static constexpr int L0_SPARSE_BUCKET_SIZE = 512;
    static constexpr int L1_SPARSE_BUCKET_SIZE = 1024; //only applicable if the flag is set to True
    //we put a flag for level 1 because creating a sparse array for level1 does not benefit size as much
    bool m_enableSparseL1Rank; // Member flag to control L1 rank strategy

public:
    FlatMultiaryWaveletTree(const int *input_param, int size, bool enableL1SparseRank = false)
        : m_enableSparseL1Rank(enableL1SparseRank) {
        originalSize = size;
        if (originalSize == 0) {
            treeSize = 0;
        } else {
            treeSize = 1;
            while (treeSize < originalSize) {
                treeSize *= arity;
            }
        }

        int *paddedInput = nullptr;
        if (treeSize > 0) {
            paddedInput = new int[treeSize];
            for (int i = 0; i < originalSize; i++) {
                paddedInput[i] = input_param[i];
            }
            for (int i = originalSize; i < treeSize; i++) {
                paddedInput[i] = i;
            }
        }
        int totalLevels = calculateLevel(); 
        numLevels = totalLevels > 0 ? totalLevels - 1 : 0;

        dataArr = nullptr;
        if (numLevels > 0 && treeSize > 0) {
            dataArr = new uint8_t*[numLevels];
            for (int lvl = 0; lvl < numLevels; lvl++) {
                dataArr[lvl] = new uint8_t[treeSize];
                for(int i = 0; i < treeSize; i++){
                    dataArr[lvl][i] = 0;
                }
            }
            if (paddedInput) {
                build_node(0, 0, paddedInput, treeSize, 0, treeSize - 1);
            }
        }
        rankMatrixArr = nullptr;
        if (numLevels > 0 && treeSize > 0) { // check treeSize for rank allocation
            rankMatrixArr = new void*[numLevels];
            for (int k = 0; k < numLevels; ++k) {
                rankMatrixArr[k] = nullptr;
            }

            int numBuckets0 = (treeSize + L0_SPARSE_BUCKET_SIZE - 1) / L0_SPARSE_BUCKET_SIZE;
            if (dataArr != nullptr && dataArr[0] != nullptr) {
                rankMatrixArr[0] = new uint16_t[static_cast<size_t>(numBuckets0) * arity];
                uint16_t* row0 = static_cast<uint16_t*>(rankMatrixArr[0]);
                vector<int> cumulative_l0(arity, 0);
                for (int b = 0; b < numBuckets0; b++) {
                    for (int sym = 0; sym < arity; sym++) {
                        row0[static_cast<size_t>(b) * arity + sym] = static_cast<uint16_t>(cumulative_l0[sym]);
                    }
                    int start = b * L0_SPARSE_BUCKET_SIZE;
                    int end = std::min((b + 1) * L0_SPARSE_BUCKET_SIZE, treeSize);
                    for (int i = start; i < end; i++) {
                        int d = dataArr[0][i];
                        cumulative_l0[d]++;
                    }
                }
            }

            if (numLevels > 1) {
                if (m_enableSparseL1Rank) { // flag to check if level 1 sparse rank is enabled
                    int nodeSize_L1 = numElementsLvl(1);
                    int numMajorSegments_L1 = (nodeSize_L1 > 0) ? (treeSize / nodeSize_L1) : 0;
                    int numSparseBucketsPerMajorSegment = 0;
                    if (nodeSize_L1 > 0) {
                         numSparseBucketsPerMajorSegment = (nodeSize_L1 + L1_SPARSE_BUCKET_SIZE - 1) / L1_SPARSE_BUCKET_SIZE;
                    }

                    if (numMajorSegments_L1 > 0 && numSparseBucketsPerMajorSegment > 0 && dataArr != nullptr && dataArr[1] != nullptr) {
                        rankMatrixArr[1] = new uint16_t[static_cast<size_t>(numMajorSegments_L1) * numSparseBucketsPerMajorSegment * arity];
                        uint16_t* row1_sparse = static_cast<uint16_t*>(rankMatrixArr[1]);
                        for (int majorSegIdx = 0; majorSegIdx < numMajorSegments_L1; ++majorSegIdx) {
                            std::vector<int> cumulative_in_major_segment(arity, 0);
                            int majorSegmentDataStartIdx = majorSegIdx * nodeSize_L1;
                            for (int sparseBlockIdxInMajorSeg = 0; sparseBlockIdxInMajorSeg < numSparseBucketsPerMajorSegment; ++sparseBlockIdxInMajorSeg) {
                                size_t sparseRankBlockEntryBase = (static_cast<size_t>(majorSegIdx) * numSparseBucketsPerMajorSegment + sparseBlockIdxInMajorSeg) * arity;
                                for (int sym = 0; sym < arity; ++sym) {
                                    row1_sparse[sparseRankBlockEntryBase + sym] = static_cast<uint16_t>(cumulative_in_major_segment[sym]);
                                }
                                int currentSparseBlockDataStart = majorSegmentDataStartIdx + sparseBlockIdxInMajorSeg * L1_SPARSE_BUCKET_SIZE;
                                int currentSparseBlockDataEnd = majorSegmentDataStartIdx + (sparseBlockIdxInMajorSeg + 1) * L1_SPARSE_BUCKET_SIZE;
                                currentSparseBlockDataEnd = std::min(currentSparseBlockDataEnd, majorSegmentDataStartIdx + nodeSize_L1);
                                currentSparseBlockDataEnd = std::min(currentSparseBlockDataEnd, treeSize);
                                for (int i = currentSparseBlockDataStart; i < currentSparseBlockDataEnd; ++i) {
                                    int d = dataArr[1][i];
                                    cumulative_in_major_segment[d]++;
                                }
                            }
                        }
                    }
                } else {
                    // if the flag is disabled do a normal dense array for level 1
                    if (dataArr != nullptr && dataArr[1] != nullptr) { // only applicable if dataArr[1] exists
                        rankMatrixArr[1] = new uint8_t[treeSize];
                        uint8_t* row1_dense = static_cast<uint8_t*>(rankMatrixArr[1]);
                        int nodeSize_L1_dense = numElementsLvl(1);
                        int numMajorSegments_L1_dense = (nodeSize_L1_dense > 0) ? (treeSize / nodeSize_L1_dense) : 0;

                        for (int majorSegIdx = 0; majorSegIdx < numMajorSegments_L1_dense; ++majorSegIdx) {
                            int segment_start_idx = majorSegIdx * nodeSize_L1_dense;
                            std::vector<int> count_in_segment(arity, 0); // counts within this major segment
                            for (int j = 0; j < nodeSize_L1_dense; ++j) {
                                int current_global_idx = segment_start_idx + j;
                                if (current_global_idx < treeSize) { 
                                    int d_symbol = dataArr[1][current_global_idx];
                                    count_in_segment[d_symbol]++;
                                    int r = count_in_segment[d_symbol] - 1; // rank of d_symbol within this segment
                                    row1_dense[current_global_idx] = static_cast<uint8_t>(r);
                                }
                            }
                        }
                    }
                }
            }
        }

        if (paddedInput != nullptr) {
            delete[] paddedInput;
        }
    }

    ~FlatMultiaryWaveletTree() {
        if (dataArr != nullptr) {
            for (int i = 0; i < numLevels; i++) {
                delete[] dataArr[i];
            }
            delete[] dataArr;
        }

        if (rankMatrixArr != nullptr) {
            if (numLevels > 0 && rankMatrixArr[0] != nullptr) {
                delete[] static_cast<uint16_t*>(rankMatrixArr[0]);
            }
            if (numLevels > 1 && rankMatrixArr[1] != nullptr) {
                if (m_enableSparseL1Rank) { // destructor depends on whether sparsity was enabled
                    delete[] static_cast<uint16_t*>(rankMatrixArr[1]); 
                } else {
                    delete[] static_cast<uint8_t*>(rankMatrixArr[1]); // for dense array in level 1
                }
            }
            delete[] rankMatrixArr;
        }
    }


    int access(int i) {
        if (i >= originalSize || i < 0 || treeSize == 0) {
            return -1;
        }
        return accessFlat(i);
    }

    uint64_t size() {
        uint64_t total_size = 0;
        if (numLevels > 0 && treeSize > 0) { 
            total_size += static_cast<uint64_t>(numLevels) * treeSize * sizeof(uint8_t);
        }

        if (numLevels > 0 && rankMatrixArr != nullptr && rankMatrixArr[0] != nullptr) { 
            int numBuckets0 = (treeSize + L0_SPARSE_BUCKET_SIZE - 1) / L0_SPARSE_BUCKET_SIZE;
            total_size += static_cast<uint64_t>(numBuckets0) * arity * sizeof(uint16_t);
        }
        if (numLevels > 1 && rankMatrixArr != nullptr && rankMatrixArr[1] != nullptr) {
            if (m_enableSparseL1Rank) {
                // for sparse array in L1 calc
                int nodeSize_L1 = numElementsLvl(1);
                int numMajorSegments_L1 = (nodeSize_L1 > 0) ? (treeSize / nodeSize_L1) : 0;
                int numSparseBucketsPerMajorSegment = 0;
                if (nodeSize_L1 > 0) {
                    numSparseBucketsPerMajorSegment = (nodeSize_L1 + L1_SPARSE_BUCKET_SIZE - 1) / L1_SPARSE_BUCKET_SIZE;
                }
                total_size += static_cast<uint64_t>(numMajorSegments_L1) * numSparseBucketsPerMajorSegment * arity * sizeof(uint16_t);
            } else {
               //dense array L1 calc
                total_size += static_cast<uint64_t>(treeSize) * sizeof(uint8_t); 
            }
        }
        return total_size;
    }

private:
    int calculateLevel() { 
        if (treeSize <= 0) return 0; 
        if (treeSize == 1) return 1;
        int local_lvl = 0;
        long long n_val = treeSize; 
        while (n_val > 1) {
            n_val /= arity;
            local_lvl++;
        }
        return local_lvl + 1;
    }
    int numElementsLvl(int lvl_idx) { 
        if (treeSize == 0) return 0;
        long long divisor = 1;
        for (int i = 0; i < lvl_idx; ++i) {
            if (static_cast<double>(divisor) * arity > std::numeric_limits<long long>::max()) {
                return 0; // divisor overflow - 0 elements per node at this depth
            }
            divisor *= arity;
        }
        if (divisor == 0) return treeSize; // Should not happen if arity >= 1
        return treeSize / divisor;
    }

    void build_node(int level_idx, int node_idx_param, const int *data_param, int data_size_param, int lo, int hi) {
        if (data_size_param <= 0 || treeSize == 0) return;
        if (level_idx == numLevels) return; 
        if (dataArr == nullptr || dataArr[level_idx] == nullptr) return; // sanity check

        int current_node_elements_total = numElementsLvl(level_idx); 
        if (current_node_elements_total == 0 && treeSize > 0 && level_idx < numLevels) return; // if segment is 0

        int node_physical_start_idx = node_idx_param * current_node_elements_total; 

        int full_range_for_node = hi - lo + 1;
        if (full_range_for_node <= 0) full_range_for_node = 1; 
        
        int val_range_per_bucket = (full_range_for_node + arity - 1) / arity;
        if (val_range_per_bucket == 0) val_range_per_bucket = 1;

        vector<uint8_t> symbols_for_current_data(data_size_param);
        vector<vector<int>> child_actual_data(arity); 
        vector<int> child_value_counts(arity, 0);     

        for (int i = 0; i < data_size_param; ++i) {
            int current_value = data_param[i];
            int bucket_symbol = (current_value - lo) / val_range_per_bucket;
            
            if (bucket_symbol >= arity) bucket_symbol = arity - 1;
            if (bucket_symbol < 0) bucket_symbol = 0; 

            symbols_for_current_data[i] = static_cast<uint8_t>(bucket_symbol);
            child_actual_data[bucket_symbol].push_back(current_value);
            child_value_counts[bucket_symbol]++;
        }

        for (int i = 0; i < data_size_param; ++i) {
            if (node_physical_start_idx + i < treeSize) { 
                 dataArr[level_idx][node_physical_start_idx + i] = symbols_for_current_data[i];
            }
        }
        
        for (int b = 0; b < arity; ++b) {
            int count_for_this_child = child_value_counts[b];
            if (count_for_this_child > 0) {
                int child_conceptual_node_index = node_idx_param * arity + b;
                int child_lo = lo + b * val_range_per_bucket;
                int child_hi = lo + (b + 1) * val_range_per_bucket - 1;
                if (child_hi >= hi && b < arity -1) child_hi = std::max(lo, hi-1); // child_hi should not go below child_lo 
                if (b == arity -1) child_hi = hi; 
                if (child_lo > child_hi) child_lo = child_hi; //lo <= hi for child

                build_node(level_idx + 1, child_conceptual_node_index, child_actual_data[b].data(), count_for_this_child, child_lo, child_hi);
            }
        }
    }
    
    int accessFlat(int i_param) {
        int current_lvl = 0;
        int current_nodeIdx = 0;
        int current_pos = i_param;

        while (current_lvl < numLevels) {
            if (dataArr == nullptr || dataArr[current_lvl] == nullptr) return -1; // not needed here - caught by constructor

            int nodeSizeForLvl = numElementsLvl(current_lvl);
            if (nodeSizeForLvl == 0) return -1;

            int node_start = current_nodeIdx * nodeSizeForLvl;
            int globalIndex = node_start + current_pos;

            if (globalIndex < 0 || globalIndex >= treeSize) return -1;

            int bucket = dataArr[current_lvl][globalIndex];
            int newPos;

            if (current_lvl == 0) {
                if (rankMatrixArr != nullptr && rankMatrixArr[0] != nullptr) {
                    int bucketIdx = globalIndex / L0_SPARSE_BUCKET_SIZE;
                    int bucketStart = bucketIdx * L0_SPARSE_BUCKET_SIZE;
                    uint16_t* row0 = static_cast<uint16_t*>(rankMatrixArr[0]);
                    int d_sym = dataArr[current_lvl][globalIndex];
                    int baseCount = row0[static_cast<size_t>(bucketIdx) * arity + d_sym];
                    int localCount = 0;
                    for (int j = bucketStart; j < globalIndex; j++) {
                        if (dataArr[current_lvl][j] == d_sym) {
                            localCount++;
                        }
                    }
                    newPos = baseCount + localCount;
                } else { newPos = current_pos; } 
            } else if (current_lvl == 1) {
                if (m_enableSparseL1Rank) { // flag enabled for sparse L1
                    if (rankMatrixArr != nullptr && rankMatrixArr[1] != nullptr) {
                        // sparse L1 lookup
                        uint16_t* row1_sparse_ptr = static_cast<uint16_t*>(rankMatrixArr[1]);
                        int nodeSize_L1 = numElementsLvl(1);
                        int numSparseBucketsPerMajorSegment = 0;
                        if (nodeSize_L1 > 0) {
                            numSparseBucketsPerMajorSegment = (nodeSize_L1 + L1_SPARSE_BUCKET_SIZE - 1) / L1_SPARSE_BUCKET_SIZE;
                        }
                        int majorSegmentSelector = current_nodeIdx; 
                        int relativeIndexInSeg = current_pos;
                        int sparseBlockIdxInMajorSeg = relativeIndexInSeg / L1_SPARSE_BUCKET_SIZE;
                        
                        size_t sparseRankTableEntryBase = (static_cast<size_t>(majorSegmentSelector) * numSparseBucketsPerMajorSegment + sparseBlockIdxInMajorSeg) * arity;
                        int symbolAtIndex = dataArr[current_lvl][globalIndex];
                        
                        int baseCount = row1_sparse_ptr[sparseRankTableEntryBase + symbolAtIndex];
                        int localCount = 0;
                        int scanStartAbsolute = (majorSegmentSelector * nodeSize_L1) + (sparseBlockIdxInMajorSeg * L1_SPARSE_BUCKET_SIZE);
                        for (int j = scanStartAbsolute; j < globalIndex; ++j) {
                            if (dataArr[current_lvl][j] == symbolAtIndex) {
                                localCount++;
                            }
                        }
                        newPos = baseCount + localCount;
                    } else {
                        newPos = current_pos; // check -> fallback
                    }
                } else { // if sparse flag was disabled then do regular access
                    if (rankMatrixArr != nullptr && rankMatrixArr[1] != nullptr) {
                        // dense lookup at L1w
                        newPos = static_cast<int>(static_cast<uint8_t*>(rankMatrixArr[1])[globalIndex]);
                    } else {
                        newPos = current_pos; // check -> fallback
                    }
                }
            } else { 
                newPos = current_pos; // just in case there are more than 2 levels
            }
            current_pos = newPos;
            current_nodeIdx = current_nodeIdx * arity + bucket;
            current_lvl++;
        }
        return current_nodeIdx;
    }

}; 
#endif
