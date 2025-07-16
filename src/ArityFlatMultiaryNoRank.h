#ifndef FLAT_MULTIARY_WAVELET_TREE_H
#define FLAT_MULTIARY_WAVELET_TREE_H

#include <cmath>
#include <cstdint>
#include <cstdlib>  
#include <iostream>
#include <limits>    
#include <vector>
#include <algorithm> // this to use min

class FlatMultiaryWaveletTree {
public:
    const int arity;
    const int dataBits;
    const uint32_t dataMask; 
    int treeSize;
    int numLevels;
    uint8_t** levelData;

    FlatMultiaryWaveletTree(const uint64_t* input, int size, int arity_)
        : arity(arity_),
          dataBits(static_cast<int>(std::ceil(std::log2(arity_)))),
          dataMask((1U << static_cast<int>(std::ceil(std::log2(arity_)))) - 1)
    {
        uint64_t tempTreeSize = 1;
        while (tempTreeSize < static_cast<uint64_t>(size)) {
            if (tempTreeSize > UINT64_MAX / static_cast<uint64_t>(arity)) {
                std::cerr << "Error: treeSize calculation would overflow uint64_t.\n";
                exit(1);
            }
            tempTreeSize *= arity;
        }
        if (tempTreeSize > static_cast<uint64_t>(std::numeric_limits<int>::max())) {
            std::cerr << "Error: Final treeSize exceeds INT_MAX. Consider changing treeSize to uint64_t.\n";
            exit(1);
        }
        treeSize = static_cast<int>(tempTreeSize);

        int* paddedInput = new int[treeSize];
        for (int i = 0; i < size; i++) paddedInput[i] = static_cast<int>(input[i]);
        for (int i = size; i < treeSize; i++) paddedInput[i] = 0;

        numLevels = calculateLevels(treeSize);

        levelData = new uint8_t*[numLevels]; // Changed to uint8_t**
        for (int lvl = 0; lvl < numLevels; lvl++) {
            levelData[lvl] = new uint8_t[treeSize](); // Initialize to zero
        }

        buildNode(0, 0, paddedInput, treeSize);
        delete[] paddedInput;
    }

    ~FlatMultiaryWaveletTree() {
        for (int i = 0; i < numLevels; i++) {
            delete[] levelData[i];
        }
        delete[] levelData;
    }

    int access(int i) const {
        int pos = i;
        int node = 0; 
        for (int lvl = 0; lvl < numLevels; ++lvl) {
            int nodeSize = numElementsAtLevel(lvl);
            if (nodeSize == 0) {
                 std::cerr << "Error! nodeSize = 0 for access at lvl=" << lvl << std::endl;
                 return -1;
            }

            // calculate global start index of the current node block for specific lvl
            uint64_t node_block_start_u64 = static_cast<uint64_t>(node) * nodeSize;
            if (node_block_start_u64 >= static_cast<uint64_t>(treeSize)) {
                std::cerr << "Error!! OOB\n";
                return -1;
            }
            int node_block_start = static_cast<int>(node_block_start_u64);

            // Get bucket ID for the element at `pos` within this node's block
            // levelData[lvl] stores bucket IDs
            uint64_t current_element_global_idx_u64 = node_block_start_u64 + pos;
            if (current_element_global_idx_u64 >= static_cast<uint64_t>(treeSize)) {
                 std::cerr << "Error!! OOB\n";
                 return -1;
            }
            uint8_t bucket = levelData[lvl][static_cast<int>(current_element_global_idx_u64)];

            // do linear Scan for Rank
            // Find rank = number of occurrences of 'bucket' in the prefix of the current node's block
            // up to and including the current position 'pos'
            int rank = 0;
            // scan = iterate from the start of the current node's block up to 'pos'
            // to count elements that have the same 'bucket' ID.
            for (int k = 0; k <= pos; ++k) { // k is local index within the node's block
                if (levelData[lvl][node_block_start + k] == bucket) {
                    rank++;
                }
            }
            // rank is 0-indexed
            // subtract 1 because rank count is 1-indexed
            int new_pos = rank - 1; // new position for the next level's access
            pos = new_pos; // Update position for next level
            node = node * arity + bucket; // Update node index to traverse down tree
        }
        // At the last level, node variable = accumulated nodeIdx * arity + bucket
        // becomes the actual value so return
        return node;
    }

    uint64_t size() const {
        return static_cast<uint64_t>(numLevels) * treeSize * sizeof(uint8_t); // Changed sizeof(uint32_t) to sizeof(uint8_t)
    }

private:
    static uint64_t powi(int base, int exp) {
        uint64_t result = 1;
        for (int i = 0; i < exp; ++i) {
            if (result > UINT64_MAX / static_cast<uint64_t>(base)) {
                return UINT64_MAX;
            }
            result *= base;
        }
        return result;
    }

    int calculateLevels(int currentTreeSize) const {
        if (currentTreeSize <= 1) return 0;
        int lvl = 0;
        uint64_t n = currentTreeSize;
        while (n > 1) {
            n /= arity;
            lvl++;
        }
        return lvl;
    }

    int numElementsAtLevel(int lvl) const {
        uint64_t divisor = powi(arity, lvl);
        if (divisor == 0) {
            std::cerr << "Error!\n";
            return 0;
        }
        return static_cast<int>(static_cast<uint64_t>(treeSize) / divisor);
    }

    void buildNode(int lvl, int nodeIdx, const int* data, int len) {
        if (lvl >= numLevels || len == 0) return;

        int nodeSize = numElementsAtLevel(lvl);
        uint64_t startIdx_u64 = static_cast<uint64_t>(nodeIdx) * nodeSize;
        if (startIdx_u64 >= static_cast<uint64_t>(treeSize)) {
            std::cerr << "Error!\n";
            return;
        }
        int startIdx = static_cast<int>(startIdx_u64);

        int minVal = data[0];
        for (int i = 1; i < len; i++) {
            minVal = std::min(minVal, data[i]);
        }

        std::vector<int> bucketArr(len);
        std::vector<int> counts(arity, 0);

        for (int i = 0; i < len; i++) {
            int norm = data[i] - minVal;
            uint64_t scaled_norm = static_cast<uint64_t>(norm) * arity;
            if (len == 0) { 
                 std::cerr << "Error!\n";
                 exit(1);
            }
            int bucket = std::min(static_cast<int>(scaled_norm / len), arity - 1);

            if (bucket < 0 || bucket >= arity) {
                std::cerr << "Error! Calculated bucket " << bucket << " out of range [0, " << arity-1 << "] at level " << lvl << ", index " << i << std::endl;
                exit(1);
            }

            bucketArr[i] = bucket;
            counts[bucket]++;
        }

        std::vector<std::vector<int>> children(arity);
        for (int b = 0; b < arity; b++) children[b].reserve(counts[b]);

        // rankCount to data to children
        std::vector<int> rankCount(arity, 0); // track ranks for child data distribution

        for (int i = 0; i < len; i++) {
            int b = bucketArr[i];
            if (b < 0 || b >= arity) {
                std::cerr << "Critical Error: Invalid bucket value in packing loop: " << b << std::endl;
                exit(1);
            }
            // used to find pos in children for next recursion level
            int r = rankCount[b]++; // rank is implicit in order of children[b]

            children[b].push_back(data[i]); // Add element to the correct child's vector
            //pack
            uint64_t global_data_idx_u64 = static_cast<uint64_t>(startIdx) + i;
            if (global_data_idx_u64 >= static_cast<uint64_t>(treeSize)) {
                std::cerr << "Error!" << global_data_idx_u64 << ", treeSize=" << treeSize << ") at level " << lvl << std::endl;
                exit(1);
            }
            levelData[lvl][global_data_idx_u64] = static_cast<uint8_t>(b); // Store only the bucket
        }

        for (int b = 0; b < arity; b++) {
            if (!children[b].empty()) {
                buildNode(lvl + 1, nodeIdx * arity + b, children[b].data(), children[b].size());
            }
        }
    }
};

#endif
