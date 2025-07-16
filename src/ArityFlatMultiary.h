#ifndef FLAT_MULTIARY_WAVELET_TREE_H
#define FLAT_MULTIARY_WAVELET_TREE_H

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <algorithm> // this to use min

class FlatMultiaryWaveletTree {
public:
    const int arity;
    const int dataBits;
    const int rankBits;
    const uint32_t dataMask;

    int treeSize;
    int numLevels;
    uint32_t** levelData;

    FlatMultiaryWaveletTree(const uint64_t* input, int size, int arity_)
        : arity(arity_),
          dataBits(static_cast<int>(std::ceil(std::log2(arity_)))),
          rankBits(32 - static_cast<int>(std::ceil(std::log2(arity_)))),
          dataMask((1U << static_cast<int>(std::ceil(std::log2(arity_)))) - 1)
    {
        uint64_t tempTreeSize = 1;
        while (tempTreeSize < static_cast<uint64_t>(size)) {
            // to prevent overflow of treesize
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

        levelData = new uint32_t*[numLevels];
        for (int lvl = 0; lvl < numLevels; lvl++) {
            levelData[lvl] = new uint32_t[treeSize](); // Initialize to zero
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
            // int might cause overflow here
            int blockSize = static_cast<int>(static_cast<uint64_t>(treeSize) / powi(arity, lvl));
            if (blockSize == 0) { // in case of any bad block size
                 std::cerr << "Error: blockSize is 0 in access, lvl=" << lvl << ", arity=" << arity << ", treeSize=" << treeSize << std::endl;
                 return -1;
            }

            uint64_t idx_u64 = static_cast<uint64_t>(node) * blockSize + pos;
            if (idx_u64 >= static_cast<uint64_t>(treeSize)) { // Check against allocated size
                std::cerr << "Error: Access index out of bounds (idx_u64=" << idx_u64 << ", treeSize=" << treeSize << ") at level " << lvl << std::endl;
                return -1;
            }
            int idx = static_cast<int>(idx_u64);
            uint32_t packed = levelData[lvl][idx];
            int bucket = packed & dataMask;
            int rank = packed >> dataBits; // get rank i.e. the new position for next level

            pos = rank; // Update position for the next level
            node = node * arity + bucket; // Update node index to traverse down tree
        }
        return node;
    }

    uint64_t size() const {
        return static_cast<uint64_t>(numLevels) * treeSize * sizeof(uint32_t);
    }

private:
    // 'result' must be uint64_t to prevent overflow while multiplication
    static uint64_t powi(int base, int exp) {
        uint64_t result = 1;
        for (int i = 0; i < exp; ++i) {
            if (result > UINT64_MAX / static_cast<uint64_t>(base)) {
                return UINT64_MAX; // safeguard for overflow - indicate overflow or throw
            }
            result *= base;
        }
        return result;
    }

    // calculateLevels = number of levels (height) of the tree
    // uint64_t to handle large n
    int calculateLevels(int currentTreeSize) const {
        if (currentTreeSize <= 1) return 0;
        int lvl = 0;
        uint64_t n = currentTreeSize; // Changed to uint64_t
        while (n > 1) {
            n /= arity;
            lvl++;
        }
        return lvl;
    }

    // numElementsAtLevel = size of segment at a given level
    int numElementsAtLevel(int lvl) const {
        // Cast treeSize to uint64_t before division to ensure correct arithmetic
        uint64_t divisor = powi(arity, lvl);
        if (divisor == 0) { 
            std::cerr << "Error in calculating number of elements at level!\n";
            return 0; // Prevent division by zero
        }
        return static_cast<int>(static_cast<uint64_t>(treeSize) / divisor);
    }

    void buildNode(int lvl, int nodeIdx, const int* data, int len) {
        if (lvl >= numLevels || len == 0) return;

        int nodeSize = numElementsAtLevel(lvl);
        uint64_t startIdx_u64 = static_cast<uint64_t>(nodeIdx) * nodeSize;
        // overflow handle- check bounds against treeSize before casting to int index
        if (startIdx_u64 >= static_cast<uint64_t>(treeSize)) {
            std::cerr << "Error! out of bounds in buildNode - startIdx_u64=" << startIdx_u64 << ", treeSize=" << treeSize << ") at level " << lvl << ", node " << nodeIdx << std::endl;
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
                 std::cerr << "Error! len = 0 for bucket.\n";
                 exit(1);
            }
            int bucket = std::min(static_cast<int>(scaled_norm / len), arity - 1);            
            if (bucket < 0 || bucket >= arity) {
                std::cerr << "Error! Buclet range exceeded!! bucket =" << bucket << " is out of expected range [0, " << arity-1 << "] at level " << lvl << ", index " << i << std::endl;
                exit(1);
            }

            bucketArr[i] = bucket;
            counts[bucket]++;
        }
        std::vector<std::vector<int>> children(arity);
        for (int b = 0; b < arity; b++) children[b].reserve(counts[b]); // Reserve memory for efficiency

        std::vector<int> rankCount(arity, 0); // track rank within each bucket for specific level

        for (int i = 0; i < len; i++) {
            int b = bucketArr[i]; // get current bucket for current element
            if (b < 0 || b >= arity) { 
                std::cerr << "Critical Error! " << b << std::endl;
                exit(1);
            }
            int r = rankCount[b]++; // get current rank and increment for next
            children[b].push_back(data[i]); // add this element to the correct child's vector

            // pack rank and bucket into the levelData array
            uint64_t global_data_idx_u64 = static_cast<uint64_t>(startIdx) + i;
            if (global_data_idx_u64 >= static_cast<uint64_t>(treeSize)) {
                std::cerr << "Error: Packed data index out of bounds (global_data_idx_u64=" << global_data_idx_u64 << ", treeSize=" << treeSize << ") at level " << lvl << std::endl;
                exit(1);
            }
            levelData[lvl][global_data_idx_u64] = (static_cast<uint32_t>(r) << dataBits) | static_cast<uint32_t>(b);
        }

        // recursive call to buildNode for each child that has data
        for (int b = 0; b < arity; b++) {
            if (!children[b].empty()) {
                buildNode(lvl + 1, nodeIdx * arity + b, children[b].data(), children[b].size());
            }
        }
    }
};

#endif
