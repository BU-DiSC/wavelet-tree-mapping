#ifndef FLAT_MULTIARY_WAVELET_TREE_H
#define FLAT_MULTIARY_WAVELET_TREE_H

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <vector>

using namespace std;

class FlatMultiaryWaveletTree {
public:
    const int arity = 256;
    int originalSize;
    int treeSize;
    int numLevels;

    enum CompressionMethod { RAW = 0, RLE = 1 };

    struct Run {
        uint8_t symbol;
        uint16_t length;
    };

    struct BlockData {
        CompressionMethod method;
        vector<uint8_t> rawData;
        vector<Run> rleRuns;
    };

    struct HybridLevelData {
        vector<BlockData> blocks;
    };

    vector<HybridLevelData> hybridLevels;

    uint16_t* rank16Level;  // level 0 rank array
    uint8_t* rank8Level;    // level 1 rank array

    FlatMultiaryWaveletTree(const int* input, int size) {
        originalSize = size;
        int paddedSize = 1;
        while (paddedSize < size) paddedSize *= arity;
        treeSize = paddedSize;

        // Allocate flat rank arrays
        rank16Level = new uint16_t[treeSize]();
        rank8Level = new uint8_t[treeSize]();

        int* paddedInput = new int[treeSize];
        for (int i = 0; i < size; i++) paddedInput[i] = input[i];
        for (int i = size; i < treeSize; i++) paddedInput[i] = i;

        numLevels = calculateLevel() - 1;
        hybridLevels.resize(numLevels);

        build_node(0, 0, paddedInput, treeSize, 0, originalSize - 1);
        delete[] paddedInput;
    }

    ~FlatMultiaryWaveletTree() {
        delete[] rank16Level;
        delete[] rank8Level;
    }

    int access(int i) {
        if (i >= originalSize || i < 0) return -1;
        return accessFlat(i);
    }

    uint64_t size() const {
        uint64_t total = 0;
        for (const auto& level : hybridLevels) {
            for (const auto& block : level.blocks) {
                if (block.method == RAW)
                    total += block.rawData.size();
                else
                    total += block.rleRuns.size() * sizeof(Run);
            }
        }
        total += treeSize * sizeof(uint16_t);  // rank16
        total += treeSize * sizeof(uint8_t);   // rank8
        return total;
    }

private:
    static const size_t BLOCK_SIZE = 128;

    int calculateLevel() {
        int lvl = 0, n = treeSize;
        while (n > 1) {
            n /= arity;
            lvl++;
        }
        return lvl + 1;
    }

    int numElementsLvl(int lvl) {
        return treeSize / (1 << (8 * lvl));
    }

    void build_node(int level, int node_index, const int* data, int data_size, int lo, int hi) {
        if (data_size <= 0 || level == numLevels) return;

        int nodeSize = numElementsLvl(level);
        int nodeStart = node_index * nodeSize;

        int range = hi - lo + 1;
        int bucketSize = max(1, range / arity);

        vector<int> bucketArr(data_size), count(arity, 0);
        for (int i = 0; i < data_size; i++) {
            int bucket = (data[i] - lo) / bucketSize;
            if (bucket >= arity) bucket = arity - 1;
            bucketArr[i] = bucket;
            count[bucket]++;
        }

        vector<vector<int>> childData(arity);
        for (int b = 0; b < arity; b++) childData[b].resize(count[b]);
        vector<int> childCounters(arity, 0);
        vector<int> globalCount(arity, 0);

        for (int i = 0; i < data_size; i++) {
            int b = bucketArr[i];
            int globalIdx = nodeStart + i;
            int rankVal = globalCount[b]++;

            if (level == 0) rank16Level[globalIdx] = rankVal;
            if (level == 1) rank8Level[globalIdx] = static_cast<uint8_t>(rankVal);

            childData[b][childCounters[b]++] = data[i];
        }

        HybridLevelData& levelData = hybridLevels[level];
        for (int i = 0; i < data_size; i += BLOCK_SIZE) {
            size_t end = min(static_cast<size_t>(data_size), i + BLOCK_SIZE);
            vector<uint8_t> block(end - i);
            for (size_t j = i; j < end; ++j)
                block[j - i] = static_cast<uint8_t>(bucketArr[j]);

            vector<Run> runs;
            for (size_t j = 0; j < block.size();) {
                size_t len = 1;
                while (j + len < block.size() && block[j + len] == block[j]) ++len;
                runs.push_back({block[j], static_cast<uint16_t>(len)});
                j += len;
            }

            BlockData blk;
            if (runs.size() * sizeof(Run) < block.size()) {
                blk.method = RLE;
                blk.rleRuns = move(runs);
            } else {
                blk.method = RAW;
                blk.rawData = move(block);
            }
            levelData.blocks.push_back(move(blk));
        }

        for (int b = 0; b < arity; b++) {
            if (!childData[b].empty()) {
                int childIdx = node_index * arity + b;
                int new_lo = lo + b * bucketSize;
                int new_hi = lo + (b + 1) * bucketSize - 1;
                build_node(level + 1, childIdx, childData[b].data(), childData[b].size(), new_lo, new_hi);
            }
        }
    }

    int accessFlat(int i) {
        int lvl = 0, nodeIdx = 0, pos = i;

        while (lvl < numLevels) {
            int nodeSize = numElementsLvl(lvl);
            int globalIdx = nodeIdx * nodeSize + pos;

            int blockIndex = globalIdx / BLOCK_SIZE;
            int offset = globalIdx % BLOCK_SIZE;

            const BlockData& blk = hybridLevels[lvl].blocks[blockIndex];
            int bucket;
            if (blk.method == RAW) {
                bucket = blk.rawData[offset];
            } else {
                size_t runOffset = 0;
                for (const Run& run : blk.rleRuns) {
                    if (runOffset + run.length > offset) {
                        bucket = run.symbol;
                        break;
                    }
                    runOffset += run.length;
                }
            }

            if (lvl == 0) {
                pos = rank16Level[globalIdx];
            } else if (lvl == 1) {
                pos = rank8Level[globalIdx];
            }

            nodeIdx = nodeIdx * arity + bucket;
            lvl++;
        }
        return nodeIdx;
    }
};

#endif
