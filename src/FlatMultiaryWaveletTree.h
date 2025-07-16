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
    const int arity = 256;  // 256-way tree
    int originalSize;
    int treeSize;
    int numLevels;

    //we will put 4 8-bit data in the 32 bits
    uint32_t **dataArr;

    //if it is of type void we can static cast it to whatever type we want
    //alternatively we could have stored different arrays for every level
    //doc did not specify any performance drawback of this method
    void **rankMatrixArr;

    FlatMultiaryWaveletTree(const int *input, int size) {
        originalSize = size;
        int paddedSize = 1;
        while (paddedSize < size) {
            paddedSize *= arity;
        }
        treeSize = paddedSize;
        // padded array will now be next power of 256
        int *paddedInput = new int[treeSize];
        for (int i = 0; i < size; i++) {
            paddedInput[i] = input[i];
        }
        // pad by using the values after size - these will be in sorted order
        for (int i = size; i < treeSize; i++) {
            paddedInput[i] = i;
        }
        int totalLevels = calculateLevel();
        numLevels = totalLevels - 1;
        // bpnumLevels = 4;

        //for 256 way tree max value of elements is 255 so bits will always be
        //max 8 bits
        dataArr = new uint32_t*[numLevels];
        for (int lvl = 0; lvl < numLevels; lvl++) {
            int numWords = (treeSize + 3)/4;
            dataArr[lvl] = new uint32_t[numWords];
            for(int i=0;i<numWords;i++){
                dataArr[lvl][i] = 0;
            }
        }
    
        build_node(0, 0, paddedInput, treeSize, 0, originalSize - 1);
        int rankLvl = numLevels;
        if (numLevels>2){
            rankLvl = 2;
        }
        rankMatrixArr = new void*[rankLvl];
        //there is at least one level
        if (numLevels > 0){
            //considering max 16M elements, 1st row will have max 16M/256 = 65K
            //65K elements is 2^16 and therefore requires 16 bits of storage.
            int numWords0 = (treeSize+1)/2;
            rankMatrixArr[0] = new uint32_t[numWords0];
            uint32_t* row0 = static_cast<uint32_t*>(rankMatrixArr[0]);
            for(int i=0;i<numWords0;i++){
                row0[i]=0;
            }
            int nodeSize = numElementsLvl(0);
            //this means the number of elements in each block or bucket for that lvl
            //as per our logic for 16M elements, it will be 65k nodeSize at 0th level
            int numNodes = treeSize / nodeSize; 
            //this gives number of blocks/buckets
            for(int node = 0; node<numNodes; node++){
                int block_start = node * nodeSize;
                vector<int> count(arity, 0);
                for(int j=0; j<nodeSize; j++){
                    int idx = block_start+j;
                    int d = unpackData(dataArr[0],idx);
                    count[d]++;
                    uint16_t r = count[d]-1;
                    packRank0(row0, idx, r);
                }
            }
        }
        //if we have more than one level
        //literally the same code as level 0
        //todo maybe we can make this concise but not sure yet
        if (numLevels > 1){
            int numWords1 = (treeSize+3)/4;
            rankMatrixArr[1] = new uint32_t[numWords1];
            uint32_t* row1 = static_cast<uint32_t*>(rankMatrixArr[1]);
            for(int i=0;i<numWords1;i++){
                row1[i]=0;
            }
            int nodeSize = numElementsLvl(1);
            //here this value (nodeSize) should be 256
            int numNodes = treeSize/nodeSize;
            for (int node = 0; node<numNodes; node++){
                int block_start = node*nodeSize;
                vector<int> count(arity,0);
                for (int j=0; j<nodeSize; j++){
                    int idx = block_start+j;
                    int d=unpackData(dataArr[1],idx);
                    count[d]++;
                    uint8_t r = static_cast<uint8_t>(count[d]-1);
                    packRank1(row1,idx,r);
                }
            }
        }
        delete[] paddedInput;
    }

    ~FlatMultiaryWaveletTree() {
        for (int i = 0; i < numLevels; i++) {
            int numWords = (treeSize+3)/4;
            delete[] dataArr[i];
        }
        delete[] dataArr;
        if (numLevels>0){
            delete[] static_cast<uint32_t*>(rankMatrixArr[0]);
        }
        if (numLevels>1){
            delete[] static_cast<uint32_t*>(rankMatrixArr[1]);
        }
        delete[] rankMatrixArr;

    }

    int access(int i) {
        if (i >= originalSize || i < 0) {
            return -1;
        }
        return accessFlat(i);
    }

    uint64_t size() {
        uint64_t total_size = 0;
        for (int lvl = 0; lvl < numLevels; lvl++) {
            total_size += ((treeSize + 3) / 4) * sizeof(uint32_t);
        }
        //this is for the data
        total_size += numLevels * sizeof(uint32_t*);
        //now for rank level 0
        if (numLevels>0){
            total_size+=((treeSize+1)/2)*sizeof(uint32_t);
        }
        //this is for rank level 1
        if(numLevels>1){
            total_size+=((treeSize+3)/4)*sizeof(uint32_t);
        }
        total_size += numLevels * sizeof(uint32_t*);
        return total_size;
    }
   
   private:
     uint8_t unpackData(uint32_t* arr, uint64_t index) {
        uint64_t wordIdx = index / 4;
        int slot = index % 4;
        uint32_t word = arr[wordIdx];
        return static_cast<uint8_t>((word >> (8 * slot)) & 0xFF);
    }

    void packData(uint32_t* arr, uint64_t index, uint8_t value) {
        uint64_t wordIdx = index / 4;
        int slot = index % 4;
        uint32_t mask = ~(0xFFu << (8 * slot));
        arr[wordIdx] = (arr[wordIdx] & mask) | (static_cast<uint32_t>(value) << (8 * slot));
    }

    void packRank0(uint32_t* arr, uint64_t index, uint16_t rankVal) {
        uint64_t wordIdx = index / 2;
        int slot = index % 2;
        uint32_t mask = ~(0xFFFFu << (16 * slot));
        arr[wordIdx] = (arr[wordIdx] & mask) | (static_cast<uint32_t>(rankVal) << (16 * slot));
    }

    void packRank1(uint32_t* arr, uint64_t index, uint8_t rankVal) {
        uint64_t wordIdx = index / 4;
        int slot = index % 4;
        uint32_t mask = ~(0xFFu << (8 * slot));
        arr[wordIdx] = (arr[wordIdx] & mask) | (static_cast<uint32_t>(rankVal) << (8 * slot));
    }

    int calculateLevel() {
        int lvl = 0;
        int n = treeSize;
        while (n > 1) {
            n /= arity;
            lvl++;
        }
        return lvl + 1;
    }

    int numElementsLvl(int lvl) { 
        //calculates number of elements in each bucket/segment for that level
        // int bits_to_shift = log(arity)/log(2);                                                            
        return treeSize / (1 << (8 * lvl));
        // because 256^lvl = 2^(8*lvl)
    }

    // level: current level (0 <= level < numLevels)
    // node_index: index of the current node at that level
    // data: pointer to the node's data (length data_size)
    // data_size: number of elements in this node
    // Revised build_node using your bucket‐determination method for non‐base levels.
    //For both the base (leaf) and non-base levels, we determine buckets by first
    // finding the minimum value in the block, then "normalizing" each value, and then
    // partitioning the normalized range evenly among 'arity' buckets.
    // Revised build_node: Now uses the global range [lo, hi]
// Revised build_node (global-range version)
void build_node(int level, int node_index, const int *data, int data_size, int lo, int hi) {
    if (data_size <= 0)
        return;

    // Stop recursion if we have reached the desired number of levels.
    if (level == numLevels)
        return;  // Do not build further; leaves are implicit.

    int nodeSize = numElementsLvl(level);
    int node_start = node_index * nodeSize;

    int range = hi - lo + 1;
    int bucketSize = range / arity;  
    // We need this for unsorted data
    vector<int> bucketArr(data_size);
    vector<int> count(arity, 0);
    for (int i = 0; i < data_size; i++) {
        int key = data[i];
        int bucket = (key - lo) / bucketSize;
        if (bucket >= arity)
            bucket = arity - 1;
        bucketArr[i] = bucket;
        count[bucket]++;
    }

    // child arrays have same size
    vector<vector<int>> childData(arity);
    for (int b = 0; b < arity; b++) {
        childData[b].resize(count[b]);
    }
    vector<int> childCounters(arity, 0);
    for (int i = 0; i < data_size; i++) {
        int b = bucketArr[i];
        childData[b][childCounters[b]++] = data[i];
    }

    //dataArr has the bucket values
    for (int i = 0; i < data_size; i++) {
        uint8_t b = static_cast<uint8_t>(bucketArr[i]);
        packData(dataArr[level], node_start + i, b);
    }

    // rceursion for child node with global range
    for (int b = 0; b < arity; b++) {
        int childCount = count[b];
        if (childCount > 0) {
            int child_node_index = node_index * arity + b;
            int new_lo = lo + b * bucketSize;
            int new_hi = lo + (b + 1) * bucketSize - 1;
            build_node(level + 1, child_node_index, childData[b].data(), childCount, new_lo, new_hi);
        }
    }
}

    int accessFlat(int i) {
        int lvl = 0;
        int nodeIdx = 0;
        int pos = i;
        while (lvl < numLevels) {
            int nodeSize = numElementsLvl(lvl);
            int node_start = nodeIdx * nodeSize;
            int globalIndex = node_start + pos;
            if (lvl == 0) {
                uint32_t* r0 = static_cast<uint32_t*>(rankMatrixArr[0]);
                __builtin_prefetch(&r0[globalIndex / 2], 0, 3);
            } else if (lvl == 1) {
                uint32_t* r1 = static_cast<uint32_t*>(rankMatrixArr[1]);
                __builtin_prefetch(&r1[globalIndex / 4], 0, 3);
            }
            int bucket = unpackData(dataArr[lvl], globalIndex);
            int newPos;
            if (lvl == 0) {
                uint32_t* r0 = static_cast<uint32_t*>(rankMatrixArr[0]);
                uint64_t wordIdx = globalIndex / 2;
                int slot = globalIndex % 2;
                newPos = (r0[wordIdx] >> (16 * slot)) & 0xFFFF;
            } else if (lvl == 1) {
                uint32_t* r1 = static_cast<uint32_t*>(rankMatrixArr[1]);
                uint64_t wordIdx = globalIndex / 4;
                int slot = globalIndex % 4;
                newPos = (r1[wordIdx] >> (8 * slot)) & 0xFF;
            } else {
                newPos = pos;
            }
            pos = newPos;
            nodeIdx = nodeIdx * arity + bucket;
            lvl++;
        }
        return nodeIdx;
    }
};

#endif
