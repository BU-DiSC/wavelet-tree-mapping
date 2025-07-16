# wavelet-tree-mapping

**Exploring Wavelet Trees as Space-Efficient Physical-to-Sorted Mapping for Learned Indexes**  

## Overview

This project implements and evaluates **Integer Wavelet Trees (IWTs)** and higher-fanout **T-way IWTs** as **space-efficient mappings** between the sorted and physical order of data in **learned indexes**.

- Evaluate trade-off between compression and access latency.
- Integrate IWTs with learned indexes like **RadixSpline**.
- Compare against classical indexes like **B+-trees** and **QuIT**.

---

## Repository Structure

```
wavelet-tree-mapping/
│
├── src/                    # header files for all wavelet tree variants in the paper
├── external/               # 3rd-party dependencies
├── workloads/              # place .bin workload files in this repo
├── results/                # output and logs will go here
├── permute_test.cpp           # wavelet tree tests independently
├── li_wavelet_test.cpp        # Learned index + wavelet trees
├── li_lsi_test.cpp            # Learned index + LSI
├── CMakeLists.txt             # Main CMake file
├── external/CMakeLists.txt    # dependency CMake file
└── README.md                  # Project guide
```

---

## Dependencies

- **CMake >= 3.4**
- **C++20**
- Git

**Automatically fetched dependencies:**
- [spdlog](https://github.com/gabime/spdlog)
- [cxxopts](https://github.com/jarro2783/cxxopts)
- [Roaring Bitmap](https://github.com/RoaringBitmap/CRoaring)
- [RadixSpline](https://github.com/BU-DiSC/RadixSpline)
- [LearnedSecondaryIndex](https://github.com/BU-DiSC/LearnedSecondaryIndex)

---

## 🔨 Build Instructions

```bash
git clone --recurse-submodules https://github.com/BU-Di-SC/wavelet-tree-mapping.git
cd wavelet-tree-mapping

mkdir build
cd build
cmake ..

make <binary_name>
```

---

## 🗃️ Generate Workloads

Our tests use `.bin` format for workloads - it can be easily tweaked to accept csv as well.  
We use the [BoDS framework](https://github.com/BU-DiSC/bods)

1. To replicate, use steps from [BoDS](https://github.com/BU-DiSC/bods) to generate bin files with varying degrees of sortedness.
2. Place the `.bin` files in a `workloads/` folder at the project root.

Example:
```
workloads/
 ├── n16777216_k0l0.bin      # Fully sorted
 ├── n16777216_k25l25.bin    # Partially sorted
 ├── ...
```

---

Example binaries:
- `permute_flat_multiary` — different variants of multiary trees can be selected later
- `li_wavelet_test` — learned index with wavelet
- `li_lsi_test` — learned index with LSI

---

## B+-Tree Baselines

B+-tree and QuIT experiments, refer:  
[Quick Insertion Tree (QuIT)](https://github.com/BU-DiSC/quick-insertion-tree/tree/main)

---

