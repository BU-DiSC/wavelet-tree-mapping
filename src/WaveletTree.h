#include <spdlog/spdlog.h>

#include <algorithm>
#include <iostream>
#include <vector>

#include "permute_vector.h"
// #include "bitvector_container.h"

using namespace std;
using namespace permute;
template <typename VECTOR_TYPE>
class WaveletTree {
    /* data */
    VECTOR_TYPE* T;
    int sigma;

   public:
    WaveletTree(vector<int> S, int _sigma) : sigma(_sigma) {
        int n = depth();
        T = new VECTOR_TYPE[n];
        build(S, 0, (int)S.size());
        spdlog::info("build complete... now building rank-select support");

        build_rs();
        spdlog::info("completed building rank-select support");
    }

    void build(vector<int>& S, const int a, const int b, int level = 0) {
        if (b - a <= 1) return;

        // int mid = (a + b) / 2;
        int mid = a + (b - a) / 2;
        // int mid = l + (r - l) / 2;
        T[level].resize(b - a);
        for (int i = a; i < b; i++) {
            if (S[i] >= mid) T[level].set(i);
        }
        // T[level].build_rs();

        auto pivot = stable_partition(S.begin() + a, S.begin() + b,
                                      [mid](int c) { return c < mid; });
        int z = distance(S.begin(), pivot);

        build(S, a, z, level + 1);
        build(S, z, b, level + 1);
        // T[level].optimize();
    }
    void build_rs() {
        for (int i = 0; i < depth(); i++) {
            T[i].build_rs();
        }
    }

    /* access to the value of S[i] */
    int access(int i, int level = 0, int l = 0, int r = -1) {
        if (r < 0) r = sigma;
        if (r - l <= 1) return l;

        // int mid = (l + r) / 2;
        int mid = l + (r - l) / 2;
        int idx = l + i;
        if (!T[level].get(idx)) {
            return access(T[level].rank0cached(idx) - T[level].rank0cached(l),
                          level + 1, l, mid);
        } else {
            return access(T[level].rank1cached(idx) - T[level].rank1cached(l),
                          level + 1, mid, r);
        }
    }

    int depth() const { return ceil(log2(sigma)); }

    int height() {
        int d = depth();
        return d;
    }

    uint64_t size() {
        uint64_t s = 0;
        uint64_t s_1 = 0;
        for (int i = 0; i < depth(); i++) {
            uint64_t size_before_optimize = T[i].get_size();
            s_1 += size_before_optimize;
            // T[i].optimize();
            s += T[i].get_size();
            // std::cout << size_before_optimize << ", " << T[i].get_size()
            //           << std::endl;
        }
        // std::cout << "s_1 = " << s_1 << std::endl;
        return s + sizeof(sigma);
    }

    ~WaveletTree() { delete[] T; }
};
