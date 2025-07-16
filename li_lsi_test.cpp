#include <chrono>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <vector>
#include <algorithm> 
#include <learned_secondary_index.hpp>

using namespace learned_secondary_index;
using namespace std;

bool file_exists_and_nonempty(const std::string &filename) {
    std::ifstream infile(filename, std::ios::binary);
    if (!infile.good()) return false;
    infile.seekg(0, std::ios::end);
    return infile.tellg() > 0;
}

std::vector<uint64_t> generateRandomLookups(const std::vector<uint64_t>& unsorted_keys, int num_lookups, int seed) {
    std::vector<uint64_t> lookups;
    std::mt19937 rng(seed);
    std::uniform_int_distribution<size_t> dist(0, unsorted_keys.size() - 1);
    for (int i = 0; i < num_lookups; ++i) {
        lookups.push_back(unsorted_keys[dist(rng)]);
    }
    return lookups;
}


double time_access(LearnedSecondaryIndex<uint64_t>& lsi, const vector<uint64_t>& unsorted_keys, const vector<uint64_t>& lookups, int num_trials) {
    using namespace std::chrono;
    long long total_ns = 0;
    int total_accesses = 0;

    for (int trial = 0; trial < num_trials; ++trial) {
        for (const auto& key : lookups) {
            auto start = high_resolution_clock::now();
            auto iter = lsi.lookup<false>(unsorted_keys.begin(), unsorted_keys.end(), key);
            if (iter != lsi.end()) {
                volatile auto result = *iter;  // prevent optimization
            }
            auto end = high_resolution_clock::now();
            total_ns += duration_cast<nanoseconds>(end - start).count();
            ++total_accesses;
        }
    }
    return static_cast<double>(total_ns) / total_accesses;
}


int main(int argc, char **argv) {
    if (argc < 6) {
        std::cout << "Usage: " << argv[0] << " <unsorted_data.bin> <num_lookups> <num_trials> <seed> <output_file>" << endl;
        return 1;
    }

    std::string unsorted_path = argv[1];
    int num_lookups = atoi(argv[2]);
    int num_trials = atoi(argv[3]);
    int seed = atoi(argv[4]);
    std::string output_file = argv[5];

    try {
        // Load unsorted data
        std::vector<uint64_t> unsorted_keys;
        {
            std::ifstream ifs(unsorted_path, ios::binary);
            if (!ifs) throw runtime_error("Failed to open unsorted data file: " + unsorted_path);
            uint32_t val32;
            while (ifs.read(reinterpret_cast<char*>(&val32), sizeof(uint32_t))) {
                unsorted_keys.push_back(static_cast<uint64_t>(val32));
            }
        }

        // Build LSI
        using LSIType = LearnedSecondaryIndex<uint64_t, learned_hashing::RadixSplineHash<uint64_t, 18, 16>, 0, false>;
        LSIType lsi(unsorted_keys.begin(), unsorted_keys.end());
        std::cout << "Loaded LearnedSecondaryIndex on " << unsorted_keys.size() << " keys." << endl;

        // Generate random lookup keys
        auto lookups = generateRandomLookups(unsorted_keys, num_lookups, seed);

        // Time access
        double avg_ns = time_access(lsi, unsorted_keys, lookups, num_trials);

        size_t total_bytes = lsi.byte_size();

        bool file_nonempty = file_exists_and_nonempty(output_file);
        ofstream outfile(output_file, ios::app);
        if (!outfile.is_open()) throw runtime_error("Could not open output file.");

        if (!file_nonempty) {
            outfile << "bin_path,num_lookups,num_trials,seed,access_avg_ns,total_bytes\n";
        }

        outfile << unsorted_path << "," 
                << num_lookups << ","
                << num_trials << ","
                << seed << ","
                << avg_ns << ","
                << total_bytes << "\n";

        outfile.close();

        cout << "Results appended to " << output_file << endl;
        cout << "Access average: " << avg_ns << " ns" << endl;
        cout << "Total index bytes: " << total_bytes << endl;

    } catch (exception &e) {
        std::cerr << "Error: " << e.what() << endl;
        return 1;
    }

    return 0;
}
