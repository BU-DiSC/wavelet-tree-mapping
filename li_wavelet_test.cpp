#include <spdlog/spdlog.h>

#include <chrono>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <vector>
#include <iomanip> //CSV formatting
#include <learned_hashing.hpp>
// #include "src/ArityFlatMultiary.h"
#include "src/FlatMultiaryWaveletTree.h"
// #include "src/FlatMultiaryWaveletTree_RankX.h"
// #include "src/FlatMultiaryWaveletTree_HRLE.h"
#include "src/WaveletTree.h"
#include "src/permute_roaring.h" 
#include "src/utils.h"
#include "rs/builder.h"
#include "rs/radix_spline.h"


using namespace std;
const size_t max_error = 16;
using RSHash = learned_hashing::RadixSplineHash<uint64_t, 18,max_error>; //18 is fanout and 16 is max error bits

// template <typename F>
// unsigned long long time_function(F f) {
//     auto start = std::chrono::high_resolution_clock::now();
//     f();
//     auto stop = std::chrono::high_resolution_clock::now();
//     return std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count();
// }


template <typename T>
bool validate_access_function(T &wt,
                              const std::vector<uint64_t> &original_permutation) {
    spdlog::info("validating");

    for (int i = 0; i < original_permutation.size(); i++) {
        int result = access(wt, i);
        if (result != original_permutation[i]) {
            spdlog::error("FAIL! at index {}: Expected {}, Got {}", i,
                          original_permutation[i], result);
            return false;
        }
    }

    spdlog::info("Success!");
    return true;
}

template <typename T>
int access(WaveletTree<T> &wt, int x) {
    return wt.access(x);
}

int access(FlatMultiaryWaveletTree &wt, int x) { return wt.access(x); }

int access(vector<int> p, int x) { return p.at(x); }

int access(int *p, int x) { return p[x]; }

template <typename T>
int binary_search(int start, int end, int val, const vector<uint64_t> &input_stream, T &wt,
                  long long &duration, long &num_access_calls,
                  bool regenerate_random = false) {
    int midpt_pos_in_sorted = start + (end - start) / 2;

    // start timer for calculating access time
    auto start_time = std::chrono::high_resolution_clock::now();
    int midpt_pos_in_input = access(wt, midpt_pos_in_sorted);
    // stop timer
    auto stop_time = std::chrono::high_resolution_clock::now();
    auto d = std::chrono::duration_cast<std::chrono::nanoseconds>(stop_time -
                                                                  start_time)
                 .count();
    duration += d;
#ifdef PRINT_ACCESS_DURATION
    std::cout << "duration = " << duration << ", time = " << d << std::endl;
#endif
    num_access_calls++;

    uint64_t midpt_val = input_stream[midpt_pos_in_input];

    if (val < midpt_val) {
        return binary_search(start, midpt_pos_in_sorted, val, input_stream, wt,
                             duration, num_access_calls);
    } else if (val > midpt_val) {
        return binary_search(midpt_pos_in_sorted, end, val, input_stream, wt,
                             duration, num_access_calls);
    } else {
        return midpt_pos_in_input;
    }
}

int binary_search(int start, int end, int val, int *input_stream,
                  int input_stream_size, int *p,
                  bool regenerate_random = false) {
    int midpt_pos_in_sorted = start + (end - start) / 2;

    int midpt_pos_in_input = access(p, midpt_pos_in_sorted);

    int midpt_val = input_stream[midpt_pos_in_input];

    if (val < midpt_val) {
        return binary_search(start, midpt_pos_in_sorted, val, input_stream,
                             input_stream_size, p);
    } else if (val > midpt_val) {
        return binary_search(midpt_pos_in_sorted, end, val, input_stream,
                             input_stream_size, p);
    } else {
        return midpt_pos_in_input;
    }
}

int scan(int val, vector<int> input_stream, vector<int> p) {
    for (int i = 0; i < p.size(); i++) {
        int x = access(p, i);
        if (input_stream[x] == val) {
            return x;
        }
    }
    return -1;
}

std::vector<uint64_t> generateRandomNumbers(int n, bool shuffle = true,
                                       int seed = -1) {
    std::vector<uint64_t> numbers;
    for (uint64_t i = 0; i < n; i++) {
        numbers.push_back(i);
    }

    if (shuffle) {
        spdlog::info("shuffling input");
        if (seed >= 0) {
            std::mt19937 g(seed);
            std::shuffle(numbers.begin(), numbers.end(), g);
        } else {
            std::random_device rd;
            std::mt19937 g(rd());
            std::shuffle(numbers.begin(), numbers.end(), g);
        }
    }

    return numbers;
}

std::vector<uint64_t> generateRandomLookups(int n, int num_lookups, int seed) {
    std::vector<uint64_t> lookups;
    std::mt19937 rng(seed);  // Random number generator with seed
    std::uniform_int_distribution<uint64_t> dist(0,
                                            n - 1);  // Distribution from 1 to n

    std::cout << "Random numbers: ";
    for (int i = 0; i < num_lookups; ++i) {
        lookups.push_back(dist(rng));
    }
    return lookups;
}

template <typename T, typename Index>
std::pair<unsigned long long, unsigned long long>
do_experiment(WaveletTree<T> &wt,
              const std::vector<uint64_t> &input_stream,
              const std::vector<uint64_t> &lookups,
              Index &rsh,
              const std::string &results_file,
              int n,
              int num_trials,
              long long &duration_wt,
              long &num_access_calls) {
    spdlog::info("Executing compare_binsearch_wt on with n = {}, num_lookups = {}, num_trials = {}",
                 n, lookups.size(), num_trials);

    unsigned long long time_taken_wt = 0;
    unsigned long long total_access_time = 0;
    unsigned long long total_search_bound_time = 0;
    duration_wt = 0;
    num_access_calls = 0;

    for (int t = 0; t < num_trials; ++t) {
        unsigned long long time_taken_for_trial = 0;
        for (auto q : lookups) {
            // 1) search bound
            // total_search_ns += time_function([&](){
            auto search_bound_start = std::chrono::high_resolution_clock::now();
            auto pred = rsh(q);
            auto search_bound_stop = std::chrono::high_resolution_clock::now();
            auto search_bound_duration = std::chrono::duration_cast<std::chrono::nanoseconds>(
                search_bound_stop - search_bound_start).count();
            total_search_bound_time += search_bound_duration;
            auto begin = std::max(static_cast<long>(pred - max_error), 0L);
            auto end = std::min(static_cast<long>(pred + max_error + 1), static_cast<long>(input_stream.size()));
            // 2) full binary search + wavelet access timing
            auto start = std::chrono::high_resolution_clock::now();
            time_function([&]() {
                int x = binary_search(begin, end, q,
                                input_stream, wt,
                                duration_wt, num_access_calls);
                x+=1;//dummy
            });
            auto stop = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count();
            time_taken_for_trial += duration;
        }
        total_access_time += (time_taken_for_trial / lookups.size());
    }
    total_access_time /= static_cast<double>(num_trials);
    double avg_wt_access_time = static_cast<double>(duration_wt) / num_access_calls;
    double avg_total_access_time = static_cast<double>(total_access_time);
    return {static_cast<unsigned long long>(avg_wt_access_time),
        static_cast<unsigned long long>(avg_total_access_time)};
    spdlog::info("Avg. Time taken per wavelet access call: {} nanoseconds", avg_wt_access_time);
    spdlog::info("Avg. Total time per access (wavelet + search bound): {} nanoseconds", avg_total_access_time);
    return {avg_wt_access_time, avg_total_access_time};
}

template<typename Index>
std::pair<unsigned long long, unsigned long long> do_experiment(
    FlatMultiaryWaveletTree &wt,
    const std::vector<uint64_t> &input_stream,
    const std::vector<uint64_t> &lookups,
    Index &rsh,
    const std::string &results_file,
    int n,
    int num_trials,
    long long &duration_wt,
    long &num_access_calls) {
    spdlog::info(
        "Executing compare_binsearch_wt on with n = {}, num_lookups = {}, "
        "and num_trials = {}",
        n, lookups.size(), num_trials);

    unsigned long long time_taken_wt = 0; // Time for wavelet tree access
    unsigned long long total_access_time = 0; // Total time (wavelet access + search bound)
    unsigned long long total_search_bound_time = 0;
    duration_wt = 0;
    num_access_calls = 0;

    for (int i = 0; i < num_trials; i++) {
        unsigned long long time_taken_for_trial = 0;
        for (int j = 0; j < lookups.size(); j++) {
            auto search_bound_start = std::chrono::high_resolution_clock::now();
                auto pred = rsh(lookups[j]);
            auto search_bound_stop = std::chrono::high_resolution_clock::now();
            auto search_bound_duration = std::chrono::duration_cast<std::chrono::nanoseconds>(
                search_bound_stop - search_bound_start).count();
            total_search_bound_time += search_bound_duration;

            auto begin = std::max(static_cast<long>(pred - max_error), 0L);
            auto end = std::min(static_cast<long>(pred + max_error + 1), static_cast<long>(input_stream.size()));   
            
            auto start = std::chrono::high_resolution_clock::now();
            time_function([&]() {
                int x = binary_search(begin, end, lookups[j],
                                        input_stream,
                                        wt,
                                        duration_wt,
                                        num_access_calls);
                x+=1;
            });
            auto stop = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count();
            time_taken_for_trial += duration;
        }
        total_access_time += (time_taken_for_trial / lookups.size());
    }
    total_access_time /= static_cast<double>(num_trials);
    double avg_wt_access_time = static_cast<double>(duration_wt) / num_access_calls;
    double avg_total_access_time = static_cast<double>(total_access_time);
    return {static_cast<unsigned long long>(avg_wt_access_time),
        static_cast<unsigned long long>(avg_total_access_time)};
    spdlog::info("Avg. Time taken per wavelet access call: {} nanoseconds", avg_wt_access_time);
    spdlog::info("Avg. Total time per access (wavelet + search bound): {} nanoseconds", avg_total_access_time);
    return {avg_wt_access_time, avg_total_access_time};
}

std::vector<uint64_t> read_numbers_from_file(const std::string& filename) {
    std::vector<uint64_t> numbers;
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        exit(1); 
    }

    std::string line;
    uint64_t value;
    while (std::getline(file, line)) {
        try {
            value = std::stoull(line);
            numbers.push_back(value);
        } catch (const std::exception& e) {
            std::cerr << "Could not parse line '" << line << "' from file " << filename << "Error = " << e.what() << std::endl;
        }
    }
    file.close();
    spdlog::info("Done reading from file{}!", filename);
    return numbers;
}

std::vector<uint64_t> read_numbers_from_bin_file(const std::string &filename) {
    std::vector<uint64_t> numbers;
    std::ifstream file(filename, std::ios::binary | std::ios::in);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open binary file " << filename << std::endl;
        exit(1);
    }

    // Read the file as `int` values (4 bytes per element)
    int value;
    while (file.read(reinterpret_cast<char *>(&value), sizeof(int))) {
        // Convert `int` to `uint64_t` and store in the vector
        numbers.push_back(static_cast<uint64_t>(value));
    }
    file.close();

    spdlog::info("Done reading from binary file: {}!", filename);
    return numbers;
}

int main(int argc, char **argv) {
    int n_actual;
    int n = atoi(argv[1]);
    cout << n;
    int num_trials = atoi(argv[2]);
    string results_file = argv[3];
    int num_lookups = atoi(argv[4]);
    int seed = atoi(argv[5]);
    bool random_input = atoi(argv[6]);
    int arity = atoi(argv[7]);
    string input_file_path = "";
    if (argc > 8) {
        input_file_path = argv[8];
    }

    spdlog::info("Writing results to file: {}", results_file);
    try {
        vector<uint64_t> input_stream;
        if (!input_file_path.empty()) {
            spdlog::info("Reading input from file: {}", input_file_path);
            input_stream = read_numbers_from_bin_file(input_file_path);
            n_actual = input_stream.size();
            if (n_actual == 0) {
                std::cerr << "Input file is empty or could not be read\n";
                return 1;
            }
        } else {
            spdlog::info("Generating random input numbers with size: {}", n);
            input_stream = generateRandomNumbers(n, random_input, seed);
            n_actual = input_stream.size();
        }
        vector<uint64_t> p = sorted_indices(input_stream);
        vector<int> p_int(p.begin(), p.end());
        vector<uint64_t> lookups = generateRandomLookups(n_actual, num_lookups, seed);
        auto start_construction = std::chrono::high_resolution_clock::now();
        // FlatMultiaryWaveletTree wt3(p.data(), n, arity);
        // FlatMultiaryWaveletTree wt3(p_int.data(), n);
        WaveletTree<permute::PermuteRoaring> wt3(p_int, n);
        auto end_construction = std::chrono::high_resolution_clock::now();
        auto construction_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end_construction - start_construction).count();
        spdlog::info("Wavelet‐tree build time = {} ns", construction_time);

        std::vector<uint64_t> sorted_p = p;
        std::sort(sorted_p.begin(), sorted_p.end());
        // rs::Builder<uint64_t> rsb(sorted_p.front(), sorted_p.back());
        RSHash rsh;
        rsh.train(sorted_p.begin(), sorted_p.end(), sorted_p.size());
        // for (const auto &val : sorted_p) rsb.AddKey(val);
        // rs::RadixSpline<uint64_t> rs = rsb.Finalize();
        
        long long duration_wt = 0;
        long num_access_calls = 0;
#ifdef PROFILE
        cout << "enter something to continue" << endl;
        cin >> c;
#endif
        auto [time_wavelet_access, time_total_access] = do_experiment(wt3, input_stream, lookups, rsh,
                                  results_file, n, num_trials,
                                  duration_wt, num_access_calls);
        uint64_t size = wt3.size();
#ifdef VALIDATE
        // if (!validate_access_function(wt3, p)) {
        //     std::cerr
        //         << "Validation multiary failed! Check the implementation of "
        //            "access or the wavelet tree."
        //         << std::endl;
        //     return 1;
        // }
#endif
        std::ofstream csv_file(results_file, std::ios::app);
        if (!csv_file.is_open()) {
            std::cerr << "Error: Could not open run.csv for writing.\n";
            return 1;
        }
        std::ifstream check_file(results_file);
        check_file.seekg(0, std::ios::end); // go to end of file
        if (check_file.tellg() == 0) { // If the file is empty
            csv_file << "filename,num_elements,arity,num_lookups,num_trials,size_bytes,time_wavelet_access_ns,time_total_access_ns,k,l\n";
        }

        if (check_file.tellg() == 0) { // If the file is empty
            csv_file << "filename,num_elements,arity,num_lookups,num_trials,size_bytes,time_wavelet_access_ns,time_total_access_ns,k,l\n";
        }
        check_file.close();
        //getting the k and l values from the input file name
        std::string filename = input_file_path.substr(input_file_path.find_last_of("/") + 1); 
        size_t k_pos = filename.find("k");
        size_t l_pos = filename.find("l");
        std::string k_value = "NA"; // Default value if k is not found
        std::string l_value = "NA"; // Default value if l is not found
        if (k_pos != std::string::npos && l_pos != std::string::npos) {
            k_value = filename.substr(k_pos + 1, l_pos - k_pos - 1); // Extract value between 'k' and 'l'
            size_t dot_pos = filename.find(".", l_pos);
            if (dot_pos != std::string::npos) {
                l_value = filename.substr(l_pos + 1, dot_pos - l_pos - 1); // Extract value after 'l' till '.'
            }
        }
        csv_file << std::fixed << std::setprecision(2); // Format floating-point numbers
        csv_file << input_file_path << ","              // File name
                 << input_stream.size() << ","          // Number of elements
                 << arity << ","                        // Arity
                 << num_lookups << ","                  // Number of lookups
                 << num_trials << ","                   // Number of trials
                 << size << ","                         // Size of wavelet tree (bytes)
                 << time_wavelet_access << ","          // Time for wavelet tree access (ns)
                 << time_total_access << ","            // Total time per access (ns)
                 << construction_time   << ","          // Construction time (ns)    
                 << k_value << ","                      // Value of k
                 << l_value << "\n"; 
        csv_file.close();
        spdlog::info("Results written to csv file");
        } catch (exception &e) {
        std::cout << "Exception caught: ";
        cout << e.what() << endl;
        return 1;
    }

    return 0;
};
