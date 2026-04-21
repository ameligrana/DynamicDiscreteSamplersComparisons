

#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <vector>
#include <chrono>
#include <cmath>
#include <fstream>
#include <string>
#include <algorithm>
#include <numeric>
#include "XoshiroCpp.hpp"
#include "QuickBucket.hpp"

using namespace std;

double median(std::vector<double> v) {
    if (v.empty()) throw std::domain_error("median of empty vector");
    std::sort(v.begin(), v.end());
    size_t n = v.size();
    if (n % 2 == 1) {
        // odd number of elements
        return v[n/2];
    } else {
        // even number of elements: average the two middle
        return 0.5 * (v[n/2 - 1] + v[n/2]);
    }
}

void append_column_in_place(const std::string& file_path, const std::vector<double>& vec) {
    std::ifstream fin(file_path);
    if (!fin.is_open()) {
        throw std::runtime_error("Could not open input file: " + file_path);
    }

    std::vector<std::string> lines;
    std::string line;

    // Read and modify header
    if (std::getline(fin, line)) {
        lines.push_back(line + ",BUS");
    } else {
        throw std::runtime_error("Empty input file: " + file_path);
    }

    // Read each row and append corresponding vector value
    size_t i = 0;
    while (std::getline(fin, line)) {
        if (i < vec.size()) {
            lines.push_back(line + ',' + std::to_string(vec[i]));
        } else {
            lines.push_back(line + ",");
        }
        ++i;
    }

    fin.close();

    if (i < vec.size()) {
        std::cerr << "Warning: vector has " << (vec.size() - i) << " extra elements\n";
    }

    // Now write everything back to the same file
    std::ofstream fout(file_path);
    if (!fout.is_open()) {
        throw std::runtime_error("Could not open output file: " + file_path);
    }

    for (const auto& l : lines) {
        fout << l << '\n';
    }
}

// Function to set up the BucketMethod with random weights
BucketMethod setup_sampler(int n, XoshiroCpp::Xoroshiro128Plus& rng) {
    vector<Element> elements;
    elements.reserve(n);

    std::normal_distribution<double> weight_dist(0.0, 1.0);

    for (int i = 0; i < n; ++i) {
        Element e(i, i, std::abs(weight_dist(rng)));
        elements.emplace_back(e);
    }

    std::uniform_int_distribution<size_t> seed_choice(0, 100000);

    return BucketMethod(seed_choice(rng), n, elements);
}

vector<int> benchmark_sample_static(XoshiroCpp::Xoroshiro128Plus& rng, BucketMethod& bucket, int n) {
    vector<int> samples;
    samples.reserve(n);

    for (int i = 0; i < n; ++i) {
        samples.push_back(bucket.random_sample_value());
    }

    return samples;
}

vector<int> benchmark_sample_dynamic_fixed(XoshiroCpp::Xoroshiro128Plus& rng, BucketMethod& bucket, int n) {
    vector<int> samples;
    samples.reserve(n);

    std::normal_distribution<double> weight_dist(0.0, 1.0);
    uniform_int_distribution<int> int_dist(0, n - 1);

    for (int i = 0; i < n; ++i) {
    	samples.push_back(bucket.random_sample_value());
        int j = int_dist(rng);
        Element e(j, j, std::abs(weight_dist(rng)));
        bucket.update_weight(j, e);
    }
    
    return samples;
}

std::vector<int> benchmark_sample_dynamic_variable(
    XoshiroCpp::Xoroshiro128Plus& rng,
    BucketMethod& bucket,
    int n,
    const std::vector<int>& insertion_order
) {

    std::vector<int> samples;
   samples.reserve(9*n);

    std::normal_distribution<double> weight_dist(0.0, 1.0);

    for (size_t i = 0; i < 9*n; ++i) {
        samples.push_back(bucket.random_sample_value());
        int index = insertion_order[i];
        double new_weight = std::abs(weight_dist(rng));
        Element e(index, index, new_weight);
        bucket.insert_element(e);
    }

    return samples;
}

std::vector<int> benchmark_sample_dynamic_decreasing(
    XoshiroCpp::Xoroshiro128Plus& rng,
    BucketMethod& bucket,
    const std::vector<int>& removal_order
) {
    int n = static_cast<int>(removal_order.size());
    int steps = n - n / 10;
    std::vector<int> samples;
    samples.reserve(steps);

    for (int i = 0; i < steps; ++i) {
        samples.push_back(bucket.random_sample_value());
        bucket.delete_element(removal_order[i]);
    }

    return samples;
}

double numerical_check() {
    std::vector<Element> elements;
    elements.reserve(3);
    elements.emplace_back(0, 0, 0.1);
    elements.emplace_back(1, 1, 0.9);
    elements.emplace_back(2, 2, 9.0*1.0e15);
    std::random_device rd;
     std::mt19937 gen(rd());
    std::uniform_int_distribution<> distrib(1, 100);
    BucketMethod bucket(distrib(gen), 3, elements);
    bucket.update_weight(2, Element(2, 2, 0.0));
    int count = 0;
    for (int i = 0; i < 1000000; ++i) {
        if (bucket.random_sample_value() == 0) {
            count += 1;
        }
    }
    return double(count) / 1000000.0;
}

#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include "XoshiroCpp.hpp"
#include "QuickBucket.hpp"

using namespace std;

// --------------------------------------------------------------------------------
// decaying_weights_sampling (normalized), using BucketMethod interface:
//
//   • Build k[i] = (2.0 + (1/(100*n))*(i+1))^1000 for i = 0..n-1
//   • Construct a BucketMethod with those initial weights (Element(i, i, k[i])).
//   • Perform t “decay” rounds: for each i=0..n-1, divide k[i] by
//       (2.0 + (1/(100*n))*(i+1))
//     and call bucket.update_weight(i, Element(i, i, k[i])).
//   • Draw 1,000,000 samples via bucket.random_sample_value(), tally raw counts.
//   • Normalize counts → probabilities (each count / total samples).
//   • Return a std::vector<double> of length n with those normalized frequencies.
// --------------------------------------------------------------------------------
std::vector<double> decaying_weights_sampling(int n, int t) {
    // 1) Build initial weight vector k[0..n-1].
    std::vector<double> k(n);
    for (int i = 0; i < n; ++i) {
        double base = 2.0 + (1.0 / (100.0 * double(n))) * double(i + 1);
        k[i] = std::pow(base, 1000.0);
    }

    // 2) Set up Xoroshiro128Plus RNG and BucketMethod
    XoshiroCpp::Xoroshiro128Plus rng(std::chrono::high_resolution_clock::now()
                                         .time_since_epoch()
                                         .count());
    // Build vector<Element> for BucketMethod constructor
    std::vector<Element> elements;
    elements.reserve(n);
    for (int i = 0; i < n; ++i) {
        // Element(index, value, weight)
        elements.emplace_back(i, i, k[i]);
    }
    // Choose a random seed for BucketMethod
    std::uniform_int_distribution<uint64_t> seed_dist(0, UINT64_MAX);
    uint64_t bucket_seed = seed_dist(rng);
    BucketMethod bucket(bucket_seed, n, elements);

    // 3) Perform t “decay” rounds
    for (int round = 0; round < t; ++round) {
        for (int i = 0; i < n; ++i) {
            double divisor = 2.0 + (1.0 / (100.0 * double(n))) * double(i + 1);
            k[i] /= divisor;
            // Update the i-th element’s weight
            bucket.update_weight(i, Element(i, i, k[i]));
        }
    }

    // 4) Draw 1,000,000 samples and tally raw counts
    constexpr int N_SAMPLES = 1000000;
    std::vector<std::size_t> counts(n, 0);
    for (int rep = 0; rep < N_SAMPLES; ++rep) {
        int sampled_value = bucket.random_sample_value();
        if (sampled_value >= 0 && sampled_value < n) {
            ++counts[sampled_value];
        }
    }

    std::vector<double> normalized(n, 0.0);
    for (int i = 0; i < n; ++i) {
        normalized[i] = double(counts[i]) / N_SAMPLES;
    }
    return normalized;
}

int main() {

    const int n = 100;
    const int t_max = 50; // after that it stops to work

    // Open output CSV. Each row: 1000 normalized values for a given t.
    std::ofstream fout("../../../data/bus_numerical.csv");
    if (!fout.is_open()) {
        std::cerr << "Error: could not open the file for writing.\n";
        return 1;
    }
    fout << std::fixed << std::setprecision(8);

    for (int t = 1; t <= t_max; ++t) {
        auto normalized = decaying_weights_sampling(n, t);

        // Write one CSV row
        for (int i = 0; i < n; ++i) {
            fout << normalized[i];
            if (i + 1 < n) fout << ',';
        }
        fout << '\n';
        fout.flush();  // ensure each row is pushed to disk immediately

        std::cerr << "Completed t = " << t << " / " << t_max << "\n";
    }

    fout.close();

    constexpr uint64_t seed = 42;
    XoshiroCpp::Xoroshiro128Plus rng(seed);

    int repetitions = 50;

    std::vector<double> static_times;
    std::vector<double> dynamic_fixed_times;
    std::vector<double> dynamic_variable_times;
    std::vector<double> dynamic_decreasing_times;

    for (int exp = 3; exp <= 7; ++exp) {
        size_t size = static_cast<size_t>(std::pow(10, exp));
        int steps = static_cast<int>(size - size / 10);

        std::vector<double> static_ns;
        std::vector<double> dynamic_fixed_ns;
        std::vector<double> dynamic_variable_ns;
        std::vector<double> dynamic_decreasing_ns;

        std::vector<int> insertion_order(9 * size);
        std::iota(insertion_order.begin(), insertion_order.end(), static_cast<int>(size));
        std::shuffle(insertion_order.begin(), insertion_order.end(), rng);

        std::vector<int> removal_order(size);
        std::iota(removal_order.begin(), removal_order.end(), 0);
        std::shuffle(removal_order.begin(), removal_order.end(), rng);

        if (exp == 6) {repetitions /= 10;}

        auto sampler1 = setup_sampler(size, rng);

        for (int rep = 1; rep <= repetitions; ++rep) {

            // Benchmark static sampling
            auto fixed_start = std::chrono::high_resolution_clock::now();
            auto fixed_samples = benchmark_sample_static(rng, sampler1, size);
            auto fixed_end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double, std::nano> fixed_time = fixed_end - fixed_start;
            static_ns.push_back(fixed_time.count() / size);

        }

        auto sampler2 = setup_sampler(size, rng);

        for (int rep = 1; rep <= repetitions; ++rep) {

            // Benchmark dynamic fixed sampling
            auto variable_start = std::chrono::high_resolution_clock::now();
            auto variable_samples = benchmark_sample_dynamic_fixed(rng, sampler2, size);
            auto variable_end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double, std::nano> variable_time = variable_end - variable_start;
            dynamic_fixed_ns.push_back(variable_time.count() / size);

        }
        
        for (int rep = 1; rep <= repetitions; ++rep) {
            auto sampler3 = setup_sampler(size, rng);

            // Benchmark dynamic increasing sampling
            auto variable2_start = std::chrono::high_resolution_clock::now();
            auto variable2_samples = benchmark_sample_dynamic_variable(rng, sampler3, size, insertion_order);
            auto variable2_end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double, std::nano> variable2_time = variable2_end - variable2_start;
            dynamic_variable_ns.push_back(variable2_time.count() / (9*size));
        
        }

        for (int rep = 1; rep <= repetitions; ++rep) {
            auto sampler4 = setup_sampler(size, rng);

            // Benchmark dynamic decreasing sampling
            auto decreasing_start = std::chrono::high_resolution_clock::now();
            auto decreasing_samples = benchmark_sample_dynamic_decreasing(rng, sampler4, removal_order);
            auto decreasing_end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double, std::nano> decreasing_time = decreasing_end - decreasing_start;
            dynamic_decreasing_ns.push_back(decreasing_time.count() / steps);
        }

        static_times.push_back(median(static_ns));
        dynamic_fixed_times.push_back(median(dynamic_fixed_ns));
        dynamic_variable_times.push_back(median(dynamic_variable_ns));
        dynamic_decreasing_times.push_back(median(dynamic_decreasing_ns));

    }

    const std::string data_dir = "../../../data/";

    try {
        append_column_in_place(data_dir + "static.csv", static_times);
        append_column_in_place(data_dir + "dynamic_fixed.csv", dynamic_fixed_times);
        append_column_in_place(data_dir + "dynamic_variable.csv", dynamic_variable_times);
        append_column_in_place(data_dir + "dynamic_decreasing.csv", dynamic_decreasing_times);
        std::cout << "All four files written successfully.\n";
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}