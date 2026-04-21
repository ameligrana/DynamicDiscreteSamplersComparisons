
#include <iostream>
#include <random>
#include <vector>
#include <chrono>
#include <cmath>
#include <fstream>
#include <string>
#include <iomanip>
#include <algorithm>
#include <numeric>
#include "DynamicProposalArrayStar.hpp"

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
        lines.push_back(line + ",PROPOSAL_ARRAY");
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

std::vector<double> generate_weights(size_t size, std::mt19937& rng) {
    std::normal_distribution<double> weight_dist(0.0, 1.0);
    std::vector<double> weights(size);

    for (size_t i = 0; i < size; ++i) {
        weights[i] = std::abs(weight_dist(rng));
    }

    return weights;
}

// Setup function for the sampler
sampling::DynamicProposalArrayStar setup_sampler(const std::vector<double>& weights) {
    return sampling::DynamicProposalArrayStar(weights);
}

sampling::DynamicProposalArrayStar setup_sampler(size_t size, std::mt19937& rng) {
    return setup_sampler(generate_weights(size, rng));
}

// Fixed sampling benchmark
std::vector<size_t> benchmark_sample_static(sampling::DynamicProposalArrayStar& sampler, std::mt19937& rng, size_t n) {
    std::vector<size_t> samples;
    samples.reserve(n);

    for (size_t i = 0; i < n; ++i) {
        samples.push_back(sampler.sample(rng));
    }

    return samples;
}

// Increasing-range sampling benchmark
std::vector<size_t> benchmark_sample_dynamic_fixed(sampling::DynamicProposalArrayStar& sampler, std::mt19937& rng, size_t n) {
    std::vector<size_t> samples;
    samples.reserve(n);

    std::uniform_int_distribution<size_t> index_dist(0, n - 1);
    std::normal_distribution<double> weight_dist(0.0, 1.0);

    for (size_t i = 0; i < n; ++i) {
        samples.push_back(sampler.sample(rng));
        size_t index = index_dist(rng);
        double new_weight = std::abs(weight_dist(rng));
        sampler.update(index, new_weight);
    }

    return samples;
}

std::vector<size_t> benchmark_sample_dynamic_variable(
    sampling::DynamicProposalArrayStar& sampler,
    std::mt19937& rng,
    size_t n,
    const std::vector<size_t>& insertion_order
) {
    std::vector<size_t> samples;
    samples.reserve(9*n);

    std::normal_distribution<double> weight_dist(0.0, 1.0);

    for (size_t j = sampler.size_indices(); j < 10 * n; ++j) {
        sampler.push(0.0);
    }

    for (size_t i = 0; i < 9*n; ++i) {
        samples.push_back(sampler.sample(rng));
        sampler.update(insertion_order[i], std::abs(weight_dist(rng)));
    }

    return samples;
}

std::vector<size_t> benchmark_sample_dynamic_decreasing(
    sampling::DynamicProposalArrayStar& sampler,
    std::mt19937& rng,
    const std::vector<size_t>& removal_order
) {
    const size_t n = removal_order.size();
    const size_t steps = n - n / 10;
    std::vector<size_t> samples;
    samples.reserve(steps);

    for (size_t i = 0; i < steps; ++i) {
        samples.push_back(sampler.sample(rng));
        sampler.update(removal_order[i], 0.0);
    }

    return samples;
}

double numerical_check(std::mt19937& rng) {
    std::vector<double> weights;
    weights.reserve(3);
    weights.emplace_back(0.1);
    weights.emplace_back(0.9);
    weights.emplace_back(100000.0);
    sampling::DynamicProposalArrayStar sampler(weights);
    sampler.update(2, 0.0);
    int c = 0;
    for (int i = 0; i < 1000; ++i) {
        if (sampler.sample(rng) == 0) {
            c += 1;
        }
    }
    return double(c) / 1000.0;
}

std::vector<double> decaying_weights_sampling(std::size_t n, std::size_t t) {
    std::vector<double> k(n);
    for (std::size_t i = 0; i < n; ++i) {
        double base = 2.0 + (1.0 / (100.0 * double(n))) * double(i + 1);
        k[i] = std::pow(base, 1000.0);
    }

    sampling::DynamicProposalArrayStar sampler(k);

    for (std::size_t round = 0; round < t; ++round) {
        for (std::size_t i = 0; i < n; ++i) {
            double divisor = 2.0 + (1.0 / (100.0 * double(n))) * double(i + 1);
            k[i] /= divisor;
            sampler.update(i, k[i]);
        }
    }

    std::mt19937 rng(std::random_device{}());
    std::vector<std::size_t> counts(n, 0);
    constexpr std::size_t N_SAMPLES = 1000000;
    for (std::size_t rep = 0; rep < N_SAMPLES; ++rep) {
        std::size_t idx = sampler.sample(rng);
        if (idx < n) {
            ++counts[idx];
        }
    }

    std::vector<double> normalized(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        normalized[i] = double(counts[i]) / double(N_SAMPLES);
    }
    return normalized;
}

int main() {

    const std::size_t n = 100;
    const std::size_t t_max = 51; // it slows down considerably afterwards, so that it seems stuck

    // Open an output CSV file. Each row will be the 1000 normalized values for a given t.
    std::ofstream fout("../../data/proposal_array_numerical.csv");
    if (!fout.is_open()) {
        std::cerr << "Error: could not open decay_normalized.csv for writing.\n";
        return 1;
    }

    // We choose to write each double with 8 digits after the decimal point.
    fout << std::fixed << std::setprecision(8);

    // Loop t from 1..500
    for (std::size_t t = 1; t <= t_max; ++t) {
        // Get the normalized 1000‐length vector for this (n=1000, current t)
        auto normalized = decaying_weights_sampling(n, t);

        // Write them as a single CSV row: normalized[0],normalized[1],...,normalized[999]\n
        for (std::size_t i = 0; i < normalized.size(); ++i) {
            fout << normalized[i];
            if (i + 1 < normalized.size()) {
                fout << ',';
            }
        }
        fout << '\n';

        fout.flush();
    }

    fout.close();

    std::mt19937 rng(42);

    int repetitions = 50;

    std::vector<double> static_times;
    std::vector<double> dynamic_fixed_times;
    std::vector<double> dynamic_variable_times;
    std::vector<double> dynamic_decreasing_times;

    for (int exp = 3; exp <= 7; ++exp) {
        size_t size = static_cast<size_t>(std::pow(10, exp));
        size_t steps = size - size / 10;

        std::vector<double> static_ns;
        std::vector<double> dynamic_fixed_ns;
        std::vector<double> dynamic_variable_ns;
        std::vector<double> dynamic_decreasing_ns;

        std::vector<size_t> insertion_order(9 * size);
        std::iota(insertion_order.begin(), insertion_order.end(), size);
        std::shuffle(insertion_order.begin(), insertion_order.end(), rng);

        std::vector<size_t> removal_order(size);
        std::iota(removal_order.begin(), removal_order.end(), 0);
        std::shuffle(removal_order.begin(), removal_order.end(), rng);

        if (exp == 6) {repetitions /= 10;}

        auto sampler1 = setup_sampler(size, rng);

        for (int rep = 1; rep <= repetitions; ++rep) {

            // Benchmark static sampling
            auto fixed_start = std::chrono::high_resolution_clock::now();
            auto fixed_samples = benchmark_sample_static(sampler1, rng, size);
            auto fixed_end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double, std::nano> fixed_time = fixed_end - fixed_start;
            static_ns.push_back(fixed_time.count() / size);
        }

        auto sampler2 = setup_sampler(size, rng);

        for (int rep = 1; rep <= repetitions; ++rep) {

            // Benchmark dynamic fixed sampling
            auto variable_start = std::chrono::high_resolution_clock::now();
            auto variable_samples = benchmark_sample_dynamic_fixed(sampler2, rng, size);
            auto variable_end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double, std::nano> variable_time = variable_end - variable_start;
            dynamic_fixed_ns.push_back(variable_time.count() / size);
        }
        
        if (exp <= 6) {
            // Stop the increasing-range proposal-array benchmark at 10^6 because it ran out of memory at 10^7.
            for (int rep = 1; rep <= repetitions; ++rep) {
                auto sampler3 = setup_sampler(size, rng);
  
                // Benchmark dynamic increasing sampling
                auto variable2_start = std::chrono::high_resolution_clock::now();
                auto variable2_samples = benchmark_sample_dynamic_variable(sampler3, rng, size, insertion_order);
                auto variable2_end = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double, std::nano> variable2_time = variable2_end - variable2_start;
                dynamic_variable_ns.push_back(variable2_time.count() / (9*size));
            }
        }

        for (int rep = 1; rep <= repetitions; ++rep) {
            auto sampler4 = setup_sampler(size, rng);

            // Benchmark dynamic decreasing sampling
            auto decreasing_start = std::chrono::high_resolution_clock::now();
            auto decreasing_samples = benchmark_sample_dynamic_decreasing(sampler4, rng, removal_order);
            auto decreasing_end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double, std::nano> decreasing_time = decreasing_end - decreasing_start;
            dynamic_decreasing_ns.push_back(decreasing_time.count() / steps);
        }

        static_times.push_back(median(static_ns));
        dynamic_fixed_times.push_back(median(dynamic_fixed_ns));
        dynamic_decreasing_times.push_back(median(dynamic_decreasing_ns));

        if (exp <= 6) {
            dynamic_variable_times.push_back(median(dynamic_variable_ns));
        }

    }

    const std::string data_dir = "../../data/";

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