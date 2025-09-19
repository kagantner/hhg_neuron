#include <iostream>
#include <vector>
#include <chrono>
#include <iomanip>
#include "NeuronSegment.h"
#include "NeuronSegment_Golden.h"
#include "HHLookupTables.h"

void test_performance() {
    std::cout << "Running Performance Benchmark..." << std::endl;

    HHLookupTables::initialize();

    // Use a large number of steps to get a meaningful measurement
    const int num_steps = 10'000'000;
    const double dt = 0.01;
    const double seg_length = 0.01;
    const double seg_diameter = 0.01;
    const double I_stim = 1.0; // A constant stimulus for simplicity

    // --- Benchmark Golden Model ---
    auto start_golden = std::chrono::high_resolution_clock::now();
    {
        NeuronSegment_Golden golden_segment(seg_length, seg_diameter);
        for (int i = 0; i < num_steps; ++i) {
            golden_segment.update(I_stim, dt);
        }
    }
    auto end_golden = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration_golden = end_golden - start_golden;

    // --- Benchmark Optimized (LUT) Model ---
    auto start_optimized = std::chrono::high_resolution_clock::now();
    {
        NeuronSegment optimized_segment(seg_length, seg_diameter);
        for (int i = 0; i < num_steps; ++i) {
            optimized_segment.update(I_stim, dt);
        }
    }
    auto end_optimized = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration_optimized = end_optimized - start_optimized;

    double speedup = duration_golden.count() / duration_optimized.count();

    std::cout << "  [INFO] Benchmark results (" << num_steps << " steps):" << std::endl;
    std::cout << "    Golden Model execution time:    " << duration_golden.count() << " ms" << std::endl;
    std::cout << "    Optimized Model execution time:  " << duration_optimized.count() << " ms" << std::endl;
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "  [PASS] Optimized model is " << speedup << "x faster." << std::endl;
}