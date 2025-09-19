#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>
#include <iomanip> // For std::setprecision
#include "NeuronSegment.h"
#include "NeuronSegment_Golden.h"
#include "HHLookupTables.h"

void test_simulation_accuracy() {
    std::cout << "Running Full Simulation Accuracy Test..." << std::endl;

    HHLookupTables::initialize();

    const double T_total = 50.0;
    const double dt = 0.01;
    const int num_steps = static_cast<int>(T_total / dt);
    const double seg_length = 0.01;
    const double seg_diameter = 0.01;
    const double surface_area = M_PI * seg_diameter * seg_length;
    const double stim_density = 10.0;
    const double I_stim_total = stim_density * surface_area;
    const double stim_start = 1.0;
    const double stim_end = 2.0;

    double max_diff = 0.0;

    NeuronSegment optimized_segment(seg_length, seg_diameter);
    NeuronSegment_Golden golden_segment(seg_length, seg_diameter);

    for (int i = 0; i < num_steps; ++i) {
        double t = i * dt;
        double current = 0.0;

        if (t >= stim_start && t < stim_end) {
            current = I_stim_total;
        }

        optimized_segment.update(current, dt);
        golden_segment.update(current, dt);

        double v_opt = optimized_segment.get_V();
        double v_gld = golden_segment.get_V();
        double diff = std::abs(v_opt - v_gld);

        if (diff > max_diff) {
            max_diff = diff;
        }
    }

    std::cout << std::fixed << std::setprecision(10);
    std::cout << "  [INFO] Accuracy report:" << std::endl;
    std::cout << "    Max difference between optimized and golden models: " << max_diff << " mV" << std::endl;

    // We accept that there will be a small, bounded error.
    // The key is that it doesn't grow uncontrollably.
    // Let's assert that the maximum error is still reasonably small.
    const double max_allowed_error = 0.1; // 0.1 mV is a physiologically reasonable error bound.
    assert(max_diff < max_allowed_error);

    std::cout << "  [PASS] Optimized model simulation tracks golden model within acceptable error." << std::endl;
}
