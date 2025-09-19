#include "../src/SynapseModel.h"
#include <iostream>
#include <cmath>
#include <stdexcept>
#include <sstream>

// Helper to check for approximate equality of floats
void check_approx(float actual, float expected, float tolerance = 1e-6) {
    if (std::abs(actual - expected) > tolerance) {
        std::stringstream ss;
        ss << "Assertion failed: actual=" << actual << ", expected=" << expected;
        throw std::runtime_error(ss.str());
    }
}

void test_lut_synapse() {
    std::cout << "Running LUT Synapse Test..." << std::endl;

    // 1. Setup the manager and build the LUTs
    SynapseManager sm;
    sm.build_luts();

    // 2. Basic checks on the generated AMPA LUT (type_id = 0)
    if (sm.luts.empty()) throw std::runtime_error("No LUTs were built.");
    const auto& ampa_lut = sm.luts[0];
    if (ampa_lut.time_points == 0) throw std::runtime_error("AMPA LUT has zero time points.");
    if (ampa_lut.data.empty()) throw std::runtime_error("AMPA LUT data is empty.");

    // The ODE solver produces a normalized fraction of open channels directly,
    // so we just check that the LUT was populated.
    std::cout << "  [PASS] LUT generated correctly." << std::endl;

    // 3. Test the interpolation logic
    // We will pick a time, calculate the expected value by hand, and compare.

    // Test case 1: A time that falls exactly on a grid point
    float t1 = 2.0f; // ms
    float expected1 = sm.luts[0].data[static_cast<int>(t1 / ampa_lut.time_step_ms)];
    float actual1 = sm.get_conductance(0, t1, 0); // voltage doesn't matter for 1D
    check_approx(actual1, expected1);
    std::cout << "  [PASS] Interpolation at grid point." << std::endl;

    // Test case 2: A time that falls between two grid points
    float t2 = 2.025f; // Exactly halfway between 2.00 and 2.05
    float f_index = t2 / ampa_lut.time_step_ms;
    int index_0 = static_cast<int>(f_index);
    float val_0 = ampa_lut.data[index_0];
    float val_1 = ampa_lut.data[index_0 + 1];
    float expected2 = val_0 + 0.5f * (val_1 - val_0); // Manual interpolation
    float actual2 = sm.get_conductance(0, t2, 0);
    check_approx(actual2, expected2);
    std::cout << "  [PASS] Interpolation between grid points." << std::endl;

    // Test case 3: Time before the start (should be 0)
    float actual3 = sm.get_conductance(0, -1.0f, 0);
    check_approx(actual3, 0.0f);
    std::cout << "  [PASS] Interpolation before start time." << std::endl;

    // Test case 4: Time after the end (should be 0)
    float actual4 = sm.get_conductance(0, ampa_lut.duration_ms + 1.0f, 0);
    check_approx(actual4, 0.0f);
    std::cout << "  [PASS] Interpolation after end time." << std::endl;

    // 4. Check reversal potentials
    check_approx(sm.luts[0].E_rev_mv, 0.0f); // AMPA
    check_approx(sm.luts[1].E_rev_mv, 0.0f); // NMDA
    std::cout << "  [PASS] Reversal potentials are correct." << std::endl;

    // 5. Test 2D bilinear interpolation for NMDA (type_id = 1)
    const auto& nmda_lut = sm.luts[1];
    float t_2d = 50.25f; // Halfway between 50.0 and 50.5
    float v_2d = -30.5f; // Halfway between -31.0 and -30.0

    // Manually calculate expected value
    int t_idx0 = static_cast<int>(t_2d / nmda_lut.time_step_ms);
    int v_idx0 = static_cast<int>((v_2d - nmda_lut.v_min_mv) / nmda_lut.v_step_mv);

    int idx00 = v_idx0 * nmda_lut.time_points + t_idx0;
    float c00 = nmda_lut.data[idx00];
    float c10 = nmda_lut.data[idx00 + 1];
    float c01 = nmda_lut.data[idx00 + nmda_lut.time_points];
    float c11 = nmda_lut.data[idx00 + nmda_lut.time_points + 1];

    float tx0 = c00 + 0.5f * (c10 - c00);
    float tx1 = c01 + 0.5f * (c11 - c01);
    float expected_2d = tx0 + 0.5f * (tx1 - tx0);

    float actual_2d = sm.get_conductance(1, t_2d, v_2d);
    check_approx(actual_2d, expected_2d);
    std::cout << "  [PASS] 2D Bilinear interpolation." << std::endl;


    std::cout << "  ...test passed." << std::endl;
}
