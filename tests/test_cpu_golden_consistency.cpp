#include "NeuronCPU.h"
#include "NeuronSegment_Golden.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <iomanip>

// A simple assertion function for testing
void assert_true(bool condition, const std::string& message) {
    if (!condition) {
        throw std::runtime_error("Assertion failed: " + message);
    }
}

// Test to verify that the NeuronCPU (on-the-fly) implementation produces results
// consistent with the original "Golden" single-segment implementation.
void test_cpu_vs_golden_consistency() {
    std::cout << "Running Test: CPU vs Golden Consistency..." << std::endl;

    // --- 1. Simulation Parameters ---
    double sim_time_ms = 50.0;
    double dt = 0.01;
    int num_steps = static_cast<int>(sim_time_ms / dt);

    // --- 2. Neuron Parameters (single segment) ---
    int num_segments = 1;
    double length_cm = 0.01;
    double diameter_cm = 0.01;
    double Ra_ohm_cm = 150.0;
    double surface_area = M_PI * diameter_cm * length_cm;

    // --- 3. Stimulus ---
    double stim_current_uA = 1.0;
    double stim_start_ms = 10.0;
    double stim_end_ms = 40.0;

    // --- 4. Initialize Models ---
    NeuronCPU cpu_neuron(num_segments, length_cm, diameter_cm, Ra_ohm_cm);
    NeuronSegment_Golden golden_neuron(length_cm, diameter_cm);

    // --- 5. Simulation Loop ---
    double current_time_ms = 0.0;
    for (int i = 0; i < num_steps; ++i) {
        bool stim_on = (current_time_ms >= stim_start_ms && current_time_ms <= stim_end_ms);
        double injected_current = stim_on ? stim_current_uA : 0.0;

        // Update CPU model
        cpu_neuron.set_injected_current(0, injected_current);
        cpu_neuron.update(dt);

        // Update Golden model
        golden_neuron.update(injected_current, dt);

        current_time_ms += dt;
    }

    // --- 6. Compare Final Voltages ---
    double v_cpu = cpu_neuron.get_segment_V(0);
    double v_golden = golden_neuron.get_V();
    double tolerance = 1e-5;

    std::cout << std::fixed << std::setprecision(8);
    std::cout << "  Final Vm (CPU):    " << v_cpu << " mV" << std::endl;
    std::cout << "  Final Vm (Golden): " << v_golden << " mV" << std::endl;
    std::cout << "  Difference:        " << std::abs(v_cpu - v_golden) << " mV" << std::endl;

    assert_true(std::abs(v_cpu - v_golden) < tolerance, "Final voltages do not match between CPU and Golden models.");

    std::cout << "  Test Passed." << std::endl;
}