#include "Neuron.h"
#include "SynapseModel.h"
#include "HHLookupTables.h"
#include <iostream>
#include <vector>
#include <iomanip> // For std::fixed and std::setprecision

int main() {
    // --- 1. Initialize shared resources ---
    HHLookupTables hh_luts;
    hh_luts.initialize();
    SynapseManager synapse_manager;
    synapse_manager.build_luts(); // Pre-compute all LUTs

    // --- 2. Simulation Parameters ---
    double sim_time_ms = 100.0;
    double dt = 0.01;
    int num_steps = static_cast<int>(sim_time_ms / dt);

    // --- 3. Neuron Parameters ---
    int num_segments = 1;
    double length_cm = 0.1;
    double diameter_cm = 0.02;
    double Ra_ohm_cm = 150.0;

    // --- 4. Create Neuron ---
    Neuron neuron(num_segments, length_cm, diameter_cm, Ra_ohm_cm, hh_luts);

    // --- 5. Simulation Setup ---
    double stimulus_current = 0.1; // uA
    double stimulus_start_ms = 10.0;
    double stimulus_end_ms = 11.0;

    // Set output precision
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "Time(ms)\tVm_first(mV)" << std::endl;

    // --- 6. Simulation Loop ---
    for (int i = 0; i < num_steps; ++i) {
        double current_time_ms = neuron.get_time();

        // Apply stimulus
        if (stimulus_start_ms <= current_time_ms && current_time_ms < stimulus_end_ms) {
            neuron.set_injected_current(0, stimulus_current);
        }

        // Update the neuron state, providing the synapse manager
        neuron.update(dt, synapse_manager);

        // Print membrane potential of first, middle, and last segments
        if (i % 10 == 0) { // Print every 0.1 ms
            std::cout << current_time_ms << "\t"
                      << neuron.get_segment_V(0) << std::endl;
        }
    }

    return 0;
}
