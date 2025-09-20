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
    double sim_time_ms = 50.0; // Shorter simulation to see synapse effect clearly
    double dt = 0.01;
    int num_steps = static_cast<int>(sim_time_ms / dt);

    // --- 3. Neuron Parameters ---
    int num_segments = 1; // Use a single segment for a simple demonstration
    double length_cm = 0.01;
    double diameter_cm = 0.01;
    double Ra_ohm_cm = 150.0;

    // --- 4. Create Neuron and Add Synapse ---
    Neuron neuron(num_segments, length_cm, diameter_cm, Ra_ohm_cm, hh_luts);

    SynapseInstance ampa_synapse;
    ampa_synapse.g_max = 0.5f; // Strong synapse for clear effect (uS)
    ampa_synapse.type_id = 0; // 0 = AMPA, based on the order in build_luts
    ampa_synapse.target_segment = 0;
    ampa_synapse.last_spike_time = -1e9; // Last spike was long ago
    neuron.add_synapse(ampa_synapse);

    // Set output precision
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "Time(ms)\tSegment_0_Vm(mV)" << std::endl;

    // --- 5. Simulation Loop ---
    for (int i = 0; i < num_steps; ++i) {
        // Trigger the synapse at t=10ms
        if (std::abs(neuron.get_time() - 10.0) < dt/2.0) {
            neuron.trigger_synapse(0);
        }

        // Update the neuron state, providing the synapse manager
        neuron.update(dt, synapse_manager);

        // Print membrane potential of the first segment
        if (i % 10 == 0) { // Print every 0.1 ms
            std::cout << neuron.get_time() << "\t"
                      << neuron.get_segment_V(0) << std::endl;
        }
    }

    return 0;
}
