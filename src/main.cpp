#include <iostream>
#include <vector>
#include <cmath> // For M_PI
#include "HHLookupTables.h"
#include "Neuron.h"

int main() {
    // --- 1. Initialize shared resources ---
    HHLookupTables::initialize();

    // --- 2. Define simulation parameters ---
    const double T_total = 50.0; // Total simulation time (ms)
    const double dt = 0.01;      // Time step (ms)
    const int num_steps = static_cast<int>(T_total / dt);

    // --- 3. Define neuron properties ---
    const int num_segments = 20;
    const double seg_length = 0.01;   // Segment length (cm)
    const double seg_diameter = 0.01; // Segment diameter (cm)
    const double Ra = 35.4;           // Axial resistivity (Ohm-cm)

    // --- 4. Create the Neuron ---
    Neuron neuron(num_segments, seg_length, seg_diameter, Ra);

    // --- 5. Define the stimulus ---
    const double stim_density = 10.0; // Stimulus current density (uA/cm^2)
    const double surface_area = M_PI * seg_diameter * seg_length;
    const double I_stim_total = stim_density * surface_area; // Total stimulus current (uA)

    const double stim_start = 1.0;    // Stimulus start time (ms)
    const double stim_end = 2.0;      // Stimulus end time (ms)
    const int stim_segment = 0;       // Inject current into the first segment

    // --- 6. Run the simulation ---
    std::cout << "Time(ms)\tV_seg0(mV)\tV_seg10(mV)\tV_seg19(mV)\n";

    for (int i = 0; i < num_steps; ++i) {
        double t = i * dt;

        // Apply the stimulus current during the specified window
        if (t >= stim_start && t < stim_end) {
            neuron.set_injected_current(stim_segment, I_stim_total);
        }

        // Update the entire neuron for one time step
        neuron.update(dt);

        // Print the voltage of a few segments every so often
        if (i % 20 == 0) { // Print every 20 steps (every 0.2 ms)
             std::cout << t << "\t"
                       << neuron.get_segment_V(0) << "\t"
                       << neuron.get_segment_V(num_segments / 2) << "\t"
                       << neuron.get_segment_V(num_segments - 1) << "\n";
        }
    }

    return 0;
}
