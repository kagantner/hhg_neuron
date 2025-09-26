#include "NeuronModel.h"
#include "NeuronLUT.h"
#include "NeuronCPU.h"
#include "NeuronGPU.h"
#include "HHLookupTables.h"

#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include <stdexcept>
#include <iomanip>

void print_usage() {
    std::cerr << "Usage: ./neuron_model <backend>" << std::endl;
    std::cerr << "Available backends: lut, cpu, gpu" << std::endl;
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        print_usage();
        return 1;
    }

    std::string backend = argv[1];

    // --- 1. Simulation Parameters ---
    double sim_time_ms = 100.0;
    double dt = 0.01;
    int num_steps = static_cast<int>(sim_time_ms / dt);

    // --- 2. Neuron Parameters ---
    int num_segments = 50;
    double length_cm = 0.01;      // Length of each segment
    double diameter_cm = 0.002;   // Diameter of the axon
    double Ra_ohm_cm = 150.0;     // Axial resistivity

    // --- 3. Stimulus Parameters ---
    double stim_current_uA = 1.0; // Stimulus current in microamperes
    double stim_start_ms = 10.0;
    double stim_end_ms = 11.0;
    int stim_segment = 0; // Stimulate the first segment

    // --- 4. Select and Initialize Backend ---
    std::unique_ptr<NeuronModel> neuron;

    if (backend == "lut") {
        std::cout << "Using LUT backend." << std::endl;
        // The LUT backend requires pre-initialized lookup tables
        HHLookupTables hh_luts;
        hh_luts.initialize();
        neuron = std::make_unique<NeuronLUT>(num_segments, length_cm, diameter_cm, Ra_ohm_cm, hh_luts);
    } else if (backend == "cpu") {
        std::cout << "Using CPU (on-the-fly) backend." << std::endl;
        neuron = std::make_unique<NeuronCPU>(num_segments, length_cm, diameter_cm, Ra_ohm_cm);
    } else if (backend == "gpu") {
        std::cout << "Using GPU backend." << std::endl;
        try {
            neuron = std::make_unique<NeuronGPU>(num_segments, length_cm, diameter_cm, Ra_ohm_cm);
        } catch (const std::runtime_error& e) {
            std::cerr << "Error initializing GPU backend: " << e.what() << std::endl;
            std::cerr << "Please ensure you have a compatible NVIDIA GPU and CUDA drivers installed." << std::endl;
            return 1;
        }
    } else {
        std::cerr << "Error: Invalid backend '" << backend << "'." << std::endl;
        print_usage();
        return 1;
    }

    // --- 5. Simulation Loop ---
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "Time(ms)\tVm_first(mV)\tVm_mid(mV)\tVm_last(mV)" << std::endl;

    double current_time_ms = 0.0;
    for (int i = 0; i < num_steps; ++i) {
        // Apply stimulus current
        if (current_time_ms >= stim_start_ms && current_time_ms <= stim_end_ms) {
            neuron->set_injected_current(stim_segment, stim_current_uA);
        }

        // Update the neuron state
        neuron->update(dt);

        // Print output at specified intervals
        if (i % 20 == 0) { // Print every 0.2 ms
            std::cout << current_time_ms << "\t"
                      << neuron->get_segment_V(0) << "\t"
                      << neuron->get_segment_V(num_segments / 2) << "\t"
                      << neuron->get_segment_V(num_segments - 1) << std::endl;
        }

        current_time_ms += dt;
    }

    return 0;
}