#include "Neuron.h"
#include <cmath>
#include <stdexcept>
#include <algorithm> // For std::fill

// Constructor
Neuron::Neuron(int num_segments, double length, double diameter, double Ra) {
    if (num_segments <= 0) {
        throw std::invalid_argument("Number of segments must be positive.");
    }

    segments.reserve(num_segments);
    for (int i = 0; i < num_segments; ++i) {
        segments.emplace_back(length, diameter);
    }

    injected_currents.resize(num_segments, 0.0);

    if (num_segments > 1) {
        double cross_area = M_PI * (diameter / 2.0) * (diameter / 2.0);
        double R_axial_kohm = (Ra * length / cross_area) / 1000.0;
        g_a = 1.0 / R_axial_kohm;
    } else {
        g_a = 0;
    }
}

// Update all segments
void Neuron::update(double dt) {
    int n_seg = segments.size();
    if (n_seg == 0) return;

    // For a single segment, there's no axial current.
    if (n_seg == 1) {
        segments[0].update(injected_currents[0], dt);
        std::fill(injected_currents.begin(), injected_currents.end(), 0.0);
        return;
    }

    // --- Jacobi-style update to prevent data dependencies within a time step ---

    // 1. Store all voltages from the previous time step (t).
    std::vector<double> prev_V(n_seg);
    for (int i = 0; i < n_seg; ++i) {
        prev_V[i] = segments[i].get_V();
    }

    // 2. Calculate all updates based on the state at time (t) and store in a temporary copy.
    std::vector<NeuronSegment> next_segments = segments;
    for (int i = 0; i < n_seg; ++i) {
        // Sealed-end boundary condition
        double v_left = (i > 0) ? prev_V[i - 1] : prev_V[i];
        double v_right = (i < n_seg - 1) ? prev_V[i + 1] : prev_V[i];

        double I_axial = g_a * (v_left - prev_V[i]) + g_a * (v_right - prev_V[i]);
        double total_current = injected_currents[i] + I_axial;

        // Update the copy, not the original
        next_segments[i].update(total_current, dt);
    }

    // 3. Atomically swap the updated state back.
    segments = next_segments;

    // 4. Reset injected currents for the next time step.
    std::fill(injected_currents.begin(), injected_currents.end(), 0.0);
}

// Set the injected current (in uA) for a specific segment
void Neuron::set_injected_current(int segment_index, double current_uA) {
    if (segment_index >= 0 && static_cast<size_t>(segment_index) < segments.size()) {
        injected_currents[segment_index] = current_uA;
    }
}

// Get the membrane potential of a specific segment
double Neuron::get_segment_V(int segment_index) const {
    if (segment_index >= 0 && static_cast<size_t>(segment_index) < segments.size()) {
        return segments[segment_index].get_V();
    }
    throw std::out_of_range("Segment index out of range.");
}

// Get the total number of segments
int Neuron::get_num_segments() const {
    return segments.size();
}
