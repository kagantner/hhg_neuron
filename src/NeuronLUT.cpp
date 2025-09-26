#include "NeuronLUT.h"
#include "HHLookupTables.h"
#include <cmath>
#include <stdexcept>
#include <algorithm> // For std::fill

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

NeuronLUT::NeuronLUT(int num_segments, double length, double diameter, double Ra, const HHLookupTables& luts) {
    if (num_segments <= 0) {
        throw std::invalid_argument("Number of segments must be positive.");
    }
    segments.reserve(num_segments);
    for (int i = 0; i < num_segments; ++i) {
        segments.emplace_back(length, diameter, luts);
    }
    injected_currents_uA.resize(num_segments, 0.0);

    if (num_segments > 1) {
        double cross_area_cm2 = M_PI * (diameter / 2.0) * (diameter / 2.0);
        double R_axial_kohm = (Ra * length / cross_area_cm2) / 1000.0;
        g_a = (R_axial_kohm > 0) ? 1.0 / R_axial_kohm : 0.0;
    } else {
        g_a = 0;
    }
}

void NeuronLUT::update(double dt) {
    int n_seg = segments.size();
    if (n_seg == 0) return;

    // Store previous voltages to calculate axial currents
    std::vector<double> prev_V(n_seg);
    for (int i = 0; i < n_seg; ++i) {
        prev_V[i] = segments[i].get_V();
    }

    // Create a copy of segments to hold the next state
    std::vector<NeuronSegment> next_segments = segments;

    // This loop calculates the new state for each segment
    for (int i = 0; i < n_seg; ++i) {
        double I_axial = 0.0;
        if (g_a > 0) {
            double v_left = (i > 0) ? prev_V[i - 1] : prev_V[i];
            double v_right = (i < n_seg - 1) ? prev_V[i + 1] : prev_V[i];
            I_axial = g_a * (v_left - prev_V[i]) + g_a * (v_right - prev_V[i]);
        }

        double total_current = injected_currents_uA[i] + I_axial;
        next_segments[i].update(total_current, dt);
    }

    segments = next_segments;
    std::fill(injected_currents_uA.begin(), injected_currents_uA.end(), 0.0);
}

void NeuronLUT::set_injected_current(int segment_index, double current_uA) {
    if (segment_index >= 0 && static_cast<size_t>(segment_index) < segments.size()) {
        injected_currents_uA[segment_index] = current_uA;
    } else {
        throw std::out_of_range("Segment index out of range.");
    }
}

double NeuronLUT::get_segment_V(int segment_index) const {
    if (segment_index >= 0 && static_cast<size_t>(segment_index) < segments.size()) {
        return segments[segment_index].get_V();
    }
    throw std::out_of_range("Segment index out of range.");
}

std::vector<double> NeuronLUT::get_all_segment_V() const {
    std::vector<double> voltages;
    voltages.reserve(segments.size());
    for (const auto& seg : segments) {
        voltages.push_back(seg.get_V());
    }
    return voltages;
}

int NeuronLUT::get_num_segments() const {
    return segments.size();
}
