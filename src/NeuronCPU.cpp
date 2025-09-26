#include "NeuronCPU.h"
#include <cmath>
#include <vector>
#include <stdexcept>
#include <algorithm> // For std::fill

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

NeuronCPU::NeuronCPU(int num_segments, double length, double diameter, double Ra) {
    if (num_segments <= 0) {
        throw std::invalid_argument("Number of segments must be positive.");
    }
    segments.resize(num_segments);
    injected_currents_uA.resize(num_segments, 0.0);

    segment_surface_area = M_PI * diameter * length; // cm^2

    // Calculate axial conductance between segments
    if (num_segments > 1) {
        double cross_area_cm2 = M_PI * (diameter / 2.0) * (diameter / 2.0);
        // Ra is in Ohm-cm, R_axial is in kOhm
        double R_axial_kohm = (Ra * length / cross_area_cm2) / 1000.0;
        if (R_axial_kohm > 0) {
            g_a = 1.0 / R_axial_kohm; // mS
        } else {
            g_a = 0;
        }
    } else {
        g_a = 0;
    }
}

void NeuronCPU::update(double dt) {
    int n_seg = segments.size();
    if (n_seg == 0) return;

    // Store previous voltages to calculate axial currents
    std::vector<double> prev_V(n_seg);
    for (int i = 0; i < n_seg; ++i) {
        prev_V[i] = segments[i].get_V();
    }

    // This loop calculates the new state for each segment
    for (int i = 0; i < n_seg; ++i) {
        // Calculate axial current from neighbors (in uA)
        double I_axial = 0.0;
        if (g_a > 0) {
            // Voltage of the left neighbor. Boundary condition: no-flux (current = 0), so V_left = V_me
            double v_left = (i > 0) ? prev_V[i - 1] : prev_V[i];
            // Voltage of the right neighbor. Boundary condition: no-flux, so V_right = V_me
            double v_right = (i < n_seg - 1) ? prev_V[i + 1] : prev_V[i];
            // I_axial is in uA because g_a is mS (1/kOhm) and V is mV. uA = mS * mV.
            I_axial = g_a * (v_left - prev_V[i]) + g_a * (v_right - prev_V[i]);
        }

        // Total current is the sum of injected and axial currents
        double total_current_uA = injected_currents_uA[i] + I_axial;

        // Convert total current to current density for the HHModel (uA/cm^2)
        double I_stim_density = (segment_surface_area > 0) ? total_current_uA / segment_surface_area : 0.0;

        // Update the segment's state
        segments[i].update(I_stim_density, dt);
    }

    // Reset injected currents for the next time step
    std::fill(injected_currents_uA.begin(), injected_currents_uA.end(), 0.0);
}

void NeuronCPU::set_injected_current(int segment_index, double current_uA) {
    if (segment_index >= 0 && static_cast<size_t>(segment_index) < segments.size()) {
        injected_currents_uA[segment_index] = current_uA;
    } else {
        throw std::out_of_range("Segment index out of range.");
    }
}

double NeuronCPU::get_segment_V(int segment_index) const {
    if (segment_index >= 0 && static_cast<size_t>(segment_index) < segments.size()) {
        return segments[segment_index].get_V();
    }
    throw std::out_of_range("Segment index out of range.");
}

int NeuronCPU::get_num_segments() const {
    return segments.size();
}

std::vector<double> NeuronCPU::get_all_segment_V() const {
    std::vector<double> voltages;
    voltages.reserve(segments.size());
    for (const auto& seg : segments) {
        voltages.push_back(seg.get_V());
    }
    return voltages;
}