#include "Neuron.h"
#include "SynapseModel.h"
#include "HHLookupTables.h"
#include <cmath>
#include <stdexcept>
#include <algorithm> // For std::fill

Neuron::Neuron(int num_segments, double length, double diameter, double Ra, const HHLookupTables& luts) : time_ms(0.0) {
    if (num_segments <= 0) {
        throw std::invalid_argument("Number of segments must be positive.");
    }
    // Use emplace_back since NeuronSegment no longer has a default constructor
    segments.reserve(num_segments);
    for (int i = 0; i < num_segments; ++i) {
        segments.emplace_back(length, diameter, luts);
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

void Neuron::update(double dt, const SynapseManager& synapse_manager) {
    time_ms += dt;
    int n_seg = segments.size();
    if (n_seg == 0) return;

    std::vector<double> prev_V(n_seg);
    for (int i = 0; i < n_seg; ++i) {
        prev_V[i] = segments[i].get_V();
    }

    std::vector<double> synaptic_currents(n_seg, 0.0);
    for (const auto& syn : synapses) {
        float time_since_spike = static_cast<float>(time_ms - syn.last_spike_time);
        if (time_since_spike < 0) continue;

        float g_norm = synapse_manager.get_conductance(syn.type_id, time_since_spike, prev_V[syn.target_segment]);

        // Get the reversal potential from the LUT for this synapse type
        const auto& lut = synapse_manager.luts[syn.type_id];
        float E_rev = lut.E_rev_mv;

        float current_nA = syn.g_max * g_norm * (prev_V[syn.target_segment] - E_rev);
        synaptic_currents[syn.target_segment] += current_nA / 1000.0; // Convert nA to uA
    }

    std::vector<NeuronSegment> next_segments = segments;
    for (int i = 0; i < n_seg; ++i) {
        double I_axial = 0.0;
        if (n_seg > 1) {
            double v_left = (i > 0) ? prev_V[i - 1] : prev_V[i];
            double v_right = (i < n_seg - 1) ? prev_V[i + 1] : prev_V[i];
            I_axial = g_a * (v_left - prev_V[i]) + g_a * (v_right - prev_V[i]);
        }
        double total_current = injected_currents[i] + synaptic_currents[i] + I_axial;
        next_segments[i].update(total_current, dt);
    }

    segments = next_segments;
    std::fill(injected_currents.begin(), injected_currents.end(), 0.0);
}

void Neuron::add_synapse(const SynapseInstance& synapse) {
    if (synapse.target_segment < 0 || static_cast<size_t>(synapse.target_segment) >= segments.size()) {
        throw std::out_of_range("Synapse target segment index out of range.");
    }
    synapses.push_back(synapse);
}

void Neuron::trigger_synapse(int synapse_idx) {
    if (synapse_idx < 0 || static_cast<size_t>(synapse_idx) >= synapses.size()) {
        throw std::out_of_range("Synapse index out of range.");
    }
    synapses[synapse_idx].last_spike_time = time_ms;
}

void Neuron::set_injected_current(int segment_index, double current_uA) {
    if (segment_index >= 0 && static_cast<size_t>(segment_index) < segments.size()) {
        injected_currents[segment_index] = current_uA;
    }
}

double Neuron::get_segment_V(int segment_index) const {
    if (segment_index >= 0 && static_cast<size_t>(segment_index) < segments.size()) {
        return segments[segment_index].get_V();
    }
    throw std::out_of_range("Segment index out of range.");
}

int Neuron::get_num_segments() const {
    return segments.size();
}

double Neuron::get_time() const {
    return time_ms;
}
