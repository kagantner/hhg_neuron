#include "SynapseModel.h"
#include <cmath> // For std::floor
#include <stdexcept>
#include <algorithm> // For std::max_element

void SynapseManager::build_luts() {
    // --- Build AMPA LUT (1D) ---
    // This implements the ODE solver technique described in the documentation.
    // A simple 2-state (Closed <-> Open) Markov model is used.
    // dO/dt = k_on * NT * (1-O) - k_off * O

    SynapticLUT ampa_lut;
    ampa_lut.E_rev_mv = 0.0f; // Excitatory
    ampa_lut.dimensions = 1;
    ampa_lut.duration_ms = 20.0f;
    ampa_lut.time_step_ms = 0.05f;
    ampa_lut.time_points = static_cast<size_t>(ampa_lut.duration_ms / ampa_lut.time_step_ms);
    ampa_lut.data.resize(ampa_lut.time_points);

    // NOTE: Rate constants are placeholders, as they were not in the spec.
    const float k_on = 0.5f;  // (ms*mM)^-1
    const float k_off = 0.1f; // ms^-1
    const float nt_duration = 1.0f; // ms

    float open_fraction = 0.0f;
    for (size_t i = 0; i < ampa_lut.time_points; ++i) {
        float t = i * ampa_lut.time_step_ms;
        float nt_conc = (t < nt_duration) ? 1.0f : 0.0f; // 1mM pulse of NT

        float dO = (k_on * nt_conc * (1.0f - open_fraction) - k_off * open_fraction) * ampa_lut.time_step_ms;
        open_fraction += dO;

        ampa_lut.data[i] = open_fraction;
    }

    luts.push_back(ampa_lut);

    // --- Build NMDA LUT (2D) ---
    SynapticLUT nmda_lut;
    nmda_lut.E_rev_mv = 0.0f; // Excitatory
    nmda_lut.dimensions = 2;
    nmda_lut.duration_ms = 250.0f;
    nmda_lut.time_step_ms = 0.5f;
    nmda_lut.time_points = static_cast<size_t>(nmda_lut.duration_ms / nmda_lut.time_step_ms);
    nmda_lut.v_min_mv = -90.0f;
    nmda_lut.v_step_mv = 1.0f;
    nmda_lut.voltage_points = static_cast<size_t>((30.0f - -90.0f) / nmda_lut.v_step_mv);
    nmda_lut.data.resize(nmda_lut.time_points * nmda_lut.voltage_points);

    // NOTE: NMDA model is also a placeholder 2-state model
    const float nmda_k_on = 0.1f;
    const float nmda_k_off = 0.01f;

    for (size_t v_idx = 0; v_idx < nmda_lut.voltage_points; ++v_idx) {
        float V = nmda_lut.v_min_mv + v_idx * nmda_lut.v_step_mv;

        // Mg2+ block is voltage-dependent
        // Using a standard sigmoidal formula as a placeholder for calculate_mg_block(v)
        float mg_block = 1.0f / (1.0f + std::exp(-0.062f * V));

        open_fraction = 0.0f; // Reset for each voltage trace
        for (size_t t_idx = 0; t_idx < nmda_lut.time_points; ++t_idx) {
            float t = t_idx * nmda_lut.time_step_ms;
            float nt_conc = (t < nt_duration) ? 1.0f : 0.0f;

            // ODE includes the voltage-dependent block on the off-rate
            float dO = (nmda_k_on * nt_conc * (1.0f - open_fraction) - (nmda_k_off * mg_block) * open_fraction) * nmda_lut.time_step_ms;
            open_fraction += dO;

            size_t index = v_idx * nmda_lut.time_points + t_idx;
            nmda_lut.data[index] = open_fraction;
        }
    }
    luts.push_back(nmda_lut);
}

float SynapseManager::get_conductance(int type_id, float time_since_spike, float voltage_mv) const {
    if (type_id < 0 || static_cast<size_t>(type_id) >= luts.size()) {
        throw std::out_of_range("Synapse type_id is out of range.");
    }

    const auto& lut = luts[type_id];

    // Event is too old to be in the table, so its contribution is zero.
    if (time_since_spike < 0 || time_since_spike >= lut.duration_ms) {
        return 0.0f;
    }

    if (lut.dimensions == 1) {
        return interpolate_1d(lut, time_since_spike);
    } else {
        return interpolate_2d(lut, time_since_spike, voltage_mv);
    }
}

float SynapseManager::interpolate_1d(const SynapticLUT& lut, float time_since_spike) const {
    // Calculate the floating-point index into the LUT data
    float f_index = time_since_spike / lut.time_step_ms;
    int index_0 = static_cast<int>(std::floor(f_index));
    float frac = f_index - index_0;

    // Bounds check: ensure the indices are valid
    if (static_cast<size_t>(index_0 + 1) >= lut.data.size()) {
        return lut.data.back();
    }

    // Get the two points for interpolation
    float val_0 = lut.data[index_0];
    float val_1 = lut.data[index_0 + 1];

    // Perform linear interpolation
    return val_0 + frac * (val_1 - val_0);
}

float SynapseManager::interpolate_2d(const SynapticLUT& lut, float time_since_spike, float voltage_mv) const {
    // Find time index and fraction
    float ft_index = time_since_spike / lut.time_step_ms;
    int t_idx0 = static_cast<int>(std::floor(ft_index));
    float t_frac = ft_index - t_idx0;

    // Find voltage index and fraction
    float fv_index = (voltage_mv - lut.v_min_mv) / lut.v_step_mv;
    int v_idx0 = static_cast<int>(std::floor(fv_index));
    float v_frac = fv_index - v_idx0;

    // Bounds checking for indices
    if (t_idx0 < 0 || static_cast<size_t>(t_idx0 + 1) >= lut.time_points) return 0.0f;
    if (v_idx0 < 0 || static_cast<size_t>(v_idx0 + 1) >= lut.voltage_points) return 0.0f;

    // Get the 4 corner points from the 2D table (flattened to 1D)
    size_t idx00 = static_cast<size_t>(v_idx0) * lut.time_points + static_cast<size_t>(t_idx0);
    size_t idx10 = idx00 + 1;
    size_t idx01 = (static_cast<size_t>(v_idx0) + 1) * lut.time_points + static_cast<size_t>(t_idx0);
    size_t idx11 = idx01 + 1;

    float c00 = lut.data[idx00];
    float c10 = lut.data[idx10];
    float c01 = lut.data[idx01];
    float c11 = lut.data[idx11];

    // Interpolate along time axis for both voltage levels
    float tx0 = c00 + t_frac * (c10 - c00);
    float tx1 = c01 + t_frac * (c11 - c01);

    // Interpolate along voltage axis to get the final value
    return tx0 + v_frac * (tx1 - tx0);
}
