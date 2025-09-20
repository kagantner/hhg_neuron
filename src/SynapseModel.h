#ifndef SYNAPSE_MODEL_H
#define SYNAPSE_MODEL_H

#include <vector>
#include <string>

// A shared Look-Up Table for a specific type of synapse (e.g., AMPA, NMDA)
struct SynapticLUT {
    int dimensions = 1;

    // Time dimension
    float duration_ms = 0.0f;
    float time_step_ms = 0.0f;
    size_t time_points = 0;

    // Optional voltage dimension (for 2D tables)
    float v_min_mv = 0.0f;
    float v_step_mv = 0.0f;
    size_t voltage_points = 0;

    // The pre-computed, normalized conductance data
    std::vector<float> data;
    float E_rev_mv = 0.0f; // Reversal potential for this synapse type
};

// An instance of a single synapse in the neuron model
struct SynapseInstance {
    float g_max;           // Max conductance (unique strength of this synapse) in uS
    double last_spike_time;  // Time of the last presynaptic event in ms
    int type_id;           // Index to select the correct SynapticLUT from the manager
    int target_segment;    // The segment of the postsynaptic neuron it connects to
};

// Manages the creation and storage of all synaptic look-up tables
class SynapseManager {
public:
    // Builds all LUTs based on their underlying models.
    // Must be called once before the simulation starts.
    void build_luts();

    // Gets the interpolated conductance for a given synapse type.
    // time_since_spike: in ms
    // voltage_mv: current postsynaptic membrane potential in mV
    // Returns: normalized conductance (0.0 to 1.0)
    float get_conductance(int type_id, float time_since_spike, float voltage_mv) const;

    // Stores all the pre-computed LUTs
    std::vector<SynapticLUT> luts;

private:
    // Helper for 1D linear interpolation
    float interpolate_1d(const SynapticLUT& lut, float time_since_spike) const;

    // Helper for 2D bilinear interpolation
    float interpolate_2d(const SynapticLUT& lut, float time_since_spike, float voltage_mv) const;
};

#endif // SYNAPSE_MODEL_H
