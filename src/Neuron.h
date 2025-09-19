#ifndef NEURON_H
#define NEURON_H

#include "NeuronSegment.h"
#include "SynapseModel.h"
#include <vector>

class SynapseManager; // Forward declaration

class Neuron {
public:
    Neuron(int num_segments, double length, double diameter, double Ra);

    // Update all segments, calculating synaptic currents using the provided manager
    void update(double dt, const SynapseManager& synapse_manager);

    // Add a synapse instance to the neuron
    void add_synapse(const SynapseInstance& synapse);

    // Trigger a synapse by its index in the neuron's synapse vector
    void trigger_synapse(int synapse_idx);

    void set_injected_current(int segment_index, double current);
    double get_segment_V(int segment_index) const;
    int get_num_segments() const;
    double get_time() const;

private:
    std::vector<NeuronSegment> segments;
    std::vector<double> injected_currents;
    std::vector<SynapseInstance> synapses;
    double time_ms; // Current simulation time

    double g_a; // Axial conductance in mS
};

#endif // NEURON_H
