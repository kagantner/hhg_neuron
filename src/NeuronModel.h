#ifndef NEURON_MODEL_H
#define NEURON_MODEL_H

#include <vector>

// Abstract base class for all neuron simulation models.
// This defines a common interface for different backends (LUT, CPU, GPU).
class NeuronModel {
public:
    // Virtual destructor to ensure proper cleanup of derived classes.
    virtual ~NeuronModel() = default;

    // Update the neuron state for a single time step.
    // dt: The time step in ms.
    virtual void update(double dt) = 0;

    // Apply an external injected current to a specific segment.
    // segment_index: The index of the segment to apply the current to.
    // current_uA: The total current in microamperes.
    virtual void set_injected_current(int segment_index, double current_uA) = 0;

    // Get the membrane potential of a specific segment.
    virtual double get_segment_V(int segment_index) const = 0;

    // Get the membrane potentials of all segments.
    virtual std::vector<double> get_all_segment_V() const = 0;

    // Get the total number of segments in the neuron.
    virtual int get_num_segments() const = 0;
};

#endif // NEURON_MODEL_H