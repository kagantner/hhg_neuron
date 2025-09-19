#ifndef NEURON_H
#define NEURON_H

#include "NeuronSegment.h"
#include <vector>

class Neuron {
public:
    // Constructor: Creates a neuron with a specified number of segments and properties
    // num_segments: The number of compartments in the neuron model
    // length: The length of each segment in cm
    // diameter: The diameter of each segment in cm
    // Ra: The axial resistivity of the axoplasm in Ohm-cm
    Neuron(int num_segments, double length, double diameter, double Ra);

    // Update all segments in the neuron for a single time step
    void update(double dt);

    // Inject a current into a specific segment for the next time step
    void set_injected_current(int segment_index, double current);

    // Get the membrane potential of a specific segment
    double get_segment_V(int segment_index) const;

    // Get the total number of segments
    int get_num_segments() const;

private:
    std::vector<NeuronSegment> segments;
    std::vector<double> injected_currents;

    // Coupling conductance between segments (mS)
    // Note: The units need to be consistent. Hodgkin-Huxley conductances are in mS/cm^2.
    // The axial current is in uA. Let's adjust this.
    double g_a; // Axial conductance in mS
};

#endif // NEURON_H
