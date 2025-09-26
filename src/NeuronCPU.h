#ifndef NEURON_CPU_H
#define NEURON_CPU_H

#include "NeuronModel.h"
#include "HHModel.h"
#include <vector>
#include <stdexcept>

// Simulates a multi-compartment neuron on the CPU, using the HHModel for on-the-fly calculations.
class NeuronCPU : public NeuronModel {
public:
    // Constructor: Sets up the neuron's geometry and initializes all segments.
    NeuronCPU(int num_segments, double length, double diameter, double Ra);

    // --- Implement NeuronModel interface ---
    void update(double dt) override;
    void set_injected_current(int segment_index, double current_uA) override;
    double get_segment_V(int segment_index) const override;
    std::vector<double> get_all_segment_V() const override;
    int get_num_segments() const override;

private:
    std::vector<HHModel> segments;
    std::vector<double> injected_currents_uA; // Total injected current per segment

    double segment_surface_area; // in cm^2
    double g_a; // Axial conductance in mS
};

#endif // NEURON_CPU_H