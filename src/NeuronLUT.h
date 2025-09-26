#ifndef NEURON_LUT_H
#define NEURON_LUT_H

#include "NeuronModel.h"
#include "NeuronSegment.h"
#include <vector>

class HHLookupTables; // Forward declaration

// The original Neuron class, refactored to be a concrete implementation of NeuronModel.
// This version uses pre-calculated lookup tables for performance.
class NeuronLUT : public NeuronModel {
public:
    NeuronLUT(int num_segments, double length, double diameter, double Ra, const HHLookupTables& luts);

    // --- Implement NeuronModel interface ---
    void update(double dt) override;
    void set_injected_current(int segment_index, double current_uA) override;
    double get_segment_V(int segment_index) const override;
    std::vector<double> get_all_segment_V() const override;
    int get_num_segments() const override;

private:
    std::vector<NeuronSegment> segments;
    std::vector<double> injected_currents_uA;
    double g_a; // Axial conductance in mS
};

#endif // NEURON_LUT_H
