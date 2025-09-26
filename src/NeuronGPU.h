#ifndef NEURON_GPU_H
#define NEURON_GPU_H

#include "NeuronModel.h"
#include <vector>
#include <stdexcept>

struct NeuronGPU_Data; // Forward declaration

// Host-side class to manage a multi-compartment neuron simulation on the GPU.
class NeuronGPU : public NeuronModel {
public:
    // Constructor: Allocates memory on the GPU and initializes the neuron state.
    NeuronGPU(int num_segments, double length, double diameter, double Ra);

    // Destructor: Frees all allocated GPU memory.
    ~NeuronGPU();

    // --- Implement NeuronModel interface ---
    void update(double dt) override;
    void set_injected_current(int segment_index, double current_uA) override;
    double get_segment_V(int segment_index) const override;
    std::vector<double> get_all_segment_V() const override;
    int get_num_segments() const override;

private:
    NeuronGPU_Data* d_data; // Pointer to implementation (PIMPL)
    int num_segments;
    std::vector<double> injected_currents_uA; // Host-side buffer
};

#endif // NEURON_GPU_H