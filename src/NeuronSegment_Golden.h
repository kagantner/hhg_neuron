#ifndef NEURON_SEGMENT_GOLDEN_H
#define NEURON_SEGMENT_GOLDEN_H

// This is the "Golden" reference version of the NeuronSegment.
// It calculates alpha and beta values on the fly instead of using lookup tables.
// It is used as a ground truth to verify the accuracy of the optimized version.

class NeuronSegment_Golden {
public:
    // Constructor takes segment geometry to handle total currents correctly
    NeuronSegment_Golden(double length, double diameter);

    // Update the neuron state for a single time step
    // I_total_inj: Total injected current in uA
    // dt: Time step in ms
    void update(double I_total_inj, double dt);

    // Get the current membrane potential in mV
    double get_V() const;

private:
    // State variables
    double V_m; // Membrane potential (mV)
    double m, h, n; // Gating variables

    // Geometric properties
    double surface_area; // in cm^2

    // Biophysical constants (per unit area)
    static constexpr double C_m_density = 1.0;      // Membrane capacitance (uF/cm^2)
    static constexpr double g_Na_density = 120.0;   // Max Na conductance (mS/cm^2)
    static constexpr double g_K_density = 36.0;     // Max K conductance (mS/cm^2)
    static constexpr double g_L_density = 0.3;      // Leak conductance (mS/cm^2)
    static constexpr double E_Na = 50.0;            // Na reversal potential (mV)
    static constexpr double E_K = -77.0;            // K reversal potential (mV)
    static constexpr double E_L = -54.387;          // Leak reversal potential (mV)

    // Private functions to calculate alpha and beta for gating variables
    double alpha_m(double V);
    double beta_m(double V);
    double alpha_h(double V);
    double beta_h(double V);
    double alpha_n(double V);
    double beta_n(double V);
};

#endif // NEURON_SEGMENT_GOLDEN_H
