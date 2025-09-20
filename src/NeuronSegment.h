#ifndef NEURON_SEGMENT_H
#define NEURON_SEGMENT_H

#include <vector>
#include <cmath>

// Forward declaration to avoid including the full header
class HHLookupTables;

class NeuronSegment {
public:
    // Default constructor is deleted because a LUT reference is required.
    NeuronSegment() = delete;

    // Constructor: Initialize the neuron with geometry and a reference to the LUTs
    NeuronSegment(double length, double diameter, const HHLookupTables& luts);

    // Update the neuron state for a single time step
    void update(double I_total_inj, double dt);

    // Get the current membrane potential in mV
    double get_V() const;

private:
    // State variables
    double V_m; // Membrane potential (mV)
    double m, h, n; // Gating variables

    // Geometric properties
    double surface_area; // in cm^2

    // Pointer to the shared lookup tables
    const HHLookupTables* luts;

    // Biophysical constants (per unit area)
    static constexpr double C_m_density = 1.0;      // Membrane capacitance (uF/cm^2)
    static constexpr double g_Na_density = 120.0;   // Max Na conductance (mS/cm^2)
    static constexpr double g_K_density = 36.0;     // Max K conductance (mS/cm^2)
    static constexpr double g_L_density = 0.3;      // Leak conductance (mS/cm^2)
    static constexpr double E_Na = 50.0;            // Na reversal potential (mV)
    static constexpr double E_K = -77.0;            // K reversal potential (mV)
    static constexpr double E_L = -54.387;          // Leak reversal potential (mV)
};

#endif // NEURON_SEGMENT_H
