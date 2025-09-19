#ifndef NEURON_SEGMENT_H
#define NEURON_SEGMENT_H

#include "HHLookupTables.h"

class NeuronSegment {
public:
    // Constructor now takes segment geometry to handle total currents correctly
    NeuronSegment(double length, double diameter);
    NeuronSegment(); // Default constructor for vector resizing

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
};

#endif // NEURON_SEGMENT_H
