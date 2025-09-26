#ifndef HH_MODEL_H
#define HH_MODEL_H

// This class implements the core Hodgkin-Huxley model equations.
// It calculates alpha and beta values on the fly, without lookup tables.
// This model is designed to be reusable for both CPU and GPU simulations.

class HHModel {
public:
    // Default constructor
    HHModel();

    // Initialize the model to a specific voltage
    void initialize(double V_m_init);

    // Update the model state for a single time step
    // I_stim: Stimulus current in uA/cm^2 (current density)
    // dt: Time step in ms
    void update(double I_stim, double dt);

    // Get the current membrane potential in mV
    double get_V() const;

    // Get the gating variables
    double get_m() const;
    double get_h() const;
    double get_n() const;

private:
    // State variables
    double V_m; // Membrane potential (mV)
    double m, h, n; // Gating variables

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

#endif // HH_MODEL_H