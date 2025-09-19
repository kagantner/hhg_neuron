#include "NeuronSegment_Golden.h"
#include <cmath>
#include <algorithm>

// Constructor: Initialize the neuron with specific geometry and to a resting state
NeuronSegment_Golden::NeuronSegment_Golden(double length, double diameter) {
    V_m = -65.0;
    surface_area = M_PI * diameter * length;

    // Initialize gating variables to their steady-state values at resting potential
    m = alpha_m(V_m) / (alpha_m(V_m) + beta_m(V_m));
    h = alpha_h(V_m) / (alpha_h(V_m) + beta_h(V_m));
    n = alpha_n(V_m) / (alpha_n(V_m) + beta_n(V_m));
}

double NeuronSegment_Golden::get_V() const {
    return V_m;
}

// Update the neuron state using on-the-fly calculations
void NeuronSegment_Golden::update(double I_total_inj, double dt) {
    if (surface_area <= 0) return;

    // --- Step 1: Update membrane potential using Forward Euler ---
    double I_Na = g_Na_density * surface_area * m * m * m * h * (V_m - E_Na);
    double I_K = g_K_density * surface_area * n * n * n * n * (V_m - E_K);
    double I_L = g_L_density * surface_area * (V_m - E_L);

    double C_total = C_m_density * surface_area;
    double dV = (I_total_inj - I_Na - I_K - I_L) / C_total;
    V_m += dV * dt;

    // --- Step 2: Update gating variables using Rush-Larsen ---
    // Get alpha and beta values by calculating them directly
    double am = alpha_m(V_m);
    double bm = beta_m(V_m);
    double ah = alpha_h(V_m);
    double bh = beta_h(V_m);
    double an = alpha_n(V_m);
    double bn = beta_n(V_m);

    // Update m gate
    double m_inf = am / (am + bm);
    double tau_m = 1.0 / (am + bm);
    m = m_inf + (m - m_inf) * exp(-dt / tau_m);

    // Update h gate
    double h_inf = ah / (ah + bh);
    double tau_h = 1.0 / (ah + bh);
    h = h_inf + (h - h_inf) * exp(-dt / tau_h);

    // Update n gate
    double n_inf = an / (an + bn);
    double tau_n = 1.0 / (an + bn);
    n = n_inf + (n - n_inf) * exp(-dt / tau_n);
}


// --- Private kinetic functions, copied from HHLookupTables ---

double NeuronSegment_Golden::alpha_m(double V) {
    if (std::abs(V + 40.0) < 1e-5) return 1.0;
    return 0.1 * (V + 40.0) / (1.0 - exp(-(V + 40.0) / 10.0));
}

double NeuronSegment_Golden::beta_m(double V) {
    return 4.0 * exp(-(V + 65.0) / 18.0);
}

double NeuronSegment_Golden::alpha_h(double V) {
    return 0.07 * exp(-(V + 65.0) / 20.0);
}

double NeuronSegment_Golden::beta_h(double V) {
    return 1.0 / (1.0 + exp(-(V + 35.0) / 10.0));
}

double NeuronSegment_Golden::alpha_n(double V) {
    if (std::abs(V + 55.0) < 1e-5) return 0.1;
    return 0.01 * (V + 55.0) / (1.0 - exp(-(V + 55.0) / 10.0));
}

double NeuronSegment_Golden::beta_n(double V) {
    return 0.125 * exp(-(V + 65.0) / 80.0);
}
