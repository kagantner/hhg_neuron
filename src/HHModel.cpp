#include "HHModel.h"
#include <cmath>

// Default constructor: Initialize the model to a resting state
HHModel::HHModel() {
    initialize(-65.0);
}

// Initialize the model to a specific voltage
void HHModel::initialize(double V_m_init) {
    V_m = V_m_init;
    // Initialize gating variables to their steady-state values at the given voltage
    m = alpha_m(V_m) / (alpha_m(V_m) + beta_m(V_m));
    h = alpha_h(V_m) / (alpha_h(V_m) + beta_h(V_m));
    n = alpha_n(V_m) / (alpha_n(V_m) + beta_n(V_m));
}

double HHModel::get_V() const { return V_m; }
double HHModel::get_m() const { return m; }
double HHModel::get_h() const { return h; }
double HHModel::get_n() const { return n; }

// Update the model state using on-the-fly calculations
void HHModel::update(double I_stim_density, double dt) {
    // --- Step 1: Update membrane potential using Forward Euler ---
    // All currents are now densities (per cm^2)
    double I_Na = g_Na_density * m * m * m * h * (V_m - E_Na);
    double I_K = g_K_density * n * n * n * n * (V_m - E_K);
    double I_L = g_L_density * (V_m - E_L);

    double dV = (I_stim_density - I_Na - I_K - I_L) / C_m_density;
    V_m += dV * dt;

    // --- Step 2: Update gating variables using Rush-Larsen ---
    double am = alpha_m(V_m);
    double bm = beta_m(V_m);
    double ah = alpha_h(V_m);
    double bh = beta_h(V_m);
    double an = alpha_n(V_m);
    double bn = beta_n(V_m);

    double m_inf = am / (am + bm);
    double tau_m = 1.0 / (am + bm);
    m = m_inf + (m - m_inf) * exp(-dt / tau_m);

    double h_inf = ah / (ah + bh);
    double tau_h = 1.0 / (ah + bh);
    h = h_inf + (h - h_inf) * exp(-dt / tau_h);

    double n_inf = an / (an + bn);
    double tau_n = 1.0 / (an + bn);
    n = n_inf + (n - n_inf) * exp(-dt / tau_n);
}


// --- Private kinetic functions, same as the "Golden" version ---

double HHModel::alpha_m(double V) {
    double V_shifted = V + 40.0;
    if (std::abs(V_shifted) < 1e-5) return 1.0; // Avoid division by zero
    return 0.1 * V_shifted / (1.0 - exp(-V_shifted / 10.0));
}

double HHModel::beta_m(double V) {
    return 4.0 * exp(-(V + 65.0) / 18.0);
}

double HHModel::alpha_h(double V) {
    return 0.07 * exp(-(V + 65.0) / 20.0);
}

double HHModel::beta_h(double V) {
    return 1.0 / (1.0 + exp(-(V + 35.0) / 10.0));
}

double HHModel::alpha_n(double V) {
    double V_shifted = V + 55.0;
    if (std::abs(V_shifted) < 1e-5) return 0.1; // Avoid division by zero
    return 0.01 * V_shifted / (1.0 - exp(-V_shifted / 10.0));
}

double HHModel::beta_n(double V) {
    return 0.125 * exp(-(V + 65.0) / 80.0);
}