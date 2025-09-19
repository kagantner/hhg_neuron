#include "NeuronSegment.h"
#include "HHLookupTables.h"
#include <cmath>
#include <algorithm>

// Default constructor for vector pre-allocation
NeuronSegment::NeuronSegment() : V_m(-65.0), m(0), h(0), n(0), surface_area(0) {}

// Constructor: Initialize the neuron with specific geometry and to a resting state
NeuronSegment::NeuronSegment(double length, double diameter) {
    V_m = -65.0;
    surface_area = M_PI * diameter * length;

    if (HHLookupTables::alpha_m_lut.empty()) {
        HHLookupTables::initialize();
    }

    // Use direct calculation for initial state to match golden model perfectly
    auto temp_alpha_m = [](double V) { if (std::abs(V + 40.0) < 1e-5) return 1.0; return 0.1 * (V + 40.0) / (1.0 - exp(-(V + 40.0) / 10.0)); };
    auto temp_beta_m = [](double V) { return 4.0 * exp(-(V + 65.0) / 18.0); };
    auto temp_alpha_h = [](double V) { return 0.07 * exp(-(V + 65.0) / 20.0); };
    auto temp_beta_h = [](double V) { return 1.0 / (1.0 + exp(-(V + 35.0) / 10.0)); };
    auto temp_alpha_n = [](double V) { if (std::abs(V + 55.0) < 1e-5) return 0.1; return 0.01 * (V + 55.0) / (1.0 - exp(-(V + 55.0) / 10.0)); };
    auto temp_beta_n = [](double V) { return 0.125 * exp(-(V + 65.0) / 80.0); };

    double initial_V = -65.0;
    m = temp_alpha_m(initial_V) / (temp_alpha_m(initial_V) + temp_beta_m(initial_V));
    h = temp_alpha_h(initial_V) / (temp_alpha_h(initial_V) + temp_beta_h(initial_V));
    n = temp_alpha_n(initial_V) / (temp_alpha_n(initial_V) + temp_beta_n(initial_V));
}

double NeuronSegment::get_V() const {
    return V_m;
}

// Helper function for linear interpolation on the LUT
inline double interpolate_lut(double V, const std::vector<double>& lut) {
    double pos = (V - HHLookupTables::V_MIN) / HHLookupTables::V_STEP;
    int idx = static_cast<int>(pos);
    double frac = pos - idx;

    idx = std::max(0, std::min(HHLookupTables::LUT_SIZE - 2, idx));

    double val1 = lut[idx];
    double val2 = lut[idx + 1];

    return val1 + frac * (val2 - val1);
}

// Update the neuron state using Rush-Larsen and interpolated LUTs
void NeuronSegment::update(double I_total_inj, double dt) {
    if (surface_area <= 0) return;

    double I_Na = g_Na_density * surface_area * m * m * m * h * (V_m - E_Na);
    double I_K = g_K_density * surface_area * n * n * n * n * (V_m - E_K);
    double I_L = g_L_density * surface_area * (V_m - E_L);

    double C_total = C_m_density * surface_area;
    double dV = (I_total_inj - I_Na - I_K - I_L) / C_total;
    V_m += dV * dt;

    double am = interpolate_lut(V_m, HHLookupTables::alpha_m_lut);
    double bm = interpolate_lut(V_m, HHLookupTables::beta_m_lut);
    double ah = interpolate_lut(V_m, HHLookupTables::alpha_h_lut);
    double bh = interpolate_lut(V_m, HHLookupTables::beta_h_lut);
    double an = interpolate_lut(V_m, HHLookupTables::alpha_n_lut);
    double bn = interpolate_lut(V_m, HHLookupTables::beta_n_lut);

    // Update gates with Rush-Larsen, matching the golden model's structure
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
