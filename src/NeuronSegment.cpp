#include "NeuronSegment.h"
#include "HHLookupTables.h"
#include <cmath>
#include <algorithm>

// Constructor: Initialize the neuron with specific geometry and to a resting state
NeuronSegment::NeuronSegment(double length, double diameter, const HHLookupTables& luts)
    : V_m(-65.0), m(0), h(0), n(0), surface_area(0), luts(&luts) {
    V_m = -65.0;
    surface_area = M_PI * diameter * length;

    // Use direct calculation for initial state to match golden model perfectly.
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

namespace {
    // Use a fast, accurate Padé approximant for exp(-x) for small x > 0
    // This is equivalent to the Crank-Nicolson method for this ODE
    inline double fast_exp_approx(double x) {
        return (1.0 - 0.5 * x) / (1.0 + 0.5 * x);
    }
}

// Update the neuron state using Rush-Larsen and interpolated LUTs
void NeuronSegment::update(double I_total_inj, double dt) {
    if (surface_area <= 0) return;

    // We store the original V_m, update the potential, then use the stored
    // value for the kinetic calculations. This improves accuracy and stability.
    const double V_at_start_of_step = V_m;

    // --- Step 1: Update gating variables using Rush-Larsen ---
    // Get alpha and beta values by interpolating from the LUTs using the voltage
    // from the beginning of the time step.

    // Pre-calculate LUT interpolation parameters
    const double pos = (V_at_start_of_step - luts->V_MIN) / luts->V_STEP;
    const int    idx_base = static_cast<int>(pos);
    const double frac = pos - idx_base;
    const int    idx = std::max(0, std::min(static_cast<int>(HHLookupTables::LUT_SIZE) - 2, idx_base));

    // Interpolate all alpha/beta values
    const double am = luts->alpha_m_lut[idx] + frac * (luts->alpha_m_lut[idx+1] - luts->alpha_m_lut[idx]);
    const double bm = luts->beta_m_lut[idx]  + frac * (luts->beta_m_lut[idx+1]  - luts->beta_m_lut[idx]);
    const double ah = luts->alpha_h_lut[idx] + frac * (luts->alpha_h_lut[idx+1] - luts->alpha_h_lut[idx]);
    const double bh = luts->beta_h_lut[idx]  + frac * (luts->beta_h_lut[idx+1]  - luts->beta_h_lut[idx]);
    const double an = luts->alpha_n_lut[idx] + frac * (luts->alpha_n_lut[idx+1] - luts->alpha_n_lut[idx]);
    const double bn = luts->beta_n_lut[idx]  + frac * (luts->beta_n_lut[idx+1]  - luts->beta_n_lut[idx]);

    // Update m gate
    const double sum_m = am + bm;
    const double m_inf = am / sum_m;
    m = m_inf + (m - m_inf) * fast_exp_approx(dt * sum_m);

    // Update h gate
    const double sum_h = ah + bh;
    const double h_inf = ah / sum_h;
    h = h_inf + (h - h_inf) * fast_exp_approx(dt * sum_h);

    // Update n gate
    const double sum_n = an + bn;
    const double n_inf = an / sum_n;
    n = n_inf + (n - n_inf) * fast_exp_approx(dt * sum_n);

    // --- Step 2: Update membrane potential using Forward Euler ---
    // Calculate ionic currents using the *updated* gating variables and the
    // voltage from the *start* of the step.
    const double I_Na = g_Na_density * surface_area * m * m * m * h * (V_at_start_of_step - E_Na);
    const double I_K = g_K_density * surface_area * n * n * n * n * (V_at_start_of_step - E_K);
    const double I_L = g_L_density * surface_area * (V_at_start_of_step - E_L);

    // Update voltage
    const double C_total = C_m_density * surface_area;
    const double dV = (I_total_inj - I_Na - I_K - I_L) / C_total;
    V_m += dV * dt;
}
