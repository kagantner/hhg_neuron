#include "HHLookupTables.h"
#include <cmath>
#include <vector>

// --- Private kinetic functions, used only for initialization ---

double HHLookupTables::alpha_m(double V) {
    // A small offset to avoid division by zero at V = -40
    if (std::abs(V + 40.0) < 1e-5) {
        return 1.0;
    }
    return 0.1 * (V + 40.0) / (1.0 - exp(-(V + 40.0) / 10.0));
}

double HHLookupTables::beta_m(double V) {
    return 4.0 * exp(-(V + 65.0) / 18.0);
}

double HHLookupTables::alpha_h(double V) {
    return 0.07 * exp(-(V + 65.0) / 20.0);
}

double HHLookupTables::beta_h(double V) {
    return 1.0 / (1.0 + exp(-(V + 35.0) / 10.0));
}

double HHLookupTables::alpha_n(double V) {
    // A small offset to avoid division by zero at V = -55
    if (std::abs(V + 55.0) < 1e-5) {
        return 0.1;
    }
    return 0.01 * (V + 55.0) / (1.0 - exp(-(V + 55.0) / 10.0));
}

double HHLookupTables::beta_n(double V) {
    return 0.125 * exp(-(V + 65.0) / 80.0);
}

// --- Public static method to initialize the tables ---

void HHLookupTables::initialize() {
    // Resize vectors to the correct size
    alpha_m_lut.resize(LUT_SIZE);
    beta_m_lut.resize(LUT_SIZE);
    alpha_h_lut.resize(LUT_SIZE);
    beta_h_lut.resize(LUT_SIZE);
    alpha_n_lut.resize(LUT_SIZE);
    beta_n_lut.resize(LUT_SIZE);

    // Populate the tables
    for (int i = 0; i < LUT_SIZE; ++i) {
        double v = V_MIN + i * V_STEP;
        alpha_m_lut[i] = alpha_m(v);
        beta_m_lut[i]  = beta_m(v);
        alpha_h_lut[i] = alpha_h(v);
        beta_h_lut[i]  = beta_h(v);
        alpha_n_lut[i] = alpha_n(v);
        beta_n_lut[i]  = beta_n(v);
    }
}
