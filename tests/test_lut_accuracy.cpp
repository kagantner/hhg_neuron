#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>
#include "HHLookupTables.h"

// To avoid making the golden model's methods public, we can
// re-declare the pure alpha/beta functions here for the test.
namespace GoldenFunctions {
    double alpha_m(double V) {
        if (std::abs(V + 40.0) < 1e-5) return 1.0;
        return 0.1 * (V + 40.0) / (1.0 - exp(-(V + 40.0) / 10.0));
    }
    double beta_m(double V) {
        return 4.0 * exp(-(V + 65.0) / 18.0);
    }
    double alpha_h(double V) {
        return 0.07 * exp(-(V + 65.0) / 20.0);
    }
    double beta_h(double V) {
        return 1.0 / (1.0 + exp(-(V + 35.0) / 10.0));
    }
    double alpha_n(double V) {
        if (std::abs(V + 55.0) < 1e-5) return 0.1;
        return 0.01 * (V + 55.0) / (1.0 - exp(-(V + 55.0) / 10.0));
    }
    double beta_n(double V) {
        return 0.125 * exp(-(V + 65.0) / 80.0);
    }
}

void test_lut_accuracy() {
    std::cout << "Running LUT Accuracy Test..." << std::endl;

    // Ensure the tables are initialized
    HHLookupTables::initialize();

    const double tolerance = 1e-9;

    for (int i = 0; i < HHLookupTables::LUT_SIZE; ++i) {
        double v = HHLookupTables::V_MIN + i * HHLookupTables::V_STEP;

        // Test alpha_m
        double golden_am = GoldenFunctions::alpha_m(v);
        double lut_am = HHLookupTables::alpha_m_lut[i];
        assert(std::abs(golden_am - lut_am) < tolerance);

        // Test beta_m
        double golden_bm = GoldenFunctions::beta_m(v);
        double lut_bm = HHLookupTables::beta_m_lut[i];
        assert(std::abs(golden_bm - lut_bm) < tolerance);

        // Test alpha_h
        double golden_ah = GoldenFunctions::alpha_h(v);
        double lut_ah = HHLookupTables::alpha_h_lut[i];
        assert(std::abs(golden_ah - lut_ah) < tolerance);

        // Test beta_h
        double golden_bh = GoldenFunctions::beta_h(v);
        double lut_bh = HHLookupTables::beta_h_lut[i];
        assert(std::abs(golden_bh - lut_bh) < tolerance);

        // Test alpha_n
        double golden_an = GoldenFunctions::alpha_n(v);
        double lut_an = HHLookupTables::alpha_n_lut[i];
        assert(std::abs(golden_an - lut_an) < tolerance);

        // Test beta_n
        double golden_bn = GoldenFunctions::beta_n(v);
        double lut_bn = HHLookupTables::beta_n_lut[i];
        assert(std::abs(golden_bn - lut_bn) < tolerance);
    }

    std::cout << "  [PASS] LUT values match golden function calculations." << std::endl;
}
