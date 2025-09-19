#ifndef HH_LOOKUP_TABLES_H
#define HH_LOOKUP_TABLES_H

#include <vector>

class HHLookupTables {
public:
    // Initializes the static lookup tables. Must be called once before simulation.
    static void initialize();

    // Lookup table parameters
    static constexpr double V_MIN = -100.0; // mV
    static constexpr double V_MAX = 50.0;   // mV
    static constexpr double V_STEP = 0.01;  // mV resolution
    static constexpr int LUT_SIZE = static_cast<int>((V_MAX - V_MIN) / V_STEP) + 1;

    // Static lookup tables, accessible to NeuronSegment
    static std::vector<double> alpha_m_lut;
    static std::vector<double> beta_m_lut;
    static std::vector<double> alpha_h_lut;
    static std::vector<double> beta_h_lut;
    static std::vector<double> alpha_n_lut;
    static std::vector<double> beta_n_lut;

private:
    // The original kinetic functions, now private and used only for table generation.
    static double alpha_m(double V);
    static double beta_m(double V);
    static double alpha_h(double V);
    static double beta_h(double V);
    static double alpha_n(double V);
    static double beta_n(double V);
};

#endif // HH_LOOKUP_TABLES_H
