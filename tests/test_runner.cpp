#include <iostream>

// Forward declarations of our test functions
void test_lut_accuracy();
void test_simulation_accuracy();
void test_performance();
void test_lut_synapse();

int main() {
    std::cout << "Test runner started..." << std::endl;
    std::cout << "------------------------" << std::endl;

    // Call all test suites
    try {
        test_lut_accuracy();
        test_simulation_accuracy();
        test_performance();
        test_lut_synapse();

        std::cout << "------------------------" << std::endl;
        std::cout << "All tests passed." << std::endl;
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "A test failed: " << e.what() << std::endl;
        return 1;
    }
}
