#include <iostream>

// Forward declarations of our test functions
void test_lut_accuracy();
void test_simulation_accuracy();
void test_performance();

int main() {
    std::cout << "Test runner started..." << std::endl;
    std::cout << "------------------------" << std::endl;

    // Call all test suites
    test_lut_accuracy();
    test_simulation_accuracy();
    test_performance();

    std::cout << "------------------------" << std::endl;
    std::cout << "All tests passed." << std::endl;

    return 0;
}
