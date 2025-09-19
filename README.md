# Hodgkin-Huxley Neuron Model Simulator

This project is a C++ implementation of a multi-compartment Hodgkin-Huxley neuron model. It simulates the propagation of action potentials along a neuron's axon, which is modeled as a series of connected cylindrical segments.

## Code Structure

The simulator is built around a few key classes:

*   `Neuron`: This class represents the entire neuron, composed of multiple segments. It manages the connections between segments, calculates the axial currents flowing between them, and orchestrates the overall simulation update loop.
*   `NeuronSegment`: This class represents a single, isopotential compartment of the neuron. It implements the core Hodgkin-Huxley equations to model the ionic currents (Na+, K+, leak) that determine the segment's membrane potential. For efficiency, it uses pre-calculated lookup tables for the voltage-dependent gating variable rate functions.
*   `HHLookupTables`: This is a helper class that pre-calculates and stores the alpha and beta rate functions for the Hodgkin-Huxley model's gating variables (m, h, n). This avoids costly calculations within the main simulation loop, significantly improving performance.
*   `main.cpp`: This is the entry point of the application. It sets up the simulation parameters (e.g., simulation time, time step), defines the neuron's geometry, creates the `Neuron` object, applies a stimulus current, and runs the simulation, printing the results to the console.

## Building the Code

The project uses a `Makefile` for easy compilation.

*   **To build the main simulation executable (`neuron_model`):**
    ```bash
    make
    ```
    This will compile all necessary source files and place the executable in the root directory.

*   **To build the test runner (`test_runner`):**
    ```bash
    make test
    ```

*   **To clean up all build artifacts:**
    ```bash
    make clean
    ```

## Running the Simulation

After building the project, you can run the simulation:

```bash
./neuron_model
```

The program will output a tab-separated table to the console. The columns represent:
1.  Time (in ms)
2.  Membrane potential of the first segment (in mV)
3.  Membrane potential of the middle segment (in mV)
4.  Membrane potential of the last segment (in mV)

You can redirect this output to a file for analysis, for example:
```bash
./neuron_model > simulation_output.tsv
```

## Running the Tests

The project includes a test suite to verify the correctness of the model. To run the tests:

```bash
make test && ./test_runner
```

This command will first build the `test_runner` executable (if it hasn't been built already) and then execute it. The results of the tests will be printed to the console.
