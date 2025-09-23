import torch
import time
from src.gpu_simulation.py_neuron import PyNeuron

def main():
    # --- 1. Simulation Parameters ---
    sim_time_ms = 100.0
    dt = 0.01
    num_steps = int(sim_time_ms / dt)

    # --- 2. Neuron Parameters ---
    num_segments = 1
    length_cm = 0.1
    diameter_cm = 0.02
    Ra_ohm_cm = 150.0

    # --- 3. Setup Device (GPU if available) ---
    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    print(f"Using device: {device}")

    # --- 4. Create Neuron ---
    neuron = PyNeuron(
        num_segments=num_segments,
        length_cm=length_cm,
        diameter_cm=diameter_cm,
        Ra_ohm_cm=Ra_ohm_cm,
        device=device
    )

    # --- 5. Simulation Setup ---
    # Prepare to store results
    v_history = []
    # Setup stimulus
    stimulus_current = 0.1 # uA
    stimulus_start_ms = 10.0
    stimulus_end_ms = 11.0

    # --- 6. Simulation Loop ---
    print("Starting simulation...")
    start_time = time.time()

    for i in range(num_steps):
        current_time_ms = i * dt

        # Apply stimulus
        injected_current = torch.zeros(num_segments, device=device, dtype=torch.float64)
        if stimulus_start_ms <= current_time_ms < stimulus_end_ms:
            # Stimulate the first segment
            injected_current[0] = stimulus_current

        # Update neuron state
        neuron.update(dt, injected_current)

        # Record membrane potential of first, middle, and last segments
        if i % 10 == 0: # Record every 0.1 ms
            v_history.append((
                current_time_ms,
                neuron.V_m[0].item(),
            ))

    end_time = time.time()
    print(f"Simulation finished in {end_time - start_time:.4f} seconds.")

    # --- 7. Print Results ---
    print("\nTime(ms)\tVm_first(mV)")
    for row in v_history[:20]: # Print first 20 recorded steps
        print(f"{row[0]:.2f}\t\t{row[1]:.4f}")

if __name__ == "__main__":
    main()
