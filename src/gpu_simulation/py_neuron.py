import torch
import numpy as np
from .hh_lookup_tables import HHLookupTables

class PyNeuron:
    # --- Constants from NeuronSegment.h ---
    E_Na = 50.0   # Sodium reversal potential (mV)
    E_K = -77.0   # Potassium reversal potential (mV)
    E_L = -54.4  # Leak reversal potential (mV)
    g_Na_density = 120.0 # Sodium conductance density (mS/cm^2)
    g_K_density = 36.0   # Potassium conductance density (mS/cm^2)
    g_L_density = 0.3    # Leak conductance density (mS/cm^2)
    C_m_density = 1.0    # Membrane capacitance (uF/cm^2)

    def __init__(self, num_segments, length_cm, diameter_cm, Ra_ohm_cm, device='cpu'):
        self.device = device
        self.num_segments = num_segments

        # --- Geometric Properties ---
        self.surface_area = torch.tensor(np.pi * diameter_cm * length_cm, device=self.device, dtype=torch.float64)
        if num_segments > 1:
            cross_area = np.pi * (diameter_cm / 2.0)**2
            R_axial_kohm = (Ra_ohm_cm * length_cm / cross_area) / 1000.0
            self.g_a = 1.0 / R_axial_kohm
        else:
            self.g_a = 0

        # --- Lookup Tables ---
        self.luts = HHLookupTables(device=self.device)

        # --- State Variables (as tensors) ---
        self.V_m = torch.full((num_segments,), -65.0, device=self.device, dtype=torch.float64)
        self.m = torch.full((num_segments,), 0.0, device=self.device, dtype=torch.float64)
        self.h = torch.full((num_segments,), 0.0, device=self.device, dtype=torch.float64)
        self.n = torch.full((num_segments,), 0.0, device=self.device, dtype=torch.float64)
        self._initialize_gating_variables()

    def _initialize_gating_variables(self):
        # Initialize to resting state, similar to C++ constructor
        V_rest = -65.0
        am = self.luts._alpha_m(torch.tensor(V_rest, device=self.device))
        bm = self.luts._beta_m(torch.tensor(V_rest, device=self.device))
        ah = self.luts._alpha_h(torch.tensor(V_rest, device=self.device))
        bh = self.luts._beta_h(torch.tensor(V_rest, device=self.device))
        an = self.luts._alpha_n(torch.tensor(V_rest, device=self.device))
        bn = self.luts._beta_n(torch.tensor(V_rest, device=self.device))

        self.m.fill_(am / (am + bm))
        self.h.fill_(ah / (ah + bh))
        self.n.fill_(an / (an + bn))

    def _fast_exp_approx(self, x):
        # Padé approximant for exp(-x)
        return (1.0 - 0.5 * x) / (1.0 + 0.5 * x)

    def update(self, dt, I_injected):
        # Ensure I_injected is a tensor on the correct device
        if not isinstance(I_injected, torch.Tensor):
            I_injected = torch.tensor(I_injected, device=self.device, dtype=torch.float64)
        if I_injected.ndim == 0:
            I_injected = I_injected.expand(self.num_segments)

        V_at_start_of_step = self.V_m

        # --- Step 1: Update gating variables (Rush-Larsen) ---
        # Interpolate alpha and beta values from LUTs
        pos = (V_at_start_of_step - self.luts.V_MIN) / self.luts.V_STEP
        idx_base = pos.long()
        frac = (pos - idx_base).unsqueeze(-1) # Add a dimension for broadcasting

        # Clamp indices to be within LUT bounds
        idx = torch.clamp(idx_base, 0, self.luts.LUT_SIZE - 2)

        # Gather from LUTs
        alpha_m = self.luts.alpha_m_lut[idx]
        alpha_m_plus_1 = self.luts.alpha_m_lut[idx+1]
        am = alpha_m + frac.squeeze(-1) * (alpha_m_plus_1 - alpha_m)

        beta_m = self.luts.beta_m_lut[idx]
        beta_m_plus_1 = self.luts.beta_m_lut[idx+1]
        bm = beta_m + frac.squeeze(-1) * (beta_m_plus_1 - beta_m)

        alpha_h = self.luts.alpha_h_lut[idx]
        alpha_h_plus_1 = self.luts.alpha_h_lut[idx+1]
        ah = alpha_h + frac.squeeze(-1) * (alpha_h_plus_1 - alpha_h)

        beta_h = self.luts.beta_h_lut[idx]
        beta_h_plus_1 = self.luts.beta_h_lut[idx+1]
        bh = beta_h + frac.squeeze(-1) * (beta_h_plus_1 - beta_h)

        alpha_n = self.luts.alpha_n_lut[idx]
        alpha_n_plus_1 = self.luts.alpha_n_lut[idx+1]
        an = alpha_n + frac.squeeze(-1) * (alpha_n_plus_1 - alpha_n)

        beta_n = self.luts.beta_n_lut[idx]
        beta_n_plus_1 = self.luts.beta_n_lut[idx+1]
        bn = beta_n + frac.squeeze(-1) * (beta_n_plus_1 - beta_n)

        # Update m, h, n
        sum_m = am + bm
        m_inf = am / sum_m
        self.m = m_inf + (self.m - m_inf) * self._fast_exp_approx(dt * sum_m)

        sum_h = ah + bh
        h_inf = ah / sum_h
        self.h = h_inf + (self.h - h_inf) * self._fast_exp_approx(dt * sum_h)

        sum_n = an + bn
        n_inf = an / sum_n
        self.n = n_inf + (self.n - n_inf) * self._fast_exp_approx(dt * sum_n)

        # --- Step 2: Calculate currents ---
        # Ionic currents
        I_Na = self.g_Na_density * self.surface_area * self.m**3 * self.h * (V_at_start_of_step - self.E_Na)
        I_K  = self.g_K_density  * self.surface_area * self.n**4 * (V_at_start_of_step - self.E_K)
        I_L  = self.g_L_density  * self.surface_area * (V_at_start_of_step - self.E_L)

        # Axial currents
        I_axial = torch.zeros_like(self.V_m)
        if self.num_segments > 1:
            v_left = torch.roll(V_at_start_of_step, 1)
            v_left[0] = V_at_start_of_step[0] # No-flux boundary

            v_right = torch.roll(V_at_start_of_step, -1)
            v_right[-1] = V_at_start_of_step[-1] # No-flux boundary

            I_axial = self.g_a * (v_left - V_at_start_of_step) + self.g_a * (v_right - V_at_start_of_step)

        I_total = I_injected + I_axial - (I_Na + I_K + I_L)

        # --- Step 3: Update membrane potential (Forward Euler) ---
        C_total = self.C_m_density * self.surface_area
        dV = I_total / C_total
        self.V_m += dV * dt
