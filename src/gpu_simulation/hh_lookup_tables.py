import torch
import numpy as np

class HHLookupTables:
    V_MIN = -100.0
    V_MAX = 50.0
    V_STEP = 0.01
    LUT_SIZE = int((V_MAX - V_MIN) / V_STEP) + 1

    def __init__(self, device='cpu', dtype=torch.float64):
        self.device = device
        self.dtype = dtype
        self.alpha_m_lut = None
        self.beta_m_lut = None
        self.alpha_h_lut = None
        self.beta_h_lut = None
        self.alpha_n_lut = None
        self.beta_n_lut = None
        self.initialize()

    def _alpha_m(self, V):
        V_shifted = V + 40.0
        # Use torch.where to handle the singularity
        return torch.where(torch.abs(V_shifted) < 1e-5,
                           torch.tensor(1.0, device=self.device),
                           0.1 * V_shifted / (1.0 - torch.exp(-V_shifted / 10.0)))

    def _beta_m(self, V):
        return 4.0 * torch.exp(-(V + 65.0) / 18.0)

    def _alpha_h(self, V):
        return 0.07 * torch.exp(-(V + 65.0) / 20.0)

    def _beta_h(self, V):
        return 1.0 / (1.0 + torch.exp(-(V + 35.0) / 10.0))

    def _alpha_n(self, V):
        V_shifted = V + 55.0
        # Use torch.where to handle the singularity
        return torch.where(torch.abs(V_shifted) < 1e-5,
                           torch.tensor(0.1, device=self.device),
                           0.01 * V_shifted / (1.0 - torch.exp(-V_shifted / 10.0)))

    def _beta_n(self, V):
        return 0.125 * torch.exp(-(V + 65.0) / 80.0)

    def initialize(self):
        v_range = torch.linspace(self.V_MIN, self.V_MAX, self.LUT_SIZE, device=self.device, dtype=self.dtype)

        self.alpha_m_lut = self._alpha_m(v_range)
        self.beta_m_lut = self._beta_m(v_range)
        self.alpha_h_lut = self._alpha_h(v_range)
        self.beta_h_lut = self._beta_h(v_range)
        self.alpha_n_lut = self._alpha_n(v_range)
        self.beta_n_lut = self._beta_n(v_range)

if __name__ == '__main__':
    # Example of how to use it
    luts = HHLookupTables(device='cuda' if torch.cuda.is_available() else 'cpu')
    print(f"LUTs created on device: {luts.device}")
    print(f"LUT size: {luts.LUT_SIZE}")
    print(f"Alpha_m LUT shape: {luts.alpha_m_lut.shape}")
    # Print first 5 values of alpha_m
    print("First 5 alpha_m values:")
    print(luts.alpha_m_lut[:5])
