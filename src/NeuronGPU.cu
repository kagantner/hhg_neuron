#include "NeuronGPU.h"
#include <cuda_runtime.h>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <vector>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// --- CUDA Error Handling Macro ---
#define CUDA_CHECK(err) { \
    cudaError_t e = (err); \
    if (e != cudaSuccess) { \
        throw std::runtime_error(std::string("CUDA Error: ") + cudaGetErrorString(e)); \
    } \
}

// --- Device-side (GPU) HH Model Functions ---

// Note: These functions are marked with __device__ to indicate they run on the GPU.
// They are direct translations of the HHModel class methods.

__device__ double alpha_m_gpu(double V) {
    double V_shifted = V + 40.0;
    if (abs(V_shifted) < 1e-5) return 1.0;
    return 0.1 * V_shifted / (1.0 - exp(-V_shifted / 10.0));
}

__device__ double beta_m_gpu(double V) {
    return 4.0 * exp(-(V + 65.0) / 18.0);
}

__device__ double alpha_h_gpu(double V) {
    return 0.07 * exp(-(V + 65.0) / 20.0);
}

__device__ double beta_h_gpu(double V) {
    return 1.0 / (1.0 + exp(-(V + 35.0) / 10.0));
}

__device__ double alpha_n_gpu(double V) {
    double V_shifted = V + 55.0;
    if (abs(V_shifted) < 1e-5) return 0.1;
    return 0.01 * V_shifted / (1.0 - exp(-V_shifted / 10.0));
}

__device__ double beta_n_gpu(double V) {
    return 0.125 * exp(-(V + 65.0) / 80.0);
}

// --- Main CUDA Kernel ---

__global__ void update_kernel(
    double* V_m, double* m, double* h, double* n,
    const double* I_inj,
    int n_seg, double dt, double g_a, double surface_area) {

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n_seg) return;

    // Biophysical constants (must be defined within device code)
    const double C_m_density = 1.0;
    const double g_Na_density = 120.0;
    const double g_K_density = 36.0;
    const double g_L_density = 0.3;
    const double E_Na = 50.0;
    const double E_K = -77.0;
    const double E_L = -54.387;

    // --- Step 1: Calculate axial current ---
    double v_prev = V_m[i]; // Store current voltage
    double v_left = (i > 0) ? V_m[i - 1] : v_prev;
    double v_right = (i < n_seg - 1) ? V_m[i + 1] : v_prev;
    double I_axial = g_a * (v_left - v_prev) + g_a * (v_right - v_prev);

    double total_current_uA = I_inj[i] + I_axial;
    double I_stim_density = (surface_area > 0) ? total_current_uA / surface_area : 0.0;

    // --- Step 2: Update membrane potential ---
    double I_Na = g_Na_density * m[i] * m[i] * m[i] * h[i] * (v_prev - E_Na);
    double I_K = g_K_density * n[i] * n[i] * n[i] * n[i] * (v_prev - E_K);
    double I_L = g_L_density * (v_prev - E_L);

    double dV = (I_stim_density - I_Na - I_K - I_L) / C_m_density;
    double new_V = v_prev + dV * dt;

    // --- Step 3: Update gating variables ---
    double am = alpha_m_gpu(new_V);
    double bm = beta_m_gpu(new_V);
    double ah = alpha_h_gpu(new_V);
    double bh = beta_h_gpu(new_V);
    double an = alpha_n_gpu(new_V);
    double bn = beta_n_gpu(new_V);

    m[i] = am / (am + bm) + (m[i] - am / (am + bm)) * exp(-dt * (am + bm));
    h[i] = ah / (ah + bh) + (h[i] - ah / (ah + bh)) * exp(-dt * (ah + bh));
    n[i] = an / (an + bn) + (n[i] - an / (an + bn)) * exp(-dt * (an + bn));

    // --- Step 4: Store new voltage ---
    V_m[i] = new_V;
}

// --- Host-side Class Implementation ---

// Definition of the struct holding GPU pointers
struct NeuronGPU_Data {
    double* d_V_m;
    double* d_m;
    double* d_h;
    double* d_n;
    double* d_I_inj;
    double g_a;
    double surface_area;
};

NeuronGPU::NeuronGPU(int num_segments, double length, double diameter, double Ra) {
    this->num_segments = num_segments;
    d_data = new NeuronGPU_Data();
    injected_currents_uA.resize(num_segments, 0.0);

    // Calculate geometry and conductance
    d_data->surface_area = M_PI * diameter * length;
    if (num_segments > 1) {
        double cross_area = M_PI * (diameter / 2.0) * (diameter / 2.0);
        double R_axial_kohm = (Ra * length / cross_area) / 1000.0;
        d_data->g_a = (R_axial_kohm > 0) ? 1.0 / R_axial_kohm : 0.0;
    } else {
        d_data->g_a = 0;
    }

    // Allocate GPU memory
    size_t size = num_segments * sizeof(double);
    CUDA_CHECK(cudaMalloc(&d_data->d_V_m, size));
    CUDA_CHECK(cudaMalloc(&d_data->d_m, size));
    CUDA_CHECK(cudaMalloc(&d_data->d_h, size));
    CUDA_CHECK(cudaMalloc(&d_data->d_n, size));
    CUDA_CHECK(cudaMalloc(&d_data->d_I_inj, size));

    // Initialize state on host
    std::vector<double> h_V_m(num_segments, -65.0);
    std::vector<double> h_m(num_segments);
    std::vector<double> h_h(num_segments);
    std::vector<double> h_n(num_segments);

    double V_init = -65.0;
    double m_init = alpha_m_gpu(V_init) / (alpha_m_gpu(V_init) + beta_m_gpu(V_init));
    double h_init = alpha_h_gpu(V_init) / (alpha_h_gpu(V_init) + beta_h_gpu(V_init));
    double n_init = alpha_n_gpu(V_init) / (alpha_n_gpu(V_init) + beta_n_gpu(V_init));

    std::fill(h_m.begin(), h_m.end(), m_init);
    std::fill(h_h.begin(), h_h.end(), h_init);
    std::fill(h_n.begin(), h_n.end(), n_init);

    // Copy initial state to GPU
    CUDA_CHECK(cudaMemcpy(d_data->d_V_m, h_V_m.data(), size, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_data->d_m, h_m.data(), size, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_data->d_h, h_h.data(), size, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_data->d_n, h_n.data(), size, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemset(d_data->d_I_inj, 0, size));
}

NeuronGPU::~NeuronGPU() {
    cudaFree(d_data->d_V_m);
    cudaFree(d_data->d_m);
    cudaFree(d_data->d_h);
    cudaFree(d_data->d_n);
    cudaFree(d_data->d_I_inj);
    delete d_data;
}

void NeuronGPU::update(double dt) {
    // Copy injected currents to GPU
    size_t size = num_segments * sizeof(double);
    CUDA_CHECK(cudaMemcpy(d_data->d_I_inj, injected_currents_uA.data(), size, cudaMemcpyHostToDevice));

    // Configure and launch kernel
    int threadsPerBlock = 256;
    int blocksPerGrid = (num_segments + threadsPerBlock - 1) / threadsPerBlock;

    update_kernel<<<blocksPerGrid, threadsPerBlock>>>(
        d_data->d_V_m, d_data->d_m, d_data->d_h, d_data->d_n,
        d_data->d_I_inj, num_segments, dt, d_data->g_a, d_data->surface_area
    );
    CUDA_CHECK(cudaGetLastError()); // Check for kernel launch errors
    CUDA_CHECK(cudaDeviceSynchronize()); // Wait for kernel to finish

    // Reset host-side injected currents
    std::fill(injected_currents_uA.begin(), injected_currents_uA.end(), 0.0);
}

void NeuronGPU::set_injected_current(int segment_index, double current_uA) {
    if (segment_index >= 0 && segment_index < num_segments) {
        injected_currents_uA[segment_index] = current_uA;
    } else {
        throw std::out_of_range("Segment index out of range.");
    }
}

double NeuronGPU::get_segment_V(int segment_index) const {
    if (segment_index < 0 || segment_index >= num_segments) {
        throw std::out_of_range("Segment index out of range.");
    }
    double V_val;
    // Copy single value from device to host
    CUDA_CHECK(cudaMemcpy(&V_val, d_data->d_V_m + segment_index, sizeof(double), cudaMemcpyDeviceToHost));
    return V_val;
}

std::vector<double> NeuronGPU::get_all_segment_V() const {
    std::vector<double> all_V(num_segments);
    CUDA_CHECK(cudaMemcpy(all_V.data(), d_data->d_V_m, num_segments * sizeof(double), cudaMemcpyDeviceToHost));
    return all_V;
}

int NeuronGPU::get_num_segments() const {
    return num_segments;
}