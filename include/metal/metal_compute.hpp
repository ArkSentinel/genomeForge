#pragma once

#include <vector>
#include <cstdint>
#include <string>

namespace genome::metal {

struct MetalDeviceInfo {
    std::string device_name;
    bool has_unified_memory{true};
    uint64_t max_buffer_size{0};
};

struct GPUCellParticle {
    float x;
    float y;
    float vx;
    float vy;
    float radius;
    float pad[3];
};

class MetalComputeEngine {
public:
    MetalComputeEngine();
    ~MetalComputeEngine();

    bool initialize();
    bool is_available() const { return initialized_; }
    const MetalDeviceInfo& get_device_info() const { return device_info_; }

    // Ejecuta multiplicación matricial por lotes en GPU Metal usando memoria unificada (Zero-Copy)
    // X: [N, K], W: [M, K] -> Out: [N, M]
    bool run_batch_gemm(const float* X, int N, int K, const float* W, int M, float* Out, float* elapsed_us = nullptr);

    // Simula dinámica física espacial intercelular en placa de Petri en GPU Metal
    bool run_cell_physics(GPUCellParticle* particles, uint32_t count, float dt, float dish_radius, float* elapsed_us = nullptr);

private:
    struct Impl;
    Impl* pimpl_{nullptr};
    bool initialized_{false};
    MetalDeviceInfo device_info_;
};

} // namespace genome::metal
