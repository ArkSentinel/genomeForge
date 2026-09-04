#pragma once

#include <string>
#include <vector>
#include <cstdint>

namespace genome::system {

struct CpuCoreMetrics {
    int core_id{0};
    bool is_performance_core{true};
    float usage_percent{0.0f};
};

struct HardwareTelemetry {
    // CPU Metrics
    float total_cpu_usage{0.0f};
    float p_core_avg_usage{0.0f};
    float e_core_avg_usage{0.0f};
    int p_core_count{0};
    int e_core_count{0};
    std::vector<CpuCoreMetrics> core_metrics;

    // Memory Metrics (Unified Memory Architecture)
    uint64_t total_memory_bytes{0};
    uint64_t used_memory_bytes{0};
    uint64_t free_memory_bytes{0};
    uint64_t wired_memory_bytes{0};
    float memory_used_percent{0.0f};

    // Platform Info
    std::string chip_model{"Apple Silicon"};
    std::string os_version{""};
};

class TelemetryCollector {
public:
    TelemetryCollector();
    ~TelemetryCollector();

    // Actualiza las métricas instantáneas calculando diferencias desde la última muestra
    void update();

    // Obtiene la instantánea actual
    const HardwareTelemetry& get_metrics() const { return current_metrics_; }

private:
    void init_platform_info();
    void sample_cpu();
    void sample_memory();

    HardwareTelemetry current_metrics_;

    // Estados previos de CPU ticks para calcular delta
    uint32_t num_cpu_cores_{0};
    std::vector<uint64_t> prev_cpu_ticks_;
    bool has_prev_ticks_{false};
};

} // namespace genome::system
