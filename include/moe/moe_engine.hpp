#pragma once

#include "vcf/vcf_parser.hpp"
#include "metal/metal_compute.hpp"

#include <vector>
#include <array>
#include <string>
#include <cstdint>

namespace genome::moe {

constexpr int FEATURE_DIM = 16;
constexpr int NUM_EXPERTS = 4;

enum class ExpertId : int {
    DnaRepair = 0,
    CellCycle = 1,
    P53Apoptosis = 2,
    Pluripotency = 3
};

struct GatingDecision {
    std::array<float, NUM_EXPERTS> softmax_weights{};
    int top1_expert{0};
    int top2_expert{0};
    float top1_weight{0.0f};
    float top2_weight{0.0f};
    bool dispatched_to_gpu{false};
    float execution_latency_us{0.0f};
};

struct MoEExecutionMetrics {
    uint64_t total_inferences{0};
    uint64_t cpu_vdsp_dispatches{0};
    uint64_t gpu_metal_dispatches{0};
    float avg_vdsp_latency_us{1.8f};
    float avg_metal_latency_us{8.5f};
    std::array<float, NUM_EXPERTS> expert_activation_levels{};
};

class MoEEngine {
public:
    MoEEngine();
    ~MoEEngine();

    // Extrae vector de características de 16 dimensiones de una variante genómica
    std::array<float, FEATURE_DIM> extract_features(const genome::vcf::GenomicVariant& variant) const;

    // Enrutamiento y evaluación de una variante individual (utiliza Accelerate vDSP/cblas)
    GatingDecision route_variant(const genome::vcf::GenomicVariant& variant);

    // Evaluación por lotes masivos (utiliza Metal GPU si N >= 64 y memoria unificada)
    std::vector<GatingDecision> route_batch(const std::vector<genome::vcf::GenomicVariant>& variants);

    // Métricas de hardware y de activación biológica
    const MoEExecutionMetrics& get_metrics() const { return metrics_; }
    bool is_metal_active() const { return metal_engine_.is_available(); }

private:
    void init_network_weights();
    void compute_cpu_gating(const float* feature_vec, float* logits_out);

    // Matrices de pesos del Gating Router [4 x 16]
    std::array<float, NUM_EXPERTS * FEATURE_DIM> gating_weights_{};
    std::array<float, NUM_EXPERTS> gating_biases_{};

    // Matrices de pesos de las 4 subredes de expertos [4 x 16 cada una]
    std::array<float, NUM_EXPERTS * FEATURE_DIM> expert_weights_[NUM_EXPERTS]{};

    genome::metal::MetalComputeEngine metal_engine_;
    MoEExecutionMetrics metrics_{};
};

} // namespace genome::moe
