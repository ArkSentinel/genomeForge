#include "moe/moe_engine.hpp"

#include <Accelerate/Accelerate.h>
#include <cmath>
#include <chrono>
#include <algorithm>
#include <cstring>
#include <iostream>

namespace genome::moe {

MoEEngine::MoEEngine() {
    init_network_weights();
}

MoEEngine::~MoEEngine() = default;

void MoEEngine::init_network_weights() {
    // Inicialización heurística orientada a vías biológicas para los 4 expertos:
    // Expert 0: Daño ADN / Reparación (reacciona a variantes patogénicas y de alta profundidad)
    // Expert 1: Ciclo Celular / Senescencia (reacciona a genes supresores como CDKN2A)
    // Expert 2: Apoptosis p53 (reacciona a TP53 y alta patogenicidad)
    // Expert 3: Pluripotencia (reacciona a MYC y factores de transcripción)

    // Pesos del Gating Router [4 x 16]
    for (int e = 0; e < NUM_EXPERTS; ++e) {
        gating_biases_[e] = 0.1f * static_cast<float>(e);
        for (int f = 0; f < FEATURE_DIM; ++f) {
            float val = 0.05f * static_cast<float>((e * 3 + f * 7) % 11 - 5);
            gating_weights_[e * FEATURE_DIM + f] = val;
        }
    }

    // Calibrar pesos específicos del router para genes diana:
    // f=5: marcador TP53, f=6: marcador CDKN2A, f=7: marcador MYC, f=2: pathogenicity_score
    gating_weights_[static_cast<int>(ExpertId::P53Apoptosis) * FEATURE_DIM + 2] = 2.5f; // Pathogenic
    gating_weights_[static_cast<int>(ExpertId::P53Apoptosis) * FEATURE_DIM + 5] = 4.0f; // TP53 gene

    gating_weights_[static_cast<int>(ExpertId::CellCycle) * FEATURE_DIM + 6] = 3.5f;    // CDKN2A gene
    gating_weights_[static_cast<int>(ExpertId::CellCycle) * FEATURE_DIM + 2] = 1.8f;    // Pathogenic

    gating_weights_[static_cast<int>(ExpertId::Pluripotency) * FEATURE_DIM + 7] = 3.8f; // MYC gene
    // Inicializar matrices de pesos para cada subred de expertos [4 x 16]
    for (int e = 0; e < NUM_EXPERTS; ++e) {
        for (int f = 0; f < FEATURE_DIM; ++f) {
            expert_weights_[e][f] = 0.05f * static_cast<float>((e * 4 + f * 3) % 7 + 1);
        }
    }

    metrics_.expert_activation_levels = {0.75f, 0.60f, 0.15f, 0.82f};
}

std::array<float, FEATURE_DIM> MoEEngine::extract_features(const genome::vcf::GenomicVariant& v) const {
    std::array<float, FEATURE_DIM> feat{};

    feat[0] = static_cast<float>(v.pos) / 250000000.0f;           // Coordenada genómica normalizada
    feat[1] = static_cast<float>(v.chrom_id) / 24.0f;             // Cromosoma normalizado (1..24)
    feat[2] = v.pathogenicity_score;                              // Puntuación patogénica [0.0 - 1.0]
    feat[3] = v.allele_frequency;                                 // Frecuencia alélica (AF)
    feat[4] = std::min(static_cast<float>(v.depth) / 150.0f, 1.0f); // Cobertura (DP) normalizada

    std::string_view gene(v.gene_symbol);
    feat[5] = (gene == "TP53") ? 1.0f : 0.0f;
    feat[6] = (gene == "CDKN2A") ? 1.0f : 0.0f;
    feat[7] = (gene == "MYC") ? 1.0f : 0.0f;
    feat[8] = (gene == "PTEN") ? 1.0f : 0.0f;
    feat[9] = (gene == "KRAS" || gene == "EGFR") ? 1.0f : 0.0f;

    feat[10] = (v.effect == genome::vcf::VariantEffect::YamanakaCollision) ? 1.0f : 0.0f;
    feat[11] = (v.effect == genome::vcf::VariantEffect::Pathogenic) ? 1.0f : 0.0f;
    feat[12] = (v.effect == genome::vcf::VariantEffect::Benign) ? 1.0f : 0.0f;

    // Razón de longitudes de alelos (InDel vs SNV)
    feat[13] = (std::strlen(v.ref) == 1 && std::strlen(v.alt) == 1) ? 0.0f : 1.0f;
    feat[14] = static_cast<float>(v.flags & 0xFF) / 255.0f;
    feat[15] = 1.0f; // Bias término de entrada

    return feat;
}

void MoEEngine::compute_cpu_gating(const float* feature_vec, float* logits_out) {
    // Uso del framework Apple Accelerate (BLAS cblas_sgemv optimizado para núcleos ARM64/AMX)
    cblas_sgemv(
        CblasRowMajor,
        CblasNoTrans,
        NUM_EXPERTS,
        FEATURE_DIM,
        1.0f,
        gating_weights_.data(),
        FEATURE_DIM,
        feature_vec,
        1,
        0.0f,
        logits_out,
        1
    );

    // Sumar sesgos
    for (int i = 0; i < NUM_EXPERTS; ++i) {
        logits_out[i] += gating_biases_[i];
    }
}

GatingDecision MoEEngine::route_variant(const genome::vcf::GenomicVariant& variant) {
    auto start_time = std::chrono::high_resolution_clock::now();

    auto features = extract_features(variant);
    float logits[NUM_EXPERTS];
    compute_cpu_gating(features.data(), logits);

    // Aplicar Softmax
    float max_l = *std::max_element(logits, logits + NUM_EXPERTS);
    float sum_exp = 0.0f;
    float exps[NUM_EXPERTS];
    for (int i = 0; i < NUM_EXPERTS; ++i) {
        exps[i] = std::exp(logits[i] - max_l);
        sum_exp += exps[i];
    }

    GatingDecision decision{};
    for (int i = 0; i < NUM_EXPERTS; ++i) {
        decision.softmax_weights[i] = exps[i] / (sum_exp + 1e-6f);
    }

    // Selección Top-2 Experts
    int best1 = 0, best2 = 1;
    if (decision.softmax_weights[best2] > decision.softmax_weights[best1]) {
        std::swap(best1, best2);
    }
    for (int i = 2; i < NUM_EXPERTS; ++i) {
        if (decision.softmax_weights[i] > decision.softmax_weights[best1]) {
            best2 = best1;
            best1 = i;
        } else if (decision.softmax_weights[i] > decision.softmax_weights[best2]) {
            best2 = i;
        }
    }

    decision.top1_expert = best1;
    decision.top2_expert = best2;
    decision.top1_weight = decision.softmax_weights[best1];
    decision.top2_weight = decision.softmax_weights[best2];
    decision.dispatched_to_gpu = false;

    // Computar activación de las subredes de expertos Top-2 seleccionadas
    float expert_response1 = 0.0f;
    float expert_response2 = 0.0f;
    for (int f = 0; f < FEATURE_DIM; ++f) {
        expert_response1 += expert_weights_[best1][f] * features[f];
        expert_response2 += expert_weights_[best2][f] * features[f];
    }
    (void)expert_response1;
    (void)expert_response2;

    auto end_time = std::chrono::high_resolution_clock::now();
    decision.execution_latency_us = std::chrono::duration<float, std::micro>(end_time - start_time).count();

    // Actualizar métricas acumuladas
    metrics_.total_inferences++;
    metrics_.cpu_vdsp_dispatches++;
    metrics_.avg_vdsp_latency_us = 0.95f * metrics_.avg_vdsp_latency_us + 0.05f * decision.execution_latency_us;

    for (int i = 0; i < NUM_EXPERTS; ++i) {
        metrics_.expert_activation_levels[i] = 
            0.90f * metrics_.expert_activation_levels[i] + 0.10f * decision.softmax_weights[i];
    }

    return decision;
}

std::vector<GatingDecision> MoEEngine::route_batch(const std::vector<genome::vcf::GenomicVariant>& variants) {
    size_t N = variants.size();
    std::vector<GatingDecision> decisions(N);
    if (N == 0) return decisions;

    // Aplanar matriz de características [N x FEATURE_DIM]
    std::vector<float> X(N * FEATURE_DIM);
    for (size_t i = 0; i < N; ++i) {
        auto feat = extract_features(variants[i]);
        std::memcpy(&X[i * FEATURE_DIM], feat.data(), FEATURE_DIM * sizeof(float));
    }

    std::vector<float> logits(N * NUM_EXPERTS);
    bool use_metal = (N >= 64 && metal_engine_.is_available());
    float elapsed_us = 0.0f;

    if (use_metal) {
        // Despachar a GPU Metal mediante buffers en memoria unificada (Zero-Copy)
        metal_engine_.run_batch_gemm(
            X.data(), static_cast<int>(N), FEATURE_DIM,
            gating_weights_.data(), NUM_EXPERTS,
            logits.data(),
            &elapsed_us
        );
        metrics_.gpu_metal_dispatches++;
        metrics_.avg_metal_latency_us = 0.90f * metrics_.avg_metal_latency_us + 0.10f * (elapsed_us / N);
    } else {
        auto start = std::chrono::high_resolution_clock::now();
        // Cómputo CPU vectorizado con Accelerate cblas_sgemm
        cblas_sgemm(
            CblasRowMajor, CblasNoTrans, CblasTrans,
            static_cast<int>(N), NUM_EXPERTS, FEATURE_DIM,
            1.0f,
            X.data(), FEATURE_DIM,
            gating_weights_.data(), FEATURE_DIM,
            0.0f,
            logits.data(), NUM_EXPERTS
        );
        auto end = std::chrono::high_resolution_clock::now();
        elapsed_us = std::chrono::duration<float, std::micro>(end - start).count();
        metrics_.cpu_vdsp_dispatches++;
        metrics_.avg_vdsp_latency_us = 0.90f * metrics_.avg_vdsp_latency_us + 0.10f * (elapsed_us / N);
    }

    // Evaluar Softmax y Top-2 por fila
    for (size_t row = 0; row < N; ++row) {
        float* row_logits = &logits[row * NUM_EXPERTS];
        for (int e = 0; e < NUM_EXPERTS; ++e) {
            row_logits[e] += gating_biases_[e];
        }

        float max_l = *std::max_element(row_logits, row_logits + NUM_EXPERTS);
        float sum_exp = 0.0f;
        float exps[NUM_EXPERTS];
        for (int e = 0; e < NUM_EXPERTS; ++e) {
            exps[e] = std::exp(row_logits[e] - max_l);
            sum_exp += exps[e];
        }

        auto& dec = decisions[row];
        for (int e = 0; e < NUM_EXPERTS; ++e) {
            dec.softmax_weights[e] = exps[e] / (sum_exp + 1e-6f);
        }

        int b1 = 0, b2 = 1;
        if (dec.softmax_weights[b2] > dec.softmax_weights[b1]) std::swap(b1, b2);
        for (int e = 2; e < NUM_EXPERTS; ++e) {
            if (dec.softmax_weights[e] > dec.softmax_weights[b1]) {
                b2 = b1; b1 = e;
            } else if (dec.softmax_weights[e] > dec.softmax_weights[b2]) {
                b2 = e;
            }
        }

        dec.top1_expert = b1;
        dec.top2_expert = b2;
        dec.top1_weight = dec.softmax_weights[b1];
        dec.top2_weight = dec.softmax_weights[b2];
        dec.dispatched_to_gpu = use_metal;
        dec.execution_latency_us = elapsed_us / static_cast<float>(N);
    }

    metrics_.total_inferences += N;
    return decisions;
}

} // namespace genome::moe
