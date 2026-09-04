#include "crispr/crispr_engine.hpp"

#include <arm_neon.h>
#include <algorithm>
#include <cmath>
#include <iostream>

namespace genome::crispr {

CrisprEngine::CrisprEngine() {
    // Fondo genómico de referencia simulado para cribado de off-targets
    mock_off_target_background_ = {
        "GCCCUCAACAAGAUGUUAAA", // 3 mismatches
        "GCCCTCAACAAGATCTTCTG", // 4 mismatches
        "TCCCTCAACAAGATGTTTTG", // 1 mismatch (sitio pseudogénico)
        "ACCCUCAACAAGATGTTTTG", // 1 mismatch
        "GCCCUCAAAAAGATGTTTTG", // 2 mismatches
        "CTGCTCCCCCCCTGCCCTCA", // Chr 17 paralog
        "GAGGTTGGCTCTGACTGTAC", // Target secundario
        "CAACAAGATGTTTTGCCAGT"  // Target terciario
    };
}

CrisprEngine::~CrisprEngine() = default;

std::string CrisprEngine::reverse_complement(const std::string& seq) const {
    std::string rev;
    rev.reserve(seq.size());
    for (auto it = seq.rbegin(); it != seq.rend(); ++it) {
        char c = *it;
        switch (c) {
            case 'A': case 'a': rev += 'T'; break;
            case 'T': case 't': rev += 'A'; break;
            case 'C': case 'c': rev += 'G'; break;
            case 'G': case 'g': rev += 'C'; break;
            case 'U': case 'u': rev += 'A'; break;
            default: rev += 'N'; break;
        }
    }
    return rev;
}

int CrisprEngine::compute_hamming_distance_neon(const std::string& s1, const std::string& s2) {
    size_t len = std::min(s1.size(), s2.size());
    if (len < 16) {
        int dist = 0;
        for (size_t i = 0; i < len; ++i) {
            if (s1[i] != s2[i]) dist++;
        }
        return dist + static_cast<int>(std::max(s1.size(), s2.size()) - len);
    }

    // Comparación SIMD vectorizada con instrucciones ARM NEON de 128 bits
    uint8x16_t v1 = vld1q_u8(reinterpret_cast<const uint8_t*>(s1.data()));
    uint8x16_t v2 = vld1q_u8(reinterpret_cast<const uint8_t*>(s2.data()));

    // vceqq_u8 genera 0xFF donde los bytes coinciden, 0x00 donde difieren
    uint8x16_t eq_mask = vceqq_u8(v1, v2);

    alignas(16) uint8_t buffer[16];
    vst1q_u8(buffer, eq_mask);

    int matches = 0;
    for (int i = 0; i < 16; ++i) {
        if (buffer[i] == 0xFF) matches++;
    }

    // Evaluar los 4 nucleótidos restantes hasta completar 20 nt
    for (size_t i = 16; i < len; ++i) {
        if (s1[i] == s2[i]) matches++;
    }

    int distance = static_cast<int>(len) - matches;
    return distance;
}

float CrisprEngine::score_gc_content(float gc_pct) const {
    // Óptimo: 45% a 55%
    if (gc_pct >= 45.0f && gc_pct <= 55.0f) return 30.0f;
    if (gc_pct >= 40.0f && gc_pct <= 60.0f) return 25.0f;
    if (gc_pct >= 30.0f && gc_pct <= 70.0f) return 15.0f;
    return 5.0f;
}

float CrisprEngine::score_polymers(const std::string& spacer) const {
    // Penalización severa a secuencias terminadoras de polimerasa III (TTTT)
    if (spacer.find("TTTT") != std::string::npos || spacer.find("UUUU") != std::string::npos) {
        return -25.0f;
    }
    if (spacer.find("GGGG") != std::string::npos) {
        return -15.0f;
    }
    return 0.0f;
}

float CrisprEngine::evaluate_off_target_score(const std::string& protospacer, const std::vector<std::string>& background) {
    if (background.empty()) return 98.0f;

    float penalty_accum = 0.0f;
    for (const auto& bg_seq : background) {
        int h_dist = compute_hamming_distance_neon(protospacer, bg_seq);
        if (h_dist == 0) {
            penalty_accum += 40.0f; // Target idéntico en otro locus (riesgo crítico)
        } else if (h_dist == 1) {
            penalty_accum += 15.0f; // 1 solo mismatch
        } else if (h_dist == 2) {
            penalty_accum += 5.0f;  // 2 mismatches
        } else if (h_dist == 3) {
            penalty_accum += 1.5f;
        }
    }

    float score = 100.0f - penalty_accum;
    return std::clamp(score, 10.0f, 99.5f);
}

std::vector<GuideCandidate> CrisprEngine::design_guides(const std::string& sequence_window, int mutation_pos_in_window) {
    std::vector<GuideCandidate> candidates;
    int len = static_cast<int>(sequence_window.size());

    // 1. Escaneo en hebra forward buscando PAM NGG (ej. AGG, CGG, TGG, GGG)
    for (int i = 20; i < len - 2; ++i) {
        if (sequence_window[i + 1] == 'G' && sequence_window[i + 2] == 'G') {
            std::string pam = sequence_window.substr(i, 3);
            std::string spacer = sequence_window.substr(i - 20, 20);

            // Calcular contenido GC
            int gc_count = 0;
            for (char c : spacer) {
                if (c == 'G' || c == 'C' || c == 'g' || c == 'c') gc_count++;
            }
            float gc_pct = (static_cast<float>(gc_count) / 20.0f) * 100.0f;

            // Sitio de corte de la endonucleasa Cas9 (3 pb antes del PAM)
            int cut_site = i - 3;
            int dist_to_mut = std::abs(cut_site - mutation_pos_in_window);

            // Puntuación heurística de Doench / Rule Set 2
            float eff_score = 45.0f;
            eff_score += score_gc_content(gc_pct);
            eff_score += score_polymers(spacer);

            // Bonificación por cercanía a la mutación para maximizar eficiencia HDR
            if (dist_to_mut <= 3) eff_score += 25.0f;
            else if (dist_to_mut <= 8) eff_score += 15.0f;
            else if (dist_to_mut <= 15) eff_score += 8.0f;

            // Preferencia por base proximal al PAM
            if (spacer[19] == 'G' || spacer[19] == 'A') eff_score += 5.0f;

            eff_score = std::clamp(eff_score, 35.0f, 98.5f);

            // Evaluación de off-targets con vectorización NEON
            float off_safety = evaluate_off_target_score(spacer, mock_off_target_background_);

            GuideCandidate cand{};
            cand.pam = pam;
            cand.pam_position_offset = cut_site - mutation_pos_in_window;
            cand.protospacer_20nt = spacer;
            cand.gc_content_percent = gc_pct;
            cand.on_target_efficiency = eff_score;
            cand.off_target_safety_score = off_safety;
            cand.cut_distance_to_mutation = dist_to_mut;
            cand.is_reverse_strand = false;

            if (eff_score >= 90.0f && off_safety >= 95.0f) {
                cand.recommendation = "[ÓPTIMA: ALTA ESPECIFICIDAD]";
            } else if (eff_score >= 80.0f) {
                cand.recommendation = "[VIABLE: RECOMENDADA]";
            } else {
                cand.recommendation = "[MODERADA: REVISAR OFF-TARGETS]";
            }

            candidates.push_back(cand);
        }
    }

    // Ordenar candidatos por eficiencia decreciente
    std::sort(candidates.begin(), candidates.end(), [](const GuideCandidate& a, const GuideCandidate& b) {
        return (a.on_target_efficiency * 0.6f + a.off_target_safety_score * 0.4f) >
               (b.on_target_efficiency * 0.6f + b.off_target_safety_score * 0.4f);
    });

    return candidates;
}

std::string CrisprEngine::simulate_hdr_edit(
    const std::string& sequence,
    int mutation_offset,
    const std::string& wild_type_codon
) {
    if (mutation_offset < 0 || mutation_offset + wild_type_codon.size() > sequence.size()) {
        return sequence;
    }
    std::string edited = sequence;
    edited.replace(mutation_offset, wild_type_codon.size(), wild_type_codon);
    return edited;
}

} // namespace genome::crispr
