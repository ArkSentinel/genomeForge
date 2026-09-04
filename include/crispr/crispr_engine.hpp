#pragma once

#include <string>
#include <vector>
#include <cstdint>

namespace genome::crispr {

struct GuideCandidate {
    std::string pam;                  // Motivo PAM detectado (ej. "CGG", "TGG")
    int pam_position_offset;          // Posición relativa al locus de la mutación
    std::string protospacer_20nt;     // Secuencia espaciadora de 20 nucleótidos
    float gc_content_percent;         // Porcentaje de contenido GC (óptimo: 40-60%)
    float on_target_efficiency;       // Puntuación de corte on-target [0.0 - 100.0%]
    float off_target_safety_score;    // Puntuación de especificidad frente a off-targets [0.0 - 100.0]
    int cut_distance_to_mutation;     // Distancia en pb desde el corte DSB hasta la mutación
    bool is_reverse_strand;           // Hebra inversa (PAM CCN)
    std::string recommendation;       // "[ÓPTIMA]", "[VIABLE]", "[BAJA ESPECIFICIDAD]"
};

class CrisprEngine {
public:
    CrisprEngine();
    ~CrisprEngine();

    // Escanea una ventana genómica en busca de motivos PAM (NGG y CCN) y diseña guías de 20 nt
    std::vector<GuideCandidate> design_guides(const std::string& sequence_window, int mutation_pos_in_window);

    // Algoritmo vectorizado (ARM NEON) análogo a Cas-OFFinder para calcular distancia de Hamming
    int compute_hamming_distance_neon(const std::string& s1, const std::string& s2);

    // Evalúa la especificidad de una guía frente a una lista de secuencias genómicas homólogas
    float evaluate_off_target_score(const std::string& protospacer, const std::vector<std::string>& reference_background);

    // Simula la edición por HDR (Homology-Directed Repair) sustituyendo el codón diana
    std::string simulate_hdr_edit(const std::string& sequence, int mutation_offset, const std::string& wild_type_codon);

private:
    float score_gc_content(float gc_pct) const;
    float score_polymers(const std::string& spacer) const;
    std::string reverse_complement(const std::string& seq) const;

    std::vector<std::string> mock_off_target_background_;
};

} // namespace genome::crispr
