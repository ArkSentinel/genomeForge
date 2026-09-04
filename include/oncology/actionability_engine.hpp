#pragma once

#include "vcf/vcf_parser.hpp"

#include <string>
#include <vector>
#include <cstdint>

namespace genome::oncology {

enum class EvidenceLevel : uint8_t {
    Level1_FDA_Approved = 1,   // FDA/EMA approved standard-of-care biomarker
    Level2_StandardOfCare = 2, // Supported by NCCN / ESMO clinical guidelines
    Level3_ClinicalTrials = 3, // Compelling biological rationale; active Phase II/III trials
    Level4_Preclinical = 4,    // In-silico or preclinical laboratory model evidence
    Resistance = 5             // Biomarker conferring secondary or primary resistance
};

struct TargetedDrug {
    std::string drug_name;
    std::string drug_class;
    std::string target_gene;
    std::string indication;
    std::string mechanism_of_action;
    std::string clinical_trial_or_approval;
    EvidenceLevel evidence_level{EvidenceLevel::Level1_FDA_Approved};
    float synthetic_lethality_score{0.0f};  // 0.0 to 1.0
    float tumor_regression_efficacy{0.85f}; // In-silico cell death rate (0.0 to 1.0)
    std::string recommended_dose{"Standard dosing regimen"};
};

struct ActionableBiomarker {
    std::string gene_symbol;
    std::string variant_desc;
    uint32_t genomic_pos{0};
    std::string alt_allele;
    std::string consequence;
    bool is_resistance_marker{false};
    float pathogenicity_score{0.0f};
    std::vector<TargetedDrug> matched_drugs;
};

struct ActionabilityAnalysis {
    std::vector<ActionableBiomarker> detected_biomarkers;
    std::vector<TargetedDrug> all_recommended_drugs;
    std::vector<std::string> resistance_warnings;
    size_t level1_count{0};
    size_t level2_count{0};
    size_t level3_count{0};
    float total_actionability_score{0.0f};
    bool high_tmb{false};
    float tumor_mutational_burden{0.0f}; // mut/Mb estimate
};

class ActionabilityEngine {
public:
    ActionabilityEngine();
    ~ActionabilityEngine();

    // Analiza las variantes detectadas en el VCF y genera recomendaciones terapéuticas
    ActionabilityAnalysis analyze_variants(const std::vector<vcf::GenomicVariant>& variants);

    // Formateo de etiquetas HUD
    static std::string evidence_level_to_string(EvidenceLevel level);
    static std::string evidence_level_badge(EvidenceLevel level);

private:
    std::vector<TargetedDrug> drug_database_;
    void init_drug_database();
};

} // namespace genome::oncology
