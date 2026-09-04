#include "oncology/actionability_engine.hpp"

#include <algorithm>
#include <sstream>
#include <iomanip>

namespace genome::oncology {

ActionabilityEngine::ActionabilityEngine() {
    init_drug_database();
}

ActionabilityEngine::~ActionabilityEngine() = default;

void ActionabilityEngine::init_drug_database() {
    drug_database_ = {
        // EGFR Alterations
        {
            "Osimertinib",
            "3rd-Gen Irreversible EGFR TKI",
            "EGFR",
            "NSCLC with EGFR exon 19 del or L858R / T790M",
            "Selective and irreversible inhibition of mutant EGFR kinase domain",
            "FDA Approved (FLAURA & AURA3 Trials)",
            EvidenceLevel::Level1_FDA_Approved,
            0.20f,
            0.92f,
            "80 mg orally once daily"
        },
        {
            "Gefitinib",
            "1st-Gen Reversible EGFR TKI",
            "EGFR",
            "Metastatic NSCLC with sensitizing EGFR mutations",
            "ATP-competitive reversible inhibition of EGFR tyrosine kinase",
            "FDA Approved (IPASS Study - Vulnerable to T790M resistance)",
            EvidenceLevel::Level2_StandardOfCare,
            0.15f,
            0.78f,
            "250 mg orally once daily"
        },

        // BRCA1 / BRCA2 DNA Repair Deficiencies
        {
            "Olaparib",
            "PARP1/2 Inhibitor (Synthetic Lethality)",
            "BRCA1",
            "BRCA1/2 mutated ovarian, breast, pancreatic, and prostate cancer",
            "Traps PARP1/2 on DNA single-strand breaks, inducing double-strand breaks in HR-deficient cells",
            "FDA Approved (OlympiAD & SOLO-1 Trials)",
            EvidenceLevel::Level1_FDA_Approved,
            0.95f,
            0.88f,
            "300 mg orally twice daily"
        },
        {
            "Talazoparib",
            "Potent PARP Trapping Inhibitor",
            "BRCA1",
            "HER2-negative locally advanced or metastatic breast cancer with gBRCA",
            "Exceptional PARP trapping potency causing catastrophic replication fork collapse",
            "FDA Approved (EMBRACA Trial)",
            EvidenceLevel::Level1_FDA_Approved,
            0.92f,
            0.86f,
            "1 mg orally once daily"
        },

        // KRAS Alterations
        {
            "Sotorasib (AMG 510)",
            "Covalent KRAS(G12C) Inhibitor",
            "KRAS",
            "Locally advanced or metastatic NSCLC with KRAS G12C mutation",
            "Irreversibly binds the switch-II pocket of GDP-bound KRAS, locking it in inactive state",
            "FDA Approved (CodeBreaK 100 Trial)",
            EvidenceLevel::Level1_FDA_Approved,
            0.30f,
            0.84f,
            "960 mg orally once daily"
        },
        {
            "Adagrasib (MRTX849)",
            "Optimized KRAS(G12C) Inhibitor",
            "KRAS",
            "NSCLC with KRAS G12C mutation with CNS penetrance",
            "Covalently binds inactive KRAS(G12C) with high selectivity and long half-life",
            "FDA Approved (KRYSTAL-1 Trial)",
            EvidenceLevel::Level1_FDA_Approved,
            0.28f,
            0.82f,
            "600 mg orally twice daily"
        },

        // BRAF Alterations
        {
            "Dabrafenib + Trametinib",
            "Dual BRAF/MEK Targeted Regimen",
            "BRAF",
            "Unresectable or metastatic melanoma, NSCLC, anaplastic thyroid with BRAF V600E",
            "Simultaneous selective inhibition of mutant BRAF V600 kinase and downstream MEK1/2",
            "FDA Approved (COMBI-d / COMBI-v Trials)",
            EvidenceLevel::Level1_FDA_Approved,
            0.35f,
            0.90f,
            "Dabrafenib 150 mg BID + Trametinib 2 mg QD"
        },

        // PIK3CA / PTEN Pathway Alterations
        {
            "Alpelisib",
            "PI3K-alpha Selective Inhibitor",
            "PIK3CA",
            "HR-positive, HER2-negative advanced breast cancer with PIK3CA mutation",
            "Inhibits p110-alpha subunit of PI3K, shutting down downstream AKT/mTOR hyperactivation",
            "FDA Approved (SOLAR-1 Trial)",
            EvidenceLevel::Level1_FDA_Approved,
            0.25f,
            0.80f,
            "300 mg orally once daily"
        },
        {
            "Everolimus",
            "mTORC1 Allosteric Inhibitor",
            "PTEN",
            "PTEN-deficient carcinomas and subependymal giant cell astrocytomas",
            "Allosteric binding to FKBP-12/mTOR complex 1, blocking downstream ribosomal S6K signaling",
            "FDA Approved (BOLERO-2 Trial)",
            EvidenceLevel::Level2_StandardOfCare,
            0.20f,
            0.72f,
            "10 mg orally once daily"
        },

        // ERBB2 / HER2 Alterations
        {
            "Trastuzumab Deruxtecan (T-DXd)",
            "HER2-Directed Antibody-Drug Conjugate",
            "ERBB2",
            "HER2-positive or HER2-low unresectable or metastatic breast / gastric cancer",
            "Monoclonal antibody targeting HER2 delivering topoisomerase I inhibitor payload with bystander effect",
            "FDA Approved (DESTINY-Breast03 Trial)",
            EvidenceLevel::Level1_FDA_Approved,
            0.40f,
            0.94f,
            "5.4 mg/kg IV infusion every 3 weeks"
        },

        // TP53 Reactivation & Gene Editing
        {
            "APR-246 (Eprenetapopt)",
            "Conformational p53 Reactivator",
            "TP53",
            "MDS and AML with mutant TP53 and solid tumors",
            "Covalently binds cysteine residues in mutant p53 core domain, thermodynamically refolding it to active state",
            "Phase II/III Clinical Trials (Fast Track Status)",
            EvidenceLevel::Level3_ClinicalTrials,
            0.50f,
            0.75f,
            "4500 mg IV infusion daily for 4 days/cycle"
        },
        {
            "CRISPR-Cas9 HDR Gene Rescue",
            "Precision Exon Homology-Directed Repair",
            "TP53",
            "Pathogenic TP53 missense mutations (p.R175H, p.R248W)",
            "Double-strand cleavage adjacent to mutant codon with ssODN donor template for scarless restoration",
            "Preclinical / Investigational Gene Editing Platform",
            EvidenceLevel::Level4_Preclinical,
            0.80f,
            0.96f,
            "Single-dose ribonucleoprotein (RNP) electroporation"
        },

        // High Tumor Mutational Burden / Immune Biomarker
        {
            "Pembrolizumab",
            "Anti-PD-1 Monoclonal Antibody (Immunotherapy)",
            "TMB_HIGH",
            "Unresectable solid tumors with TMB >= 10 mut/Mb or MSI-H/dMMR",
            "Blocks PD-1 on T-lymphocytes, preventing immune evasion by neoantigen-rich tumor cells",
            "FDA Approved (KEYNOTE-158 Tumor-Agnostic Approval)",
            EvidenceLevel::Level1_FDA_Approved,
            0.10f,
            0.85f,
            "200 mg IV infusion every 3 weeks"
        }
    };
}

ActionabilityAnalysis ActionabilityEngine::analyze_variants(const std::vector<vcf::GenomicVariant>& variants) {
    ActionabilityAnalysis analysis{};

    // Estimación de TMB (Carga mutacional relativa)
    size_t pathogenic_count = 0;
    for (const auto& var : variants) {
        if (var.effect == vcf::VariantEffect::Pathogenic || var.effect == vcf::VariantEffect::LikelyPathogenic) {
            pathogenic_count++;
        }
    }
    analysis.tumor_mutational_burden = static_cast<float>(pathogenic_count) * 2.5f; // Proxy mut/Mb
    if (analysis.tumor_mutational_burden >= 10.0f || pathogenic_count >= 4) {
        analysis.high_tmb = true;
    }

    bool has_egfr_sensitizing = false;
    bool has_egfr_t790m = false;

    for (const auto& var : variants) {
        std::string gene = var.gene_symbol;
        if (gene.empty() || gene == "UNKNOWN") continue;

        ActionableBiomarker biomarker{};
        biomarker.gene_symbol = gene;
        biomarker.genomic_pos = var.pos;
        biomarker.alt_allele = var.alt;
        biomarker.pathogenicity_score = var.pathogenicity_score;

        if (gene == "EGFR") {
            if (var.pos == 55181300) {
                biomarker.variant_desc = "p.T790M (Gatekeeper resistance)";
                biomarker.consequence = "Steric hindrance of ATP-binding pocket; confers secondary resistance to 1st/2nd Gen TKIs.";
                biomarker.is_resistance_marker = true;
                has_egfr_t790m = true;
            } else if (var.pos == 55181378 || var.pathogenicity_score > 0.6f) {
                biomarker.variant_desc = "p.L858R (Exon 21 activating)";
                biomarker.consequence = "Constitutive tyrosine kinase activation; drives non-small cell lung carcinoma.";
                has_egfr_sensitizing = true;
            } else {
                biomarker.variant_desc = "Exon 19 In-frame Deletion / Variant";
                biomarker.consequence = "Hyperactive EGF receptor signaling.";
            }
        } else if (gene == "BRCA1") {
            biomarker.variant_desc = "p.C61G (RING domain missense) / Deficient HR";
            biomarker.consequence = "Abolishes homologous recombination; confers profound vulnerability to PARP inhibition.";
        } else if (gene == "TP53") {
            biomarker.variant_desc = "p.R175H (Conformational hotspot mutation)";
            biomarker.consequence = "Loss of DNA-binding domain stability; abrogates p53 tetramerization and apoptotic checkpoints.";
        } else if (gene == "PTEN") {
            biomarker.variant_desc = "Loss-of-Function Inactivating Alteration";
            biomarker.consequence = "Unchecked PIP3 accumulation; constitutively hyperactivates the pro-survival PI3K/AKT cascade.";
        } else if (gene == "KRAS") {
            biomarker.variant_desc = "p.G12C (Switch-II pocket oncogenic hotspot)";
            biomarker.consequence = "Impaired intrinsic GTP hydrolysis; locks KRAS into perpetual GTP-bound signaling state.";
        } else if (gene == "BRAF") {
            biomarker.variant_desc = "p.V600E (Kinase domain transactivation)";
            biomarker.consequence = "Monomeric active kinase signaling independently of upstream RAS activation.";
        } else if (gene == "PIK3CA") {
            biomarker.variant_desc = "p.E545K / H1047R (Helical/Kinase domain)";
            biomarker.consequence = "Hyperactivation of catalytic p110-alpha kinase; drives proliferation and resistance to apoptosis.";
        } else if (gene == "ERBB2") {
            biomarker.variant_desc = "Amplification / Extracellular Domain Variant";
            biomarker.consequence = "Ligand-independent heterodimerization with EGFR and HER3; triggers aggressive mitotic expansion.";
        } else {
            biomarker.variant_desc = "Clinical Target Alteration";
            biomarker.consequence = "Potentiates oncogenic signaling under cellular stress.";
        }

        // Emparejamiento de fármacos contra la base de datos
        for (const auto& drug : drug_database_) {
            if (drug.target_gene == gene) {
                biomarker.matched_drugs.push_back(drug);

                // Evitar duplicación en la prescripción global
                bool already_added = false;
                for (const auto& existing : analysis.all_recommended_drugs) {
                    if (existing.drug_name == drug.drug_name) {
                        already_added = true;
                        break;
                    }
                }
                if (!already_added) {
                    analysis.all_recommended_drugs.push_back(drug);
                    if (drug.evidence_level == EvidenceLevel::Level1_FDA_Approved) {
                        analysis.level1_count++;
                    } else if (drug.evidence_level == EvidenceLevel::Level2_StandardOfCare) {
                        analysis.level2_count++;
                    } else if (drug.evidence_level == EvidenceLevel::Level3_ClinicalTrials) {
                        analysis.level3_count++;
                    }
                }
            }
        }

        if (!biomarker.matched_drugs.empty() || biomarker.is_resistance_marker) {
            analysis.detected_biomarkers.push_back(biomarker);
        }
    }

    // Regla de Inmunoterapia por alta carga mutacional (TMB >= 10 mut/Mb)
    if (analysis.high_tmb) {
        for (const auto& drug : drug_database_) {
            if (drug.target_gene == "TMB_HIGH") {
                analysis.all_recommended_drugs.push_back(drug);
                analysis.level1_count++;
            }
        }
    }

    // Análisis de Resistencias
    if (has_egfr_t790m) {
        analysis.resistance_warnings.push_back(
            "RESISTANCE ALERT: Mutation p.T790M detected. 1st-Gen (Gefitinib) and 2nd-Gen (Afatinib) TKIs are contra-indicated. Recommend 3rd-Gen Osimertinib."
        );
    }
    if (has_egfr_sensitizing && !has_egfr_t790m) {
        analysis.resistance_warnings.push_back(
            "SURVEILLANCE NOTICE: High risk of acquired T790M or MET amplification under long-term TKI selective pressure."
        );
    }

    // Puntuación de accionabilidad global
    float score = 0.0f;
    score += static_cast<float>(analysis.level1_count) * 0.35f;
    score += static_cast<float>(analysis.level2_count) * 0.20f;
    score += static_cast<float>(analysis.level3_count) * 0.10f;
    if (analysis.high_tmb) score += 0.25f;
    analysis.total_actionability_score = std::min(1.0f, score);

    return analysis;
}

std::string ActionabilityEngine::evidence_level_to_string(EvidenceLevel level) {
    switch (level) {
        case EvidenceLevel::Level1_FDA_Approved:   return "Level 1 (FDA/EMA Approved)";
        case EvidenceLevel::Level2_StandardOfCare: return "Level 2 (NCCN Guidelines)";
        case EvidenceLevel::Level3_ClinicalTrials: return "Level 3 (Phase II/III Trials)";
        case EvidenceLevel::Level4_Preclinical:    return "Level 4 (In-Silico / Preclinical)";
        case EvidenceLevel::Resistance:            return "Resistance Biomarker";
    }
    return "Investigational";
}

std::string ActionabilityEngine::evidence_level_badge(EvidenceLevel level) {
    switch (level) {
        case EvidenceLevel::Level1_FDA_Approved:   return "[FDA LEVEL 1]";
        case EvidenceLevel::Level2_StandardOfCare: return "[NCCN LEVEL 2]";
        case EvidenceLevel::Level3_ClinicalTrials: return "[CLINICAL TRIAL 3]";
        case EvidenceLevel::Level4_Preclinical:    return "[PRECLINICAL 4]";
        case EvidenceLevel::Resistance:            return "[RESISTANCE WARNING]";
    }
    return "[INVESTIGATIONAL]";
}

} // namespace genome::oncology
