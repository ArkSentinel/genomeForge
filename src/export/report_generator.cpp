#include "export/report_generator.hpp"

#include <fstream>
#include <sstream>
#include <iomanip>
#include <chrono>
#include <ctime>

namespace genome::exp {

static std::string get_current_iso_timestamp() {
    auto now = std::chrono::system_clock::now();
    std::time_t now_c = std::chrono::system_clock::to_time_t(now);
    std::tm tm{};
    localtime_r(&now_c, &tm);
    char buf[32];
    std::strftime(buf, sizeof(buf), "%Y-%m-%d %H:%M:%S", &tm);
    return std::string(buf);
}

bool ReportGenerator::generate_molecular_tumor_board_report(
    const std::string& filepath,
    const std::string& sample_id,
    const vcf::VcfParser& parser,
    const oncology::ActionabilityAnalysis& onco_analysis,
    const std::vector<crispr::GuideCandidate>& crispr_candidates,
    const twin::CellTwinSimulator& cell_twin,
    const yamanaka::YamanakaSimulator& yamanaka_sim
) {
    std::ofstream out(filepath);
    if (!out.is_open()) return false;

    const auto& y_state = yamanaka_sim.get_state();
    const auto& metrics = cell_twin.get_metrics();
    std::string risk_str = (y_state.r_score > 0.5f ? "HIGH RISK (Malignant Drift)" : 
                           (y_state.r_score > 0.25f ? "MODERATE RISK" : "STABLE / LOW RISK"));

    out << "# MOLECULAR TUMOR BOARD CLINICAL SUMMARY REPORT\n";
    out << "## Precision Oncology & Targeted Therapeutics Evaluation\n\n";
    out << "================================================================================\n";
    out << "Fecha del Informe:          " << get_current_iso_timestamp() << "\n";
    out << "ID de Muestra / Paciente:   " << sample_id << "\n";
    out << "Genoma de Referencia:       GRCh38 / hg38 (1000 Genomes Standard)\n";
    out << "Motor de Análisis:          GenomeForge-TUI v2.2 (Apple Silicon AMX/Metal)\n";
    out << "Archivo VCF Origen:         " << parser.get_filepath() << "\n";
    out << "================================================================================\n\n";

    // 1. DIAGNÓSTICO GENÓMICO GLOBAL
    out << "### 1. DIAGNÓSTICO GENÓMICO Y BIOMARCADORES DE RIESGO\n\n";
    out << "| Parámetro Clínico | Valor Observado | Interpretación Diagnóstica |\n";
    out << "| :--- | :---: | :--- |\n";
    out << "| **Variantes Totales Analizadas** | " << parser.total_variants() << " | Cobertura en loci oncológicos diana |\n";
    out << "| **Variantes Patogénicas (Drivers)** | " << parser.pathogenic_count() << " | Alteraciones conductoras confirmadas |\n";
    out << "| **Colisiones de Reprogramación** | " << parser.yamanaka_collision_count() << " | Conflictos con factores de pluripotencia |\n";
    out << "| **Carga Mutacional Tumoral (TMB)** | " << std::fixed << std::setprecision(1) << onco_analysis.tumor_mutational_burden << " mut/Mb | " 
        << (onco_analysis.high_tmb ? "ALTA (Candidato a Inmunoterapia Anti-PD-1)" : "BAJA / MODERADA") << " |\n";
    out << "| **Puntaje de Riesgo Oncogénico (R-Score)** | " << std::fixed << std::setprecision(4) << y_state.r_score << " | "
        << risk_str << " |\n";
    out << "| **Integridad Funcional de TP53** | " << std::fixed << std::setprecision(2) << (y_state.p53_activity * 100.0f) << "% | "
        << (y_state.p53_activity < 0.3f ? "COLAPSO TUMORAL (Pérdida de apoptosis)" : "VIGILANCIA ACTIVA") << " |\n\n";

    // 2. ALTERACIONES CONDUCTORAS DETECTADAS
    out << "### 2. ALTERACIONES MOLECULARES CONDUCTORAS (ONCOGENIC DRIVERS)\n\n";
    if (onco_analysis.detected_biomarkers.empty()) {
        out << "*No se identificaron mutaciones patogénicas de alta prioridad en los loci analizados. Perfil genómico compatible con control sano.*\n\n";
    } else {
        out << "| Gen | Coordenada GRCh38 | Alteración / Codón | Consecuencia Fisiopatológica | Estado |\n";
        out << "| :--- | :---: | :--- | :--- | :--- |\n";
        for (const auto& bio : onco_analysis.detected_biomarkers) {
            out << "| **" << bio.gene_symbol << "** | " << bio.genomic_pos << " | " << bio.variant_desc 
                << " | " << bio.consequence << " | " 
                << (bio.is_resistance_marker ? "[RESISTENCIA SECUNDARIA]" : "[ONCOGÉNICO ACCIONABLE]") << " |\n";
        }
        out << "\n";
    }

    // 3. RECOMENDACIONES DE TERAPIA DIRIGIDA
    out << "### 3. RECOMENDACIONES DE TERAPIA DIRIGIDA (ACTIONABLE ONCOLOGY PRESCRIPTION)\n\n";
    if (onco_analysis.all_recommended_drugs.empty()) {
        out << "*No se requieren terapias dirigidas citotóxicas para este perfil. Reprogramación fisiológica celular recomendada.*\n\n";
    } else {
        out << "| Fármaco | Nivel Evidencia | Clase Terapéutica | Gen Diana | Mecanismo de Acción / Ensayo |\n";
        out << "| :--- | :---: | :--- | :---: | :--- |\n";
        for (const auto& drug : onco_analysis.all_recommended_drugs) {
            out << "| **" << drug.drug_name << "** | " << oncology::ActionabilityEngine::evidence_level_badge(drug.evidence_level)
                << " | " << drug.drug_class << " | " << drug.target_gene << " | " << drug.mechanism_of_action 
                << " (" << drug.clinical_trial_or_approval << ") |\n";
        }
        out << "\n";
    }

    // 4. ADVERTENCIAS DE RESISTENCIA Y CONTRAINDICACIONES
    out << "### 4. ADVERTENCIAS FARMACOGENÓMICAS Y MONITOREO CLONAL\n\n";
    if (onco_analysis.resistance_warnings.empty()) {
        out << "[+] No se detectaron biomarcadores conocidos de resistencia primaria ni contraindicaciones farmacogenómicas activas.\n\n";
    } else {
        for (const auto& warn : onco_analysis.resistance_warnings) {
            out << "[!] " << warn << "\n";
        }
        out << "\n";
    }

    // 5. INTERVENCIÓN POR EDICIÓN GÉNICA CRISPR-CAS9 (HDR RESCUE)
    out << "### 5. PLAN DE RESCATE POR EDICIÓN GÉNICA CRISPR-CAS9 (HDR IN-SILICO)\n\n";
    if (crispr_candidates.empty()) {
        out << "*No se requirió diseño de guías CRISPR para esta muestra.*\n\n";
    } else {
        const auto& best = crispr_candidates[0];
        out << "* **Locus Diana:** " << parser.get_active_locus() << "\n";
        out << "* **Secuencia de Guía (Protospacer 20 nt):** `5'-" << best.protospacer_20nt << "-3'`\n";
        out << "* **Motivo PAM (NGG):** `" << best.pam << "` (Offset relativo: " << best.pam_position_offset << ")\n";
        out << "* **Score de Eficiencia On-Target (Rule Set 2):** " << std::fixed << std::setprecision(1) << best.on_target_efficiency << " / 100.0\n";
        out << "* **Molde de Reparación Donante (ssODN HDR):** `5'-TGCAGCTGTGGGTTGATTCCACACCCCCGCCCGGCACCCGCGTCCGCGCCATGGC-3'`\n";
        out << "* **Objetivo Terapéutico:** Sustituir codón mutado restaurando el marco de lectura salvaje y la funcionalidad de apoptosis.\n\n";
    }

    // 6. PROYECCIÓN FARMACODINÁMICA DEL GEMELO DIGITAL CELULAR
    out << "### 6. PRONÓSTICO FARMACODINÁMICO DEL GEMELO DIGITAL (IN-SILICO SCREEN)\n\n";
    out << "* **Fármaco Activo en Ensayo:** " << metrics.active_drug_name << "\n";
    out << "* **Dosis Simulada:** " << std::fixed << std::setprecision(2) << metrics.drug_dose_uM << " uM\n";
    out << "* **Tasa de Regresión Tumoral Calculada:** " << std::fixed << std::setprecision(1) << metrics.tumor_regression_percent << "%\n";
    out << "* **Población Celular:** " << metrics.total_cells << " células totales (" << metrics.transformed_count << " clones tumorales, " << metrics.somatic_count << " somáticas)\n";
    out << "* **Clones con Resistencia Emergente:** " << metrics.drug_resistant_cells << " células (" 
        << (metrics.drug_resistant_cells > 0 ? "Requiere terapia combinada dual" : "Sensibilidad completa") << ")\n";
    out << "* **Confluencia Final del Tejido:** " << std::fixed << std::setprecision(1) << (metrics.confluence * 100.0f) << "%\n\n";

    out << "---\n";
    out << "*Informe generado automáticamente por GenomeForge-TUI. Destinado a evaluación preclínica e investigación traslacional.*\n";

    return true;
}

bool ReportGenerator::generate_clinical_report_json(
    const std::string& filepath,
    const std::string& sample_id,
    const vcf::VcfParser& parser,
    const oncology::ActionabilityAnalysis& onco_analysis,
    const std::vector<crispr::GuideCandidate>& crispr_candidates,
    const twin::CellTwinSimulator& cell_twin,
    const yamanaka::YamanakaSimulator& yamanaka_sim
) {
    std::ofstream out(filepath);
    if (!out.is_open()) return false;

    const auto& y_state = yamanaka_sim.get_state();
    const auto& metrics = cell_twin.get_metrics();

    out << "{\n";
    out << "  \"report_type\": \"MOLECULAR_TUMOR_BOARD_SUMMARY\",\n";
    out << "  \"timestamp\": \"" << get_current_iso_timestamp() << "\",\n";
    out << "  \"sample_id\": \"" << sample_id << "\",\n";
    out << "  \"reference_genome\": \"GRCh38\",\n";
    out << "  \"vcf_source\": \"" << parser.get_filepath() << "\",\n";

    out << "  \"genomics\": {\n";
    out << "    \"total_variants\": " << parser.total_variants() << ",\n";
    out << "    \"pathogenic_drivers\": " << parser.pathogenic_count() << ",\n";
    out << "    \"yamanaka_collisions\": " << parser.yamanaka_collision_count() << ",\n";
    out << "    \"tmb_mut_per_mb\": " << std::fixed << std::setprecision(2) << onco_analysis.tumor_mutational_burden << ",\n";
    out << "    \"high_tmb\": " << (onco_analysis.high_tmb ? "true" : "false") << ",\n";
    out << "    \"r_score\": " << std::fixed << std::setprecision(4) << y_state.r_score << ",\n";
    out << "    \"p53_integrity\": " << std::fixed << std::setprecision(4) << y_state.p53_activity << "\n";
    out << "  },\n";

    out << "  \"actionable_biomarkers\": [\n";
    for (size_t i = 0; i < onco_analysis.detected_biomarkers.size(); ++i) {
        const auto& b = onco_analysis.detected_biomarkers[i];
        out << "    {\n";
        out << "      \"gene\": \"" << b.gene_symbol << "\",\n";
        out << "      \"variant\": \"" << b.variant_desc << "\",\n";
        out << "      \"position\": " << b.genomic_pos << ",\n";
        out << "      \"is_resistance\": " << (b.is_resistance_marker ? "true" : "false") << "\n";
        out << "    }" << (i + 1 < onco_analysis.detected_biomarkers.size() ? "," : "") << "\n";
    }
    out << "  ],\n";

    out << "  \"recommended_therapies\": [\n";
    for (size_t i = 0; i < onco_analysis.all_recommended_drugs.size(); ++i) {
        const auto& d = onco_analysis.all_recommended_drugs[i];
        out << "    {\n";
        out << "      \"drug_name\": \"" << d.drug_name << "\",\n";
        out << "      \"drug_class\": \"" << d.drug_class << "\",\n";
        out << "      \"target\": \"" << d.target_gene << "\",\n";
        out << "      \"evidence_level\": \"" << oncology::ActionabilityEngine::evidence_level_to_string(d.evidence_level) << "\",\n";
        out << "      \"efficacy\": " << std::fixed << std::setprecision(2) << d.tumor_regression_efficacy << "\n";
        out << "    }" << (i + 1 < onco_analysis.all_recommended_drugs.size() ? "," : "") << "\n";
    }
    out << "  ],\n";

    out << "  \"digital_twin_pharmacodynamics\": {\n";
    out << "    \"active_drug\": \"" << metrics.active_drug_name << "\",\n";
    out << "    \"dose_uM\": " << std::fixed << std::setprecision(2) << metrics.drug_dose_uM << ",\n";
    out << "    \"tumor_regression_percent\": " << std::fixed << std::setprecision(1) << metrics.tumor_regression_percent << ",\n";
    out << "    \"resistant_cells_count\": " << metrics.drug_resistant_cells << "\n";
    out << "  }\n";
    out << "}\n";

    return true;
}

} // namespace genome::exp
