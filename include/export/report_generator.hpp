#pragma once

#include "vcf/vcf_parser.hpp"
#include "oncology/actionability_engine.hpp"
#include "crispr/crispr_engine.hpp"
#include "digital_twin/cell_twin.hpp"
#include "yamanaka/yamanaka_sim.hpp"

#include <string>

namespace genome::exp {

class ReportGenerator {
public:
    // Genera el reporte clínico completo de Comité Molecular de Tumores en Markdown
    static bool generate_molecular_tumor_board_report(
        const std::string& filepath,
        const std::string& sample_id,
        const vcf::VcfParser& parser,
        const oncology::ActionabilityAnalysis& onco_analysis,
        const std::vector<crispr::GuideCandidate>& crispr_candidates,
        const twin::CellTwinSimulator& cell_twin,
        const yamanaka::YamanakaSimulator& yamanaka_sim
    );

    // Genera el reporte clínico estructurado en JSON para integración FHIR/AI
    static bool generate_clinical_report_json(
        const std::string& filepath,
        const std::string& sample_id,
        const vcf::VcfParser& parser,
        const oncology::ActionabilityAnalysis& onco_analysis,
        const std::vector<crispr::GuideCandidate>& crispr_candidates,
        const twin::CellTwinSimulator& cell_twin,
        const yamanaka::YamanakaSimulator& yamanaka_sim
    );
};

} // namespace genome::exp
