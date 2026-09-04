#pragma once

#include "system/telemetry.hpp"
#include "vcf/vcf_parser.hpp"
#include "yamanaka/yamanaka_sim.hpp"
#include "moe/moe_engine.hpp"
#include "crispr/crispr_engine.hpp"
#include "digital_twin/cell_twin.hpp"

#include <string>

namespace genome::ipc {

class StateInspector {
public:
    // Serializa el estado completo en vivo a JSON y lo escribe de forma atómica en disco
    static bool dump_live_state_json(
        const std::string& filepath,
        int active_tab,
        const std::string& active_locus,
        const genome::system::HardwareTelemetry& hw_metrics,
        const genome::vcf::VcfParser& vcf_parser,
        const genome::yamanaka::YamanakaSimulator& yamanaka_sim,
        const genome::moe::MoEEngine& moe_engine,
        const std::vector<genome::crispr::GuideCandidate>& crispr_candidates,
        bool crispr_applied,
        const genome::twin::CellTwinSimulator& cell_twin,
        size_t sim_speed_multiplier = 1,
        float petri_zoom = 1.0f,
        float petri_pan_x = 0.0f,
        float petri_pan_y = 0.0f
    );

    // Lee el archivo de estado en vivo e imprime un reporte estructurado de diagnóstico en stdout
    static bool print_inspection_report(const std::string& filepath);
};

} // namespace genome::ipc
