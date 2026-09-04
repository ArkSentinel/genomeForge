#include "ipc/state_inspector.hpp"

#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <chrono>
#include <filesystem>

namespace genome::ipc {

bool StateInspector::dump_live_state_json(
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
    size_t sim_speed_multiplier,
    float petri_zoom,
    float petri_pan_x,
    float petri_pan_y
) {
    // Asegurar que el directorio contenedor exista
    std::filesystem::path out_p(filepath);
    if (out_p.has_parent_path()) {
        std::error_code ec;
        std::filesystem::create_directories(out_p.parent_path(), ec);
    }

    std::string temp_path = filepath + ".tmp";
    std::ofstream out(temp_path);
    if (!out.is_open()) return false;

    auto now = std::chrono::system_clock::now();
    auto time_t_now = std::chrono::system_clock::to_time_t(now);
    std::tm tm_buf{};
    localtime_r(&time_t_now, &tm_buf);
    char time_str[32];
    std::strftime(time_str, sizeof(time_str), "%Y-%m-%dT%H:%M:%S", &tm_buf);

    const auto& y_state = yamanaka_sim.get_state();
    const auto& y_factors = yamanaka_sim.get_factors();
    const auto& moe_m = moe_engine.get_metrics();
    const auto& twin_m = cell_twin.get_metrics();

    std::string risk_str = (y_state.r_score > 0.5f) ? "HIGH RISK (Malignant Drift)" :
                           ((y_state.r_score > 0.25f) ? "MODERATE RISK" : "STABLE / LOW RISK");

    out << "{\n"
        << "  \"schema_version\": \"2.0.0\",\n"
        << "  \"timestamp\": \"" << time_str << "\",\n"
        << "  \"system_state\": {\n"
        << "    \"active_tab_index\": " << active_tab << ",\n"
        << "    \"active_locus\": \"" << active_locus << "\",\n"
        << "    \"chip_model\": \"" << hw_metrics.chip_model << "\",\n"
        << "    \"cpu_total_usage\": " << std::fixed << std::setprecision(1) << hw_metrics.total_cpu_usage << ",\n"
        << "    \"p_core_avg_usage\": " << hw_metrics.p_core_avg_usage << ",\n"
        << "    \"e_core_avg_usage\": " << hw_metrics.e_core_avg_usage << ",\n"
        << "    \"memory_used_gb\": " << std::setprecision(2) << (hw_metrics.used_memory_bytes / 1e9) << ",\n"
        << "    \"memory_total_gb\": " << std::setprecision(2) << (hw_metrics.total_memory_bytes / 1e9) << ",\n"
        << "    \"memory_used_percent\": " << std::setprecision(1) << hw_metrics.memory_used_percent << ",\n"
        << "    \"realtime_speed_multiplier\": " << sim_speed_multiplier << "\n"
        << "  },\n"
        << "  \"genomics\": {\n"
        << "    \"total_variants\": " << vcf_parser.total_variants() << ",\n"
        << "    \"pathogenic_variants\": " << vcf_parser.pathogenic_count() << ",\n"
        << "    \"yamanaka_collisions\": " << vcf_parser.yamanaka_collision_count() << ",\n"
        << "    \"average_af\": " << std::setprecision(4) << vcf_parser.average_allele_frequency() << ",\n"
        << "    \"max_pathogenicity_score\": " << std::setprecision(4) << vcf_parser.max_pathogenicity_score() << "\n"
        << "  },\n"
        << "  \"yamanaka_dynamics\": {\n"
        << "    \"pluripotency\": " << std::setprecision(4) << y_state.pluripotency << ",\n"
        << "    \"senescence\": " << std::setprecision(4) << y_state.senescence << ",\n"
        << "    \"genomic_instability\": " << std::setprecision(4) << y_state.instability << ",\n"
        << "    \"p53_level\": " << std::setprecision(4) << y_state.p53_activity << ",\n"
        << "    \"p21_level\": " << std::setprecision(4) << y_state.p21_activity << ",\n"
        << "    \"r_score\": " << std::setprecision(4) << y_state.r_score << ",\n"
        << "    \"risk_category\": \"" << risk_str << "\",\n"
        << "    \"oct4\": " << std::setprecision(3) << y_factors.oct4 << ",\n"
        << "    \"sox2\": " << std::setprecision(3) << y_factors.sox2 << ",\n"
        << "    \"klf4\": " << std::setprecision(3) << y_factors.klf4 << ",\n"
        << "    \"cmyc\": " << std::setprecision(3) << y_factors.cmyc << "\n"
        << "  },\n"
        << "  \"moe_routing\": {\n"
        << "    \"cpu_accelerate_processed\": " << moe_m.cpu_vdsp_dispatches << ",\n"
        << "    \"metal_gpu_batches_processed\": " << moe_m.gpu_metal_dispatches << ",\n"
        << "    \"last_gpu_latency_us\": " << std::setprecision(2) << moe_m.avg_metal_latency_us << "\n"
        << "  },\n"
        << "  \"crispr_studio\": {\n"
        << "    \"candidates_count\": " << crispr_candidates.size() << ",\n"
        << "    \"crispr_applied\": " << (crispr_applied ? "true" : "false") << ",\n"
        << "    \"best_guide_spacer\": \"" << (crispr_candidates.empty() ? "NONE" : crispr_candidates[0].protospacer_20nt) << "\",\n"
        << "    \"best_guide_pam\": \"" << (crispr_candidates.empty() ? "NONE" : crispr_candidates[0].pam) << "\",\n"
        << "    \"best_guide_efficiency\": " << std::setprecision(2) << (crispr_candidates.empty() ? 0.0f : crispr_candidates[0].on_target_efficiency) << "\n"
        << "  },\n"
        << "  \"digital_twin\": {\n"
        << "    \"total_cells\": " << twin_m.total_cells << ",\n"
        << "    \"somatic_count\": " << twin_m.somatic_count << ",\n"
        << "    \"intermediate_count\": " << twin_m.intermediate_count << ",\n"
        << "    \"ipsc_count\": " << twin_m.ipsc_count << ",\n"
        << "    \"transformed_count\": " << twin_m.transformed_count << ",\n"
        << "    \"apoptotic_count\": " << twin_m.apoptotic_count << ",\n"
        << "    \"confluence_percent\": " << std::setprecision(2) << (twin_m.confluence * 100.0f) << ",\n"
        << "    \"mitotic_index\": " << std::setprecision(2) << twin_m.mitotic_index << ",\n"
        << "    \"ipsc_efficiency\": " << std::setprecision(2) << twin_m.ipsc_efficiency << ",\n"
        << "    \"average_r_score\": " << std::setprecision(4) << twin_m.average_r_score << ",\n"
        << "    \"doxycycline_active\": " << (twin_m.dox_active ? "true" : "false") << ",\n"
        << "    \"crispr_rescued\": " << (twin_m.crispr_applied ? "true" : "false") << ",\n"
        << "    \"time_warp_scale\": \"" << cell_twin.get_time_scale_label() << "\",\n"
        << "    \"time_multiplier\": " << std::setprecision(1) << twin_m.time_multiplier << ",\n"
        << "    \"simulated_seconds\": " << std::setprecision(2) << twin_m.simulated_seconds << ",\n"
        << "    \"simulated_time\": \"" << twin_m.simulated_time_str << "\",\n"
        << "    \"petri_zoom_level\": " << std::fixed << std::setprecision(2) << petri_zoom << ",\n"
        << "    \"petri_pan_x\": " << std::setprecision(2) << petri_pan_x << ",\n"
        << "    \"petri_pan_y\": " << std::setprecision(2) << petri_pan_y << "\n"
        << "  }\n"
        << "}\n";

    out.close();
    std::rename(temp_path.c_str(), filepath.c_str());
    return true;
}

bool StateInspector::print_inspection_report(const std::string& filepath) {
    std::ifstream in(filepath);
    if (!in.is_open()) {
        std::vector<std::string> fallbacks = {
            "../" + filepath,
            "data/" + filepath,
            "../data/live_state.json"
        };
        for (const auto& alt : fallbacks) {
            in.open(alt);
            if (in.is_open()) break;
        }
    }

    if (!in.is_open()) {
        std::cerr << "[!] Error: No active simulation telemetry found at " << filepath << std::endl;
        std::cerr << "    Run GenomeForge-TUI first to generate live telemetry state." << std::endl;
        return false;
    }

    std::stringstream buffer;
    buffer << in.rdbuf();
    std::string json_content = buffer.str();

    std::cout << "\n";
    std::cout << "╔══════════════════════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║ ⌬ GENOMEFORGE-TUI :: LIVE STATE TELEMETRY REPORT (AI & CLI INSPECTOR)           ║\n";
    std::cout << "╚══════════════════════════════════════════════════════════════════════════════════╝\n";
    std::cout << json_content << "\n";
    return true;
}

} // namespace genome::ipc
