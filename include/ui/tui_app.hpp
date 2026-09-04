#pragma once

#include "system/telemetry.hpp"
#include "vcf/vcf_parser.hpp"
#include "yamanaka/yamanaka_sim.hpp"
#include "moe/moe_engine.hpp"
#include "crispr/crispr_engine.hpp"

#include <string>
#include <vector>
#include <memory>
#include <chrono>

namespace genome::ui {

enum class AppTab : int {
    Dashboard = 0,
    VcfExplorer = 1,
    YamanakaDynamics = 2,
    MoeRouter = 3,
    CrisprStudio = 4,
    SystemLogs = 5
};

class GenomeForgeApp {
public:
    GenomeForgeApp(const std::string& vcf_path);
    ~GenomeForgeApp();

    // Inicia el loop de renderizado interactivo en la terminal
    int run();

private:
    void init_simulation_state();
    void add_log(const std::string& level, const std::string& message);
    void tick_simulation();

    std::string vcf_path_;
    genome::vcf::VcfParser vcf_parser_;
    genome::system::TelemetryCollector telemetry_;

    // Motores bioinformáticos y matemáticos reales
    genome::yamanaka::YamanakaSimulator yamanaka_sim_;
    genome::moe::MoEEngine moe_engine_;
    genome::crispr::CrisprEngine crispr_engine_;

    // Pestaña activa
    AppTab current_tab_{AppTab::Dashboard};

    // Estado interactivo de la UI
    int selected_variant_index_{0};
    bool yamanaka_running_{true};

    // Estado de CRISPR
    std::vector<genome::crispr::GuideCandidate> crispr_candidates_;
    int selected_guide_index_{0};
    bool crispr_applied_{false};
    std::string target_sequence_window_;

    // Logs interactivos
    struct LogEntry {
        std::string timestamp;
        std::string level;
        std::string message;
    };
    std::vector<LogEntry> logs_;
};

} // namespace genome::ui
