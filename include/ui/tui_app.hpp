#pragma once

#include "system/telemetry.hpp"
#include "vcf/vcf_parser.hpp"
#include "vcf/locus_catalog.hpp"
#include "yamanaka/yamanaka_sim.hpp"
#include "moe/moe_engine.hpp"
#include "crispr/crispr_engine.hpp"
#include "digital_twin/cell_twin.hpp"
#include "oncology/actionability_engine.hpp"
#include "export/report_generator.hpp"

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
    SystemLogs = 5,
    DigitalTwin = 6,
    OncoTherapy = 7
};

class GenomeForgeApp {
public:
    GenomeForgeApp(const std::string& vcf_path);
    ~GenomeForgeApp();

    // Inicia el loop de renderizado interactivo en la terminal
    int run();

    // Ejecuta la simulación en modo headless (sin interfaz) para auditoría o pruebas
    int run_headless(size_t max_ticks, bool export_csv = true, const std::string& snapshot_path = "data/live_state.json");

    // Exporta datasets en formato CSV
    void export_datasets();

    // Exporta reporte de Comité Molecular de Tumores
    void export_molecular_tumor_board_report(const std::string& filepath = "data/molecular_tumor_board_report.md");

    // Control de velocidad de simulación en tiempo real (1x a 32x)
    void increase_sim_speed();
    void decrease_sim_speed();
    size_t get_sim_speed_multiplier() const { return sim_speed_multiplier_; }

    // Control de zoom y pan óptico de la placa de Petri
    void zoom_in_petri();
    void zoom_out_petri();
    void reset_petri_view();
    void pan_petri(float dx, float dy);
    float get_petri_zoom() const { return petri_zoom_; }
    float get_petri_pan_x() const { return petri_pan_x_; }
    float get_petri_pan_y() const { return petri_pan_y_; }

    // Controles farmacológicos in-silico
    void administer_targeted_drug(genome::twin::DrugTherapy drug);
    void washout_drugs();

private:
    void init_simulation_state();
    void refresh_onco_analysis();
    void add_log(const std::string& level, const std::string& message);
    void tick_simulation();
    void switch_locus_by_index(size_t index);
    void cycle_locus();
    void toggle_sample_mode();
    void save_live_state();

    bool healthy_reference_mode_{false};
    size_t tick_counter_{0};
    size_t sim_speed_multiplier_{1};
    float petri_zoom_{1.0f};
    float petri_pan_x_{0.0f};
    float petri_pan_y_{0.0f};
    std::string vcf_path_;
    genome::vcf::VcfParser vcf_parser_;
    genome::system::TelemetryCollector telemetry_;

    // Motores bioinformáticos y matemáticos reales
    genome::yamanaka::YamanakaSimulator yamanaka_sim_;
    genome::moe::MoEEngine moe_engine_;
    genome::crispr::CrisprEngine crispr_engine_;
    genome::twin::CellTwinSimulator cell_twin_;
    genome::oncology::ActionabilityEngine onco_engine_;
    genome::oncology::ActionabilityAnalysis onco_analysis_;

    // Pestaña activa
    AppTab current_tab_{AppTab::Dashboard};

    // Estado interactivo de la UI
    int selected_variant_index_{0};
    size_t active_locus_index_{0};
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
