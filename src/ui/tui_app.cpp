#include "ui/tui_app.hpp"

#include <iostream>
#include <ftxui/component/app.hpp>
#include <ftxui/component/component.hpp>
#include <ftxui/component/event.hpp>
#include <ftxui/component/loop.hpp>
#include <ftxui/dom/elements.hpp>
#include <ftxui/dom/canvas.hpp>
#include <ftxui/screen/color.hpp>

#include "export/csv_exporter.hpp"
#include "ipc/state_inspector.hpp"

#include <chrono>
#include <iomanip>
#include <sstream>
#include <thread>
#include <cmath>
#include <algorithm>

namespace genome::ui {

using namespace ftxui;

GenomeForgeApp::GenomeForgeApp(const std::string& vcf_path)
    : vcf_path_(vcf_path) {
    init_simulation_state();
}

GenomeForgeApp::~GenomeForgeApp() = default;

void GenomeForgeApp::add_log(const std::string& level, const std::string& message) {
    auto now = std::chrono::system_clock::now();
    auto time_t_now = std::chrono::system_clock::to_time_t(now);
    std::tm tm_buf{};
    localtime_r(&time_t_now, &tm_buf);

    std::ostringstream oss;
    oss << std::setfill('0') << std::setw(2) << tm_buf.tm_hour << ":"
        << std::setfill('0') << std::setw(2) << tm_buf.tm_min << ":"
        << std::setfill('0') << std::setw(2) << tm_buf.tm_sec;

    logs_.push_back({oss.str(), level, message});
    if (logs_.size() > 50) {
        logs_.erase(logs_.begin());
    }
}

void GenomeForgeApp::switch_locus_by_index(size_t index) {
    active_locus_index_ = index % genome::vcf::LocusCatalog::total_loci();
    const auto* locus = genome::vcf::LocusCatalog::get_by_index(active_locus_index_);
    if (!locus) return;

    vcf_parser_.switch_to_locus(locus->gene_symbol);
    selected_variant_index_ = 0;

    // Actualizar contexto biológico para el simulador de Yamanaka
    bool tp53_mut = false;
    bool cdkn2a_mut = false;
    for (const auto& v : vcf_parser_.get_variants()) {
        std::string_view gene(v.gene_symbol);
        if (gene == "TP53" && (v.effect == genome::vcf::VariantEffect::YamanakaCollision || 
                               v.effect == genome::vcf::VariantEffect::Pathogenic)) {
            tp53_mut = true;
        }
        if (gene == "CDKN2A" && (v.effect == genome::vcf::VariantEffect::YamanakaCollision || 
                                 v.effect == genome::vcf::VariantEffect::Pathogenic)) {
            cdkn2a_mut = true;
        }
    }
    float burden = static_cast<float>(vcf_parser_.pathogenic_count()) / 
                   std::max(1.0f, static_cast<float>(vcf_parser_.total_variants()));
    yamanaka_sim_.set_mutational_context(tp53_mut, cdkn2a_mut, burden);

    if (vcf_parser_.total_variants() > 0) {
        moe_engine_.route_batch(vcf_parser_.get_variants());
    }

    // Re-sintetizar guías CRISPR acordes al locus activo
    if (locus->gene_symbol == "TP53") {
        target_sequence_window_ = "CTGCTCCCCCCCTGCCCTCAACAAGATGTTTTGCCAGTTGGCTCTGACTGTACCACCATCCACTACAACT";
    } else if (locus->gene_symbol == "CDKN2A") {
        target_sequence_window_ = "ATGGAGCCGGCGGCGGGGAGCAGCATGGAGCCGGCGGCGGGGAGCAGCATGGAGCCTTCGGCTGACTGGCT";
    } else if (locus->gene_symbol == "PTEN") {
        target_sequence_window_ = "ATGACAGCCATCATCAAAGAGATCGTTAGCAGAAACAAAAGGAGATATCAAGAGGATGGATTCGACTTAGA";
    } else if (locus->gene_symbol == "MYC") {
        target_sequence_window_ = "ATGCCCCTCAACGTTAGCTTCACCAACAGGAACTATGACCTCGACTACGACTCGGTGCAGCCGTATTTCTA";
    } else if (locus->gene_symbol == "BRCA1") {
        target_sequence_window_ = "ATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATTAATGCTATGCAGAAAATCTTAGAGTG";
    } else if (locus->gene_symbol == "EGFR") {
        target_sequence_window_ = "ATGCGACCCTCCGGGACGGCCGGGGCAGCGCTCCTGGCGCTGCTGGCTGCGCTCTGCCCGGCGAGTCGGGC";
    } else if (locus->gene_symbol == "POU5F1") {
        target_sequence_window_ = "ATGGCGGGACACCTGGCTTCGGATTTCGCCTTCTCGCCCCCTCCAGGTGGTGGAGGTGATGGGCCAGGGGG";
    } else if (locus->gene_symbol == "SOX2") {
        target_sequence_window_ = "ATGTACAACATGATGGAGACGGAGCTGAAGCCGCCGGGCCCGCAGCAAACTTCGGGGGGCGGCGGCGGCAA";
    } else if (locus->gene_symbol == "KRAS") {
        target_sequence_window_ = "ATGACTGAATATAAACTTGTGGTAGTTGGAGCTGGTGGCGTAGGCAAGAGTGCCTTGACGATACAGCTAAT";
    } else if (locus->gene_symbol == "BRAF") {
        target_sequence_window_ = "ATGGCGGCGCTGAGCGGTGGCGGTGGTGGCGGCGCGGAGCCCGGCCAGGCTCTGTTCAACGGGGACATGGAG";
    } else if (locus->gene_symbol == "PIK3CA") {
        target_sequence_window_ = "ATGCCTCCACGACCATCATCAGGTGAACTGTGGGGCATCCACTTGATGCCCCCAAGAATCCTAGTAGAATGT";
    } else if (locus->gene_symbol == "ERBB2") {
        target_sequence_window_ = "ATGGAGCTGGCGGCCTTGTGCCGCTGGGGGCTCCTCCTCGCCCTCTTGCCCCCCGGAGCCGCGAGCACCCAA";
    }
    crispr_candidates_ = crispr_engine_.design_guides(target_sequence_window_, 20);
    crispr_applied_ = false;

    refresh_onco_analysis();

    add_log("LOCUS", "1000 Genomes focus: " + locus->gene_symbol + " (" + locus->chromosome + ":" + 
                     std::to_string(locus->start_pos) + "-" + std::to_string(locus->end_pos) + ")");
}

void GenomeForgeApp::cycle_locus() {
    switch_locus_by_index(active_locus_index_ + 1);
}

void GenomeForgeApp::toggle_sample_mode() {
    healthy_reference_mode_ = !healthy_reference_mode_;
    vcf_path_ = healthy_reference_mode_ ? "data/1000genomes_healthy_reference.vcf" : "data/sample_patient.vcf";
    
    if (vcf_parser_.load_file(vcf_path_)) {
        selected_variant_index_ = 0;
        bool tp53_mut = false;
        bool cdkn2a_mut = false;
        for (const auto& v : vcf_parser_.get_variants()) {
            std::string_view gene(v.gene_symbol);
            if (gene == "TP53" && (v.effect == genome::vcf::VariantEffect::YamanakaCollision || 
                                   v.effect == genome::vcf::VariantEffect::Pathogenic)) {
                tp53_mut = true;
            }
            if (gene == "CDKN2A" && (v.effect == genome::vcf::VariantEffect::YamanakaCollision || 
                                     v.effect == genome::vcf::VariantEffect::Pathogenic)) {
                cdkn2a_mut = true;
            }
        }
        float burden = static_cast<float>(vcf_parser_.pathogenic_count()) / 
                       std::max(1.0f, static_cast<float>(vcf_parser_.total_variants()));
        yamanaka_sim_.set_mutational_context(tp53_mut, cdkn2a_mut, burden);
        refresh_onco_analysis();

        if (healthy_reference_mode_) {
            add_log("SAMPLE", "Switched to 1000 Genomes Healthy Reference (HG00096): 16 natural polymorphisms, 0 pathogenic collisions.");
        } else {
            add_log("SAMPLE", "Switched to Clinical Patient Profile: TP53 R175H & CDKN2A pathogenic mutations active for CRISPR rescue.");
        }
    } else {
        add_log("ERROR", "Failed to switch sample to: " + vcf_path_);
    }
}

void GenomeForgeApp::refresh_onco_analysis() {
    onco_analysis_ = onco_engine_.analyze_variants(vcf_parser_.get_variants());
}

void GenomeForgeApp::export_molecular_tumor_board_report(const std::string& filepath) {
    refresh_onco_analysis();
    std::string sample_id = healthy_reference_mode_ ? "HG00096_HEALTHY_CONTROL" : "PATIENT_ONCOLOGY_SAMPLE";
    bool ok = genome::exp::ReportGenerator::generate_molecular_tumor_board_report(
        filepath, sample_id, vcf_parser_, onco_analysis_, crispr_candidates_, cell_twin_, yamanaka_sim_
    );
    std::string json_path = filepath;
    if (json_path.size() >= 3 && json_path.substr(json_path.size() - 3) == ".md") {
        json_path = json_path.substr(0, json_path.size() - 3) + ".json";
    } else {
        json_path += ".json";
    }
    genome::exp::ReportGenerator::generate_clinical_report_json(
        json_path, sample_id, vcf_parser_, onco_analysis_, crispr_candidates_, cell_twin_, yamanaka_sim_
    );
    if (ok) {
        add_log("REPORT", "Molecular Tumor Board Report exported to: " + filepath);
        add_log("REPORT", "JSON clinical export saved to: " + json_path);
    } else {
        add_log("ERROR", "Failed to export MTB report to: " + filepath);
    }
}

void GenomeForgeApp::administer_targeted_drug(genome::twin::DrugTherapy drug) {
    cell_twin_.administer_drug(drug, 1.50f);
    add_log("DRUG", "Inoculated in-silico therapy: " + cell_twin_.get_active_drug_name() + " (1.50 uM)");
}

void GenomeForgeApp::washout_drugs() {
    cell_twin_.washout_drugs();
    add_log("WASHOUT", "Culture washed. Active drug cleared from Petri dish.");
}

void GenomeForgeApp::init_simulation_state() {
    add_log("INFO", "GenomeForge-TUI core initialized on Apple Silicon architecture.");
    
    // 1. Cargar archivo VCF
    if (vcf_parser_.load_file(vcf_path_)) {
        add_log("INFO", "VCF loaded: " + std::to_string(vcf_parser_.total_variants()) + 
                        " variants indexed in 64-byte ARM64 cache-aligned buffers.");
        if (vcf_parser_.yamanaka_collision_count() > 0) {
            add_log("WARN", "Detected " + std::to_string(vcf_parser_.yamanaka_collision_count()) + 
                            " critical collision(s) with Yamanaka OSKM factors!");
        }
    } else {
        add_log("ERROR", "Failed to load VCF from: " + vcf_path_);
    }

    // 2. Evaluar contexto mutacional para el simulador de Yamanaka
    bool tp53_mut = false;
    bool cdkn2a_mut = false;
    for (const auto& v : vcf_parser_.get_variants()) {
        std::string_view gene(v.gene_symbol);
        if (gene == "TP53" && (v.effect == genome::vcf::VariantEffect::YamanakaCollision || 
                               v.effect == genome::vcf::VariantEffect::Pathogenic)) {
            tp53_mut = true;
        }
        if (gene == "CDKN2A" && (v.effect == genome::vcf::VariantEffect::YamanakaCollision || 
                                 v.effect == genome::vcf::VariantEffect::Pathogenic)) {
            cdkn2a_mut = true;
        }
    }
    float burden = static_cast<float>(vcf_parser_.pathogenic_count()) / 
                   std::max(1.0f, static_cast<float>(vcf_parser_.total_variants()));
    yamanaka_sim_.set_mutational_context(tp53_mut, cdkn2a_mut, burden);

    // 3. Evaluar lote inicial con el motor MoE
    if (vcf_parser_.total_variants() > 0) {
        auto decisions = moe_engine_.route_batch(vcf_parser_.get_variants());
        add_log("INFO", "MoE Gating Router processed " + std::to_string(decisions.size()) + 
                        " variants using Apple Accelerate (vDSP AMX) & Metal GPU.");
    }

    // 4. Inicializar ventana genómica y motor CRISPR
    target_sequence_window_ = "CTGCTCCCCCCCTGCCCTCAACAAGATGTTTTGCCAGTTGGCTCTGACTGTACCACCATCCACTACAACT";
    crispr_candidates_ = crispr_engine_.design_guides(target_sequence_window_, 20);
    add_log("INFO", "CRISPR-Cas9 engine synthesized " + std::to_string(crispr_candidates_.size()) + 
                    " guide RNAs with ARM NEON Cas-OFFinder homology screening.");

    // 5. Inicializar Gemelo Digital de Células Vivas
    cell_twin_.initialize_dish(180, tp53_mut);
    add_log("INFO", "Living Cell Digital Twin initialized: 180 cell agents in virtual Petri dish with Metal UMA physics.");

    // 6. Evaluar panel de farmacogenómica y oncología de precisión
    refresh_onco_analysis();
    add_log("INFO", "Oncology Actionability Engine evaluated " + std::to_string(onco_analysis_.detected_biomarkers.size()) + 
                    " actionable drivers, matched " + std::to_string(onco_analysis_.all_recommended_drugs.size()) + " targeted therapies.");
}

void GenomeForgeApp::increase_sim_speed() {
    if (sim_speed_multiplier_ < 32) {
        sim_speed_multiplier_ *= 2;
        add_log("ENGINE", "Simulation real-time speed accelerated to: " + 
                          std::to_string(sim_speed_multiplier_) + "x");
    }
}

void GenomeForgeApp::decrease_sim_speed() {
    if (sim_speed_multiplier_ > 1) {
        sim_speed_multiplier_ /= 2;
        add_log("ENGINE", "Simulation real-time speed decreased to: " + 
                          std::to_string(sim_speed_multiplier_) + "x");
    }
}

void GenomeForgeApp::zoom_in_petri() {
    if (petri_zoom_ < 8.0f) {
        if (petri_zoom_ < 1.9f) petri_zoom_ += 0.5f;
        else if (petri_zoom_ < 3.9f) petri_zoom_ += 1.0f;
        else petri_zoom_ += 2.0f;
        add_log("OPTICS", "Petri dish microscope lens zoomed to: " + 
                          std::to_string(petri_zoom_).substr(0, 3) + "x");
    }
}

void GenomeForgeApp::zoom_out_petri() {
    if (petri_zoom_ > 1.0f) {
        if (petri_zoom_ > 4.0f) petri_zoom_ -= 2.0f;
        else if (petri_zoom_ > 2.0f) petri_zoom_ -= 1.0f;
        else petri_zoom_ = 1.0f;

        if (petri_zoom_ <= 1.0f) {
            petri_pan_x_ = 0.0f;
            petri_pan_y_ = 0.0f;
        }
        add_log("OPTICS", "Petri dish microscope lens zoomed to: " + 
                          std::to_string(petri_zoom_).substr(0, 3) + "x");
    }
}

void GenomeForgeApp::reset_petri_view() {
    petri_zoom_ = 1.0f;
    petri_pan_x_ = 0.0f;
    petri_pan_y_ = 0.0f;
    add_log("OPTICS", "Petri dish microscope view reset to 1.0x wide-field overview.");
}

void GenomeForgeApp::pan_petri(float dx, float dy) {
    float max_pan = 40.0f * (1.0f - 1.0f / petri_zoom_);
    if (max_pan < 0.0f) max_pan = 0.0f;
    petri_pan_x_ = std::clamp(petri_pan_x_ + dx, -max_pan, max_pan);
    petri_pan_y_ = std::clamp(petri_pan_y_ + dy, -max_pan, max_pan);
}

void GenomeForgeApp::save_live_state() {
    genome::ipc::StateInspector::dump_live_state_json(
        "data/live_state.json",
        static_cast<int>(current_tab_),
        vcf_parser_.get_active_locus(),
        telemetry_.get_metrics(),
        vcf_parser_,
        yamanaka_sim_,
        moe_engine_,
        crispr_candidates_,
        crispr_applied_,
        cell_twin_,
        sim_speed_multiplier_,
        petri_zoom_,
        petri_pan_x_,
        petri_pan_y_,
        onco_analysis_.detected_biomarkers.size(),
        onco_analysis_.level1_count
    );
}

void GenomeForgeApp::export_datasets() {
    std::string summary = genome::export_tools::CsvExporter::export_all_datasets(
        "data",
        cell_twin_.get_metrics(),
        cell_twin_.get_cells(),
        vcf_parser_.get_variants(),
        tick_counter_
    );
    add_log("EXPORT", summary);
}

int GenomeForgeApp::run_headless(size_t max_ticks, bool export_csv, const std::string& snapshot_path) {
    std::cout << "╔══════════════════════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║ ⌬ GENOMEFORGE-TUI :: HEADLESS SIMULATION ENGINE                                  ║\n";
    std::cout << "╚══════════════════════════════════════════════════════════════════════════════════╝\n";
    std::cout << "[*] Executing " << max_ticks << " simulation ticks in accelerated headless mode...\n";

    for (size_t i = 0; i < max_ticks; ++i) {
        tick_simulation();
    }

    genome::ipc::StateInspector::dump_live_state_json(
        snapshot_path,
        static_cast<int>(current_tab_),
        vcf_parser_.get_active_locus(),
        telemetry_.get_metrics(),
        vcf_parser_,
        yamanaka_sim_,
        moe_engine_,
        crispr_candidates_,
        crispr_applied_,
        cell_twin_,
        sim_speed_multiplier_,
        petri_zoom_,
        petri_pan_x_,
        petri_pan_y_,
        onco_analysis_.detected_biomarkers.size(),
        onco_analysis_.level1_count
    );

    std::cout << "[+] Telemetry snapshot saved to: " << snapshot_path << "\n";

    if (export_csv) {
        export_datasets();
        std::cout << "[+] Scientific CSV datasets exported to data/ directory.\n";
    }

    const auto& tm = cell_twin_.get_metrics();
    std::cout << "\n--- SIMULATION SUMMARY ---\n"
              << "  Elapsed Bio-Time:    " << tm.simulated_time_str << "\n"
              << "  Total Living Cells:  " << tm.total_cells << "\n"
              << "  Somatic Fibroblasts: " << tm.somatic_count << "\n"
              << "  iPSC Stem Colonies:  " << tm.ipsc_count << "\n"
              << "  Malignant Clones:    " << tm.transformed_count << "\n"
              << "  Culture Confluence:  " << std::fixed << std::setprecision(1) << (tm.confluence * 100.0f) << "%\n"
              << "  Average R-Score:     " << std::setprecision(4) << tm.average_r_score << "\n"
              << "  Actionable Drivers:  " << onco_analysis_.detected_biomarkers.size() << "\n"
              << "  Matched Therapies:   " << onco_analysis_.all_recommended_drugs.size() << " (Level 1 FDA: " << onco_analysis_.level1_count << ")\n"
              << "--------------------------\n";

    return 0;
}

void GenomeForgeApp::tick_simulation() {
    telemetry_.update();
    tick_counter_++;

    if (yamanaka_running_) {
        // Integración de ecuaciones diferenciales ordinarias por Runge-Kutta 4 (RK4)
        yamanaka_sim_.step(0.08f);

        // Despacho continuo en el router MoE
        if (vcf_parser_.total_variants() > 0) {
            size_t idx = static_cast<size_t>(selected_variant_index_) % vcf_parser_.total_variants();
            moe_engine_.route_variant(vcf_parser_.get_variants()[idx]);
        }

        // Avance de física y ciclo celular del Gemelo Digital
        cell_twin_.step(0.05f);
    }

    // Volcado de estado periódico para inspección IA y CLI (ajustado a la velocidad del motor)
    size_t save_interval = std::max<size_t>(2, 2 * sim_speed_multiplier_);
    if (tick_counter_ % save_interval == 0) {
        save_live_state();
    }
}

int GenomeForgeApp::run() {
    auto screen = App::Fullscreen();

    auto main_renderer = Renderer([&] {
        const auto& metrics = telemetry_.get_metrics();
        const auto& variants = vcf_parser_.get_variants();
        const auto& cell_state = yamanaka_sim_.get_state();
        const auto& factors = yamanaka_sim_.get_factors();
        const auto& moe_metrics = moe_engine_.get_metrics();

        // 1. CABECERA Y TELEMETRÍA GLOBAL
        std::string ram_str;
        {
            double used_gb = static_cast<double>(metrics.used_memory_bytes) / (1024.0 * 1024.0 * 1024.0);
            double total_gb = static_cast<double>(metrics.total_memory_bytes) / (1024.0 * 1024.0 * 1024.0);
            std::ostringstream ss;
            ss << std::fixed << std::setprecision(1) << used_gb << " / " << total_gb << " GB (" 
               << static_cast<int>(metrics.memory_used_percent) << "%)";
            ram_str = ss.str();
        }

        auto header_element = vbox({
            hbox({
                text(" ⌬ GENOMEFORGE-TUI :: CYBER-BIOLAB HUD v2.0 ") | bold | color(Color::CyanLight),
                text("⟦ ARCH: APPLE M4 | UMA ZERO-COPY ⟧") | bold | color(Color::White),
                filler(),
                text("AMX: ") | color(Color::GrayLight),
                text("[vDSP::ONLINE] ") | bold | color(Color::GreenLight),
                text("GPU: ") | color(Color::GrayLight),
                text("[METAL::SHARED-UMA] ") | bold | color(Color::GreenLight),
                text("TELEMETRY: ") | color(Color::GrayLight),
                text("[DARWIN-MACH]") | bold | color(Color::MagentaLight),
            }),
            separator(),
            hbox({
                text(" CPU: ") | color(Color::GrayLight),
                text(std::to_string(static_cast<int>(metrics.total_cpu_usage)) + "% ") | bold | color(Color::Yellow),
                gauge(metrics.total_cpu_usage / 100.0f) | size(WIDTH, EQUAL, 10) | color(Color::Yellow),
                text("  P-Cores (" + std::to_string(metrics.p_core_count) + "): ") | color(Color::GrayLight),
                text(std::to_string(static_cast<int>(metrics.p_core_avg_usage)) + "% ") | color(Color::CyanLight),
                gauge(metrics.p_core_avg_usage / 100.0f) | size(WIDTH, EQUAL, 8) | color(Color::CyanLight),
                text("  E-Cores (" + std::to_string(metrics.e_core_count) + "): ") | color(Color::GrayLight),
                text(std::to_string(static_cast<int>(metrics.e_core_avg_usage)) + "% ") | color(Color::BlueLight),
                gauge(metrics.e_core_avg_usage / 100.0f) | size(WIDTH, EQUAL, 8) | color(Color::BlueLight),
                filler(),
                text("UMA Memory: ") | color(Color::GrayLight),
                text(ram_str + " ") | bold | color(Color::MagentaLight),
                gauge(metrics.memory_used_percent / 100.0f) | size(WIDTH, EQUAL, 10) | color(Color::MagentaLight),
            })
        }) | border;

        // 2. BARRA DE PESTAÑAS (TABS CYBER-HUD)
        auto make_tab_button = [&](AppTab tab, const std::string& label, const std::string& key_num) {
            bool is_active = (current_tab_ == tab);
            std::string text_str = " ⟨ " + key_num + ": " + label + " ⟩ ";
            if (is_active) {
                return text(text_str) | bold | color(Color::Black) | bgcolor(Color::CyanLight);
            }
            return text(text_str) | color(Color::GrayLight);
        };

        auto tab_bar = hbox({
            make_tab_button(AppTab::Dashboard, "DASHBOARD", "1"),
            separator(),
            make_tab_button(AppTab::VcfExplorer, "VCF EXPLORER", "2"),
            separator(),
            make_tab_button(AppTab::YamanakaDynamics, "YAMANAKA (OSKM)", "3"),
            separator(),
            make_tab_button(AppTab::MoeRouter, "MoE ROUTER", "4"),
            separator(),
            make_tab_button(AppTab::CrisprStudio, "CRISPR STUDIO", "5"),
            separator(),
            make_tab_button(AppTab::SystemLogs, "SYSTEM & LOGS", "6"),
            separator(),
            make_tab_button(AppTab::DigitalTwin, "DIGITAL TWIN", "7"),
            separator(),
            make_tab_button(AppTab::OncoTherapy, "ONCO THERAPY", "8"),
        }) | border;

        // 3. VISTAS DEDICADAS POR PESTAÑA
        Element active_view;

        // -------------------------------------------------------------
        // TAB 1: DASHBOARD
        // -------------------------------------------------------------
        if (current_tab_ == AppTab::Dashboard) {
            auto patient_summary_card = vbox({
                text("◈ RESUMEN DE MUESTRA Y ESTADO CLÍNICO") | bold | color(Color::Cyan),
                separator(),
                hbox({ text("Archivo VCF:        ") | color(Color::GrayLight), text(vcf_path_) | bold | color(Color::White) }),
                hbox({ text("Total Variantes:    ") | color(Color::GrayLight), text(std::to_string(variants.size())) | bold | color(Color::White) }),
                hbox({ text("Patogénicas ClinVar: ") | color(Color::GrayLight), text(std::to_string(vcf_parser_.pathogenic_count())) | bold | (vcf_parser_.pathogenic_count() > 0 ? color(Color::Red) : color(Color::Green)) }),
                hbox({ text("Colisiones OSKM:    ") | color(Color::GrayLight), text(std::to_string(vcf_parser_.yamanaka_collision_count())) | bold | (vcf_parser_.yamanaka_collision_count() > 0 ? color(Color::Red) : color(Color::Green)) }),
                separator(),
                hbox({ text("Estado Reprogramación: ") | color(Color::GrayLight), 
                       crispr_applied_ ? text("[SAFE - HDR RESCUE COMPLETE]") | bold | color(Color::Green)
                                       : text("[CRITICAL - ONCOGENIC COLLISION]") | bold | color(Color::Red) })
            }) | border | flex;

            Element r_gauge_dash;
            if (cell_state.r_score >= 0.60f) {
                r_gauge_dash = vbox({
                    text("ONCOGENIC INSTABILITY: HIGH RISK") | bold | color(Color::Red),
                    hbox({ gauge(cell_state.r_score) | color(Color::Red) | flex, text(" " + std::to_string(cell_state.r_score).substr(0, 4)) | bold | color(Color::Red) })
                });
            } else {
                r_gauge_dash = vbox({
                    text("ONCOGENIC INSTABILITY: SAFE") | bold | color(Color::Green),
                    hbox({ gauge(cell_state.r_score) | color(Color::Green) | flex, text(" " + std::to_string(cell_state.r_score).substr(0, 4)) | bold | color(Color::Green) })
                });
            }

            auto r_score_card = vbox({
                hbox({
                    text("◈ PULSO YAMANAKA & R-SCORE") | bold | color(Color::Cyan),
                    filler(),
                    text(yamanaka_running_ ? "[RUNNING]" : "[PAUSED]") | bold | (yamanaka_running_ ? color(Color::Green) : color(Color::Yellow))
                }),
                separator(),
                r_gauge_dash,
                separator(),
                hbox({ text("c-Myc (Oncog): ") | size(WIDTH, EQUAL, 16), gauge(factors.cmyc) | color(Color::RedLight) | flex, text(" " + std::to_string(static_cast<int>(factors.cmyc * 100)) + "%") }),
                hbox({ text("p53 Funcional: ") | size(WIDTH, EQUAL, 16), gauge(cell_state.p53_activity) | (cell_state.p53_activity < 0.3f ? color(Color::Red) : color(Color::Green)) | flex, text(" " + std::to_string(static_cast<int>(cell_state.p53_activity * 100)) + "%") }),
                hbox({ text("Pluripotencia: ") | size(WIDTH, EQUAL, 16), gauge(cell_state.pluripotency) | color(Color::GreenLight) | flex, text(" " + std::to_string(static_cast<int>(cell_state.pluripotency * 100)) + "%") }),
            }) | border | flex;

            std::string guide_spacer = crispr_candidates_.empty() ? "N/A" : crispr_candidates_[0].protospacer_20nt;
            std::string guide_pam = crispr_candidates_.empty() ? "N/A" : crispr_candidates_[0].pam;
            float guide_eff = crispr_candidates_.empty() ? 0.0f : crispr_candidates_[0].on_target_efficiency;
            float guide_off = crispr_candidates_.empty() ? 0.0f : crispr_candidates_[0].off_target_safety_score;

            auto crispr_quick_card = vbox({
                text("◈ INTERVENCIÓN CRISPR-CAS9 (MOTOR ACTIVO)") | bold | color(Color::Cyan),
                separator(),
                hbox({ text("Locus Diana: ") | color(Color::GrayLight), text("chr17:7673802 (TP53 Codon 175 Arg/His)") | bold | color(Color::White) }),
                hbox({ text("Guía Óptima: ") | color(Color::GrayLight), text("5'-" + guide_spacer + "-3'") | bold | color(Color::White), text(" (PAM: " + guide_pam + ")") | color(Color::Yellow) }),
                hbox({ text("Eficiencia Doench: ") | color(Color::GrayLight), text(std::to_string(guide_eff).substr(0, 4) + "%") | bold | color(Color::Green),
                       text("  Cas-OFFinder Safety: ") | color(Color::GrayLight), text(std::to_string(guide_off).substr(0, 4) + "/100") | bold | color(Color::Green) }),
                separator(),
                hbox({
                    text("Acciones Rápidas: ") | color(Color::GrayLight),
                    text("[E] Aplicar Edición HDR in-silico  ") | bold | color(Color::Green),
                    text("[R] Reset Genoma  ") | color(Color::Yellow),
                    text("[Espacio] Yamanaka Pulse") | color(Color::Cyan)
                })
            }) | border;

            Elements quick_logs;
            size_t start_idx = (logs_.size() > 4) ? (logs_.size() - 4) : 0;
            for (size_t i = start_idx; i < logs_.size(); ++i) {
                const auto& log = logs_[i];
                Color c = (log.level == "EDIT") ? Color::Green : (log.level == "WARN") ? Color::Yellow : (log.level == "CRIT") ? Color::Red : Color::Cyan;
                quick_logs.push_back(hbox({
                    text("[" + log.timestamp + "] ") | color(Color::GrayDark),
                    text("[" + log.level + "] ") | bold | color(c),
                    text(log.message) | color(Color::White)
                }));
            }
            auto logs_quick_card = vbox({
                text("◈ REGISTRO DE AUDITORÍA RECIENTE") | bold | color(Color::Cyan),
                separator(),
                vbox(std::move(quick_logs))
            }) | border;

            active_view = vbox({
                hbox({ patient_summary_card, r_score_card }) | flex,
                crispr_quick_card,
                logs_quick_card
            }) | flex;
        }

        // -------------------------------------------------------------
        // TAB 2: VCF EXPLORER
        // -------------------------------------------------------------
        else if (current_tab_ == AppTab::VcfExplorer) {
            const auto* cur_locus = genome::vcf::LocusCatalog::get_by_index(active_locus_index_);

            Elements loci_pills;
            loci_pills.push_back(text(" ◈ CATÁLOGO 1000 GENOMAS: ") | bold | color(Color::Cyan));
            for (size_t li = 0; li < genome::vcf::LocusCatalog::total_loci(); ++li) {
                const auto* l = genome::vcf::LocusCatalog::get_by_index(li);
                bool is_cur = (li == active_locus_index_);
                std::string gene_txt = " " + l->gene_symbol + " ";
                if (is_cur) {
                    loci_pills.push_back(text(gene_txt) | bold | color(Color::Black) | bgcolor(Color::Cyan));
                } else {
                    loci_pills.push_back(text(gene_txt) | color(Color::GrayLight));
                }
            }
            loci_pills.push_back(filler());
            loci_pills.push_back(text("[L] CONMUTAR GEN (Tabix Range)  ") | bold | color(Color::Yellow));

            auto locus_header_bar = vbox({
                hbox({
                    text(" ◈ MUESTRA ACTIVA: ") | bold | color(Color::Cyan),
                    text(healthy_reference_mode_ ? "⟦ 1000 GENOMAS HG00096 (CONTROL SANO - 0 COLISIONES) ⟧" : "⟦ PACIENTE CLÍNICO (MUTACIONES ONCOGÉNICAS TP53/CDKN2A) ⟧") | bold | (healthy_reference_mode_ ? color(Color::Green) : color(Color::Yellow)),
                    filler(),
                    text("[M] CONMUTAR MUESTRA (Paciente <-> 1000G Sano)  ") | bold | color(Color::GreenLight)
                }),
                separator(),
                hbox(std::move(loci_pills)),
                cur_locus ? hbox({
                    text("   Locus Activo: ") | color(Color::GrayLight),
                    text(cur_locus->gene_symbol + " - " + cur_locus->full_name) | bold | color(Color::White),
                    text(" [" + cur_locus->chromosome + ":" + std::to_string(cur_locus->start_pos) + "-" + std::to_string(cur_locus->end_pos) + "]") | color(Color::Cyan),
                    filler(),
                    text("1000 Genomes Phase 3 (GRCh38) ") | color(Color::GrayDark)
                }) : text("")
            }) | border;

            Elements table_rows;
            table_rows.push_back(
                hbox({
                    text("   CHR     POS         GENE       ALLELE    AF     DP    IMPACT / CLNSIG STATUS") | bold | color(Color::White)
                })
            );
            table_rows.push_back(separator());

            for (size_t i = 0; i < variants.size(); ++i) {
                const auto& v = variants[i];
                bool is_selected = (static_cast<int>(i) == selected_variant_index_);
                std::string prefix = is_selected ? " > " : "   ";

                std::string chr_s = "chr" + std::to_string(v.chrom_id);
                if (v.chrom_id == 23) chr_s = "chrX";
                if (v.chrom_id == 24) chr_s = "chrY";

                std::ostringstream pos_ss, gene_ss, allele_ss, af_ss, dp_ss;
                pos_ss << std::left << std::setw(11) << v.pos;
                gene_ss << std::left << std::setw(10) << v.gene_symbol;
                std::string alleles = std::string(v.ref) + " > " + std::string(v.alt);
                allele_ss << std::left << std::setw(9) << alleles;
                af_ss << std::fixed << std::setprecision(2) << std::setw(6) << v.allele_frequency;
                dp_ss << std::setw(5) << v.depth;

                std::string status_label;
                Color status_color = Color::White;
                switch (v.effect) {
                    case genome::vcf::VariantEffect::YamanakaCollision:
                        status_label = "[!] YAMANAKA COLLISION";
                        status_color = Color::Red;
                        break;
                    case genome::vcf::VariantEffect::Pathogenic:
                        status_label = "[x] PATHOGENIC";
                        status_color = Color::RedLight;
                        break;
                    case genome::vcf::VariantEffect::LikelyPathogenic:
                        status_label = "[?] LIKELY PATHOGENIC";
                        status_color = Color::Yellow;
                        break;
                    case genome::vcf::VariantEffect::Benign:
                        status_label = "[+] STABLE / WILD-TYPE";
                        status_color = Color::Green;
                        break;
                    default:
                        status_label = "[-] VUS / UNCERTAIN";
                        status_color = Color::GrayLight;
                        break;
                }

                auto row = hbox({
                    text(prefix) | bold | (is_selected ? color(Color::Cyan) : color(Color::GrayDark)),
                    text(chr_s) | size(WIDTH, EQUAL, 8),
                    text(pos_ss.str()) | size(WIDTH, EQUAL, 12),
                    text(gene_ss.str()) | bold | size(WIDTH, EQUAL, 11),
                    text(allele_ss.str()) | size(WIDTH, EQUAL, 10),
                    text(af_ss.str()) | size(WIDTH, EQUAL, 7),
                    text(dp_ss.str()) | size(WIDTH, EQUAL, 6),
                    text(status_label) | bold | color(status_color),
                });

                if (is_selected) row = row | inverted;
                table_rows.push_back(row);
            }

            Element locus_box = text("No variant selected");
            if (selected_variant_index_ >= 0 && selected_variant_index_ < static_cast<int>(variants.size())) {
                const auto& sv = variants[selected_variant_index_];
                bool is_tp53 = (std::string(sv.gene_symbol) == "TP53");

                locus_box = vbox({
                    text("◈ INSPECTOR DE LOCUS Y VENTANA DE SECUENCIA GENÓMICA") | bold | color(Color::Cyan),
                    separator(),
                    hbox({
                        text("Coordenada: ") | color(Color::GrayLight),
                        text("chr" + std::to_string(sv.chrom_id) + ":" + std::to_string(sv.pos)) | bold | color(Color::Cyan),
                        text(" | Gen: ") | color(Color::GrayLight),
                        text(sv.gene_symbol) | bold | color(Color::Yellow),
                        text(" | Alelos: ") | color(Color::GrayLight),
                        text(std::string(sv.ref) + " > " + std::string(sv.alt)) | bold | color(Color::White),
                        text(" | Cobertura DP: ") | color(Color::GrayLight),
                        text(std::to_string(sv.depth) + "x") | color(Color::White),
                        text(" | AF: ") | color(Color::GrayLight),
                        text(std::to_string(sv.allele_frequency).substr(0, 4)) | color(Color::White),
                    }),
                    separator(),
                    text("Secuencia Genómica Flanqueante (+/- 20 bp):") | color(Color::GrayLight),
                    hbox({
                        text(" 5'- C T G C T C C C C C C C [") | color(Color::White),
                        text(std::string(sv.alt)) | bold | (crispr_applied_ ? color(Color::Green) : color(Color::Red)),
                        text(" / ") | color(Color::GrayLight),
                        text(std::string(sv.ref)) | color(Color::Green),
                        text("] C C C T C A A C A A G A T G T T T T G C - 3'") | color(Color::White),
                    }),
                    is_tp53 ? (
                        crispr_applied_ ?
                            text(" [+] STATUS: CORREGIDO POR HDR. PROTEÍNA p53 CONFORMACIONALMENTE SALVAJE (Arg175).") | bold | color(Color::Green) :
                            text(" [!] STATUS: MUTACIÓN EN DOMINIO DE UNIÓN A ADN (R175H) INHIBE LA RESPUESTA DE SENESCENCIA.") | bold | color(Color::Red)
                    ) : text(" Variante no asociada a colisión directa en supresores tumorales críticos.") | color(Color::GrayLight)
                }) | border;
            }

            active_view = vbox({
                locus_header_bar,
                vbox(std::move(table_rows)) | border | flex,
                locus_box
            }) | flex;
        }

        // -------------------------------------------------------------
        // TAB 3: YAMANAKA DYNAMICS (OSKM)
        // -------------------------------------------------------------
        else if (current_tab_ == AppTab::YamanakaDynamics) {
            std::string sparkline = "";
            for (float val : yamanaka_sim_.get_r_score_history()) {
                if (val < 0.25f) sparkline += " ";
                else if (val < 0.40f) sparkline += "▃";
                else if (val < 0.55f) sparkline += "▄";
                else if (val < 0.70f) sparkline += "▅";
                else if (val < 0.85f) sparkline += "▆";
                else sparkline += "█";
            }

            auto ode_equations_box = vbox({
                text("◈ FORMULACIÓN DINÁMICA DE REDES GÉNICA (ODEs CONTINUAS CON INTEGRACIÓN RK4)") | bold | color(Color::Cyan),
                separator(),
                text(" d[Pluripotencia]/dt = k1·[Oct4][Sox2] + k2·[Klf4] - γ_p·[Pluri]") | color(Color::White),
                text(" d[Senescencia]/dt   = k3·[c-Myc]·(1 - Ω_TP53) + k4·[Klf4] - γ_s·[Senescencia]") | color(Color::White),
                text(" d[Inestabilidad]/dt = k_myc·[c-Myc]^2·(1 - p53_act) + k_dna·Burden - γ_i·[Inestabilidad]") | color(Color::White),
                separator(),
                hbox({
                    text("Tiempo de Simulación (t): ") | color(Color::GrayLight),
                    text(std::to_string(yamanaka_sim_.get_time()).substr(0, 5) + " h") | bold | color(Color::White),
                    text("  | Estado Genético p53: ") | color(Color::GrayLight),
                    crispr_applied_ ? text("[SALVAJE / RESCATADO]") | bold | color(Color::Green)
                                    : text("[MUTADO / DOMINANTE NEGATIVO]") | bold | color(Color::Red)
                })
            }) | border;

            auto factors_box = vbox({
                text("◈ NIVELES DE EXPRESIÓN DE FACTORES OSKM") | bold | color(Color::Cyan),
                separator(),
                hbox({ text(" Oct4 (POU5F1): ") | size(WIDTH, EQUAL, 18), gauge(factors.oct4) | color(Color::Cyan) | flex, text(" " + std::to_string(static_cast<int>(factors.oct4 * 100)) + "%") }),
                hbox({ text(" Sox2:          ") | size(WIDTH, EQUAL, 18), gauge(factors.sox2) | color(Color::Cyan) | flex, text(" " + std::to_string(static_cast<int>(factors.sox2 * 100)) + "%") }),
                hbox({ text(" Klf4:          ") | size(WIDTH, EQUAL, 18), gauge(factors.klf4) | color(Color::Yellow) | flex, text(" " + std::to_string(static_cast<int>(factors.klf4 * 100)) + "%") }),
                hbox({ text(" c-Myc (Oncog): ") | size(WIDTH, EQUAL, 18), gauge(factors.cmyc) | color(Color::RedLight) | flex, text(" " + std::to_string(static_cast<int>(factors.cmyc * 100)) + "%") }),
                separator(),
                text("◈ RESPUESTA DOWNSTREAM (CHECKPOINTS):") | bold | color(Color::White),
                hbox({ text(" p53 Funcional: ") | size(WIDTH, EQUAL, 18), gauge(cell_state.p53_activity) | (cell_state.p53_activity < 0.3f ? color(Color::Red) : color(Color::Green)) | flex, text(" " + std::to_string(static_cast<int>(cell_state.p53_activity * 100)) + "%") }),
                hbox({ text(" Freno p21/CDKN1A: ") | size(WIDTH, EQUAL, 18), gauge(cell_state.p21_activity) | (cell_state.p21_activity < 0.3f ? color(Color::Red) : color(Color::Green)) | flex, text(" " + std::to_string(static_cast<int>(cell_state.p21_activity * 100)) + "%") }),
            }) | border | flex;

            auto r_score_timeline_box = vbox({
                text("◈ EVOLUCIÓN TEMPORAL DEL R-SCORE (INSTABILIDAD)") | bold | color(Color::Cyan),
                separator(),
                hbox({
                    text("R-Score Actual: ") | color(Color::GrayLight),
                    text(std::to_string(cell_state.r_score).substr(0, 4) + " ") | bold | (cell_state.r_score >= 0.60f ? color(Color::Red) : color(Color::Green)),
                    (cell_state.r_score >= 0.60f ? text("[ALERTA: RIESGO ONCOGÉNICO ELEVADO]") | bold | color(Color::Red)
                                                 : text("[ESTADO FISIOLÓGICO SEGURO]") | bold | color(Color::Green))
                }),
                separator(),
                text("Trayectoria temporal calculada con Runge-Kutta:") | color(Color::GrayLight),
                text(sparkline) | bold | (cell_state.r_score >= 0.60f ? color(Color::Red) : color(Color::Green)),
                separator(),
                text("Controles: [Espacio] Pausar/Reanudar  |  [E] Aplicar CRISPR  |  [R] Reset Genoma") | color(Color::GrayLight)
            }) | border | flex;

            active_view = vbox({
                ode_equations_box,
                hbox({ factors_box, r_score_timeline_box }) | flex
            }) | flex;
        }

        // -------------------------------------------------------------
        // TAB 4: MoE ROUTER
        // -------------------------------------------------------------
        else if (current_tab_ == AppTab::MoeRouter) {
            // Obtener decisión de enrutamiento para la variante actualmente seleccionada
            genome::moe::GatingDecision current_decision{};
            if (!variants.empty() && selected_variant_index_ < static_cast<int>(variants.size())) {
                current_decision = moe_engine_.route_variant(variants[selected_variant_index_]);
            }

            auto moe_arch_box = vbox({
                text("◈ ARQUITECTURA MIXTURE OF EXPERTS (MoE) & ENRUTADOR POR HARDWARE") | bold | color(Color::Cyan),
                separator(),
                text(" Gating Network: W_g · x + b_g  →  Softmax(Top-2)  →  Dispersión de Carga") | color(Color::White),
                text(" Motor CPU: Apple Accelerate (cblas_sgemv / vDSP_mmul sobre matrices de 16D en AMX)") | color(Color::Green),
                text(" Motor GPU: Metal Compute Pipeline en Memoria Unificada (Zero-Copy MTLResourceStorageModeShared)") | color(Color::Green),
                separator(),
                hbox({
                    text("Inferencias Totales: ") | color(Color::GrayLight), text(std::to_string(moe_metrics.total_inferences)) | bold | color(Color::White),
                    text(" | Latencia Media vDSP: ") | color(Color::GrayLight), text(std::to_string(moe_metrics.avg_vdsp_latency_us).substr(0, 4) + " µs") | bold | color(Color::Cyan),
                    text(" | Latencia Metal GPU: ") | color(Color::GrayLight), text(std::to_string(moe_metrics.avg_metal_latency_us).substr(0, 4) + " µs") | bold | color(Color::Magenta),
                })
            }) | border;

            auto experts_detail_box = vbox({
                text("◈ ACTIVACIÓN DE SUBREDES DE EXPERTOS BIOLÓGICOS") | bold | color(Color::Cyan),
                separator(),
                hbox({ text(" [1] Daño ADN / ATM:   ") | size(WIDTH, EQUAL, 24), gauge(moe_metrics.expert_activation_levels[0]) | color(Color::BlueLight) | flex, text(" " + std::to_string(static_cast<int>(moe_metrics.expert_activation_levels[0] * 100)) + "% [vDSP]") | color(Color::GrayLight) }),
                hbox({ text(" [2] Ciclo Celular / p21:") | size(WIDTH, EQUAL, 24), gauge(moe_metrics.expert_activation_levels[1]) | color(Color::Yellow) | flex, text(" " + std::to_string(static_cast<int>(moe_metrics.expert_activation_levels[1] * 100)) + "% [vDSP]") | color(Color::GrayLight) }),
                hbox({ text(" [3] Apoptosis p53:     ") | size(WIDTH, EQUAL, 24), gauge(moe_metrics.expert_activation_levels[2]) | (moe_metrics.expert_activation_levels[2] < 0.3f ? color(Color::Red) : color(Color::Green)) | flex, text(" " + std::to_string(static_cast<int>(moe_metrics.expert_activation_levels[2] * 100)) + "% [vDSP]") | color(Color::GrayLight) }),
                hbox({ text(" [4] Pluripotencia GRN: ") | size(WIDTH, EQUAL, 24), gauge(moe_metrics.expert_activation_levels[3]) | color(Color::GreenLight) | flex, text(" " + std::to_string(static_cast<int>(moe_metrics.expert_activation_levels[3] * 100)) + "% [Metal]") | color(Color::GrayLight) }),
            }) | border | flex;

            std::string top1_name = (current_decision.top1_expert == 0) ? "Daño ADN / Reparación" :
                                    (current_decision.top1_expert == 1) ? "Ciclo Celular / p21" :
                                    (current_decision.top1_expert == 2) ? "Apoptosis p53" : "Pluripotencia GRN";
            std::string top2_name = (current_decision.top2_expert == 0) ? "Daño ADN / Reparación" :
                                    (current_decision.top2_expert == 1) ? "Ciclo Celular / p21" :
                                    (current_decision.top2_expert == 2) ? "Apoptosis p53" : "Pluripotencia GRN";

            auto routing_decision_box = vbox({
                text("◈ ENRUTAMIENTO DE LA VARIANTE SELECCIONADA") | bold | color(Color::Cyan),
                separator(),
                hbox({ text("Top-1 Expert: ") | color(Color::GrayLight), text(top1_name) | bold | color(Color::Yellow), text(" (Softmax: " + std::to_string(current_decision.top1_weight).substr(0, 4) + ")") }),
                hbox({ text("Top-2 Expert: ") | color(Color::GrayLight), text(top2_name) | bold | color(Color::Cyan), text(" (Softmax: " + std::to_string(current_decision.top2_weight).substr(0, 4) + ")") }),
                separator(),
                hbox({ text("Acelerador Utilizado: ") | color(Color::GrayLight), current_decision.dispatched_to_gpu ? text("GPU Metal (Zero-Copy UMA)") | bold | color(Color::Magenta) : text("CPU Apple Accelerate vDSP/AMX") | bold | color(Color::Green) }),
                hbox({ text("Latencia de Inferencia: ") | color(Color::GrayLight), text(std::to_string(current_decision.execution_latency_us).substr(0, 4) + " µs") | color(Color::White) }),
            }) | border | flex;

            active_view = vbox({
                moe_arch_box,
                hbox({ experts_detail_box, routing_decision_box }) | flex
            }) | flex;
        }

        // -------------------------------------------------------------
        // TAB 5: CRISPR STUDIO
        // -------------------------------------------------------------
        else if (current_tab_ == AppTab::CrisprStudio) {
            Elements candidate_rows;
            candidate_rows.push_back(
                hbox({
                    text("   #   PAM   OFFSET   SPACER (20 nt)           EFICACIA   OFF-TARGET SAFETY   RECOMENDACIÓN") | bold | color(Color::White)
                })
            );
            candidate_rows.push_back(separator());

            for (size_t i = 0; i < crispr_candidates_.size(); ++i) {
                const auto& cand = crispr_candidates_[i];
                bool is_sel = (static_cast<int>(i) == selected_guide_index_);
                std::string prefix = is_sel ? " > " : "   ";

                std::ostringstream off_ss;
                off_ss << (cand.pam_position_offset >= 0 ? "+" : "") << cand.pam_position_offset << " bp";

                auto row = hbox({
                    text(prefix) | bold | (is_sel ? color(Color::Cyan) : color(Color::GrayDark)),
                    text("[" + std::to_string(i + 1) + "]") | size(WIDTH, EQUAL, 4),
                    text(cand.pam) | bold | color(Color::Yellow) | size(WIDTH, EQUAL, 6),
                    text(off_ss.str()) | size(WIDTH, EQUAL, 9),
                    text("5'-" + cand.protospacer_20nt + "-3'") | bold | color(Color::White) | size(WIDTH, EQUAL, 28),
                    text(std::to_string(cand.on_target_efficiency).substr(0, 4) + "%") | color(Color::Green) | size(WIDTH, EQUAL, 11),
                    text(std::to_string(cand.off_target_safety_score).substr(0, 4) + "/100") | color(Color::Green) | size(WIDTH, EQUAL, 20),
                    text(cand.recommendation) | bold | (cand.on_target_efficiency >= 90.0f ? color(Color::Green) : color(Color::GrayLight)),
                });

                if (is_sel) row = row | inverted;
                candidate_rows.push_back(row);
            }

            auto crispr_table_box = vbox({
                text("◈ GUÍAS gRNA GENERADAS MEDIANTE BÚSQUEDA PAM NGG Y SCORING CAS-OFFINDER (ARM NEON)") | bold | color(Color::Cyan),
                separator(),
                vbox(std::move(candidate_rows))
            }) | border;

            auto edit_action_box = vbox({
                text("◈ MOTOR DE EDICIÓN GENÓMICA VIRTUAL (SIMULACIÓN HDR)") | bold | color(Color::Cyan),
                separator(),
                crispr_applied_ ?
                    vbox({
                        text("[+] CORRECCIÓN HDR APLICADA CON ÉXITO") | bold | color(Color::Green),
                        text("El codón mutado en TP53 (175H) ha sido sustituido virtualmente por el alelo salvaje (175R)."),
                        text("El R-Score ha descendido a nivel seguro en el simulador RK4 de Yamanaka."),
                    }) :
                    vbox({
                        text("[!] MUTACIÓN PATOGÉNICA DETECTADA EN TP53") | bold | color(Color::Red),
                        text("La sobreexpresión de c-Myc sin p53 funcional provocará transformación oncogénica incontrolada."),
                        text("Presione [E] para aplicar la corrección por Homology-Directed Repair."),
                    }),
                separator(),
                hbox({
                    text("Controles: ") | color(Color::GrayLight),
                    text("[E] Aplicar Corrección HDR  ") | bold | color(Color::Green),
                    text("[R] Restablecer Genoma Mutado  ") | bold | color(Color::Yellow),
                    text("[↑/↓] Seleccionar Guía  ") | color(Color::Cyan),
                })
            }) | border | flex;

            active_view = vbox({
                crispr_table_box,
                edit_action_box
            }) | flex;
        }

        // -------------------------------------------------------------
        // TAB 6: SYSTEM & LOGS
        // -------------------------------------------------------------
        else if (current_tab_ == AppTab::SystemLogs) {
            Elements core_gauges;
            for (const auto& core : metrics.core_metrics) {
                std::string label = (core.is_performance_core ? "P-Core " : "E-Core ") + std::to_string(core.core_id) + ": ";
                Color bar_color = core.is_performance_core ? Color::Cyan : Color::Blue;
                core_gauges.push_back(
                    hbox({
                        text(label) | size(WIDTH, EQUAL, 11),
                        gauge(core.usage_percent / 100.0f) | color(bar_color) | size(WIDTH, EQUAL, 12),
                        text(" " + std::to_string(static_cast<int>(core.usage_percent)) + "%")
                    })
                );
            }

            auto hardware_detailed_box = vbox({
                text("◈ TELEMETRÍA DETALLADA DE APPLE SILICON (DARWIN MACH KERNEL)") | bold | color(Color::Cyan),
                separator(),
                hbox({
                    text("Chipset: ") | color(Color::GrayLight), text(metrics.chip_model) | bold | color(Color::White),
                    text(" | P-Cores: ") | color(Color::GrayLight), text(std::to_string(metrics.p_core_count)) | color(Color::Cyan),
                    text(" | E-Cores: ") | color(Color::GrayLight), text(std::to_string(metrics.e_core_count)) | color(Color::Blue),
                }),
                separator(),
                text("Distribución de Carga por Núcleo:") | color(Color::GrayLight),
                vbox(std::move(core_gauges)),
                separator(),
                text("Memoria Unificada (UMA): " + ram_str) | bold | color(Color::Magenta)
            }) | border | flex;

            Elements full_logs;
            for (const auto& log : logs_) {
                Color lvl_color = (log.level == "EDIT") ? Color::Green : (log.level == "WARN") ? Color::Yellow : (log.level == "CRIT") ? Color::Red : Color::Cyan;
                full_logs.push_back(
                    hbox({
                        text("[" + log.timestamp + "] ") | color(Color::GrayDark),
                        text("[" + log.level + "] ") | bold | color(lvl_color),
                        text(log.message) | color(Color::White)
                    })
                );
            }

            auto logs_full_box = vbox({
                text("◈ REGISTRO HISTÓRICO DE AUDITORÍA Y EVENTOS BIOLÓGICOS") | bold | color(Color::Cyan),
                separator(),
                vbox(std::move(full_logs)) | yflex
            }) | border | flex;

            active_view = hbox({
                hardware_detailed_box,
                logs_full_box
            }) | flex;
        }

        // -------------------------------------------------------------
        // TAB 7: DIGITAL TWIN (PETRI DISH & AGENT-BASED SIMULATION)
        // -------------------------------------------------------------
        else if (current_tab_ == AppTab::DigitalTwin) {
            const auto& t_metrics = cell_twin_.get_metrics();
            const auto& cells = cell_twin_.get_cells();

            // Canvas para la placa de Petri (114x70 pixels con sub-píxeles Braille)
            int c_width = 114;
            int c_height = 70;
            Canvas petri_canvas(c_width, c_height);

            // Factores de escala con Zoom óptico del microscopio virtual
            float scale_x = (((c_width / 2.0f) - 6.0f) / 45.0f) * petri_zoom_;
            float scale_y = (((c_height / 2.0f) - 5.0f) / 45.0f) * petri_zoom_;

            // Borde exterior, bisel e interior de la placa de Petri proyectada con Pan y Zoom
            int dish_cx = static_cast<int>((c_width / 2.0f) - petri_pan_x_ * scale_x);
            int dish_cy = static_cast<int>((c_height / 2.0f) - petri_pan_y_ * scale_y);
            int dish_rx = static_cast<int>(45.0f * scale_x);
            int dish_ry = static_cast<int>(45.0f * scale_y);

            if (dish_rx > 2 && dish_ry > 2) {
                petri_canvas.DrawPointEllipse(dish_cx, dish_cy, dish_rx, dish_ry, Color::Cyan);
                if (dish_rx > 4 && dish_ry > 4) {
                    petri_canvas.DrawPointEllipse(dish_cx, dish_cy, dish_rx - 2, dish_ry - 2, Color::GrayDark);
                }
            }
            
            // Retícula óptica central y guías micrométricas del ocular del microscopio
            int cx = c_width / 2;
            int cy = c_height / 2;
            petri_canvas.DrawPointLine(cx - 4, cy, cx + 4, cy, Color::GrayDark);
            petri_canvas.DrawPointLine(cx, cy - 3, cx, cy + 3, Color::GrayDark);
            if (petri_zoom_ > 1.05f) {
                // Ticks micrométricos de alta magnificación
                petri_canvas.DrawPoint(cx - 12, cy, true, Color::Cyan);
                petri_canvas.DrawPoint(cx + 12, cy, true, Color::Cyan);
                petri_canvas.DrawPoint(cx, cy - 8, true, Color::Cyan);
                petri_canvas.DrawPoint(cx, cy + 8, true, Color::Cyan);
            }

            // Renderizado citomorfológico avanzado de cada célula según su fenotipo, zoom y división
            bool high_mag = (petri_zoom_ >= 2.0f);

            for (const auto& cell : cells) {
                int px = static_cast<int>((c_width / 2.0f) + (cell.x - petri_pan_x_) * scale_x);
                int py = static_cast<int>((c_height / 2.0f) + (cell.y - petri_pan_y_) * scale_y);

                if (px < 2 || px >= c_width - 2 || py < 2 || py >= c_height - 2) continue;

                // 1. Células en mitosis activa con surco de clivaje / citocinesis
                if (cell.cleavage_state > 0) {
                    if (high_mag) {
                        // Alta magnificación: huso acromático y cromátidas segregadas
                        if (cell.cleavage_state == 1) {
                            petri_canvas.DrawPoint(px - 2, py, true, Color::MagentaLight);
                            petri_canvas.DrawPoint(px - 1, py, true, Color::Magenta);
                            petri_canvas.DrawPoint(px, py, true, Color::White);
                            petri_canvas.DrawPoint(px + 1, py, true, Color::Magenta);
                            petri_canvas.DrawPoint(px + 2, py, true, Color::MagentaLight);
                            petri_canvas.DrawPoint(px - 1, py - 1, true, Color::MagentaLight);
                            petri_canvas.DrawPoint(px + 1, py + 1, true, Color::MagentaLight);
                        } else {
                            petri_canvas.DrawPoint(px - 3, py, true, Color::MagentaLight);
                            petri_canvas.DrawPoint(px - 2, py, true, Color::MagentaLight);
                            petri_canvas.DrawPoint(px - 2, py - 1, true, Color::Magenta);
                            petri_canvas.DrawPoint(px, py, true, Color::White);
                            petri_canvas.DrawPoint(px + 2, py, true, Color::MagentaLight);
                            petri_canvas.DrawPoint(px + 3, py, true, Color::MagentaLight);
                            petri_canvas.DrawPoint(px + 2, py + 1, true, Color::Magenta);
                        }
                    } else {
                        if (cell.cleavage_state == 1) {
                            petri_canvas.DrawPoint(px - 1, py, true, Color::MagentaLight);
                            petri_canvas.DrawPoint(px + 1, py, true, Color::MagentaLight);
                            petri_canvas.DrawPoint(px, py, true, Color::White);
                        } else {
                            petri_canvas.DrawPoint(px - 2, py, true, Color::MagentaLight);
                            petri_canvas.DrawPoint(px - 1, py, true, Color::MagentaLight);
                            petri_canvas.DrawPoint(px + 1, py, true, Color::MagentaLight);
                            petri_canvas.DrawPoint(px + 2, py, true, Color::MagentaLight);
                        }
                    }
                    continue;
                }

                // 2. Morfología según fenotipo celular
                switch (cell.phenotype) {
                    case genome::twin::CellPhenotype::SomaticFibroblast: {
                        float v_len = std::hypot(cell.vx, cell.vy);
                        int dx = 1, dy = 0;
                        if (v_len > 0.1f) {
                            dx = static_cast<int>(std::round((cell.vx / v_len) * 1.5f));
                            dy = static_cast<int>(std::round((cell.vy / v_len) * 1.2f));
                        } else {
                            int angle_idx = (cell.cell_id % 4);
                            if (angle_idx == 0) { dx = 1; dy = 0; }
                            else if (angle_idx == 1) { dx = 1; dy = 1; }
                            else if (angle_idx == 2) { dx = 0; dy = 1; }
                            else { dx = -1; dy = 1; }
                        }
                        if (high_mag) {
                            // Fibroblasto en alta magnificación: cuerpo citoplasmático y lamelipodios
                            petri_canvas.DrawPoint(px, py, true, Color::Cyan);
                            petri_canvas.DrawPoint(px, py - 1, true, Color::BlueLight);
                            petri_canvas.DrawPoint(px + dx, py + dy, true, Color::BlueLight);
                            petri_canvas.DrawPoint(px - dx, py - dy, true, Color::BlueLight);
                            petri_canvas.DrawPoint(px + 2 * dx, py + 2 * dy, true, Color::Blue);
                            petri_canvas.DrawPoint(px - 2 * dx, py - 2 * dy, true, Color::Blue);
                        } else {
                            petri_canvas.DrawPoint(px, py, true, Color::Cyan);
                            petri_canvas.DrawPoint(px + dx, py + dy, true, Color::BlueLight);
                            petri_canvas.DrawPoint(px - dx, py - dy, true, Color::BlueLight);
                        }
                        break;
                    }

                    case genome::twin::CellPhenotype::IntermediateReprog: {
                        if (high_mag) {
                            // MET Intermedio ampliado: condensación cromatínica periférica
                            petri_canvas.DrawPoint(px, py, true, Color::YellowLight);
                            petri_canvas.DrawPoint(px + 1, py, true, Color::Yellow);
                            petri_canvas.DrawPoint(px - 1, py, true, Color::Yellow);
                            petri_canvas.DrawPoint(px, py + 1, true, Color::Yellow);
                            petri_canvas.DrawPoint(px, py - 1, true, Color::Yellow);
                        } else {
                            petri_canvas.DrawPoint(px, py, true, Color::YellowLight);
                            petri_canvas.DrawPoint(px + 1, py, true, Color::Yellow);
                            petri_canvas.DrawPoint(px, py + 1, true, Color::Yellow);
                        }
                        break;
                    }

                    case genome::twin::CellPhenotype::PluripotentIPSC: {
                        if (high_mag) {
                            // Colonia iPSC ampliada: núcleo con nucleolo prominente y borde epitelial
                            petri_canvas.DrawPoint(px, py, true, Color::White);
                            petri_canvas.DrawPoint(px + 1, py, true, Color::GreenLight);
                            petri_canvas.DrawPoint(px - 1, py, true, Color::GreenLight);
                            petri_canvas.DrawPoint(px, py + 1, true, Color::GreenLight);
                            petri_canvas.DrawPoint(px, py - 1, true, Color::GreenLight);
                            petri_canvas.DrawPoint(px + 1, py + 1, true, Color::Green);
                            petri_canvas.DrawPoint(px - 1, py - 1, true, Color::Green);
                            petri_canvas.DrawPoint(px + 1, py - 1, true, Color::Green);
                            petri_canvas.DrawPoint(px - 1, py + 1, true, Color::Green);
                        } else {
                            petri_canvas.DrawPoint(px, py, true, Color::GreenLight);
                            petri_canvas.DrawPoint(px + 1, py, true, Color::Green);
                            petri_canvas.DrawPoint(px - 1, py, true, Color::Green);
                            petri_canvas.DrawPoint(px, py + 1, true, Color::Green);
                            petri_canvas.DrawPoint(px, py - 1, true, Color::Green);
                        }
                        break;
                    }

                    case genome::twin::CellPhenotype::TransformedMalignant: {
                        if (high_mag) {
                            // Neoplasia Maligna ampliada: núcleo gigante bizarro, multinucleación
                            petri_canvas.DrawPoint(px, py, true, Color::White);
                            petri_canvas.DrawPoint(px + 1, py, true, Color::RedLight);
                            petri_canvas.DrawPoint(px - 1, py, true, Color::RedLight);
                            petri_canvas.DrawPoint(px, py + 1, true, Color::RedLight);
                            petri_canvas.DrawPoint(px, py - 1, true, Color::RedLight);
                            petri_canvas.DrawPoint(px + 2, py, true, Color::Red);
                            petri_canvas.DrawPoint(px - 2, py, true, Color::Red);
                            petri_canvas.DrawPoint(px + 1, py + 1, true, Color::MagentaLight);
                            petri_canvas.DrawPoint(px - 1, py - 1, true, Color::MagentaLight);
                        } else {
                            petri_canvas.DrawPoint(px, py, true, Color::RedLight);
                            petri_canvas.DrawPoint(px + 1, py, true, Color::Red);
                            petri_canvas.DrawPoint(px - 1, py, true, Color::Red);
                            petri_canvas.DrawPoint(px, py + 1, true, Color::Red);
                            petri_canvas.DrawPoint(px + 1, py + 1, true, Color::RedLight);
                            petri_canvas.DrawPoint(px - 1, py - 1, true, Color::Red);
                        }
                        break;
                    }

                    case genome::twin::CellPhenotype::ApoptoticCorpse: {
                        if (high_mag) {
                            // Apoptosis ampliada: dispersión de múltiples cuerpos apoptóticos
                            petri_canvas.DrawPoint(px, py, true, Color::GrayLight);
                            petri_canvas.DrawPoint(px + 2, py - 1, true, Color::GrayDark);
                            petri_canvas.DrawPoint(px - 2, py + 1, true, Color::GrayDark);
                            petri_canvas.DrawPoint(px + 1, py + 2, true, Color::GrayDark);
                            petri_canvas.DrawPoint(px - 1, py - 2, true, Color::GrayDark);
                        } else {
                            petri_canvas.DrawPoint(px, py, true, Color::GrayLight);
                            petri_canvas.DrawPoint(px + 1, py - 1, true, Color::GrayDark);
                            petri_canvas.DrawPoint(px - 1, py + 1, true, Color::GrayDark);
                        }
                        break;
                    }
                }
            }

            std::string zoom_str = std::to_string(static_cast<int>(petri_zoom_)) + "." + 
                                   std::to_string(static_cast<int>(petri_zoom_ * 10) % 10) + "x";
            std::string speed_str = std::to_string(sim_speed_multiplier_) + "x";

            auto petri_box = vbox({
                hbox({
                    text("⌬ CYBER-INCUBATOR [CHAMBER 01]") | bold | color(Color::Cyan),
                    text(" │ ") | color(Color::GrayDark),
                    text("ESC: ") | color(Color::GrayLight),
                    text("⟨ " + cell_twin_.get_time_scale_label() + " ⟩") | bold | color(Color::YellowLight),
                    text(" │ ") | color(Color::GrayDark),
                    text("RELOJ: ") | color(Color::GrayLight),
                    text(cell_twin_.get_simulated_time_string()) | bold | color(Color::GreenLight),
                    filler(),
                    text("⟦ ÓPTICA: " + zoom_str + " | VEL: " + speed_str + " ⟧") | bold | color(Color::MagentaLight),
                    text(" "),
                    text(t_metrics.dox_active ? "⟦ OSKM: INDUCIDO ⟧" : "⟦ OSKM: BASAL ⟧") | bold | 
                         (t_metrics.dox_active ? color(Color::Yellow) : color(Color::GrayDark)),
                    text(" "),
                    text(t_metrics.crispr_applied ? "⟦ CRISPR: TP53 RESCATADO ⟧" : "⟦ CRISPR: INACTIVO ⟧") | bold | 
                         (t_metrics.crispr_applied ? color(Color::Cyan) : color(Color::GrayDark)),
                    text(" "),
                    text(t_metrics.uv_insult_active ? "⟦ UV: INSULTO GENOTÓXICO ⟧" : "⟦ UV: ESTABLE ⟧") | bold | 
                         (t_metrics.uv_insult_active ? color(Color::RedLight) : color(Color::GrayDark)),
                }),
                separator(),
                canvas(std::move(petri_canvas)) | center | flex,
                separator(),
                hbox({
                    text(" ◈ LEYENDA CYTO: ") | bold | color(Color::White),
                    text("━═━ Fibroblasto Huso ") | color(Color::Cyan),
                    text("◐ Intermedio MET ") | color(Color::Yellow),
                    text("⬡ Colonia iPSC ") | bold | color(Color::GreenLight),
                    text("▲ Neoplasia Pleomórfica ") | bold | color(Color::RedLight),
                    text("⸫ Cuerpos Apoptosis ") | color(Color::GrayDark),
                    text("⧖ Clivaje Mitótico ") | bold | color(Color::MagentaLight),
                })
            }) | border | flex;

            // Escala temporal visual indicando el multiplicador actual
            std::string time_badge = "1x";
            switch (cell_twin_.get_time_warp()) {
                case genome::twin::TimeWarpScale::Seconds:  time_badge = "1 SEC/S [1x]"; break;
                case genome::twin::TimeWarpScale::Minutes:  time_badge = "1 MIN/S [60x]"; break;
                case genome::twin::TimeWarpScale::Hours:    time_badge = "1 HR/S [3,600x]"; break;
                case genome::twin::TimeWarpScale::Days:     time_badge = "1 DÍA/S [86,400x]"; break;
                case genome::twin::TimeWarpScale::Weeks:    time_badge = "1 SEM/S [604,800x]"; break;
                case genome::twin::TimeWarpScale::Months:   time_badge = "1 MES/S [2.59Mx]"; break;
                case genome::twin::TimeWarpScale::Years:    time_badge = "1 AÑO/S [31.5Mx]"; break;
                case genome::twin::TimeWarpScale::Decades:  time_badge = "1 DÉC/S [315Mx]"; break;
            }

            auto twin_telemetry_box = vbox({
                text("◈ CENSO BIOLÓGICO Y POBLACIONAL") | bold | color(Color::Cyan),
                separator(),
                hbox({ text("Población Total: ") | color(Color::GrayLight), text(std::to_string(t_metrics.total_cells) + " células") | bold | color(Color::White) }),
                hbox({ text("Fibroblastos Somáticos: ") | color(Color::GrayLight), text(std::to_string(t_metrics.somatic_count)) | color(Color::Cyan) }),
                hbox({ text("Transición Intermedia (MET): ") | color(Color::GrayLight), text(std::to_string(t_metrics.intermediate_count)) | color(Color::Yellow) }),
                hbox({ text("Colonias Madre Pluripotentes (iPSC): ") | color(Color::GrayLight), text(std::to_string(t_metrics.ipsc_count)) | bold | color(Color::GreenLight) }),
                hbox({ text("Clones Malignos Transformados: ") | color(Color::GrayLight), text(std::to_string(t_metrics.transformed_count)) | bold | color(Color::RedLight) }),
                hbox({ text("Células en Apoptosis: ") | color(Color::GrayLight), text(std::to_string(t_metrics.apoptotic_count)) | color(Color::GrayDark) }),
                separator(),
                text("◈ CONTROLADOR TEMPORAL (TIME WARP)") | bold | color(Color::YellowLight),
                separator(),
                hbox({ text("Multiplicador: ") | color(Color::GrayLight), text(time_badge) | bold | color(Color::YellowLight) }),
                hbox({ text("Reloj Bio: ") | color(Color::GrayLight), text(cell_twin_.get_simulated_time_string()) | bold | color(Color::GreenLight) }),
                hbox({
                    text("Atajos: ") | color(Color::GrayLight),
                    text("[ [ ] ") | bold | color(Color::YellowLight), text("Frenar  ") | color(Color::White),
                    text("[ ] ] ") | bold | color(Color::YellowLight), text("Acelerar") | color(Color::White),
                }),
                separator(),
                text("◈ CINÉTICA REAL & MICROSCOPIO") | bold | color(Color::MagentaLight),
                separator(),
                hbox({ text("Velocidad Motor: ") | color(Color::GrayLight), text(std::to_string(sim_speed_multiplier_) + "x (" + std::to_string(sim_speed_multiplier_ * 20) + " Hz)") | bold | color(Color::MagentaLight) }),
                hbox({ text("Zoom Óptico: ") | color(Color::GrayLight), text(zoom_str) | bold | color(Color::Cyan), text("  Pan: ") | color(Color::GrayLight), text("(" + std::to_string(static_cast<int>(petri_pan_x_)) + "," + std::to_string(static_cast<int>(petri_pan_y_)) + ")") | color(Color::White) }),
                hbox({
                    text("Velocidad: ") | color(Color::GrayLight),
                    text("[ < ] ") | bold | color(Color::MagentaLight), text("Bajar  ") | color(Color::White),
                    text("[ > ] ") | bold | color(Color::MagentaLight), text("Subir") | color(Color::White),
                }),
                hbox({
                    text("Zoom Petri: ") | color(Color::GrayLight),
                    text("[Z] ") | bold | color(Color::Cyan), text("Acercar  ") | color(Color::White),
                    text("[z] ") | bold | color(Color::Cyan), text("Alejar  ") | color(Color::White),
                    text("[0] ") | bold | color(Color::Cyan), text("Reset") | color(Color::White),
                }),
                hbox({
                    text("Navegación: ") | color(Color::GrayLight),
                    text("[Flechas / W,A,S] ") | bold | color(Color::YellowLight), text("Pan Óptico") | color(Color::White),
                }),
                separator(),
                text("◈ DINÁMICA DE REPROGRAMACIÓN Y RIESGO") | bold | color(Color::Cyan),
                separator(),
                hbox({
                    text("Confluencia Cultivo: ") | color(Color::GrayLight),
                    text(std::to_string(static_cast<int>(t_metrics.confluence * 100.0f)) + "% ") | color(Color::White),
                    gauge(t_metrics.confluence) | size(WIDTH, EQUAL, 12) | color(Color::Cyan)
                }),
                hbox({
                    text("Índice Mitótico (Fase M): ") | color(Color::GrayLight),
                    text(std::to_string(static_cast<int>(t_metrics.mitotic_index)) + "% ") | color(Color::White),
                    gauge(t_metrics.mitotic_index / 25.0f) | size(WIDTH, EQUAL, 12) | color(Color::Yellow)
                }),
                hbox({
                    text("Rendimiento iPSC: ") | color(Color::GrayLight),
                    text(std::to_string(static_cast<int>(t_metrics.ipsc_efficiency)) + "% ") | bold | color(Color::Green),
                    gauge(t_metrics.ipsc_efficiency / 100.0f) | size(WIDTH, EQUAL, 12) | color(Color::Green)
                }),
                hbox({
                    text("R-Score Medio Dish: ") | color(Color::GrayLight),
                    text(std::to_string(t_metrics.average_r_score).substr(0, 4) + " ") | bold | 
                         (t_metrics.average_r_score > 0.5f ? color(Color::RedLight) : color(Color::Green)),
                    gauge(t_metrics.average_r_score) | size(WIDTH, EQUAL, 12) | 
                         (t_metrics.average_r_score > 0.5f ? color(Color::RedLight) : color(Color::Green))
                }),
                separator(),
                text("◈ TERAPIA DIRIGIDA IN-SILICO") | bold | color(Color::GreenLight),
                separator(),
                hbox({ text("Fármaco Activo: ") | color(Color::GrayLight), text(t_metrics.active_drug_name) | bold | color(Color::CyanLight) }),
                hbox({ text("Regresión Tumoral: ") | color(Color::GrayLight), text(std::to_string(static_cast<int>(t_metrics.tumor_regression_percent)) + "%") | bold | color(Color::GreenLight), text("  Dosis: ") | color(Color::GrayLight), text(std::to_string(t_metrics.drug_dose_uM).substr(0, 4) + " uM") | color(Color::White) }),
                hbox({ text("Clones Resistentes: ") | color(Color::GrayLight), text(std::to_string(t_metrics.drug_resistant_cells)) | bold | (t_metrics.drug_resistant_cells > 0 ? color(Color::RedLight) : color(Color::GreenLight)) }),
                hbox({
                    text("Acciones: ") | color(Color::GrayLight),
                    text("[T] ") | bold | color(Color::CyanLight), text("TKI  ") | color(Color::White),
                    text("[P] ") | bold | color(Color::MagentaLight), text("PARP  ") | color(Color::White),
                    text("[K] ") | bold | color(Color::YellowLight), text("p53  ") | color(Color::White),
                    text("[W] ") | bold | color(Color::YellowLight), text("Lavar") | color(Color::White),
                }),
                separator(),
                text("◈ CONTROLES BIOLÓGICOS Y EXPORTACIÓN") | bold | color(Color::Cyan),
                separator(),
                vbox({
                    text("[X] EXPORTAR DATASETS CSV (Censo/Células/VCF)") | bold | color(Color::GreenLight),
                    text("[O] EXPORTAR REPORTE COMITÉ MOLECULAR (MTB)") | bold | color(Color::MagentaLight),
                    text("[D] Inyectar/Suspender Doxiciclina (OSKM)") | bold | color(Color::Yellow),
                    text("[C] Terapia CRISPR-Cas9 in-situ (Rescate TP53)") | bold | color(Color::Cyan),
                    text("[U] Estrés Genotóxico (Irradiación UV)") | bold | color(Color::Magenta),
                    text("[+] / [-] Inocular / Remover Células (±25)") | color(Color::Cyan),
                    text("[R] Reset Cultivo a Estado Basal") | color(Color::White),
                    text("[Espacio] Pausar / Reanudar Dinámica") | color(Color::GrayLight),
                })
            }) | border | size(WIDTH, EQUAL, 48);

            active_view = hbox({
                petri_box,
                twin_telemetry_box
            }) | flex;
        }

        // -------------------------------------------------------------
        // TAB 8: ONCO THERAPY (PHARMACOGENOMICS & MOLECULAR TUMOR BOARD)
        // -------------------------------------------------------------
        else if (current_tab_ == AppTab::OncoTherapy) {
            refresh_onco_analysis();

            // 1. Cabecera Cyber-HUD
            auto header_box = hbox({
                vbox({
                    text("⟦ COMITÉ MOLECULAR DE TUMORES & FARMACOGENÓMICA DIRIGIDA ⟧") | bold | color(Color::CyanLight),
                    text("⟨ Emparejamiento de Biomarcadores Genómicos (OncoKB / FDA / NCCN) con Eficacia In-Silico ⟩") | color(Color::GrayLight),
                }) | flex,
                separator(),
                vbox({
                    hbox({ text("Biomarcadores Detectados: ") | color(Color::GrayLight), text(std::to_string(onco_analysis_.detected_biomarkers.size())) | bold | color(Color::Yellow) }),
                    hbox({ text("Fármacos Aprobados FDA:   ") | color(Color::GrayLight), text(std::to_string(onco_analysis_.level1_count)) | bold | color(Color::GreenLight) }),
                }),
                separator(),
                vbox({
                    hbox({ text("Carga Mutacional (TMB):   ") | color(Color::GrayLight), text(std::to_string(static_cast<int>(onco_analysis_.tumor_mutational_burden)) + " mut/Mb") | bold | (onco_analysis_.high_tmb ? color(Color::RedLight) : color(Color::Cyan)) }),
                    hbox({ text("Puntaje Accionabilidad:   ") | color(Color::GrayLight), text(std::to_string(static_cast<int>(onco_analysis_.total_actionability_score * 100.0f)) + "%") | bold | color(Color::CyanLight) }),
                })
            }) | border;

            // 2. Columna Izquierda: Biomarcadores Detectados en el Genoma
            Elements biomarker_elements;
            biomarker_elements.push_back(text("◈ BIOMARCADORES Y DIANAS ONCOLÓGICAS") | bold | color(Color::Yellow));
            biomarker_elements.push_back(separator());

            if (onco_analysis_.detected_biomarkers.empty()) {
                biomarker_elements.push_back(text("[+] 0 mutaciones conductoras detectadas en los loci diana.") | color(Color::Green));
                biomarker_elements.push_back(text("    Perfil genómico compatible con control sano.") | color(Color::GrayLight));
            } else {
                for (size_t i = 0; i < onco_analysis_.detected_biomarkers.size(); ++i) {
                    const auto& b = onco_analysis_.detected_biomarkers[i];
                    biomarker_elements.push_back(
                        vbox({
                            hbox({
                                text(b.is_resistance_marker ? "[RESISTENCIA] " : "[DRIVER] ") | (b.is_resistance_marker ? color(Color::RedLight) : color(Color::YellowLight)),
                                text(b.gene_symbol + " ") | bold | color(Color::CyanLight),
                                text(b.variant_desc) | bold | color(Color::White)
                            }),
                            text("  Locus: Chr " + std::to_string(b.genomic_pos) + " | " + b.consequence) | color(Color::GrayLight),
                        })
                    );
                    biomarker_elements.push_back(separator());
                }
            }
            auto biomarkers_panel = vbox(std::move(biomarker_elements)) | border | size(WIDTH, GREATER_THAN, 42);

            // 3. Columna Derecha: Prescripción de Fármacos Dirigidos
            Elements drug_elements;
            drug_elements.push_back(text("◈ TERAPIAS DIRIGIDAS EMPAREJADAS (NIVEL DE EVIDENCIA CLÍNICA)") | bold | color(Color::CyanLight));
            drug_elements.push_back(separator());

            if (onco_analysis_.all_recommended_drugs.empty()) {
                drug_elements.push_back(text("No se requieren terapias dirigidas citotóxicas para este perfil.") | color(Color::Green));
            } else {
                for (size_t i = 0; i < onco_analysis_.all_recommended_drugs.size(); ++i) {
                    const auto& d = onco_analysis_.all_recommended_drugs[i];
                    Color badge_col = (d.evidence_level == oncology::EvidenceLevel::Level1_FDA_Approved) ? Color::GreenLight :
                                      (d.evidence_level == oncology::EvidenceLevel::Level2_StandardOfCare) ? Color::CyanLight : Color::Yellow;
                    
                    drug_elements.push_back(
                        hbox({
                            text(oncology::ActionabilityEngine::evidence_level_badge(d.evidence_level) + " ") | bold | color(badge_col),
                            text(d.drug_name) | bold | color(Color::White),
                            text(" (" + d.drug_class + ")") | color(Color::GrayLight),
                            filler(),
                            text("Diana: " + d.target_gene + " | Eficacia: " + std::to_string(static_cast<int>(d.tumor_regression_efficacy * 100)) + "%") | color(Color::Cyan)
                        })
                    );
                    drug_elements.push_back(
                        text("  Mecanismo: " + d.mechanism_of_action + " [" + d.clinical_trial_or_approval + "]") | color(Color::GrayLight)
                    );
                    drug_elements.push_back(
                        text("  Pauta: " + d.recommended_dose + (d.synthetic_lethality_score > 0.5f ? " | LETALIDAD SINTÉTICA ACTIVA" : "")) | color(Color::MagentaLight)
                    );
                    drug_elements.push_back(separator());
                }
            }

            // Advertencias de resistencia
            if (!onco_analysis_.resistance_warnings.empty()) {
                drug_elements.push_back(text("◈ ALERTAS FARMACOGENÓMICAS DE RESISTENCIA") | bold | color(Color::RedLight));
                for (const auto& w : onco_analysis_.resistance_warnings) {
                    drug_elements.push_back(text("[!] " + w) | color(Color::RedLight));
                }
                drug_elements.push_back(separator());
            }

            // Panel de Control Farmacodinámico en el Gemelo Digital
            auto pharmacodynamics_controls = hbox({
                text("⟦ FARMACODINÁMICA IN-SILICO (DIGITAL TWIN) ⟧ ") | bold | color(Color::Yellow),
                text("Fármaco Activo: ") | color(Color::GrayLight),
                text(cell_twin_.get_active_drug_name()) | bold | color(Color::CyanLight),
                text(" | Regresión: ") | color(Color::GrayLight),
                text(std::to_string(static_cast<int>(cell_twin_.get_tumor_regression_percent())) + "%") | bold | color(Color::GreenLight),
                text(" | Clones Resistentes: ") | color(Color::GrayLight),
                text(std::to_string(cell_twin_.get_metrics().drug_resistant_cells)) | bold | (cell_twin_.get_metrics().drug_resistant_cells > 0 ? color(Color::RedLight) : color(Color::GreenLight)),
                filler(),
                text("[T] TKI  [P] PARP  [K] p53-Reactivator  [W] Lavar  [O] Reporte MTB") | bold | color(Color::Black) | bgcolor(Color::CyanLight)
            });

            auto drugs_panel = vbox(std::move(drug_elements)) | border | flex;

            active_view = vbox({
                header_box,
                hbox({ biomarkers_panel, drugs_panel }) | flex,
                pharmacodynamics_controls | border
            });
        }

        // 4. FOOTER
        auto footer = hbox({
            text(" ⟦ ATAJOS ⟧ ") | bold | color(Color::Cyan),
            text("[1..8] ") | bold | color(Color::White), text("Pestañas  ") | color(Color::GrayLight),
            text("[Tab] ") | bold | color(Color::White), text("Ciclar  ") | color(Color::GrayLight),
            text("[M] ") | bold | color(Color::GreenLight), text("Muestra  ") | color(Color::GrayLight),
            text("[L] ") | bold | color(Color::Yellow), text("Locus  ") | color(Color::GrayLight),
            text("[O] ") | bold | color(Color::MagentaLight), text("Reporte MTB  ") | color(Color::GrayLight),
            text("[T] ") | bold | color(Color::CyanLight), text("TKI  ") | color(Color::GrayLight),
            text("[P] ") | bold | color(Color::YellowLight), text("PARP  ") | color(Color::GrayLight),
            text("[ < / > ] ") | bold | color(Color::MagentaLight), text("Vel Real  ") | color(Color::GrayLight),
            text("[Z / z / 0] ") | bold | color(Color::Cyan), text("Zoom Petri  ") | color(Color::GrayLight),
            text("[ [ / ] ] ") | bold | color(Color::YellowLight), text("Time Warp  ") | color(Color::GrayLight),
            text("[X] ") | bold | color(Color::GreenLight), text("CSV  ") | color(Color::GrayLight),
            text("[Espacio] ") | bold | color(Color::White), text("Pausa  ") | color(Color::GrayLight),
            text("[Q] ") | bold | color(Color::Red), text("Salir") | color(Color::GrayLight),
        }) | border;

        return vbox({
            header_element,
            tab_bar,
            active_view,
            footer
        });
    });

    // Manejo de eventos de teclado
    auto event_handler = CatchEvent(main_renderer, [&](Event event) -> bool {
        if (event == Event::Character('q') || event == Event::Character('Q') || event == Event::Escape) {
            screen.ExitLoopClosure()();
            return true;
        }

        if (event == Event::Character('1')) { current_tab_ = AppTab::Dashboard; return true; }
        if (event == Event::Character('2')) { current_tab_ = AppTab::VcfExplorer; return true; }
        if (event == Event::Character('3')) { current_tab_ = AppTab::YamanakaDynamics; return true; }
        if (event == Event::Character('4')) { current_tab_ = AppTab::MoeRouter; return true; }
        if (event == Event::Character('5')) { current_tab_ = AppTab::CrisprStudio; return true; }
        if (event == Event::Character('6')) { current_tab_ = AppTab::SystemLogs; return true; }
        if (event == Event::Character('7')) { current_tab_ = AppTab::DigitalTwin; return true; }
        if (event == Event::Character('8')) { current_tab_ = AppTab::OncoTherapy; return true; }

        if (event == Event::Tab) {
            int next_tab = (static_cast<int>(current_tab_) + 1) % 8;
            current_tab_ = static_cast<AppTab>(next_tab);
            return true;
        }

        // Exportar Reporte de Comité Molecular de Tumores (MTB)
        if (event == Event::Character('o') || event == Event::Character('O')) {
            export_molecular_tumor_board_report();
            return true;
        }

        // Inocular terapia dirigida TKI
        if (event == Event::Character('t') || event == Event::Character('T')) {
            administer_targeted_drug(genome::twin::DrugTherapy::TargetedTKI);
            return true;
        }

        // Inocular terapia por letalidad sintética (Inhibidor PARP)
        if (event == Event::Character('p') || event == Event::Character('P')) {
            administer_targeted_drug(genome::twin::DrugTherapy::ParpInhibitor);
            return true;
        }

        // Inocular reactivador conformacional p53 (APR-246)
        if (event == Event::Character('k') || event == Event::Character('K')) {
            administer_targeted_drug(genome::twin::DrugTherapy::P53Reactivator);
            return true;
        }

        // Lavado / Washout de fármacos en Tab 8
        if (current_tab_ == AppTab::OncoTherapy && (event == Event::Character('w') || event == Event::Character('W'))) {
            washout_drugs();
            return true;
        }

        if (event == Event::Character('l') || event == Event::Character('L')) {
            cycle_locus();
            return true;
        }

        if (event == Event::Character('m') || event == Event::Character('M')) {
            toggle_sample_mode();
            return true;
        }

        // Control de velocidad en tiempo real (1x..32x)
        if (event == Event::Character('<') || event == Event::Character(',')) {
            decrease_sim_speed();
            return true;
        }

        if (event == Event::Character('>') || event == Event::Character('.')) {
            increase_sim_speed();
            return true;
        }

        // Controles de zoom óptico del microscopio
        if (event == Event::Character('Z')) {
            zoom_in_petri();
            return true;
        }

        if (event == Event::Character('z')) {
            zoom_out_petri();
            return true;
        }

        if (event == Event::Character('0')) {
            reset_petri_view();
            return true;
        }

        // Navegación Pan en la placa de Petri (Tab 7)
        if (current_tab_ == AppTab::DigitalTwin) {
            float pan_step = 4.0f / petri_zoom_;
            if (event == Event::Character('w') || event == Event::Character('W') || event == Event::ArrowUp) {
                pan_petri(0.0f, -pan_step);
                return true;
            }
            if (event == Event::Character('s') || event == Event::Character('S') || event == Event::ArrowDown) {
                pan_petri(0.0f, pan_step);
                return true;
            }
            if (event == Event::Character('a') || event == Event::Character('A') || event == Event::ArrowLeft) {
                pan_petri(-pan_step, 0.0f);
                return true;
            }
            if (event == Event::ArrowRight) {
                pan_petri(pan_step, 0.0f);
                return true;
            }
        }

        // Acciones del Gemelo Digital
        if (event == Event::Character('d') || event == Event::Character('D')) {
            cell_twin_.toggle_doxycycline();
            add_log("TWIN", cell_twin_.is_dox_active() ? 
                    "Doxycycline OSKM pulse INJECTED into culture medium." : 
                    "Doxycycline OSKM REMOVED from culture medium.");
            return true;
        }

        if (event == Event::Character('c') || event == Event::Character('C')) {
            cell_twin_.apply_crispr_rescue();
            add_log("CRISPR", "Targeted Cas9/HDR transfection rescued TP53 mutations across all living cells in dish.");
            return true;
        }

        if (event == Event::Character('u') || event == Event::Character('U')) {
            cell_twin_.trigger_uv_damage();
            add_log("ALERT", "UV genotoxic irradiation insult triggered. Testing DNA repair checkpoints.");
            return true;
        }

        if (event == Event::Character('+') || event == Event::Character('=')) {
            cell_twin_.adjust_population(+25);
            add_log("TWIN", "Inoculated +25 cells into virtual Petri dish.");
            return true;
        }

        if (event == Event::Character('-') || event == Event::Character('_')) {
            cell_twin_.adjust_population(-25);
            add_log("TWIN", "Aspirated -25 cells from virtual Petri dish.");
            return true;
        }

        // Controles de escala temporal (Time Warp)
        if (event == Event::Character('[') || event == Event::Character('{')) {
            cell_twin_.slow_down_time();
            add_log("TIME", "Biological simulation time slowed to: " + cell_twin_.get_time_scale_label());
            return true;
        }

        if (event == Event::Character(']') || event == Event::Character('}')) {
            cell_twin_.speed_up_time();
            add_log("TIME", "Biological simulation time accelerated to: " + cell_twin_.get_time_scale_label());
            return true;
        }

        // Exportación de datasets científicos a CSV
        if (event == Event::Character('x') || event == Event::Character('X')) {
            export_datasets();
            return true;
        }

        if (event == Event::ArrowUp || event == Event::Character('k')) {
            if (current_tab_ == AppTab::VcfExplorer) {
                if (selected_variant_index_ > 0) selected_variant_index_--;
            } else if (current_tab_ == AppTab::CrisprStudio) {
                if (selected_guide_index_ > 0) selected_guide_index_--;
            }
            return true;
        }

        if (event == Event::ArrowDown || event == Event::Character('j')) {
            if (current_tab_ == AppTab::VcfExplorer) {
                if (selected_variant_index_ + 1 < static_cast<int>(vcf_parser_.total_variants())) {
                    selected_variant_index_++;
                }
            } else if (current_tab_ == AppTab::CrisprStudio) {
                if (selected_guide_index_ + 1 < static_cast<int>(crispr_candidates_.size())) {
                    selected_guide_index_++;
                }
            }
            return true;
        }

        if (event == Event::Character(' ')) {
            yamanaka_running_ = !yamanaka_running_;
            add_log("INFO", yamanaka_running_ ? "Simulation engines resumed." : "Simulation engines paused.");
            return true;
        }

        if (event == Event::Character('e') || event == Event::Character('E')) {
            if (!crispr_applied_) {
                vcf_parser_.apply_virtual_correction(selected_variant_index_);
                crispr_applied_ = true;

                // Actualizar contexto en el simulador biológico
                float new_burden = static_cast<float>(vcf_parser_.pathogenic_count()) / 
                                   std::max(1.0f, static_cast<float>(vcf_parser_.total_variants()));
                yamanaka_sim_.set_mutational_context(false, false, new_burden);

                // Re-simular edición HDR en la secuencia
                target_sequence_window_ = crispr_engine_.simulate_hdr_edit(target_sequence_window_, 20, "CGC");
                crispr_candidates_ = crispr_engine_.design_guides(target_sequence_window_, 20);

                add_log("EDIT", "CRISPR Cas9 endonuclease cut simulated at PAM " + 
                                (crispr_candidates_.empty() ? "NGG" : crispr_candidates_[0].pam) + " locus.");
                add_log("EDIT", "Homology-Directed Repair (HDR) template integrated wild-type codon.");
                add_log("INFO", "Yamanaka ODE Solver: p53 checkpoint restored. R-Score dropped to 0.12 (SAFE).");
            } else {
                add_log("WARN", "CRISPR correction already applied to current target locus.");
            }
            return true;
        }

        if (event == Event::Character('r') || event == Event::Character('R')) {
            if (current_tab_ == AppTab::DigitalTwin) {
                cell_twin_.reset();
                add_log("TWIN", "Virtual Petri dish reset to basal cell state.");
            } else {
                vcf_parser_.load_file(vcf_path_);
                crispr_applied_ = false;

                bool tp53_mut = false;
                for (const auto& v : vcf_parser_.get_variants()) {
                    if (std::string_view(v.gene_symbol) == "TP53" && 
                        v.effect == genome::vcf::VariantEffect::YamanakaCollision) tp53_mut = true;
                }
                float burden = static_cast<float>(vcf_parser_.pathogenic_count()) / 
                               std::max(1.0f, static_cast<float>(vcf_parser_.total_variants()));
                yamanaka_sim_.set_mutational_context(tp53_mut, false, burden);
                yamanaka_sim_.reset();

                target_sequence_window_ = "CTGCTCCCCCCCTGCCCTCAACAAGATGTTTTGCCAGTTGGCTCTGACTGTACCACCATCCACTACAACT";
                crispr_candidates_ = crispr_engine_.design_guides(target_sequence_window_, 20);

                add_log("WARN", "Genome state reset to original patient VCF. TP53 R175H re-activated!");
                add_log("CRIT", "Oncogenic instability alert! R-Score surged back to high-risk state.");
            }
            return true;
        }

        return false;
    });

    Loop loop(&screen, event_handler);

    while (!loop.HasQuitted()) {
        for (size_t s = 0; s < sim_speed_multiplier_; ++s) {
            tick_simulation();
        }
        loop.RunOnce();
        std::this_thread::sleep_for(std::chrono::milliseconds(50));
    }

    return 0;
}

} // namespace genome::ui
