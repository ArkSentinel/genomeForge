#pragma once

#include <vector>
#include <cstdint>
#include <random>
#include <string>

namespace genome::twin {

enum class CellCyclePhase : uint8_t {
    G0_Quiescent = 0,
    G1_Growth = 1,
    S_Synthesis = 2,
    G2_PreMitotic = 3,
    M_Mitosis = 4,
    Apoptotic = 5,
    Senescent = 6
};

enum class CellPhenotype : uint8_t {
    SomaticFibroblast = 0,   // Normal differentiated cell (Blue)
    IntermediateReprog = 1,  // Mesenchymal-to-Epithelial Transition MET (Yellow)
    PluripotentIPSC = 2,     // Reprogrammed stem cell colony (Green)
    TransformedMalignant = 3,// Oncogenic escape clone (Red)
    ApoptoticCorpse = 4      // Cell undergoing apoptosis (Gray)
};

enum class TimeWarpScale : uint8_t {
    Seconds = 0,    // 1x (Realtime)
    Minutes = 1,    // 60x
    Hours = 2,      // 3,600x
    Days = 3,       // 86,400x
    Weeks = 4,      // 604,800x
    Months = 5,     // 2,592,000x
    Years = 6,      // 31,536,000x
    Decades = 7     // 315,360,000x
};

enum class DrugTherapy : uint8_t {
    None = 0,
    TargetedTKI = 1,        // e.g. Osimertinib (EGFR)
    ParpInhibitor = 2,      // e.g. Olaparib (BRCA1 synthetic lethality)
    P53Reactivator = 3,     // e.g. APR-246 (TP53 conformation refolder)
    KrasInhibitor = 4,      // e.g. Sotorasib (KRAS G12C)
    CytotoxicChemo = 5      // Broad-spectrum anti-mitotic chemotherapy
};

// Alineación estricta a 64 bytes para línea de caché ARM64
struct alignas(64) CellAgent {
    float x{0.0f};           // Coordenada X en la placa de Petri
    float y{0.0f};           // Coordenada Y en la placa de Petri
    float vx{0.0f};          // Velocidad mecánica Vx
    float vy{0.0f};          // Velocidad mecánica Vy
    float radius{2.2f};      // Radio morfológico celular
    float mass{1.0f};        // Masa inercial
    CellCyclePhase phase{CellCyclePhase::G1_Growth};
    CellPhenotype phenotype{CellPhenotype::SomaticFibroblast};
    uint8_t cleavage_state{0}; // 0 = unilocular, 1 = constricción anafásica, 2 = telofase
    uint8_t pad0{0};         // Alineamiento a múltiplo de 4 bytes
    float cycle_clock{0.0f}; // Progreso en la fase actual (0.0 a 1.0)
    float P_score{0.05f};    // Nivel de Pluripotencia (Oct4/Sox2/Nanog)
    float S_score{0.95f};    // Nivel de Identidad Somática
    float p53_level{0.80f};  // Proteína p53 activa
    float r_score{0.10f};    // Índice de Riesgo Oncogénico
    uint8_t has_tp53_mutation{0}; // 1 si porta variante patogénica TP53
    uint8_t is_transformed{0};    // 1 si ha escapado a control proliferativo
    uint16_t generation{0};       // Generación mitótica
    uint32_t cell_id{0};          // Identificador único
    uint32_t parent_id{0};        // ID progenitor
    uint8_t drug_resistant{0};    // 1 si adquirió resistencia secundaria al fármaco
    uint8_t pad1[3]{0};           // Relleno final para totalizar exactamente 64 bytes
};
static_assert(sizeof(CellAgent) == 64, "CellAgent must be exactly 64 bytes");

struct PopulationMetrics {
    size_t total_cells{0};
    size_t somatic_count{0};
    size_t intermediate_count{0};
    size_t ipsc_count{0};
    size_t transformed_count{0};
    size_t apoptotic_count{0};
    float mitotic_index{0.0f};       // Porcentaje en fase M
    float average_r_score{0.0f};     // Riesgo oncogénico poblacional
    float confluence{0.0f};          // Confluencia del cultivo (0.0 a 1.0)
    float ipsc_efficiency{0.0f};     // Rendimiento de reprogramación (%)
    bool dox_active{false};          // Pulso de doxiciclina / Yamanaka activo
    bool crispr_applied{false};      // Edición génica administrada
    bool uv_insult_active{false};    // Estrés por radiación genotóxica
    DrugTherapy active_drug{DrugTherapy::None};
    std::string active_drug_name{"None (Vehicle Control)"};
    float drug_dose_uM{0.0f};
    float tumor_regression_percent{0.0f};
    size_t drug_resistant_cells{0};
    TimeWarpScale time_warp{TimeWarpScale::Seconds};
    double time_multiplier{1.0};
    double simulated_seconds{0.0};
    std::string simulated_time_str{"T+ 000y:00m:00d:00h:00m:00s"};
};

class CellTwinSimulator {
public:
    CellTwinSimulator();
    ~CellTwinSimulator();

    void initialize_dish(size_t initial_count = 160, bool patient_has_tp53_mut = true);
    void reset();

    // Avanza un paso temporal de la biofísica y del ciclo celular
    void step(float base_dt = 0.05f);

    // Controles de escala temporal biológica (Time Warp)
    void speed_up_time();
    void slow_down_time();
    void set_time_warp(TimeWarpScale scale);
    TimeWarpScale get_time_warp() const { return time_warp_; }
    double get_time_multiplier() const;
    std::string get_time_scale_label() const;
    const std::string& get_simulated_time_string() const { return metrics_.simulated_time_str; }

    // Controles biológicos interactivos
    void toggle_doxycycline();
    void set_doxycycline(bool active) { dox_active_ = active; }
    void apply_crispr_rescue();
    void trigger_uv_damage();
    void adjust_population(int delta);

    // Controles farmacológicos in-silico
    void administer_drug(DrugTherapy drug, float dose_uM = 1.0f);
    void washout_drugs();
    DrugTherapy get_active_drug() const { return active_drug_; }
    const std::string& get_active_drug_name() const { return active_drug_name_; }
    float get_drug_concentration() const { return drug_dose_uM_; }
    float get_tumor_regression_percent() const { return metrics_.tumor_regression_percent; }

    const std::vector<CellAgent>& get_cells() const { return cells_; }
    const PopulationMetrics& get_metrics() const { return metrics_; }
    float get_dish_radius() const { return dish_radius_; }
    bool is_dox_active() const { return dox_active_; }
    bool is_crispr_applied() const { return crispr_applied_; }

private:
    void update_cell_cycle(CellAgent& cell, float effective_dt, std::vector<CellAgent>& new_daughters);
    void update_physics_cpu(float dt);
    void compute_population_metrics();
    void update_time_string();

    float dish_radius_{45.0f};       // Radio del cultivo en la placa virtual
    bool dox_active_{false};         // Presencia de inductor Yamanaka
    bool crispr_applied_{false};     // Si se corrigió TP53
    bool uv_insult_active_{false};   // Daño inducido en ADN
    float uv_timer_{0.0f};

    // Farmacodinámica
    DrugTherapy active_drug_{DrugTherapy::None};
    std::string active_drug_name_{"None (Vehicle Control)"};
    float drug_dose_uM_{0.0f};
    float initial_transformed_baseline_{0.0f};

    TimeWarpScale time_warp_{TimeWarpScale::Seconds};
    double simulated_seconds_{0.0};

    uint32_t next_cell_id_{1};
    std::vector<CellAgent> cells_;
    PopulationMetrics metrics_;
    std::mt19937 rng_;
};

} // namespace genome::twin
