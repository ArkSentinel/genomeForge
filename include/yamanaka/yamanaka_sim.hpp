#pragma once

#include <vector>
#include <cstdint>
#include <string>

namespace genome::yamanaka {

struct YamanakaFactors {
    float oct4{0.85f};
    float sox2{0.80f};
    float klf4{0.90f};
    float cmyc{0.95f};
};

struct CellularState {
    float pluripotency{0.10f};       // Nivel de pluripotencia celular [0.0 - 1.0]
    float senescence{0.30f};          // Barrera de senescencia fisiológica
    float instability{0.20f};         // Inestabilidad genómica acumulada
    float p53_activity{1.0f};         // Actividad de p53 (1.0 = salvaje, 0.0 = mutado/inactivo)
    float p21_activity{1.0f};         // Freno ciclo celular p21/CDKN1A
    float r_score{0.15f};             // Oncogenic Risk Score [0.0 - 1.0]
};

class YamanakaSimulator {
public:
    YamanakaSimulator();
    ~YamanakaSimulator();

    // Configura el estado mutacional del paciente (procedente del VCF)
    void set_mutational_context(bool tp53_mutated, bool cdkn2a_mutated, float pathogenic_burden);

    // Ejecuta un paso temporal dt mediante Runge-Kutta de 4º Orden (RK4)
    void step(float dt);

    // Reinicia la simulación a condiciones basales
    void reset();

    // Control de factores exógenos
    void set_factors(const YamanakaFactors& factors) { factors_ = factors; }
    const YamanakaFactors& get_factors() const { return factors_; }

    const CellularState& get_state() const { return state_; }
    float get_time() const { return current_time_; }

    // Historial del R-Score para representación gráfica
    const std::vector<float>& get_r_score_history() const { return r_history_; }

private:
    // Derivadas del sistema dinámico
    struct StateDerivatives {
        float dP_dt;
        float dS_dt;
        float dI_dt;
        float dp53_dt;
        float dp21_dt;
    };

    StateDerivatives compute_derivatives(float t, const CellularState& s);
    void update_r_score();

    CellularState state_{};
    YamanakaFactors factors_{};
    float current_time_{0.0f};

    // Parámetros mutacionales del paciente
    bool tp53_mutated_{false};
    bool cdkn2a_mutated_{false};
    float pathogenic_burden_{0.0f};

    // Historial circular de R-Score
    std::vector<float> r_history_;
};

} // namespace genome::yamanaka
