#include "yamanaka/yamanaka_sim.hpp"

#include <cmath>
#include <algorithm>

namespace genome::yamanaka {

YamanakaSimulator::YamanakaSimulator() {
    reset();
}

YamanakaSimulator::~YamanakaSimulator() = default;

void YamanakaSimulator::reset() {
    current_time_ = 0.0f;
    state_.pluripotency = 0.08f;
    state_.senescence = 0.35f;
    state_.instability = tp53_mutated_ ? 0.75f : 0.08f;
    state_.p53_activity = tp53_mutated_ ? 0.08f : 0.95f;
    state_.p21_activity = cdkn2a_mutated_ ? 0.15f : 0.88f;

    r_history_.clear();
    update_r_score();
    r_history_.resize(40, state_.r_score);
}

void YamanakaSimulator::set_mutational_context(bool tp53_mutated, bool cdkn2a_mutated, float pathogenic_burden) {
    tp53_mutated_ = tp53_mutated;
    cdkn2a_mutated_ = cdkn2a_mutated;
    pathogenic_burden_ = std::clamp(pathogenic_burden, 0.0f, 1.0f);

    if (tp53_mutated_) {
        state_.p53_activity = std::min(state_.p53_activity, 0.10f);
    } else {
        state_.p53_activity = std::max(state_.p53_activity, 0.90f);
    }

    if (cdkn2a_mutated_) {
        state_.p21_activity = std::min(state_.p21_activity, 0.20f);
    } else {
        state_.p21_activity = std::max(state_.p21_activity, 0.85f);
    }

    update_r_score();
}

YamanakaSimulator::StateDerivatives YamanakaSimulator::compute_derivatives(float t, const CellularState& s) {
    // Oscilación o pulsos circadianos en factores exógenos
    float pulse = 1.0f + 0.03f * std::sin(t * 1.5f);
    float oct4 = factors_.oct4 * pulse;
    float sox2 = factors_.sox2 * pulse;
    float klf4 = factors_.klf4 * pulse;
    float cmyc = factors_.cmyc * pulse;

    // 1. Pluripotencia: activación cooperativa Oct4-Sox2 + Klf4 - degradación
    const float k_p1 = 0.45f;
    const float k_p2 = 0.30f;
    const float gamma_p = 0.35f;
    float dP_dt = (k_p1 * oct4 * sox2 + k_p2 * klf4) * (1.0f - s.pluripotency) - gamma_p * s.pluripotency;

    // 2. Senescencia: barrera protectora activada por c-Myc mediada por p53
    const float k_s1 = 0.38f;
    const float k_s2 = 0.20f;
    const float gamma_s = 0.30f;
    float dS_dt = (k_s1 * cmyc * s.p53_activity + k_s2 * klf4) * (1.0f - s.senescence) - gamma_s * s.senescence;

    // 3. Actividad de p53: diana de estrés de replicación inducido por oncogenes (c-Myc)
    // Si TP53 está mutado en el genoma, la tasa decae o permanece inerte
    float dp53_dt = 0.0f;
    if (tp53_mutated_) {
        dp53_dt = -0.15f * s.p53_activity;
    } else {
        // En células wild-type, c-Myc induce estabilización de p53
        dp53_dt = 0.25f * cmyc * (1.0f - s.p53_activity) - 0.10f * s.p53_activity;
    }

    // 4. Actividad de p21/CDKN1A (control del ciclo celular)
    float dp21_dt = 0.0f;
    if (cdkn2a_mutated_) {
        dp21_dt = -0.12f * s.p21_activity;
    } else {
        dp21_dt = 0.30f * s.p53_activity * (1.0f - s.p21_activity) - 0.15f * s.p21_activity;
    }

    // 5. Inestabilidad Genómica: amplificación por c-Myc en ausencia de freno p53
    const float k_i1 = 0.85f;
    const float k_i2 = 0.40f;
    const float gamma_i = 0.25f;
    float uncheck_factor = (1.0f - s.p53_activity);
    float dI_dt = (k_i1 * (cmyc * cmyc) * uncheck_factor + k_i2 * pathogenic_burden_) - gamma_i * s.instability;

    return {dP_dt, dS_dt, dI_dt, dp53_dt, dp21_dt};
}

void YamanakaSimulator::update_r_score() {
    float oncog_drive = factors_.cmyc / (state_.p53_activity + 0.08f) * 0.15f;
    float raw_r = 0.45f * state_.instability + 0.35f * oncog_drive + 0.20f * pathogenic_burden_;
    state_.r_score = std::clamp(raw_r, 0.05f, 0.98f);
}

void YamanakaSimulator::step(float dt) {
    // Implementación de Runge-Kutta de 4º Orden (RK4)
    float t = current_time_;

    // k1
    auto k1 = compute_derivatives(t, state_);

    // k2
    CellularState s_k2 = state_;
    s_k2.pluripotency = std::clamp(state_.pluripotency + 0.5f * dt * k1.dP_dt, 0.0f, 1.0f);
    s_k2.senescence   = std::clamp(state_.senescence   + 0.5f * dt * k1.dS_dt, 0.0f, 1.0f);
    s_k2.instability  = std::clamp(state_.instability  + 0.5f * dt * k1.dI_dt, 0.0f, 1.0f);
    s_k2.p53_activity = std::clamp(state_.p53_activity + 0.5f * dt * k1.dp53_dt, 0.0f, 1.0f);
    s_k2.p21_activity = std::clamp(state_.p21_activity + 0.5f * dt * k1.dp21_dt, 0.0f, 1.0f);
    auto k2 = compute_derivatives(t + 0.5f * dt, s_k2);

    // k3
    CellularState s_k3 = state_;
    s_k3.pluripotency = std::clamp(state_.pluripotency + 0.5f * dt * k2.dP_dt, 0.0f, 1.0f);
    s_k3.senescence   = std::clamp(state_.senescence   + 0.5f * dt * k2.dS_dt, 0.0f, 1.0f);
    s_k3.instability  = std::clamp(state_.instability  + 0.5f * dt * k2.dI_dt, 0.0f, 1.0f);
    s_k3.p53_activity = std::clamp(state_.p53_activity + 0.5f * dt * k2.dp53_dt, 0.0f, 1.0f);
    s_k3.p21_activity = std::clamp(state_.p21_activity + 0.5f * dt * k2.dp21_dt, 0.0f, 1.0f);
    auto k3 = compute_derivatives(t + 0.5f * dt, s_k3);

    // k4
    CellularState s_k4 = state_;
    s_k4.pluripotency = std::clamp(state_.pluripotency + dt * k3.dP_dt, 0.0f, 1.0f);
    s_k4.senescence   = std::clamp(state_.senescence   + dt * k3.dS_dt, 0.0f, 1.0f);
    s_k4.instability  = std::clamp(state_.instability  + dt * k3.dI_dt, 0.0f, 1.0f);
    s_k4.p53_activity = std::clamp(state_.p53_activity + dt * k3.dp53_dt, 0.0f, 1.0f);
    s_k4.p21_activity = std::clamp(state_.p21_activity + dt * k3.dp21_dt, 0.0f, 1.0f);
    auto k4 = compute_derivatives(t + dt, s_k4);

    // Actualización integrada RK4
    state_.pluripotency = std::clamp(state_.pluripotency + (dt / 6.0f) * (k1.dP_dt + 2.0f * k2.dP_dt + 2.0f * k3.dP_dt + k4.dP_dt), 0.0f, 1.0f);
    state_.senescence   = std::clamp(state_.senescence   + (dt / 6.0f) * (k1.dS_dt + 2.0f * k2.dS_dt + 2.0f * k3.dS_dt + k4.dS_dt), 0.0f, 1.0f);
    state_.instability  = std::clamp(state_.instability  + (dt / 6.0f) * (k1.dI_dt + 2.0f * k2.dI_dt + 2.0f * k3.dI_dt + k4.dI_dt), 0.0f, 1.0f);
    state_.p53_activity = std::clamp(state_.p53_activity + (dt / 6.0f) * (k1.dp53_dt + 2.0f * k2.dp53_dt + 2.0f * k3.dp53_dt + k4.dp53_dt), 0.0f, 1.0f);
    state_.p21_activity = std::clamp(state_.p21_activity + (dt / 6.0f) * (k1.dp21_dt + 2.0f * k2.dp21_dt + 2.0f * k3.dp21_dt + k4.dp21_dt), 0.0f, 1.0f);

    current_time_ += dt;
    update_r_score();

    // Guardar en historial circular
    r_history_.push_back(state_.r_score);
    if (r_history_.size() > 60) {
        r_history_.erase(r_history_.begin());
    }
}

} // namespace genome::yamanaka
