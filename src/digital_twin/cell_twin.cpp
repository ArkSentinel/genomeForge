#include "digital_twin/cell_twin.hpp"
#include "metal/metal_compute.hpp"

#include <cmath>
#include <algorithm>
#include <iostream>

namespace genome::twin {

CellTwinSimulator::CellTwinSimulator()
    : rng_(1337) {
    initialize_dish(180, true);
}

CellTwinSimulator::~CellTwinSimulator() = default;

void CellTwinSimulator::initialize_dish(size_t initial_count, bool patient_has_tp53_mut) {
    cells_.clear();
    next_cell_id_ = 1;
    dox_active_ = false;
    crispr_applied_ = false;
    uv_insult_active_ = false;
    uv_timer_ = 0.0f;

    std::uniform_real_distribution<float> dist_r(0.0f, dish_radius_ * 0.75f);
    std::uniform_real_distribution<float> dist_theta(0.0f, 2.0f * M_PI);
    std::uniform_real_distribution<float> dist_phase(0.0f, 1.0f);
    std::uniform_real_distribution<float> dist_vel(-0.5f, 0.5f);
    std::bernoulli_distribution dist_mut(patient_has_tp53_mut ? 0.65 : 0.0);

    cells_.reserve(initial_count * 2);
    for (size_t i = 0; i < initial_count; ++i) {
        CellAgent c{};
        float r = std::sqrt(dist_phase(rng_)) * (dish_radius_ * 0.75f);
        float theta = dist_theta(rng_);
        c.x = r * std::cos(theta);
        c.y = r * std::sin(theta);
        c.vx = dist_vel(rng_);
        c.vy = dist_vel(rng_);
        c.radius = 2.2f;
        c.mass = 1.0f;
        c.phase = static_cast<CellCyclePhase>(static_cast<int>(dist_phase(rng_) * 4.0f));
        c.phenotype = CellPhenotype::SomaticFibroblast;
        c.cycle_clock = dist_phase(rng_);
        c.P_score = 0.05f;
        c.S_score = 0.95f;
        c.cell_id = next_cell_id_++;
        c.generation = 0;

        if (dist_mut(rng_)) {
            c.has_tp53_mutation = 1;
            c.p53_level = 0.15f; // Pérdida de función p53
            c.r_score = 0.60f;
        } else {
            c.has_tp53_mutation = 0;
            c.p53_level = 0.85f; // p53 silvestre activa
            c.r_score = 0.08f;
        }

        cells_.push_back(c);
    }

    compute_population_metrics();
}

void CellTwinSimulator::reset() {
    initialize_dish(180, true);
}

void CellTwinSimulator::toggle_doxycycline() {
    dox_active_ = !dox_active_;
}

void CellTwinSimulator::apply_crispr_rescue() {
    crispr_applied_ = true;
    for (auto& c : cells_) {
        if (c.has_tp53_mutation) {
            c.has_tp53_mutation = 0;
            c.p53_level = 0.88f; // Restaurar p53 funcional
            c.r_score = 0.12f;
            if (c.phenotype == CellPhenotype::TransformedMalignant) {
                c.is_transformed = 0;
                // Si la célula transformada ya reprogramó hacia pluripotencia, se estabiliza como iPSC
                if (c.P_score > 0.60f) {
                    c.phenotype = CellPhenotype::PluripotentIPSC;
                } else {
                    c.phenotype = CellPhenotype::IntermediateReprog;
                }
            }
        }
    }
}

void CellTwinSimulator::trigger_uv_damage() {
    uv_insult_active_ = true;
    uv_timer_ = 6.0f; // 6 segundos de estrés genotóxico
}

void CellTwinSimulator::adjust_population(int delta) {
    if (delta < 0 && cells_.size() > 20) {
        size_t to_remove = std::min(static_cast<size_t>(-delta), cells_.size() - 20);
        cells_.erase(cells_.end() - to_remove, cells_.end());
    } else if (delta > 0 && cells_.size() < 3000) {
        std::uniform_real_distribution<float> dist_r(0.0f, dish_radius_ * 0.7f);
        std::uniform_real_distribution<float> dist_theta(0.0f, 2.0f * M_PI);
        for (int i = 0; i < delta; ++i) {
            CellAgent c{};
            float r = dist_r(rng_);
            float theta = dist_theta(rng_);
            c.x = r * std::cos(theta);
            c.y = r * std::sin(theta);
            c.radius = 2.2f;
            c.mass = 1.0f;
            c.phase = CellCyclePhase::G1_Growth;
            c.phenotype = CellPhenotype::SomaticFibroblast;
            c.cycle_clock = 0.1f;
            c.cell_id = next_cell_id_++;
            c.p53_level = 0.85f;
            cells_.push_back(c);
        }
    }
    compute_population_metrics();
}

void CellTwinSimulator::update_physics_cpu(float dt) {
    size_t count = cells_.size();
    float k_repulsion = 14.0f;

    for (size_t i = 0; i < count; ++i) {
        auto& c1 = cells_[i];
        if (c1.phase == CellCyclePhase::Apoptotic) continue;

        float fx = 0.0f;
        float fy = 0.0f;

        for (size_t j = 0; j < count; ++j) {
            if (i == j) continue;
            auto& c2 = cells_[j];
            float dx = c1.x - c2.x;
            float dy = c1.y - c2.y;
            float dist_sq = dx * dx + dy * dy;
            float min_dist = c1.radius + c2.radius;

            if (dist_sq < min_dist * min_dist && dist_sq > 0.0001f) {
                float dist = std::sqrt(dist_sq);
                float overlap = min_dist - dist;
                float force = overlap * k_repulsion;
                fx += (dx / dist) * force;
                fy += (dy / dist) * force;
            }

            // Adhesión clonal en colonias iPSC (Green)
            if (c1.phenotype == CellPhenotype::PluripotentIPSC && c2.phenotype == CellPhenotype::PluripotentIPSC) {
                if (dist_sq < (min_dist * 2.5f) * (min_dist * 2.5f) && dist_sq > min_dist * min_dist) {
                    float dist = std::sqrt(dist_sq);
                    float attraction = 1.5f;
                    fx -= (dx / dist) * attraction;
                    fy -= (dy / dist) * attraction;
                }
            }
        }

        // Fricción del medio de cultivo
        fx -= c1.vx * 2.2f;
        fy -= c1.vy * 2.2f;

        c1.vx += fx * dt;
        c1.vy += fy * dt;
        c1.x += c1.vx * dt;
        c1.y += c1.vy * dt;

        // Borde circular de la placa de Petri
        float r_sq = c1.x * c1.x + c1.y * c1.y;
        float max_r = dish_radius_ - c1.radius;
        if (r_sq > max_r * max_r) {
            float r = std::sqrt(r_sq);
            c1.x = (c1.x / r) * max_r;
            c1.y = (c1.y / r) * max_r;
            c1.vx = -c1.vx * 0.4f;
            c1.vy = -c1.vy * 0.4f;
        }
    }
}

void CellTwinSimulator::update_cell_cycle(CellAgent& cell, float dt, std::vector<CellAgent>& new_daughters) {
    if (cell.phase == CellCyclePhase::Apoptotic) {
        cell.radius *= 0.985f; // Encogimiento apoptótico
        cell.cycle_clock += dt * 0.5f;
        return;
    }

    // Estrés genotóxico UV
    if (uv_insult_active_) {
        if (!cell.has_tp53_mutation) {
            // p53 funcional induce arresto o apoptosis protectora
            cell.p53_level = std::min(1.0f, cell.p53_level + dt * 0.4f);
            if (cell.p53_level > 0.92f) {
                cell.phase = CellCyclePhase::Apoptotic;
                cell.phenotype = CellPhenotype::ApoptoticCorpse;
                return;
            }
        } else {
            // Célula con TP53 mutado no responde a daño UV y acumula mutaciones
            cell.r_score = std::min(1.0f, cell.r_score + dt * 0.3f);
        }
    }

    // Cinética de Reprogramación con Factores de Yamanaka (Doxiciclina)
    if (dox_active_) {
        // Avance en la reprogramación celular
        cell.P_score = std::min(1.0f, cell.P_score + dt * 0.08f);
        cell.S_score = std::max(0.0f, cell.S_score - dt * 0.08f);

        if (cell.P_score > 0.35f && cell.P_score < 0.75f) {
            cell.phenotype = CellPhenotype::IntermediateReprog;
        } else if (cell.P_score >= 0.75f && !cell.is_transformed) {
            cell.phenotype = CellPhenotype::PluripotentIPSC;
        }

        // Colisión oncogénica ante sobreexpresión de c-Myc si TP53 está mutado
        if (cell.has_tp53_mutation) {
            cell.r_score = std::min(1.0f, cell.r_score + dt * 0.15f);
            if (cell.r_score > 0.70f) {
                cell.is_transformed = 1;
                cell.phenotype = CellPhenotype::TransformedMalignant;
            }
        } else {
            // En células normales, p53 previene el escape oncogénico
            cell.r_score = std::max(0.05f, cell.r_score - dt * 0.05f);
        }
    } else {
        // Sin doxiciclina, las células en estado intermedio regresan lentamente a identidad diferenciada
        if (cell.phenotype == CellPhenotype::IntermediateReprog) {
            cell.P_score = std::max(0.05f, cell.P_score - dt * 0.04f);
            cell.S_score = std::min(0.95f, cell.S_score + dt * 0.04f);
            if (cell.P_score < 0.20f) {
                cell.phenotype = CellPhenotype::SomaticFibroblast;
            }
        }
    }

    // Progresión del ciclo celular
    float speed_factor = 0.12f;
    if (cell.is_transformed) {
        speed_factor = 0.35f; // Las células cancerosas proliferan 3 veces más rápido
    } else if (cell.phenotype == CellPhenotype::PluripotentIPSC) {
        speed_factor = 0.20f; // iPSCs tienen fase G1 corta
    }

    cell.cycle_clock += dt * speed_factor;

    // Transición de fases G1 -> S -> G2 -> M
    if (cell.cycle_clock >= 1.0f) {
        cell.cycle_clock = 0.0f;
        switch (cell.phase) {
            case CellCyclePhase::G1_Growth:
                cell.phase = CellCyclePhase::S_Synthesis;
                break;
            case CellCyclePhase::S_Synthesis:
                cell.phase = CellCyclePhase::G2_PreMitotic;
                break;
            case CellCyclePhase::G2_PreMitotic:
                cell.phase = CellCyclePhase::M_Mitosis;
                break;
            case CellCyclePhase::M_Mitosis: {
                // Evento de Mitosis: División celular
                cell.phase = CellCyclePhase::G1_Growth;
                cell.generation++;

                // Solo divide si la placa no está saturada (>1500 células)
                if (cells_.size() < 1200) {
                    CellAgent daughter = cell;
                    daughter.cell_id = next_cell_id_++;
                    daughter.parent_id = cell.cell_id;
                    daughter.cycle_clock = 0.05f;

                    // Pequeña dispersión espacial
                    float angle = (cell.cell_id % 360) * (M_PI / 180.0f);
                    float offset = cell.radius * 1.5f;
                    daughter.x = cell.x + std::cos(angle) * offset;
                    daughter.y = cell.y + std::sin(angle) * offset;
                    cell.x -= std::cos(angle) * offset;
                    cell.y -= std::sin(angle) * offset;

                    new_daughters.push_back(daughter);
                }
                break;
            }
            default:
                cell.phase = CellCyclePhase::G1_Growth;
                break;
        }
    }
}

void CellTwinSimulator::step(float dt) {
    if (uv_insult_active_) {
        uv_timer_ -= dt;
        if (uv_timer_ <= 0.0f) {
            uv_insult_active_ = false;
        }
    }

    // 1. Biofísica intercelular
    update_physics_cpu(dt);

    // 2. Ciclo celular y proliferación
    std::vector<CellAgent> new_daughters;
    for (auto& cell : cells_) {
        update_cell_cycle(cell, dt, new_daughters);
    }

    if (!new_daughters.empty()) {
        cells_.insert(cells_.end(), new_daughters.begin(), new_daughters.end());
    }

    // 3. Limpieza de células apoptóticas completamente degradadas
    cells_.erase(
        std::remove_if(cells_.begin(), cells_.end(), [](const CellAgent& c) {
            return c.phase == CellCyclePhase::Apoptotic && (c.cycle_clock > 1.5f || c.radius < 0.6f);
        }),
        cells_.end()
    );

    compute_population_metrics();
}

void CellTwinSimulator::compute_population_metrics() {
    metrics_.total_cells = cells_.size();
    metrics_.somatic_count = 0;
    metrics_.intermediate_count = 0;
    metrics_.ipsc_count = 0;
    metrics_.transformed_count = 0;
    metrics_.apoptotic_count = 0;

    size_t mitotic_count = 0;
    float sum_r_score = 0.0f;

    for (const auto& c : cells_) {
        switch (c.phenotype) {
            case CellPhenotype::SomaticFibroblast: metrics_.somatic_count++; break;
            case CellPhenotype::IntermediateReprog: metrics_.intermediate_count++; break;
            case CellPhenotype::PluripotentIPSC: metrics_.ipsc_count++; break;
            case CellPhenotype::TransformedMalignant: metrics_.transformed_count++; break;
            case CellPhenotype::ApoptoticCorpse: metrics_.apoptotic_count++; break;
        }

        if (c.phase == CellCyclePhase::M_Mitosis) {
            mitotic_count++;
        }

        sum_r_score += c.r_score;
    }

    if (!cells_.empty()) {
        metrics_.mitotic_index = (static_cast<float>(mitotic_count) / static_cast<float>(cells_.size())) * 100.0f;
        metrics_.average_r_score = sum_r_score / static_cast<float>(cells_.size());
        float dish_area = M_PI * dish_radius_ * dish_radius_;
        float cells_area = cells_.size() * (M_PI * 2.2f * 2.2f);
        metrics_.confluence = std::min(1.0f, cells_area / dish_area);
        metrics_.ipsc_efficiency = (static_cast<float>(metrics_.ipsc_count) / static_cast<float>(cells_.size())) * 100.0f;
    } else {
        metrics_.mitotic_index = 0.0f;
        metrics_.average_r_score = 0.0f;
        metrics_.confluence = 0.0f;
        metrics_.ipsc_efficiency = 0.0f;
    }

    metrics_.dox_active = dox_active_;
    metrics_.crispr_applied = crispr_applied_;
    metrics_.uv_insult_active = uv_insult_active_;
}

} // namespace genome::twin
