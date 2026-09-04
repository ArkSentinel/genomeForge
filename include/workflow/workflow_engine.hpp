#pragma once

#include <string>
#include <vector>
#include <memory>
#include <functional>
#include <cstdint>

namespace genome::workflow {

enum class StepType : uint8_t {
    LoadVcfSample = 0,
    RunGenomicDiagnosis = 1,
    ScreenTargetedTherapies = 2,
    InoculateDrug = 3,
    SimulateBioTime = 4,
    EvaluateResistance = 5,
    ApplyCrisprRescue = 6,
    ExportClinicalReports = 7
};

enum class StepStatus : uint8_t {
    Pending = 0,
    Running = 1,
    Completed = 2,
    Failed = 3,
    Skipped = 4
};

struct WorkflowStep {
    std::string id;
    std::string title;
    std::string description;
    std::string details;
    StepType type{StepType::LoadVcfSample};
    StepStatus status{StepStatus::Pending};
    
    // Parámetros de ejecución clínica
    std::string param_file;
    std::string param_drug;
    float param_dose_uM{1.50f};
    size_t param_ticks{30};
    std::string param_locus{"TP53"};
    
    // Resultados observados tras ejecución
    std::string result_summary;
};

class ClinicalWorkflow {
public:
    ClinicalWorkflow(std::string id, std::string name, std::string category, std::string description);

    const std::string& get_id() const { return id_; }
    const std::string& get_name() const { return name_; }
    const std::string& get_category() const { return category_; }
    const std::string& get_description() const { return description_; }

    void add_step(const WorkflowStep& step);
    const std::vector<WorkflowStep>& get_steps() const { return steps_; }
    std::vector<WorkflowStep>& get_steps_mut() { return steps_; }

    size_t get_current_step_index() const { return current_step_index_; }
    void set_current_step_index(size_t index) { current_step_index_ = index; }

    bool is_completed() const;
    bool is_running() const { return is_running_; }
    void set_running(bool running) { is_running_ = running; }

    float get_progress_percent() const;
    float get_progress_ratio() const { return get_progress_percent() / 100.0f; }
    void reset();

    WorkflowStep* get_current_step();

private:
    std::string id_;
    std::string name_;
    std::string category_;
    std::string description_;
    std::vector<WorkflowStep> steps_;
    size_t current_step_index_{0};
    bool is_running_{false};
};

class WorkflowEngine {
public:
    WorkflowEngine();

    // Catálogo de flujos clínicos
    const std::vector<ClinicalWorkflow>& get_workflows() const { return workflows_; }
    ClinicalWorkflow* get_active_workflow();
    size_t get_active_workflow_index() const { return active_workflow_index_; }
    size_t get_active_index() const { return active_workflow_index_; }
    void select_workflow(size_t index);

    // Reinicia el flujo activo
    void reset_active_workflow();

    // Genera la lista inicial de flujos clínicos predefinidos
    void init_default_protocols();

    // Helper para descripción textual de estados y tipos
    static std::string step_type_to_string(StepType type);
    static std::string step_status_to_string(StepStatus status);

private:
    std::vector<ClinicalWorkflow> workflows_;
    size_t active_workflow_index_{0};
};

} // namespace genome::workflow
