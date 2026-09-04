#include "workflow/workflow_engine.hpp"
#include <algorithm>

namespace genome::workflow {

ClinicalWorkflow::ClinicalWorkflow(std::string id, std::string name, std::string category, std::string description)
    : id_(std::move(id)), name_(std::move(name)), category_(std::move(category)), description_(std::move(description)) {}

void ClinicalWorkflow::add_step(const WorkflowStep& step) {
    steps_.push_back(step);
}

bool ClinicalWorkflow::is_completed() const {
    if (steps_.empty()) return false;
    for (const auto& step : steps_) {
        if (step.status != StepStatus::Completed && step.status != StepStatus::Skipped) {
            return false;
        }
    }
    return true;
}

float ClinicalWorkflow::get_progress_percent() const {
    if (steps_.empty()) return 0.0f;
    size_t completed = 0;
    for (const auto& step : steps_) {
        if (step.status == StepStatus::Completed || step.status == StepStatus::Skipped) {
            completed++;
        }
    }
    return (static_cast<float>(completed) / static_cast<float>(steps_.size())) * 100.0f;
}

void ClinicalWorkflow::reset() {
    current_step_index_ = 0;
    is_running_ = false;
    for (auto& step : steps_) {
        step.status = StepStatus::Pending;
        step.result_summary.clear();
    }
}

WorkflowStep* ClinicalWorkflow::get_current_step() {
    if (current_step_index_ < steps_.size()) {
        return &steps_[current_step_index_];
    }
    return nullptr;
}

WorkflowEngine::WorkflowEngine() {
    init_default_protocols();
}

ClinicalWorkflow* WorkflowEngine::get_active_workflow() {
    if (active_workflow_index_ < workflows_.size()) {
        return &workflows_[active_workflow_index_];
    }
    return nullptr;
}

void WorkflowEngine::select_workflow(size_t index) {
    if (index < workflows_.size()) {
        active_workflow_index_ = index;
    }
}

void WorkflowEngine::reset_active_workflow() {
    auto* wf = get_active_workflow();
    if (wf) wf->reset();
}

void WorkflowEngine::init_default_protocols() {
    workflows_.clear();

    // -------------------------------------------------------------
    // PROTOCOLO 1: ONCOLOGÍA DE PRECISIÓN INTEGRAL
    // -------------------------------------------------------------
    {
        ClinicalWorkflow wf(
            "WF_ONCO_FULL",
            "Oncologia de Precision Integral",
            "ONCOLOGIA CLINICA",
            "Protocolo multiobjetivo para evaluacion de drivers tumorales, prescripcion FDA y respuesta in-silico."
        );

        WorkflowStep s1;
        s1.id = "S1_LOAD";
        s1.title = "Ingesta de Muestra Oncologica";
        s1.description = "Cargar panel multiobjetivo de secuenciacion genomica de alta cobertura.";
        s1.type = StepType::LoadVcfSample;
        s1.param_file = "data/sample_patient_oncology.vcf";
        wf.add_step(s1);

        WorkflowStep s2;
        s2.id = "S2_DIAG";
        s2.title = "Diagnostico Molecular y Drivers";
        s2.description = "Identificar alteraciones patogenicas en EGFR, BRCA1, KRAS, BRAF, PIK3CA y TP53.";
        s2.type = StepType::RunGenomicDiagnosis;
        wf.add_step(s2);

        WorkflowStep s3;
        s3.id = "S3_SCREEN";
        s3.title = "Prescripcion de Terapia Dirigida";
        s3.description = "Emparejar alteraciones contra base de datos FDA Nivel 1 y NCCN Nivel 2.";
        s3.type = StepType::ScreenTargetedTherapies;
        wf.add_step(s3);

        WorkflowStep s4;
        s4.id = "S4_DRUG";
        s4.title = "Inoculacion Farmacologica In-Silico";
        s4.description = "Administrar Osimertinib (1.50 uM) en la placa de cultivo del Gemelo Digital.";
        s4.type = StepType::InoculateDrug;
        s4.param_drug = "TargetedTKI";
        s4.param_dose_uM = 1.50f;
        wf.add_step(s4);

        WorkflowStep s5;
        s5.id = "S5_SIM";
        s5.title = "Simulacion de Regresion Tumoral";
        s5.description = "Monitorear citotoxicidad selectiva y disminucion de carga neoplásica (40 ticks).";
        s5.type = StepType::SimulateBioTime;
        s5.param_ticks = 40;
        wf.add_step(s5);

        WorkflowStep s6;
        s6.id = "S6_RESIST";
        s6.title = "Auditoria de Resistencia Clonal";
        s6.description = "Evaluar confluencia de clones con mutaciones adaptativas de escape farmacologico.";
        s6.type = StepType::EvaluateResistance;
        wf.add_step(s6);

        WorkflowStep s7;
        s7.id = "S7_REPORT";
        s7.title = "Emision de Reporte Clinico MTB";
        s7.description = "Generar dictamen para Comite Molecular de Tumores en formatos Markdown y JSON.";
        s7.type = StepType::ExportClinicalReports;
        s7.param_file = "data/molecular_tumor_board_report.md";
        wf.add_step(s7);

        workflows_.push_back(std::move(wf));
    }

    // -------------------------------------------------------------
    // PROTOCOLO 2: RESCATE DE RESISTENCIA TKI (T790M)
    // -------------------------------------------------------------
    {
        ClinicalWorkflow wf(
            "WF_RESIST_RESCUE",
            "Rescate de Resistencia TKI (T790M)",
            "FARMACOGENOMICA",
            "Deteccion de mutacion gatekeeper en EGFR, contraindicacion de Gefitinib y transicion a Osimertinib."
        );

        WorkflowStep s1;
        s1.id = "S1_LOAD_RESIST";
        s1.title = "Ingesta de Panel con Resistencia";
        s1.description = "Cargar VCF de paciente en progresion post-TKI de primera generacion.";
        s1.type = StepType::LoadVcfSample;
        s1.param_file = "data/sample_patient_resistance.vcf";
        wf.add_step(s1);

        WorkflowStep s2;
        s2.id = "S2_GATEKEEPER";
        s2.title = "Deteccion de Mutacion Gatekeeper";
        s2.description = "Identificar p.T790M en EGFR e impedimento esterico en el bolsillo de anclaje de ATP.";
        s2.type = StepType::RunGenomicDiagnosis;
        wf.add_step(s2);

        WorkflowStep s3;
        s3.id = "S3_TKI_SWITCH";
        s3.title = "Transicion a TKI de 3ra Generacion";
        s3.description = "Inocular Osimertinib (2.00 uM) para superar la resistencia esterica de T790M.";
        s3.type = StepType::InoculateDrug;
        s3.param_drug = "TargetedTKI";
        s3.param_dose_uM = 2.00f;
        wf.add_step(s3);

        WorkflowStep s4;
        s4.id = "S4_CRISPR";
        s4.title = "Rescate Genico CRISPR HDR de TP53";
        s4.description = "Electroporacion in-silico de ssODN para restaurar el guardian del genoma.";
        s4.type = StepType::ApplyCrisprRescue;
        s4.param_locus = "TP53";
        wf.add_step(s4);

        WorkflowStep s5;
        s5.id = "S5_SIM";
        s5.title = "Monitoreo de Cinética de Respuesta";
        s5.description = "Simular viabilidad celular y reduccion del R-Score poblacional (35 ticks).";
        s5.type = StepType::SimulateBioTime;
        s5.param_ticks = 35;
        wf.add_step(s5);

        WorkflowStep s6;
        s6.id = "S6_EXPORT";
        s6.title = "Informe de Farmacogenomica y Alerta";
        s6.description = "Exportar informe con contraindicaciones clinicas registradas.";
        s6.type = StepType::ExportClinicalReports;
        s6.param_file = "data/resistance_tumor_board_report.md";
        wf.add_step(s6);

        workflows_.push_back(std::move(wf));
    }

    // -------------------------------------------------------------
    // PROTOCOLO 3: LETALIDAD SINTÉTICA (BRCA1 / PARP)
    // -------------------------------------------------------------
    {
        ClinicalWorkflow wf(
            "WF_PARP_SYNTHETIC",
            "Letalidad Sintetica (BRCA1 / PARP)",
            "TERAPIA BLANCO",
            "Explotacion de deficiencia en recombinacion homologa mediante inhibidores de PARP1/2."
        );

        WorkflowStep s1;
        s1.id = "S1_LOAD";
        s1.title = "Ingesta de Muestra con Deficiencia HR";
        s1.description = "Cargar perfil tumoral con variante patogenica p.C61G en gen BRCA1.";
        s1.type = StepType::LoadVcfSample;
        s1.param_file = "data/sample_patient_oncology.vcf";
        wf.add_step(s1);

        WorkflowStep s2;
        s2.id = "S2_HR_SCREEN";
        s2.title = "Screening de Deficiencia de Recombinacion";
        s2.description = "Verificar vulnerabilidad sintetica por colapso de reparacion de roturas de doble hebra.";
        s2.type = StepType::RunGenomicDiagnosis;
        wf.add_step(s2);

        WorkflowStep s3;
        s3.id = "S3_INOCULATE_PARP";
        s3.title = "Infusion de Inhibidor PARP (Olaparib)";
        s3.description = "Inocular Olaparib (1.50 uM) para inducir atrapamiento de PARP1/2 sobre el ADN.";
        s3.type = StepType::InoculateDrug;
        s3.param_drug = "ParpInhibitor";
        s3.param_dose_uM = 1.50f;
        wf.add_step(s3);

        WorkflowStep s4;
        s4.id = "S4_COLLAPSE";
        s4.title = "Simulacion de Colapso Replicativo";
        s4.description = "Monitorear muerte celular masiva en clones deficientes en BRCA1 (45 ticks).";
        s4.type = StepType::SimulateBioTime;
        s4.param_ticks = 45;
        wf.add_step(s4);

        WorkflowStep s5;
        s5.id = "S5_REPORT";
        s5.title = "Exportar Informe de Letalidad Sintetica";
        s5.description = "Generar documentacion clinica de respuesta terapeutica.";
        s5.type = StepType::ExportClinicalReports;
        s5.param_file = "data/molecular_tumor_board_report.md";
        wf.add_step(s5);

        workflows_.push_back(std::move(wf));
    }

    // -------------------------------------------------------------
    // PROTOCOLO 4: REPROGRAMACIÓN Y RESCATE OSKM
    // -------------------------------------------------------------
    {
        ClinicalWorkflow wf(
            "WF_OSKM_REPROGRAM",
            "Reprogramacion y Rescate OSKM",
            "MEDICINA REGENERATIVA",
            "Induccion de factores de Yamanaka, deteccion de colisiones en supresores y rescate CRISPR."
        );

        WorkflowStep s1;
        s1.id = "S1_LOAD";
        s1.title = "Ingesta de Muestra Somatica Basal";
        s1.description = "Cargar genoma fibroblástico del paciente portador de variantes.";
        s1.type = StepType::LoadVcfSample;
        s1.param_file = "data/sample_patient.vcf";
        wf.add_step(s1);

        WorkflowStep s2;
        s2.id = "S2_DIAG";
        s2.title = "Evaluacion de Colisiones con Factores OSKM";
        s2.description = "Identificar conflictos de reprogramacion en locus TP53 y CDKN2A.";
        s2.type = StepType::RunGenomicDiagnosis;
        wf.add_step(s2);

        WorkflowStep s3;
        s3.id = "S3_CRISPR";
        s3.title = "Edicion Genica CRISPR HDR de Precision";
        s3.description = "Sustituir codon mutante p.R175H restituyendo apoptosis celular previa a la induccion.";
        s3.type = StepType::ApplyCrisprRescue;
        s3.param_locus = "TP53";
        wf.add_step(s3);

        WorkflowStep s4;
        s4.id = "S4_SIM";
        s4.title = "Simulacion de Transicion hacia Pluripotencia";
        s4.description = "Resolver cinetica ODE de Yamanaka y evaluar emergencia de colonias iPSC (50 ticks).";
        s4.type = StepType::SimulateBioTime;
        s4.param_ticks = 50;
        wf.add_step(s4);

        WorkflowStep s5;
        s5.id = "S5_EXPORT";
        s5.title = "Dictamen de Seguridad de Celulas Madre";
        s5.description = "Emitir dictamen de estabilidad genomica y R-Score seguro (< 0.15).";
        s5.type = StepType::ExportClinicalReports;
        s5.param_file = "data/molecular_tumor_board_report.md";
        wf.add_step(s5);

        workflows_.push_back(std::move(wf));
    }

    // -------------------------------------------------------------
    // PROTOCOLO 5: DIAGNÓSTICO GENÓMICO RÁPIDO
    // -------------------------------------------------------------
    {
        ClinicalWorkflow wf(
            "WF_EXPRESS_SCREEN",
            "Cribado Diagnostico Rapido (Express)",
            "DIAGNOSTICO AMBULATORIO",
            "Analisis directo de carga mutacional, calculo de R-Score y prescripcion inicial inmediata."
        );

        WorkflowStep s1;
        s1.id = "S1_LOAD";
        s1.title = "Carga Express de Muestra VCF";
        s1.description = "Mapeo ultra-rapido de coordenadas GRCh38 en memoria ARM64 unificada.";
        s1.type = StepType::LoadVcfSample;
        s1.param_file = "data/sample_patient_oncology.vcf";
        wf.add_step(s1);

        WorkflowStep s2;
        s2.id = "S2_TMB";
        s2.title = "Cuantificacion de TMB y Riesgo";
        s2.description = "Calcular Carga Mutacional Tumoral (mut/Mb) e idoneidad para inmunoterapia.";
        s2.type = StepType::RunGenomicDiagnosis;
        wf.add_step(s2);

        WorkflowStep s3;
        s3.id = "S3_TARGETS";
        s3.title = "Matching Instantaneo de Dianas";
        s3.description = "Clasificar terapias disponibles por niveles de evidencia FDA y NCCN.";
        s3.type = StepType::ScreenTargetedTherapies;
        wf.add_step(s3);

        WorkflowStep s4;
        s4.id = "S4_OUT";
        s4.title = "Generacion Instantanea de Informe Clinico";
        s4.description = "Exportar resumen molecular para discusion clinica.";
        s4.type = StepType::ExportClinicalReports;
        s4.param_file = "data/molecular_tumor_board_report.md";
        wf.add_step(s4);

        workflows_.push_back(std::move(wf));
    }
}

std::string WorkflowEngine::step_type_to_string(StepType type) {
    switch (type) {
        case StepType::LoadVcfSample:           return "Carga de Muestra VCF";
        case StepType::RunGenomicDiagnosis:     return "Diagnostico Molecular";
        case StepType::ScreenTargetedTherapies: return "Cribado de Terapias Dirigidas";
        case StepType::InoculateDrug:           return "Inoculacion Farmacologica";
        case StepType::SimulateBioTime:         return "Simulacion en Gemelo Digital";
        case StepType::EvaluateResistance:      return "Auditoria de Resistencia";
        case StepType::ApplyCrisprRescue:       return "Rescate Genico CRISPR HDR";
        case StepType::ExportClinicalReports:   return "Exportacion de Reportes";
    }
}

std::string WorkflowEngine::step_status_to_string(StepStatus status) {
    switch (status) {
        case StepStatus::Pending:   return "[ PENDIENTE ]";
        case StepStatus::Running:   return "[ EN CURSO ]";
        case StepStatus::Completed: return "[ COMPLETADO ]";
        case StepStatus::Failed:    return "[ FALLIDO ]";
        case StepStatus::Skipped:   return "[ OMITIDO ]";
    }
}

} // namespace genome::workflow
