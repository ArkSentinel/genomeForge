#pragma once

#include <string>
#include <vector>
#include <cstdint>
#include <memory>

struct sqlite3;

namespace genome::storage {

struct SimulationRecord {
    int64_t id{0};
    std::string timestamp;
    std::string sample_name;
    std::string vcf_path;
    size_t total_variants{0};
    size_t pathogenic_variants{0};
    float r_score{0.0f};
    float confluence{0.0f};
    size_t living_cells{0};
    std::string active_drug{"None"};
    float tumor_regression_pct{0.0f};
    size_t resistant_cells{0};
    std::string notes;
};

struct ReportRecord {
    int64_t id{0};
    std::string timestamp;
    std::string sample_id;
    float tmb_score{0.0f};
    float r_score{0.0f};
    size_t biomarkers_count{0};
    size_t level1_count{0};
    std::string report_md_path;
    std::string report_json_path;
};

struct WorkflowRunRecord {
    int64_t id{0};
    std::string timestamp;
    std::string workflow_id;
    std::string workflow_name;
    size_t total_steps{0};
    size_t completed_steps{0};
    std::string status;
    std::string summary_findings;
};

class SQLiteStorage {
public:
    SQLiteStorage();
    explicit SQLiteStorage(const std::string& db_path);
    ~SQLiteStorage();

    // Apertura y cierre de base de datos
    bool open(const std::string& db_path = "data/genome_forge.db");
    void close();
    bool is_open() const { return db_ != nullptr; }
    const std::string& get_db_path() const { return db_path_; }

    // Inicializacion de esquema relacional
    bool init_schema();

    // Inserciones estructuradas
    int64_t record_simulation(const SimulationRecord& rec);
    int64_t record_clinical_report(const ReportRecord& rec);
    int64_t record_workflow_run(const WorkflowRunRecord& rec);

    // Consultas
    std::vector<SimulationRecord> get_recent_simulations(size_t limit = 10);
    std::vector<ReportRecord> get_recent_reports(size_t limit = 10);
    std::vector<WorkflowRunRecord> get_recent_workflows(size_t limit = 10);

    // Salida informativa en terminal para CLI / Inspector
    void print_database_summary() const;

    // Helper para fecha/hora ISO8601 actual
    static std::string current_iso_timestamp();

private:
    sqlite3* db_{nullptr};
    std::string db_path_{"data/genome_forge.db"};
};

} // namespace genome::storage
