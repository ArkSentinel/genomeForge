#include "storage/sqlite_storage.hpp"

#include <sqlite3.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <chrono>
#include <ctime>
#include <filesystem>

namespace genome::storage {

std::string SQLiteStorage::current_iso_timestamp() {
    auto now = std::chrono::system_clock::now();
    std::time_t now_c = std::chrono::system_clock::to_time_t(now);
    std::tm tm_buf{};
#if defined(_WIN32)
    gmtime_s(&tm_buf, &now_c);
#else
    gmtime_r(&now_c, &tm_buf);
#endif
    char buf[64];
    std::strftime(buf, sizeof(buf), "%Y-%m-%dT%H:%M:%SZ", &tm_buf);
    return std::string(buf);
}

SQLiteStorage::SQLiteStorage() {
    open("data/genome_forge.db");
}

SQLiteStorage::SQLiteStorage(const std::string& db_path) {
    open(db_path);
}

SQLiteStorage::~SQLiteStorage() {
    close();
}

bool SQLiteStorage::open(const std::string& db_path) {
    close();
    db_path_ = db_path;

    // Crear directorio data/ si no existe
    try {
        std::filesystem::path p(db_path_);
        if (p.has_parent_path()) {
            std::filesystem::create_directories(p.parent_path());
        }
    } catch (...) {}

    int rc = sqlite3_open_v2(db_path_.c_str(), &db_, 
                             SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE, nullptr);
    if (rc != SQLITE_OK) {
        std::cerr << "[!] Error abriendo SQLite DB en " << db_path_ << ": " 
                  << sqlite3_errmsg(db_) << "\n";
        db_ = nullptr;
        return false;
    }

    return init_schema();
}

void SQLiteStorage::close() {
    if (db_) {
        sqlite3_close_v2(db_);
        db_ = nullptr;
    }
}

bool SQLiteStorage::init_schema() {
    if (!db_) return false;

    const char* schema_sql = R"(
        CREATE TABLE IF NOT EXISTS simulations (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            timestamp TEXT NOT NULL,
            sample_name TEXT NOT NULL,
            vcf_path TEXT NOT NULL,
            total_variants INTEGER NOT NULL,
            pathogenic_variants INTEGER NOT NULL,
            r_score REAL NOT NULL,
            confluence REAL NOT NULL,
            living_cells INTEGER NOT NULL,
            active_drug TEXT NOT NULL,
            tumor_regression_pct REAL NOT NULL,
            resistant_cells INTEGER NOT NULL,
            notes TEXT
        );

        CREATE TABLE IF NOT EXISTS clinical_reports (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            timestamp TEXT NOT NULL,
            sample_id TEXT NOT NULL,
            tmb_score REAL NOT NULL,
            r_score REAL NOT NULL,
            biomarkers_count INTEGER NOT NULL,
            level1_count INTEGER NOT NULL,
            report_md_path TEXT NOT NULL,
            report_json_path TEXT NOT NULL
        );

        CREATE TABLE IF NOT EXISTS workflow_runs (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            timestamp TEXT NOT NULL,
            workflow_id TEXT NOT NULL,
            workflow_name TEXT NOT NULL,
            total_steps INTEGER NOT NULL,
            completed_steps INTEGER NOT NULL,
            status TEXT NOT NULL,
            summary_findings TEXT
        );

        CREATE INDEX IF NOT EXISTS idx_simulations_timestamp ON simulations(timestamp);
        CREATE INDEX IF NOT EXISTS idx_reports_sample ON clinical_reports(sample_id);
        CREATE INDEX IF NOT EXISTS idx_workflow_name ON workflow_runs(workflow_name);
    )";

    char* errmsg = nullptr;
    int rc = sqlite3_exec(db_, schema_sql, nullptr, nullptr, &errmsg);
    if (rc != SQLITE_OK) {
        std::cerr << "[!] Error creando esquema SQLite: " 
                  << (errmsg ? errmsg : "Desconocido") << "\n";
        sqlite3_free(errmsg);
        return false;
    }

    return true;
}

int64_t SQLiteStorage::record_simulation(const SimulationRecord& rec) {
    if (!db_) return -1;

    const char* sql = R"(
        INSERT INTO simulations (
            timestamp, sample_name, vcf_path, total_variants, pathogenic_variants,
            r_score, confluence, living_cells, active_drug, tumor_regression_pct,
            resistant_cells, notes
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?);
    )";

    sqlite3_stmt* stmt = nullptr;
    if (sqlite3_prepare_v2(db_, sql, -1, &stmt, nullptr) != SQLITE_OK) {
        return -1;
    }

    std::string ts = rec.timestamp.empty() ? current_iso_timestamp() : rec.timestamp;
    sqlite3_bind_text(stmt, 1, ts.c_str(), -1, SQLITE_TRANSIENT);
    sqlite3_bind_text(stmt, 2, rec.sample_name.c_str(), -1, SQLITE_TRANSIENT);
    sqlite3_bind_text(stmt, 3, rec.vcf_path.c_str(), -1, SQLITE_TRANSIENT);
    sqlite3_bind_int64(stmt, 4, static_cast<sqlite3_int64>(rec.total_variants));
    sqlite3_bind_int64(stmt, 5, static_cast<sqlite3_int64>(rec.pathogenic_variants));
    sqlite3_bind_double(stmt, 6, rec.r_score);
    sqlite3_bind_double(stmt, 7, rec.confluence);
    sqlite3_bind_int64(stmt, 8, static_cast<sqlite3_int64>(rec.living_cells));
    sqlite3_bind_text(stmt, 9, rec.active_drug.c_str(), -1, SQLITE_TRANSIENT);
    sqlite3_bind_double(stmt, 10, rec.tumor_regression_pct);
    sqlite3_bind_int64(stmt, 11, static_cast<sqlite3_int64>(rec.resistant_cells));
    sqlite3_bind_text(stmt, 12, rec.notes.c_str(), -1, SQLITE_TRANSIENT);

    int rc = sqlite3_step(stmt);
    sqlite3_finalize(stmt);

    if (rc == SQLITE_DONE) {
        return sqlite3_last_insert_rowid(db_);
    }
    return -1;
}

int64_t SQLiteStorage::record_clinical_report(const ReportRecord& rec) {
    if (!db_) return -1;

    const char* sql = R"(
        INSERT INTO clinical_reports (
            timestamp, sample_id, tmb_score, r_score, biomarkers_count,
            level1_count, report_md_path, report_json_path
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?);
    )";

    sqlite3_stmt* stmt = nullptr;
    if (sqlite3_prepare_v2(db_, sql, -1, &stmt, nullptr) != SQLITE_OK) {
        return -1;
    }

    std::string ts = rec.timestamp.empty() ? current_iso_timestamp() : rec.timestamp;
    sqlite3_bind_text(stmt, 1, ts.c_str(), -1, SQLITE_TRANSIENT);
    sqlite3_bind_text(stmt, 2, rec.sample_id.c_str(), -1, SQLITE_TRANSIENT);
    sqlite3_bind_double(stmt, 3, rec.tmb_score);
    sqlite3_bind_double(stmt, 4, rec.r_score);
    sqlite3_bind_int64(stmt, 5, static_cast<sqlite3_int64>(rec.biomarkers_count));
    sqlite3_bind_int64(stmt, 6, static_cast<sqlite3_int64>(rec.level1_count));
    sqlite3_bind_text(stmt, 7, rec.report_md_path.c_str(), -1, SQLITE_TRANSIENT);
    sqlite3_bind_text(stmt, 8, rec.report_json_path.c_str(), -1, SQLITE_TRANSIENT);

    int rc = sqlite3_step(stmt);
    sqlite3_finalize(stmt);

    if (rc == SQLITE_DONE) {
        return sqlite3_last_insert_rowid(db_);
    }
    return -1;
}

int64_t SQLiteStorage::record_workflow_run(const WorkflowRunRecord& rec) {
    if (!db_) return -1;

    const char* sql = R"(
        INSERT INTO workflow_runs (
            timestamp, workflow_id, workflow_name, total_steps, completed_steps,
            status, summary_findings
        ) VALUES (?, ?, ?, ?, ?, ?, ?);
    )";

    sqlite3_stmt* stmt = nullptr;
    if (sqlite3_prepare_v2(db_, sql, -1, &stmt, nullptr) != SQLITE_OK) {
        return -1;
    }

    std::string ts = rec.timestamp.empty() ? current_iso_timestamp() : rec.timestamp;
    sqlite3_bind_text(stmt, 1, ts.c_str(), -1, SQLITE_TRANSIENT);
    sqlite3_bind_text(stmt, 2, rec.workflow_id.c_str(), -1, SQLITE_TRANSIENT);
    sqlite3_bind_text(stmt, 3, rec.workflow_name.c_str(), -1, SQLITE_TRANSIENT);
    sqlite3_bind_int64(stmt, 4, static_cast<sqlite3_int64>(rec.total_steps));
    sqlite3_bind_int64(stmt, 5, static_cast<sqlite3_int64>(rec.completed_steps));
    sqlite3_bind_text(stmt, 6, rec.status.c_str(), -1, SQLITE_TRANSIENT);
    sqlite3_bind_text(stmt, 7, rec.summary_findings.c_str(), -1, SQLITE_TRANSIENT);

    int rc = sqlite3_step(stmt);
    sqlite3_finalize(stmt);

    if (rc == SQLITE_DONE) {
        return sqlite3_last_insert_rowid(db_);
    }
    return -1;
}

std::vector<SimulationRecord> SQLiteStorage::get_recent_simulations(size_t limit) {
    std::vector<SimulationRecord> list;
    if (!db_) return list;

    const char* sql = "SELECT id, timestamp, sample_name, vcf_path, total_variants, pathogenic_variants, "
                      "r_score, confluence, living_cells, active_drug, tumor_regression_pct, "
                      "resistant_cells, notes FROM simulations ORDER BY id DESC LIMIT ?;";

    sqlite3_stmt* stmt = nullptr;
    if (sqlite3_prepare_v2(db_, sql, -1, &stmt, nullptr) != SQLITE_OK) {
        return list;
    }

    sqlite3_bind_int64(stmt, 1, static_cast<sqlite3_int64>(limit));

    while (sqlite3_step(stmt) == SQLITE_ROW) {
        SimulationRecord rec;
        rec.id = sqlite3_column_int64(stmt, 0);
        rec.timestamp = reinterpret_cast<const char*>(sqlite3_column_text(stmt, 1));
        rec.sample_name = reinterpret_cast<const char*>(sqlite3_column_text(stmt, 2));
        rec.vcf_path = reinterpret_cast<const char*>(sqlite3_column_text(stmt, 3));
        rec.total_variants = static_cast<size_t>(sqlite3_column_int64(stmt, 4));
        rec.pathogenic_variants = static_cast<size_t>(sqlite3_column_int64(stmt, 5));
        rec.r_score = static_cast<float>(sqlite3_column_double(stmt, 6));
        rec.confluence = static_cast<float>(sqlite3_column_double(stmt, 7));
        rec.living_cells = static_cast<size_t>(sqlite3_column_int64(stmt, 8));
        rec.active_drug = reinterpret_cast<const char*>(sqlite3_column_text(stmt, 9));
        rec.tumor_regression_pct = static_cast<float>(sqlite3_column_double(stmt, 10));
        rec.resistant_cells = static_cast<size_t>(sqlite3_column_int64(stmt, 11));
        const unsigned char* notes = sqlite3_column_text(stmt, 12);
        if (notes) rec.notes = reinterpret_cast<const char*>(notes);

        list.push_back(std::move(rec));
    }
    sqlite3_finalize(stmt);
    return list;
}

std::vector<ReportRecord> SQLiteStorage::get_recent_reports(size_t limit) {
    std::vector<ReportRecord> list;
    if (!db_) return list;

    const char* sql = "SELECT id, timestamp, sample_id, tmb_score, r_score, biomarkers_count, "
                      "level1_count, report_md_path, report_json_path FROM clinical_reports "
                      "ORDER BY id DESC LIMIT ?;";

    sqlite3_stmt* stmt = nullptr;
    if (sqlite3_prepare_v2(db_, sql, -1, &stmt, nullptr) != SQLITE_OK) {
        return list;
    }

    sqlite3_bind_int64(stmt, 1, static_cast<sqlite3_int64>(limit));

    while (sqlite3_step(stmt) == SQLITE_ROW) {
        ReportRecord rec;
        rec.id = sqlite3_column_int64(stmt, 0);
        rec.timestamp = reinterpret_cast<const char*>(sqlite3_column_text(stmt, 1));
        rec.sample_id = reinterpret_cast<const char*>(sqlite3_column_text(stmt, 2));
        rec.tmb_score = static_cast<float>(sqlite3_column_double(stmt, 3));
        rec.r_score = static_cast<float>(sqlite3_column_double(stmt, 4));
        rec.biomarkers_count = static_cast<size_t>(sqlite3_column_int64(stmt, 5));
        rec.level1_count = static_cast<size_t>(sqlite3_column_int64(stmt, 6));
        rec.report_md_path = reinterpret_cast<const char*>(sqlite3_column_text(stmt, 7));
        rec.report_json_path = reinterpret_cast<const char*>(sqlite3_column_text(stmt, 8));

        list.push_back(std::move(rec));
    }
    sqlite3_finalize(stmt);
    return list;
}

std::vector<WorkflowRunRecord> SQLiteStorage::get_recent_workflows(size_t limit) {
    std::vector<WorkflowRunRecord> list;
    if (!db_) return list;

    const char* sql = "SELECT id, timestamp, workflow_id, workflow_name, total_steps, "
                      "completed_steps, status, summary_findings FROM workflow_runs "
                      "ORDER BY id DESC LIMIT ?;";

    sqlite3_stmt* stmt = nullptr;
    if (sqlite3_prepare_v2(db_, sql, -1, &stmt, nullptr) != SQLITE_OK) {
        return list;
    }

    sqlite3_bind_int64(stmt, 1, static_cast<sqlite3_int64>(limit));

    while (sqlite3_step(stmt) == SQLITE_ROW) {
        WorkflowRunRecord rec;
        rec.id = sqlite3_column_int64(stmt, 0);
        rec.timestamp = reinterpret_cast<const char*>(sqlite3_column_text(stmt, 1));
        rec.workflow_id = reinterpret_cast<const char*>(sqlite3_column_text(stmt, 2));
        rec.workflow_name = reinterpret_cast<const char*>(sqlite3_column_text(stmt, 3));
        rec.total_steps = static_cast<size_t>(sqlite3_column_int64(stmt, 4));
        rec.completed_steps = static_cast<size_t>(sqlite3_column_int64(stmt, 5));
        rec.status = reinterpret_cast<const char*>(sqlite3_column_text(stmt, 6));
        const unsigned char* sum = sqlite3_column_text(stmt, 7);
        if (sum) rec.summary_findings = reinterpret_cast<const char*>(sum);

        list.push_back(std::move(rec));
    }
    sqlite3_finalize(stmt);
    return list;
}

void SQLiteStorage::print_database_summary() const {
    std::cout << "================================================================================\n"
              << " ⌬ GENOMEFORGE :: SQLITE3 DATABASE REPOSITORY & AUDIT TRAIL\n"
              << " ⟦ SQLite v" << sqlite3_libversion() << " | Base de Datos: " << db_path_ << " ⟧\n"
              << "================================================================================\n";

    if (!db_) {
        std::cout << "[!] Base de datos no disponible o no abierta.\n";
        return;
    }

    auto count_table = [&](const char* table) -> int64_t {
        std::string q = "SELECT COUNT(*) FROM " + std::string(table) + ";";
        sqlite3_stmt* stmt = nullptr;
        int64_t count = 0;
        if (sqlite3_prepare_v2(db_, q.c_str(), -1, &stmt, nullptr) == SQLITE_OK) {
            if (sqlite3_step(stmt) == SQLITE_ROW) {
                count = sqlite3_column_int64(stmt, 0);
            }
            sqlite3_finalize(stmt);
        }
        return count;
    };

    int64_t n_sims = count_table("simulations");
    int64_t n_reports = count_table("clinical_reports");
    int64_t n_wfs = count_table("workflow_runs");

    std::uintmax_t file_size = 0;
    try {
        std::error_code ec;
        file_size = std::filesystem::file_size(db_path_, ec);
    } catch (...) {}

    std::cout << " Metricas del Repositorio:\n"
              << "  - Tamano archivo DB:       " << file_size << " bytes (" << (file_size / 1024) << " KB)\n"
              << "  - Corridas de Simulacion:  " << n_sims << "\n"
              << "  - Reportes MTB Emitidos:   " << n_reports << "\n"
              << "  - Protocolos de Workflow:  " << n_wfs << "\n\n";

    // Ultimos Workflows Registrados
    const_cast<SQLiteStorage*>(this)->get_recent_workflows(5);
    auto wfs = const_cast<SQLiteStorage*>(this)->get_recent_workflows(5);
    if (!wfs.empty()) {
        std::cout << " Ultimas Ejecuciones de Workflows Clinicos:\n";
        for (const auto& w : wfs) {
            std::cout << "  [ID " << w.id << "] " << w.timestamp << " | " << w.workflow_name 
                      << " (" << w.completed_steps << "/" << w.total_steps << " pasos) | "
                      << "Estado: " << w.status << "\n"
                      << "       -> " << w.summary_findings << "\n";
        }
        std::cout << "\n";
    }

    // Ultimos Reportes Clinicos
    auto reps = const_cast<SQLiteStorage*>(this)->get_recent_reports(5);
    if (!reps.empty()) {
        std::cout << " Ultimos Informes de Comite Molecular de Tumores (MTB):\n";
        for (const auto& r : reps) {
            std::cout << "  [ID " << r.id << "] " << r.timestamp << " | Muestra: " << r.sample_id
                      << " | TMB: " << r.tmb_score << " mut/Mb | R-Score: " << r.r_score
                      << " | Dianas FDA L1: " << r.level1_count << "\n"
                      << "       -> Markdown: " << r.report_md_path << "\n";
        }
        std::cout << "\n";
    }

    // Ultimas Simulaciones
    auto sims = const_cast<SQLiteStorage*>(this)->get_recent_simulations(5);
    if (!sims.empty()) {
        std::cout << " Ultimas Corridas de Simulacion Celular:\n";
        for (const auto& s : sims) {
            std::cout << "  [ID " << s.id << "] " << s.timestamp << " | Muestra: " << s.sample_name
                      << " | Celulas: " << s.living_cells << " | Confluencia: " 
                      << static_cast<int>(s.confluence * 100) << "% | Farmaco: " << s.active_drug
                      << " | Regresion: " << static_cast<int>(s.tumor_regression_pct) << "%\n";
        }
        std::cout << "\n";
    }

    std::cout << "================================================================================\n";
}

} // namespace genome::storage
