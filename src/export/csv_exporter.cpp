#include "export/csv_exporter.hpp"

#include <fstream>
#include <sstream>
#include <iomanip>
#include <chrono>

namespace genome::export_tools {

std::string CsvExporter::get_timestamp_string() {
    auto now = std::chrono::system_clock::now();
    auto time_t_now = std::chrono::system_clock::to_time_t(now);
    std::tm tm_buf{};
    localtime_r(&time_t_now, &tm_buf);

    char buf[32];
    std::strftime(buf, sizeof(buf), "%Y%m%d_%H%M%S", &tm_buf);
    return std::string(buf);
}

bool CsvExporter::export_population_census(
    const std::string& filepath,
    const genome::twin::PopulationMetrics& metrics,
    size_t tick_step
) {
    bool file_exists = false;
    {
        std::ifstream check(filepath);
        file_exists = check.is_open();
    }

    std::ofstream out(filepath, std::ios::app);
    if (!out.is_open()) return false;

    if (!file_exists) {
        out << "Tick,SimulatedSeconds,SimulatedTime,TotalCells,SomaticFibroblasts,"
            << "IntermediatesMET,PluripotentiPSC,MalignantTransformed,ApoptoticCorpses,"
            << "ConfluencePercent,MitoticIndexPercent,iPSCEfficiencyPercent,AverageRScore,"
            << "DoxycyclineActive,CRISPRRescued,UVInsultActive\n";
    }

    out << tick_step << ","
        << std::fixed << std::setprecision(2) << metrics.simulated_seconds << ","
        << "\"" << metrics.simulated_time_str << "\","
        << metrics.total_cells << ","
        << metrics.somatic_count << ","
        << metrics.intermediate_count << ","
        << metrics.ipsc_count << ","
        << metrics.transformed_count << ","
        << metrics.apoptotic_count << ","
        << std::setprecision(2) << (metrics.confluence * 100.0f) << ","
        << std::setprecision(2) << metrics.mitotic_index << ","
        << std::setprecision(2) << metrics.ipsc_efficiency << ","
        << std::setprecision(4) << metrics.average_r_score << ","
        << (metrics.dox_active ? 1 : 0) << ","
        << (metrics.crispr_applied ? 1 : 0) << ","
        << (metrics.uv_insult_active ? 1 : 0) << "\n";

    return true;
}

bool CsvExporter::export_cell_agents(
    const std::string& filepath,
    const std::vector<genome::twin::CellAgent>& cells
) {
    std::ofstream out(filepath);
    if (!out.is_open()) return false;

    out << "CellID,ParentID,Generation,Phase,Phenotype,CleavageState,X,Y,Vx,Vy,Radius,"
        << "CycleClock,PluripotencyScore,SomaticScore,p53Level,RScore,HasTP53Mutation,IsTransformed\n";

    for (const auto& c : cells) {
        std::string phase_str = "G1";
        switch (c.phase) {
            case genome::twin::CellCyclePhase::G0_Quiescent: phase_str = "G0"; break;
            case genome::twin::CellCyclePhase::G1_Growth: phase_str = "G1"; break;
            case genome::twin::CellCyclePhase::S_Synthesis: phase_str = "S"; break;
            case genome::twin::CellCyclePhase::G2_PreMitotic: phase_str = "G2"; break;
            case genome::twin::CellCyclePhase::M_Mitosis: phase_str = "M"; break;
            case genome::twin::CellCyclePhase::Apoptotic: phase_str = "Apoptosis"; break;
            case genome::twin::CellCyclePhase::Senescent: phase_str = "Senescence"; break;
        }

        std::string pheno_str = "SomaticFibroblast";
        switch (c.phenotype) {
            case genome::twin::CellPhenotype::SomaticFibroblast: pheno_str = "SomaticFibroblast"; break;
            case genome::twin::CellPhenotype::IntermediateReprog: pheno_str = "IntermediateReprog"; break;
            case genome::twin::CellPhenotype::PluripotentIPSC: pheno_str = "PluripotentIPSC"; break;
            case genome::twin::CellPhenotype::TransformedMalignant: pheno_str = "TransformedMalignant"; break;
            case genome::twin::CellPhenotype::ApoptoticCorpse: pheno_str = "ApoptoticCorpse"; break;
        }

        out << c.cell_id << ","
            << c.parent_id << ","
            << c.generation << ","
            << phase_str << ","
            << pheno_str << ","
            << static_cast<int>(c.cleavage_state) << ","
            << std::fixed << std::setprecision(3) << c.x << ","
            << std::setprecision(3) << c.y << ","
            << std::setprecision(3) << c.vx << ","
            << std::setprecision(3) << c.vy << ","
            << std::setprecision(2) << c.radius << ","
            << std::setprecision(3) << c.cycle_clock << ","
            << std::setprecision(3) << c.P_score << ","
            << std::setprecision(3) << c.S_score << ","
            << std::setprecision(3) << c.p53_level << ","
            << std::setprecision(4) << c.r_score << ","
            << static_cast<int>(c.has_tp53_mutation) << ","
            << static_cast<int>(c.is_transformed) << "\n";
    }

    return true;
}

bool CsvExporter::export_genomic_variants(
    const std::string& filepath,
    const std::vector<genome::vcf::GenomicVariant>& variants
) {
    std::ofstream out(filepath);
    if (!out.is_open()) return false;

    out << "Chromosome,Position,GeneSymbol,RefAllele,AltAllele,Effect,PathogenicityScore,AlleleFrequency,ReadDepth\n";

    for (const auto& v : variants) {
        std::string chr_s = "chr" + std::to_string(v.chrom_id);
        if (v.chrom_id == 23) chr_s = "chrX";
        if (v.chrom_id == 24) chr_s = "chrY";

        std::string effect_str = "Uncertain";
        switch (v.effect) {
            case genome::vcf::VariantEffect::Benign: effect_str = "Benign"; break;
            case genome::vcf::VariantEffect::LikelyBenign: effect_str = "Likely_benign"; break;
            case genome::vcf::VariantEffect::Uncertain: effect_str = "VUS"; break;
            case genome::vcf::VariantEffect::LikelyPathogenic: effect_str = "Likely_pathogenic"; break;
            case genome::vcf::VariantEffect::Pathogenic: effect_str = "Pathogenic"; break;
            case genome::vcf::VariantEffect::YamanakaCollision: effect_str = "YamanakaCollision"; break;
        }

        out << chr_s << ","
            << v.pos << ","
            << v.gene_symbol << ","
            << v.ref << ","
            << v.alt << ","
            << effect_str << ","
            << std::fixed << std::setprecision(4) << v.pathogenicity_score << ","
            << std::setprecision(4) << v.allele_frequency << ","
            << v.depth << "\n";
    }

    return true;
}

std::string CsvExporter::export_all_datasets(
    const std::string& output_dir,
    const genome::twin::PopulationMetrics& metrics,
    const std::vector<genome::twin::CellAgent>& cells,
    const std::vector<genome::vcf::GenomicVariant>& variants,
    size_t tick_step
) {
    std::string ts = get_timestamp_string();
    std::string f_census = output_dir + "/export_cell_census_" + ts + ".csv";
    std::string f_agents = output_dir + "/export_cell_agents_" + ts + ".csv";
    std::string f_vars   = output_dir + "/export_genomic_variants_" + ts + ".csv";

    export_population_census(f_census, metrics, tick_step);
    export_cell_agents(f_agents, cells);
    export_genomic_variants(f_vars, variants);

    return "Datasets exported successfully: census (" + std::to_string(metrics.total_cells) + " cells), agents (" +
           std::to_string(cells.size()) + " records), variants (" + std::to_string(variants.size()) + " rows).";
}

} // namespace genome::export_tools
