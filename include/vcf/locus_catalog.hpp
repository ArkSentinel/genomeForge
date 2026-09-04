#pragma once

#include <string>
#include <vector>
#include <cstdint>
#include <optional>

namespace genome::vcf {

struct GenomicLocus {
    std::string gene_symbol;
    std::string full_name;
    std::string chromosome;
    uint16_t chrom_id{0};
    uint32_t start_pos{0};
    uint32_t end_pos{0};
    std::string clinical_significance;
    std::string yamanaka_role;
    std::string s3_remote_url;
};

class LocusCatalog {
public:
    // Retorna la lista de loci clínicos de referencia en GRCh38
    static const std::vector<GenomicLocus>& get_all_loci();

    // Búsqueda por símbolo genético
    static const GenomicLocus* find_by_gene(const std::string& gene);

    // Búsqueda por índice
    static const GenomicLocus* get_by_index(size_t index);
    static size_t total_loci();

    // Genera el contenido VCF v4.2 para un locus específico con anotaciones 1000 Genomas
    static std::string generate_slice_vcf(const GenomicLocus& locus, const std::string& sample_id = "HG00096");

    // Guarda el slice VCF a disco
    static bool export_slice_to_file(const GenomicLocus& locus, const std::string& filepath, const std::string& sample_id = "HG00096");
};

} // namespace genome::vcf
