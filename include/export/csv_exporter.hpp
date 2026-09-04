#pragma once

#include "digital_twin/cell_twin.hpp"
#include "vcf/vcf_parser.hpp"

#include <string>
#include <vector>

namespace genome::export_tools {

class CsvExporter {
public:
    // Exporta el censo biológico poblacional a un archivo CSV
    static bool export_population_census(
        const std::string& filepath,
        const genome::twin::PopulationMetrics& metrics,
        size_t tick_step
    );

    // Exporta el detalle microscópico de cada célula individual en la placa
    static bool export_cell_agents(
        const std::string& filepath,
        const std::vector<genome::twin::CellAgent>& cells
    );

    // Exporta las variantes genéticas analizadas con anotaciones patológicas
    static bool export_genomic_variants(
        const std::string& filepath,
        const std::vector<genome::vcf::GenomicVariant>& variants
    );

    // Exporta los 3 datasets de forma atómica a la carpeta data/ y retorna un resumen
    static std::string export_all_datasets(
        const std::string& output_dir,
        const genome::twin::PopulationMetrics& metrics,
        const std::vector<genome::twin::CellAgent>& cells,
        const std::vector<genome::vcf::GenomicVariant>& variants,
        size_t tick_step
    );

    static std::string get_timestamp_string();
};

} // namespace genome::export_tools
