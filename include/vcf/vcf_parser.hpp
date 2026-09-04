#pragma once

#include <string>
#include <vector>
#include <cstdint>
#include <string_view>

namespace genome::vcf {

enum class VariantEffect : uint8_t {
    Benign = 0,
    LikelyBenign = 1,
    Uncertain = 2,
    LikelyPathogenic = 3,
    Pathogenic = 4,
    YamanakaCollision = 5
};

struct alignas(64) GenomicVariant {
    uint32_t pos{0};                // Coordenada genómica (1-indexed) (offset 0..3)
    uint16_t chrom_id{0};           // ID numérico de cromosoma (offset 4..5)
    VariantEffect effect{VariantEffect::Uncertain}; // (offset 6)
    uint8_t pad0{0};                // Alineación de 1 byte (offset 7)
    char ref[8]{0};                 // Alelo de referencia (offset 8..15)
    char alt[8]{0};                 // Alelo alternativo (offset 16..23)
    char gene_symbol[16]{0};        // Símbolo del gen (offset 24..39)
    float pathogenicity_score{0.0f};// Puntuación estimada (offset 40..43)
    float allele_frequency{0.0f};   // Frecuencia alélica observada (offset 44..47)
    uint32_t depth{0};              // Profundidad de cobertura (offset 48..51)
    uint32_t flags{0};              // Banderas auxiliares (offset 52..55)
    uint8_t pad1[8]{0};             // Relleno para completar exactamente 64 bytes (offset 56..63)
};
static_assert(sizeof(GenomicVariant) == 64, "GenomicVariant must match exact 64-byte cache line");

class VcfParser {
public:
    VcfParser();
    ~VcfParser();

    bool load_file(const std::string& filepath);

    const std::vector<GenomicVariant>& get_variants() const { return variants_; }
    size_t total_variants() const { return variants_.size(); }
    size_t pathogenic_count() const { return pathogenic_count_; }
    size_t yamanaka_collision_count() const { return yamanaka_collision_count_; }

    const std::string& get_filepath() const { return filepath_; }
    bool has_loaded() const { return loaded_; }

    // Simula corrección virtual de una variante (in-silico CRISPR edit)
    bool apply_virtual_correction(size_t index);

private:
    void parse_line(const std::string& line);
    uint16_t parse_chrom(std::string_view chrom_str);

    std::string filepath_;
    std::vector<GenomicVariant> variants_;
    size_t pathogenic_count_{0};
    size_t yamanaka_collision_count_{0};
    bool loaded_{false};
};

} // namespace genome::vcf
