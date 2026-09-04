#include "vcf/vcf_parser.hpp"

#include <fstream>
#include <sstream>
#include <iostream>
#include <cstring>
#include <algorithm>

namespace genome::vcf {

VcfParser::VcfParser() = default;
VcfParser::~VcfParser() = default;

uint16_t VcfParser::parse_chrom(std::string_view chrom_str) {
    if (chrom_str.starts_with("chr") || chrom_str.starts_with("CHR")) {
        chrom_str.remove_prefix(3);
    }
    if (chrom_str == "X" || chrom_str == "x") return 23;
    if (chrom_str == "Y" || chrom_str == "y") return 24;
    if (chrom_str == "M" || chrom_str == "MT" || chrom_str == "mt") return 25;
    
    try {
        int val = std::stoi(std::string(chrom_str));
        if (val >= 1 && val <= 22) return static_cast<uint16_t>(val);
    } catch (...) {}
    return 0;
}

void VcfParser::parse_line(const std::string& line) {
    if (line.empty() || line[0] == '#') return;

    std::stringstream ss(line);
    std::string chrom_s, pos_s, id_s, ref_s, alt_s, qual_s, filter_s, info_s;

    if (!(ss >> chrom_s >> pos_s >> id_s >> ref_s >> alt_s >> qual_s >> filter_s >> info_s)) {
        return;
    }

    GenomicVariant var{};
    var.chrom_id = parse_chrom(chrom_s);
    try {
        var.pos = static_cast<uint32_t>(std::stoul(pos_s));
    } catch (...) {
        var.pos = 0;
    }

    std::strncpy(var.ref, ref_s.c_str(), sizeof(var.ref) - 1);
    std::strncpy(var.alt, alt_s.c_str(), sizeof(var.alt) - 1);

    // Parsing de campos INFO (GENE, IMPACT, CLNSIG, AF, DP)
    std::string gene_found = "UNKNOWN";
    std::string clnsig_found = "";
    std::string impact_found = "";
    float af_val = 0.0f;
    uint32_t dp_val = 0;

    std::stringstream info_stream(info_s);
    std::string token;
    while (std::getline(info_stream, token, ';')) {
        auto eq_pos = token.find('=');
        if (eq_pos != std::string::npos) {
            std::string key = token.substr(0, eq_pos);
            std::string val = token.substr(eq_pos + 1);

            if (key == "GENE") {
                gene_found = val;
            } else if (key == "CLNSIG") {
                clnsig_found = val;
            } else if (key == "IMPACT") {
                impact_found = val;
            } else if (key == "AF") {
                try { af_val = std::stof(val); } catch (...) {}
            } else if (key == "DP") {
                try { dp_val = static_cast<uint32_t>(std::stoul(val)); } catch (...) {}
            }
        }
    }

    std::strncpy(var.gene_symbol, gene_found.c_str(), sizeof(var.gene_symbol) - 1);
    var.allele_frequency = af_val;
    var.depth = dp_val;

    // Clasificación y detección de colisión con Yamanaka
    bool is_pathogenic = (clnsig_found.find("Pathogenic") != std::string::npos || 
                          clnsig_found.find("pathogenic") != std::string::npos);
    bool is_critical_tumor_suppressor = (gene_found == "TP53" || gene_found == "CDKN2A" || 
                                         gene_found == "PTEN" || gene_found == "RB1");

    if (is_pathogenic && is_critical_tumor_suppressor) {
        var.effect = VariantEffect::YamanakaCollision;
        var.pathogenicity_score = 0.95f;
        yamanaka_collision_count_++;
        pathogenic_count_++;
    } else if (is_pathogenic) {
        var.effect = VariantEffect::Pathogenic;
        var.pathogenicity_score = 0.85f;
        pathogenic_count_++;
    } else if (clnsig_found.find("Likely_pathogenic") != std::string::npos) {
        var.effect = VariantEffect::LikelyPathogenic;
        var.pathogenicity_score = 0.65f;
        pathogenic_count_++;
    } else if (clnsig_found.find("Likely_benign") != std::string::npos ||
               clnsig_found.find("Benign") != std::string::npos) {
        var.effect = VariantEffect::Benign;
        var.pathogenicity_score = 0.05f;
    } else {
        var.effect = VariantEffect::Uncertain;
        var.pathogenicity_score = 0.30f;
    }

    variants_.push_back(var);
}

bool VcfParser::load_file(const std::string& filepath) {
    filepath_ = filepath;
    variants_.clear();
    pathogenic_count_ = 0;
    yamanaka_collision_count_ = 0;

    std::ifstream file(filepath);
    if (!file.is_open()) {
        loaded_ = false;
        return false;
    }

    std::string line;
    while (std::getline(file, line)) {
        parse_line(line);
    }

    loaded_ = true;
    return true;
}

bool VcfParser::apply_virtual_correction(size_t index) {
    if (index >= variants_.size()) return false;

    auto& v = variants_[index];
    if (v.effect == VariantEffect::YamanakaCollision) {
        if (yamanaka_collision_count_ > 0) yamanaka_collision_count_--;
        if (pathogenic_count_ > 0) pathogenic_count_--;
    } else if (v.effect == VariantEffect::Pathogenic || v.effect == VariantEffect::LikelyPathogenic) {
        if (pathogenic_count_ > 0) pathogenic_count_--;
    }

    // Efectuar corrección in-silico (alelo ALT reemplazado por REF)
    std::strncpy(v.alt, v.ref, sizeof(v.alt) - 1);
    v.effect = VariantEffect::Benign;
    v.pathogenicity_score = 0.02f;
    return true;
}

} // namespace genome::vcf
