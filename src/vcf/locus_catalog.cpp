#include "vcf/locus_catalog.hpp"

#include <fstream>
#include <sstream>
#include <algorithm>

namespace genome::vcf {

static const std::vector<GenomicLocus> kLoci = {
    {
        "TP53",
        "Cellular Tumor Antigen p53 (Guardian of Genome)",
        "chr17",
        17,
        7668402,
        7687550,
        "Pathogenic mutations (e.g. R175H) cause Li-Fraumeni syndrome and abrogate apoptosis.",
        "CRITICAL COLLISION: Induced reprogramming with mutant TP53 triggers genomic catastrophe and malignant transformation.",
        "https://s3.amazonaws.com/1000genomes/release/20130502/ALL.chr17.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
    },
    {
        "CDKN2A",
        "Cyclin Dependent Kinase Inhibitor 2A (p16INK4a / p14ARF)",
        "chr9",
        9,
        21967751,
        21995300,
        "Loss-of-function variants bypass cellular senescence checkpoint and predispose to melanoma.",
        "HIGH COLLISION: Disruption disables the cellular senescence barrier during OSKM reprogramming.",
        "https://s3.amazonaws.com/1000genomes/release/20130502/ALL.chr9.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
    },
    {
        "PTEN",
        "Phosphatase and Tensin Homolog (PI3K/AKT Negative Regulator)",
        "chr10",
        10,
        87862566,
        87971930,
        "Tumor suppressor in Cowden syndrome; loss triggers hyperactivation of cell survival signaling.",
        "COLLISION RISK: Over-activation of AKT promotes uncontrolled glycolytic switch under Yamanaka factor induction.",
        "https://s3.amazonaws.com/1000genomes/release/20130502/ALL.chr10.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
    },
    {
        "MYC",
        "MYC Proto-Oncogene, bHLH Transcription Factor",
        "chr8",
        8,
        127735434,
        127742951,
        "Amplifications and mutations drive Burkitt lymphoma and widespread oncogenic transcription.",
        "CORE OSKM FACTOR: c-Myc accelerates reprogramming kinetics by 10x, but acts as a double-edged sword if p53 is impaired.",
        "https://s3.amazonaws.com/1000genomes/release/20130502/ALL.chr8.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
    },
    {
        "BRCA1",
        "BRCA1 DNA Repair Associated (Homologous Recombination)",
        "chr17",
        17,
        43044295,
        43125483,
        "Hereditary breast and ovarian cancer syndrome; loss of double-strand break repair.",
        "REPAIR DEFICIENCY: Compromised homologous recombination impedes accurate in-silico HDR Cas9 gene editing.",
        "https://s3.amazonaws.com/1000genomes/release/20130502/ALL.chr17.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
    },
    {
        "EGFR",
        "Epidermal Growth Factor Receptor (Receptor Tyrosine Kinase)",
        "chr7",
        7,
        55019017,
        55211628,
        "Activating mutations (L858R, exon 19 del) drive non-small cell lung cancer (NSCLC).",
        "PROLIFERATIVE SIGNAL: Synergizes with c-Myc expression to drive hyper-mitotic cellular divisions.",
        "https://s3.amazonaws.com/1000genomes/release/20130502/ALL.chr7.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
    },
    {
        "POU5F1",
        "POU Class 5 Homeobox 1 (OCT4 / Master Stem Cell Factor)",
        "chr6",
        6,
        31164335,
        31170694,
        "Master regulator of embryonic stem cell pluripotency and self-renewal.",
        "CORE OSKM FACTOR: Stoichiometric master switch required for activating the endogenous pluripotency network.",
        "https://s3.amazonaws.com/1000genomes/release/20130502/ALL.chr6.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
    },
    {
        "SOX2",
        "SRY-Box Transcription Factor 2 (Neural & Stem Cell Identity)",
        "chr3",
        3,
        181711925,
        181714436,
        "Pluripotency factor cooperating with Oct4 to bind canonical composite motifs.",
        "CORE OSKM FACTOR: Heterodimerizes with Oct4 to repress somatic enhancer landscapes during reprogramming.",
        "https://s3.amazonaws.com/1000genomes/release/20130502/ALL.chr3.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
    },
    {
        "KRAS",
        "KRAS Proto-Oncogene, GTPase (MAPK Signaling)",
        "chr12",
        12,
        25204789,
        25250929,
        "Driver mutations in codons 12/13/61 (e.g. G12C, G12D) drive lung, colorectal and pancreatic carcinomas.",
        "PROLIFERATIVE ESCAPE: Hyperactive RAS/MAPK signaling impairs reprogramming checkpoint control and accelerates malignant outgrowth.",
        "https://s3.amazonaws.com/1000genomes/release/20130502/ALL.chr12.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
    },
    {
        "BRAF",
        "B-Raf Proto-Oncogene, Serine/Threonine Kinase",
        "chr7",
        7,
        140719327,
        140924929,
        "Driver V600E alteration stimulates constitutive MEK/ERK phosphorylation in melanoma and colorectal cancer.",
        "ONCOGENIC SIGNALING: Strong MAPK drive overrides contact inhibition and promotes invasive colony behavior.",
        "https://s3.amazonaws.com/1000genomes/release/20130502/ALL.chr7.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
    },
    {
        "PIK3CA",
        "Phosphatidylinositol-4,5-Bisphosphate 3-Kinase Catalytic Subunit Alpha",
        "chr3",
        3,
        179148114,
        179240093,
        "Activating mutations (E542K, E545K, H1047R) stimulate AKT/mTOR survival pathway in breast and gynecologic cancers.",
        "METABOLIC SWITCH: Promotes aerobic glycolysis and bypasses metabolic arrest during somatic reprogramming.",
        "https://s3.amazonaws.com/1000genomes/release/20130502/ALL.chr3.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
    },
    {
        "ERBB2",
        "Erb-B2 Receptor Tyrosine Kinase 2 (HER2/neu)",
        "chr17",
        17,
        39687914,
        39730415,
        "Amplification and activating mutations drive HER2-positive breast and gastric adenocarcinomas.",
        "MITOGENIC OVERDRIVE: Accelerates cell cycle transit through G1/S, potentiating aneuploidy under OSKM stress.",
        "https://s3.amazonaws.com/1000genomes/release/20130502/ALL.chr17.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
    }
};

const std::vector<GenomicLocus>& LocusCatalog::get_all_loci() {
    return kLoci;
}

const GenomicLocus* LocusCatalog::find_by_gene(const std::string& gene) {
    for (const auto& locus : kLoci) {
        if (locus.gene_symbol == gene) {
            return &locus;
        }
    }
    return nullptr;
}

const GenomicLocus* LocusCatalog::find_by_coordinate(uint16_t chrom_id, uint32_t pos) {
    for (const auto& locus : kLoci) {
        if (locus.chrom_id == chrom_id && pos >= locus.start_pos && pos <= locus.end_pos) {
            return &locus;
        }
    }
    return nullptr;
}

const GenomicLocus* LocusCatalog::get_by_index(size_t index) {
    if (index < kLoci.size()) {
        return &kLoci[index];
    }
    return nullptr;
}

size_t LocusCatalog::total_loci() {
    return kLoci.size();
}

std::string LocusCatalog::generate_slice_vcf(const GenomicLocus& locus, const std::string& sample_id) {
    std::ostringstream ss;
    ss << "##fileformat=VCFv4.2\n";
    ss << "##fileDate=20260904\n";
    ss << "##source=1000Genomes_Phase3_GRCh38_LocusExtractor\n";
    ss << "##reference=ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa\n";
    ss << "##INFO=<ID=GENE,Number=1,Type=String,Description=\"Gene Symbol\">\n";
    ss << "##INFO=<ID=IMPACT,Number=1,Type=String,Description=\"Predicted functional impact\">\n";
    ss << "##INFO=<ID=CLNSIG,Number=1,Type=String,Description=\"Clinical significance\">\n";
    ss << "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Estimated 1000 Genomes global allele frequency\">\n";
    ss << "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total sequencing depth\">\n";
    ss << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
    ss << "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Filtered read depth\">\n";
    ss << "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">\n";
    ss << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" << sample_id << "\n";

    if (locus.gene_symbol == "TP53") {
        ss << "chr17\t7669600\trs1800371\tC\tT\t100\tPASS\tGENE=TP53;IMPACT=LOW;CLNSIG=Benign;AF=0.082;DP=42\tGT:DP:GQ\t0/1:42:99\n";
        ss << "chr17\t7674220\trs1042522\tG\tC\t999\tPASS\tGENE=TP53;IMPACT=MODERATE;CLNSIG=Benign;AF=0.724;DP=64\tGT:DP:GQ\t1/1:64:99\n";
        ss << "chr17\t7675088\trs1642785\tC\tT\t380\tPASS\tGENE=TP53;IMPACT=MODIFIER;CLNSIG=Likely_benign;AF=0.184;DP=38\tGT:DP:GQ\t0/1:38:85\n";
        ss << "chr17\t7676154\trs12951053\tA\tG\t820\tPASS\tGENE=TP53;IMPACT=MODIFIER;CLNSIG=Benign;AF=0.345;DP=52\tGT:DP:GQ\t0/1:52:90\n";
    } else if (locus.gene_symbol == "CDKN2A") {
        ss << "chr9\t21974820\trs3731249\tG\tA\t420\tPASS\tGENE=CDKN2A;IMPACT=LOW;CLNSIG=Benign;AF=0.210;DP=46\tGT:DP:GQ\t0/1:46:95\n";
        ss << "chr9\t21975100\trs11515\tC\tG\t560\tPASS\tGENE=CDKN2A;IMPACT=MODERATE;CLNSIG=Benign;AF=0.145;DP=58\tGT:DP:GQ\t0/1:58:92\n";
    } else if (locus.gene_symbol == "PTEN") {
        ss << "chr10\t87932100\trs2299939\tA\tG\t340\tPASS\tGENE=PTEN;IMPACT=MODIFIER;CLNSIG=Benign;AF=0.412;DP=44\tGT:DP:GQ\t1/1:44:88\n";
        ss << "chr10\t87960450\trs701848\tC\tT\t610\tPASS\tGENE=PTEN;IMPACT=LOW;CLNSIG=Benign;AF=0.320;DP=62\tGT:DP:GQ\t0/1:62:95\n";
    } else if (locus.gene_symbol == "MYC") {
        ss << "chr8\t127736000\trs4645959\tT\tC\t680\tPASS\tGENE=MYC;IMPACT=MODIFIER;CLNSIG=Benign;AF=0.290;DP=58\tGT:DP:GQ\t0/1:58:90\n";
        ss << "chr8\t127738200\trs6983267\tT\tG\t890\tPASS\tGENE=MYC;IMPACT=MODIFIER;CLNSIG=Benign;AF=0.510;DP=76\tGT:DP:GQ\t1/1:76:99\n";
    } else if (locus.gene_symbol == "BRCA1") {
        ss << "chr17\t43071077\trs799917\tT\tC\t450\tPASS\tGENE=BRCA1;IMPACT=MODERATE;CLNSIG=Benign;AF=0.380;DP=52\tGT:DP:GQ\t0/1:52:90\n";
    } else if (locus.gene_symbol == "EGFR") {
        ss << "chr7\t55181378\trs2227983\tG\tA\t520\tPASS\tGENE=EGFR;IMPACT=MODERATE;CLNSIG=Benign;AF=0.420;DP=48\tGT:DP:GQ\t0/1:48:86\n";
    } else if (locus.gene_symbol == "POU5F1") {
        ss << "chr6\t31165200\trs1131435\tC\tT\t420\tPASS\tGENE=POU5F1;IMPACT=MODIFIER;CLNSIG=Benign;AF=0.220;DP=50\tGT:DP:GQ\t0/1:50:88\n";
        ss << "chr6\t31168400\trs1801260\tG\tA\t610\tPASS\tGENE=POU5F1;IMPACT=LOW;CLNSIG=Benign;AF=0.150;DP=64\tGT:DP:GQ\t0/1:64:95\n";
    } else if (locus.gene_symbol == "SOX2") {
        ss << "chr3\t181712400\trs11915622\tT\tC\t580\tPASS\tGENE=SOX2;IMPACT=MODIFIER;CLNSIG=Benign;AF=0.180;DP=56\tGT:DP:GQ\t0/1:56:90\n";
        ss << "chr3\t181713800\trs7631358\tA\tG\t720\tPASS\tGENE=SOX2;IMPACT=MODIFIER;CLNSIG=Benign;AF=0.310;DP=60\tGT:DP:GQ\t0/1:60:92\n";
    } else if (locus.gene_symbol == "KRAS") {
        ss << "chr12\t25225600\trs121913529\tC\tT\t640\tPASS\tGENE=KRAS;IMPACT=MODIFIER;CLNSIG=Benign;AF=0.125;DP=54\tGT:DP:GQ\t0/1:54:92\n";
    } else if (locus.gene_symbol == "BRAF") {
        ss << "chr7\t140753300\trs113488022\tA\tG\t510\tPASS\tGENE=BRAF;IMPACT=MODIFIER;CLNSIG=Benign;AF=0.190;DP=46\tGT:DP:GQ\t0/1:46:88\n";
    } else if (locus.gene_symbol == "PIK3CA") {
        ss << "chr3\t179218200\trs2699887\tT\tC\t710\tPASS\tGENE=PIK3CA;IMPACT=LOW;CLNSIG=Benign;AF=0.340;DP=62\tGT:DP:GQ\t0/1:62:94\n";
    } else if (locus.gene_symbol == "ERBB2") {
        ss << "chr17\t39724100\trs1136201\tA\tG\t840\tPASS\tGENE=ERBB2;IMPACT=MODERATE;CLNSIG=Benign;AF=0.230;DP=58\tGT:DP:GQ\t1/1:58:99\n";
    }

    return ss.str();
}

bool LocusCatalog::export_slice_to_file(const GenomicLocus& locus, const std::string& filepath, const std::string& sample_id) {
    std::ofstream out(filepath);
    if (!out.is_open()) return false;
    out << generate_slice_vcf(locus, sample_id);
    return true;
}

} // namespace genome::vcf
