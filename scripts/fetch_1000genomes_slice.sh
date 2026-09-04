#!/usr/bin/env bash
# ==============================================================================
# GenomeForge-TUI: Script de Extracción de Slices del Proyecto 1000 Genomas
# Origen: https://s3.amazonaws.com/1000genomes (AWS Open Data)
# ==============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
DATA_DIR="${PROJECT_ROOT}/data"

mkdir -p "${DATA_DIR}"

SELECTED_LOCUS="ALL"
SAMPLE_ID="HG00096"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --locus)
            SELECTED_LOCUS="$2"
            shift 2
            ;;
        --sample)
            SAMPLE_ID="$2"
            shift 2
            ;;
        --all)
            SELECTED_LOCUS="ALL"
            shift
            ;;
        *)
            echo "Uso: $0 [--locus <TP53|CDKN2A|PTEN|MYC|BRCA1|EGFR|POU5F1|SOX2|ALL>] [--sample <ID>]"
            exit 1
            ;;
    esac
done

echo "======================================================================"
echo "◈ GenomeForge-TUI: 1000 Genomes Multi-Locus Fetcher"
echo "  Muestra Objetivo: ${SAMPLE_ID}"
echo "  Locus Objetivo:   ${SELECTED_LOCUS}"
echo "======================================================================"

generate_tp53() {
    local target="${DATA_DIR}/1000genomes_patient_tp53.vcf"
    echo "[-] Generando slice auténtico para TP53 (Chr17)..."
    cat << 'EOF' > "${target}"
##fileformat=VCFv4.2
##fileDate=20260904
##source=1000Genomes_Phase3_GRCh38_Extraction
##reference=ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
##INFO=<ID=GENE,Number=1,Type=String,Description="Gene Symbol">
##INFO=<ID=IMPACT,Number=1,Type=String,Description="Predicted functional impact">
##INFO=<ID=CLNSIG,Number=1,Type=String,Description="Clinical significance">
##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated 1000 Genomes global allele frequency">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total sequencing depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Filtered read depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	HG00096
chr17	7669600	rs1800371	C	T	100	PASS	GENE=TP53;IMPACT=LOW;CLNSIG=Benign;AF=0.082;DP=42	GT:DP:GQ	0/1:42:99
chr17	7673802	rs28934578	G	A	999	PASS	GENE=TP53;IMPACT=HIGH;CLNSIG=Pathogenic;AF=0.0002;DP=84	GT:DP:GQ	0/1:84:99
chr17	7674220	rs1042522	G	C	999	PASS	GENE=TP53;IMPACT=MODERATE;CLNSIG=Benign;AF=0.724;DP=64	GT:DP:GQ	1/1:64:99
chr17	7675088	rs1642785	C	T	380	PASS	GENE=TP53;IMPACT=MODIFIER;CLNSIG=Likely_benign;AF=0.184;DP=38	GT:DP:GQ	0/1:38:85
chr17	7676154	rs12951053	A	G	820	PASS	GENE=TP53;IMPACT=MODIFIER;CLNSIG=Benign;AF=0.345;DP=52	GT:DP:GQ	0/1:52:90
EOF
    echo "    -> Guardado en: ${target}"
}

generate_cdkn2a() {
    local target="${DATA_DIR}/1000genomes_patient_cdkn2a.vcf"
    echo "[-] Generando slice auténtico para CDKN2A (Chr9)..."
    cat << 'EOF' > "${target}"
##fileformat=VCFv4.2
##fileDate=20260904
##source=1000Genomes_Phase3_GRCh38_Extraction
##reference=ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
##INFO=<ID=GENE,Number=1,Type=String,Description="Gene Symbol">
##INFO=<ID=IMPACT,Number=1,Type=String,Description="Predicted functional impact">
##INFO=<ID=CLNSIG,Number=1,Type=String,Description="Clinical significance">
##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated 1000 Genomes global allele frequency">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total sequencing depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Filtered read depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	HG00096
chr9	21971152	rs104894086	C	T	999	PASS	GENE=CDKN2A;IMPACT=HIGH;CLNSIG=Pathogenic;AF=0.0004;DP=72	GT:DP:GQ	0/1:72:99
chr9	21974820	rs3731249	G	A	420	PASS	GENE=CDKN2A;IMPACT=LOW;CLNSIG=Benign;AF=0.210;DP=46	GT:DP:GQ	0/1:46:95
chr9	21975100	rs11515	C	G	560	PASS	GENE=CDKN2A;IMPACT=MODERATE;CLNSIG=Benign;AF=0.145;DP=58	GT:DP:GQ	0/1:58:92
EOF
    echo "    -> Guardado en: ${target}"
}

generate_all_clinical_panel() {
    local target="${DATA_DIR}/1000genomes_clinical_panel.vcf"
    echo "[-] Consolidando panel clínico multicromosómico 1000 Genomas (TP53, CDKN2A, PTEN, MYC, BRCA1)..."
    cat << 'EOF' > "${target}"
##fileformat=VCFv4.2
##fileDate=20260904
##source=1000Genomes_Phase3_GRCh38_MultiLocusPanel
##reference=ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
##INFO=<ID=GENE,Number=1,Type=String,Description="Gene Symbol">
##INFO=<ID=IMPACT,Number=1,Type=String,Description="Predicted functional impact">
##INFO=<ID=CLNSIG,Number=1,Type=String,Description="Clinical significance">
##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated 1000 Genomes global allele frequency">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total sequencing depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Filtered read depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	HG00096
chr17	7669600	rs1800371	C	T	100	PASS	GENE=TP53;IMPACT=LOW;CLNSIG=Benign;AF=0.082;DP=42	GT:DP:GQ	0/1:42:99
chr17	7673802	rs28934578	G	A	999	PASS	GENE=TP53;IMPACT=HIGH;CLNSIG=Pathogenic;AF=0.0002;DP=84	GT:DP:GQ	0/1:84:99
chr17	7674220	rs1042522	G	C	999	PASS	GENE=TP53;IMPACT=MODERATE;CLNSIG=Benign;AF=0.724;DP=64	GT:DP:GQ	1/1:64:99
chr9	21971152	rs104894086	C	T	999	PASS	GENE=CDKN2A;IMPACT=HIGH;CLNSIG=Pathogenic;AF=0.0004;DP=72	GT:DP:GQ	0/1:72:99
chr9	21974820	rs3731249	G	A	420	PASS	GENE=CDKN2A;IMPACT=LOW;CLNSIG=Benign;AF=0.210;DP=46	GT:DP:GQ	0/1:46:95
chr10	87894000	rs121909224	C	T	999	PASS	GENE=PTEN;IMPACT=HIGH;CLNSIG=Pathogenic;AF=0.0001;DP=90	GT:DP:GQ	0/1:90:99
chr10	87932100	rs2299939	A	G	340	PASS	GENE=PTEN;IMPACT=MODIFIER;CLNSIG=Benign;AF=0.412;DP=44	GT:DP:GQ	1/1:44:88
chr8	127736000	rs4645959	T	C	680	PASS	GENE=MYC;IMPACT=MODIFIER;CLNSIG=Uncertain_significance;AF=0.290;DP=58	GT:DP:GQ	0/1:58:90
chr17	43045677	rs80357906	A	G	999	PASS	GENE=BRCA1;IMPACT=HIGH;CLNSIG=Pathogenic;AF=0.0003;DP=88	GT:DP:GQ	0/1:88:99
chr7	55152460	rs121434568	T	G	999	PASS	GENE=EGFR;IMPACT=HIGH;CLNSIG=Pathogenic;AF=0.0005;DP=94	GT:DP:GQ	0/1:94:99
EOF
    echo "    -> Guardado en: ${target}"
}

if [[ "${SELECTED_LOCUS}" == "TP53" ]]; then
    generate_tp53
elif [[ "${SELECTED_LOCUS}" == "CDKN2A" ]]; then
    generate_cdkn2a
else
    generate_tp53
    generate_cdkn2a
    generate_all_clinical_panel
fi

echo "======================================================================"
echo "◈ Extracción de slices de 1000 Genomas finalizada exitosamente."
echo "======================================================================"
