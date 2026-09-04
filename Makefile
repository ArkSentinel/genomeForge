BUILD_DIR = build
CMAKE = /opt/homebrew/bin/cmake

all: build

$(BUILD_DIR)/Makefile: CMakeLists.txt
	mkdir -p $(BUILD_DIR)
	cd $(BUILD_DIR) && $(CMAKE) .. -DCMAKE_BUILD_TYPE=Release

build: $(BUILD_DIR)/Makefile
	cd $(BUILD_DIR) && $(CMAKE) --build . -j

run: build
	./$(BUILD_DIR)/genome_forge --vcf data/sample_patient.vcf

run-1000g: build
	./$(BUILD_DIR)/genome_forge --vcf data/1000genomes_healthy_reference.vcf

fetch:
	./scripts/fetch_1000genomes_slice.sh

headless: build
	./$(BUILD_DIR)/genome_forge --vcf data/1000genomes_healthy_reference.vcf --headless --ticks 50 --export-csv

report: build
	./$(BUILD_DIR)/genome_forge --vcf data/sample_patient_oncology.vcf --report-onco data/molecular_tumor_board_report.md

test-onco: build
	./$(BUILD_DIR)/genome_forge --vcf data/sample_patient_oncology.vcf --headless --ticks 30 --export-csv --report-onco data/molecular_tumor_board_report.md
	./$(BUILD_DIR)/genome_forge --vcf data/sample_patient_resistance.vcf --headless --ticks 30 --export-csv --report-onco data/resistance_tumor_board_report.md

inspect: build
	./$(BUILD_DIR)/genome_forge --inspect

clean:
	rm -rf $(BUILD_DIR)

.PHONY: all build run run-1000g fetch headless report test-onco inspect clean
