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

headless: build
	./$(BUILD_DIR)/genome_forge --vcf data/1000genomes_clinical_panel.vcf --headless --ticks 50 --export-csv

inspect: build
	./$(BUILD_DIR)/genome_forge --inspect

clean:
	rm -rf $(BUILD_DIR)

.PHONY: all build run headless inspect clean
