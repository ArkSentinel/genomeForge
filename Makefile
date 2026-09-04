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

clean:
	rm -rf $(BUILD_DIR)

.PHONY: all build run clean
