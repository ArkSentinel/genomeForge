cmake_minimum_required(VERSION 3.20)
project(GenomeForgeTUI VERSION 1.0.0 LANGUAGES CXX OBJCXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_OBJCXX_STANDARD 20)
set(CMAKE_OBJCXX_STANDARD_REQUIRED ON)

# Opciones de compilador y optimizaciones ARM64 para Apple Silicon
if(APPLE)
    add_compile_definitions(ACCELERATE_NEW_LAPACK)
    add_compile_options(
        -O3
        -Wall
        -Wextra
        -Wno-unused-parameter
    )
endif()

# Deshabilitar ejemplos y documentación de FTXUI
set(FTXUI_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(FTXUI_BUILD_DOCS OFF CACHE BOOL "" FORCE)

# Incorporar FTXUI
add_subdirectory(external/ftxui)

# Binario principal con todos los motores bioinformáticos y matemáticos
add_executable(genome_forge
    src/main.cpp
    src/system/telemetry.cpp
    src/vcf/vcf_parser.cpp
    src/vcf/locus_catalog.cpp
    src/yamanaka/yamanaka_sim.cpp
    src/moe/moe_engine.cpp
    src/metal/metal_compute.mm
    src/crispr/crispr_engine.cpp
    src/digital_twin/cell_twin.cpp
    src/export/csv_exporter.cpp
    src/ipc/state_inspector.cpp
    src/ui/tui_app.cpp
)

target_include_directories(genome_forge PRIVATE
    ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(genome_forge PRIVATE
    ftxui::screen
    ftxui::dom
    ftxui::component
    "-framework Accelerate"
    "-framework Metal"
    "-framework Foundation"
)
