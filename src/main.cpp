#include "ui/tui_app.hpp"
#include "ipc/state_inspector.hpp"

#include <iostream>
#include <string>
#include <vector>

void print_usage(const char* prog_name) {
    std::cout << "================================================================================\n"
              << " ⌬ GENOMEFORGE-TUI :: CLI & HUD ENGINE CONTROLLER\n"
              << " ⟦ Apple Silicon M-Series | UMA Zero-Copy | vDSP Accelerate | Metal Compute ⟧\n"
              << "================================================================================\n\n"
              << "Uso:\n"
              << "  " << prog_name << " [opciones]\n\n"
              << "Opciones Principales:\n"
              << "  --vcf <archivo.vcf>    Ruta al archivo VCF de entrada (default: data/sample_patient.vcf)\n"
              << "  --headless             Ejecuta la simulación en modo headless (sin interfaz TUI)\n"
              << "  --ticks <N>            Número de iteraciones a ejecutar en modo headless (default: 30)\n"
              << "  --snapshot <path>      Ruta del archivo JSON de telemetría (default: data/live_state.json)\n"
              << "  --inspect [path]       Inspecciona el estado de la TUI en vivo y muestra reporte HUD\n"
              << "  --export-csv           Fuerza la exportación de datasets CSV al finalizar\n"
              << "  --help, -h             Muestra esta ayuda detallada\n\n"
              << "Controles en TUI Interactiva:\n"
              << "  [1..7] o [Tab]         Navegar entre las 7 pestañas del Cyber-HUD\n"
              << "  [ [ / ] ]              Time Warp biológico: desacelerar o acelerar la simulación\n"
              << "  [X]                    Exportar datasets científicos a CSV (data/export_*.csv)\n"
              << "  [L]                    Conmutar locus genómico clínico (TP53, OCT4, SOX2, KLF4, MYC)\n"
              << "  [D]                    Inyectar/retirar pulso de Doxiciclina (OSKM inducido)\n"
              << "  [C]                    Terapia génica CRISPR-Cas9 in-situ (Rescate de TP53)\n"
              << "  [U]                    Insulto genotóxico por irradiación UV\n"
              << "  [+ / -]                Inocular / aspirar células del cultivo virtual (±25)\n"
              << "  [E]                    Simulación HDR en el explorador de variantes\n"
              << "  [R]                    Reiniciar cultivo / genoma al estado inicial\n"
              << "  [Espacio]              Pausar / reanudar dinámica celular y ODEs\n"
              << "  [Q / Esc]              Salir de la aplicación\n\n";
}

int main(int argc, char* argv[]) {
    std::string vcf_path = "data/sample_patient.vcf";
    bool headless = false;
    bool inspect = false;
    std::string inspect_path = "data/live_state.json";
    std::string snapshot_path = "data/live_state.json";
    size_t ticks = 30;
    bool export_csv = false;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            print_usage(argv[0]);
            return 0;
        } else if (arg == "--vcf" && i + 1 < argc) {
            vcf_path = argv[++i];
        } else if (arg == "--headless") {
            headless = true;
        } else if (arg == "--ticks" && i + 1 < argc) {
            ticks = std::stoull(argv[++i]);
        } else if (arg == "--snapshot" && i + 1 < argc) {
            snapshot_path = argv[++i];
        } else if (arg == "--inspect") {
            inspect = true;
            if (i + 1 < argc && argv[i + 1][0] != '-') {
                inspect_path = argv[++i];
            }
        } else if (arg == "--export-csv") {
            export_csv = true;
        }
    }

    // Modo de inspección de estado para AI Assistant / CLI
    if (inspect) {
        if (!genome::ipc::StateInspector::print_inspection_report(inspect_path)) {
            return 1;
        }
        return 0;
    }

    try {
        genome::ui::GenomeForgeApp app(vcf_path);

        // Modo Headless para tests, scripting y telemetría automatizada
        if (headless) {
            std::cout << "[GenomeForge] Iniciando ejecución HEADLESS (" << ticks << " ticks)...\n";
            int ret = app.run_headless(ticks, export_csv, snapshot_path);
            std::cout << "[GenomeForge] Ejecución HEADLESS completada exitosamente. Código: " << ret << "\n";
            return ret;
        }

        // Modo TUI Interactivo estándar
        return app.run();
    } catch (const std::exception& e) {
        std::cerr << "[CRITICAL ERROR] Fallo durante la ejecución de GenomeForge: " << e.what() << std::endl;
        return 1;
    }
}
