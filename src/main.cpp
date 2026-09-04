#include "ui/tui_app.hpp"

#include <iostream>
#include <string>
#include <vector>

void print_usage(const char* prog_name) {
    std::cout << "GenomeForge-TUI: Herramienta de Modelado Genómico y TUI para Apple Silicon\n\n"
              << "Uso:\n"
              << "  " << prog_name << " [opciones]\n\n"
              << "Opciones:\n"
              << "  --vcf <archivo.vcf>   Ruta al archivo VCF de entrada (por defecto: data/sample_patient.vcf)\n"
              << "  --help, -h            Muestra esta ayuda\n\n"
              << "Controles en TUI:\n"
              << "  [↑/↓ o j/k]           Navegar por las variantes genómicas\n"
              << "  [E]                   Aplicar edición CRISPR virtual (in-silico HDR repair)\n"
              << "  [R]                   Reiniciar genoma al estado basal con mutaciones\n"
              << "  [Espacio]             Pausar / reanudar pulso dinámico de Yamanaka\n"
              << "  [Q o Esc]             Salir de la aplicación\n";
}

int main(int argc, char* argv[]) {
    std::string vcf_path = "data/sample_patient.vcf";

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            print_usage(argv[0]);
            return 0;
        } else if (arg == "--vcf" && i + 1 < argc) {
            vcf_path = argv[++i];
        }
    }

    try {
        genome::ui::GenomeForgeApp app(vcf_path);
        return app.run();
    } catch (const std::exception& e) {
        std::cerr << "Error crítico durante la ejecución de GenomeForge: " << e.what() << std::endl;
        return 1;
    }
}
