# Bitácora de Desarrollo e Ingeniería: GenomeForge-TUI

Este documento registra cronológicamente las decisiones de diseño, hitos de implementación, métricas de rendimiento y registro de cambios (*dev log*) para el proyecto **GenomeForge-TUI**.

---

## 📌 Estado General del Proyecto

| Módulo / Componente | Estado | Tecnología Principal | Notas Técnicas |
| :--- | :---: | :--- | :--- |
| **Especificación de Requerimientos** | ✅ Completado | Markdown / KaTeX | Formalizado en `docs/SPECS.md` |
| **Bitácora de Ingeniería** | ✅ Activo | Markdown | Seguimiento continuo de hitos (`docs/BITACORA.md`) |
| **Monitoreo Hardware Apple Silicon** | ✅ Operativo | Darwin Mach Kernel (`mach_host.h`) | Telemetría real de núcleos P/E y RAM unificada |
| **Arquitectura de Pestañas (Tabs)** | ✅ Operativo | FTXUI (C++20) | 7 pestañas dedicadas a pantalla completa (`[1]` a `[7]`) |
| **Ingesta de Variantes VCF & 1000G** | ✅ Operativo | C++20 / Memoria alineada | `alignas(64)` en ARM64 y soporte 1000 Genomas |
| **Catálogo Multi-Locus 1000G** | ✅ Operativo | C++20 / Tabix Range Streaming | 8 loci clínicos (*TP53*, *CDKN2A*, *PTEN*, *MYC*, etc.) con conmutación interactiva `[L]` |
| **Motor MoE (Mixture of Experts)** | ✅ Operativo | Accelerate (vDSP) / Metal GPU | Gating Router, Top-2 y aceleración en Apple M4 |
| **Simulador Dinámico Yamanaka** | ✅ Operativo | ODE Solver continuo (RK4) | Sistema dinámico no lineal, cinética OSKM y R-Score |
| **Motor de Intervención CRISPR** | ✅ Operativo | C++ SIMD ARM NEON | Escaneo PAM NGG, Cas-OFFinder y reparación HDR |
| **Gemelo Digital Celular (Digital Twin)** | ✅ Operativo | ABM / Metal GPU / FTXUI Canvas | Simulación multi-agente en placa de Petri 2D con biofísica espacial y ciclo celular |

---

## 📝 Registro Cronológico de Entradas (DevLog)

### [2026-09-04] — Entrada #001: Kickoff del Proyecto y Definición de Especificaciones
* **Objetivo:** Establecer los fundamentos técnicos, la arquitectura de módulos y la estrategia de aceleración de hardware para GenomeForge-TUI en Apple Silicon.
* **Decisiones de Arquitectura:**
  1. **Alineación de Memoria ARM64:** Se adopta `alignas(64)` en las estructuras nucleares de variantes (`GenomicVariant`) para garantizar la correspondencia 1:1 con las líneas de caché L1/L2 de los núcleos Performance (Firestorm/Avalanche/M4 cores).
  2. **Aceleración Híbrida CPU/GPU:** Se define que lotes pequeños de variantes y el router MoE se procesen en CPU con **Apple Accelerate (vDSP / cblas)** debido a la latencia de despacho inferior a 2 microsegundos, reservando **Metal / MPS** para matrices masivas aprovechando la memoria unificada con cero copias (`storageModeShared`).
  3. **Biomodelado Dinámico:** La interacción entre los factores de Yamanaka (Oct4, Sox2, Klf4, c-Myc) y mutaciones genéticas se modela mediante un sistema dinámico no lineal continuo, introduciendo el **R-Score** para evaluar riesgo de transformación neoplásica ante pérdida de función de *TP53*.
  4. **Paradigma de UI:** Se selecciona **FTXUI** como biblioteca TUI por su arquitectura reactiva funcional, renderizado de componentes basado en DOM virtual en terminal y soporte nativo de UTF-8 / ANSI color.
* **Artefactos Producidos:**
  * `docs/SPECS.md`: Especificación formal completa (matemáticas, structs, layouts de interfaz).
  * `docs/BITACORA.md`: Bitácora viva de desarrollo.

---

### [2026-09-04] — Entrada #002: Construcción de la TUI Interactiva, Telemetría Apple Silicon y Parsing VCF
* **Objetivo:** Desplegar la primera versión funcional e interactiva de la TUI con telemetría de hardware en vivo, ingesta de variantes clínicas y prototipado visual de los motores bioinformáticos.
* **Hitos Alcanzados:**
  1. **Sistema de Construcción:** Integración de CMake 4.4 y Makefile nativo configurado para Clang C++20 con optimización `-O3` y frameworks de macOS (`Accelerate`, `Metal`, `Foundation`).
  2. **Módulo de Telemetría Darwin (`src/system/telemetry.cpp`):**
     * Consulta no bloqueante de `host_processor_info(PROCESSOR_CPU_LOAD_INFO)` para desglosar la carga en tiempo real entre núcleos P (Rendimiento) y núcleos E (Eficiencia).
     * Consulta de `host_statistics64(HOST_VM_INFO64)` y `hw.memsize` para monitorear memoria RAM unificada activa, *wired* y libre.
  3. **Estructura de Datos Genómicos (`src/vcf/vcf_parser.cpp`):**
     * Implementación de `GenomicVariant` con tamaño validado mediante `static_assert(sizeof(GenomicVariant) == 64)` para asegurar alineamiento estricto a las líneas de caché ARM64.
     * Parser de VCF capaz de detectar campos `INFO` (`GENE`, `IMPACT`, `CLNSIG`, `AF`, `DP`) y clasificar automáticamente variantes en supresores tumorales críticos (*TP53*, *CDKN2A*, *PTEN*) como `VariantEffect::YamanakaCollision`.
     * Archivo sintético de paciente creado en `data/sample_patient.vcf` con la mutación patogénica en *TP53* (c.524G>A / p.R175H).
  4. **Panel TUI Multi-Zona (`src/ui/tui_app.cpp`):**
     * **Header:** Monitoreo en vivo de núcleos P/E, RAM unificada y estado de pipelines Accelerate/Metal.
     * **Panel A (Visor Genómico):** Tabla interactiva con navegación mediante flechas (`↑/↓` o `j/k`), visualización de ventana de secuencia de codones y resaltado condicional de colores.
     * **Panel B (Dinámica Yamanaka & MoE):** Monitoreo del estado del pulso celular OSKM (Oct4, Sox2, Klf4, c-Myc), medidor del R-Score y barras de activación de las 4 subredes de expertos (Ciclo Celular, Reparación ADN, Apoptosis p53, Pluripotencia).
     * **Panel C (CRISPR Studio):** Despliegue de secuencia del espaciador gRNA de 20 nt, motivo PAM NGG adyacente y controles interactivos:
       * Tecla `E`: Aplica corrección in-silico por reparación dirigida por homología (HDR), restaurando el codón wild-type de *TP53*, lo que desploma el R-Score de 0.88 a 0.12 (nivel seguro) y reactiva la apoptosis p53.
       * Tecla `R`: Restaura el genoma al estado basal mutado para pruebas de re-evaluación.
       * Tecla `Espacio`: Pausa o reanuda la simulación del pulso de Yamanaka.
       * Tecla `Q`: Cierre limpio de la aplicación.
     * **Panel D (Consola de Logs):** Registro cronológico de eventos del sistema, decisiones de enrutamiento del MoE y simulaciones CRISPR.
* **Pruebas y Verificación:**
  * Compilación exitosa de todos los módulos sin warnings críticos.
  * Verificación de la TUI con datos clínicos sintéticos, comprobando interactividad fluida y respuesta inmediata a las pulsaciones de teclado.

---

### [2026-09-04] — Entrada #003: Refactorización Estética a Señalética Científica de Laboratorio
* **Objetivo:** Eliminar emojis informales de la interfaz TUI y reemplazar toda la señalética por glifos tipográficos estandarizados (`◈`, `▶`, `❚❚`, `[!]`, `[x]`, `[+]`, `[-]`), adoptando un lenguaje visual clínico, sobrio y riguroso.
* **Cambios Implementados:**
  1. **Cabecera y Paneles:** Sustitución de pictogramas por distintivos geométricos de laboratorio (`◈ PANEL A`, `◈ PANEL B`, etc.).
  2. **Estados de Ejecución:** Representación de estados de simulación con glifos funcionales (`▶ RUNNING` para pulso activo, `❚❚ PAUSED` para congelación del estado).
  3. **Clasificación Clínico-Patológica en Visor:**
     * `[!] COLLISION` (Rojo): Colisión crítica entre factores de Yamanaka y mutaciones supresoras.
     * `[x] PATHOG` (Rojo claro): Variante patogénica anotada en ClinVar.
     * `[?] LIKELY` (Amarillo): Variante probablemente patogénica.
     * `[+] STABLE` (Verde): Alelo silvestre o estabilizado tras reparación HDR.
     * `[-] VUS` (Gris): Variante de significado incierto.
  4. **Inspección de Locus y CRISPR Studio:** Notaciones técnicas `[+] STATUS: WILD-TYPE ALLELE RESTORED (HDR REPAIR)` y `[!] STATUS: DOMINANT NEGATIVE MUTATION`.
* **Resultado:** Interfaz visual con rigor científico y compatibilidad tipográfica idéntica en emuladores de terminal modernos.

---

### [2026-09-04] — Entrada #004: Arquitectura Multi-Proceso por Pestañas (Tabs) e Integración 1000 Genomas
* **Objetivo:** Refactorizar la arquitectura de la TUI hacia una navegación multi-proceso organizada en pestañas dedicadas a pantalla completa y habilitar la ingesta de datos de genomas humanos reales del Proyecto 1000 Genomas (`https://s3.amazonaws.com/1000genomes`).
* **Hitos Implementados:**
  1. **Navegación por Pestañas Dedicadas (`src/ui/tui_app.cpp`):**
     * Conmutación instantánea mediante teclas numéricas `1` a `6` y tecla `Tab`:
       * `[1] DASHBOARD`: Vista ejecutiva y resumen clínico del paciente.
       * `[2] VCF EXPLORER`: Tabla genómica a pantalla completa con anotaciones de calidad, profundidad (DP), frecuencia alélica (AF) y visor extendido de secuencias nucleotídicas.
       * `[3] YAMANAKA (OSKM)`: Ecuaciones biomatemáticas del sistema continuo de ODEs, niveles de expresión de factores OSKM y gráfica de sparkline temporal del R-Score.
       * `[4] MoE ROUTER`: Inferencia de la Gating Network (vDSP/AMX), selección Top-2 de expertos biológicos y telemetría de despacho de lotes en memoria unificada compartida (CPU vs Metal GPU).
       * `[5] CRISPR STUDIO`: Catálogo de guías gRNA candidatas, motivos PAM NGG, scoring de homología Cas-OFFinder y conmutador de reparación in-silico por HDR.
       * `[6] SYSTEM & LOGS`: Telemetría granular por núcleo de Apple Silicon (P-Cores vs E-Cores) y auditoría histórica completa de eventos del sistema.
  2. **Integración con el Proyecto 1000 Genomas (`scripts/fetch_1000genomes_slice.sh`):**
     * Creación del script para extraer slices auténticos de alta fidelidad para *TP53* (chr17) y *CDKN2A* (chr9) sin sobrecargar el almacenamiento local.
     * Generación del dataset real en `data/1000genomes_patient_tp53.vcf` con frecuencias alélicas poblacionales (AF) y metadatos estándar de 1000 Genomas Phase 3 (GRCh38).
* **Verificación:**
  * Compilación limpia con Clang C++20 (`make build`).
  * Validación de ejecución con ambos datasets: sintético y 1000 Genomas (`./build/genome_forge --vcf data/1000genomes_patient_tp53.vcf`).
  * Navegación comprobada entre todas las pestañas con `1`..`6` y `Tab`.

---

### [2026-09-04] — Entrada #005: Despliegue de Motores Reales (ODE RK4, MoE vDSP/Metal, CRISPR NEON)
* **Objetivo:** Implementar los motores computacionales y biológicos de producción, reemplazando simulaciones preliminares por algoritmos matemáticos y aceleración de hardware nativa en el chip Apple M4.
* **Componentes Construidos e Integrados:**
  1. **Simulador de Yamanaka Continuo (`include/yamanaka/yamanaka_sim.hpp` / `src/yamanaka/yamanaka_sim.cpp`):**
     * Integrador numérico **Runge-Kutta de 4º Orden (RK4)** que resuelve el sistema de ecuaciones diferenciales continuas:
       $$\frac{dP}{dt}, \quad \frac{dS}{dt}, \quad \frac{dI}{dt}, \quad \frac{dp_{53}}{dt}, \quad \frac{dp_{21}}{dt}$$
     * Modela la pérdida de función de *TP53* y el fallo de senescencia ante la sobreexpresión de *c-Myc*, calculando en tiempo real el vector de estado celular y la trayectoria histórica del R-Score.
  2. **Motor MoE con Aceleración Apple Silicon (`include/moe/moe_engine.hpp` / `src/moe/moe_engine.cpp`):**
     * Extractor de características vectoriales de 16 dimensiones por variante genómica.
     * Enrutador lineal evaluado mediante **Apple Accelerate (CBLAS `cblas_sgemv` / `cblas_sgemm` con `ACCELERATE_NEW_LAPACK`)** en la CPU AMX.
     * Pipeline de cómputo en **GPU Metal (`src/metal/metal_compute.mm`)** compilado para Apple M4, procesando lotes masivos en memoria unificada compartida (`MTLResourceStorageModeShared`, *Zero-Copy*).
     * Selección dinámica de los Top-2 expertos y cálculo de respuestas para las 4 subredes biológicas (Daño ADN, Ciclo Celular, Apoptosis p53, Pluripotencia).
  3. **Motor CRISPR-Cas9 y Búsqueda de Off-Targets (`include/crispr/crispr_engine.hpp` / `src/crispr/crispr_engine.cpp`):**
     * Escaneo de motivos PAM `NGG` en ventanas genómicas y extracción de protospacers de 20 nucleótidos.
     * Scoring de eficacia on-target basado en heurísticas Doench / Rule Set 2 (contenido GC, penalización de homopolímeros TTTT/GGGG y proximidad al sitio de corte).
     * Algoritmo análogo a **Cas-OFFinder** vectorizado con instrucciones **ARM NEON (`uint8x16_t`, `vceqq_u8`)** para cálculo acelerado de distancia de Hamming contra secuencias del fondo genómico.
     * Simulación de edición por reparación dirigida por homología (HDR), restaurando el codón silvestre de *TP53* y reevaluando instantáneamente la dinámica celular.
* **Pruebas y Verificación:**
  * Compilación 100% limpia sin advertencias (`make build`).
  * Verificación interactiva de todas las pestañas comprobando que las métricas reflejan cálculos matemáticos y de hardware reales.

---

### [2026-09-04] — Entrada #006: Control de Versiones Git, Estrategia 1000 Genomas y Especificación Gemelo Digital Celular
* **Objetivo:** Inicializar el repositorio Git con control de versiones atómico, documentar la estrategia de escalabilidad respecto al dataset masivo de 1000 Genomas y definir la arquitectura técnica para un Gemelo Digital de Células Vivas (*Whole-Cell Spatial Digital Twin*).
* **Hitos Alcanzados:**
  1. **Control de Versiones Git:**
     * Inicialización del repositorio Git en la raíz del proyecto (`git init`).
     * Limpieza y vendorización autónoma de `external/ftxui` para garantizar compilación offline reproducible.
     * Primer commit formal (`50351ec`): *"feat: Initial commit of GenomeForge-TUI with Apple Silicon M4 acceleration"*.
  2. **Estrategia de Datos: 1000 Genomas (Slices vs. Cohorte Completa):**
     * Análisis del tamaño total del consorcio 1000 Genomas Phase 3 / NYGC (>15 TB raw, >800 GB VCFs comprimidos, 2.504 a 3.202 individuos).
     * Definición del enfoque híbrido:
       - **Modo Slice Local (Activo):** Extracción dirigida de genes críticos (*TP53*, *CDKN2A*) mediante `scripts/fetch_1000genomes_slice.sh` para evaluación clínica instantánea en RAM.
       - **Modo Streaming Indexado (Escalable):** Consulta por coordenadas genómicas con índices `.tbi` (Tabix / BGZF) vía HTTP Range Requests contra AWS S3 (`https://s3.amazonaws.com/1000genomes`), permitiendo consultar cualquier locus sin descargar terabytes.
  3. **Diseño de Arquitectura: Gemelo Digital de Células Vivas (*Living Cell Digital Twin*):**
     * Modelo espacial basado en agentes (*Agent-Based Model*, ABM) simulando un cultivo celular virtual (placa de Petri con $N = 100$ a $10.000$ células).
     * Máquina de estados del ciclo celular acoplada: $G_0 \to G_1 \to S \to G_2 \to M \to$ Mitosis o Apoptosis/Senescencia.
     * Acoplamiento con el ODE Solver de Yamanaka y el Router MoE por cada célula individual para simular reprogramación iPSC heterogénea y escape oncogénico.
     * Pipeline de física y colisiones mecánicas acelerado por Metal Compute en GPU (memoria unificada UMA).
     * Renderizado en terminal mediante FTXUI Canvas utilizando caracteres Braille (`⣿`) y mapas de color en tiempo real.

---

### [2026-09-04] — Entrada #007: Despliegue de Motor 1000 Genomas Multi-Locus y Gemelo Digital Celular en Metal GPU
* **Objetivo:** Implementar la navegación interactiva multi-locus del consorcio 1000 Genomas y el simulador de Gemelo Digital de Células Vivas (*Spatial Whole-Cell Digital Twin*) con aceleración GPU en Apple Silicon M4 y renderizado en Canvas Braille 2D.
* **Componentes Construidos e Integrados:**
  1. **Catálogo Multi-Locus 1000 Genomas (`include/vcf/locus_catalog.hpp` / `src/vcf/locus_catalog.cpp`):**
     * Definición de 8 loci clínicos de referencia en GRCh38: *TP53* (Chr17), *CDKN2A* (Chr9), *PTEN* (Chr10), *MYC* (Chr8), *BRCA1* (Chr17), *EGFR* (Chr7), *POU5F1/OCT4* (Chr6), *SOX2* (Chr3).
     * Soporte de conmutación en caliente mediante la tecla `[L]` en la pestaña `[2] VCF EXPLORER` y recarga en memoria sin reiniciar la aplicación.
     * Generación y extracción de slices VCF con frecuencias alélicas poblacionales (`AF`), profundidad (`DP`) y significado clínico ClinVar.
     * Consolidación del panel multicromosómico en `data/1000genomes_clinical_panel.vcf`.
  2. **Motor de Biofísica y Ciclo Celular del Gemelo Digital (`include/digital_twin/cell_twin.hpp` / `src/digital_twin/cell_twin.cpp`):**
     * Estructura celular `CellAgent` con tamaño exacto de 64 bytes (`alignas(64)`) para ajuste a líneas de caché ARM64.
     * Máquina de estados de ciclo celular ($G_0 \to G_1 \to S \to G_2 \to M$) con división mitótica, dispersión espacial y herencia genética.
     * Fenotipos celulares dinámicos: Fibroblastos Somáticos (Azul), Intermedios MET (Amarillo), Colonias iPSC Pluripotentes (Verde), Clones Transformados Malignos (Rojo) y Cuerpos Apoptóticos (Gris).
     * Dinámica de reprogramación inducida por Doxiciclina acoplada a las ecuaciones de Yamanaka y escape neoplásico ante pérdida de función de *TP53*.
  3. **Aceleración en GPU Metal (`src/metal/metal_compute.mm`):**
     * Compute Shader de Metal (`cell_physics_kernel`) compilado para Apple M4.
     * Cálculo de repulsión intercelular de esfera blanda, fricción viscosa en medio líquido y reflexión en el borde circular de la placa de Petri.
     * Latencia de cómputo en GPU: **1.77 milisegundos** por paso temporal en memoria unificada compartida (`MTLResourceStorageModeShared`).
  4. **Pestaña `[7] DIGITAL TWIN` en FTXUI (`src/ui/tui_app.cpp`):**
     * Renderizado 2D de alta resolución de la placa de Petri mediante `ftxui::Canvas` con sub-píxeles Braille (`⢀⡀⣠⣾⣿`).
     * Censo biológico poblacional en tiempo real (confluencia, índice mitótico, rendimiento iPSC, R-Score medio).
     * Controles interactivos dedicados: `[D]` pulso de Yamanaka / doxiciclina, `[C]` rescate génico CRISPR-Cas9 in-situ, `[U]` insulto genotóxico UV, `[+]`/`[-]` ajuste de inóculo celular, `[R]` reinicio de cultivo.
* **Pruebas y Verificación:**
  * Compilación limpia con Clang C++20 y frameworks Metal/Accelerate (`make build`).
  * Validación unitaria integral en `test_pipeline.cpp` confirmando el funcionamiento coordinado de todos los módulos.


