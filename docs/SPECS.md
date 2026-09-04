# GenomeForge-TUI: Especificación Técnica y Diseño Funcional

**Versión:** 1.0.0  
**Fecha:** 2026-09-04  
**Plataforma Objetivo:** Apple Silicon (macOS Sonoma 14+ / Sequoia 15+, ARM64 M1–M4)  
**Lenguaje Estándar:** C++20 / C++23 (Apple Clang)  

---

## 1. Visión General del Sistema

**GenomeForge-TUI** es una herramienta interactiva de alto rendimiento en terminal diseñada para:
1. **Ingestar y mapear variantes genéticas locales (VCF)** en estructuras de memoria optimizadas para la jerarquía de caché de procesadores Apple Silicon (ARM64).
2. **Simular redes de regulación génica (GRN)** mediante una arquitectura de **Mezcla de Expertos (MoE)**, ejecutando el enrutamiento y la multiplicación matricial en CPU (vDSP/AMX) y GPU (Metal) sin copias redundantes gracias a la Memoria Unificada (*Unified Memory Architecture*).
3. **Modelar la dinámica de reprogramación celular por factores de Yamanaka (OSKM: Oct4, Sox2, Klf4, c-Myc)** evaluando colisiones con mutaciones patogénicas (ej. *TP53*, *CDKN2A*, *PTEN*) y cuantificando un **R-Score (Oncogenic Risk Score)**.
4. **Diseñar intervenciones CRISPR-Cas9 automatizadas** con búsqueda de motivos PAM (`NGG`), cálculo de homología de corte y filtrado de *off-targets* vectorizado con ARM NEON.
5. **Ofrecer una interfaz de usuario reactiva en terminal (TUI con FTXUI)** que muestra telemetría del chip en tiempo real, visor genómico interactivo y logs de inferencia biológica.

```
┌────────────────────────────────────────────────────────────────────────┐
│                        GenomeForge-TUI Core                            │
└───────────────────────────────────┬────────────────────────────────────┘
                                    │
       ┌────────────────────────────┼───────────────────────────┐
       ▼                            ▼                           ▼
┌──────────────┐             ┌──────────────┐            ┌──────────────┐
│  FTXUI TUI   │             │ Inferencia   │            │  Motores     │
│  Multi-Panel │             │ Accelerate   │            │  Biológicos  │
│  Live Render │             │ Metal (GPU)  │            │  OSKM/CRISPR │
└──────┬───────┘             └──────┬───────┘            └──────┬───────┘
       │                            │                           │
       └────────────────────────────┼───────────────────────────┘
                                    ▼
       ┌────────────────────────────────────────────────────────┐
       │     Unified Memory Buffer (Zero-Copy ARM64 Cache)      │
       │       Variantes VCF alineadas a 64 bytes (L1/L2)       │
       └────────────────────────────────────────────────────────┘
```

---

## 2. Requisitos de Hardware y Plataforma

* **Arquitectura:** ARM64 (Apple Silicon M1, M2, M3, M4).
* **Jerarquía de Caché:** Optimización para líneas de caché de 64 bytes (o 128 bytes en clústeres P-Core). Alineación forzada con `alignas(64)` para evitar *false sharing* y *cache line splitting*.
* **Memoria Unificada (UMA):** Uso de buffers compartidos (`MTL::ResourceStorageModeShared`) permitiendo que la CPU cargue y preprocese datos genómicos y la GPU ejecute kernels de inferencia sin copias de memoria intermedias.
* **Telemetría de Hardware:** Lectura directa del kernel Darwin vía `mach/mach_host.h`, `processor_info` y `sysctlbyname` para monitorear:
  * Carga de núcleos de rendimiento (P-Cores) y eficiencia (E-Cores).
  * Uso de memoria RAM compartida (activa, inactiva, *wired* y libre).
  * Presión de memoria del sistema.

---

## 3. Módulos del Sistema

### Módulo A: Ingesta de VCF y Mapeo en Memoria ARM64

#### Responsabilidad
Parsear archivos Variant Call Format (VCF v4.2+) de forma no bloqueante, extrayendo variantes de un solo nucleótido (SNVs) e inserciones/deleciones (InDels) en genes diana clave para senescencia, oncogénesis y reprogramación.

#### Estructuras de Datos
```cpp
namespace genome::vcf {

enum class VariantEffect : uint8_t {
    Benign = 0,
    LikelyBenign = 1,
    Uncertain = 2,
    LikelyPathogenic = 3,
    Pathogenic = 4,
    YamanakaCollision = 5
};

struct alignas(64) GenomicVariant {
    uint32_t pos;                   // Posición genómica (1-indexed) (offset 0..3)
    uint16_t chrom_id;              // ID del cromosoma (1..24 para 1-22, X, Y) (offset 4..5)
    VariantEffect effect;           // Clasificación funcional (offset 6)
    uint8_t pad0;                   // Alineación de 1 byte (offset 7)
    char ref[8];                    // Alelo de referencia (offset 8..15)
    char alt[8];                    // Alelo alternativo (offset 16..23)
    char gene_symbol[16];           // Símbolo del gen (ej. "TP53", "CDKN2A") (offset 24..39)
    float pathogenicity_score;      // Puntuación de patogenicidad [0.0, 1.0] (offset 40..43)
    float allele_frequency;         // Frecuencia alélica observada (AF) (offset 44..47)
    uint32_t flags;                 // Banderas de filtrado / flags de calidad (offset 48..51)
    uint32_t depth;                 // Profundidad de lectura DP (offset 52..55)
    uint8_t pad1[8];                // Relleno explícito a 64 bytes (offset 56..63)
};
static_assert(sizeof(GenomicVariant) == 64, "GenomicVariant must fit in a 64-byte cache line");

} // namespace genome::vcf
```

---

### Módulo B: Motor MoE (Mixture of Experts) con Aceleración Hardware

#### Arquitectura del MoE
El motor divide el procesamiento de variantes en $E$ redes de expertos independientes especializadas por vía biológica:
1. **Experto 1: Ciclo Celular y Senescencia** (control p16INK4a, p21, RB1).
2. **Experto 2: Respuesta a Daño de ADN / Apoptosis** (vía TP53, BAX, PUMA).
3. **Experto 3: Estabilidad Genómica y Telomérica** (TERT, ATM, ATR, BRCA1).
4. **Experto 4: Pluripotencia y Red de Transcripción** (Oct4, Nanog, Sox2 regulon).

```
               [ Vector de Características de la Variante (x ∈ R^D) ]
                                      │
                                      ▼
                        ┌───────────────────────────┐
                        │   Gating Router Network   │
                        │     W_g · x + b_g         │
                        │    (Accelerate / vDSP)    │
                        └─────────────┬─────────────┘
                                      │
                          g = Softmax(TopK(W_g · x))
                                      │
               ┌──────────────────────┼──────────────────────┐
               ▼                      ▼                      ▼
         [ Experto 1 ]          [ Experto 2 ]          [ Experto 3 ]
         Ciclo Celular          Daño ADN/p53           Pluripotencia
               │                      │                      │
               └──────────────────────┼──────────────────────┘
                                      ▼
                      [ Respuesta Ponderada: y = Σ g_i E_i(x) ]
```

#### Enrutamiento Matemático
Para un vector de características de variante $\mathbf{x} \in \mathbb{R}^D$:
$$\mathbf{h} = \mathbf{W}_g \mathbf{x} + \mathbf{b}_g$$
$$g_i = \frac{\exp(h_i)}{\sum_{j=1}^E \exp(h_j)}$$
En la implementación de CPU, la multiplicación $\mathbf{W}_g \mathbf{x}$ se ejecuta mediante `cblas_sgemv` o `vDSP_mmul` del framework `Accelerate`. Cuando el lote de variantes excede $N > 256$, la ejecución se despacha al kernel de Metal en la GPU.

---

### Módulo C: Simulación de Yamanaka (OSKM) y Cálculo del R-Score

#### Formulación Biomatemática (Sistema Dinámico de ODEs)
La reprogramación celular mediada por los factores de Yamanaka se modela a través de un sistema de ecuaciones diferenciales no lineales acopladas:

$$\frac{d[\text{Pluri}]}{dt} = k_1 [\text{Oct4}][\text{Sox2}] + k_2 [\text{Klf4}] - \gamma_p [\text{Pluri}]$$

$$\frac{d[\text{Senescencia}]}{dt} = k_3 [\text{c-Myc}] \cdot (1 - \Omega_{\text{TP53}}) + k_4 [\text{Klf4}] - \gamma_s [\text{Senescencia}]$$

$$\frac{d[\text{Inestabilidad}]}{dt} = k_{\text{myc}} [\text{c-Myc}]^2 \cdot \Omega_{\text{TP53}} + k_{\text{dna}} (1 - \Omega_{\text{Reparación}}) - \gamma_i [\text{Inestabilidad}]$$

Donde:
* $\Omega_{\text{TP53}} \in [0, 1]$ representa la pérdida de función de p53 (0 = normal/salvaje, 1 = mutación patogénica completa como R175H).
* Si p53 está mutado ($\Omega_{\text{TP53}} \approx 1.0$), el pulso de *c-Myc* no logra activar la senescencia protectora o la apoptosis, provocando una amplificación exponencial del término de Inestabilidad.

#### Métrica de Riesgo Oncogénico (R-Score)
$$R(t) = w_1 \cdot [\text{Inestabilidad}](t) + w_2 \cdot \frac{[\text{c-Myc}](t)}{[\text{p53}_{\text{funcional}}](t) + \epsilon} + w_3 \cdot \sum_{v \in \text{VCF}} P(v)$$

* **$R < 0.25$:** Estado Seguro / Reprogramación Fisiológica (Verde).
* **$0.25 \le R < 0.60$:** Advertencia / Estrés Replicativo Elevado (Amarillo).
* **$R \ge 0.60$:** Peligro Crítico / Riesgo de Transformación Neoplásica (Rojo).

---

### Módulo D: Motor CRISPR-Cas9 y Diseño de Intervenciones

#### Mecanismo de Búsqueda y Edición
1. **Identificación del Locus Mutado:** Coordenadas $chr:pos$ de la mutación diana.
2. **Escaneo PAM:** Detección de motivos `NGG` en la hebra *forward* (ventana $[-25, +25]$ bp) y `CCN` en la hebra *reverse complement*.
3. **Selección de Protospacer:** Extracción de la secuencia de 20 nucleótidos adyacente a cada PAM detectado:
   $$5'\text{-}[N_{20}]\text{-NGG-}3'$$
4. **Scoring de Eficiencia y Especificidad:**
   * Distancia al sitio de corte (corte a 3 pb corriente arriba del PAM).
   * Contenido GC óptimo ($40\% - 60\%$).
   * Penalización de *homopolímeros* ($T_4, G_4$).
5. **Filtrado de Off-Targets:** Cálculo de distancias de Hamming con secuencias homólogas genómicas usando comparaciones vectorizadas NEON (`veorq_u8`, `vcntq_u8`).
6. **Edición Virtual In-Silico:** El motor permite al usuario simular la corrección del alelo (`ALT -> REF`), recalculando instantáneamente el R-Score para certificar la neutralización del riesgo.

---

### Módulo E: Interfaz de Usuario de Terminal (FTXUI)

#### Arquitectura Visual Multi-Panel
```
┌────────────────────────────────────────────────────────────────────────┐
│ GenomeForge-TUI v1.0 [Apple Silicon M-Series | Metal GPU & vDSP AMX]   │
│ CPU: [P-Cores: 18%  E-Cores: 8%]  Unified RAM: 4.2 / 16.0 GB  UMA: OK  │
├───────────────────────────────────┬────────────────────────────────────┤
│ PANEL A: VISOR GENÓMICO EN VIVO   │ PANEL B: MOE & YAMANAKA DYNAMICS   │
│ Chromosome 17 - TP53 (17p13.1)    │ Reprogramming Status: HIGH RISK    │
│                                   │                                    │
│ 5'-CTGCTCCCCCCCTGCCCTCAACAAGA-3'  │ R-Score: 0.782 [██████████░░] ALERTA│
│              ▲                    │                                    │
│   Pos: 7,673,802  G > A [R175H]   │ Expert Activation:                 │
│   Effect: PATHOGENIC MUTATION     │  [1] DNA Repair:    [███████░░░] 72%│
│                                   │  [2] Cell Cycle:    [█████░░░░░] 54%│
│ Status: Yamanaka Collision Alert! │  [3] p53 Apoptosis: [█░░░░░░░░░] 12%│
├───────────────────────────────────┴────────────────────────────────────┤
│ PANEL C: CRISPR-CAS9 REPAIR STUDIO & INTERVENTION PROPOSALS            │
│ Target PAM: TGG (+12 bp) | Spacer: 5'-GAGGTTGGCTCTGACTGTAC-3' (20 nt)  │
│ Cutting Efficiency: 89.4% | Off-Target Risk: MINIMAL (Score: 98.2/100) │
│ Action: [Press 'E' to Apply In-Silico CRISPR Edit]                     │
├────────────────────────────────────────────────────────────────────────┤
│ PANEL D: LOGS & INFERENCE TELEMETRY                                    │
│ [INFO] VCF Loaded: 14,280 variants indexed in Unified Memory.          │
│ [WARN] Yamanaka simulation pulse triggered c-Myc explosion at t=42.    │
│ [CRIT] Unchecked proliferation detected due to p53 loss-of-function.   │
└────────────────────────────────────────────────────────────────────────┘
```

#### Arquitectura de Navegación por Pestañas (Process Tabs)
Para maximizar el detalle analítico y permitir la inspección profunda de cada subsistema sin saturar la pantalla, la interfaz se estructura en 7 pestañas independientes a pantalla completa:

* **Pestaña `[1] DASHBOARD`:** Resumen ejecutivo de la muestra del paciente, estado global de la simulación, nivel de alerta del R-Score y acceso rápido a corrección.
* **Pestaña `[2] VCF EXPLORER`:** Tabla genómica interactiva exhaustiva con catálogo de loci 1000 Genomas, selector en caliente (`[L]`), profundidad de lectura (DP), frecuencias alélicas poblacionales (AF) y visualizador de secuencia flanqueante.
* **Pestaña `[3] YAMANAKA (OSKM)`:** Simulador biomatemático de ecuaciones diferenciales ordinarias (ODEs continuas), cinética de factores OSKM (Oct4, Sox2, Klf4, c-Myc), eje de control downstream (p53 / p21) y gráfica temporal del R-Score.
* **Pestaña `[4] MoE ROUTER`:** Motor de Mezcla de Expertos, vector de pesos de la Gating Network, distribución Softmax sobre las 4 subredes biológicas (Daño ADN, Ciclo Celular, Apoptosis, Pluripotencia) y telemetría de latencia en CPU (Accelerate vDSP/AMX) vs GPU (Metal Shared UMA).
* **Pestaña `[5] CRISPR STUDIO`:** Catálogo de ARNs guía (gRNA) candidatos para motivos PAM NGG, scoring de especificidad, estimación de corte de doble hebra (DSB) y conmutador de reparación in-silico por HDR.
* **Pestaña `[6] SYSTEM & LOGS`:** Auditoría completa de eventos del sistema y telemetría granular de hardware Apple Silicon por núcleo individual (P-Cores vs E-Cores) y páginas de memoria unificada.
* **Pestaña `[7] DIGITAL TWIN`:** Gemelo digital multicelular espacial (*Agent-Based Model*) con placa de Petri virtual en 2D, física de colisiones acelerada por GPU Metal en memoria unificada, máquina de estados de ciclo celular ($G_0 \to M$), colonias iPSC (verdes), clones transformados por mutación en *TP53* (rojos) y renderizado sub-pixel con caracteres Braille.

#### Atajos de Teclado Globales
* `1` a `7`: Salto directo instantáneo a cualquier pestaña de proceso.
* `Tab`: Ciclar ordenadamente entre las 7 pestañas.
* `L`: Conmutar dinámicamente el gen/locus activo del consorcio 1000 Genomas (*TP53*, *CDKN2A*, *PTEN*, *MYC*, *BRCA1*, *EGFR*, *OCT4*, *SOX2*).
* `D`: Inyectar o suspender Doxiciclina en el medio de cultivo (inducción de factores de Yamanaka OSKM).
* `C`: Administrar terapia génica dirigida CRISPR-Cas9 in-situ sobre las células del Gemelo Digital.
* `U`: Disparar estrés genotóxico por radiación UV para poner a prueba los puntos de control del ciclo celular.
* `+` / `-`: Inocular o aspirar células en el cultivo de la placa de Petri (+25 / -25 células).
* `↑` / `↓` o `k` / `j`: Navegar filas en tablas genómicas y catálogos de guías.
* `Espacio`: Pausar o reanudar el pulso dinámico de las simulaciones.
* `E`: Ejecutar intervención virtual CRISPR-Cas9 (HDR repair) sobre la variante VCF seleccionada.
* `R`: Restablecer el cultivo celular (si se está en Tab 7) o restablecer el genoma original del paciente.
* `Q` o `Esc`: Cierre limpio de la aplicación.

---

## 4. Gemelo Digital Celular Espacial (*Whole-Cell Spatial Digital Twin*)

### Biofísica Intercelular y Dinámica de Población
* **Modelo Basado en Agentes (ABM):** Población de $N = 100$ a $2.000$ células en una placa de Petri 2D de radio $R = 45.0$ mm.
* **Estructura `CellAgent` (64 bytes, `alignas(64)`):** Coordenadas $(x,y)$, velocidades $(v_x,v_y)$, radio celular, masa inercial, reloj de ciclo celular, niveles de pluripotencia ($P$), identidad somática ($S$), supresión tumoral ($p_{53}$), riesgo oncogénico ($R$-Score) y mutaciones genómicas.
* **Máquina de Estados de Ciclo Celular:**
  $$G_0 \text{ (Quiescencia)} \longrightarrow G_1 \text{ (Crecimiento)} \longrightarrow S \text{ (Síntesis ADN)} \longrightarrow G_2 \longrightarrow M \text{ (Mitosis)}$$
* **Eventos Fenotípicos:**
  * **Fibroblasto Somático (Azul):** Célula diferenciada quiescente o en ciclo normal.
  * **Intermedio MET (Amarillo):** Célula en transición mesénquima-epitelio durante reprogramación.
  * **Colonia Pluripotente iPSC (Verde):** Agrupaciones compactas de células madre con alta adhesión intercelular.
  * **Escape Oncogénico / Neoplásico (Rojo):** Células con mutación patogénica en *TP53* sometidas a inducción con *c-Myc*, proliferando a velocidad $3\times$.
  * **Apoptosis (Gris):** Muerte celular programada inducida por $p_{53}$ funcional ante estrés genotóxico UV.
* **Aceleración Metal GPU:** Kernel `cell_physics_kernel` en Objective-C++ y MSL, resolviendo colisiones elásticas de esfera blanda a 60 FPS en memoria unificada compartida (`storageModeShared`).

---

## 5. Estándar de Validación con Datos Reales (Proyecto 1000 Genomas)

El sistema soporta la ingesta y benchmarking frente a conjuntos de variantes auténticas provenientes del **1000 Genomes Project** (`https://s3.amazonaws.com/1000genomes`):
* **Genoma de Referencia:** GRCh38 / hg38.
* **Muestras Clave:** Muestra HG00096 / NA12878.
* **Catálogo Multi-Locus:**
  * `TP53` (Chr17: 7,668,402 - 7,687,550 bp) — Guardián del genoma y colisión con Yamanaka.
  * `CDKN2A` (Chr9: 21,967,751 - 21,995,300 bp) — Barrera de senescencia celular.
  * `PTEN` (Chr10: 87,862,566 - 87,971,930 bp) — Supresor tumoral vía PI3K/AKT.
  * `MYC` (Chr8: 127,735,434 - 127,742,951 bp) — Proto-oncogén y factor de reprogramación OSKM.
  * `BRCA1` (Chr17: 43,044,295 - 43,125,483 bp) — Reparación por recombinación homóloga.
  * `EGFR` (Chr7: 55,019,017 - 55,211,628 bp) — Receptor tirosina quinasa de proliferación.
  * `POU5F1` / `OCT4` (Chr6: 31,164,335 - 31,170,694 bp) — Factor maestro de pluripotencia.
  * `SOX2` (Chr3: 181,711,925 - 181,714,436 bp) — Factor de transcripción stem cell.
* **Script de Automatización:** `scripts/fetch_1000genomes_slice.sh`.

---

## 6. Pipeline de Construcción y Rendimiento

* **Compilador:** `clang++` con `-std=c++20 -O3` (optimizado para Apple Silicon ARM64).
* **Frameworks de macOS:**
  * `-framework Accelerate` (vDSP, CBLAS, LAPACK).
  * `-framework Metal` (cómputo GPU y memoria unificada Zero-Copy UMA).
  * `-framework Foundation` (servicios del sistema).
* **Dependencia TUI:** FTXUI integrada de forma autónoma (`external/ftxui`).
