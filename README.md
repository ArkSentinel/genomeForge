# GenomeForge-TUI 🧬

> **Herramienta bioinformática interactiva en terminal (TUI) desarrollada en C++20 con aceleración de hardware para Apple Silicon (M1–M4: CPU AMX/vDSP y GPU Metal) para simulación de redes de regulación génica (MoE), dinámica de reprogramación celular (Yamanaka OSKM), evaluación de inestabilidad genómica y diseño de intervenciones CRISPR-Cas9.**

---

## 📋 Tabla de Contenidos
* [Arquitectura del Sistema](#-arquitectura-del-sistema)
* [Características Principales](#-características-principales)
* [Requisitos del Sistema](#-requisitos-del-sistema)
* [Compilación y Ejecución](#-compilación-y-ejecución)
* [Navegación y Controles por Pestañas](#-navegación-y-controles-por-pestañas)
* [Motores Científicos y Matemáticos](#-motores-científicos-y-matemáticos)
* [Integración con el Proyecto 1000 Genomas](#-integración-con-el-proyecto-1000-genomas)
* [Documentación Adicional](#-documentación-adicional)

---

## 🏗 Arquitectura del Sistema

```
┌────────────────────────────────────────────────────────────────────────────────────────┐
│                                 GenomeForge-TUI (C++20)                                │
│          [1] DASHBOARD | [2] VCF | [3] OSKM | [4] MoE | [5] CRISPR | [6] SYSTEM        │
└───────────────────────────────────────────┬────────────────────────────────────────────┘
                                            │
       ┌────────────────────────────────────┼────────────────────────────────────┐
       ▼                                    ▼                                    ▼
┌──────────────┐                    ┌──────────────┐                     ┌──────────────┐
│ Yamanaka RK4 │                    │  Motor MoE   │                     │ CRISPR-Cas9  │
│  Continuous  │                    │  Accelerate  │                     │  ARM NEON    │
│  ODE Solver  │                    │  Metal (M4)  │                     │ Cas-OFFinder │
└──────┬───────┘                    └──────┬───────┘                     └──────┬───────┘
       │                                   │                                    │
       └───────────────────────────────────┼────────────────────────────────────┘
                                           ▼
       ┌────────────────────────────────────────────────────────────────────────┐
       │     Unified Memory Buffer (Zero-Copy MTLResourceStorageModeShared)     │
       │               Variantes VCF alineadas a 64 bytes (L1/L2)               │
       └────────────────────────────────────────────────────────────────────────┘
```

---

## ⚡ Características Principales

1. **Aceleración por Hardware en Apple Silicon:**
   * **CBLAS & vDSP:** Inferencia ultrarrápida del *Gating Router* en CPU AMX (~10–40 microsegundos).
   * **Metal Compute (GPU):** Pipelines de cómputo en Shading Language (MSL) ejecutándose en la GPU (Apple M4) con buffers en **Memoria Unificada Compartida (Zero-Copy sin transferencias por bus)**.
   * **ARM NEON SIMD:** Comparación de secuencias y cálculo de homología CRISPR vectorizado con instrucciones de 128 bits (`uint8x16_t`, `vceqq_u8`).
2. **Telemetría en Vivo del Kernel Darwin:**
   * Consulta directa de `mach/mach_host.h`, `host_processor_info` y `host_statistics64` para desglosar la carga en tiempo real entre núcleos de Rendimiento (**P-Cores**) y núcleos de Eficiencia (**E-Cores**), además del uso de memoria RAM unificada.
3. **Ingesta de Variantes VCF Optimizada para Caché:**
   * Estructura `GenomicVariant` con alineación estricta de 64 bytes (`alignas(64)`), garantizando correspondencia 1:1 con las líneas de caché de procesadores ARM64.
   * Clasificación de patogenicidad ClinVar y detección de colisiones de reprogramación en supresores tumorales (*TP53*, *CDKN2A*, *PTEN*).
4. **Simulador Continuo de Yamanaka (RK4):**
   * Resuelve numéricamente un sistema acoplado de 5 ecuaciones diferenciales ordinarias no lineales (ODEs) con el método Runge-Kutta de 4º Orden.
   * Cuantifica en tiempo real el **R-Score (Oncogenic Risk Score)** ante sobreexpresión de *c-Myc* en células con pérdida de función de *p53*.
5. **Motor CRISPR-Cas9 y Cas-OFFinder:**
   * Detección de motivos PAM `NGG`, extracción de protospacers de 20 nt y cálculo de eficacia con heurísticas Doench / Rule Set 2.
   * Cribado de *off-targets* mediante distancia de Hamming vectorizada en ARM NEON.
   * Simulación in-silico de reparación dirigida por homología (HDR), restaurando el genoma salvaje y recalculando la dinámica celular al instante.
6. **Interfaz TUI Reactiva en 6 Pestañas (FTXUI):**
   * Señalética científica sobria y de laboratorio (`◈`, `▶`, `❚❚`, `[!]`, `[x]`, `[+]`, `[-]`).
   * Navegación por pestañas independientes a pantalla completa (`[1]` a `[6]`).

---

## 💻 Requisitos del Sistema

* **Hardware:** Mac con procesador Apple Silicon (M1, M2, M3, M4 o variantes Pro/Max/Ultra).
* **Sistema Operativo:** macOS 14 (Sonoma) o macOS 15 (Sequoia) o superior.
* **Herramientas de Compilación:**
  * Xcode Command Line Tools (`clang++` compatible con C++20).
  * CMake 3.20+ (o superior).
  * Make (`/usr/bin/make`).
* **Frameworks de macOS Vinculados:**
  * `-framework Accelerate`
  * `-framework Metal`
  * `-framework Foundation`

---

## 🚀 Compilación y Ejecución

### 1. Compilación
Desde la raíz del repositorio:
```bash
make build
```

### 2. Ejecución
Para iniciar con el archivo del paciente sintético:
```bash
make run
```

Para iniciar con el genoma auténtico del **Proyecto 1000 Genomas**:
```bash
./build/genome_forge --vcf data/1000genomes_patient_tp53.vcf
```

Para consultar las opciones de línea de comandos:
```bash
./build/genome_forge --help
```

---

## 🎮 Navegación y Controles por Pestañas

| Atajo | Acción |
| :---: | :--- |
| **`1`** | Pestaña **`[1] DASHBOARD`**: Resumen clínico del paciente, monitor R-Score, estado del pulso y auditoría. |
| **`2`** | Pestaña **`[2] VCF EXPLORER`**: Tabla genómica a pantalla completa, frecuencias alélicas, cobertura DP y ventana de secuencia nucleotídica. |
| **`3`** | Pestaña **`[3] YAMANAKA (OSKM)`**: Ecuaciones diferenciales continuas ODE, niveles de factores OSKM y gráfica de sparkline temporal del R-Score. |
| **`4`** | Pestaña **`[4] MoE ROUTER`**: Arquitectura Mixture of Experts, pesos de la Gating Network, selección Top-2 y telemetría vDSP (CPU) vs Metal (GPU). |
| **`5`** | Pestaña **`[5] CRISPR STUDIO`**: Catálogo de ARNs guía (gRNA), análisis PAM NGG, scoring Cas-OFFinder y control de corte HDR. |
| **`6`** | Pestaña **`[6] SYSTEM & LOGS`**: Telemetría granular por núcleo de CPU (P-Cores vs E-Cores), páginas de memoria unificada y log de auditoría. |
| **`Tab`** | Ciclar ordenadamente a la siguiente pestaña. |
| **`↑` / `↓` o `k` / `j`** | Navegar por la tabla de variantes (pestaña 2) o por la lista de guías CRISPR (pestaña 5). |
| **`E`** | **Aplicar corrección CRISPR virtual (HDR repair):** Sustituye el codón mutado en *TP53*, reevalúa la dinámica celular y desploma el R-Score de 0.88 a 0.12 (nivel seguro). |
| **`R`** | **Reiniciar genoma:** Reinyecta la mutación patogénica para re-evaluar la alerta oncogénica. |
| **`Espacio`** | Pausar o reanudar el pulso dinámico de los factores de Yamanaka. |
| **`Q` o `Esc`** | Salir limpiamente de la aplicación. |

---

## 🔬 Motores Científicos y Matemáticos

### Sistema Dinámico de Reprogramación (Yamanaka ODEs)
Modelado a través de un sistema diferencial de 5 ecuaciones acopladas resueltas mediante **Runge-Kutta 4 (RK4)**:
$$\frac{dP}{dt} = k_{p1} [\text{Oct4}][\text{Sox2}] + k_{p2} [\text{Klf4}] - \gamma_p P$$
$$\frac{dS}{dt} = k_{s1} [\text{c-Myc}] \cdot p_{53} + k_{s2} [\text{Klf4}] - \gamma_s S$$
$$\frac{dI}{dt} = k_{i1} [\text{c-Myc}]^2 \cdot (1 - p_{53}) + k_{i2} \cdot \text{Burden} - \gamma_i I$$

### Métrica de Riesgo Oncogénico (R-Score)
$$R(t) = w_1 \cdot I(t) + w_2 \cdot \frac{[\text{c-Myc}](t)}{p_{53}(t) + \epsilon} + w_3 \cdot \text{Burden}$$
* $R < 0.25$: Estado seguro / reprogramación fisiológica viable.
* $0.25 \le R < 0.60$: Estrés replicativo moderado.
* $R \ge 0.60$: Alerta neoplásica crítica / pérdida de apoptosis.

---

## 🧬 Integración con el Proyecto 1000 Genomas

Para calibrar el sistema con datos de secuenciación humana real sin saturar el almacenamiento con terabytes de datos brutos, se incluye un extractor automatizado de *slices* indexados de alta fidelidad:

```bash
# Descargar/generar slice auténtico de TP53 (Chr 17) y CDKN2A (Chr 9)
./scripts/fetch_1000genomes_slice.sh
```

El dataset resultante en `data/1000genomes_patient_tp53.vcf` contiene polimorfismos reales con metadatos estándar de 1000 Genomas Phase 3 (GRCh38).

---

## 📚 Documentación Adicional

* [`docs/SPECS.md`](docs/SPECS.md): Especificación técnica exhaustiva, diseño de structs ARM64, ecuaciones matemáticas y arquitectura de software.
* [`docs/BITACORA.md`](docs/BITACORA.md): Bitácora viva de desarrollo con el registro cronológico de decisiones de ingeniería, hitos y benchmarks de hardware.
