# GenomeForge-TUI

> Herramienta bioinformática interactiva en terminal (TUI) desarrollada en C++20 con aceleración de hardware para Apple Silicon (M1-M4: AMX/vDSP en CPU y Metal Shading Language en GPU). Integra simulación de redes de regulación génica (MoE), dinámica de reprogramación celular continua (factores de Yamanaka OSKM), evaluación de riesgo oncogénico (R-Score), diseño de intervenciones genómicas CRISPR-Cas9, Gemelo Digital Celular espacial (Whole-Cell Spatial Digital Twin) y una plataforma clínica de Oncología de Precisión y Farmacogenómica con generación de reportes de Comité Molecular de Tumores (MTB).

---

## Tabla de Contenidos
* [1. Visión General del Sistema](#1-visión-general-del-sistema)
* [2. Arquitectura de Hardware Apple Silicon](#2-arquitectura-de-hardware-apple-silicon)
* [3. Navegación y Procesos del Cyber-HUD (8 Pestañas)](#3-navegación-y-procesos-del-cyber-hud-8-pestañas)
* [4. Fundamentos Matemáticos y Modelado Biológico](#4-fundamentos-matemáticos-y-modelado-biológico)
  * [Dinámica No Lineal de Yamanaka (ODE Solver RK4)](#dinámica-no-lineal-de-yamanaka-ode-solver-rk4)
  * [Métrica de Riesgo Oncogénico (R-Score)](#métrica-de-riesgo-oncogénico-r-score)
  * [Red de Enrutamiento MoE y Subredes de Expertos](#red-de-enrutamiento-moe-y-subredes-de-expertos)
  * [Motor de Intervención CRISPR-Cas9 y Homología HDR](#motor-de-intervención-crispr-cas9-y-homología-hdr)
* [5. Gemelo Digital Celular Espacial (Living Cell Digital Twin)](#5-gemelo-digital-celular-espacial-living-cell-digital-twin)
  * [Modelo Basado en Agentes (ABM)](#modelo-basado-en-agentes-abm)
  * [Controlador Temporal Multi-Escala (Time Warp)](#controlador-temporal-multi-escala-time-warp)
  * [Visualización Citomorfológica en Canvas Braille](#visualización-citomorfológica-en-canvas-braille)
* [6. Motor de Oncología de Precisión y Farmacogenómica](#6-motor-de-oncología-de-precisión-y-farmacogenómica)
  * [Niveles de Evidencia Clínica y Base Curada](#niveles-de-evidencia-clínica-y-base-curada)
  * [Farmacodinamia In-Silico en el Gemelo Digital](#farmacodinamia-in-silico-en-el-gemelo-digital)
  * [Reportes Clínicos de Comité Molecular de Tumores (MTB)](#reportes-clínicos-de-comité-molecular-de-tumores-mtb)
* [7. Interoperabilidad: Modo Headless, Inspector IA, Datasets CSV y Reportes](#7-interoperabilidad-modo-headless-inspector-ia-datasets-csv-y-reportes)
* [8. Integración Multi-Locus con el Proyecto 1000 Genomas (12 Loci Clínicos)](#8-integración-multi-locus-con-el-proyecto-1000-genomas-12-loci-clínicos)
* [9. Requisitos, Compilación y Referencia CLI](#9-requisitos-compilación-y-referencia-cli)
* [10. Tabla de Atajos de Teclado](#10-tabla-de-atajos-de-teclado)
* [11. Documentación y Referencias Técnicas](#11-documentación-y-referencias-técnicas)

---

## 1. Visión General del Sistema

GenomeForge-TUI es una estación de trabajo computacional para medicina de precisión y biología sintética ejecutada íntegramente en la terminal. Su propósito es simular el impacto biológico de la reprogramación celular mediante factores de Yamanaka (*Oct4*, *Sox2*, *Klf4*, *c-Myc* / OSKM) sobre genomas portadores de variantes patogénicas en supresores tumorales y proto-oncogenes críticos (*TP53*, *CDKN2A*, *PTEN*, *KRAS*, *BRAF*, *PIK3CA*, *ERBB2*, *EGFR*, *BRCA1*), cuantificar el riesgo de escape neoplásico, proyectar la respuesta terapéutica a fármacos oncológicos dirigidos (TKIs, inhibidores PARP, reactivadores de p53, anticuerpos conjugados e inmunoterapias) y permitir el rescate genómico in-silico mediante edición CRISPR-Cas9 por recombinación homóloga (HDR).

El sistema fusiona bioinformática estructural, oncología traslacional, física de partículas intercelulares, cinéticas bioquímicas de ecuaciones diferenciales continuas y computación de alto rendimiento aprovechando la arquitectura de memoria unificada (UMA) de los chips de Apple.

---

## 2. Arquitectura de Hardware Apple Silicon

GenomeForge-TUI está diseñado específicamente para exprimir las capacidades arquitectónicas del silicio de Apple:

```text
+-------------------------------------------------------------------------------------------------------+
|                                      GenomeForge-TUI Core (C++20)                                     |
|           < 1: DASHBOARD > | < 2: VCF > | < 3: OSKM > | < 4: MoE > | < 5: CRISPR >                    |
|           < 6: SYSTEM > | < 7: DIGITAL TWIN > | < 8: ONCO THERAPY > | < 9: WORKFLOWS >                 |
+---------------------------------------------------┬---------------------------------------------------+
                                                    │
             +--------------------------------------+--------------------------------------+
             |                                      |                                      |
             v                                      v                                      v
    +------------------+                  +------------------+                   +------------------+
    |   Yamanaka RK4   |                  |   Router MoE     |                   |  Biofísica 2D    |
    |   Continuo ODE   |                  |  Accelerate vDSP |                   |    Metal GPU     |
    |   Sistema 5x5    |                  |  Apple Matrix AMX|                   |  Kernel MSL M4   |
    +--------+---------+                  +--------+---------+                   +--------+---------+
             |                                      |                                      |
             +--------------------------------------+--------------------------------------+
                                                    │
                                                    v
+-------------------------------------------------------------------------------------------------------+
|                             Unified Memory Architecture (Zero-Copy Buffer)                            |
|                    GenomicVariant (alignas(64))  |  CellAgent (alignas(64))                           |
|                                  MTLResourceStorageModeShared                                         |
+-------------------------------------------------------------------------------------------------------+
```

1. **Alineación a Líneas de Caché ARM64 (`alignas(64)`):**
   * Tanto las variantes genómicas (`GenomicVariant`) como los agentes celulares (`CellAgent`) están estrictamente alineados y rellenados a 64 bytes (`static_assert(sizeof(...) == 64)`).
   * Esto elimina el *cache line splitting* y el *false sharing*, permitiendo que los núcleos Performance (Firestorm/Avalanche/M4 cores) carguen los datos mediante instrucciones vectorizadas de un solo ciclo (`LDP`/`STP`).
2. **Memoria Unificada Cero-Copias (Zero-Copy UMA):**
   * Los buffers de partículas celulares se asignan utilizando `MTLResourceStorageModeShared`. La CPU prepara y actualiza los parámetros biológicos en memoria y la GPU Metal procesa los kernels de colisión mecánica sobre el mismo espacio de direccionamiento físico sin transferencias a través de un bus PCIe.
3. **Aceleración Vectorial con Apple Accelerate:**
   * La red de compuerta (Gating Network) del modelo Mixture of Experts (MoE) ejecuta la multiplicación matriz-vector a través de `cblas_sgemv` y transformaciones vectoriales `vDSP`, logrando latencias de inferencia de 1.8 microsegundos en la unidad matricial Apple AMX.
4. **Vectorización SIMD ARM NEON (128-bit):**
   * El cálculo de especificidad de guías frente a homólogos genómicos (análogo a Cas-OFFinder) está vectorizado con registros `uint8x16_t` e intrínsecos `vceqq_u8`, evaluando 16 nucleótidos por instrucción de ciclo.
5. **Telemetría no bloqueante del Kernel Mach (`mach_host.h`):**
   * Muestreo continuo sin impacto en el rendimiento de núcleos P (Performance) vs E (Efficiency), paginación y memoria *wired*.

---

## 3. Navegación y Procesos del Cyber-HUD (9 Pestañas + Soporte Total de Ratón)

La interfaz gráfica en terminal (construida con FTXUI en C++20) sigue una estética de Cyber-HUD científico, organizada en 9 pestañas a pantalla completa con soporte integral de clics de ratón, rueda de desplazamiento y selector modal de muestras:

* **< 1: DASHBOARD >:**
  * Resumen clínico del paciente, diagnóstico genómico global, indicador holográfico del R-Score poblacional, confluencia del cultivo, estado de los inductores biológicos (Doxiciclina / CRISPR) y consola de eventos.
* **< 2: VCF EXPLORER >:**
  * Tabla exhaustiva de variantes del paciente y del consorcio 1000 Genomas. Coordenadas cromosómicas, alelos REF/ALT, profundidad de secuenciación (DP), frecuencia alélica poblacional (AF) y visualizador de la ventana flanqueante (+/- 20 pares de bases) con codones diana. Soporta selección de variantes por clic de ratón y conmutación en caliente de locus mediante la tecla `[L]`.
* **< 3: YAMANAKA DYNAMICS >:**
  * Monitor de cinéticas de reprogramación continua mediante ecuaciones diferenciales ordinarias (ODEs). Expresión de factores OSKM (Oct4, Sox2, Klf4, c-Myc), barrera de senescencia, inestabilidad genómica acumulada, puntos de control de ciclo p53/p21 y gráfica sparkline histórica del R-Score.
* **< 4: MoE ROUTER >:**
  * Arquitectura Mixture of Experts. Extracción de embedding genómico de 16 dimensiones, distribución de pesos Softmax en la Gating Network y activación proporcional sobre las 4 subredes biológicas. Telemetría comparada de latencia AMX/vDSP frente a Metal GPU.
* **< 5: CRISPR STUDIO >:**
  * Entorno de edición genómica dirigida. Escaneo de motivos PAM NGG en hebras forward y reverse, extracción de protospacers de 20 nt, scoring de corte on-target (heurística Rule Set 2 / Doench) y distancia de corte a la mutación. Selección interactiva de candidatos gRNA por clic y control `[E]` para ejecutar reparación in-silico por homología (HDR).
* **< 6: SYSTEM & LOGS >:**
  * Diagnóstico en vivo de la plataforma Apple Silicon. Monitoreo individualizado por núcleo de CPU (P-Cores vs E-Cores), deslose de memoria unificada (activa, cableada, libre) y registro de auditoría del sistema.
* **< 7: DIGITAL TWIN >:**
  * Placa de Petri virtual en 2D interactiva. Renderizado citomorfológico de alta fidelidad en sub-píxeles Braille, censo celular en vivo, dinámica clonal, controlador de velocidad temporal multi-escala (*Time Warp*), panel de farmacodinamia activa y botones de acción clicables (Doxiciclina, CRISPR, UV, Sembrar, Aspirar, Zoom in/out/reset).
* **< 8: ONCO THERAPY >:**
  * Plataforma de Oncología de Precisión y Farmacogenómica. Matching automatizado de alteraciones patogénicas frente a fármacos dirigidos aprobados por la FDA y guías NCCN. Monitoreo de carga mutacional tumoral (TMB) para elegibilidad de inmunoterapia anti-PD-1, detección de mutaciones de resistencia secundaria (*EGFR* T790M contra Gefitinib), botones interactivos de inoculación rápida por clic (`[T]` TKI, `[P]` PARP, `[K]` p53, `[W]` Lavar) y exportador directo del reporte de Comité Molecular de Tumores (`[O]`).
* **< 9: WORKFLOWS (ESTUDIO CLÍNICO VISUAL) >:**
  * Pipeline gráfico guiado sin código para personal médico y clínico. Presenta 5 protocolos estandarizados con diagrama de flujo visual (cajas con estado `[DONE]`, `[ACTIVO]`, `[ESPERA]`, `[FALLIDO]` unidas por flechas de flujo), controles de ejecución con clic o atajo (`[Enter]` Flujo Completo, `[S]` Paso a Paso, `[R]` Reiniciar, `[F]` Cambiar Muestra) y tarjeta de hallazgos clínicos acumulados en tiempo real.
* **Selector Modal de Muestras VCF (`[F]`):**
  * Ventana modal interactiva sobrepuesta en cualquier momento de la simulación. Escanea el almacenamiento local en busca de archivos `.vcf`, permitiendo al operador seleccionar y cargar perfiles de pacientes en caliente con un clic o con las flechas de navegación y Enter.

---

## 4. Fundamentos Matemáticos y Modelado Biológico

### Dinámica No Lineal de Yamanaka (ODE Solver RK4)

El proceso de reprogramación celular hacia pluripotencia inducida (iPSC) y su conflicto con la integridad genómica se formula como un sistema continuo de 5 ecuaciones diferenciales ordinarias acopladas:

$$\frac{dP}{dt} = k_{p1} \cdot [\text{Oct4}][\text{Sox2}] + k_{p2} \cdot [\text{Klf4}] - \gamma_p \cdot P$$

$$\frac{dS}{dt} = k_{s1} \cdot [\text{c-Myc}] \cdot (1 - \Omega_{\text{TP53}}) + k_{s2} \cdot [\text{Klf4}] - \gamma_s \cdot S$$

$$\frac{dI}{dt} = k_{\text{myc}} \cdot [\text{c-Myc}]^2 \cdot (1 - p_{53}) + k_{\text{dna}} \cdot \text{Burden} - \gamma_i \cdot I$$

$$\frac{dp_{53}}{dt} = k_{\text{act}} \cdot I \cdot \Omega_{\text{TP53}} - \gamma_{\text{deg}} \cdot p_{53}$$

$$\frac{dp_{21}}{dt} = k_{21} \cdot p_{53} - \gamma_{21} \cdot p_{21}$$

Donde:
* $P(t)$: Nivel de Pluripotencia ($0.0 \to 1.0$).
* $S(t)$: Nivel de Senescencia celular defensiva.
* $I(t)$: Inestabilidad genómica acumulada.
* $\Omega_{\text{TP53}}$: Flag binario de integridad funcional del gen *TP53* ($1.0$ salvaje / restaurado, $0.0$ mutado p.R175H).
* $\text{Burden}$: Carga mutacional patogénica relativa extraída del VCF.

El sistema se resuelve numéricamente en cada ciclo temporal mediante integración de Runge-Kutta de 4º Orden (RK4):

$$y_{n+1} = y_n + \frac{\Delta t}{6} (k_1 + 2k_2 + 2k_3 + k_4)$$

---

### Métrica de Riesgo Oncogénico (R-Score)

Durante la reprogramación celular, la sobreexpresión ectópica del oncogén *c-Myc* en un fondo genético con pérdida de función de *TP53* impide el arresto fisiológico en $G_1/S$, conduciendo a transformación maligna. Esta dinámica se cuantifica mediante el **R-Score**:

$$R(t) = w_1 \cdot I(t) + w_2 \cdot \frac{[\text{c-Myc}](t)}{p_{53}(t) + \epsilon} + w_3 \cdot \text{Burden}$$

* **Interpretación Clínica del R-Score:**
  * $R < 0.25$: **STABLE / LOW RISK** (Reprogramación fisiológica viable, arresto p53 activo).
  * $0.25 \le R \le 0.50$: **MODERATE RISK** (Estrés replicativo transitorio, balance homeostático).
  * $R > 0.50$: **HIGH RISK / MALIGNANT DRIFT** (Colisión crítica: proliferación incontrolada, aneuploidía y teratogénesis).

---

### Red de Enrutamiento MoE y Subredes de Expertos

Cada variante genómica se proyecta en un vector de características de 16 dimensiones:

$$x = \phi(\text{pos}, \text{chrom}, \text{effect}, \text{AF}, \text{DP}, \dots) \in \mathbb{R}^{16}$$

La compuerta del MoE (Gating Network) calcula la afinidad sobre las 4 subredes biológicas:

$$H = W_{\text{gate}} \cdot x + b_{\text{gate}}$$

$$G(x) = \text{Softmax}(H) \in \mathbb{R}^4$$

Las 4 subredes de expertos corresponden a las vías moleculares clave:
1. **Experto 0 (DNA Repair):** Vías de escisión de nucleótidos (NER) y recombinación homóloga.
2. **Experto 1 (Cell Cycle Checkpoints):** Regulación ciclinas/CDK, arresto $G_1/S$ y $G_2/M$.
3. **Experto 2 (p53 Apoptosis):** Señalización pro-apoptótica (Bax, Puma, Caspasa-3).
4. **Experto 3 (Stem Pluripotency):** Red de transcripción Nanog, Oct4 y remodelación de cromatina.

El motor selecciona dinámicamente las dos subredes con mayor ponderación (Top-2 Gating):

$$y = \sum_{i \in \text{Top-2}} G(x)_i \cdot E_i(x)$$

---

### Motor de Intervención CRISPR-Cas9 y Homología HDR

1. **Búsqueda PAM:** Escanea la secuencia diana identificando motivos canónicos `5'-NGG-3'` (hebra codificante) y motivos `5'-CCN-3'` (hebra molde / reversa).
2. **Scoring de Eficacia On-Target:** Evalúa el espaciador de 20 nt mediante la heurística Doench / Rule Set 2:
   $$\text{Score} = f(\text{Contenido GC}, \text{Polímeros T}, \text{Afinidad PAM adyacente})$$
3. **Distancia de Hamming ARM NEON (Cas-OFFinder):** Búsqueda de homología en ventanas genómicas para predecir cortes no deseados (*off-targets*).
4. **Reparación HDR in-silico (`[E]`):**
   * Reemplaza el codón mutado (ej. *TP53* c.524G>A; p.Arg175His -> CGC salvaje).
   * Restaura $\Omega_{\text{TP53}} = 1.0$ en el ODE de Yamanaka.
   * El R-Score se desploma instantáneamente a niveles seguros ($< 0.15$), estabilizando el cultivo.

---

## 5. Gemelo Digital Celular Espacial (Living Cell Digital Twin)

### Modelo Basado en Agentes (ABM)

El módulo `CellTwinSimulator` simula un cultivo en vitro 2D dentro de una placa de Petri virtual de radio $R = 45.0\text{ mm}$ que aloja hasta $3.000$ células simultáneas.

Cada célula individual (`CellAgent`) contiene una máquina de estados finitos del ciclo celular:

$$G_0 \text{ (Quiescencia)} \longrightarrow G_1 \text{ (Crecimiento)} \longrightarrow S \text{ (Síntesis ADN)} \longrightarrow G_2 \text{ (Premitótico)} \longrightarrow M \text{ (Mitosis)}$$

Al completar la fase $M$, la célula madre genera una célula hija por gemación mitótica en una orientación elástica con separación de masas, o sufre apoptosis / senescencia según su nivel de daño en ADN y actividad de p53.

```text
                                +-------------------+
                                | G0 (Quiescencia)  |
                                +---------+---------+
                                          |
                                          v
      +-------------------------> G1 (Crecimiento)  <------------------------+
      |                                   |                                  |
      |                                   v                                  |
      |                            S (Síntesis)                              |
      |                                   |                                  |
      |                                   v                                  |
      |                           G2 (Premitótico)                           |
      |                                   |                                  |
      |                                   v                                  |
      |                              M (Mitosis)                             |
      |                             /           \                            |
      |                            /             \                           |
      |      [cleavage_state = 1] v               v [p53 > 0.92 / Daño UV]   |
      |             Anafase (Dumbbell)       Apoptosis (Corpse)              |
      |                           |               |                          |
      |      [cleavage_state = 2] v               v                          |
      |             Telofase (Doblete)       Degradación / Lisis             |
      |                           |                                          |
      +---------------------------+------------------------------------------+
                               Mitosis completada
                            (Célula Madre + Hija)
```

---

### Controlador Temporal Multi-Escala (Time Warp)

La simulación permite alterar dinámicamente la escala del tiempo a través de 8 órdenes de magnitud discretos mediante las teclas `[` (frenar) y `]` (acelerar):

| Escala (`TimeWarpScale`) | Multiplicador | Dominio Biológico Modelado |
| :--- | :---: | :--- |
| **Seconds** | $1.0\times$ | Tiempo real, difusión física molecular y rebote elástico. |
| **Minutes** | $60.0\times$ | Señalización paracrina, fosforilación y vías quinasas. |
| **Hours** | $3,600.0\times$ | Transición de fases del ciclo celular ($G_1 \to S \to G_2 \to M$). |
| **Days** | $86,400.0\times$ | Cinética de reprogramación de Yamanaka (inducción OSKM). |
| **Weeks** | $604,800.0\times$ | Crecimiento y expansión de colonias pluripotentes iPSC. |
| **Months** | $2,592,000.0\times$ | Deriva genómica, selección clonal y acumulación mutacional. |
| **Years** | $31,536,000.0\times$ | Senescencia replicativa celular acumulada (acortamiento telómeros). |
| **Decades** | $315,360,000.0\times$ | Reloj epigenético de metilación e inestabilidad a largo plazo. |

* **Reloj Biológico Formateado:** Registra el tiempo simulado exacto: `T+ 000y:00m:04d:18h:32m:00s`.
* **Estabilizador Biofísico:** La cinemática espacial y el ciclo celular incorporan un factor adaptativo `bio_speed_boost` para asegurar transiciones morfológicas fluidas sin sobrecarga de cómputo en fotogramas de gran escala.

---

### Visualización Citomorfológica en Canvas Braille y Microscopio Virtual de Alta Magnificación

Sobre el `Canvas` de FTXUI ($114 \times 70$ sub-píxeles Braille), la placa de Petri se modela con un microscopio virtual con zoom óptico continuo/discreto ($1.0\times$ a $8.0\times$) y navegación pan:

* **Microscopio Virtual (`[Z / z / 0]`):**
  * Campo amplio ($1.0\times$): Visión panorámica completa del cultivo y bordes de la placa de Petri.
  * Alta magnificación ($\ge 2.0\times$): Resolución ultraestructural sub-celular detallada.
  * Paneo de cámara (`[Flechas]` / `[W, A, S]`): Desplaza el campo de visión sobre cualquier sector del cultivo con fijación automática de márgenes (*camera clamp*).
* **Morfología por Fenotipo Celular:**
  * **Fibroblasto Somático (`SomaticFibroblast`):** Morfología fusiforme alargada bipolar orientada según el vector de velocidad $(v_x, v_y)$ con núcleo oval central cian; en alta magnificación proyecta lamelipodios polares extendidos.
  * **Intermedio MET (`IntermediateReprog`):** Célula redondeada en transición mesénquima-epitelio en color amarillo con condensación cromatínica periférica.
  * **Colonia iPSC Pluripotente (`PluripotentIPSC`):** Colonias epiteliales compactas tipo cobblestone con núcleo esmeralda y nucleolo blanco prominente, delimitadas por un halo verde cohesivo.
  * **Clon Neoplásico Maligno (`TransformedMalignant`):** Pleomorfismo celular severo con núcleos gigantes hipercromáticos e irregulares en rojo carmesí, multinucleación atípica y pseudópodos invasivos.
  * **Cuerpo Apoptótico (`ApoptoticCorpse`):** Picnosis y blebbing celular disperso en múltiples vesículas y fragmentos picnóticos en gris pizarra.
  * **Surco de División Mitótica (`cleavage_state`):** Constricción en reloj de arena (dumbbell) con huso acromático visible en anafase tardía y dobletes celulares con anillo de telofase/midbody previos a la citocinesis en magenta y blanco.

---

### Acelerador Cinético en Tiempo Real (1x..32x)

A diferencia del **Time Warp** biológico (que altera el diferencial $dt$ de segundos a décadas), el **Multiplicador de Velocidad en Tiempo Real** acelera el throughput físico y de integración de la simulación mediante sub-stepping por fotograma:

| Multiplicador Real | Frecuencia de Cómputo | Capacidad de Simulación |
| :---: | :---: | :--- |
| **$1\times$** | $20\text{ Hz}$ | Velocidad estándar en tiempo real (~20 fotogramas de interfaz por segundo). |
| **$2\times$** | $40\text{ Hz}$ | $2$ pasos de física celular y solver ODE por fotograma. |
| **$4\times$** | $80\text{ Hz}$ | Aceleración media para observar migración celular y confluencia rápida. |
| **$8\times$** | $160\text{ Hz}$ | Alta velocidad para dinámicas de división y reprogramación activa. |
| **$16\times$** | $320\text{ Hz}$ | Dinámica acelerada para ensayos de apoptosis masiva por irradiación UV. |
| **$32\times$** | $640\text{ Hz}$ | Máximo rendimiento cinético en Apple Silicon para simulación ultra-rápida. |

* **Atajos:** `[ < ]` / `[ , ]` para decrementar la velocidad y `[ > ]` / `[ . ]` para incrementarla.

## 6. Interoperabilidad: Modo Headless, Inspector IA y Datasets CSV

GenomeForge-TUI está diseñado para operar tanto de forma interactiva en la terminal como de forma autónoma integrada con agentes de inteligencia artificial, pipelines bioinformáticos y scripts de análisis.

### 1. Ejecución Headless
Permite ejecutar simulaciones aceleradas sin entorno de interfaz gráfica:
```bash
./build/genome_forge --vcf data/1000genomes_clinical_panel.vcf --headless --ticks 50 --export-csv
# o utilizando el atajo del Makefile:
make headless
```

### 2. Inspector de Estado en Vivo (IPC para Agentes IA)
La aplicación serializa atómicamente su estado global en `data/live_state.json` cada 2 ticks de cómputo (utilizando escritura a archivo temporal `.tmp` y renombrado atómico para garantizar cero corrupciones por concurrencia).

Cualquier agente IA o script de terminal puede inspeccionar el estado actual ejecutando:
```bash
make inspect
# o directamente:
./build/genome_forge --inspect data/live_state.json
```

El informe JSON proporciona:
* Índices de pestaña y locus clínico activo (*TP53*, *CDKN2A*, etc.).
* Telemetría de hardware (chips M-Series, uso de CPU, RAM unificada).
* Variables continuas del solver ODE ($P, S, I, p_{53}, p_{21}, R\text{-Score}$).
* Estadísticas del Router MoE y latencia de aceleradores (vDSP / Metal GPU).
* Candidatos CRISPR seleccionados y scoring de eficacia.
* Censo biológico del gemelo digital (células por fenotipo, confluencia, índice mitótico, escala temporal activa).

### 3. Exportador Científico a CSV
Al presionar la tecla `[X]` en la TUI o activar `--export-csv`, se generan en la carpeta `data/`:
1. `export_cell_census_<timestamp>.csv`: Confluencia, índice mitótico, recuento poblacional y R-Score medio en función del tiempo.
2. `export_cell_agents_<timestamp>.csv`: Censo microscópico unicelular con coordenadas $(X, Y)$, velocidades $(V_x, V_y)$, fases del ciclo ($G_0, G_1, S, G_2, M$), fenotipo, estado de clivaje y nivel de $p_{53}$.
3. `export_genomic_variants_<timestamp>.csv`: Catálogo de variantes del panel con alelos REF/ALT, profundidad $DP$, frecuencia $AF$ y scoring ClinVar.

---

## 6. Motor de Oncología de Precisión y Farmacogenómica

GenomeForge-TUI incorpora un motor de oncología de precisión (`ActionabilityEngine`) y farmacodinamia computacional que traduce variantes genómicas en prescripciones terapéuticas dirigidas y pronósticos de respuesta celular.

```text
+-----------------------------------------------------------------------------------------------+
|                       ONCOLOGY ACTIONABILITY & PHARMACOGENOMICS PIPELINE                      |
|                                                                                               |
|  [VCF Patient Data] ──> [LocusCatalog 12 Genes] ──> [Oncogenic Driver Matcher]                |
|                                                               │                               |
|                                                               v                               |
|                     ┌──────────────────────────────────────────────────┐                      |
|                     │     ActionabilityEngine Knowledge Base (FDA/NCCN)│                      |
|                     │  - Level 1: FDA Approved Targeted Therapies      │                      |
|                     │  - Level 2: NCCN Guideline Standards of Care     │                      |
|                     │  - Level 3: Phase II/III Clinical Trials         │                      |
|                     │  - Level 4: Preclinical CRISPR-Cas9 HDR Rescue   │                      |
|                     │  - Resistance: Secondary Gatekeeper Markers      │                      |
|                     └────────────────────────┬─────────────────────────┘                      |
|                                              │                                                |
|                      ┌───────────────────────┴───────────────────────┐                        |
|                      v                                               v                        |
|        ┌───────────────────────────┐                   ┌───────────────────────────┐          |
|        │  In-Silico Pharmacodynamics│                   │ Molecular Tumor Board (MTB)│         |
|        │  - Drug Inoculation       │                   │ - Markdown Clinical Report│          |
|        │  - Tumor Regression Rate  │                   │ - JSON Telemetry Interop  │          |
|        │  - Clonal Resistance Emergence│               │ - CRISPR HDR Rescue Plan  │          |
|        └───────────────────────────┘                   └───────────────────────────┘          |
+-----------------------------------------------------------------------------------------------+
```

### Niveles de Evidencia Clínica y Base Curada

El motor clasifica cada recomendación bajo una taxonomía de evidencia estandarizada:

* **Nivel 1 (Aprobado FDA / EMA):**
  * **Osimertinib:** Inhibidor irreversible de tirosina quinasa (TKI) de 3ª generación para mutaciones sensibilizantes de *EGFR* (p.L858R, deleciones del exón 19).
  * **Olaparib / Talazoparib:** Inhibidores de PARP1/2 por letalidad sintética en células con deficiencia de recombinación homóloga por mutaciones patogénicas en *BRCA1* o *BRCA2*.
  * **Sotorasib / Adagrasib:** Inhibidores covalentes del bolsillo switch-II de *KRAS* p.G12C en cáncer de pulmón no microcítico.
  * **Dabrafenib + Trametinib:** Régimen dual de inhibición de cinasas sobre *BRAF* p.V600E y downstream *MEK1/2*.
  * **Alpelisib:** Inhibidor selectivo de la subunidad catalítica p110-alfa de PI3K para mutaciones activadoras en *PIK3CA*.
  * **Trastuzumab Deruxtecan (T-DXd):** Anticuerpo conjugado con carga de inhibidor de topoisomerasa I dirigido contra sobreexpresión/amplificación de *ERBB2* (HER2).
  * **Pembrolizumab:** Inmunoterapia anti-PD-1 activada automáticamente ante alta carga mutacional tumoral ($\text{TMB} \ge 10\text{ mut/Mb}$ o $\text{Drivers} \ge 4$).
* **Nivel 2 (Guías NCCN / Estándar de Cuidado):**
  * **Gefitinib:** TKI reversible de 1ª generación para *EGFR* (contraindicado formalmente en presencia de resistencia secundaria).
  * **Everolimus:** Inhibidor alostérico de mTORC1 en tumores con pérdida de función del supresor *PTEN*.
* **Nivel 3 (Ensayos Clínicos Fase II/III):**
  * **APR-246 (Eprenetapopt):** Reactivador conformacional de p53 que se une covalentemente a cisteínas del dominio de unión a ADN en mutaciones hotspot conformacionales (*TP53* p.R175H).
* **Nivel 4 (Preclínico / In-Silico):**
  * **CRISPR-Cas9 HDR Gene Rescue:** Edición genética de precisión mediante rotura de doble cadena adyacente y plantilla donante de hebra simple (ssODN) para corregir el codón mutado.
* **Alertas de Resistencia Farmacogenómica:**
  * Detección específica de la mutación *gatekeeper* *EGFR* p.T790M. El sistema emite una alerta roja inmediata, contraindica el uso de Gefitinib o inhibidores reversibles y prioriza Osimertinib.

### Farmacodinamia In-Silico en el Gemelo Digital

La pestaña `< 8: ONCO THERAPY >` y los controles de teclado permiten inocular fármacos en el gemelo digital celular (`CellTwinSimulator`):

* **Cinéticas de Lisis Selectiva:** Los fármacos activos ejercen una presión de eliminación dependiente de la dosis sobre clones neoplásicos transformados.
* **Letalidad Sintética:** Las células portadoras de alelos mutados en *BRCA1* experimentan colapso replicativo letal ante la presencia de inhibidores PARP.
* **Restauración de Checkpoints:** La inoculación de APR-246 restaura gradualmente los niveles de $p_{53}$, reactivando la apoptosis funcional en clones desregulados.
* **Emergencia Estocástica de Resistencia Clonal:** Cada división mitótica bajo presión selectiva tiene una probabilidad del 2.5% de generar un clon resistente (`drug_resistant = true`), visualizado en la placa de Petri y cuantificado en la métrica de regresión tumoral.
* **Lavado de Fármacos (`[W]`):** Retira los compuestos activos y restaura la cinética basal del cultivo.

### Reportes Clínicos de Comité Molecular de Tumores (MTB)

Al pulsar la tecla `[O]` en la TUI o mediante el parámetro CLI `--report-onco <ruta.md>`, GenomeForge genera dos artefactos complementarios:

1. **Informe Clínico en Markdown (`.md`):**
   * Metadatos del paciente y versión del motor de referencia GRCh38.
   * Diagnóstico genómico resumido (TMB, R-Score, estado de *TP53*).
   * Tabla de alteraciones moleculares conductoras con coordenadas y codones.
   * Prescripción de terapias dirigidas clasificadas por nivel de evidencia FDA/NCCN.
   * Advertencias farmacogenómicas de resistencia y protocolos de monitoreo clonal.
   * Plan de rescate genético por edición CRISPR-Cas9 HDR (guía protospacer 20 nt, PAM y molde ssODN).
   * Pronóstico farmacodinámico proyectado por el Gemelo Digital.
2. **Exportación Estructurada en JSON (`.json`):**
   * Esquema JSON completo para interoperabilidad con agentes de IA, registros clínicos electrónicos (EHR) y pipelines NGS automatizados.

---

## 7. Interoperabilidad: Modo Headless, Inspector IA, Datasets CSV y Reportes

GenomeForge-TUI está diseñado para integrarse de forma nativa en flujos de trabajo autónomos y pipelines de bioinformática:

### 1. Ejecución Headless para Pipelines Automatizados
```bash
./build/genome_forge --vcf data/sample_patient_oncology.vcf --headless --ticks 50 --export-csv --report-onco data/molecular_tumor_board_report.md
# o utilizando el atajo del Makefile:
make test-onco
```

### 2. Inspector de Estado en Vivo (IPC para Agentes IA)
La aplicación serializa atómicamente su estado global en `data/live_state.json` cada 2 ticks de cómputo (utilizando escritura a archivo temporal `.tmp` y renombrado atómico para garantizar cero corrupciones por concurrencia).

Cualquier agente IA o script de terminal puede inspeccionar el estado actual ejecutando:
```bash
make inspect
# o directamente:
./build/genome_forge --inspect data/live_state.json
```

El informe JSON proporciona:
* Índices de pestaña y locus clínico activo (*TP53*, *CDKN2A*, *EGFR*, etc.).
* Telemetría de hardware (chips M-Series, uso de CPU, RAM unificada).
* Variables continuas del solver ODE ($P, S, I, p_{53}, p_{21}, R\text{-Score}$).
* Estadísticas del Router MoE y latencia de aceleradores (vDSP / Metal GPU).
* Candidatos CRISPR seleccionados y scoring de eficacia.
* Censo biológico del gemelo digital (células por fenotipo, confluencia, índice mitótico, escala temporal activa).
* Módulo `oncology_actionability` con fármacos activos, dosis en $\mu\text{M}$, tasa de regresión tumoral y clones resistentes.

### 3. Exportador Científico a CSV
Al presionar la tecla `[X]` en la TUI o activar `--export-csv`, se generan en la carpeta `data/`:
1. `export_cell_census_<timestamp>.csv`: Confluencia, índice mitótico, recuento poblacional y R-Score medio en función del tiempo.
2. `export_cell_agents_<timestamp>.csv`: Censo microscópico unicelular con coordenadas $(X, Y)$, velocidades $(V_x, V_y)$, fases del ciclo ($G_0, G_1, S, G_2, M$), fenotipo, estado de clivaje y nivel de $p_{53}$.
3. `export_genomic_variants_<timestamp>.csv`: Catálogo de variantes del panel con alelos REF/ALT, profundidad $DP$, frecuencia $AF$ y scoring ClinVar.

### 4. Automatización y Verificación de Workflows Clínicos
Para auditar y validar los 5 protocolos clínicos de extremo a extremo sin interfaz gráfica:
```bash
make test-workflows
# o directamente:
./build/genome_forge --test-workflows
```
El motor ejecuta de forma secuencial la ingesta de muestras, diagnósticos moleculares, cribado de fármacos FDA Nivel 1, inoculación en el gemelo digital, simulación biofísica, verificación de resistencia clonal y emisión de reportes MTB en Markdown y JSON.

---

## 8. Integración Multi-Locus con el Proyecto 1000 Genomas (12 Loci Clínicos)

Para validar las simulaciones frente a datos genómicos reales de población humana y paneles de oncología de precisión, GenomeForge-TUI incorpora un catálogo multi-locus ampliado a 12 genes clínicos en GRCh38 / hg38:

| Locus | Coordenada Genómica (GRCh38) | Relevancia Biológica y Farmacogenómica |
| :--- | :--- | :--- |
| **TP53** | Chr17: 7,668,402 - 7,687,550 | Guardián del genoma; mutaciones hotspot (R175H) restaurables vía APR-246 y CRISPR. |
| **CDKN2A** | Chr9: 21,967,751 - 21,995,300 | Barrera fisiológica de senescencia celular (p16INK4a / p14ARF). |
| **PTEN** | Chr10: 87,862,566 - 87,971,930 | Fosfatasa supresora de proliferación; pérdida diana de Everolimus (mTORC1). |
| **MYC** | Chr8: 127,735,434 - 127,742,951 | Proto-oncogén y factor de reprogramación OSKM (c-Myc). |
| **BRCA1** | Chr17: 43,044,295 - 43,125,483 | Reparación de ADN por HR; diana de letalidad sintética con inhibidores PARP (Olaparib). |
| **EGFR** | Chr7: 55,019,017 - 55,211,628 | Receptor tirosina quinasa diana de Osimertinib (L858R) y monitoreo de T790M. |
| **POU5F1** | Chr6: 31,164,335 - 31,170,694 | Factor de transcripción maestro nuclear de pluripotencia (Oct4). |
| **SOX2** | Chr3: 181,711,925 - 181,714,436 | Factor central del mantenimiento de células madre pluripotentes. |
| **KRAS** | Chr12: 25,204,789 - 25,250,929 | GTPasa oncogénica; diana de inhibidores covalentes específicos G12C (Sotorasib, Adagrasib). |
| **BRAF** | Chr7: 140,719,327 - 140,924,929 | Serina/treonina quinasa diana de terapia combinada Dabrafenib + Trametinib (V600E). |
| **PIK3CA** | Chr3: 179,148,114 - 179,240,093 | Subunidad catalítica p110-alfa diana de Alpelisib en tumores con activación de PI3K. |
| **ERBB2** | Chr17: 39,687,914 - 39,730,415 | Tirosina quinasa HER2 diana de anticuerpos conjugados (Trastuzumab Deruxtecan). |

* **Mapeo Inteligente por Coordenadas Genómicas (GRCh38):** Los archivos VCF brutos del consorcio 1000 Genomes no contienen etiquetas informativas `GENE=` en la columna `INFO`. `VcfParser` y `LocusCatalog` resuelven la pertenencia a los loci clínicos mediante resolución de intervalos por coordenadas cromosómicas en GRCh38.
* **Clasificación Poblacional y Controles Sanos ($AF \ge 0.01$):** Los donantes del 1000 Genomes Project son controles poblacionales sanos. El motor bioinformático clasifica automáticamente las variantes con frecuencia alélica común ($AF \ge 0.01$) sin anotación clínica de riesgo como polimorfismos naturales benignos (`[+] STABLE / WILD-TYPE`), resultando en $0$ colisiones de Yamanaka y un R-Score basal verde (< 0.10).
* **Conmutador de Muestras (`[M]`):** Permite alternar en tiempo real entre:
  * **Modo Paciente Clínico Oncológico (`data/sample_patient_oncology.vcf`):** Panel multiobjetivo con variantes patogénicas en *EGFR*, *BRCA1*, *KRAS*, *BRAF*, *PTEN*, *PIK3CA*, *ERBB2* y *TP53* para evaluar terapias dirigidas.
  * **Modo Paciente Resistencia TKI (`data/sample_patient_resistance.vcf`):** Perfil con mutación sensibilizante *EGFR* L858R más mutación *gatekeeper* T790M para probar contraindicaciones y cambio a Osimertinib.
  * **Modo 1000 Genomas Control Sano (`data/1000genomes_healthy_reference.vcf`):** Perfil de control poblacional sano (`HG00096`) con polimorfismos naturales benignos y estabilidad genómica garantizada.
* **Conmutación en Caliente de Locus (`[L]`):** Permite ciclar en tiempo real entre los 12 genes clínicos en la pestaña `[2] VCF EXPLORER`, reconfigurando automáticamente el contexto mutacional, el solver ODE y el diseño de guías CRISPR.

---

## 9. Requisitos, Compilación y Referencia CLI

### Requisitos del Sistema
* **Hardware:** Mac con procesador Apple Silicon (M1, M2, M3, M4 o variantes Pro/Max/Ultra).
* **Sistema Operativo:** macOS 14 (Sonoma) o macOS 15 (Sequoia) o superior.
* **Compilador:** Clang compatible con C++20 (incluido en Xcode Command Line Tools).
* **Gestor de Construcción:** CMake 3.20+ y Make.
* **Frameworks Nativos:** `Accelerate`, `Metal`, `Foundation`.

### Compilación Rápida
```bash
# Compilar binario optimizado (-O3) para Apple Silicon ARM64
make build

# Generar reporte clínico de Comité Molecular de Tumores (MTB)
make report

# Ejecutar suite de pruebas oncológicas y resistencia en modo headless
make test-onco

# Inspeccionar telemetría viva en formato JSON
make inspect

# Limpiar artefactos de compilación
make clean
```

### Opciones de Línea de Comandos (CLI)
```text
Uso:
  ./build/genome_forge [opciones]

Opciones Principales:
  --vcf <archivo.vcf>    Ruta al archivo VCF de entrada (default: data/sample_patient.vcf)
  --headless             Ejecuta la simulación en modo headless (sin interfaz TUI)
  --ticks <N>            Número de iteraciones a ejecutar en modo headless (default: 30)
  --snapshot <path>      Ruta del archivo JSON de telemetría (default: data/live_state.json)
  --inspect [path]       Inspecciona el estado de la TUI en vivo y muestra reporte HUD
  --export-csv           Fuerza la exportación de datasets CSV al finalizar
  --report-onco [path]   Genera el reporte clínico de Comité Molecular de Tumores (.md y .json)
  --help, -h             Muestra la ayuda detallada del sistema
```

---

## 10. Tabla de Atajos de Teclado

| Tecla / Atajo | Función y Efecto Biológico |
| :---: | :--- |
| **`1` .. `8`** | Salto directo a cualquiera de las 8 pestañas del Cyber-HUD (incluyendo Onco Therapy). |
| **`Tab`** | Ciclar secuencialmente a la siguiente pestaña (módulo 8). |
| **`<` / `,`** | **Disminuir Velocidad Real:** Reduce los sub-pasos del motor físico/ODE ($32\times \to 1\times$, $640\text{ Hz} \to 20\text{ Hz}$). |
| **`>` / `.`** | **Aumentar Velocidad Real:** Acelera los sub-pasos del motor físico/ODE ($1\times \to 32\times$, $20\text{ Hz} \to 640\text{ Hz}$). |
| **`Z`** | **Zoom In Microscopio:** Incrementa el aumento óptico de la placa de Petri ($1.0\times \to 8.0\times$). |
| **`z`** | **Zoom Out Microscopio:** Reduce el aumento óptico de la placa de Petri hacia el campo amplio. |
| **`0`** | **Reset Microscopio:** Restablece el zoom a $1.0\times$ y centra las coordenadas pan en $(0,0)$. |
| **`↑ ↓ ← →` / `W A S`** | **Pan Óptico (Tab 7):** Paneo espacial bidireccional sobre las células de la placa de Petri. |
| **`[`** | **Frenar Time Warp:** Reduce la escala temporal biológica (ej. Horas -> Minutos -> Segundos). |
| **`]`** | **Acelerar Time Warp:** Incrementa la escala temporal biológica (ej. Días -> Semanas -> Años -> Décadas). |
| **`O`** | **Reporte MTB:** Exporta el informe clínico de Comité Molecular de Tumores a Markdown y JSON. |
| **`T`** | **Inocular TKI:** Administra Osimertinib (1.50 uM) contra clones oncogénicos dirigidos por EGFR. |
| **`P`** | **Inocular PARP:** Administra Olaparib (1.50 uM) por letalidad sintética en mutaciones de BRCA1. |
| **`K`** | **Inocular p53:** Administra APR-246 (1.50 uM) para replegamiento conformacional de p53 mutante. |
| **`W`** | **Washout Farmacológico:** Suspende y lava cualquier terapia activa (control vehículo 0.0 uM). |
| **`X`** | **Exportar Datasets CSV:** Guarda censo, células individuales y variantes en `data/`. |
| **`L`** | **Conmutar Locus 1000G:** Alterna entre los 12 loci clínicos (*TP53*, *CDKN2A*, *PTEN*, *MYC*, *BRCA1*, *EGFR*, *OCT4*, *SOX2*, *KRAS*, *BRAF*, *PIK3CA*, *ERBB2*). |
| **`M`** | **Conmutar Muestra:** Alterna entre Muestra Oncológica, Control de Resistencia y Control Sano 1000 Genomas. |
| **`D`** | **Doxiciclina (OSKM):** Inyecta o suspende el inductor de reprogramación en el cultivo celular. |
| **`C`** | **CRISPR In-Situ:** Administra terapia génica en la placa de Petri para corregir mutaciones de *TP53*. |
| **`U`** | **Estrés Genotóxico:** Dispara un pulso de radiación UV para poner a prueba los puntos de control. |
| **`+` / `-`** | Inocular / aspirar células del cultivo virtual de la placa de Petri (±25 células). |
| **`↑` / `↓` o `k` / `j`** | Navegar verticalmente por listas de variantes VCF o candidatos gRNA CRISPR (Tabs 2 y 5). |
| **`E`** | **Edición CRISPR (HDR):** Aplica corrección in-silico por homología sobre la variante VCF seleccionada. |
| **`R`** | **Reset Basal:** Restaura el cultivo celular (en Tab 7) o reinyecta la mutación en el genoma original. |
| **`Espacio`** | Pausar o reanudar el avance temporal de los motores biológicos y la biofísica. |
| **`Q` o `Esc`** | Salida limpia del programa. |

---

## 11. Documentación y Referencias Técnicas

* [`docs/SPECS.md`](docs/SPECS.md): Especificación técnica formal, structs alineados a caché ARM64, derivaciones matemáticas de ODEs, kernels de Metal y arquitectura del sistema.
* [`docs/BITACORA.md`](docs/BITACORA.md): Bitácora cronológica viva de ingeniería con todos los registros de desarrollo (*devlog*), decisiones técnicas y benchmarks de rendimiento.



