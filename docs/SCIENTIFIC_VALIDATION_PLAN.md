# GenomeForge — Plan de recuperación y validación científica

**Estado:** propuesta técnica, pendiente de implementación y revisión científica.  
**Alcance:** investigación computacional y educación; **no** diagnóstico, prescripción, selección terapéutica ni predicción clínica validada.  
**Base de auditoría:** inspección estática del repositorio (noviembre de 2026). No se compiló, ejecutó ni validó experimentalmente el software.

## 1. Resumen ejecutivo

GenomeForge contiene software sustancial: parser VCF, motor de reprogramación OSKM con integración RK4, modelos celulares espaciales, diseño heurístico de guías CRISPR, reglas de oncología, workflows, exportación, almacenamiento SQLite, TUI y componentes acelerados con Apple Silicon. La existencia de código y resultados numéricos **no demuestra** fidelidad biológica ni utilidad clínica.

**Decisión propuesta:** conservar el proyecto, detener la expansión funcional y establecer una ruta de validación incremental centrada en **un fenómeno biológico bien delimitado**. Separar rigurosamente simulación ilustrativa, cálculo bioinformático verificable, hipótesis mecanística y evidencia experimental.

## 2. Hallazgos de auditoría y riesgos

| Prioridad | Área / referencia | Observación verificable | Riesgo | Acción |
|---|---|---|---|---|
| P0 | `src/yamanaka/yamanaka_sim.cpp` | El R-Score combina inestabilidad, impulso oncogénico y carga patogénica con coeficientes constantes (0.45/0.35/0.20) y límites 0.05–0.98. | Se puede interpretar erróneamente como probabilidad clínica. | Renombrar como índice interno exploratorio; retirar etiquetas de riesgo clínico hasta calibrarlo. |
| P0 | `src/digital_twin/cell_twin.cpp` | `apply_crispr_rescue()` elimina la mutación TP53 y fija p53=0.88 y r_score=0.12; cambia estados malignos mediante reglas directas. | El «rescate» favorable está impuesto, no predicho. | Identificar el flujo como demostración; sustituir por modelo explícito con incertidumbre y sin promesas terapéuticas. |
| P0 | `src/crispr/crispr_engine.cpp` | `mock_off_target_background_` es un conjunto pequeño de secuencias simuladas; la eficiencia y seguridad se puntúan con heurísticas. | Los porcentajes no representan tasas experimentales ni cribado genómico exhaustivo. | Etiquetar `demo`; no reportar «seguridad» como validada; comparar con referencias de evaluación publicadas. |
| P0 | `src/oncology/actionability_engine.cpp` | Base de fármacos y evidencia embebida en código, con puntuaciones constantes. | Evidencia desactualizada, contexto clínico insuficiente, falsa autoridad. | Desactivar lenguaje de prescripción y recomendaciones clínicas; separar base de conocimiento versionada y trazable. |
| P1 | `src/workflow/workflow_engine.cpp` | Workflows denominados «clínicos» generan informes tipo comité molecular. | Un informe sintético puede confundirse con consejo médico. | Marcar todos los resultados como simulación/no uso clínico; renombrar workflows. |
| P1 | `src/digital_twin/cell_twin.cpp` | Interacciones celulares por bucles anidados, coste aproximado O(N²). | Escalabilidad limitada; optimizar sin validar no mejora ciencia. | Perfilar y, solo después, evaluar rejilla espacial / neighbor lists. |
| P1 | `src/yamanaka/yamanaka_sim.cpp` | RK4 está implementado; parámetros cinéticos constantes sin trazabilidad experimental visible. | Precisión numérica confundida con exactitud biológica. | Verificación de convergencia, sensibilidad y calibración con datos independientes. |
| P1 | `src/vcf/` | Parsing y clasificación de variantes propios. | Interpretaciones incorrectas de VCF, genoma de referencia, alelos y anotaciones. | Pruebas de conformidad y comparación con herramientas establecidas. |
| P2 | `docs/BITACORA.md`, `README.md` | Se describen componentes como «operativos» y reportes «clínicos». | Sobreafirmación de madurez científica. | Separar estado de ingeniería de estado de validación. |

**Advertencia:** los nombres de medicamentos, dosis y supuestas eficacias de los ejemplos existentes no constituyen recomendaciones. No utilizar datos reales identificables de pacientes para pruebas iniciales.

## 3. Niveles de evidencia obligatorios

Cada módulo, métrica, gráfico, exportación y reporte debe declarar un nivel:

- **L0 — Demostración:** valores sintéticos, efectos programados, mock data.
- **L1 — Verificación computacional:** algoritmo comprobado con casos conocidos, invariantes y pruebas reproducibles.
- **L2 — Validación retrospectiva:** comparación cuantitativa con datos experimentales publicados, fuera del conjunto de calibración.
- **L3 — Validación externa:** replicación independiente, distintas condiciones y datasets.
- **L4 — Validación para uso previsto:** solo tras revisión especializada, gobernanza y evaluación regulatoria correspondiente.

**Regla:** L1 no implica L2; ningún resultado de GenomeForge puede presentarse como clínicamente válido por el solo hecho de ejecutarse.

## 4. Fases y criterios de salida

### Fase 0 — Contención y reproducibilidad (P0)
**Objetivo:** impedir afirmaciones clínicas injustificadas y congelar una referencia auditable.

- [ ] Crear inventario de módulos, entradas, salidas, datasets y dependencias.
- [ ] Añadir aviso «RESEARCH / DEMONSTRATION ONLY — NOT FOR CLINICAL USE» a TUI, CLI, CSV y Markdown exportado.
- [ ] Distinguir datos sintéticos de reales en metadatos; no tratar muestras sintéticas como pacientes.
- [ ] Desactivar expresiones «FDA-prescribed», «safe», «cure», «tumor board recommendation» en resultados automáticos sin revisión.
- [ ] Fijar versiones de dependencias, compilador y datasets; documentar build Apple Silicon.
- [ ] Capturar hashes, seeds, parámetros, versión del modelo y fecha en cada ejecución.
- [ ] Ejecutar compilación limpia y smoke tests, registrar resultados reales; no asumir que pasan.

**Salida:** cada resultado se identifica inequívocamente como simulación y se reproduce con la misma configuración.

### Fase 1 — Verificación de ingeniería (P0–P1)
**Objetivo:** asegurar que los algoritmos implementan sus propias especificaciones.

- [ ] Pruebas unitarias de VCF: multi-alélicos, indels, genotipos, contigs, coordenadas, valores faltantes, REF/ALT, cabeceras y variantes estructurales fuera de alcance.
- [ ] Contrastar parsing/anotación con herramientas de referencia y conjuntos de prueba públicos, documentando diferencias.
- [ ] RK4: problemas analíticos conocidos, estudio de convergencia al reducir dt, invariantes y estabilidad numérica.
- [ ] CRISPR: casos de PAM en ambas hebras, complemento inverso, extremos, secuencias ambiguas, mismatch y longitudes.
- [ ] Pruebas de determinismo, seeds y sensibilidad de parámetros del gemelo digital.
- [ ] Validar formatos de exportación, etiquetas de unidades, provenance y consistencia entre CLI/TUI/reportes.
- [ ] CI con compilación y tests en macOS ARM64; análisis estático, sanitizers cuando sea viable.

**Salida:** suite automatizada con resultados reproducibles, cobertura de casos límite y errores conocidos documentados.

### Fase 2 — Selección de un único modelo científico (P0)
**Objetivo:** escoger un problema que pueda contrastarse con evidencia, no un «gemelo celular universal».

**Candidato inicial:** dinámica simplificada de reprogramación OSKM o respuesta poblacional a una perturbación experimental **sin inferencias terapéuticas**.

- [ ] Formular pregunta científica concreta, observable y falsable.
- [ ] Definir variables medidas, unidades, escala temporal y condiciones iniciales.
- [ ] Identificar un dataset público experimental con licencia y metadatos adecuados.
- [ ] Separar conjunto de calibración y conjunto de evaluación antes de ajustar parámetros.
- [ ] Establecer baseline trivial y modelo de referencia de literatura.
- [ ] Documentar ecuaciones, hipótesis, fuentes de parámetros y aspectos omitidos.
- [ ] Aprobar selección con revisión de alguien con experiencia en biología experimental.

**Salida:** ficha del modelo, protocolo predefinido y dataset versionado; no avanzar sin datos comparables.

### Fase 3 — Calibración, incertidumbre y validación (P0)
- [ ] Estimar parámetros solo con datos de entrenamiento, incluyendo intervalos plausibles.
- [ ] Medir error en observaciones retenidas: MAE/RMSE u otra métrica justificada, más error del baseline.
- [ ] Realizar análisis de sensibilidad local/global e identificabilidad.
- [ ] Reportar intervalos de incertidumbre, no únicamente curvas deterministas.
- [ ] Analizar fallos por condición, lote, línea celular y perturbación cuando corresponda.
- [ ] Repetir evaluación con dataset independiente; registrar explícitamente cuándo no generaliza.
- [ ] Publicar protocolo reproducible, figuras, datos permitidos y resultados negativos.

**Salida:** informe con comparación cuantitativa frente a baseline y límites de extrapolación. Si no supera baseline, reducir o abandonar el modelo.

### Fase 4 — Rehabilitación gradual de subsistemas (P1)
**VCF:** biblioteca de parsing y anotación con conformidad demostrada; separar anotación de interpretación clínica.  
**CRISPR:** herramienta de exploración de secuencias y evaluación computacional comparativa, sin inferir seguridad real de edición.  
**Oncología:** catálogo de conocimiento con fuente, versión, fecha, tumor, variante, indicación, nivel de evidencia y revisión humana; nunca «prescribir» automáticamente.  
**Gemelo celular:** modelo específico de experimento con reglas explícitas y datos para calibración; retirar el rescate determinista de TP53 de cualquier modo de investigación.  
**R-Score:** mantener solo como variable interna no clínica o eliminarla; si se redefine, exigir constructo, validación y calibración.

**Salida:** cada subsistema muestra su nivel L0–L4 y su informe de validación.

### Fase 5 — Rendimiento, arquitectura y publicación (P2)
- [ ] Perfilar tiempo, memoria y escalado por número de células/variantes.
- [ ] Comparar CPU vs Metal/Accelerate con resultados numéricos equivalentes y benchmarks reproducibles.
- [ ] Sustituir O(N²) por índice espacial si el perfil lo justifica; medir coste/precisión.
- [ ] Separar `core`, `models`, `datasets`, `validation`, `ui`, `reports` y `experimental`.
- [ ] Documentar instalación, licencias de datos, límites, seguridad, reproducibilidad y contribución científica.
- [ ] Someter modelos a revisión técnica/biológica externa antes de afirmar capacidad predictiva.

**Salida:** release de investigación reproducible, con limitaciones explícitas y sin pretensiones clínicas.

## 5. Arquitectura objetivo

```text
Datasets versionados + procedencia
             |
      Ingesta y validación
             |
    Core computacional puro
      /        |         \
   VCF      Modelos     CRISPR
            (L0-L4)    (L0-L4)
      \        |         /
     Motor de experimentos
     (seeds, configs, runs)
             |
     Evaluación independiente
  (baselines, métricas, incertidumbre)
             |
       Reportes trazables
             |
         TUI / CLI
```

**Principios:** módulos sin efectos laterales cuando sea posible; unidades explícitas; parámetros configurables y versionados; resultados siempre vinculados a modelo, dataset y nivel de evidencia. Separar lógica de presentación, no usar etiquetas clínicas para variables exploratorias.

## 6. Estructura documental propuesta

```text
docs/
  SCIENTIFIC_VALIDATION_PLAN.md   # este plan
  MODEL_REGISTRY.md               # inventario, supuestos, nivel L0–L4
  DATA_PROVENANCE.md              # fuente/licencia/versión de datasets
  VALIDATION_PROTOCOL.md          # criterios y métricas predefinidos
  KNOWN_LIMITATIONS.md            # errores, no-usos y resultados negativos
validation/
  fixtures/
  baselines/
  experiments/
  reports/
```

No crear archivos vacíos solo por completar la estructura; incorporarlos al ejecutar cada fase.

## 7. Primer backlog concreto

1. **GF-001 (P0):** rotular todos los outputs y la UI como investigación/no uso clínico.
2. **GF-002 (P0):** inventariar y etiquetar valores mock, constantes y resultados forzados.
3. **GF-003 (P0):** quitar el carácter clínico del R-Score y de las recomendaciones automáticas.
4. **GF-004 (P0):** crear harness reproducible (config, seed, hash, versionado).
5. **GF-005 (P1):** test suite de VCF y comparaciones de referencia.
6. **GF-006 (P1):** pruebas de convergencia RK4 y sensibilidad.
7. **GF-007 (P1):** batería CRISPR para ambas hebras y casos límite; documentar fondo off-target simulado.
8. **GF-008 (P0):** seleccionar una pregunta biológica y dataset externo.
9. **GF-009 (P0):** baseline, protocolo y métricas de validación.
10. **GF-010 (P1):** informe reproducible con resultados positivos y negativos.

## 8. Criterios de no aceptación

No considerar un módulo «validado» por: compilar, pasar tests funcionales, generar un informe, producir una curva plausible, ejecutar más rápido en GPU, usar datos VCF reales sin ground truth, o citar estudios sin reproducir sus observaciones. No permitir que puntuaciones heurísticas se presenten como probabilidades clínicas ni que un efecto programado se describa como descubrimiento.

## 9. Preguntas abiertas

- ¿Cuál es la **pregunta biológica única** de la primera validación?
- ¿Qué mediciones publicadas existen para calibrar y evaluar el modelo?
- ¿Cuáles parámetros están respaldados por literatura y cuáles fueron elegidos heurísticamente?
- ¿Qué variantes y referencias genómicas soporta realmente el parser?
- ¿Qué resultados deben quedar como demos educativas en lugar de investigación?
- ¿Quién realizará revisión biológica independiente?

---

**Definición de éxito:** GenomeForge será valioso si sus afirmaciones son proporcionales a su evidencia. Un simulador pequeño, reproducible y cuantitativamente validado es preferible a una plataforma amplia que genere resultados clínicamente convincentes pero no comprobados.
