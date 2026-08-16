#!/usr/bin/env bash
# Executa o pipeline completo na ordem correta. Ver README.md para detalhes.
set -e
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
RAIZ="$(cd "$(dirname "$0")" && pwd)"
mkdir -p "$RAIZ/dados"

echo "== Etapa 1/5: pipeline principal em R (~7 min) =="
( cd "$RAIZ/01_pipeline_R" && Rscript metabric_pipeline.R )

echo "== Etapa 2/5: analises complementares (~40 min) =="
( cd "$RAIZ/02_analises_complementares" \
  && python3 analise_complementar.py \
  && python3 parte_b_farmaco.py \
  && python3 parte_c_rf.py )

echo "== Etapa 3/5: benchmark de algoritmos (~30 min) =="
( cd "$RAIZ/03_benchmark_ml" \
  && python3 ml_tarefa1_subtipo.py \
  && python3 ml_tarefa2_sobrevida.py \
  && python3 ml_tarefa2_final.py )

echo "== Etapa 4/5: treino final e exportacao (~5 min) =="
( cd "$RAIZ/04_treino_e_exportacao" && python3 treina_exporta.py )

echo "== Etapa 5/5: figuras e tabelas (~1 min) =="
( cd "$RAIZ/05_figuras" && python3 gera_figuras.py )

echo "== Concluido. Figuras em 05_figuras/, tabelas em 06_tabelas/. =="
