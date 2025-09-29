#!/bin/bash

# Anotación funcional con eggnog

source /home/alejandro.vazquez/Desktop/A.VAZQUEZ_TFM_2025_MBIF/TRANSCRIPTOMICS/CODE/config.env

#Primera fase de eggnog

#echo -e  "Comenzando proceso de anotación funcional... Generando alineaciones crudas y genes ortólogos... "

#emapper.py \
# -i "${SALMON_DIR}/SALMON_COMMON/all_proteins.faa" \
# --itype proteins \
# -o eggnog_annotation \
# --output_dir "${EGGNOG_DIR}" \
# --data_dir "${EGGNOG_DB}" \
# --cpu 2 \
# --override \

#echo -e "primera fase de eggnog finalizada con éxito. Comenzando la segunda fase... "

# Segunda fase de eggnog

echo -e "procesando los hits crudos y generando tablas las tablas finales de anotación... "

emapper.py \
  --annotate_hits_table "${EGGNOG_DIR}/eggnog_annotation.emapper.seed_orthologs" \
  --data_dir "${EGGNOG_DB}" \
  --output eggnog_annotation \
  --output_dir "${EGGNOG_DIR}" \
  --cpu 2 \
  --override


echo -e "Transformación de alineamientos crudos en anotaciones funcionales completas realizada."

