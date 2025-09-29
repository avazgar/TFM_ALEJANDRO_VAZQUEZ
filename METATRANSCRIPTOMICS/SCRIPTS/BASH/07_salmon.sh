#!/bin/bash

# === CONFIGURACIÓN ===
source /home/alejandro.vazquez/Desktop/A.VAZQUEZ_TFM_2025_MBIF/TRANSCRIPTOMICS/CODE/config.env

PRODIGAL_DIR="${PRODIGAL_DIR}"
SALMON_DIR="${SALMON_DIR}"
mkdir -p "${SALMON_DIR}"

echo "Iniciando cuantificación de abundancias con Salmon..."

#for SAMPLE_DIR in "${PRODIGAL_DIR}"/*; do
#    base=$(basename "${SAMPLE_DIR}")
#    GENES="${SAMPLE_DIR}/genes.fna"

 #   if [ -f "${GENES}" ]; then
 #       echo "Indexando genes de ${base}..."
 #       INDEX="${SAMPLE_DIR}/salmon_index"
 #       salmon index -t "${GENES}" -i "${INDEX}" --threads "${THREADS}"

  #      echo "Cuantificando abundancia para ${base}..."
  #      R1="${TRIM_DIR}/${base}_1_final.fq.gz"
  #      R2="${TRIM_DIR}/${base}_2_final.fq.gz"
  #      OUT_DIR="${SALMON_DIR}/${base}"
  #      mkdir -p "${OUT_DIR}"

  #      salmon quant -i "${INDEX}" -l A \
  #          -1 "${R1}" -2 "${R2}" \
  #          -o "${OUT_DIR}" \
  #          --validateMappings \
  #          --threads "${THREADS}"
  #  else
  #      echo "genes.fna no encontrado para ${base}, saltando..."
  #  fi
#done

#echo "Generando matriz unificada de abundancias (TPM + counts)..."

# Crear tablas TPM y counts
#MATRIX_DIR="${SALMON_DIR}/matrices"
#mkdir -p "${MATRIX_DIR}"


salmon quantmerge \
    --quants $(find ${SALMON_DIR} -maxdepth 1 -type d ! -name matrices -not -path ${SALMON_DIR}) \
    --column tpm \
    --output ${SALMON_DIR}/matrices/salmon_TPM.tsv

salmon quantmerge \
    --quants $(find ${SALMON_DIR} -maxdepth 1 -type d ! -name matrices -not -path ${SALMON_DIR}) \
    --column numreads \
    --output ${SALMON_DIR}/matrices/salmon_counts.tsv


echo "Cuantificación finalizada. Matrices en: ${MATRIX_DIR}"

