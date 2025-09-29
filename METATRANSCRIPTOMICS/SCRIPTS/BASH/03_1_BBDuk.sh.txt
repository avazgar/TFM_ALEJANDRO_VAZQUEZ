#!/bin/bash

# Cargar configuración
source /home/alejandro.vazquez/Desktop/A.VAZQUEZ_TFM_2025_MBIF/TRANSCRIPTOMICS/CODE/config.env

# Base de datos rRNA (SILVA fast version)
RIBO_DB="${BASE_DIR}/DATA/PROCESSED/SILVA_DATABASE/smr_v4.3_fast_db.fasta"

# Crear carpeta de salida
OUT_DIR="${BASE_DIR}/RESULTS/BBDUK"
mkdir -p "${OUT_DIR}"

# Iterar sobre archivos recortados
for R1 in "${TRIM_DIR}"/*1_subsampled.fq.gz ; do
    BASENAME=$(basename "${R1}" _1_subsampled.fq.gz )
    R2="${TRIM_DIR}/${BASENAME}_2_subsampled.fq.gz "

    echo "Procesando muestra ${BASENAME} con BBDuk..."

    bbduk.sh \
        in1="${R1}" \
        in2="${R2}" \
        out1="${OUT_DIR}/${BASENAME}_no_rRNA_R1.fq.gz" \
        out2="${OUT_DIR}/${BASENAME}_no_rRNA_R2.fq.gz" \
        ref="${RIBO_DB}" \
        k=25 \
        hdist=1 \
        stats="${OUT_DIR}/${BASENAME}_bbduk_stats.txt" \
        -Xmx5500m

    echo "Filtrado completado para ${BASENAME}"
done


