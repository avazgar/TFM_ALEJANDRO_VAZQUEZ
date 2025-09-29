#!/bin/bash

# === CONFIGURACIÓN ===
source /home/alejandro.vazquez/Desktop/A.VAZQUEZ_TFM_2025_MBIF/TRANSCRIPTOMICS/CODE/config.env

SILVA_INDEX="${INDEX_DIR}/SILVA_138.2_SSU_INDEX"
OUTPUT_DIR="${CHECK_rRNA_DIR}/BOWTIE2/rRNA_CHECK"
mkdir -p "${OUTPUT_DIR}"

echo "Iniciando chequeo de rRNA contra SILVA SSU (16S/18S)..."

# Crear archivo resumen
echo -e "Sample\tTotal_Reads\tAligned_Reads\tOverall_Alignment(%)" > "${OUTPUT_DIR}/rRNA_summary.tsv"

for R1 in "${TRIM_DIR}"/*_1_final.fq.gz; do
    base=$(basename "${R1}" _1_final.fq.gz)
    R2="${TRIM_DIR}/${base}_2_final.fq.gz"

    echo "Procesando muestra: ${base}"

    # Submuestreo 10% para agilizar (los RNAr se encuentran distribuidos de forma homogénea, los resultados obtenidos en este estudio serían fiables a pesar del subsampling)
    seqtk sample -s100 "${R1}" 0.1 > "${OUTPUT_DIR}/${base}_sub_R1.fq"
    seqtk sample -s100 "${R2}" 0.1 > "${OUTPUT_DIR}/${base}_sub_R2.fq"

    # Mapeo con Bowtie2
    bowtie2 -x "${SILVA_INDEX}" \
            -1 "${OUTPUT_DIR}/${base}_sub_R1.fq" \
            -2 "${OUTPUT_DIR}/${base}_sub_R2.fq" \
            --very-sensitive -p "${THREADS}" \
            -S "${OUTPUT_DIR}/${base}_rRNA.sam" \
            2> "${OUTPUT_DIR}/${base}_bowtie2.log"

    # Extraer estadísticas del log
    total=$(grep "reads;" "${OUTPUT_DIR}/${base}_bowtie2.log" | awk '{print $1}')
    aligned=$(grep "overall alignment rate" "${OUTPUT_DIR}/${base}_bowtie2.log" | awk '{print $1}')
    overall=$(grep "overall alignment rate" "${OUTPUT_DIR}/${base}_bowtie2.log" | awk '{print $1}')

    echo -e "${base}\t${total}\t${aligned}\t${overall}" >> "${OUTPUT_DIR}/rRNA_summary.tsv"

    # Limpieza opcional
    rm "${OUTPUT_DIR}/${base}_rRNA.sam"
done

echo "Chequeo de rRNA finalizado. Resumen en: ${OUTPUT_DIR}/rRNA_summary.tsv"


