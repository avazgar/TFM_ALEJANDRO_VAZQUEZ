#!/bin/bash

# === CONFIGURACIÓN ===
source /home/alejandro.vazquez/Desktop/A.VAZQUEZ_TFM_2025_MBIF/TRANSCRIPTOMICS/CODE/config.env

mkdir -p "${TRIM_DIR}"
mkdir -p "${TRIM_DIR}/fastp_reports"

echo "Iniciando trimming de adaptadores y poly-G..."

for R1 in "${RAW_DIR}"/*_1_SUBSAMPLED.fastq.gz; do
    base=$(basename "${R1}" _1_SUBSAMPLED.fastq.gz)
    R2="${RAW_DIR}/${base}_2_SUBSAMPLED.fastq.gz"

    if [[ -f "${R2}" ]]; then
        echo "🔹 Procesando muestra: ${base}"

        # === 1. Trim Galore! ===
        trim_galore \
            --paired \
            --quality 20 \
            --length "${MIN_LENGTH}" \
            --stringency 3 \
            --cores "${THREADS}" \
            --output "${TRIM_DIR}" \
            "${R1}" "${R2}"

        # Archivos generados por Trim Galore!
        TG_R1="${TRIM_DIR}/${base}_1_SUBSAMPLED_val_1.fq.gz"
        TG_R2="${TRIM_DIR}/${base}_2_SUBSAMPLED_val_2.fq.gz"

        # === 2. fastp (eliminación poly-G) ===
        fastp \
            --in1 "${TG_R1}" --in2 "${TG_R2}" \
            --out1 "${TRIM_DIR}/${base}_1_final.fq.gz" \
            --out2 "${TRIM_DIR}/${base}_2_final.fq.gz" \
            --trim_poly_g \
            --detect_adapter_for_pe \
            --thread "${THREADS}" \
            --html "${TRIM_DIR}/fastp_reports/${base}_fastp.html" \
            --json "${TRIM_DIR}/fastp_reports/${base}_fastp.json"

    else
        echo "Error: No se encontró el par para ${R1}"
    fi
done

echo "Trimming completo: archivos finales en ${TRIM_DIR} (sufijo _final.fq.gz)"

