#!/bin/bash

# === CONFIGURACIÓN ===
source /home/alejandro.vazquez/Desktop/A.VAZQUEZ_TFM_2025_MBIF/TRANSCRIPTOMICS/CODE/config.env

ASSEMBLY_DIR="${ASSEMBLY_DIR}"
mkdir -p "${ASSEMBLY_DIR}"

echo "Iniciando ensamblaje de novo con MEGAHIT (modo secuencial)..."

for R1 in "${TRIM_DIR}"/*_1_final.fq.gz; do
    base=$(basename "${R1}" _1_final.fq.gz)
    R2="${TRIM_DIR}/${base}_2_final.fq.gz"
    SAMPLE_OUT="${ASSEMBLY_DIR}/${base}"

    echo "Ensamblando muestra: ${base}"

    # Si existe la carpeta de salida, eliminarla
    if [ -d "${SAMPLE_OUT}" ]; then
        echo "Carpeta de salida existente, eliminándola..."
        rm -rf "${SAMPLE_OUT}"
    fi

    # Ejecutar MEGAHIT (sin --tmp-dir)
    megahit \
        -1 "${R1}" -2 "${R2}" \
        --min-contig-len 300 \
        --k-list 21,33,55,77,99,127 \
        --num-cpu-threads "${THREADS}" \
        --out-dir "${SAMPLE_OUT}" \
        --memory 0.9

    # Verificar si el ensamblaje fue exitoso
    if [ -f "${SAMPLE_OUT}/final.contigs.fa" ]; then
        if command -v seqkit &> /dev/null; then
            seqkit stats "${SAMPLE_OUT}/final.contigs.fa" > "${SAMPLE_OUT}/assembly_stats.tsv"
        fi
        echo "Ensamblaje completado: ${SAMPLE_OUT}/final.contigs.fa"
    else
        echo "Error: No se generó el archivo final.contigs.fa para ${base}"
    fi
done

echo "Todos los ensamblajes finalizados. Resultados en: ${ASSEMBLY_DIR}"
