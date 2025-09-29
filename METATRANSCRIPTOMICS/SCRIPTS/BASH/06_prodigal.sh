#!/bin/bash

# === CONFIGURACIÓN ===
source /home/alejandro.vazquez/Desktop/A.VAZQUEZ_TFM_2025_MBIF/TRANSCRIPTOMICS/CODE/config.env

ASSEMBLY_DIR="${ASSEMBLY_DIR}"
PRODIGAL_DIR="${PRODIGAL_DIR}"
mkdir -p "${PRODIGAL_DIR}"

echo "Generando resumen global de ensamblajes..."

# Archivo de resumen
SUMMARY_FILE="${ASSEMBLY_DIR}/assembly_summary.tsv"
echo -e "Sample\tNum_Contigs\tTotal_bp\tMin_bp\tAvg_bp\tMax_bp\tN50" > "${SUMMARY_FILE}"

for SAMPLE_DIR in "${ASSEMBLY_DIR}"/*; do
    base=$(basename "${SAMPLE_DIR}")
    CONTIGS="${SAMPLE_DIR}/final.contigs.fa"

    if [ -f "${CONTIGS}" ]; then
        # Estadísticas con seqkit
        stats=$(seqkit stats -T "${CONTIGS}" | tail -n1)
        # Añadir estadísticas básicas
        echo -e "${base}\t${stats}" >> "${SUMMARY_FILE}"

        # Calcular N50 con seqkit fx2tab (opcional)
        N50=$(seqkit fx2tab -nl "${CONTIGS}" | awk '{print $2}' | sort -nr | \
              awk 'BEGIN{sum=0;total=0}{len[NR]=$1;sum+=len[NR]}END{for(i=1;i<=NR;i++){total+=len[i]; if(total>=sum/2){print len[i]; exit}}}')
        # Añadir N50 al final de la línea
        sed -i "s/^${base}.*/&\t${N50}/" "${SUMMARY_FILE}"
    fi
done

echo "Resumen guardado en: ${SUMMARY_FILE}"

echo "Ejecutando Prodigal para cada ensamblaje..."

for SAMPLE_DIR in "${ASSEMBLY_DIR}"/*; do
    base=$(basename "${SAMPLE_DIR}")
    CONTIGS="${SAMPLE_DIR}/final.contigs.fa"

    if [ -f "${CONTIGS}" ]; then
        echo "🔹 Procesando muestra: ${base}"

        OUT_DIR="${PRODIGAL_DIR}/${base}"
        mkdir -p "${OUT_DIR}"

        prodigal \
            -i "${CONTIGS}" \
            -a "${OUT_DIR}/proteins.faa" \
            -d "${OUT_DIR}/genes.fna" \
            -o "${OUT_DIR}/prodigal.gff" \
            -p meta

        echo "Prodigal completado para ${base}: ${OUT_DIR}"
    else
        echo "Contigs no encontrados para ${base}, saltando..."
    fi
done

echo "Todos los ensamblajes procesados con Prodigal. Resultados en: ${PRODIGAL_DIR}"

