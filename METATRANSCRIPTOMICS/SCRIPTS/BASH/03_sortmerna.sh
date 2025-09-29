#!/bin/bash

# Cargar configuración
source /home/alejandro.vazquez/Desktop/A.VAZQUEZ_TFM_2025_MBIF/TRANSCRIPTOMICS/CODE/config.env

echo "Iniciando el filtrado de rRNA con SortMeRNA..."

# Ruta a la base de datos default de rRNA (ajusta si está en otra carpeta)
REF_DB="${BASE_DIR}/DATA/PROCESSED/sortmerna_DOCS/sortmerna_DB/smr_v4.3_fast_db.fasta"

# Comprobar si el índice ya existe
if [ ! -f "${REF_INDEX}.idx" ]; then
    echo "Índices no encontrados. Generando con indexdb_rna..."
    indexdb_rna --ref "${REF_DB}" --ref "${REF_INDEX}"
    echo "Índices generados correctamente."
else
    echo "Índices ya existen. Saltando paso de indexación."
fi

# Crear directorio de salida si no existe
SORTMERNA_OUT="/opt/SORTMERNA_TMP"
mkdir -p "${SORTMERNA_OUT}"

# Iterar sobre archivos paired-end recortados
for R1 in "${TRIM_DIR}"/*_1_val_1.fq.gz ; do
    BASENAME=$(basename "${R1}" _1_val_1.fq.gz)
    R2="${TRIM_DIR}/${BASENAME}_2_val_2.fq.gz"

    echo "Procesando muestra ${BASENAME}..."

    sortmerna \
        --ref "${REF_DB}" \
        --reads "${R1}" \
        --reads "${R2}" \
        --paired_in \
        --fastx \
        --workdir "${SORTMERNA_OUT}/${BASENAME}/workdir" \
        --other "${SORTMERNA_OUT}/${BASENAME}/${BASENAME}_no_rRNA" \
        --aligned "${SORTMERNA_OUT}/${BASENAME}/${BASENAME}_rRNA" \
        --threads 2

    # Contar las secuencias rRNA eliminadas
    RRNA_FILE="${SORTMERNA_OUT}/${BASENAME}/${BASENAME}_rRNA.fq"
    if [[ -f "${RRNA_FILE}" ]]; then
        RRNA_COUNT=$(($(wc -l < "$RRNA_FILE") / 4))
    else
        RRNA_COUNT=0
    fi

    echo "${BASENAME}: se eliminaron ${RRNA_COUNT} secuencias rRNA."
done

