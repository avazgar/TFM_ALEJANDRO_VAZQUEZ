# --- Configuración de rutas desde variables de entorno ---
source /home/alejandro.vazquez/Desktop/A.VAZQUEZ_TFM_2025_MBIF/TRANSCRIPTOMICS/CODE/config.env

# --- Definición de variables específicas ---
SALMON_OUT="${SALMON_DIR}/SALMON_COMMON"
THREADS="${THREADS}"

# --- Crear carpetas necesarias ---
mkdir -p "${SALMON_OUT}/matrices"

# --- Paso 1: Fusionar todos los genes.fna (no proteins.faa) ---
echo "Fusionando genes predichos..."
cat "${PRODIGAL_DIR}"/*/genes.fna > "${SALMON_OUT}/all_genes.fna"

# --- Paso 2: Crear índice común ---
echo "Generando índice común con Salmon..."
salmon index -t "${SALMON_OUT}/all_genes.fna" -i "${SALMON_OUT}/salmon_index" -p ${THREADS}

# --- Paso 3: Cuantificar todas las muestras ---
echo "Cuantificando muestras..."
for sample_dir in "${PRODIGAL_DIR}"/*; do
    sample_id=$(basename "$sample_dir")
    fq1="${TRIM_DIR}/${sample_id}_1_final.fq.gz"
    fq2="${TRIM_DIR}/${sample_id}_2_final.fq.gz"
    outdir="${SALMON_OUT}/${sample_id}"

    echo "   • Cuantificando $sample_id..."
    salmon quant -i "${SALMON_OUT}/salmon_index" -l A \
        -1 "$fq1" -2 "$fq2" -p "$THREADS" \
        -o "$outdir" --validateMappings
done

# --- Paso 4: Generar matrices de abundancia ---
echo "Generando matrices..."
salmon quantmerge \
    --quants $(find ${SALMON_OUT} -mindepth 1 -maxdepth 1 -type d -name 'SRR*') \
    --column tpm \
    --output ${SALMON_OUT}/matrices/salmon_TPM.tsv

salmon quantmerge \
    --quants $(find ${SALMON_OUT} -mindepth 1 -maxdepth 1 -type d -name 'SRR*') \
    --column numreads \
    --output ${SALMON_OUT}/matrices/salmon_counts.tsv

echo "Cuantificación completada. Matrices en:"
echo "   - ${SALMON_OUT}/matrices/salmon_TPM.tsv"
echo "   - ${SALMON_OUT}/matrices/salmon_counts.tsv"

