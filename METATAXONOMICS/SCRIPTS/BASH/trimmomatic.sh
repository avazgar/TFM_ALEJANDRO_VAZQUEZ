#!/bin/bash

# Crear las variables

INPUT_DIRECTORY="../rawfastq/"
THREADS=2
OUTPUT_DIRECTORY="../trimmomatic_output/"

# Crear el directorio

mkdir -p "${OUTPUT_DIRECTORY}"

# Comprobar si trimmomatic fue instalado y está disponible

command -v trimmomatic >/dev/null 2>&1 || { echo >&2 "Trimmomatic no está instalado. Abortando."; exit 1; }


for sample in "${INPUT_DIRECTORY}"/*_R1.fastq.gz; do
	base=$(basename "$sample" _R1.fastq.gz)
	rev="${INPUT_DIRECTORY}/${base}_R2.fastq.gz"


    	if [[ ! -f "$rev" ]]; then
        echo "No se encontró el archivo R2 para ${base}. Saltando muestra."
        continue
    	fi

	echo "Procesando muestra: ${base}"

	trimmomatic PE -threads "${THREADS}" -phred33 \
		"${INPUT_DIRECTORY}/${base}_R1.fastq.gz" "${rev}" \
		"${OUTPUT_DIRECTORY}/${base}_R1_paired.fastq.gz" "${OUTPUT_DIRECTORY}/${base}_R1_unpaired.fastq.gz" \
		"${OUTPUT_DIRECTORY}/${base}_R2_paired.fastq.gz" "${OUTPUT_DIRECTORY}/${base}_R2_unpaired.fastq.gz" \
		LEADING:3 TRAILING:3 SLIDINGWINDOW:50:20 MINLEN:215
done

echo "Trimmomatic finalizado"

