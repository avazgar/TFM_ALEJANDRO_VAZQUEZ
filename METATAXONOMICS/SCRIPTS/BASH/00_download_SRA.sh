#!/bin/bash

# THIS .sh FILE AIMS TO AUTOMATICALLY DOWNLOAD THE REQUIRED SRA FILES FROM NCBI REPOSITORY AND THEIR CONVERSION TO FASTQ FILES

## FILE CONTAINING THE LIST OF SRR FILES

SRR_LIST="/home/alejandro.vazquez/Desktop/A.VAZQUEZ_TFM_2025_MBIF/input/sra_list.txt" 
RAW_DIR="/home/alejandro.vazquez/Desktop/A.VAZQUEZ_TFM_2025_MBIF/rawfastq"
## CREATION OF THE DIRECTORY

mkdir -p "${RAW_DIR}"

echo "Inicio de la descarga de los archivos SRR desde ${SRR_LIST}"

while read -r SRR; do
	echo "Procesando ${SRR}..."

	# Download with prefetch
	prefetch "$SRR" --output-directory "${RAW_DIR}"

	# Conversion to FASTQ
	fastq-dump "$RAW_DIR/$SRR" \
	  --split-files \
	  --gzip \
	  --outdir "$RAW_DIR"

	echo "¡{$SRR} completado!"


done < "$SRR_LIST"

echo "¡Todos los SRR fueron descargados y convertidos a FASTQ correctamente!"
