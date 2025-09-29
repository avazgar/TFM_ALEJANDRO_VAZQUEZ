#!/bin/bash

# THIS SCRIPT AIMS TO CARRY OUT A QUALITY CONTROL OF THE PREVIOUSLY DOWNLOAD AND CONVERTED FASTQ FILES WITH FASTQC TOOL

## LOAD CONFIG.ENV

source /home/alejandro.vazquez/Desktop/A.VAZQUEZ_TFM_2025_MBIF/TRANSCRIPTOMICS/CODE/config.env

## CREATE OUTPUT DIRECTORY (IF IT DOESN'T EXIST)

mkdir -p "${QC_DIR}/FASTQC_RAW"

## APPLY FASTQC OVER THE ORIGINAL FASTQ

echo "Ejecutando fastQC en las lecturas crudas..."

fastqc "${RAW_DIR}"/*.fastq.gz -o "${QC_DIR}/FASTQC_RAW" -t "$THREADS"

## APPLY MULTIQC TO GET A QUALITY REPORT

echo "Ejecutando MultiQC..."

multiqc "${QC_DIR}/FASTQC_RAW" -o "${QC_DIR}/FASTQC_RAW"

echo "El QC inicial se ha completado con éxito. Los resultados están disponibles en ${QC_DIR}/FASTQC_RAW"
