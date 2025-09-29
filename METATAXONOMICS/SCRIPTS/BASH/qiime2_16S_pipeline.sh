#!/bin/bash

# ****************************************************************************** #
# Pipeline completo para metataxonómica utilizando qiime2
# Asignación taxonómica en el gen marcador 16S rRNA (V4-V5, 515f-926r)
# Uso de clasificador Bayesiano adaptado para nuestra región de interés
# Uso de SILVA v.138 SSURef NR99 full-length region sequences
# Pipeline adaptado a la estructura de directorios presente en 16S_analysis
# Lanzar este script desde el directorio "scripts"

# =============================================================================== #

# **Variables generales**

METADATA="../input/metadata.tsv"
MANIFEST="../input/manifest.csv"
CLASSIFIER="../input/classifier/classifier_16SV4-V5.qza"
THREADS=2
FORWARD_PRIMER="GTGYCAGCMGCCGCGGTAA"
REVERSE_PRIMER="CCGYCAATTYMTTTRAGTTT"
MAX_EE_F=2
MAX_EE_R=2
TRUNC_Q=2
MIN_OVERLAP=10


# =========================================================== #
# *******************Importación de datos********************


#echo "Importando los archivos FastQ..."

#if ! qiime tools import \
# --type SampleData[PairedEndSequencesWithQuality] \
# --input-format PairedEndFastqManifestPhred33 \
# --input-path "${MANIFEST}" \
# --output-path ../qiime2/00_import/demux-seqs_16SV4-V5.qza; then
#	echo "No es posible importar los archivos FastQ, compruebe que los archivos "manifest.csv" y "metadata.tsv" han sido generados correctamente"
#	exit 1
#fi

#echo "¡Archivos importados correctamente!"

#if ! qiime demux summarize \
# --i-data ../qiime2/00_import/demux-seqs_16SV4-V5.qza \
# --o-visualization ../qiime2/07_visualizations/demux-seqs_16SV4-V5.qzv; then

#	echo "No se encuentran los archivos necesarios para generar el reporte interactivo, compruebe las secuencias demultiplexadas .qza"
#	exit 1
#fi

#echo "generando el resumen de calidad inicial de las secuencias, abriendo reporte interactivo para su visualización..."

#qiime tools view ../qiime2/07_visualizations/demux-seqs_16SV4-V5.qzv

#echo "Identificando y eliminando las secuencias de los primers con Cutadapt..."

# =============================== 1 ================================ #
# **********************Eliminación de primers**********************

#if ! qiime cutadapt trim-paired \
#  --i-demultiplexed-sequences ../qiime2/00_import/demux-seqs_16SV4-V5.qza \
#  --p-front-f "${FORWARD_PRIMER}" \
#  --p-front-r "${REVERSE_PRIMER}" \
#  --o-trimmed-sequences ../qiime2/00_import/demux-seqs_16SV4-V5_trimmed.qza \
#  --verbose \
#  --p-cores "${THREADS}"; then
#	echo "No se han podido filtrar los primers de las secuencias"
#	exit 1
#fi

#echo "Iniciando visualización tras el recorte..."

#if ! qiime demux summarize \
# --i-data ../qiime2/00_import/demux-seqs_16SV4-V5_trimmed.qza \
# --o-visualization ../qiime2/07_visualizations/demux-seqs_16SV4-V5_trimmed.qzv; then
#	echo "No se han podido visualizar los resultados"
#	exit 1
#fi

#qiime tools view ../qiime2/07_visualizations/demux-seqs_16SV4-V5_trimmed.qzv

# ============================= 1.5 ============================== #
# **********************Denoising con DADA2**********************


echo "Introduce los valores de truncamiento (p-trim-left & p-trunc-len) en función de la calidad observada:"


while [[ -z "${TRIM_LEFT_F}" ]]; do
	read -p "trim-left-f: " TRIM_LEFT_F
done

while [[ -z "${TRIM_LEFT_R}" ]]; do
	read -p "trim-left-r: " TRIM_LEFT_R
done

while [[ -z "${TRUNC_LEN_F}" ]]; do
	read -p "trunc-len-f: " TRUNC_LEN_F
done

while [[ -z "${TRUNC_LEN_R}" ]]; do
	read -p "trunc-len-r: " TRUNC_LEN_R
done


echo "Evaluando...Este proceso puede tomar unos minutos."


if ! qiime dada2 denoise-paired \
 --i-demultiplexed-seqs ../qiime2/00_import/demux-seqs_16SV4-V5_trimmed.qza \
 --p-trim-left-f "${TRIM_LEFT_F}" \
 --p-trim-left-r "${TRIM_LEFT_R}" \
 --p-trunc-len-f "${TRUNC_LEN_F}" \
 --p-trunc-len-r "${TRUNC_LEN_R}" \
 --p-max-ee-f "${MAX_EE_F}" \
 --p-max-ee-r "${MAX_EE_R}" \
 --p-trunc-q "${TRUNC_Q}" \
 --p-min-overlap "${MIN_OVERLAP}" \
 --p-pooling-method independent \
 --p-chimera-method consensus \
 --p-min-fold-parent-over-abundance 1.0 \
 --p-n-threads "${THREADS}" \
 --o-representative-sequences ../qiime2/02_dada2/rep-seqs_16SV4-V5_dada2_tq2-fpoa1.qza \
 --o-table ../qiime2/02_dada2/table_16SV4-V5_dada2_tq2-fpoa1.qza \
 --o-denoising-stats ../qiime2/01_quality/stats_16SV4-V5_dada2_tq2-fpoa1.qza; then
	echo "Error al realizar el denoise de "demux-seqs_16SV4-V5_OPT.qza"..."
	exit 1
fi

if ! qiime metadata tabulate \
 --m-input-file ../qiime2/01_quality/stats_16SV4-V5_dada2_tq2-fpoa1.qza \
 --o-visualization ../qiime2/07_visualizations/stats-16SV4-V5-dada2_tq2-fpoa1.qzv; then
	echo "Error al generar las estadísticas del proceso de denoising, comprueba la estructura de directorios."
	exit 1
fi

if ! qiime feature-table summarize \
  --i-table ../qiime2/02_dada2/table_16SV4-V5_dada2_tq2-fpoa1.qza \
  --o-visualization ../qiime2/07_visualizations/table_16SV4-V5_dada2_tq2-fpoa1.qzv \
  --m-sample-metadata-file "${METADATA}"; then
	echo "Error al generar la tabla de ASVs, comprueba la estructura de directorios."
	exit 1
fi

if ! qiime feature-table tabulate-seqs \
  --i-data ../qiime2/02_dada2/rep-seqs_16SV4-V5_dada2_tq2-fpoa1.qza \
  --o-visualization ../qiime2/07_visualizations/rep-seqs_16SV4-V5_dada2_tq2-fpoa1.qzv; then
	echo "Error al generar las secuencias representativas, comprueba la estructura de directorios."
	exit 1
fi

echo "¡El proceso de denoise ha finalizado con éxito! Dirígete al directorio! Mostrando las estadísticas del proceso..."

qiime tools view ../qiime2/07_visualizations/stats-16SV4-V5-dada2_tq2-fpoa1.qzv

read -p "En caso de que las estadísticas sean las esperadas, Pulsa [ENTER] para continuar con el siguiente paso del pipeline, de lo contrario, interrumpe con Ctrl + C: "


# ============================== 2 ============================== #
# *********************Asignación taxonómica*********************
# *******************Clasificador customizado********************

echo "Iniciando la asignación taxonómica con clasificador bayesiano adaptado a nuestra región de estudio V4-V5 (515f-926r)..."

if ! qiime feature-classifier classify-sklearn \
  --i-classifier "${CLASSIFIER}" \
  --i-reads ../qiime2/02_dada2/rep-seqs_16SV4-V5_dada2_OPT.qza \
  --p-n-jobs "${THREADS}" \
  --o-classification ../qiime2/03_taxonomy/taxonomy_16SV4-V5_OPT.qza; then
	echo "Asignación taxonómica incompleta, revisa el archivo de referencia y el clasificador"
	exit 1
fi

if ! qiime metadata tabulate \
 --m-input-file ../qiime2/03_taxonomy/taxonomy_16SV4-V5_OPT.qza \
 --o-visualization ../qiime2/03_taxonomy/taxonomy_16SV4-V5_OPT.qzv; then
	echo "Asignación taxonómica incompleta, no se ha podido generar el archivo de taxonomía .qzv"
	exit 1
fi

echo "¡Asignación taxonómica realizada con éxito!"


# ================================= 3 ================================= #
# ************************Visualizar abundancias************************
# *******************Antes de limpieza mitocondrial*********************


echo "Colapsando abundancias a nivel 4 (orden) antes del filtrado de mitocondrias y cloroplastos..."
if ! qiime taxa collapse \
  --i-table ../qiime2/02_dada2/table_16SV4-V5_dada2.qza \
  --i-taxonomy ../qiime2/03_taxonomy/taxonomy_16SV4-V5.qza \
  --p-level 4 \
  --o-collapsed-table ../qiime2/02_dada2/table_16SV4-V5_dada2_order.qza; then
	echo "No se ha podido colapsar la tabla de abundancias a nivel de orden"
	exit 1
fi

echo "Generando resumen interactivo de la tabla de ASVs..."

if ! qiime feature-table summarize \
  --i-table ../qiime2/02_dada2/table_16SV4-V5_dada2_order.qza \
  --o-visualization ../qiime2/07_visualizations/table_16SV4-V5_dada2_order.qzv \
  --m-sample-metadata-file "${METADATA}"; then
	echo "Error al ejecutar la acción, comprueba la estructura de directorios."
	exit 1
fi

echo "Abriendo el reporte interactivo de abundancias a nivel de orden..."

qiime tools view ../qiime2/07_visualizations/table_16SV4-V5_dada2_order.qzv

read -p "Pulsa [ENTER] para continuar con el filtrado de mitocondrias: "

# ============================= 4 ============================== #
# ****************Eliminación de DNA Mitocondrial****************

if ! qiime taxa filter-table \
  --i-table ../qiime2/02_dada2/table_16SV4-V5_dada2.qza  \
  --i-taxonomy ../qiime2/03_taxonomy/taxonomy_16SV4-V5.qza \
  --p-exclude mitochondria \
  --o-filtered-table ../qiime2/04_filtering/table_no_mitochondria.qza; then
	echo "Error al generar la tabla de taxonomías filtrada."
	exit 1
fi

if ! qiime taxa filter-seqs \
  --i-sequences ../qiime2/02_dada2/rep-seqs_16SV4-V5_dada2.qza \
  --i-taxonomy ../qiime2/03_taxonomy/taxonomy_16SV4-V5.qza \
  --p-exclude mitochondria \
  --o-filtered-sequences ../qiime2/04_filtering/rep-seqs_no_mitochondria.qza; then
	echo "Error al eliminar las secuencias asignables a mitocondrias de las secuencias representativas"
	exit 1
fi
echo "¡Secuencias de mitocondrias filtradas con éxito! Iniciando la asignación taxonómica con el clasificador bayesiano sobre las secuencias filtradas..."

if ! qiime feature-classifier classify-sklearn \
  --i-classifier "${CLASSIFIER}" \
  --i-reads ../qiime2/04_filtering/rep-seqs_no_mitochondria.qza \
  --p-n-jobs "${THREADS}" \
  --o-classification ../qiime2/04_filtering/taxonomy_no_mitochondria_16SV4-V5.qza; then
	echo "Error al filtrar los IDs atribuibles a mitocondrias de los ASVs"
	exit 1
fi

# ============================= 3.5 ============================== #
# *******************Evaluación de cloroplastos*******************

read -p "¡Filtrado mitocondrial finalizado con éxito! Pulsa [ENTER] para continuar con la evaluación de cloroplastos presentes: "

echo "Evaluando si "o__Chloroplast" contiene entradas no compatibles con cianobacterias verdaderas..."

echo "Los archivos resultantes de esta fase se depositarán en un nuevo directorio...Creando nueva ruta para depositar archivos --> qiime2/03_taxonomy_chloroplast_check"

mkdir -p ../qiime2/03_taxonomy/chloroplast_check

# Exportar la taxonomía filtrada de mitocondrias

if ! qiime tools export \
  --input-path ../qiime2/04_filtering/taxonomy_no_mitochondria_16SV4-V5.qza \
  --output-path ../qiime2/03_taxonomy/chloroplast_check; then
	echo "No se ha podido exportar la taxonomía filtrada del paso anterior, el archivo .tsv no se ha generado"
	exit 1
fi

cd ../qiime2/03_taxonomy/chloroplast_check

# Filtrar las líneas con "o__Chloroplast; f__Chloroplast" del archivo resultante del filtrado de mitocondrias

echo "Evaluando si hay cloroplastos a excluir..."

awk -F'\t' '$2 ~ /o__Chloroplast/ && $2 ~ /f__Chloroplast/' taxonomy.tsv > all_chloroplast.tsv

# Eliminamos del archivo las excepciones al filtrado, es decir, cianobacterias asignadas incorrectamente a cloroplastos

grep -vE 's__uncultured_cyanobacterium|s__uncultured_marine|s__uncultured_bacterium' all_chloroplast.tsv > chloroplast_to_filter.tsv

# Extraemos los IDs de los features que queremos eliminar (cloroplastos)

cut -f 1 chloroplast_to_filter.tsv > ids_chloroplast.txt
awk 'BEGIN { print "FeatureID" } NF > 0 { print $1 }' ids_chloroplast.txt > ids_chloroplast_metadata.tsv


if [[ -s ../qiime2/03_taxonomy/chloroplast_check/ids_chloroplast_metadata.tsv ]] && \
   [[ $(wc -l < ../qiime2/03_taxonomy/chloroplast_check/ids_chloroplast_metadata.tsv) -gt 1 ]]; then

	echo "Cloroplastos detectados, aplicando filtro..."

	cd ../../../scripts/

	# Aplicar el filtrado con qiime2

	echo "Generando la tabla de frecuencias sin cloroplastos ni mitocondrias..."

	if ! qiime feature-table filter-features \
  	--i-table ../qiime2/04_filtering/table_no_mitochondria.qza \
  	--m-metadata-file ../qiime2/03_taxonomy/chloroplast_check/ids_chloroplast_metadata.tsv \
	--p-exclude-ids \
  	--o-filtered-table ../qiime2/04_filtering/table_no_mito_no_chloro.qza; then
		echo "Error al filtrar la tabla de ASVs"
		exit 1
	fi

	echo "¡tabla de frecuencias generada!"

	echo "Generando las secuencias representativas sin cloroplastos ni mitocondrias..."

	if ! qiime feature-table filter-seqs \
  	--i-data ../qiime2/04_filtering/rep-seqs_no_mitochondria.qza \
  	--m-metadata-file ../qiime2/03_taxonomy/chloroplast_check/ids_chloroplast_metadata.tsv \
  	--p-exclude-ids \
  	--o-filtered-data ../qiime2/04_filtering/rep-seqs_no_mito_no_chloro.qza; then
		echo "Error al filtrar las secuencias representativas"
		exit 1
	fi

	echo "¡Secuencias representativas generadas!"

	echo "Generando una la taxonomía final con los filtros aplicados: No cloroplastos + No mitocondrias..."

	if ! qiime feature-classifier classify-sklearn \
  	--i-classifier "${CLASSIFIER}" \
  	--i-reads ../qiime2/04_filtering/rep-seqs_no_mito_no_chloro.qza \
  	--p-n-jobs "${THREADS}" \
  	--o-classification ../qiime2/04_filtering/taxonomy_no_mito_no_chloro.qza; then
        	echo "Error al filtrar los IDs atribuibles a cloroplastos de los ASVs"
        	exit 1
	fi

	echo "¡Nueva taxonomía con los filtros de cloroplastos y mitocondrias generada con éxito!"

	# Variables de salida actualizadas

	TABLE_PATH="../qiime2/04_filtering/table_no_mito_no_chloro.qza"
	REPS_PATH="../qiime2/04_filtering/rep-seqs_no_mito_no_chloro.qza"
	TAXO_PATH="../qiime2/04_filtering/taxonomy_no_mito_no_chloro.qza"

	echo "¡Filtrado de secuencias potencialmente atribuibles a cloroplastos realizado con éxito!"

else
	echo "No se detectaron cloroplastos... Manteniendo los archivos post-mitocondriales para continuar el análisis..."

	#Actualizar variables de salida a los archivos filtrados de mitocondrias

	TABLE_PATH="../qiime2/04_filtering/table_no_mitochondria.qza"
        REPS_PATH="../qiime2/04_filtering/rep-seqs_no_mitochondria.qza"
        TAXO_PATH="../qiime2/04_filtering/taxonomy_no_mitochondria_16SV4-V5.qza"

fi

cd ../../../scripts/

# ================================= 4 ================================== #
# ************************* Arbol filogenetico *************************

echo "Construyendo el árbol filogenético..."

# Construir árbol filogenético

if ! qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences "${REPS_PATH}" \
  --o-alignment ../qiime2/05_phylogeny/aligned-rep-seqs.qza \
  --o-masked-alignment ../qiime2/05_phylogeny/masked-aligned-rep-seqs.qza \
  --o-tree ../qiime2/05_phylogeny/unrooted-tree.qza \
  --o-rooted-tree ../qiime2/05_phylogeny/rooted-tree.qza; then
	echo "No se ha podido construir el árbol filogenético, revisa los archivos necesarios."
	exit 1
fi

# Visualizar tabla para determinar profundidad de muestreo
qiime feature-table summarize \
  --i-table "${TABLE_PATH}" \
  --o-visualization ../qiime2/05_phylogeny/table_no_mito_no_chloro.qzv

qiime tools view ../qiime2/05_phylogeny/table_no_mito_no_chloro.qzv

read -p "Introduce manualmente el valor de "--p-min-frequency" en base a lo observado: " MIN_FREQUENCY

#Filtro por frecuencia mínima

qiime feature-table filter-samples \
  --i-table "${TABLE_PATH}" \
  --p-min-frequency "${MIN_FREQUENCY}" \
  --o-filtered-table ../qiime2/04_filtering/table_16SV4-V5_filtered.qza
qiime feature-table summarize \
  --i-table ../qiime2/04_filtering/table_16SV4-V5_filtered.qza \
  --o-visualization ../qiime2/04_filtering/table_16SV4-V5_filtered.qzv \
  --m-sample-metadata-file "${METADATA}"

# Actualizamos la variable TABLE_PATH

TABLE_PATH="../qiime2/04_filtering/table_16SV4-V5_filtered.qza"

# Visualizar la tabla filtrada para determinar --p-sampling-depth en el análisis de rarefacción

qiime tools view ../qiime2/04_filtering/table_16SV4-V5_filtered.qzv


# ================================= 5 ================================== #
# ***********************Diversidad y rarefacción***********************

# Ejecutar métricas básicas con rarefacción

read -p "Introduce manualmente el valor de "p-sampling-depth" en base a lo observado: " SAMPLING_DEPTH

qiime diversity alpha-rarefaction \
  --i-table "${TABLE_PATH}" \
  --i-phylogeny ../qiime2/05_phylogeny/rooted-tree.qza \
  --p-max-depth "${SAMPLING_DEPTH}" \
  --m-metadata-file "${METADATA}" \
  --o-visualization ../qiime2/06_diversity/alpha-rarefaction.qzv

echo "Abriendo el reporte interactivo de diversidad alfa..."

qiime tools view ../qiime2/06_diversity/alpha-rarefaction.qzv

# Eliminar el directorio core_metrics_results existente del análisis previo.

rm -rf ../qiime2/06_diversity/core_metrics_results

if ! qiime diversity core-metrics-phylogenetic \
  --i-phylogeny ../qiime2/05_phylogeny/rooted-tree.qza \
  --i-table "${TABLE_PATH}" \
  --p-sampling-depth "${SAMPLING_DEPTH}" \
  --m-metadata-file "${METADATA}" \
  --output-dir ../qiime2/06_diversity/core_metrics_results; then
    echo "Error al ejecutar core-metrics"
    exit 1
fi

# ================================= 6 ================================== #
# **********************Barplots de diversidad final**********************

if ! qiime taxa barplot \
  --i-table "${TABLE_PATH}" \
  --i-taxonomy "${TAXO_PATH}" \
  --m-metadata-file "${METADATA}" \
  --o-visualization ../qiime2/07_visualizations/taxonomy-barplot.qzv; then
	echo "Error al generar los Barplots de diversidad"
	exit 1
fi

# ===================================== 7 ====================================== #
# ************************16S Qiime2 Pipeline Finalizado ************************

echo "Pipeline de QIIME2 para 16S completado con éxito, los resultados y artefactos intermedios del flujo están disponibles en los correspondientes directorios del proyecto"


# ===================================== 8 ====================================== #
# ***************************Exportar archivos Rstudio**************************



echo "Exportación de archivos necesarios para el análisis estadístico en Microeco..."

mkdir -p ../qiime2/08_export/microeco_package

# Exportar la tabla de frecuencias final

if ! qiime tools export \
  --input-path "${TABLE_PATH}" \
  --output-path ../qiime2/08_export/microeco_package/; then
	echo "No se ha podido exportar la tabla de frecuencias, comprueba el error reportado en la terminal"
	exit 1
fi

# Exportar la taxonomía filtrada final

if ! qiime tools export \
  --input-path "${TAXO_PATH}" \
  --output-path ../qiime2/08_export/microeco_package/; then
	echo "No se ha podido exportar la taxonomía, comprueba el error reportado en la terminal"
	exit 1
fi

# Exportar las secuencias representativas finales

if ! qiime tools export \
  --input-path "${REPS_PATH}" \
  --output-path ../qiime2/08_export/microeco_package/; then
	echo "No se ha podido exportar las secuencias representativas, comprueba el error reportado en la terminal"
	exit 1
fi

# Copiar metadata para compatibilidad
cp ../input/metadata.tsv ../qiime2/08_export/microeco_package/


echo "¡Archivos exportados correctamente! Revisa el directorio microeco_package"


