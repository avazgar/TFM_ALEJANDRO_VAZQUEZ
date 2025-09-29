#SCRIPT DISEÑADO PARA GENERAR EL ARCHIVO MANIFEST.CSV CON EL QUE TRABAJAREMOS EN QIIME2

##Delimitamos los directorios donde se encuentran los archivos fastq, el archivo de salida y el log.file:

FASTQ_DIR="../rawfastq"
MANIFEST_FILE="../input/manifest.csv"
LOG_FILE="../logs/create_manifest.log"

##Generamos el encabezado de MANIFEST_FILE

 echo "sample-id,absolute-filepath,direction" > "${MANIFEST_FILE}" 
 
 > "${LOG_FILE}"

##Recorremos los archivos fastq forward y reverse y los vamos incorporando siguiendo el orden del encabezado generado previamente, ocultaremos el stderr:

 for forward_file in "${FASTQ_DIR}"/*_R1.fastq.gz; do

	###Verificamos que es un archivo regular:

	if [[ -f "${forward_file}" ]]; then

		###Extraemos el nombre de la muestra, situado antes de _1.fastq.gz:

		filename=$(basename "${forward_file}")
		sample_id="${filename%%_R1.fastq.gz}"

		###Obtenemos las rutas absolutas de cada archivo:

		abs_forward=$(realpath "${forward_file}")
		reverse_file="${forward_file/_R1.fastq.gz/_R2.fastq.gz}"
		abs_reverse=$(realpath "${reverse_file}")

			###Comprobamos la existencia del archivo reverse:

			if [[ -f "${reverse_file}" ]]; then

			####Si el archivo reverse existe escribimos ambas líneas en el manifiesto separadas por comas:

			echo "${sample_id},${abs_forward},forward" >> "${MANIFEST_FILE}"
			echo "${sample_id},${abs_reverse},reverse" >> "${MANIFEST_FILE}"

			####Incorporamos cambios al log:

			echo "Añadida la muestra ${sample_id} al manifest.csv" | tee -a "${LOG_FILE}"

			else

			echo "Muestra ${sample_id} omitida: archivo reverse no encontrado" | tee -a "${LOG_FILE}"

			fi

 	fi

done

echo "Archivo Manifest generado: ${MANIFEST_FILE}" | tee -a "${LOG_FILE}"
