###########################################################################################################################
# Script: 00_environment_configuration.R
# Description:
# 	- Define the working directory
#	- Define the project name
#	- Create the structure of directories
#	- Define the variables within the project
# Date: 23/09/2025
# Owner: Alejandro Vázquez García
# Context: Metataxonomic processing - ANALYSIS AND INTERPRETATION OF RESULTS - RSTUDIO

# Define the working directory (the current directory of the project)

setwd("C:/Users/aleja/Desktop/TFM_A_VAZQUEZ/METATAXONOMICS_TFM")

getwd()  # To clarify which is the current working directory and create the aforementioned structure from it.

# Ruta al directorio del proyecto

nombre_proyecto <- "FONS_CLIMATIC"
ruta_proyecto <- file.path(getwd())

# Create the structure of directories within the project

dir.create(file.path(ruta_proyecto,"scripts"))
dir.create(file.path(ruta_proyecto,"data"))
dir.create(file.path(ruta_proyecto,"data","raw"))
dir.create(file.path(ruta_proyecto,"data","processed"))
dir.create(file.path(ruta_proyecto,"documents"))
dir.create(file.path(ruta_proyecto,"figures"))

# Create variables to each location within the project

raw_path <- file.path(ruta_proyecto,"data","raw")
processed_path <- file.path(ruta_proyecto,"data","processed")
scripts_path <- file.path(ruta_proyecto,"scripts")
figures_path <- file.path(ruta_proyecto,"figures")
data_path <- file.path(ruta_proyecto,"data")
documents_path <- file.path(ruta_proyecto,"documents")

# Ensure being on the base directory:

setwd(file.path(ruta_proyecto))
