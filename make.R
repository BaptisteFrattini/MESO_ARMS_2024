#' MESO_ARMS_2024: A Research Compendium
#' 
#' @description 
#' A paragraph providing a full description of the project and describing each 
#' step of the workflow.
#' 
#' @author Baptiste Frattini \email{baptiste.frattini22@gmail.com}
#' 
#' @date 2024/01/23



## Install Dependencies (listed in DESCRIPTION) ----

devtools::install_deps(upgrade = "never")

### function to run when opening the project 4 the first time

install.packages("targets")
install.packages("here")
install.packages("renv")
renv::restore()

## Load Project Addins (R Functions and Packages) ----

## usefull command for package loading
renv::init()
renv::install()
renv::status()
renv::snapshot()


# make the pipeline
targets::tar_visnetwork()
targets::tar_make()
targets::tar_visnetwork()


## Global Variables ----

# You can list global variables here (or in a separate R script)


## Run Project ----

# List all R scripts in a sequential order and using the following form:
# source(here::here("analyses", "script_X.R"))
