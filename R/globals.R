# globals shared across markdown docs

library(Seurat)
library(readr)
library(tidyr)
library(plyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(Matrix)
library(matrixStats)
library(purrr)
library(R.utils)
library(viridis)
library(ComplexHeatmap)
library(doParallel)
library(itertools)
#### Paths ####

project_dir <- path.expand("~/Projects/zemans")
data_dir <- file.path(project_dir, "data")
results_dir <- file.path(project_dir, "results")
docs_dir <- file.path(project_dir, "docs")
db_dir <- file.path(project_dir, "dbases")

# vector of figure paths
figs_dir <-  file.path(results_dir, "Figures") %>%
  dir(pattern = "Figure_[1-4]$",
      include.dirs = TRUE,
      full.names = T)


##### Functions ####

#' When writing out excel workbooks using openxlsx::write.xlsx()
#' this function will set the class attributes for a column, which
#' enforces a column type in the resulting xlsx file. 
#' Useful for avoid gene names being clobbered to dates and 
#' setting scientific number formatting

set_xlsx_class <- function(df, col, xlsx_class){
  for(i in seq_along(col)){
    class(df[[col[i]]]) <- xlsx_class
  }
  df
}







