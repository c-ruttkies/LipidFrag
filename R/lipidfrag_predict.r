source("lipidfrag_functions.r")

# reads MetFrag generated CSV result file from LipidMaps query and predicts lipid class
# 'Identifier' column must contain LipidMaps ID
#
predict.lipidmaps.class <- function(metfrag.csv.file, lipid.class.models, sep = ",") {
  # read result file
  data <- read.csv(metfrag.csv.file, header = T, sep = sep)
  
  annotated.candidate.list <- get.fcps.candidate.list(data, lipid.class.models)
  
  annotated.candidate.list[, "FCP"]
}