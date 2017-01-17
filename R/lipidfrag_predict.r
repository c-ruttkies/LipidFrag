source("lipidfrag_functions.r")

# reads MetFrag generated CSV result file from LipidMaps query and predicts lipid class with highest FCP and
# MetFrag score

predict.lipidmaps.class <- function(metfrag.csv.file, lipid.class.models, sep = ",") {
  # read result file
  annotated.candidate.list <- annotate.lipidmaps.class(metfrag.csv.file, lipid.class.models, sep)
  # get FCP values
  fcps <- as.numeric(as.character(annotated.candidate.list[,"FCP"]))
  fcps[is.na(fcps)] <- 0.0
  annotated.candidate.list[,"FCP"] <- fcps
  # get entries with maximal FCP values
  max.indexes <- which(fcps == max(fcps))
  max.fragmenterscore.index <- which.max(annotated.candidate.list[max.indexes, "FragmenterScore"])
  
  return(annotated.candidate.list[,c("Identifier", "FCP", "LipidMapsClass")][max.indexes,][max.fragmenterscore.index,])
}



# reads MetFrag generated CSV result file from LipidMaps query and predicts lipid class for all entries in the candidate list
# 'Identifier' column must contain LipidMaps ID
#
annotate.lipidmaps.class <- function(metfrag.csv.file, lipid.class.models, sep = ",") {
  # read result file
  data <- read.csv(metfrag.csv.file, header = T, sep = sep, comment.char="")
  
  annotated.candidate.list <- get.fcps.candidate.list(data, lipid.class.models)
  
  return(annotated.candidate.list)
}


# calculate the posterior foreground class probability (FCP) for a given MetFrag result list

get.fcps.candidate.list <- function(candidate.list, lipid.class.models) {
  class.names.model <- names(lipid.class.models)
  fcps <- t(sapply(1:dim(candidate.list)[1], function(row.index) {
    row <- candidate.list[row.index,]
    identifier <- as.character(row[1, "Identifier"])
    score <- as.numeric(row[1, "FragmenterScore"])
    matching.model.index <- which(sapply(class.names.model, function(x) {any(startsWith(identifier, unlist(strsplit(x,"_"))))}))
    if(length(matching.model.index) == 0) { return(c(0, identifier)) }
    if(length(matching.model.index) == 0) stop("More than one matching model for Identifier ", identifier, ". ", "Which shall I use?", sep="")
    return(c(get.posterior.foreground(score, lipid.class.models[[matching.model.index]]), class.names.model[matching.model.index]))
  }))
  annotated.candidate.list <- cbind(candidate.list, fcps)
  colnames(annotated.candidate.list) <- c(colnames(candidate.list), c("FCP", "LipidMapsClass"))
  return(annotated.candidate.list)
} 

