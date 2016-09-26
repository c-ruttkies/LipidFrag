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
  
  return(annotated.candidate.list[,c("FCP", "LipidMapsClass")][max.indexes,][max.fragmenterscore.index,])
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
    matching.model.index <- which(sapply(class.names.model, function(x) {startsWith(identifier, x)}))
    if(length(matching.model.index) == 0) { return(c(NA, identifier)) }
    if(length(matching.model.index) == 0) stop("More than one matching model for Identifier ", identifier, ". ", "Which shall I use?", sep="")
    return(c(get.posterior.foreground(score, lipid.class.models[[matching.model.index]]), class.names.model[matching.model.index]))
  }))
  annotated.candidate.list <- cbind(candidate.list, fcps)
  colnames(annotated.candidate.list) <- c(colnames(candidate.list), c("FCP", "LipidMapsClass"))
  return(annotated.candidate.list)
} 

# get pessimistic rank of candidates matching given InChIKeys or LipidMaps identifiers
# Output: Rank, BetterCandidates, EqualCandidates, WorseCandidates, TotalCandiates, Identifier

get.pessimistic.rank <- function(candidate.list, inchikey1s = NULL, lmids = NULL) {
  if(is.null(inchikey1s) && is.null(lmids)) {
    cat("Error: InChIKey1 or LMID needed.\n")
    return(NA)
  }
  if(is.null(lmids)) {
    correct.scores <- candidate.list[as.character(candidate.list[,c("InChIKey1")]) %in% inchikey1s, c("Score", "Identifier")]
    wrong.scores <- candidate.list[!as.character(candidate.list[,c("InChIKey1")]) %in% inchikey1s, c("Score")]
  }
  else {
    inchikey1s <- as.character(candidate.list[as.character(candidate.list[,"Identifier"]) %in% lmids, "InChIKey1"])
    correct.scores <- candidate.list[as.character(candidate.list[,c("InChIKey1")]) %in% inchikey1s, c("Score", "Identifier")]
    wrong.scores <- candidate.list[!as.character(candidate.list[,c("InChIKey1")]) %in% inchikey1s, c("Score")]
  }
  if(dim(correct.scores)[1] == 0) {
    cat("Error: Identifier not found in candidate list.\n")
    return(NA)
  }
  max.correct.score <- correct.scores[which.max(correct.scores[,1]),]
  if(length(wrong.scores) == 0) {
    return.vals <- data.frame(Rank=1,BC=0,EC=0,WC=0,TC=1,ID=max.correct.score[,2])
    return(return.vals)
  }
  equal.cands <- length(which(max.correct.score[,1] == wrong.scores))
  better.cands <- length(which(max.correct.score[,1] < wrong.scores))
  worse.cands <- length(which(max.correct.score[,1] > wrong.scores))
  return.vals <- data.frame(Rank = better.cands + equal.cands + 1, BC = better.cands, 
                    EQ = equal.cands, WC = worse.cands, TC = length(wrong.scores) + 1, 
                    ID = max.correct.score[,2])
  return(return.vals)
}

# get pessimistic rank of FCP filtered candidates matching given InChIKeys or LipidMaps identifiers

get.pessimistic.rank.annotated <- function(annotated.candidate.list, inchikey1s = NULL, lmids = NULL, fcp.threshold = 0.9) {
  if(length(which(colnames(annotated.candidate.list) %in% "FCP")) == 0) {
    cat("Error: Given candidate list is not annotated with FCP.\n")
    return(NA)
  }
  # filter candidate list by FCP
  candidate.list.filtered <- annotated.candidate.list[!is.na(annotated.candidate.list[,"FCP"]),]
  candidate.list.filtered <- candidate.list.filtered[as.numeric(as.vector(candidate.list.filtered[,"FCP"])) >= fcp.threshold,]
  to.return <- get.pessimistic.rank(candidate.list.filtered, inchikey1s = inchikey1s, lmids = lmids)
  if(length(to.return) && is.na(to.return)) {return(to.return)}
  to.return <- 
    cbind(to.return, 
        as.numeric(as.vector(candidate.list.filtered[
          as.character(candidate.list.filtered[,"Identifier"]) == as.character(as.vector(to.return["ID"])[1,1]),"FCP"
      ])))
  names(to.return)[7] <- "FCP"
  return(to.return)
}
