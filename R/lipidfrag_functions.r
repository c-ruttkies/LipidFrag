require(gdata)

#######
#
# NEEDED FUNCTIONS FOR LIPIDFRAG
#
# this part contains functions needed for model generation and lipid class prediction of LipidFrag
# - gamma distribution density calculation
# - gamma distribution based model calculation
#
#######


# calculate weighted density of gamma distribution

get.weighted.gamma.density <- function(value, pi, alpha, beta) {
  return(pi * dgamma(value, alpha, 1/beta))
}


# calculate model parameters of gamma distribution

get.model.params.gamma <- function(values) {
  beta <- var(values) / mean(values)
  alpha <- mean(values) / beta
  return(c(alpha, beta))
}


# calculate posterior foreground class probability (FCP) for a given value

get.posterior.foreground <- function(value, model) {
  backDensity <- get.weighted.gamma.density(value, 0.5, model[["bg"]][1], model[["bg"]][2])
  foreDensity <- get.weighted.gamma.density(value, 0.5, model[["fg"]][1], model[["fg"]][2])
  
  if(value == 0) return(0.0)
  return(foreDensity / (foreDensity + backDensity))
}


# calculate the posterior foreground class probability (FCP) for a given MetFrag result list

get.fcps.candidate.list <- function(candidate.list, lipid.class.models) {
  class.names.model <- names(models)
  fcps <- t(apply(candidate.list, 1, function(row) {
    identifier <- as.character(row["Identifier"])
    score <- as.numeric(row["FragmenterScore"])
    matching.model.index <- which(startsWith(identifier, class.names.model))
    if(length(matching.model.index) == 0) { return(c(NA, identifier)) }
    if(length(matching.model.index) == 0) stop("More than one matching model for Identifier ", identifier, ". ", "Which shall I use?", sep="")
    return(c(get.posterior.foreground(score, lipid.class.models[[matching.model.index]]), class.names.model[matching.model.index]))
  }))
  annotated.candidate.list <- cbind(candidate.list, fcps)
  colnames(annotated.candidate.list) <- c(colnames(annotated.candidate.list), c("FCP", "LipidClass"))
  return(annotated.candidate.list)
} 

# calculate ROC value for given background and foreground values

calculate.single.roc <- function(relsFore, relsBack, thresh) {
  truePos<-table(relsFore>thresh)["TRUE"]
  trueNeg<-table(relsBack>thresh)["TRUE"]
  falseNeg<-length(relsFore)-truePos
  falsePos<-length(relsBack)-trueNeg
  if(is.na(falseNeg)){falseNeg=0}
  if(is.na(falsePos)){falsePos=0}
  if(is.na(trueNeg)){trueNeg=0}
  if(is.na(truePos)){truePos=0}
  
  spec<-truePos/(truePos+falseNeg)
  sens<-trueNeg/(trueNeg+falsePos)
  if(is.nan(sens)) sens=0
  if(is.nan(spec)) spec=0
  return(cbind(truePos,trueNeg,falseNeg,falsePos,thresh,sens,spec))
}
