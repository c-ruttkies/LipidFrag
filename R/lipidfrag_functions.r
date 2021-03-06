require(gdata)
require(MASS)

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

# calculate model parameters of gamma distribution by mle

get.model.params.gamma.mle <- function(values) {
  estimates <- fitdistr(values, "gamma", start=list(shape=1, rate=1))$estimate
  beta <- 1 / estimates["rate"]
  alpha <- estimates["shape"]
  return(c(alpha, beta))
}



# calculate posterior foreground class probability (FCP) for a given value

get.posterior.foreground <- function(value, model, detailed = F) {
  type <- model[["type"]]
  backDensity <- get.weighted.gamma.density(value, 0.5, model[["bg"]][1], model[["bg"]][2])
  foreDensity <- 0
  if(type == "mix") {
    # calculate mixture model foreground probability weighted by the number of classes: 1 / length(model[["fg"]])
    sapply(1:length(model[["fg"]]), function(index) {
      foreDensity <<- foreDensity + get.weighted.gamma.density(value, 1 / length(model[["fg"]]), model[["fg"]][[index]][1], model[["fg"]][[index]][2])
    })
    foreDensity <<- foreDensity * 0.5
  } else {
    foreDensity <- get.weighted.gamma.density(value, 0.5, model[["fg"]][1], model[["fg"]][2])
  }
  if(value == 0) {
    if(detailed) {
      return(data.frame(FCP=0,ForeProb=0,BackProb=0))
    } else {
      return(data.frame(FCP=0))
    }
  }
  if(detailed) {
    return(data.frame(FCP=(foreDensity / (foreDensity + backDensity)),ForeProb=foreDensity,BackProb=backDensity))
  } else {
    return(data.frame(FCP=(foreDensity / (foreDensity + backDensity))))
  }
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

# get pessimistic rank of candidates matching given InChIKeys or LipidMaps identifiers
# Output: Rank, BetterCandidates, EqualCandidates, WorseCandidates, TotalCandiates, Identifier

get.pessimistic.rank <- function(candidate.list, inchikey1s = NULL, lmids = NULL) {
  if(is.null(inchikey1s) && is.null(lmids)) {
    stop("Error: InChIKey1 or LMID needed.\n")
  }
  if(is.null(lmids)) {
    correct.scores <- candidate.list[as.character(candidate.list[,c("InChIKey1")]) %in% inchikey1s, c("Score", "Identifier")]
    wrong.scores <- candidate.list[!as.character(candidate.list[,c("InChIKey1")]) %in% inchikey1s, c("Score", "InChIKey1")]
  }
  else {
    inchikey1s <- as.character(candidate.list[as.character(candidate.list[,"Identifier"]) %in% lmids, "InChIKey1"])
    correct.scores <- candidate.list[as.character(candidate.list[,c("InChIKey1")]) %in% inchikey1s, c("Score", "Identifier")]
    wrong.scores <- candidate.list[!as.character(candidate.list[,c("InChIKey1")]) %in% inchikey1s, c("Score", "InChIKey1")]
    # filter by inchikey1
    wrong.scores <- sapply(unique(as.character(wrong.scores[,2])), function(x) wrong.scores[which(as.character(wrong.scores[,2]) == x)[1],1])
  }
  if(dim(correct.scores)[1] == 0) {
    print("Error: Identifier not found in candidate list.\n")
    return("NA")
  }
  max.correct.score <- correct.scores[which.max(correct.scores[,1]),]
  if(length(wrong.scores) == 0) {
    return.vals <- data.frame(Rank=1,BC=0,EC=0,WC=0,TC=1,ID=max.correct.score[,2],FPIncludeRate=0)
    return(return.vals)
  }
  equal.cands <- length(which(max.correct.score[,1] == wrong.scores))
  better.cands <- length(which(max.correct.score[,1] < wrong.scores))
  worse.cands <- length(which(max.correct.score[,1] > wrong.scores))
  return.vals <- data.frame(Rank = better.cands + equal.cands + 1, BC = better.cands, 
                            EC = equal.cands, WC = worse.cands, TC = length(wrong.scores) + 1, 
                            ID = max.correct.score[,2], FCPIncludeRate = (length(wrong.scores) / (length(wrong.scores) + 1)))
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
  
  if(is.null(inchikey1s) && is.null(lmids)) {
    cat("Error: InChIKey1 or LMID needed.\n")
    return(NA)
  }
  if(is.null(lmids)) {
    correct.scores <- candidate.list.filtered[as.character(candidate.list.filtered[,c("InChIKey1")]) %in% inchikey1s, c("Score", "Identifier")]
    wrong.scores <- candidate.list.filtered[!as.character(candidate.list.filtered[,c("InChIKey1")]) %in% inchikey1s, c("Score", "InChIKey1")]
    #unfiltered
    correct.scores.unfiltered <- annotated.candidate.list[as.character(annotated.candidate.list[,c("InChIKey1")]) %in% inchikey1s, c("Score", "Identifier")]
    wrong.scores.unfiltered <- annotated.candidate.list[!as.character(annotated.candidate.list[,c("InChIKey1")]) %in% inchikey1s, c("Score", "InChIKey1")]
  } else {
    inchikey1s <- as.character(candidate.list.filtered[as.character(candidate.list.filtered[,"Identifier"]) %in% lmids, "InChIKey1"])
    correct.scores <- candidate.list.filtered[as.character(candidate.list.filtered[,c("InChIKey1")]) %in% inchikey1s, c("Score", "Identifier")]
    wrong.scores <- candidate.list.filtered[!as.character(candidate.list.filtered[,c("InChIKey1")]) %in% inchikey1s, c("Score", "InChIKey1")]
    # filter by inchikey1
    wrong.scores <- sapply(unique(as.character(wrong.scores[,2])), function(x) wrong.scores[which(as.character(wrong.scores[,2]) == x)[1],1])
    #
    # unfiltered
    #
    inchikey1s.unfiltered <- as.character(annotated.candidate.list[as.character(annotated.candidate.list[,"Identifier"]) %in% lmids, "InChIKey1"])
    correct.scores.unfiltered <- annotated.candidate.list[as.character(annotated.candidate.list[,c("InChIKey1")]) %in% inchikey1s, c("Score", "Identifier")]
    wrong.scores.unfiltered <- annotated.candidate.list[!as.character(annotated.candidate.list[,c("InChIKey1")]) %in% inchikey1s, c("Score", "InChIKey1")]
    # filter by inchikey1
    wrong.scores.unfiltered <- sapply(unique(as.character(wrong.scores.unfiltered[,2])), function(x) wrong.scores.unfiltered[which(as.character(wrong.scores.unfiltered[,2]) == x)[1],1])
  }
  if(dim(correct.scores)[1] == 0) {
    cat("Error: Identifier not found in candidate list.\n")
    return(NA)
  }
  max.correct.score <- correct.scores[which.max(correct.scores[,1]),]
  if(length(wrong.scores) == 0) {
    return.vals <- data.frame(Rank=1,BC=0,EC=0,WC=0,TC=1,ID=as.character(max.correct.score[,2]),FPIncludeRate=0)
    return.vals <- 
      cbind(return.vals, paste(0,1,fcp.threshold),
            candidate.list.filtered[
              as.character(candidate.list.filtered[,"Identifier"]) == as.character(as.vector(return.vals["ID"])[1,1]),
              c("FCP","LipidMapsClass")
              ])
    return(return.vals)
  }
  equal.cands <- length(which(max.correct.score[,1] == wrong.scores))
  better.cands <- length(which(max.correct.score[,1] < wrong.scores))
  worse.cands <- length(which(max.correct.score[,1] > wrong.scores))
  return.vals <- data.frame(Rank = better.cands + equal.cands + 1, BC = better.cands, 
                            EC = equal.cands, WC = worse.cands, TC = length(wrong.scores) + 1, 
                            ID = max.correct.score[,2], FPIncludeRate=(length(wrong.scores) / (length(wrong.scores.unfiltered) + 1)), Values=paste(length(wrong.scores), (length(wrong.scores.unfiltered) + 1), fcp.threshold))
  
  return.vals <- 
    cbind(return.vals, 
          candidate.list.filtered[
            as.character(candidate.list.filtered[,"Identifier"]) == as.character(as.vector(return.vals["ID"])[1,1]),
            c("FCP","LipidMapsClass")
            ])
  return(return.vals)
}

# retrieve upper liebisch level for given lipidmaps id

get.liebisch.annotation <- function(lmid, detailed = F) {
  require(jsonlite)
  out <- try(fromJSON(paste("http://www.lipidmaps.org/rest/compound/lm_id/",lmid,"/name",sep="")))
  if(!detailed) {
    anno <- data.frame(Liebisch="NA")
  } else {
    anno <- data.frame(Liebisch="NA",CommonName="NA")
  }
  if(!is.null(out$name)) {
    format.ok <- length(grep("^[A-z]*-*[A-z]*[0-9]*\\([A-Za-z]*-*[0-9]+:[0-9]+.*\\)", out$name))
    if(format.ok > 0) {
      prefix <- gsub("(^[A-z]*-*[A-z]*[0-9]*)\\(.*", "\\1", out$name)
      content <- gsub("^[A-z]*-*[A-z]*[0-9]*\\((.*)\\).*", "\\1", out$name)
      content <- gsub("\\(([0-9]*,*[A-Z]*)+\\)", "", content)
      content.split <- unlist(strsplit(content, "/"))
      liebisch.content <- paste(content.split[order(as.numeric(gsub(":",".",gsub("[A-Za-z]*-*", "", content.split))))], collapse="_")
      if(detailed) {
        anno <- data.frame(Liebisch=paste(prefix,"(",liebisch.content,")",sep=""), CommonName=out$name)
      } else {
        anno <- data.frame(Liebisch=paste(prefix,"(",liebisch.content,")",sep=""))
      }
    }
  }
  return(anno)
}
