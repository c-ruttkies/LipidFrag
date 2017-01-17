source("lipidfrag_functions.r")

# helper function to generate model for single lipid class

generate.single.model <- function(fg.values, bg.values, type="sin") {
  if(type != "sin" & type != "mix") {
    stop(paste("Model type",type,"not know. Use: sin, mix"))
  }
  # calculate model parameters for foreground data single model
  if (type == "sin") {
    fg.gamma.model = get.model.params.gamma(unlist(fg.values))
  }
  else {
    # generate a mixture model for the foreground data
    fg.gamma.model <- list()
    sapply(names(fg.values), function(name) {
      fg.gamma.model[[name]]<<-get.model.params.gamma(fg.values[[name]])
    })
  }
  # calculate model parameters for background data
  bg.gamma.model = get.model.params.gamma(bg.values)
  
  model <- list()
  model[["fg"]] <- fg.gamma.model
  model[["bg"]] <- bg.gamma.model
  model[["data_fg"]] <- fg.values
  model[["data_bg"]] <- bg.values
  model[["type"]]<-type
  
  return(model)
}


# reads lipid class training file with foreground and background scores to generate the prediction models
# generates a model for one particular lipid class
#
# example:
#
# example ClassifierType
# LMSP0201_LMSP0202_mix
# LMSP0201_LMSP0202_sin
# LMGP0101_sin

generate.model <- function(filename) {
  data <- read.csv(filename, header = T, stringsAsFactor = F)
  # check for foreground class
  fg.global.indexes <- which(as.character(data[,"Type"]) == "fg")
  bg.global.indexes <- which(as.character(data[,"Type"]) == "bg")
  # get data
  fg.data <- data[fg.global.indexes,]
  bg.data <- data[bg.global.indexes,]
  # detect foreground class names
  fg.global.classifier.types <- unique(as.character(data[fg.global.indexes, "ClassifierType"]))
  models<-list()
  # generate classifier for each classifier type
  sapply(fg.global.classifier.types, function(classifier.type) {
    classifier.type.split <- unlist(strsplit(classifier.type, "_"))
    type <- classifier.type.split[length(classifier.type.split)]
    fg.values <- list()
    bg.values <- list()
    classifier.name <- ""
    sapply(1:(length(classifier.type.split) - 1), function(index) {
      fg.values[[classifier.type.split[index]]] <<- fg.data[as.character(fg.data[,"CandidateLipidClassType"]) == classifier.type.split[index], "Score"]
      # filtering values lower than 10 for foreground data due to low quality spectra
      fg.values[[classifier.type.split[index]]] <<- fg.values[[index]][fg.values[[index]] > 10]
      # build classifier name
      classifier.name <<- paste(classifier.name, classifier.type.split[index], sep="_")
      if(length(which(fg.data[fg.data[,"CandidateLipidClassType"] == classifier.type.split[index], "SpectrumLipidClassType"] %in% classifier.type.split[index] == FALSE)) != 0) {
        stop(paste("Foreground data of", classifier.type.split[index], "has scores from other lipid classes.", sep = " "))
      }
    })
    # collect background data
    bg.values[[index]] <- bg.data[bg.data[,"CandidateLipidClassType"] == paste(classifier.type.split[1:(length(classifier.type.split)-1)],collapse="_"), "Score"]
    bg.values <- unlist(bg.values)
    # filtering values lower than 1 for background data
    bg.values <- bg.values[ bg.values > 1 ]
    # check data
    if(length(bg.values) == 0) {
      stop("No background (bg) data given for", classifier.type, sep = " ")
    }
    classifier.name <- gsub("^_*", "", classifier.name)
    models[[classifier.name]] <<- generate.single.model(fg.values, bg.values, type)

  })
  # return
  return(models)
}

# plot model distribution together with the underlying data

plot.model.data <- function(model, main = "Model", ...) {
  # get boundaries
  maximum.value = max(c(model[["data_bg"]], unlist(model[["data_fg"]])))
  a<-hist(model[["data_bg"]], breaks="FD", plot=F)$density
  b<-hist(unlist(model[["data_fg"]]), breaks="FD", plot=F)$density
  # plot data
  hist(model[["data_bg"]], ylim=c(0,max(c(a,b))+0.01), xlim=c(0,maximum.value), probability=T, breaks="FD", col=rgb(0,1,0,0.5), main=main, xlab="Scores", ...)
  hist(unlist(model[["data_fg"]]), probability=T, breaks="FD", add=T, col=rgb(0,0,1,0.5))
  # plot densities of model distributions
  scores.range <- seq(0, maximum.value, 1)
  lines(scores.range, sapply(scores.range, function(x) dgamma(x, model[["bg"]][1], 1 / model[["bg"]][2])), col = rgb(0,1,0,0.5), lwd=6)
  if(model$type == "mix") {
    sapply(1:length(model[["fg"]]), function(index) {
      lines(scores.range, sapply(scores.range, function(x) dgamma(x, model[["fg"]][[index]][1], 1 / model[["fg"]][[index]][2])), col = rgb(0,0,1,0.5), lwd=6)
    })
  } else {
    lines(scores.range, sapply(scores.range, function(x) dgamma(x, model[["fg"]][1], 1 / model[["fg"]][2])), col = rgb(0,0,1,0.5), lwd=6)
  }
}


# perform ROC calculation for a specific classifier with 10-fold corss-validation

calculate.roc <- function(foreground.values, background.values, main = "", return.fcps = F, ...) {
  require(ROCR)
  
  is.foreground <- c(rep(T, length(foreground.values)), rep(F, length(background.values)))
  groups_fore <- rep(1:10, ceiling(length(foreground.values) / 10))[1:length(foreground.values)]
  groups_back <- rep(1:10, ceiling(length(background.values) / 10))[1:length(background.values)]
  
  groups_fore <- sample(groups_fore, length(groups_fore))
  groups_back <- sample(groups_back, length(groups_back))
  
  groups <- c(groups_fore, groups_back)
  
  data <- c(foreground.values, background.values)
  
  group.numbers <- sort(unique(groups))
  
  relsForeAll <- c()
  relsBackAll <- c()
  
  returnValsProbabilities <- c()
  returnValsLabels <- c()
  
  for(i in 1:10) {
    test.group.number <- i
    backTrain <- data[groups != i & !is.foreground]
    backTest <- data[groups==i & !is.foreground]
    foreTrain <- data[groups!=i & is.foreground]
    foreTest <- data[groups==i & is.foreground]
    #gamma dist bg
    paramsBack <- get.model.params.gamma(backTrain)
    betaBack <- paramsBack[2]
    alphaBack <- paramsBack[1]
    
    #gamma dist fg
    paramsFore <- get.model.params.gamma(foreTrain)
    betaFore <- paramsFore[2]
    alphaFore <- paramsFore[1]
    
    probsBack_back_mod <- sapply(backTest, function(x) get.weighted.gamma.density(x, 0.5, alphaBack, betaBack)) #spez
    probsFore_fore_mod <- sapply(foreTest, function(x) get.weighted.gamma.density(x, 0.5, alphaFore, betaFore))       #sens
    probsFore_back_mod <- sapply(foreTest, function(x) get.weighted.gamma.density(x, 0.5, alphaBack, betaBack))
    probsBack_fore_mod <- sapply(backTest, function(x) get.weighted.gamma.density(x, 0.5, alphaFore, betaFore))
    
    probsBack_back_modp <- probsBack_back_mod / (probsBack_back_mod + probsBack_fore_mod)
    probsFore_fore_modp <- probsFore_fore_mod / (probsFore_back_mod + probsFore_fore_mod)
    probsFore_back_modp <- probsFore_back_mod / (probsFore_back_mod + probsFore_fore_mod)
    probsBack_fore_modp <- probsBack_fore_mod / (probsBack_back_mod + probsBack_fore_mod)
    
    returnValsProbabilities <- c(returnValsProbabilities, probsFore_fore_modp)
    returnValsProbabilities <- c(returnValsProbabilities, probsBack_fore_modp)
    names(returnValsProbabilities) <- NULL
    returnValsLabels <- c(returnValsLabels, rep("fg", length(probsFore_fore_modp)))
    returnValsLabels <- c(returnValsLabels, rep("bg", length(probsBack_fore_modp)))
    
    relsForeAll <- c(relsForeAll, probsFore_fore_modp / (probsFore_fore_modp + probsFore_back_modp))
    relsBackAll <- c(relsBackAll, probsBack_back_modp / (probsBack_fore_modp + probsBack_back_modp))
    
  }
  
  
  data <- sapply(seq(0, 1, 0.01), function(x) calculate.single.roc(relsForeAll, relsBackAll, x))
  sorted_spez <- sort(1-data[7,])
  sorted_sens <- sort(data[6,])
  
  pred <- prediction(as.numeric(returnValsProbabilities), returnValsLabels)
  perf <- performance( pred, "tpr", "fpr" )
  auc <- performance( pred, "auc" )
  
  plot(perf, main = paste(main, " AUC: ", round(attributes(auc)$y.values[[1]], 3), sep = " "), ...)
  text(0.5, 0.5, paste("AUC:",round(attributes(auc)$y.values[[1]], 3)), cex=2)
  
  if(return.fcps) return(cbind(returnValsProbabilities, returnValsLabels))
}



calculate.roc <- function(filename, classifier.to.test, main = "", return.fcps = F, ...) {
  require(ROCR)
  
  data <- read.csv(filename, header = T, stringsAsFactors = F)
  data.to.test <- data[data[,"ClassifierType"] == classifier.to.test,]
  spectrum.ids <- unique(data.to.test[,"SpectrumID"])
  
  foreground.values <- data.to.test[data.to.test[,"Type"] == "fg", "Score"]
  background.values <- data.to.test[data.to.test[,"Type"] == "bg", "Score"]
  
  is.foreground <- c(rep(T, length(foreground.values)), rep(F, length(background.values)))
  
  unique.fold.ids <- sample(((1:length(spectrum.ids))%%10)+1)
  names(unique.fold.ids) <- spectrum.ids
  fold.ids <- unique.fold.ids[data.to.test[,"SpectrumID"]]
  
  returnValsProbabilities <- c()
  returnValsLabels <- c()

  for(i in 1:10) {
    tmp.trainings.file <- tempfile(pattern = "file", tmpdir = tempdir(), fileext = "")
    # write out training data for this test fold
    write.table(data.to.test[fold.ids != i,], tmp.trainings.file, sep=",", quote=F)
    model <- generate.model(tmp.trainings.file)
    
    test.data <- data.to.test[fold.ids == i,]
    
    returnValsProbabilities <- sapply(test.data[,"Score"], function(x) get.posterior.foreground(x, model))
    
    returnValsLabels <- test.data[,"Type"]

  }

  pred <- prediction(as.numeric(returnValsProbabilities), returnValsLabels)
  perf <- performance( pred, "tpr", "fpr" )
  auc <- performance( pred, "auc" )

  plot(perf, main = paste(main, " AUC: ", round(attributes(auc)$y.values[[1]], 3), sep = " "), ...)
  text(0.5, 0.5, paste("AUC:",round(attributes(auc)$y.values[[1]], 3)), cex=2)

  if(return.fcps) return(cbind(returnValsProbabilities, returnValsLabels))
}

