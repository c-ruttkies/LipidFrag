source("lipidfrag_functions.r")


# reads lipid class training file with foreground and background scores to generate the prediction models
# generates a model for one particular lipid class
#
# example:
#
# Type,SpectrumLipidClassType,CandidateLipidClassType,Score
# fg,LMGP01,LMGP01,245.160932791981
# fg,LMGP01,LMGP01,186.459601776604
# fg,LMGP01,LMGP01,65.7972209963903
# bg,LMGP01,LMGL03,8.498761648540373
# bg,LMGP01,LMGL03,2.473272817399843
# bg,LMGP01,LMGP06,0.055950568492718
# ...
# fg,LMGL03,LMGL03,153.78184843018255
# fg,LMGL03,LMGL03,181.36545409410016
# fg,LMGL03,LMGL03,154.19201120395428
# bg,LMGL03,LMGP06,54.35734294611269
# bg,LMGL03,LMGP06,10.243225557067095
#

generate.model <- function(filename) {
  data <- read.csv(filename, header = T)
  # check for foreground class
  fg.global.indexes <- which(as.character(data[,"Type"]) == "fg")
  bg.global.indexes <- which(as.character(data[,"Type"]) == "bg")
  # detect foreground class names
  fg.global.lipid.classes <- unique(as.character(data[fg.global.indexes, "SpectrumLipidClassType"]))
  models<-list()
  
  sapply(fg.global.lipid.classes, function(lipid.class) {
    # get fore- and background training data of current lipid class
    fg.single.data <- data[as.character(data[,"Type"]) == "fg" & as.character(data[, "SpectrumLipidClassType"]) == lipid.class, ]
    bg.single.data <- data[as.character(data[,"Type"]) == "bg" & as.character(data[, "SpectrumLipidClassType"]) == lipid.class, ]
    # check data
    ## check for foreground scores being of the same lipid class
    if(length(which(fg.single.data[, "CandidateLipidClassType"] %in% lipid.class == FALSE)) != 0) {
      stop(paste("Foreground data of", lipid.class, "has scores from other lipid classes.", sep = " "))
    }
    if(dim(bg.single.data)[1] == 0) {
      stop("No background (bg) data given for", lipid.class, sep = " ")
    }
    ## check for background scores being of other lipid classes
    if(length(which(bg.single.data[, "CandidateLipidClassType"] %in% lipid.class == TRUE)) != 0) {
      stop(paste("Background data of", lipid.class, "has scores from the correct lipid classes.", sep = " "))
    }
    # start model generation
    fg.values <- fg.single.data[, "Score"]
    bg.values <- bg.single.data[, "Score"]
    # filtering values lower than 1 for background data
    bg.values <- bg.values[ bg.values > 1 ] 
    
    # calculate model parameters for foreground data
    fg.gamma.model = get.model.params.gamma(fg.values)
    # calculate model parameters for background data
    bg.gamma.model = get.model.params.gamma(bg.values)
    
    models[[lipid.class]] <<- list()
    models[[lipid.class]][["fg"]] <<- fg.gamma.model
    models[[lipid.class]][["bg"]] <<- bg.gamma.model
    models[[lipid.class]][["data_fg"]] <<- fg.values
    models[[lipid.class]][["data_bg"]] <<- bg.values
  })
  # return
  return(models)
}


# plot model distribution together with the underlying data

plot.model.data <- function(model, main = "") {
  # get boundaries
  maximum.value = max(c(model[["data_bg"]], model[["data_fg"]]))
  a<-hist(model[["data_bg"]], breaks="FD", plot=F)$density
  b<-hist(model[["data_fg"]], breaks="FD", plot=F)$density
  # plot data
  hist(model[["data_bg"]], ylim=c(0,max(c(a,b))+0.01), xlim=c(0,maximum.value), probability=T, breaks="FD", col=rgb(0,1,0,0.5), main=main, xlab="Scores", cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
  hist(model[["data_fg"]], probability=T, breaks="FD", add=T, col=rgb(0,0,1,0.5))
  # plot densities of model distributions
  scores.range <- seq(0, maximum.value, 1)
  lines(scores.range, sapply(scores.range, function(x) dgamma(x, model[["bg"]][1], 1 / model[["bg"]][2])), col = rgb(0,1,0,0.5), lwd=6)
  lines(scores.range, sapply(scores.range, function(x) dgamma(x, model[["fg"]][1], 1 / model[["fg"]][2])), col = rgb(0,0,1,0.5), lwd=6)
}


# perform ROC calculation for a specific classifier with 10-fold corss-validation

calculate.roc <- function(foreground.values, background.values, main = "") {
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
    #exp dist bg
    paramsBack <- get.model.params.gamma(backTrain)
    betaBack <- paramsBack[2]
    alphaBack <- paramsBack[1]
    
    #norm dist fg
    paramsFore <- get.model.params.gamma(foreTrain)
    betaFore <- paramsFore[2]
    alphaFore <- paramsFore[1]
    
    probsBack_back_mod <- sapply(backTest, function(x) get.weighted.gamma.density(x, 0.5, alphaBack, betaBack)) #spez
    probsFore_fore_mod <- sapply(foreTest, function(x) get.weighted.gamma.density(x, 0.5, alphaFore, betaFore))       #sens
    probsFore_back_mod <- sapply(foreTest, function(x) get.weighted.gamma.density(x, 0.5, alphaBack, betaBack))
    probsBack_fore_mod <- sapply(backTest, function(x) get.weighted.gamma.density(x, 0.5, alphaFore, betaFore))
    
    probsBack_back_modp <- probsBack_back_mod / (probsBack_back_mod + probsBack_fore_mod) #spez
    probsFore_fore_modp <- probsFore_fore_mod / (probsFore_back_mod + probsFore_fore_mod) #sens
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
  auc <- sum(sapply(2:length(data[6,]), function(x) ((sorted_spez[x] - sorted_spez[x - 1]) / 2) * (sorted_sens[x] + sorted_sens[x - 1])))
  
  pred <- prediction(as.numeric(returnValsProbabilities), returnValsLabels)
  perf <- performance( pred, "tpr", "fpr" )
  auc <- performance( pred, "auc" )
  
  par(mar=c(5.1,5.1,4.1,2.1))
  plot(perf, main = paste(main, " AUC: ", round(attributes(auc)$y.values[[1]], 3), sep = " "), cex.axis = 2, cex.lab = 2, lwd = 4)
  text(0.5, 0.5, paste("AUC:",round(attributes(auc)$y.values[[1]], 3)), cex=2)
}

