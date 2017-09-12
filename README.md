LipidFrag
=========

MetFrag based lipid identification

##### Requirements
- Java >= 1.6
- R >= 3.1.1
- R packages: gdata, ROCR, jsonlite

MetFrag
---------
As input LipidFrag requires a MetFrag generated candidate list by using a LipidMaps database. LipidMaps identifiers encode the lipid ontology which is used by LipidFrag. When using MetFragCLI or MetFragR LipidMaps can be included by downloading the database file from the [MSBI website](https://msbi.ipb-halle.de/~cruttkie/databases/lipidmaps.csv) and setting the database:
``R
settingsObject[["MetFragDatabaseType"]]<-"LocalCSV"
settingsObject[["LocalDatabasePath"]]<-"PATH_TP_LIPIDMAPS_CSV"
```
or for MetFragCLI
``bash
MetFragDatabaseType=LocalCSV
LocalDatabasePath=PATH_TP_LIPIDMAPS_CSV
```

The generated MetFrag (CSV) output file can then be used as input for 
``R
predict.lipidmaps.class
```
function. See lipidfrag_main.r to try an example.

R scripts
---------

##### lipidfrag_train.r
- includes functions to train prediction models based on lipid standard data
- uses files model_scores_pos.txt and model_scores_neg.txt for prediction models in positive and negative mode

```R
models.pos <- generate.model("../data/model_scores_pos.txt")
models.neg <- generate.model("../data/model_scores_neg.txt")
```

##### lipidfrag_predict.r
- includes functions to predict lipid classes based on trained models and MetFrag CSV files generated by using LipidMaps as candidate source

##### lipidfrag_functions.r
- includes general functions used by lipidfrag_train.r and lipidfrag_predict.r

##### lipidfrag_main.r
- includes an example showing how to use LipidFrag
- necessary data files are included in the repository

Additionals
-----------

##### Plot model value distribution
- after training the models score values can be plotted together with trained distribution parameters
```R
# distributions for negative ionization models
sapply(names(models.neg), function(class) {plot.model.data(models.neg[[class]], main = paste("Distribution of ", class, " (neg)", sep=""))})
# distributions for positive ionization models
sapply(names(models.pos), function(class) {plot.model.data(models.pos[[class]], main = paste("Distribution of ", class, " (pos)", sep=""))})
```


##### Plot ROC curves
- after training the models ROC curves can be plotted to get quality of classifications based on given standard data

```R
# ROC curves for negative ionization models
sapply(names(models.neg), function(class) calculate.roc("../data/model_scores_pos.txt", paste(class, models.pos[[class]]$type, sep="_"), main = paste("ROC of ", class, " (pos)", sep=""), variant = 1))
# ROC curves for positive ionization models
sapply(names(models.pos), function(class) calculate.roc("../data/model_scores_neg.txt", paste(class, models.neg[[class]]$type, sep="_"), main = paste("ROC of ", class, " (neg)", sep=""), variant = 1))
```
