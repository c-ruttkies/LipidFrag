# generate models
source("lipidfrag_train.r")
source("lipidfrag_predict.r")

# initialise the prediction models
models.pos <- generate.model("../data/model_scores_pos.txt")
models.neg <- generate.model("../data/model_scores_neg.txt")

# predict the lipid (main-/sub-)class for a given MetFrag annotation
predicted.lipidmaps.class.neg <- predict.lipidmaps.class("../data/cel_results/neg/Cel_1-A,8_05_21044.914.914.0_variant_0.csv", models.neg, sep="|")
predicted.lipidmaps.class.pos <- predict.lipidmaps.class("../data/cel_results/pos/Cel_AutoMS_2-A,2_05_18615.944.944.0_variant_2.csv", models.pos, sep="|")
