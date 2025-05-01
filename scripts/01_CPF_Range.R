setwd("C:/Users/bat20/OneDrive - Brigham Young University/BYU/2024/Fall/Lab/BreastCancer")
#########################################################
#             Loading the Functions                     #
#########################################################

source("functions/readFiles.R")
source("functions/bakeFiles.R")
source("functions/mdMetrics.R")
source("functions/rocCurve.R")
library(tidymodels)
library(tidyverse)

#########################################################
#             Reading in the Data                       #
#########################################################

dataSets = c("GSE25065",
             "GSE21653",
             "GSE58644",
             "GSE62944",
             "GSE81538",
             "METABRIC",
             "GSE96058N")

genes = readLines("variables/genes.txt")

readFiles("GSE25055", columns = c("Class", genes))
readFiles(dataSets, columns = c("Class", genes))

dataSets = c("test", dataSets)

#########################################################
#             Preparing the Recipe                      #
#########################################################

set.seed(42)
split = initial_split(GSE25055, prop = 0.75)
train = training(split)
test = testing(split)

formula = Class ~ .

subFolder = "01_CPF_Range/"
if(!dir.exists("baked/")) dir.create("baked/")
if(!dir.exists(paste0("baked/", subFolder))) dir.create(paste0("baked/", subFolder))
if(!dir.exists(paste0("plots/", subFolder))) dir.create(paste0("plots/", subFolder))

recipe = recipe(formula, data = train) |>
    step_range(all_predictors()) |>
    step_mutate(Class = as.factor(Class)) |>
    prep(training = train)

train = bake(recipe, new_data = NULL)
write_tsv(train, paste0("baked/", subFolder, "train.tsv"))

bakeFiles(paste0(subFolder, dataSets), recipe)

saveRDS(recipe, "recipes/01_CPF_Range.rds")

#########################################################
#             Preparing the Model                       #
#########################################################

set.seed(42)
model = rand_forest(trees = 25,
                    mode = "classification",
                    engine = "ranger")

fitModel = model |>
            fit(formula, data = train)

saveRDS(fitModel, "models/01_CPF_Range.rds")

#########################################################
#             Obtaining Results                         #
#########################################################

dataSetVariables = mget(dataSets)

i = 1
for(dataSet in dataSetVariables) {
    setName = dataSets[i]
    mdMetrics(fitModel,
              dataSet,
              setName = setName,
              num = i)
    i = i + 1
}

i = 1
for(dataSet in dataSetVariables) {
    setName = dataSets[i]
    rocCurve(fitModel,
             dataSet,
             filename = setName,
             folder = paste0("plots/", subFolder),
             mdText = "CPF Range Standardization",
             title = paste0("GSE25055 prediciting on ", setName),
             num = i,
             plot = FALSE)
    i = i + 1
}

rm(i,
   genes,
   dataSet,
   dataSets,
   dataSetVariables,
   setName,
   subFolder,
   readFiles,
   bakeFiles,
   mdMetrics,
   rocCurve,
   model,
   formula,
   split)



