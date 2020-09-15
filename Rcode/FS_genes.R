source("functions.R")

ematWork = readRDS(file = "../Data_out/ematMerged_BatchCorr.rds")
sampleMeta = read.csv("../Data_out/sampleMeta.csv", stringsAsFactors = F, row.names = 1)

all(colnames(ematWork) == rownames(sampleMeta))

# ematSTA = ematWork[,sampleMeta[sampleMeta$class == "STA", "sample"]]
ematRN = ematWork[,sampleMeta$sample[sampleMeta$class %in% c("AR", "Normal")]]
outcomeRN = sampleMeta[colnames(ematRN), "class"] 
responseRN = ifelse(outcomeRN == "AR", 1, 0)


set.seed(309)
# inTraining <- caret::createDataPartition(responseRN, p = .8, list = FALSE)
# write.csv(inTraining, file = "../Data_out/inTrainingSamples.csv")
inTraining <- read.csv("../Data_out/inTrainingSamples.csv", stringsAsFactors = FALSE, header = TRUE)[,1]
train.GE <- t(scale(t(ematRN[,  inTraining])))
test.GE  <- t(scale(t(ematRN[, -inTraining])))


sam.resGE = sam(data = train.GE, 
                cl = responseRN[inTraining], 
                method = d.stat,
                gene.names = rownames(ematRN),
                use.dm = TRUE,
                rand = 309
)

delta = findDelta(sam.resGE, 0.05)[2,1]
sam.matGE = summary(sam.resGE, delta, bonf = TRUE)@mat.sig
sam.matGE = sam.matGE[order(-sam.matGE$R.fold),]

genesRN = sam.matGE[sam.matGE$q.value < 0.05 & 
                      (sam.matGE$R.fold >= 1.5 | sam.matGE$R.fold <= 1/1.5),]

### Find correlations
corGenesRN = sapply(genesRN, function(x) cor(train.GE[x,], responseRN[inTraining]))
#### Applying feseR

cl <-  makeCluster(detectCores() - 1, type = 'FORK', outfile = "")
registerDoParallel(cl)

genesFS.res <- combineFS(features = t(train.GE[genesRN,]), 
                     class = responseRN[inTraining],
                     univariate = 'corr', 
                     # n.percent = 0.95,
                     mincorr = max(abs(corGenesRN))*2/3,
                     multivariate = 'mcorr', 
                     maxcorr = 1,
                     wrapper = 'rfe.rf', 
                     number.cv = 5,
                     group.sizes = c(1:20),
                     extfolds = 100,
                     partition = 4/5,
                     metric = "ROC",
                     tolerance = 1)
genesFS.res
saveRDS(genesFS.res, "../Data_out/genesFS_results.rds")

impGenes = genesFS.res$opt.variables
write.csv(impGenes, "../Data_out/impGenes.csv")

plotUMAP(ematWork[impGenes, ], sampleFull, interest = "class", n = 15, pointSize = 1.3)
plotUMAP(ematRN[impGenes, ], sampleFull, interest = "class", n = 15, pointSize = 1.3)

## Model performance with selected features
ctrl <- trainControl(method = "repeatedcv",
                     number = 5,
                     repeats = 100,
                     savePredictions = "all",
                     summaryFunction = twoClassSummary,
                     classProbs = TRUE,
                     allowParallel = TRUE,
                     verbose = TRUE)

modelGE.train = train(x = t(train.GE[impGenes,]), 
                      y = factor(outcomeRN[inTraining], levels = c("Normal", "AR")),
                      method = "rf",
                      trControl = ctrl,
                      metric = "ROC")
modelGE.train
saveRDS(modelGE.train, "../Data_out/modelGE.trainSet.rds")

thres <- thresholder(modelGE.train, threshold = seq(0.3, 0.7, by = 0.01), final = TRUE)
ggplot(thres, aes(x = prob_threshold, y = Sensitivity)) + 
  geom_point() + 
  geom_point(aes(y = Specificity), col = "red")

test.predGE = predict(modelGE.train, newdata = t(test.GE[impGenes,]), type = "prob")

plotROC(test = test.predGE$Normal, target = outcomeRN[-inTraining])
confusionMatrix(as.factor(as.numeric(test.predGE$AR >= 0.57)), as.factor(responseRN[-inTraining]))


##### 
# Building the prediction model on the whole data with feature selected genes
modelGE = train(x = t(ematRN[impGenes,]), 
                y = factor(outcomeRN, levels = c("Normal", "AR")),
                method = "rf",
                trControl = ctrl,
                metric = "ROC")
modelGE

stopCluster(cl)


saveRDS(modelGE, "../Data_out/modelGE.allData.rds")
