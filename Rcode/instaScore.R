source("functions.R")


cellList = setNames(object = make.names(rownames(xcellRN)), rownames(xcellRN))
impCells = cellList[match(impCells, cellList)]


geneCellMatrixRN = rbind(ematRN[impGenes, ], xcellRN[names(impCells),])
geneCellMatrixSTA = t(scale(t(rbind(ematSTA[impGenes, ], xcellSTA[names(impCells),]))))



###

inTraining <- read.csv("../Data_out/inTrainingSamples.csv", stringsAsFactors = FALSE, header = TRUE)[,1]
trainRN = t(scale(t(geneCellMatrixRN[,inTraining])))
testRN = t(scale(t(geneCellMatrixRN[,-inTraining])))


ctrl <- trainControl(method = "repeatedcv",
                     number = 5,
                     repeats = 100,
                     savePredictions = "all",
                     summaryFunction = twoClassSummary,
                     classProbs = TRUE,
                     allowParallel = TRUE,
                     verboseIter = TRUE)

modelIS.train = train(x = t(trainRN), 
                 y = factor(outcomeRN[inTraining], levels = c("Normal", "AR")),
                 # preProcess = "scale",
                 method = "glm",
                 trControl = ctrl,
                 metric = "ROC")

modelIS.train
saveRDS(modelIS.train, "../Data_out/modelIS_train.rds")

##

modelIS.train = readRDS("../Data_out/modelIS_train.rds")

thresSC <- thresholder(modelIS.train, threshold = seq(0.5, 0.6, by = 0.005), final = TRUE)
ggplot(thresSC, aes(x = prob_threshold, y = Sensitivity)) + 
  geom_point() + 
  geom_point(aes(y = Specificity), col = "red") + 
  theme_classic() # 0.575


test.tmp = predict(modelIS.train, newdata = t(testRN), type = "prob")
plotROC(test = test.tmp$Normal, target = outcomeRN[-inTraining])
confusionMatrix(as.factor(as.numeric(test.tmp$AR > 0.5)), 
                as.factor(responseRN[-inTraining]), mode = "everything")
confusionMatrix(as.factor(as.numeric(predict(modelIS.train, newdata = t(trainRN), type = "prob")$AR > 0.5)), 
                as.factor(responseRN[inTraining]), mode = "everything")


modelIS.all = train(x = scale(t(geneCellMatrixRN)), 
                 y = factor(outcomeRN, levels = c("Normal", "AR")),
                 method = "glm",
                 trControl = ctrl,
                 metric = "ROC")

modelIS.all
saveRDS(modelIS.all, "../Data_out/modelIS_allData.rds")


## Applying to hSTA
sampleMeta$InstaScore = instaScore(geneCellMatrixSTA, model = modelIS.all)
sampleMeta$classIS = sampleMeta$class 

sampleMeta[colnames(geneCellMatrixSTA), "classIS"] = ifelse(sampleMeta$InstaScore > 0, "hSTA/mAR", "hSTA/mSTA")
table(sampleFull$classIS)

write.csv(sampleMeta, file = "../Data_out/sampleMeta_IS.csv")


plotUMAP(rbind(ematWork[impGenes, ], xcellWork[names(impCells),]), sampleFull, interest = "classByIS", n = 10, pointSize = 1.3)
