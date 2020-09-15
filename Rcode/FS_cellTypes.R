source("functions.R")

xcellWork = read.csv(file = "../Data_out/xcellScores.csv", 
                     stringsAsFactors = F, row.names = 1) %>% as.matrix()
pxcell = read.csv(file = "../Data_out/xcell_pvalues.csv", 
                  stringsAsFactors = F, row.names = 1) %>% as.matrix()
sampleMeta = read.csv("../Data_out/sampleMeta.csv", stringsAsFactors = F, row.names = 1)
inTraining <- read.csv("../Data_out/inTrainingSamples.csv", stringsAsFactors = FALSE, header = TRUE)[,1]

all(colnames(xcellWork) == rownames(sampleMeta))

xcellRN = xcellWork[,sampleMeta$sample[sampleMeta$class %in% c("AR", "Normal")]]
sigCells = rownames(pxcell)[rowMedians(pxcell[,inTraining]) < 0.2]
# xcellRN = xcellRN[rownames(xcellRN) %in% sigCells, ]

responseRN = ifelse(sampleMeta[colnames(xcellRN), "class"] == "AR", 1, 0)


sam.resCT = sam(data = xcellRN[sigCells,inTraining], 
                cl = responseRN[inTraining], 
                method = wilc.stat, 
                gene.names = rownames(xcellRN), 
                rand = 309,
                R.unlog = FALSE
                #B = 1000
)

delta = findDelta(sam.resCT, 0.05)[2,1]

sam.sum = summary(sam.resCT, 0.1, bonf = TRUE)
sam.mat = sam.sum@mat.sig
sam.matCT = sam.mat[order(-sam.mat$R.fold),]


cellsRN = rownames(sam.matCT[sam.matCT$q.value <= 0.05,])

### Find correlations
corCellsRN = sapply(cellsRN, function(x) cor(xcellRN[x,inTraining,drop=T], responseRN[inTraining]))


train.CT <- t(scale(t(xcellRN[cellsRN,  inTraining])))
test.CT  <- t(scale(t(xcellRN[cellsRN, -inTraining])))
rownames(train.CT) = make.names(rownames(train.CT))
rownames(test.CT) = make.names(rownames(test.CT))


cl <-  makeCluster(detectCores() - 1, type = 'FORK', outfile = "")
registerDoParallel(cl)

cellsFS.res <- combineFS(features = t(train.CT), 
                         class = responseRN[inTraining],
                         univariate = 'corr', 
                         # n.percent = 0.95,
                         mincorr = max(abs(corCellsRN))*2/3,
                         multivariate = 'mcorr', 
                         maxcorr = 1,
                         wrapper = 'rfe.rf', 
                         number.cv = 5,
                         group.sizes = c(1:10),
                         extfolds = 100,
                         partition = 4/5,
                         metric = "ROC",
                         tolerance = 1)
cellsFS.res
saveRDS(cellsFS.res, "../Data_out/cellsFS_results.rds")

impCells = cellsFS.res$opt.variables
write.csv(impCells, "../Data_out/impCells.csv")

plotUMAP(xcellRN[names(impCells), ], sampleFull, interest = "class", n = 15, pointSize = 1.3)

## Model performance with selected features
ctrl <- trainControl(method = "repeatedcv",
                     number = 5,
                     repeats = 100,
                     savePredictions = "all",
                     summaryFunction = twoClassSummary,
                     classProbs = TRUE,
                     allowParallel = TRUE,
                     verbose = TRUE)

modelCT.train = train(x = t(train.CT[impCells,]), 
                y = factor(outcomeRN[inTraining], levels = c("Normal", "AR")),
                method = "rf",
                trControl = ctrl,
                metric = "ROC")
modelCT.train
saveRDS(modelCT.train, "../Data_out/modelCT.trainSet.rds")

thres <- thresholder(modelCT.train, threshold = seq(0.3, 0.7, by = 0.01), final = TRUE)
ggplot(thres, aes(x = prob_threshold, y = Sensitivity)) + 
  geom_point() + 
  geom_point(aes(y = Specificity), col = "red")

test.predCT = predict(modelCT.train, newdata = t(test.CT), type = "prob")

plotROC(test = test.predCT$Normal, target = outcomeRN[-inTraining])
confusionMatrix(as.factor(as.numeric(test.predCT$AR >= 0.39)), as.factor(responseRN[-inTraining]))
confusionMatrix(as.factor(as.numeric(predict(modelCT.train, newdata = t(train.CT), type = "prob")$AR >= 0.39)), as.factor(responseRN[inTraining]))


###
# Building prediction model on the whole CT data
rownames(xcellRN) = make.names(rownames(xcellRN))
modelCT = train(x = t(xcellRN[impCells,]), 
                y = factor(outcomeRN, levels = c("Normal", "AR")),
                method = "rf",
                trControl = ctrl,
                metric = "ROC")
modelCT
saveRDS(modelCT, "../Data_out/modelCT.allData.rds")

stopCluster(cl)
