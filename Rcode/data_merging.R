source("functions.R")


sampleMeta = read.csv("../Data_out/sampleMeta.csv", stringsAsFactors = F, row.names = 1)
esetList = readRDS("../Data_out/esetListNoDup.rds")

ematList = cleanStudyData(esetList, sampleMeta)
ematList = lapply(ematList, function(x) geneId2Symbol(x))
saveRDS(ematList, file = "../Data_out/ematList.rds")

# ematList = lapply(ematList, function(x) x[, colnames(x) %in% sampleFull$sample])

geneIds = Reduce(intersect, lapply(ematList, function(x) rownames(x)))
ematList2 = foreach(studyName = names(ematList)) %do% {ematNow = ematList[[studyName]][geneIds,]}
ematListScaled = lapply(ematList2, function(emat) (emat - mean(emat)) / sd(emat)) # zero centering

# Combine all data to one matrix
ematListScaled = lapply(ematListScaled, function(emat) normalizeQuantiles(emat))
ematMerged = do.call(cbind, ematListScaled)
saveRDS(ematMerged, file = "../Data_out/ematMerged_noBatchCorr.rds")

ematMerged = readRDS("../Data_out/ematMerged_noBatchCorr.rds")

## Batch correction with Combat
mod = model.matrix(~as.factor(subclass), data = sampleMeta)
ematMerged.Combat = ComBat(ematMerged, 
                           batch = sampleMeta$batch, 
                           mod = mod, 
                           BPPARAM = MulticoreParam(workers = detectCores() - 1)
)
saveRDS(ematMerged.Combat, file = "../Data_out/ematMerged_BatchCorr.rds")




###plots
plotPCA(ematMerged, sampleMeta, interest = "study")
plotPCA(ematMerged.Combat, sampleMeta, interest = "study")
plotPCA(ematMerged.Combat, sampleMeta, interest = "class")

plotUMAP(ematMerged, sampleMeta, interest = "study", n = 10, pointSize = 0.8)
plotUMAP(ematMerged.Combat, sampleMeta, interest = "study", n = 10, pointSize = 0.8)
plotUMAP(ematMerged.Combat, sampleMeta, interest = "class", n = 10, pointSize = 0.8)
plotUMAP(ematMerged.Combat, sampleMeta, interest = "subclass", n = 10, pointSize = 0.8)


rm(geneIds, esetList, ematNow, ematList, ematList2, ematListScaled, ematMerged, ematMerged.Combat, mod)