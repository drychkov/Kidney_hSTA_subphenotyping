source("functions.R")


sarwalValidData = read.csv(file = "../Data_in/SarwalLab_SNS01_longitudinal.csv", stringsAsFactors = F)
rownames(sarwalValidData) = sarwalValidData$sample

s = intersect(sarwalValidData$sample, sampleMeta[sampleMeta$class %in% "STA", "sample"])
sarwalValidData = sarwalValidData[s,]

sarwalValidData$InstaScore = sampleMeta[s, "InstaScore"]
sarwalValidData$classIS = sampleMeta[s, "classIS"]


confusionMatrix(as.factor(ifelse(sarwalValidData[, "classIS"] == "hSTA/mAR",1,0)),
                as.factor(sarwalValidData[, "graftLoss"]), mode = "everything", positive = "1")
## Plots