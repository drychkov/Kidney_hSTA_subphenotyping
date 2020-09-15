source("functions.R")

ematWork = readRDS(file = "../Data_out/ematMerged_BatchCorr.rds")
sampleMeta = read.csv("../Data_out/sampleMeta.csv", stringsAsFactors = F, row.names = 1)

all(colnames(ematWork) == rownames(sampleMeta))

# ematSTA = ematWork[,sampleMeta[sampleMeta$class == "STA", "sample"]]
ematRN = ematWork[,sampleMeta$sample[sampleMeta$class %in% c("AR", "Normal")]]
outcomeRN = sampleMeta[colnames(ematRN), "class"] 
responseRN = ifelse(outcomeRN == "AR", 1, 0)

sam.resGE = sam(data = ematRN, #t(scale(t(ematRN))), 
                cl = responseRN, 
                method = d.stat,
                gene.names = rownames(ematRN),
                use.dm = TRUE,
                rand = 309
)

delta = findDelta(sam.resGE, 0.05)[2,1]
sam.matGE = summary(sam.resGE, delta, bonf = TRUE)@mat.sig


saveRDS(sam.resGE, "../Data_out/SAM_genes_RN_results.rds")


gRNup = sam.matGE[sam.matGE$q.value < 0.05 & sam.matGE$R.fold >= 1.5,]
gRNdown = sam.matGE[sam.matGE$q.value < 0.05 & sam.matGE$R.fold <= 1/1.5,]
gRNup = gRNup[order(-gRNup$R.fold),]
gRNdown = gRNdown[order(gRNdown$R.fold),]

write.csv(rbind(gRNup, gRNdown)[,4:7], file = "../Data_out/SAM_genesRN_FC15.csv")
write.csv(gRNup[,4:7], file = "../Data_out/SAM_genesRNup_FC15.csv")
write.csv(gRNdown[,4:7], file = "../Data_out/SAM_genesRNdown_FC15.csv")

plotHeatmap(ematRN[c(rownames(gRNup), rownames(gRNdown)),], sampleMeta, fontsize = 11)
plotPCA(ematRN[c(rownames(gRNup), rownames(gRNdown)),], sampleMeta, cex = 1.1, cexLab = 1.1)

plotUMAP(ematRN[c(rownames(gRNup), rownames(gRNdown)),], sampleMeta, n = 10)
plotUMAP(ematRN[c(rownames(gRNup), rownames(gRNdown)),], sampleMeta, n = 10, interest = "subclass")

# Identifying cluster significance ----------------------------------------

genes = c(rownames(gRNup), rownames(gRNdown))
shc_result <- shc(x = scale(t(ematRN[genes,])), 
                  metric = "euclidean", 
                  linkage = "ward", 
                  rcpp = TRUE,
                  n_sim = 100,
                  ci_emp = F,
                  alpha = 0.05)
summary(shc_result)
data.frame(result = head(shc_result$nd_type, 5),
           round(head(shc_result$p_norm, 5), 5),
           round(head(shc_result$p_emp, 5), 5))
plot(shc_result, 
     hang = -1, 
     alpha = 0.05, 
     group = sampleMeta[colnames(ematRN), "class"], 
     use_labs = FALSE,
     fwer = T,
     ci_emp = F)