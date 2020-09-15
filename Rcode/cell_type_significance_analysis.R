source("functions.R")

xcellWork = read.csv(file = "../Data_out/xcellScores.csv", 
                     stringsAsFactors = F, row.names = 1) %>% as.matrix()
pxcell = read.csv(file = "../Data_out/xcell_pvalues.csv", 
                  stringsAsFactors = F, row.names = 1) %>% as.matrix()
sampleMeta = read.csv("../Data_out/sampleMeta.csv", stringsAsFactors = F, row.names = 1)

all(colnames(xcellWork) == rownames(sampleMeta))

xcellRN = xcellWork[,sampleMeta$sample[sampleMeta$class %in% c("AR", "Normal")]]
sigCells = rownames(pxcell)[rowMedians(pxcell) < 0.2]
xcellRN = xcellRN[rownames(xcellRN) %in% sigCells, ]

responseRN = ifelse(sampleMeta[colnames(xcellRN), "class"] == "AR", 1, 0)

sam.resCT = sam(data = xcellRN,
                cl = responseRN, 
                method = wilc.stat, 
                gene.names = rownames(xcellRN), 
                rand = 309,
                R.unlog = FALSE
                #B = 1000
)

delta = findDelta(sam.resCT, 0.05)[2,1]
sam.resCT
saveRDS(sam.resCT, "../Data_out/SAM_cells_RN_results.rds")

sam.sum = summary(sam.resCT, 0.1, bonf = TRUE)
sam.mat = sam.sum@mat.sig
sam.matCT = sam.mat[order(-sam.mat$R.fold),]

write.csv(sam.mat, file = "../Data_out/sigCellsRN.csv")
# sam.mat = read.csv("../Data_out/sigCellsRN.csv", header = TRUE, stringsAsFactors = FALSE, row.names = 1)

cellsRN = rownames(sam.matCT[sam.matCT$q.value <= 0.05,])
cellsRNup = rownames(sam.matCT[sam.matCT$q.value <= 0.05 & sam.matCT$R.fold > 1,])
cellsRNdown = rownames(sam.matCT[sam.matCT$q.value <= 0.05 & sam.matCT$R.fold < 1,])

# write.csv(cellsRN, "Data/kidney/cellsRN.csv", row.names = FALSE)
# write.csv(cellsRNup, "Data/kidney/cellsRNup.csv", row.names = FALSE)
# write.csv(cellsRNdown, "Data/kidney/cellsRNdown.csv", row.names = FALSE)


plotHeatmap(xcellRN[cellsRN,], sampleMeta, rowNames = T)



# Identifying cluster significance ----------------------------------------

shc_result <- shc(x = scale(t(xcellRN[cellsRN,])), 
                  metric = "euclidean", 
                  linkage = "ward", 
                  rcpp = TRUE,
                  n_sim = 1000,
                  n_min = 50,
                  icovest = 1,
                  ci_emp = F,
                  alpha = 0.05)
plot(shc_result, 
     hang = -1, 
     alpha = 0.05, 
     group = sampleFull[colnames(xcellRN), "class"], 
     use_labs = FALSE,
     fwer = TRUE,
     ci_emp = FALSE)

