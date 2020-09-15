source("functions.R")

ematWork = readRDS(file = "../Data_out/ematMerged_BatchCorr.rds")

cellTypes = read.csv(file = '../Data_in/cellTypeNames.csv', stringsAsFactors = FALSE, header = TRUE, na.strings = NA)
cellTypes = as.list(cellTypes)
cellTypes = lapply(cellTypes, function(x) x[x != ""])
# cellList = c(unname(unlist(cellTypes)), c("ImmuneScore", "StromaScore", "MicroenvironmentScore"))
cellList = unname(unlist(cellTypes))

xcellWork = xCellAnalysis(expr = ematWork, 
                          rnaseq = FALSE, 
                          file.name = "../Data_out/xcellScores.csv", 
                          cell.types.use = cellList, # Find enrichments for cell types of interest
                          parallel.sz = detectCores() - 1,
                          parallel.type = "FORK")

pxcell = xCellSignifcanceBetaDist(xcellWork, rnaseq = FALSE)
colnames(pxcell) = colnames(xcellWork)
write.csv(pxcell, "../Data_out/xcell_pvalues.csv")


rm(cellTypes, xcellWork, pxcell)