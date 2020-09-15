source("functions")

sampleMeta = read.csv("../Data_out/sampleMetaAll.csv", stringsAsFactors = F, row.names = 1)
esetList = sapply(unique(sampleMeta$study), function(x) readRDS(paste0("../Data_in/", x, ".rds")))


######## Automatic identification #-----
doppe.res <- doppelgangR(esetList,
                         corFinder.args = list(use.ComBat = F),
                         phenoFinder.args = NULL,
                         cache.dir = NULL,
                         impute.knn.args = NULL,
                         intermediate.pruning = TRUE,
                         within.datasets.only = FALSE,
                         BPPARAM = MulticoreParam(workers = detectCores() - 1))

doppe = summary(doppe.res)

doppe[, c("set1", "GSM_set1")] = stringr::str_split_fixed(doppe$sample1, ":", n = 2)
doppe[, c("set2", "GSM_set2")] = stringr::str_split_fixed(doppe$sample2, ":", n = 2)

doppe = doppe[, c("set1", "set2", "GSM_set1", "GSM_set2", "expr.similarity")]
doppe = doppe[order(doppe$set1),]
doppe$class1 = sampleMeta[doppe$GSM_set1, "subclass"]
doppe$class2 = sampleMeta[doppe$GSM_set2, "subclass"]

doppe$isSamePlatform = sampleMeta[doppe$GSM_set1, "gpl"] == sampleMeta[doppe$GSM_set2, "gpl"]

doppe$isSameDiag = (doppe$class1 == doppe$class2) |
  ((!doppe$class1 %in% c("STA", "Normal")) & (!doppe$class2 %in% c("STA", "Normal")))

duplMeta = sampleMeta[as.vector(unlist(t(doppe[, c("GSM_set1", "GSM_set2")]))),]
duplMeta$sets = rep(1:2, nrow(doppe))

write_csv(duplMeta, "../Data_out/dupl_meta.csv")
write.csv(doppe, "../Data_out/doppe_res.csv")



