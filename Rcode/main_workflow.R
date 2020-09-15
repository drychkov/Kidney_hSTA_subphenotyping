source("functions.R")


esetList = sapply(unique(sampleMeta$study), function(x) readRDS(paste0("../Data_in/", x, ".rds")))


# Get study metadata ------------------------------------------------------
source("get_meta_GEO.R")


# Raw data downloading ----------------------------------------------------

# All data should be downloaded to ../Data_in/GEO
#
# Affymetrix datasets should be in separate folders named by a study name GSE#####
#
# Download non-normalized Illumna files:
# GSE43974_non-normalized_data.txt.gz
# GSE47097_non_normalized.txt.gz
# GSE52694_non_normalized.txt.gz
# GSE65326_non-normalized.txt.gz
# GSE69677_non-normalized.txt.gz
#
# Download Agilent data:
# txt.gz files from GSE60807 study to a separate folder 
# GSE10419_family.soft.gz

# Raw data preprocessing --------------------------------------------------
source("raw_data_preprocessing.R")


# Search for duplicates ---------------------------------------------------
source("duplicates_search.R")


# Remove duplicates from the data -----------------------------------------

# After the manual curation remove duplicates from the metadata and gene expression sets
sampleMeta = read.csv("../Data_out/sampleMetaAll.csv", stringsAsFactors = F, row.names = 1)
duplicates2Remove = read.csv("../Data_out/doppe_res.csv", stringsAsFactors = F)
duplicates2Remove = duplicates2Remove[duplicates2Remove$manual_flag,]

# Removing one outlier
sampleMeta = sampleMeta[!sampleMeta$sample == "GSM1874749",]

esetList = sapply(unique(sampleMeta$study), 
                  function(x) readRDS(paste0("../Data_out/GEO/", x, ".rds")))
for(name in names(esetList)){
  toRemove = duplicates2Remove[duplicates2Remove$set2 == name,"GSM_set2"]
  esetList[[name]] = esetList[[name]][, !colnames(esetList[[name]]) %in% toRemove]
}
saveRDS(esetList, file = "../Data_out/esetListNoDup.rds")

sampleMeta = sampleMeta[!rownames(sampleMeta) %in% duplicates2Remove$GSM_set2,]
write_csv(sampleMeta, "../Data_out/sampleMeta.csv")


# Data merging ------------------------------------------------------------
source("data_merging.R")


# DGE Analysis ------------------------------------------------------------
source("DGE_analysis.R")


# Pathway Analysis --------------------------------------------------------
source("pathway_analysis.R")


# Cell Type Analysis ------------------------------------------------------
source("cell_type_significance_analysis.R")


# Feature selection procedure on genes ------------------------------------
source("FS_genes.R")


# Feature selection procedure on cell types -------------------------------
source("FS_cellTypes.R")


# InstaScore introduction and application to hSTA -------------------------
source("instaScore.R")


# InstaScore Validation ---------------------------------------------------
source("instaScore_validation.R")





