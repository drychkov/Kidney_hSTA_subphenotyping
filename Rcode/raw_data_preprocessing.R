source("package_install_load.R")


########### Affymetrix datasets #-----

affyStudies = data.frame(stringsAsFactors = FALSE,
                         study = c("GSE11166", "GSE14328", "GSE21374", "GSE22459", 
                                   "GSE25902", "GSE30718", "GSE34437", "GSE34748", 
                                   "GSE36059", "GSE44131", "GSE48581", "GSE50058",
                                   "GSE50084", "GSE53605", "GSE53769", "GSE54888", 
                                   "GSE57387", "GSE72925", "GSE7392", "GSE76882", "GSE9493"),
                         platform = c("hgu133plus2hsentrezgprobe", "hgu133plus2hsentrezgprobe", 
                                      "hgu133plus2hsentrezgprobe", "hgu133plus2hsentrezgprobe",
                                      "hgu133plus2hsentrezgprobe", "hgu133plus2hsentrezgprobe",
                                      "hgu133plus2hsentrezgprobe", "hgu133plus2hsentrezgprobe", 
                                      "hgu133plus2hsentrezgprobe", "hugene10sthsentrezgprobe", 
                                      "hgu133plus2hsentrezgprobe", "hgu133plus2hsentrezgprobe", 
                                      "hugene10sthsentrezgprobe", "hgu133a2hsentrezgprobe", 
                                      "hugene20sthsentrezgprobe", "hugene10sthsentrezgprobe", 
                                      "huex10sthsentrezgprobe", "hgu133plus2hsentrezgprobe", 
                                      "hgu133plus2hsentrezgprobe","hthgu133pluspmhsentrezgprobe", 
                                      "hgu133plus2hsentrezgprobe"))

cl = makeCluster(spec = detectCores() - 1, type = "FORK")
doParallel::registerDoParallel(cores = cl)
foreach(studyName = affyStudies$study) %do% {
  cat("Processing ", studyName, "...", "\n", sep = "")
  
  celFilePath = file.path("../Data_in/GEO/", studyName, "*.cel.gz")
  pkgName = affyStudies[studyName, "platform"]
  library(package = pkgName, character.only = TRUE, quietly = TRUE)
  
  normSet = SCAN(celFilePath, probeSummaryPackage = pkgName)
  
  saveRDS(normSet, file = paste0("../Data_out/GEO/", studyName, ".rds"))
}
stopCluster(cl)

rm(affyStudies, cl, celFilePath, pkgName, studyName, normSet)

########### Illumina datasets 


ensembl = useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
# View(listDatasets(ensembl))
# View(listFilters(ensembl))
# q = listAttributes(ensembl)
# View(q[grep("affy|agilent|illumina", q$name),])

####### Study GSE43974 #-----
q = read.ilmn(files = "GSE43974_non-normalized_data.txt.gz", path = "../Data_in/GEO/",
              probeid = "PROBE_ID", expr = "AVG_Signal", other.columns = "Detection")
y = neqc(q)
mapping = getBM(attributes = c("illumina_humanht_12_v4", "entrezgene"), 
                mart = ensembl,
                verbose = F,
                uniqueRows = T)
mapping = getGeneProbeMappingDirect(featureDf = mapping, 
                                    geneColname = "entrezgene", 
                                    probeColname = "illumina_humanht_12_v4")

exprsByGene = calcExprsByGeneEmat(y$E, mapping)
# sum(is.na(exprsByGene))
eset <- new("ExpressionSet", exprs = exprsByGene)

esetGEO = getGEO("GSE43974")
colnames(eset) = esetGEO$geo_accession[match(colnames(eset), esetGEO$description)]

saveRDS(eset, "../Data_out/GEO/GSE43974.rds")
rm(q, y, mapping, exprsByGene, eset, esetGEO)

######## Study GSE47097 #-----
q = read.ilmn(files="GSE47097_non_normalized.txt.gz", path = "../Data_in/GEO/",
              probeid = "ID_REF", expr = "SAMPLE", other.columns = "Detection")
y = neqc(q)

esetGEO = getGEO("GSE47097")
mapping = getBM(attributes = c("illumina_humanref_8_v3", "entrezgene"), 
                mart = ensembl,
                verbose = F,
                uniqueRows = T)
mapping = getGeneProbeMappingDirect(featureDf = mapping, 
                                    geneColname = "entrezgene", 
                                    probeColname = "illumina_humanref_8_v3")
exprsByGene = calcExprsByGeneEmat(y$E, mapping)
# sum(is.na(exprsByGene))
eset <- new("ExpressionSet", exprs = exprsByGene)
colnames(eset) = esetGEO$geo_accession[match(paste0("SAMPLE", colnames(q)), esetGEO$description)]

saveRDS(eset, "../Data_out/GEO/GSE47097.rds")
rm(q, y, mapping, exprsByGene, eset, esetGEO)

######### Study GSE52694 #-----
q = read.ilmn(files = "GSE52694_non_normalized.txt.gz", path = '../Data_in/GEO/',
              probeid = "ID_REF", expr = "BORDERLINE_PRAGUE_")
y = backgroundCorrect.matrix(q$E, method = "normexp", normexp.method = "mle", offset = 16)
y = log2(y)
y = normalizeQuantiles(y)

mapping = getBM(attributes = c("illumina_humanht_12_v4", "entrezgene"), 
                mart = ensembl,
                verbose = F,
                uniqueRows = T)
mapping = getGeneProbeMappingDirect(featureDf = mapping, 
                                    geneColname = "entrezgene", 
                                    probeColname = "illumina_humanht_12_v4")
exprsByGene = calcExprsByGeneEmat(y, mapping)
eset <- new("ExpressionSet", exprs = exprsByGene)

esetGEO = getGEO("GSE52694")
colnames(eset) = esetGEO$geo_accession[match(paste0("BORDERLINE_PRAGUE_", colnames(eset)), 
                                             as.character(esetGEO$description))]

saveRDS(eset, "../Data_out/GEO/GSE52694.rds")
rm(q, y, mapping, exprsByGene, eset, esetGEO)

####### Study GSE65326 #-----
q = read.ilmn(files = "GSE65326_non-normalized.txt.gz", path = '../Data_in/GEO/',
              probeid = "ID_REF", expr = "57",  other.columns = "Detection")
y = neqc(q)
mapping = getBM(attributes = c("illumina_humanht_12_v4", "entrezgene"), 
                mart = ensembl,
                verbose = F,
                uniqueRows = T)
mapping = getGeneProbeMappingDirect(featureDf = mapping, 
                                    geneColname = "entrezgene", 
                                    probeColname = "illumina_humanht_12_v4")
exprsByGene = calcExprsByGeneEmat(y$E, mapping)
# sum(is.na(exprsByGene))
exprsByGene = exprsByGene[complete.cases(exprsByGene),]
eset <- new("ExpressionSet", exprs = exprsByGene)

esetGEO = getGEO("GSE65326")
colnames(eset) = esetGEO$geo_accession[match(paste0("57", colnames(eset)), esetGEO$title)]

saveRDS(eset, "../Data_out/GEO/GSE65326.rds")
rm(q, y, mapping, exprsByGene, eset, esetGEO)

###### Study GSE69677 #----
q = read.ilmn(files = "GSE69677_non-normalized.txt.gz", path = '../Data_in/GEO/',
              probeid = "PROBE_ID", expr = "AVG_Signal",  other.columns = "Detection")
y = neqc(q)
mapping = read.csv("../Data_in/GEO/HumanHT-12_V4_0_R2_15002873_B_WGDASL.txt.TXT", 
                   skip = 8, sep = "\t")
mapping = getGeneProbeMappingDirect(featureDf = mapping, 
                                    geneColname = "Entrez_Gene_ID", 
                                    probeColname = "Probe_Id")
exprsByGene = calcExprsByGeneEmat(y$E, mapping)
# sum(is.na(exprsByGene))
exprsByGene = exprsByGene[complete.cases(exprsByGene),]

eset <- new("ExpressionSet", exprs = exprsByGene)
esetGEO = getGEO(filename = "../Data_in/GEO/GSE69677_series_matrix.txt.gz")
eset = eset[,colnames(eset)[(colnames(eset) %in% esetGEO$description)]]
colnames(eset) = esetGEO$geo_accession[match(colnames(eset), esetGEO$description)]

saveRDS(eset, "../Data_out/GEO/GSE69677.rds")
rm(q, y, mapping, exprsByGene, eset, esetGEO)

############# Agilent datasets


###### Study GSE60807 #-----
studyName = "GSE60807"

q = read.maimages(files = dir("../Data_in/GEO/GSE60807"), 
                  path = "../Data_in/GEO/GSE60807", source = "agilent", green.only = TRUE)
colnames(q) = unname(fixGeoSampleNames(colnames(q)))

y = neqc(q, status = q$genes$ControlType, negctrl = -1, regular = 0)
rownames(y$E) = y$genes$ProbeName

mapping = getBM(attributes = c("agilent_wholegenome_4x44k_v1", "entrezgene"), 
                mart = ensembl,
                verbose = F,
                uniqueRows = T)
mapping = getGeneProbeMappingDirect(featureDf = mapping, 
                                    geneColname = "entrezgene", 
                                    probeColname = "agilent_wholegenome_4x44k_v1")

exprsByGene = calcExprsByGeneEmat(y$E, mapping)
# sum(is.na(exprsByGene))
eset <- new("ExpressionSet", exprs = exprsByGene)

saveRDS(eset, "../Data_out/GEO/GSE60807.rds")
rm(q, y, mapping, exprsByGene, eset)

######### Study GSE10419 #-------

q = getGEO(filename = "../Data_in/GEO/GSE10419_family.soft.gz")
a = foreach(name = names(q@gsms)) %do% {emat = q@gsms[name]}[[1]]
b = lapply(a, function(x) q@gsms)[[1]]
c = lapply(b, function(x) x@dataTable@table[,6])
y = do.call(cbind, c)

esetGEO = getGEO("GSE10419")[[1]]
rownames(y) = esetGEO@featureData@data$SPOT_ID

emat = backgroundCorrect.matrix(y, method = "normexp", normexp.method = 'mle', offset = 16)
emat = log2(emat)
emat = normalizeQuantiles(emat)

mapping = esetGEO@featureData@data[, c("SPOT_ID", "GENE")]
mapping = getGeneProbeMappingDirect(featureDf = mapping, 
                                    geneColname = "GENE", 
                                    probeColname = "SPOT_ID")
exprsByGene = calcExprsByGeneEmat(emat, mapping)
sum(is.na(exprsByGene))
eset <- new("ExpressionSet", exprs = exprsByGene)

saveRDS(eset, "../Data_out/GEO/GSE10419.rds")
rm(q, a, b, c, y, mapping, emat, exprsByGene, eset)


rm(ensembl)
