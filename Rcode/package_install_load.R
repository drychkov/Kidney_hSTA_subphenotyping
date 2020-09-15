# files were obtained from http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/22.0.0/entrezg.asp
brainArrayPckgs = c("hgu133plus2hsentrezgprobe", "hugene10sthsentrezgprobe", "hgu133a2hsentrezgprobe", "hugene20sthsentrezgprobe", "huex10sthsentrezgprobe", "hthgu133pluspmhsentrezgprobe")
q = brainArrayPckgs[!brainArrayPckgs %in% rownames(installed.packages())]
install.packages(pkgs = q, repos = NULL, type = 'source')


pkgs = c('GEOquery', 'GEOmetadb', 'affy', 'oligo', 'annotate', 'org.Hs.eg.db', 'foreach', 'impute', 'tidyverse', 'sva', 'devtools', 'limma', 'biomaRt', 'doParallel', 'SCAN.UPC', 'ggsci', 'reshape2', 'ROCR', 'cba', 'pheatmap',  'ggplot2',  'gridExtra',  'RColorBrewer', 'DBI', 'RSQLite', 'doppelgangR', 'BiocParallel', 'clusterProfiler', 'ggbeeswarm', 'readxl'
         #'ReactomePA', 
         # 'Rtsne',
         # 'data.table',
         # "RUVcorr",
         # "RUVnormalize",
         # "convert",
         # 'affyQCReport',
         # 'arrayQualityMetrics',
         # 'preprocessCore'
)

if (!require("BiocManager")) install.packages("BiocManager")
pkg2install = pkgs[!pkgs %in% rownames(installed.packages())]
BiocManager::install(pkg2install)
invisible(lapply(pkgs, library, character.only = TRUE))

if (!require("xCell")) BiocManager::install('dviraran/xCell')
library("xCell")
if (!require("feseR")) BiocManager::install('enriquea/feseR')
library("feseR")
if (!require("sigclust2")) BiocManager::install('pkimes/sigclust2')
library("sigclust2")
if (!require("Rclusterpp")) BiocManager::install('nolanlab/Rclusterpp')
library("Rclusterpp")


# for (pkg in pkgs) {
#   q = require(pkg, character.only = TRUE)
#   if (!q) BiocManager::install(pkg)
# }

# if (!require("pacman")) install.packages("pacman"); library(pacman)
# pacman::p_load(char = pkgs, install = T)
# rm(pkgs)

rm(brainArrayPckgs, pkgs, pkg2install)






