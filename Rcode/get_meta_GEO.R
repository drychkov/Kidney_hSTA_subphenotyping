source("functions.R")

# Download database if not existed
if(!file.exists('GEOmetadb.sqlite')) getSQLiteFile(destdir = "../Data_in/")

# Using SQL to extract data from GSE
getGSEInfo <- function(gse_list_char){
  sqlQuery <- paste("SELECT DISTINCT gse, title, pubmed_id, summary, overall_design, contributor",
                    "FROM",
                    "gse" ,
                    "WHERE",
                    "gse IN (", gse_list_char, ")", sep=" ")
  rs <- dbGetQuery(con, sqlQuery)
  return(rs)
}


# Using SQL to extract annotations from GSM
getGSMInfo <- function(gse_list_char){
  sqlQuery <- paste("SELECT DISTINCT",
                    "*",
                    "FROM",
                    "  gse_gsm join gsm on gse_gsm.gsm=gsm.gsm",
                    "WHERE",
                    "  gse_gsm.gse IN (", gse_list_char, ")",
                    sep=" ")
  rs <- dbGetQuery(con,sqlQuery)
  return(rs)
}

# Connect to database
con <- dbConnect(SQLite(),'../Data_in/GEOmetadb.sqlite')

gse_list_char = "'GSE10419', 'GSE11166', 'GSE14328', 'GSE21374', 'GSE22459', 'GSE25902', 'GSE30718', 'GSE34437', 'GSE34748', 'GSE36059', 'GSE43974', 'GSE44131', 'GSE47097', 'GSE48581', 'GSE50058', 'GSE50084', 'GSE52694', 'GSE53605', 'GSE53769', 'GSE54888', 'GSE57387', 'GSE60807', 'GSE65326', 'GSE69677', 'GSE72925', 'GSE7392', 'GSE76882', 'GSE9493'"

gsm_info <- getGSMInfo(gse_list_char)
write.csv(gsm_info, '../Data_out/gsm_info.csv')
dbDisconnect(con)

rm(con, gsm_info, gse_list_char, getGSMInfo)


# A further (semi-)manual revision for better sample annotation is required for metadata in sampleMeta.csv
# The resulting metadata should have at least these columns: study, sample, subclass, gpl, batch.

sampleMeta = read.csv('../Data_out/gsm_info.csv', stringsAsFactors = FALSE)

detect = stringr::str_detect(unname(unlist(anno[,'subclass'])), "ABMR|TCMR|Mixed|BL|AR[+]CAN|AR|BL[+]CAN|STA|Normal")
sampleMeta = sampleMeta[detect,]

sampleMeta$class = case_when(
  sampleMeta$subclass %in% c("ABMR","TCMR", "Mixed","AR+CAN", "BL", "BL+CAN", "AR") ~ "AR",
  sampleMeta$subclass == "STA" ~ "STA",
  sampleMeta$subclass == "Normal" ~ "Normal"
)
rownames(sampleMeta) = sampleMeta$sample
write.csv(sampleMeta, file = "../Data_out/sampleMetaAll.csv")


## Cell file names cleanup
dirz = list.dirs('../Data_in/GEO/')
for (onedir in dirz) {
  filez = list.celfiles(onedir)
  cwd = setwd(onedir)
  for (onefile in filez) {
    newfile = paste0(unlist(strsplit(onefile, split="[._]"))[1], ".cel.gz")
    file.rename(onefile, newfile)
  }
  setwd(cwd)
}
rm(dirz, filez, cwd, newfile)