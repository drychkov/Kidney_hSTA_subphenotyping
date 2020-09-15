source("functions.R")

genesUp = read.csv(file = "../Data_out/SAM_genesRNup_FC15.csv", stringsAsFactors = F, row.names = 1) %>% rownames()
genesDown = read.csv(file = "../Data_out/SAM_genesRNdown_FC15.csv", stringsAsFactors = F, row.names = 1) %>% rownames()


# genesUp = rownames(gRNup)
# geneIds = bitr(genesUp, fromType = "SYMBOL",
#                toType = "ENTREZID",
#                OrgDb = "org.Hs.eg.db")
# names(genesUp) = geneIds$ENTREZID
# 
# genesDown = rownames(gRNdown)
# geneIds = bitr(genesDown, fromType = "SYMBOL",
#                toType = "ENTREZID",
#                OrgDb = "org.Hs.eg.db")
# names(genesDown) = geneIds$ENTREZID


egoUp <- enrichGO(gene          = symbol2geneID(genesUp),
                  universe      = bitr(rownames(ematWork), fromType = "SYMBOL",
                                       toType = "ENTREZID",
                                       OrgDb = "org.Hs.eg.db")$ENTREZID,
                  OrgDb         = org.Hs.eg.db,
                  keyType       = "ENTREZID",
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  readable      = TRUE)


pdf(file = "../Plots/genesUp_GO_pathways.pdf")
# barplot(egoUp, showCategory = 20)
dotplot(egoUp, showCategory = 20, color = "qvalue")
dev.off()



egoDown <- enrichGO(gene          = symbol2geneID(genesDown),
                    universe      = bitr(rownames(ematWork), fromType = "SYMBOL",
                                         toType = "ENTREZID",
                                         OrgDb = "org.Hs.eg.db")$ENTREZID,
                    OrgDb         = org.Hs.eg.db,
                    keyType       = "ENTREZID",
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    readable      = TRUE)

pdf(file = "../Plots/genesDown_GO_pathways.pdf")
# barplot(egoDown, showCategory = 20)
dotplot(egoDown, showCategory = 20, color = "qvalue")
dev.off()

rm(egoUp, egoDown, genesUp, genesDown)
