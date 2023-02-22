setwd("C:/Users/86188/Desktop/work/brainnetome/project/help/macaqueIPLProteomics/src/")

library(ggplot2)
library(clusterProfiler)
library(org.Mmu.eg.db)


########################################
# Left Hemi
########################################
gene_sdi_L <- read.csv("../res/5.gredient_sdi_cellType/gene_sdi_L.tsv", header = F)
head_n = round(nrow(gene_sdi_L) * 0.25)
tail_n = round(nrow(gene_sdi_L) * 0.75)

decoupled_genelist <- gene_sdi_L$V1[1:head_n]
coupled_genelist <- gene_sdi_L$V1[tail_n:nrow(gene_sdi_L)]

## coupled genelist
coupled_go <- enrichGO(coupled_genelist, OrgDb = org.Mmu.eg.db, keyType = "SYMBOL", 
                       ont = "BP", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
dotplot(coupled_go, color = "p.adjust", showCategory=15, title="Left Hemi GO: Coupled")
ggsave("../res/5.gredient_sdi_cellType/coupled_go_L.pdf", width=9, height=5 )

coupled_kegg <- enrichKEGG (bitr(coupled_genelist, fromType = "SYMBOL", toType = "UNIPROT", OrgDb = org.Mmu.eg.db)$UNIPROT, 
                             organism="mcc", keyType = "uniprot",
                             pvalueCutoff = 0.05, qvalueCutoff = 0.05)
dotplot(coupled_kegg, color = "p.adjust", showCategory=15, title="Left Hemi KEGG: coupled")
ggsave("../res/5.gredient_sdi_cellType/coupled_kegg_L.pdf", width=6.5, height=5 )


## decoupled genelist
decoupled_go <- enrichGO(decoupled_genelist, OrgDb = org.Mmu.eg.db, keyType = "SYMBOL", 
                       ont = "BP", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
dotplot(decoupled_go, color = "p.adjust", showCategory=15, title="Left Hemi GO: Decoupled")
ggsave("../res/5.gredient_sdi_cellType/decoupled_go_L.pdf", width=9, height=5 )

decoupled_kegg <- enrichKEGG (bitr(decoupled_genelist, fromType = "SYMBOL", toType = "UNIPROT", OrgDb = org.Mmu.eg.db)$UNIPROT, 
                            organism="mcc", keyType = "uniprot",
                            pvalueCutoff = 0.05, qvalueCutoff = 0.05)
dotplot(decoupled_kegg, color = "p.adjust", showCategory=10, title="Left Hemi KEGG: Deoupled")
ggsave("../res/5.gredient_sdi_cellType/decoupled_kegg_L.pdf", width=6, height=5 )


########################################
# Right Hemi
########################################
gene_sdi_R <- read.csv("../res/5.gredient_sdi_cellType/gene_sdi_R.tsv", header = F)
head_n = round(nrow(gene_sdi_R) * 0.25)
tail_n = round(nrow(gene_sdi_R) * 0.75)

decoupled_genelist <- gene_sdi_R$V1[1:head_n]
coupled_genelist <- gene_sdi_R$V1[tail_n:nrow(gene_sdi_R)]

## coupled genelist
coupled_go <- enrichGO(coupled_genelist, OrgDb = org.Mmu.eg.db, keyType = "SYMBOL", 
                       ont = "BP", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
dotplot(coupled_go, color = "p.adjust", showCategory=15, title="Left Hemi GO: Coupled")
ggsave("../res/5.gredient_sdi_cellType/coupled_go_R.pdf", width=9, height=5 )

coupled_kegg <- enrichKEGG (bitr(coupled_genelist, fromType = "SYMBOL", toType = "UNIPROT", OrgDb = org.Mmu.eg.db)$UNIPROT, 
                            organism="mcc", keyType = "uniprot",
                            pvalueCutoff = 0.05, qvalueCutoff = 0.05)
dotplot(coupled_kegg, color = "p.adjust", showCategory=15, title="Left Hemi KEGG: coupled")
ggsave("../res/5.gredient_sdi_cellType/coupled_kegg_R.pdf", width=6.5, height=5 )


## decoupled genelist
decoupled_go <- enrichGO(decoupled_genelist, OrgDb = org.Mmu.eg.db, keyType = "SYMBOL", 
                         ont = "BP", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
dotplot(decoupled_go, color = "p.adjust", showCategory=15, title="Left Hemi GO: Decoupled")
ggsave("../res/5.gredient_sdi_cellType/decoupled_go_R.pdf", width=9, height=5 )

decoupled_kegg <- enrichKEGG (bitr(decoupled_genelist, fromType = "SYMBOL", toType = "UNIPROT", OrgDb = org.Mmu.eg.db)$UNIPROT, 
                              organism="mcc", keyType = "uniprot",
                              pvalueCutoff = 0.05, qvalueCutoff = 0.05)
dotplot(decoupled_kegg, color = "p.adjust", showCategory=10, title="Left Hemi KEGG: Deoupled")
ggsave("../res/5.gredient_sdi_cellType/decoupled_kegg_R.pdf", width=6, height=5 )


##########################
# Compare analysis
##########################
gene_sdi_L <- read.table("../res/5.gredient_sdi_cellType/gene_sdi_fc_L.tsv", header = F, sep="\t")
head_n = round(nrow(gene_sdi_L) * 0.25)
tail_n = round(nrow(gene_sdi_L) * 0.75)

decoupled_genelist_L <- as.data.frame(bitr(gene_sdi_L$V1[1:head_n], 
                                           fromType = "SYMBOL", toType = "ENTREZID", 
                                           OrgDb = org.Mmu.eg.db))
decoupled_genelist_L$Group = "Decouple"
decoupled_genelist_L$Hemi = "L"
colnames(decoupled_genelist_L)[1:2] = c("Symbol", "Entrez")

coupled_genelist_L <- as.data.frame(bitr(gene_sdi_L$V1[tail_n:nrow(gene_sdi_L)], 
                                         fromType = "SYMBOL", toType = "ENTREZID", 
                                         OrgDb = org.Mmu.eg.db))
coupled_genelist_L$Group = "Couple"
coupled_genelist_L$Hemi = "L"
colnames(coupled_genelist_L)[1:2] = c("Symbol", "Entrez")


gene_sdi_R <- read.table("../res/5.gredient_sdi_cellType/gene_sdi_fc_R.tsv", header = F, sep="\t")
head_n = round(nrow(gene_sdi_R) * 0.25)
tail_n = round(nrow(gene_sdi_R) * 0.75)

decoupled_genelist_R <- as.data.frame(bitr(gene_sdi_R$V1[1:head_n],
                                           fromType = "SYMBOL", toType = "ENTREZID", 
                                           OrgDb = org.Mmu.eg.db))
decoupled_genelist_R$Group = "Decouple"
decoupled_genelist_R$Hemi = "R"
colnames(decoupled_genelist_R)[1:2] = c("Symbol", "Entrez")

coupled_genelist_R <- as.data.frame(bitr(gene_sdi_R$V1[tail_n:nrow(gene_sdi_R)],
                                         fromType = "SYMBOL", toType = "ENTREZID", 
                                         OrgDb = org.Mmu.eg.db))
coupled_genelist_R$Group = "Couple"
coupled_genelist_R$Hemi = "R"
colnames(coupled_genelist_R)[1:2] = c("Symbol", "Entrez")


genelist <- rbind(coupled_genelist_L, coupled_genelist_R, decoupled_genelist_L, 
                  decoupled_genelist_R)

xx <- compareCluster(Entrez~Group+Hemi, data=genelist, 
                     fun="groupGO", ont="BP", level=5, readable=TRUE,
                     OrgDb=org.Mmu.eg.db)
dotplot(xx, showCategory=10)


xx <- compareCluster(Entrez~Group+Hemi, data=genelist, 
                     fun="enrichGO", ont="BP", 
                     readable=TRUE,
                     pvalueCutoff = 0.05, qvalueCutoff = 0.05,
                     OrgDb=org.Mmu.eg.db )
dotplot(xx, showCategory=10)

a <- as.data.frame(xx)



xx <- compareCluster(Entrez~Group+Hemi, data=genelist, 
                     fun="enrichKEGG", 
                     #pvalueCutoff=0.05, 
                     qvalueCutoff = 0.5,
                     organism="mcc" )
dotplot(xx, showCategory=10)

b <- as.data.frame(xx)
write.table(b, "../res/5.gredient_sdi_cellType/compare_cluster_KEGG_FC.tsv", sep="\t", quote=F, row.names = T, col.names = T)










