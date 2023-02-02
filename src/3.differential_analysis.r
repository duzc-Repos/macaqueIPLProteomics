setwd("C:/Users/86188/Desktop/work/brainnetome/project/help/macaqueIPLProteomics/src/")

library(readxl)
library(DEqMS)
library(ggplot2)
library(clusterProfiler)
library(org.Mmu.eg.db)

#######################################
# 0. Load raw proteomics data (after PD processed)
#######################################

expr.log <- read.table("../res/0.preprocessing/sample_expr_rbe.tsv")
rawinfo <- read.table("../res/0.preprocessing/sample_info.tsv ", header = T, row.names = 1)
rawinfo$Group <- gsub("-", "_", rawinfo$Group)

################################################
## differential analysis (Part1)
#  1. Pairwise model   (one vs. one)
#  2. Enrichment model (one vs. all)
################################################

design <- model.matrix(~0 + as.factor(rawinfo$Group))
colnames(design) <- sapply(strsplit(colnames(design), split=")", fixed=T), function(x){x[2]})
colnames(design) <- gsub("-", "_", colnames(design), fixed=T)
rownames(design) <- rownames(rawinfo)
fit <- lmFit(expr.log, design = design)
eb <- eBayes(fit)

## 1. pairwise model
region_combs <- combn(colnames(design), 2)
region_contrasts <- apply(region_combs, 2, function(x) {
  z <- paste(x, collapse = '-')
  makeContrasts(contrasts = z, levels = design)
})
rownames(region_contrasts) <- colnames(design)
colnames(region_contrasts) <- apply(region_combs, 2, paste, collapse="-")

eb_contrasts <- eBayes(contrasts.fit(fit, region_contrasts))
pvals_contrasts <- eb_contrasts$p.value
fdr_contrasts <- apply(pvals_contrasts, 2, p.adjust, 'BH')
pairwise_model_conuting <- data.frame("FDR" = colSums(fdr_contrasts < 0.05),
                                      "Pval-2" = colSums(pvals_contrasts < 0.01),
                                      "Pval-6" = colSums(pvals_contrasts < 1e-6) )
pairwise_count_mat = as.data.frame(matrix(0, nrow=ncol(design), ncol = ncol(design)))
idx = c("L_PF", "L_PFG", "L_PFop", "L_PGop", "L_PG", "L_Opt",
        "R_PF", "R_PFG", "R_PFop", "R_PGop", "R_PG", "R_Opt")
colnames(pairwise_count_mat) <- idx
rownames(pairwise_count_mat) <- idx
for(i in rownames(pairwise_model_conuting)){
  g = unlist(strsplit(i, split="-", fixed=T))
  pairwise_count_mat[g[1], g[2]] = pairwise_model_conuting[i, 1]
  pairwise_count_mat[g[2], g[1]] = pairwise_model_conuting[i, 1]
}
write.table(pairwise_count_mat, file="../res/3.diff/pairwise_counting.txt", quote = F, row.names = T, col.names = T)


## 2. enrichment model
eb0_list <- lapply(idx, function(x){
  cat(x, "\n")
  res <- as.factor(rawinfo$Group == x)
  m <- model.matrix(~res)
  tmplm <- lmFit(expr.log, design = m)
  eBayes(tmplm)
})
pvals0_contrasts <- sapply(eb0_list, function(x){ x$p.value[, 2, drop=FALSE] })
fdr0_contrasts <- apply(pvals0_contrasts, 2, p.adjust, 'BH')
tstats0_contrasts <- sapply(eb0_list, function(x) {x$t[, 2, drop = FALSE]})
enrichment_model_counting <- data.frame("FDR" = colSums(fdr0_contrasts < 0.05),
                                          "Pval-2" = colSums(pvals0_contrasts < 0.01),
                                          "Pval-6" = colSums(pvals0_contrasts < 1e-6))
rownames(enrichment_model_counting) <- idx

# Fecth Gene List
de_summary <- list()
pairwise_genelist <- c()
for (j in 1:nrow(pairwise_model_conuting)){
  tmp <- topTable(eb_contrasts, coef = j, number = Inf)[c(1, 4, 5)]
  fdrsig <- tmp$adj.P.Val < 0.05
  tmp$reg <- "No"
  tmp$reg[fdrsig & tmp$logFC>0] = "Up"
  tmp$reg[fdrsig & tmp$logFC<0] = "Down"
  de_summary[[rownames(pairwise_model_conuting)[j]]] <- tmp
  pairwise_genelist <- append(pairwise_genelist, rownames(tmp)[tmp$reg != "No"], length(pairwise_genelist))
}
pairwise_genelist <- unique(pairwise_genelist)

enrichment_genelist <- c()
for (j in 1:nrow(enrichment_model_counting)){
  tmp <- topTable(eb0_list[[j]], number = Inf)[c(1, 4, 5)]
  fdrsig <- tmp$adj.P.Val < 0.05
  tmp$reg <- "No"
  tmp$reg[fdrsig & tmp$logFC>0] = "Up"
  de_summary[[rownames(enrichment_model_counting)[j]]] <- tmp
  enrichment_genelist <- append(enrichment_genelist, rownames(tmp)[tmp$reg != "No"], length(enrichment_genelist))
  
  enrichment_model_counting[j, 1] <- sum(tmp$reg == "Up")
  enrichment_model_counting[j, 2] <- sum(tmp$P.Value<0.01 & tmp$logFC>0)
  enrichment_model_counting[j, 3] <- sum(tmp$P.Value<1e-6 & tmp$logFC>0)
}
enrichment_genelist <- unique(enrichment_genelist)

write.table(enrichment_model_counting, file="../res/3.diff/enrichment_counting.txt", quote = F, row.names = T, col.names = F)
combined_genelist <- unique(append(pairwise_genelist, enrichment_genelist, length(pairwise_genelist)))
write.table(combined_genelist, file="../res/3.diff/combined_genelist.txt", quote = F, row.names = F, col.names = F)

for (j in names(de_summary) ){
  file_name <- paste("../res/3.diff/", j, ".tsv", sep="")
  write.table(de_summary[[j]], file_name, quote = F, sep = "\t")
}



##############################################
# heatmap
##############################################
diff_pro_expr <- expr.log[combined_genelist, ]
colnames(diff_pro_expr) <- gsub(".", "-", gsub("X", "INS ", colnames(diff_pro_expr)), fixed = T)

anno_col = data.frame(row.names = colnames(diff_pro_expr), 
                      "Hemi"= rawinfo$Hemi,
                      "Region"=factor(rawinfo$Region, levels = c("PF", "PFop",  
                                                                "PFG", "PGop",
                                                                "PG",  "Opt")))

a <- factor(rawinfo$Group, levels = c("L_PF", "R_PF", "L_PFop", "R_PFop", 
                                 "L_PFG", "R_PFG", "L_PGop", "R_PGop",
                                 "L_PG", "R_PG", "L_Opt", "R_Opt"))
b <- factor(rawinfo$Group, levels = c("L_PF", "L_PFop", "L_PFG", "L_PGop", "L_PG", "L_Opt",
                                      "R_PF", "R_PFop", "R_PFG", "R_PGop", "R_PG",  "R_Opt"))
tmp <- pheatmap(as.matrix(diff_pro_expr[, order(b)]), 
         cluster_rows = T, cutree_rows=5,
         cluster_cols = F, 
         scale = "row", 
         treeheight_col=20, #treeheight_row=0,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "euclidean",
         show_rownames = F, show_colnames = F,
         annotation_col = anno_col,
         clustering_method = "complete",
         #annotation_colors = anno_col_color[0]
         )

cut_tree = cutree(tmp$tree_row, k=5)

pheatmap(as.matrix(diff_pro_expr[cut_tree==1, order(b)]), 
                cluster_rows = T,
                cluster_cols = F, 
                scale = "row", 
                treeheight_col=20, #treeheight_row=0,
                clustering_distance_rows = "correlation",
                clustering_distance_cols = "euclidean",
                show_rownames = F, show_colnames = F,
                annotation_col = anno_col,
                clustering_method = "complete",
                #annotation_colors = anno_col_color[0]
)

cluster1_genelist <- rownames(diff_pro_expr[cut_tree==1, ])
cluster1_go <- enrichGO(cluster1_genelist, OrgDb = org.Mmu.eg.db, keyType = "SYMBOL", 
                       ont = "BP", #universe = rownames(diff_pro_expr),
                       pvalueCutoff = 0.05, qvalueCutoff = 0.05)
dotplot(cluster1_go, color = "p.adjust", showCategory=10)
cnetplot(cluster1_go, showCategory = 10)

cluster1_kegg <- enrichKEGG (bitr(cluster1_genelist, fromType = "SYMBOL", toType = "UNIPROT", OrgDb = org.Mmu.eg.db)$UNIPROT, 
                           organism="mcc", keyType = "uniprot", 
                           pvalueCutoff = 0.05, qvalueCutoff = 0.05)
dotplot(cluster1_kegg, color = "p.adjust", showCategory=10)
cnetplot(cluster1_kegg, showCategory = 10)



pheatmap(as.matrix(diff_pro_expr[cut_tree==4, order(b)]), 
         cluster_rows = T,
         cluster_cols = F, 
         scale = "row", 
         treeheight_col=20, #treeheight_row=0,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "euclidean",
         show_rownames = F, show_colnames = F,
         annotation_col = anno_col,
         clustering_method = "complete",
         #annotation_colors = anno_col_color[0]
)

cluster2_genelist <- rownames(diff_pro_expr[cut_tree==4, ])
cluster2_go <- enrichGO(cluster2_genelist, OrgDb = org.Mmu.eg.db, keyType = "SYMBOL", 
                        ont = "BP", 
                        pvalueCutoff = 0.05, qvalueCutoff = 0.05)
dotplot(cluster2_go, color = "p.adjust", showCategory=10)
cnetplot(cluster2_go, showCategory = 10)

cluster2_kegg <- enrichKEGG (bitr(cluster2_genelist, fromType = "SYMBOL", toType = "UNIPROT", OrgDb = org.Mmu.eg.db)$UNIPROT, 
                             organism="mcc", keyType = "uniprot", 
                             pvalueCutoff = 0.05, qvalueCutoff = 0.05)
dotplot(cluster2_kegg, color = "p.adjust", showCategory=10)


pheatmap(as.matrix(diff_pro_expr[cut_tree==3, order(b)]), 
         cluster_rows = T,
         cluster_cols = F, 
         scale = "row", 
         treeheight_col=20, #treeheight_row=0,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "euclidean",
         show_rownames = F, show_colnames = F,
         annotation_col = anno_col,
         clustering_method = "complete",
         #annotation_colors = anno_col_color[0]
)

cluster3_genelist <- rownames(diff_pro_expr[cut_tree==3, ])
cluster3_go <- enrichGO(cluster3_genelist, OrgDb = org.Mmu.eg.db, keyType = "SYMBOL", 
                        ont = "BP", 
                        pvalueCutoff = 0.05, qvalueCutoff = 0.05)
dotplot(cluster2_go, color = "p.adjust", showCategory=10)
cnetplot(cluster2_go, showCategory = 10)

















