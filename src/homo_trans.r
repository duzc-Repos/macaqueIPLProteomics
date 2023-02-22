setwd("C:/Users/86188/Desktop/work/brainnetome/project/help/macaqueIPLProteomics/src/")

library(readxl)
library(homologene)

# House keeping genes
hk_genes <- read.table("../data/refdata/HK_genes_human.txt")
hk_genes <- unique(na.omit(homologene(hk_genes$V1, inTax = 9606, outTax = 9544)$`9544`))
hk_genes <- as.data.frame(hk_genes)
hk_genes$Class <- "HKG"
colnames(hk_genes) <- c("Genes", "Class")

# Synaptome
synaptome <- unique(read_excel("../data/refdata/41593_2018_195_MOESM4_ESM.xlsx", col_names = "gene", sheet = "Synaptome"))
synaptome <- unique(na.omit(homologene(synaptome$gene, inTax = 9606, outTax = 9544)$`9544`))
synaptome <- as.data.frame(synaptome)
synaptome$Class <- "Synaptome"
colnames(synaptome) <- c("Genes", "Class")

# EVs
ev_genes <- read.table("../data/refdata/EV_TOP_100.txt", sep="\t", header = T)
ev_genes <- unique(na.omit(homologene(ev_genes$GENE.SYMBOL, inTax = 9606, outTax = 9544)$`9544`))
ev_genes <- as.data.frame(ev_genes)
ev_genes$Class <- "EVs"
colnames(ev_genes) <- c("Genes", "Class")

# LRs
lr_genes <- read.table("../data/refdata/CellTalk_human_lr_pair.txt", header = T)
l_genes <- unique(na.omit(homologene(unique(lr_genes$ligand_gene_symbol), inTax = 9606, outTax = 9544)$`9544`))
l_genes <- as.data.frame(l_genes)
l_genes$Class <- "Ligand"
colnames(l_genes) <- c("Genes", "Class")

r_genes <- unique(na.omit(homologene(unique(lr_genes$receptor_gene_symbol), inTax = 9606, outTax = 9544)$`9544`))
r_genes <- as.data.frame(r_genes)
r_genes$Class <- "Receptor"
colnames(r_genes) <- c("Genes", "Class")


signaling_related_genes <- rbind(hk_genes, synaptome, ev_genes, l_genes, r_genes)
write.table(signaling_related_genes, file = "../data/refdata/signaling_related_genes.tsv", quote = F, row.names = F, col.names = T, sep="\t")

# Layers 
sheet_name = c("Layer2enriched", "Layer3enriched", "Layer4enriched", "Layer5enriched", "Layer6enriched")
layer_genes <- lapply(sheet_name, function(x){
  a <- as.data.frame(unique(na.omit(read_excel("../data/refdata/neuron_11047_mmc6.xlsx", sheet = x)$`Gene Symbol`)))
  a$Class <- strsplit(x, split = "enriched", fixed = T)[[1]][1]
  colnames(a) <- c("Genes", "Class")
  return(a)
})
layer_genes <- Reduce(function(x, y){rbind(x, y)}, layer_genes )
write.table(layer_genes, file = "../data/refdata/layer_genes.tsv", quote = F, row.names = F, col.names = T, sep="\t")




















