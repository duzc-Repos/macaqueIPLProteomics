library(CellChat)
library(homologene)
library(gprofiler2)
options(stringsAsFactors = FALSE)

CellChatDB <- CellChatDB.human
write.table(CellChatDB$interaction, file="../res/gprofiler2/cellchatdb_human.tsv", sep="\t", quote=F, row.names=F, col.names=T)
symbol_mapper_human2macaque <- gorth( CellChatDB$geneInfo$Symbol, source_organism="hsapiens", target_organism="mmulatta" )[, c("input", "ortholog_name")]
write.table(symbol_mapper_human2macaque, file=".../res/gprofiler2/human2macaque.tsv", sep="\t", quote=F, row.names=F, col.names=T)

CellChatDB <- CellChatDB.mouse
write.table(CellChatDB$interaction, file="../res/gprofiler2/cellchatdb_mouse.tsv", sep="\t", quote=F, row.names=F, col.names=T)
symbol_mapper_mouse2macaque <- gorth( CellChatDB$geneInfo$Symbol, source_organism="mmusculus", target_organism="mmulatta" )[, c("input", "ortholog_name")]
write.table(symbol_mapper_mouse2macaque, file=".../res/gprofiler2/mouse2macaque.tsv", sep="\t", quote=F, row.names=F, col.names=T)
