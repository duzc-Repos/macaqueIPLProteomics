setwd("C:/Users/86188/Desktop/work/brainnetome/project/help/macaqueIPLProteomics/src/")

library(mistyR)
library(future)

# data manipulation
library(dplyr)
library(purrr)
library(distances)

# plotting
library(ggplot2)

#plan(multisession)



##################################
# Left Hemi
##################################
sample.expr <- t(read.table("../res/4.SC_exploration/sample_expr_L.tsv", header = TRUE,  sep = "\t",  row.names=1))
colnames(sample.expr) <- gsub("-", "_", colnames(sample.expr), fixed = T)
sample.pos <- read.table("../res/4.SC_exploration/sample_info_L.tsv", header=T, sep="\t", row.names=1) %>% select(x, y, z)
sample.pos <- sample.pos[rownames(sample.expr), ]
sample.expr <- as_tibble(sample.expr)
sample.pos <- as_tibble(sample.pos)

misty.intra <- create_initial_view(sample.expr)
misty.views <- misty.intra %>% add_juxtaview(sample.pos, 10)
misty.views <- misty.views %>% add_paraview(sample.pos, 10)
summary(misty.views)

span = 100
for(i in seq(5000, ncol(sample.expr), span) ){
  start = i
  endl = min(i+span, ncol(sample.expr))
  cat(start, endl-1, "\n")
  out_folder = paste("../res/5.gredient_sdi_cellType/left_hemi/batch", round(i/span)+1, sep="")
  
  misty.views %>% run_misty(num.trees=10000, num.threads = 4,
                            target.subset = colnames(sample.expr)[start:endl],
                            results.folder = out_folder)
  
  #break
}


##################################
# Right Hemi
##################################
sample.expr <- t(read.table("../res/4.SC_exploration/sample_expr_R.tsv", header = TRUE,  sep = "\t",  row.names=1))
colnames(sample.expr) <- gsub("-", "_", colnames(sample.expr), fixed = T)
sample.pos <- read.table("../res/4.SC_exploration/sample_info_R.tsv", header=T, sep="\t", row.names=1) %>% select(x, y, z)
sample.pos <- sample.pos[rownames(sample.expr), ]
sample.expr <- as_tibble(sample.expr)
sample.pos <- as_tibble(sample.pos)

misty.intra <- create_initial_view(sample.expr)
misty.views <- misty.intra %>% add_juxtaview(sample.pos, 10)
misty.views <- misty.views %>% add_paraview(sample.pos, 10)
summary(misty.views)

span = 100
for(i in seq(4400, ncol(sample.expr), span) ){
  start = i
  endl = min(i+span, ncol(sample.expr))
  cat(start, endl-1, "\n")
  out_folder = paste("../res/5.gredient_sdi_cellType/right_hemi/batch", round(i/span)+1, sep="")
  
  misty.views %>% run_misty(num.trees=10000, num.threads = 4,
                            target.subset = colnames(sample.expr)[start:endl],
                            results.folder = out_folder)
  
  #break
}







