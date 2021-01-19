require(plyr)
require(purrr)
require(reshape2)
require(pheatmap)
library(dplyr)

cl <- read.table( "../data/proteinGroups_cellLines_ppbNorm_min50validValues_svdImputeMissing_removeBatchEffect.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE,  na.strings = c(NA, "NA", "NaN"), strip.white = T) # comment.char = "#",

tm <- read.table( "../data/proteinGroups_tumours_ppbNorm_min50validValues_svdImputeMissing_removeBatchEffect.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE,  na.strings = c(NA, "NA", "NaN"), strip.white = T) # comment.char = "#",

row.names(cl) <- sub(";.*", "", cl$Majority.protein.IDs)
row.names(tm) <- tm$Leading.razor.ID
cl <- cl[ , -1]
tm <- tm[ , -1]
####
cl.tm <- merge(cl, tm, by= "row.names", all = T )
row.names(cl.tm) <- cl.tm$Row.names
cl.tm <- cl.tm[ , -1]

metadata <- read.table( "Supplementary-Table-1-samples-linegae-metadata_FINAL.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE,  na.strings = c(NA, "NA", "NaN"), strip.white = T) # comment.char = "#",
#
cor_mat <- cor(as.matrix(cl.tm), use = "pairwise.complete.obs", method = "pearson")
##
annotation <- data.frame(
  Lineage = metadata$Lineage.Assigned[match(colnames(cl.tm), metadata$Assay.Name)], 
  Sample.type =metadata$Sample.type[match(colnames(cl.tm), metadata$Assay.Name)] , 
  Study = metadata$Study..Batch.[match(colnames(cl.tm), metadata$Assay.Name)])# , 
row.names(annotation) <- colnames(cl.tm)

annotation_colors = list( 
  Sample.type = c(CellLine ="#e0e0e0", Tumour = "#1a1a1a"),
  Lineage = c(
    Skin = "#a6cee3",
    Prostate = "black", # "#33a02c",
    Ovary = "#b15928",
    Lymph.node = "#ffff99",
    Lung = "#a6cee3",  # "#b3de69",
    Liver = "#a6cee3",# "#fccde5",
    Colorectal = "#6a3d9a",
    Kidney = "#a6cee3" , # "#35978f",
    Cervix = "#a6cee3" ,# #fb9a99",
    Breast =  "#c51b8a",
    Breast.Metastatic = "#c51b8a",
    Brain = "#a6cee3", 
    Bone = "#a6cee3" , 
    Blood = "#a6cee3" )) 

##########################
##
pheatmap(cor_mat, scale = "none", show_colnames = F, show_rownames = F, cluster_cols=T, cluster_rows=T, annotation_col = annotation, annotation_row = annotation, annotation_colors = annotation_colors)

pheatmap(cor_mat, scale = "none", cluster_rows = T, cluster_cols = T, clustering_distance_rows = 'euclidean', clustering_distance_cols = 'euclidean', clustering_method = 'ward.D2',  legend = T, drop_levels = F, na_col = "grey",  show_colnames = F,show_rownames = F, border_color=NA, cex = 1, heatmap_legend_side = "top",  annotation_colors = annotation_colors, annotation_row = annotation,  annotation_col = annotation, annotation_legend_side = "top") # ,


pheatmap(cor_mat, scale = "none", cluster_rows = T, cluster_cols = T, clustering_distance_rows = 'euclidean', clustering_distance_cols = 'euclidean', clustering_method = 'ward.D2',  legend = T, drop_levels = F, na_col = "grey",  show_colnames = F,show_rownames = F, border_color=NA, cex = 1, heatmap_legend_side = "top",  annotation_colors = annotation_colors,annotation_col = annotation, annotation_legend_side = "top") # ,

###
cor_mat[lower.tri(cor_mat,diag=TRUE)]=NA # put NA
cor_mat<-as.data.frame(as.table(cor_mat)) # as a dataframe
cor_mat<-na.omit(cor_mat) # remove NA
cor_mat<-cor_mat[with(cor_mat, order(-Freq)), ] # order by correlation

summary(cor_mat)
median(cor_mat$Freq)
min(cor_mat$Freq, na.rm = T)
max(cor_mat$Freq,na.rm = T)
#############################################
#
rmat.test <- cor_mat
annotation$name <- row.names(annotation)
rmat.test <- merge(rmat.test, annotation, by.x = "Var1", by.y = "name")
rmat.test <- merge(rmat.test, annotation, by.x = "Var2", by.y = "name")
rmat.test.backup <- rmat.test
rmat.test$lineage.type.x <- paste(rmat.test$Lineage.x, rmat.test$Sample.type.x, sep = "_")
rmat.test$lineage.type.y <- paste(rmat.test$Lineage.y, rmat.test$Sample.type.y, sep = "_")
x <- rmat.test %>% group_by(Lineage.x, Sample.type.x, Lineage.y, Sample.type.y) %>% summarize(median.cor = median(Freq))
x.x <- rmat.test %>% group_by(lineage.type.x, lineage.type.y) %>% summarize(median.cor = median(Freq))
####

library(lattice)
bwplot(Freq ~ Lineage.x  + Lineage.y  | Sample.type.x + Sample.type.y ,         data=rmat.test,
       main="Cors")
#
bwplot(Freq ~ Lineage.x  : Lineage.y  | Sample.type.x + Sample.type.y , 
       data=rmat.test,layout=c(1,1),
       main="Cors")

panel.mn <- function(x,y,box.width=.5,horiz=FALSE,...){
  panel.bwplot(x,y,box.width,horiz=horiz,...)
  y <- tapply(y,x,median,na.rm=TRUE)
  x <- seq_along(y)
  panel.segments(x0=x-box.width/2,x1=x+box.width/2,
                 y0 = y,y1=y,...)
}

bwplot(Freq ~  Sample.type.x : Sample.type.y ,       data=rmat.test,
       main="", panel=panel.mn, ylab = "Spearman correlation coefficient")
##
tapply(rmat.test$Freq, rmat.test$Sample.type.x : rmat.test$Sample.type.y,  summary)
#