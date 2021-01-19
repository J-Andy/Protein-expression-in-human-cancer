#mRNA-protein dataset merging
# This script will
# 1. map Protein groups to ENSG identifiers 
# 2. merge genentech, ccle and anger ExpressionAtlas datasets with proteomics matrix
# 3. create a separate df for each cell line to perform mRNA-protein correlations
# 4. final plots
require(reshape)
require(data.table)


### #first run script to prepare protein data: "generate_protein_quant_files.R"




################## 1. ########################################################
# load the mapping reference file - this is  now a file provided by UniProt but can ba any reference file 
##
xxy.map <- read.table( "HUMAN_9606_idmapping_selected.tab", header = F, sep = "\t", fill = T, stringsAsFactors = FALSE)
colnames(xxy.map)  <- c("UniProtKB.AC","UniProtKB.ID","GeneID..EntrezGene.","RefSeq","GI", "PDB","GO","UniRef100","UniRef90","UniRef50","UniParc","PIR","NCBI.taxon","MIM","UniGene",     "PubMed","EMBL","EMBL.CDS","Ensembl","Ensembl_TRS","Ensembl_PRO","Additional.PubMed")
uniprot.map <- xxy.map[ , c(1, 19)]
rm(xxy.map)

# # using cl data frame - prepared by "generate_protein_quant_files.R" script

# cl is the cell lines matrix
x <- cl
x <- 2^x
xx <- melt(t(x))
colnames(xx) <- c("ID", "Variable", "Value")
#
xx <- cast(xx, ID~Variable, value = "Value", mean)
row.names(xx) <- xx$ID
xx <- xx[ , -1]
xx <- as.data.frame(xx)
xx <- as.data.frame(t(xx))
table(colnames(xx))
table(colnames(cl))

##### perform the protein group to gene mapping
c.l.cl <- colnames(cl)
data.to.map <-  as.data.frame(cl) # xx 
colnames(data.to.map) <- c.l.cl
data.to.map$Majority.protein.IDs <- row.names(data.to.map)
##
data.to.map$ENSG <- "NA"
for(i in 1:nrow(data.to.map)){
  x <- data.frame(strsplit(data.to.map[ i, "Majority.protein.IDs"], split = ";"), stringsAsFactors = FALSE)
  colnames(x) <- "prot"
  # extract canonical UniProt protein ID
  x[,1] <- regmatches(x[,1],regexpr("[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}", x[,1]))
  x <- merge(x, uniprot.map, by.x = "prot", by.y = "UniProtKB.AC" )
  all.genes <- sub(" ","", unlist(strsplit(x[ ,2], ";" )))
  data.to.map[ i, "ENSG"] <- paste( unique(all.genes), collapse = ";")
}
######################### match the protein ids with the uniprot data  ###############################################
xx <- data.to.map
# remove protein groups that have no mapping to an ENSG gene IDs
xx <- xx[ xx$ENSG != "" , ]
# remove all protein groups that map to multiple ENSG gene IDs (this is quiet a lot) - the reasoning here is that we cannot establish for sure which gene is contibuting the signal to protein; all genes contribute equally or onbe gene is a majority? 
xx <- xx[ grep(";", xx$ENSG, invert = TRUE) , ]
# for genes that map to multiple proteins, in order to detemine the amount of protein that gene is producing - sum the protein quantification values
xx.Majority.protein.IDs <- aggregate(xx$Majority.protein.IDs, list(ESNG = xx$ENSG ), function(x) paste0( (x) )  )
xx <- aggregate(xx[ , 1:(ncol(xx)-2)], list(ENSG = xx$ENSG), sum, na.rm =TRUE)
#
xx <- cbind(xx.Majority.protein.IDs, xx)

###########################################################################
##################################################################################
genentech <- read.table("E-MTAB-2706-query-results.fpkms.tsv", header = T, sep = "\t", fill = T, stringsAsFactors = FALSE)# , check.names = FALSE

ccle <- read.table("E-MTAB-2770-query-results.fpkms.tsv", header = T, sep = "\t", fill = T, stringsAsFactors = FALSE, comment.char = "#")

sanger <- read.table("E-MTAB-3983-query-results.fpkms.tsv", header = T, sep = "\t", fill = T, stringsAsFactors = FALSE, comment.char = "#")
#######
MyMerge <- function(x, y){
  df <- merge(x, y, by= "Gene.ID", all.x= TRUE, all.y= TRUE)
  return(df)
}

cgs <- Reduce(MyMerge, list(ccle, genentech, sanger))
cgs.cols <- colnames(cgs)
row.names(cgs) <- cgs$Gene.ID
ind <- cgs$Gene.ID %in% xx$ESNG
rna <- cgs[ ind , ]
rna.cells <-  sapply( strsplit( colnames(rna), "\\.\\."), "[",1 )
rna.cells <- toupper(rna.cells)
rna.cells <- sub("NCI\\.", "", rna.cells)
rna.cells <- gsub("\\.", "", rna.cells)
rna.cells <- gsub("^X", "", rna.cells)
rna.cells
cgs <- cgs[ ind , -c(1, grep("GENENAM.*" , rna.cells)) ]
colnames(cgs) <- rna.cells[-c(1, grep("GENENAM.*" , rna.cells))]
cgs <- melt(t(cgs))
#
colnames(cgs) <- c("ID", "Variable", "Value")
cgs$ID <- as.character(cgs$ID)
cgs$Variable <- as.character(cgs$Variable)
cgs <- as.data.table(cgs)
cgs <- cgs[, mean(Value, na.rm = T ), by = c("ID", "Variable")]
cgs <- dcast(cgs, ID~Variable, value.var = "V1")
row.names(cgs) <- cgs$ID
cgs.cols <- cgs$ID
cgs <- as.data.frame(cgs[ , -1])
cgs <- as.data.frame(t(cgs))
colnames(cgs) <- cgs.cols
#########################################################################################
#########################################################################################
prot <- xx
#########################################################################################
#########################################################################################
row.names(prot) <- prot$ESNG
prot <- prot[ , -c(1:3)]
colnames(prot) <- c.l.cl
##############################
l.prot <- list()
l.rna <- list()
cors <- vector(mode = "numeric", length = length(colnames(prot)))
names(cors) <-   colnames(prot)
cors.on.how.many.genes <- c()
for(i in 1:length(colnames(prot))){
  print(names(prot)[i])
  rna.1 <- cgs[, colnames(cgs) %in% colnames(prot)[i], drop = FALSE]
  prot.1 <- prot[ , i, drop = FALSE]
  rna.1$prot <- row.names(rna.1)
  prot.1$prot <- row.names(prot.1)
  rna.prot <- merge(rna.1, prot.1, by="prot", all= F)
  rna.prot <- rna.prot[ complete.cases(rna.prot), ]
  print( nrow(rna.prot) )
  if(ncol(rna.prot) == 2){
    cors[i] <- NA
    cors.on.how.many.genes[i] <- NA
  } else {
    cor.1 <- cor( rna.prot[ ,2], rna.prot[,3],use="pairwise.complete.obs", method = "spearman")
    cors[i] <- cor.1
    cors.on.how.many.genes[i] <- nrow(rna.prot)
  }
  
}
names(cors)
cors <- data.frame(cors)
cors$c.l <- c.l.cl
cors$assay.names <- assay.names
##

mean(cors.on.how.many.genes, na.rm = T )
hist(cors[,1],  xlim = c(0.4,0.75), col = "#fa9fb5")#breaks = 147,
abline(v = median(cors[,1], na.rm = T), lwd = 2, col = "black", lty = 2.5)
summary(cors$cors)
table(round(cors, 2))
median(cors[,1], na.rm = T)
#


############
############
#
meta <-  cell.metadata <- read.table( "Supplementary-Table-1-samples-linegae-metadata_FINAL.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE,  na.strings = c(NA, "NA", "NaN"), strip.white = T)
cell.metadata <- read.table( "cell-lines-metadata.complete.cosmic.txt", sep ="\t", header = TRUE, stringsAsFactors = FALSE)
cell.metadata <- cell.metadata[ !duplicated(cell.metadata$my.cells), ]

lineage <- as.character(cell.metadata$Lineage[match((cors$c.l), cell.metadata$my.cells)])
lineage.stats <- table(lineage)
lineage <- gsub("large_intestine" , "colorectal", lineage)
lineage <- paste0(toupper(substr(lineage , 1, 1)), substr(lineage , 2, nchar(lineage )))
cors <- cbind(cors, lineage)

table(cors$lineage)
cors <- cors[ complete.cases(cors), ]
cors <- cors[cors$lineage %in% names(which(table(cors$lineage) > 2)), ]
cors <- droplevels(cors)
cors$batch <- factor(cors$batch)
cors$lineage <- factor(cors$lineage)
boxplot(cors$cors ~ cors$lineage, las = 2)
aggregate(cors ~ as.factor(lineage), data = cors, median) 
bymedian.cors <- with(cors, reorder(lineage, -cors, median))
op <- par(mar=c(7,4,4,1))
boxplot(cors ~ bymedian.cors, data = cors,
        ylab = "Correlation",
        main = "", varwidth = F,
        col = "#fde0dd", las = 2,outline=FALSE, ylim = c(min(cors$cors, na.rm = T), max(cors$cors, na.rm = T)),frame.plot = FALSE ) #xlab = "Cell line lineage",
stripchart(cors ~ bymedian.cors, data = cors,  vertical=T, method="jitter", pch=19,add=TRUE, col = "#737373", cex = 1) 
rm(op)#

cors$batch <- sapply(strsplit(cors$assay.names, "_"), "[" , 3)
bymedian.cors <- with(cors, reorder(batch, -cors, median))
boxplot(cors ~ bymedian.cors, data = cors,
        ylab = "Correlation",
        main = "", varwidth = F,
        col = "#fde0dd", las = 2,outline=FALSE, ylim = c(min(cors$cors, na.rm = T), max(cors$cors, na.rm = T)),frame.plot = FALSE ) #xlab = "Cell line lineage",
stripchart(cors ~ bymedian.cors, data = cors,  vertical=T, method="jitter", pch=19,add=TRUE, col = "#737373", cex = 1) 


kruskal.test(cors$cors ~ cors$lineage)
kruskal.test(cors$cors ~ as.factor(cors$batch))


aov1 <- aov(cors$cors ~ cors$lineage * as.factor(cors$batch))
summary(aov1)


glm1 <- glm(cors$cors ~ cors$lineage * cors$batch)
summary(glm1)
#

two.way <- aov(cors$cors ~ cors$lineage + cors$batch)
summary(two.way)

plot(two.way)

tukey.two.way<-TukeyHSD(two.way)
x <- tukey.two.way$`cors$lineage`
x.x <-tukey.two.way$`cors$batch` 
plot(tukey.two.way, las = 1)
boxplot(tukey.two.way)
##
tes.result <- pairwise.wilcox.test(cors$cors , cors$lineage,
                                   p.adjust.method = "BH")
tes.result <- pairwise.wilcox.test(cors$cors , as.factor(cors$batch),
                                   p.adjust.method = "BH")
tes.result
sum(tes.result$p.value < 0.05, na.rm = T)
tes.result <- tes.result[["p.value"]][tes.result$p.value < 0.05]

require(plyr)
library(corrplot)
plot.table(as.matrix(tes.result$p.value))
# , smain='Cohort(users)', highlight = TRUE, colorbar = TRUE)
M <- as.matrix(tes.result$p.value)
corrplot(M, is.cor = FALSE, typ = "lower", tl.col = "black", method = "number",  col = "red",p.mat = tes.result$p.value, sig.level = c(.05), order = "original", insig = "blank", cl.pos = "n",  number.digits = 3 )

