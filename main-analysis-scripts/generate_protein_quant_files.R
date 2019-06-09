# load required packages
# to install bioconductor packages try http:// if https:// URLs are not supported
# for example, to get the pcaMethods package run:
# source("https://bioconductor.org/biocLite.R")
# biocLite("pcaMethods")
require(data.table);require(pcaMethods);require(limma);require(reshape2); require(lattice);require(plyr)

# load utility functions
lapply(list.files(path = "./utility-functions/", pattern = "[.]R$", recursive = TRUE), function(x) source( paste("./utility-functions/" ,x, sep="")  ))

##################################################################################
##################################################################################
# process cell line data
##################################################################################
##################################################################################
#
# the initial proteinGroups.txt file containing cell line data is available in PXD013455 inside the txt-celllines.zip folder. You will need to download that file in order to follow the analysis
tmp  <-  fread( "E:/cancer-analysis/data/proteomics/cellLines/txt-7studies/proteinGroups.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE,  na.strings = c(NA, "NA", "NaN"), strip.white = T)
# perform standard clean-up steps 
tmp <- data.frame(tmp)
tmp <- tmp[ tmp[, "Reverse"] != "+", ]
tmp <- tmp[ tmp[, "Potential.contaminant"] != "+", ]  
# subset the data to iBAQ only
tmp.2 <- tmp[ ,c(2, grep("iBAQ.", colnames(tmp))) ]
tmp.2 <- tmp.2[ , -1]
# fot normalise
tmp.2 <- fot.normalise(tmp.2)
tmp <- cbind(tmp, tmp.2)

# one can keep only protein groups that were identified by unique peptides. We keep all proteions in this case. There will be further filtering steps later on that are effectiviely removing any no-unique groups. 
tmp <- tmp[  tmp$Unique.peptides >= 0 , ]
tmp <- tmp[ ,c(2, grep("ppb.", colnames(tmp))) ]
tmp[tmp == 0] <- NA
tmp[ tmp == "NaN"] <- NA
# create cell lines matrix
cl <- tmp
rm(tmp.2, tmp)
# save the intermediate file, not included in this repo due to size
# write.table( cl, "./data/proteinGroups_cellLines_ppbNorm.txt", sep = "\t", row.names = F)
#
#
# some further data wrangling is required here to remove assays that we are not interested in (such as tumour samples, or non-cancer cell lines)
assays.to.keep <- c("Majority.protein.IDs", "iBAQ.Tumour.")
tmp.tm <- cl[ , grep(paste(assays.to.keep, collapse = "|"), colnames(cl), invert = FALSE), drop = FALSE]
assays.to.rm <- c( "MCF10A", "HEK293","iBAQ.FallopianTubeEpithelialCells", "iBAQ.HEYA8", "iBAQ.Tumour.", "_kinome", "_LysC", "Chymo_", "GluC_")
cl <- cl[ , grep(paste(assays.to.rm, collapse = "|"), colnames(cl), invert = TRUE), drop = FALSE]
row.names(cl) <- cl$Majority.protein.IDs
cl <- cl[, -1]
#######################################
# calculate the % of missing values in the current matrix
sum(is.na(cl)) / (sum(is.na(cl)) + sum(!is.na(cl)))
# keep the proteins (rows) with a minimum of 50% valid values 
cl <- cl[-which(rowMeans(is.na(cl)) > 0.5), ]
# % of missing values after filtering
sum(is.na(cl)) / (sum(is.na(cl)) + sum(!is.na(cl)))
# IMPORTANT: here we log2 transform the intensity values
cl <- log2(cl)
# save the intermediate file
# write.table( cbind(Majority.protein.IDs = row.names(cl), cl), "./data/proteinGroups_cellLines_ppbNorm_min50validValues.txt", sep = "\t", row.names = F)
######################################
##
#
# impute missing data for the cell line matrix
#
##
result <- pca(cl, method="svdImpute", nPcs=3, center = TRUE, scale="uv")
## Get the estimated complete observations
cObs <- completeObs(result)
## Now plot the scores
plotPcs(result, type = "scores")
plot(sDev(result))
# slplot(result)
cl <- as.data.frame(cObs)
#
# write.table( cbind(Majority.protein.IDs = row.names(cl), cl), "./data/proteinGroups_cellLines_ppbNorm_min50validValues_svdImputeMissing.txt", sep = "\t", row.names = F)
##
#
# remove major batch effects, NOTE when it comes to differnetial analysis with limma it is recommended to include batch as covariate.
#
##
batch <- sapply(strsplit(colnames(cl), "_"), function(x) x[3] )
batch
dat.rm.batch <- removeBatchEffect(cl, batch )
cl <- dat.rm.batch
# write.table( cbind(Majority.protein.IDs = row.names(cl), cl), "./data/proteinGroups_cellLines_ppbNorm_min50validValues_svdImputeMissing_removeBatchEffect.txt", sep = "\t", row.names = F)



##################################################################################
##################################################################################
# process tumours data
##################################################################################
##################################################################################
# the proteinGroups.txt file containing tumor data is available in PXD013455 inside txt-tumours.zip folder
# you will need to download that file in order to follow this analysis
tmp  <-  fread( "E:/cancer-analysis/data/proteomics/tumours/txt/proteinGroups.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE,  na.strings = c(NA, "NA", "NaN"), strip.white = T) # comment.char = "#",
# standard clean-up steps 
tmp <- data.frame(tmp)
tmp <- tmp[ tmp[, "Reverse"] != "+", ]
tmp <- tmp[ tmp[, "Potential.contaminant"] != "+", ]  
# subset the data to iBAQ only
tmp.2 <- tmp[ ,c(2, grep("iBAQ.", colnames(tmp))) ]
tmp.2 <- tmp.2[ , -1]
# fot normalise
tmp.2 <- fot.normalise(tmp.2)
tmp <- cbind(tmp, tmp.2)
# one can keep only protein groups that were identified by unique peptides. We keep all proteions in this case. There will be further filtering steps later on that are effectiviely removing no-unique groups. 
tmp <- tmp[  tmp$Unique.peptides >= 0 , ]
tmp <- tmp[ ,c(2, grep("ppb.", colnames(tmp))) ]
tmp[tmp == 0] <- NA
tmp[ tmp == "NaN"] <- NA
# create tumours matrix
tm <- tmp
rm(tmp.2, tmp)
# save the intermediate file, not included in the repo due to size
# write.table( cl, "./data/proteinGroups_tumours_ppbNorm.txt", sep = "\t", row.names = F)
######################################
# some further data wrangling is required here to remove assays that we are not interested in, and merge with remaining tumour data 
assays.to.keep <- c("Majority.protein.IDs", "Zhang_", "L.Tumour")
tm <- tm[ , grep(paste(assays.to.keep, collapse = "|"), colnames(tm), invert = FALSE), drop = FALSE]
assays.to.rm <- c(  "_LysC")
tmp.tm <- tmp.tm[ , grep(paste(assays.to.rm, collapse = "|"), colnames(tmp.tm), invert = TRUE), drop = FALSE]
#######################
tm$Leading.razor.ID <- sub(";.*", "", tm$Majority.protein.IDs)
tmp.tm$Leading.razor.ID <- sub(";.*", "", tmp.tm$Majority.protein.IDs)
##
tm <- merge(tm, tmp.tm, by ="Leading.razor.ID", all = TRUE) 
tm <- tm[ ,colnames(tm) != "Majority.protein.IDs.y"]
tm <- tm[ ,colnames(tm) != "Majority.protein.IDs.x"]
row.names(tm) <- tm$Leading.razor.ID
tm <- tm[ , -1]
#######################################
# % of missing values in the current matrix
sum(is.na(tm)) / (sum(is.na(tm)) + sum(!is.na(tm)))
# keep the proteins (rows) with a minimum of 50% valid values
tm <- tm[-which(rowMeans(is.na(tm)) > 0.5), ]
# % of missing values after filtering
sum(is.na(tm)) / (sum(is.na(tm)) + sum(!is.na(tm)))
# IMPORTANT: here we log2 transform the intensity values
tm <- log2(tm)
# save the intermediate file
# write.table( cbind(Leading.razor.ID = row.names(tm), tm), "./data/proteinGroups_tumours_ppbNorm_min50validValues.txt", sep = "\t", row.names = F)
######################################
#
result <- pca(tm, method="svdImpute", nPcs=3, center = TRUE, scale="uv")
## Get the estimated complete observations
cObs <- completeObs(result)
## Now plot the scores
plotPcs(result, type = "scores")
plot(sDev(result))
# slplot(result)
tm <- as.data.frame(cObs)
# write.table( cbind(Leading.razor.ID = row.names(tm), tm), "./data/proteinGroups_tumours_ppbNorm_min50validValues_svdImputeMissing.txt", sep = "\t", row.names = F)
##
#
# remove major batch effects, NOTE when it comes to differnetial analysis with limma it is recommended to include batch as covariate.
#
##
batch <- sapply(strsplit(colnames(tm), "_"), function(x) x[3] )
batch
dat.rm.batch <- removeBatchEffect(tm, batch )
tm <- dat.rm.batch
# write.table( cbind(Leading.razor.ID = row.names(tm), tm), "./data/proteinGroups_tumours_ppbNorm_min50validValues_svdImputeMissing_removeBatchEffect.txt", sep = "\t", row.names = F)
# clean up env
rm(list = ls())
##############################################################################
##############################################################################
##############################################################################
##############################################################################