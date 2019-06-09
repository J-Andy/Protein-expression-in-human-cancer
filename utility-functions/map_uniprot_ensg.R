# Utility function to perform UniProt protein ID to ensmbl gene IDs mapping
# This script will
# 1. map Protein groups to ENSG identifiers 
################## 1. ########################################################

# load the mapping reference file - this is a file provided by UniProt but can ba any reference file 

# # download id mapping file for Human from UniProt
## no need to run the download part if you already have the mapping file
# uniprot_url <- "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping_selected.tab.gz"
# temp <- tempfile()
# download.file(uniprot_url, temp, method = "libcurl", mode = "wb")
# uniprot.map <- read.table(gzfile(temp), header = F, sep = "\t", fill = T, stringsAsFactors = FALSE)
# unlink(temp)
# colnames(uniprot.map) <- c("UniProtKB.AC","UniProtKB.ID","GeneID..EntrezGene.","RefSeq","GI", "PDB","GO","UniRef100","UniRef90","UniRef50","UniParc","PIR","NCBI.taxon","MIM","UniGene",     "PubMed","EMBL","EMBL.CDS","Ensembl","Ensembl_TRS","Ensembl_PRO","Additional.PubMed")
# uniprot.map <- uniprot.map[ , c(1, 19)]


# if you have the full file downloaded from uniprot
# xxy.map <- read.table( "E:/cancer-analysis/data/mRNA/expression-atlas-cellLines/HUMAN_9606_idmapping_selected.tab/HUMAN_9606_idmapping_selected.tab", header = F, sep = "\t", fill = T, stringsAsFactors = FALSE)
# colnames(xxy.map)  <- c("UniProtKB.AC","UniProtKB.ID","GeneID..EntrezGene.","RefSeq","GI", "PDB","GO","UniRef100","UniRef90","UniRef50","UniParc","PIR","NCBI.taxon","MIM","UniGene",     "PubMed","EMBL","EMBL.CDS","Ensembl","Ensembl_TRS","Ensembl_PRO","Additional.PubMed")
# uniprot.map <- xxy.map[ , c(1, 19)]
# rm(xxy.map)

map_uniprot_ensg <- function(data.to.map = TRUE, xxy.map = TRUE){

# in this case we use the mapping file that is already provided
if(xxy.map == TRUE){
xxy.map <- read.table( "HUMAN_9606_idmapping_selected.txt", header = F, sep = "\t", fill = T, stringsAsFactors = FALSE)
colnames(xxy.map)  <- c("UniProtKB.AC","Ensembl")
# a note about Accession and ID numbers in Uniprot: https://www.uniprot.org/help/difference_accession_entryname
# What is the difference between an accession number (AC) and the entry name?
#   
#   An accession number (AC) is assigned to each sequence upon inclusion into UniProtKB. Accession numbers are stable from release to release. If several UniProtKB entries are merged into one, for reasons of minimizing redundancy, the accession numbers of all relevant entries are kept. Each entry has one primary AC and optional secondary ACs.
# 
# The 'Entry name' (formerly ID) is a unique identifier, often containing biologically relevant information. It is sometimes necessary, for reasons of consistency, to change IDs (for example to ensure that related entries have similar names). Another common cause for changing an ID is when an entry is promoted from UniProtKB's TrEMBL section (with computationally-annotated records) to the Swiss-Prot section (with fully curated records). However, an accession number (AC) is always conserved, and therefore allows unambiguous citation of UniProt entries. If a UniProtKB entry contains more than one accession number, the first one (or primary accession number) should be cited. 
uniprot.map <- xxy.map
rm(xxy.map)
}


if(data.to.map == TRUE){
data.to.map <- read.table("data.to.map.example-withInt.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE,  na.strings = c(NA, "NA", "NaN"), strip.white = T)
}

##### perform the protein group to gene mapping


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
# remove all protein groups that map to multiple ENSG gene IDs (this removes a lot of proteins) - the reasoning here is that we cannot establish for sure which gene is contributing the signal to the protein abundance; all genes contribute equally or one gene is a majority? 
xx <- xx[ grep(";", xx$ENSG, invert = TRUE) , ]
# for genes that map to multiple proteins, in order to determine the amount of protein that gene is producing - sum the protein quantification values
xx.Majority.protein.IDs <- aggregate(xx$Majority.protein.IDs, list(ESNG = xx$ENSG ), function(x) paste0( (x) )  )
### select which columns to aggregate
# in the example file
xx <- aggregate(xx[ , 3:3], list(ENSG = xx$ENSG), sum, na.rm =TRUE)
# other cases
# xx <- aggregate(xx[ , 1:(ncol(xx)-2)], list(ENSG = xx$ENSG), sum, na.rm =TRUE)

#

xx <- cbind((as.character(xx.Majority.protein.IDs$x)), xx)


# check if it worked in the example file - if TURE it did
xx[ grep("ENSG00000162664", xx$ESNG), 4] == 20

}