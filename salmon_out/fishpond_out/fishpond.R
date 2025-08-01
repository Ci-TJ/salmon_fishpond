#Error in `collect()`: ! Failed to collect lazy table. Caused by error in `db_collect()`: ! Arguments in `...` must be used. ? Problematic argument: ? ..1 = Inf ? Did you misspell an argument name? Run `rlang::last_trace()` to see where the error occurred.
#https://stackoverflow.com/questions/77370659/error-failed-to-collect-lazy-table-caused-by-error-in-db-collect-using
#install.packages("devtools")
#devtools::install_version("dbplyr", version = "2.3.4")
#env NGS


cat("====================================================",date(),"====================================================","\n")
cat("====================================================",date(),"====================================================","\n")
library(tximeta)
library(fishpond)
library(data.table)
library(tidytable)
suppressPackageStartupMessages(library(SummarizedExperiment))
library(readxl)

########################################
##00 settings
########################################
species <- "mouse"
seed <- 1
keep <- TRUE

##01 Loading data
WD <- "/home/user_li/linqin_tmp/Analysis/RNA-seq/tac112/salmon_out/"
project <- "tac112"

colN <- dir(file.path(WD, "fishpond_out"))[grep("xlsx$",dir(file.path(WD, "fishpond_out")))]
coldata0 <- read_xlsx(file.path(WD, "fishpond_out",colN)) %>% filter(Condition != "other") #only contain the wildtype and tac samples
#Add in 20240623
coldata <- coldata0 
coldata$names<-coldata$Run ## Note: fishpond need colnames "names"
coldata$condition <- coldata$Condition
head(coldata)
coldata$files <- file.path(WD, "quant_out",coldata$names, "quant.sf")
all(file.exists(coldata$files))
if (all(file.exists(coldata$files))){
	coldata <- coldata
	} else {
	coldata <- data.frame(coldata)[file.exists(coldata$files),]
	}


head(coldata)
all(file.exists(coldata$files))
rownames(coldata) <- 1:nrow(coldata)

########################################
##02 Tximeta
########################################
indexdir_hsa='/home/user_li/linqin_tmp/Genome/hsa/release107/index_human_r107_salmon_k31'
gtfpath_hsa='/home/user_li/linqin_tmp/Genome/hsa/release107/Homo_sapiens.GRCh38.107.gtf'
faPath_hsa = "/home/user_li/linqin_tmp/Genome/hsa/release107/Homo_sapiens.GRCh38.all.rna.fa"

indexdir_mus='/home/user_li/linqin_tmp/Genome/mus/release107/index_mouse_r107_salmon_k31'
gtfpath_mus='/home/user_li/linqin_tmp/Genome/mus/release107/Mus_musculus.GRCm39.107.gtf'
faPath_mus = "/home/user_li/linqin_tmp/Genome/mus/release107/Mus_musculus.GRCm39.all.rna.fa"

if (species == "human"){
	index <- indexdir_hsa
	gtf <- gtfpath_hsa
	fasta <- faPath_hsa
	org = "Homo Sapiens"
	genome = "GRCh38"

	} else {
	index <- indexdir_mus
	gtf <- gtfpath_mus
	fasta <- faPath_mus
	org = "Mus musculus"
        genome = "GRCm39"
	}


makeLinkedTxome(
    indexDir = index,
    source = "Ensembl",
    organism = org,
    release = "107",
    genome = genome,
    fasta = fasta,
    gtf = gtf,
    write = FALSE
)

###############################################
##3.1 4w TAC vs healthy
###############################################
#condition for Sham_*&TAC_*
#filter the samples
cg <- "Sham-w4"
eg <- "TAC-w4"
outName <- "TACvsSham_4w_gse66630"

coldata1 <- coldata %>% drop_na(condition) %>% filter(`Series ID` == "GSE66630") %>% filter(condition %in% c(cg, eg)) #Note coldata1 !!!!!!!!!!!!!!
se <- tximeta(coldata1) # delete skipMeta=TRUE
assayNames(se)
head(rownames(se))
se
class(se)
cat("tximeta OK!", "\n")
##03 Differential expression analysis at gene level
gse <- summarizeToGene(se)

#note condition for TAC_w4 vs Sham_w4
y <- gse[,gse$condition %in% c(cg, eg)]
y$condition <- factor(y$condition, c(cg, eg)) #the ordering of 'sham' and 'tac' will determine the log2FC
y
y <- scaleInfReps(y)

## Don't filter the genes by count and number of samples, preserve all the genes for annotation
#y <- labelKeep(y)
#y <- y[mcols(y)$keep,]
if (keep){
	y <- labelKeep(y)
	y <- y[mcols(y)$keep,]
	}

cat("scaleInfReps-labelKeep OK!", "\n")
set.seed(seed)
cat("Start swish!")
y <- swish(y, x="condition", nperms=100) #x="condition"
y
head(mcols(y))
cat("swish OK!", "\n")



##04 Output data
Gstat <-data.frame(mcols(y))
annot <- Gstat %>% select(gene_id, gene_name) #select.() to select()

TPM<-data.frame((assays(gse)$abundance));TPM$gene_id<-rownames(TPM); TPM <- annot %>% inner_join(TPM, by = "gene_id")
count<-data.frame((assays(gse)$count));count$gene_id<-rownames(count); count <- annot %>% inner_join(count, by = "gene_id")

fwrite(count, paste0("Gcount_", outName, ".txt"), sep = "\t")
fwrite(TPM, paste0("Gtpm_", outName, ".txt"), sep = "\t")
fwrite(Gstat, paste0("Gstat_", outName, ".txt"), sep = "\t")

rm(se, gse, y, Gstat, annot, TPM, count) #free up memory
cat(outName, "done!", "\n")


###############################################
##3.2 5w TAC vs healthy
###############################################
#condition for Sham_*&TAC_*
#filter the samples
cg <- "Sham-w5"
eg <- "TAC-w5"
outName <- "TACvsSham_5w_gse112055"

coldata1 <- coldata %>% drop_na(condition) %>% filter(`Series ID` == "GSE112055") %>% filter(condition %in% c(cg, eg)) #Note coldata1 !!!!!!!!!!!!!!
se <- tximeta(coldata1) # delete skipMeta=TRUE
assayNames(se)
head(rownames(se))
se
class(se)
cat("tximeta OK!", "\n")
##03 Differential expression analysis at gene level
gse <- summarizeToGene(se)

#note condition for TAC_w5 vs Sham_w5
y <- gse[,gse$condition %in% c(cg, eg)]
y$condition <- factor(y$condition, c(cg, eg)) #the ordering of 'sham' and 'tac' will determine the log2FC
y
y <- scaleInfReps(y)

## Don't filter the genes by count and number of samples, preserve all the genes for annotation
#y <- labelKeep(y)
#y <- y[mcols(y)$keep,]
if (keep){
	y <- labelKeep(y)
	y <- y[mcols(y)$keep,]
	}

cat("scaleInfReps-labelKeep OK!", "\n")
set.seed(seed)
cat("Start swish!")
y <- swish(y, x="condition", nperms=100) #x="condition"
y
head(mcols(y))
cat("swish OK!", "\n")



##04 Output data
Gstat <-data.frame(mcols(y))
annot <- Gstat %>% select(gene_id, gene_name) #select.() to select()

TPM<-data.frame((assays(gse)$abundance));TPM$gene_id<-rownames(TPM); TPM <- annot %>% inner_join(TPM, by = "gene_id")
count<-data.frame((assays(gse)$count));count$gene_id<-rownames(count); count <- annot %>% inner_join(count, by = "gene_id")

fwrite(count, paste0("Gcount_", outName, ".txt"), sep = "\t")
fwrite(TPM, paste0("Gtpm_", outName, ".txt"), sep = "\t")
fwrite(Gstat, paste0("Gstat_", outName, ".txt"), sep = "\t")

rm(se, gse, y, Gstat, annot, TPM, count) #free up memory
cat(outName, "done!", "\n")

###############################################
##3.3 8w TAC vs healthy
###############################################
#condition for Sham_*&TAC_*
#filter the samples
cg <- "Sham-w8"
eg <- "TAC-w8"
outName <- "TACvsSham_8w_gse66630"

coldata1 <- coldata %>% drop_na(condition) %>% filter(`Series ID` == "GSE66630") %>% filter(condition %in% c(cg, eg)) #Note coldata1 !!!!!!!!!!!!!!
se <- tximeta(coldata1) # delete skipMeta=TRUE
assayNames(se)
head(rownames(se))
se
class(se)
cat("tximeta OK!", "\n")
##03 Differential expression analysis at gene level
gse <- summarizeToGene(se)

#note condition for TAC_w4 vs Sham_w4
y <- gse[,gse$condition %in% c(cg, eg)]
y$condition <- factor(y$condition, c(cg, eg)) #the ordering of 'sham' and 'tac' will determine the log2FC
y
y <- scaleInfReps(y)

## Don't filter the genes by count and number of samples, preserve all the genes for annotation
#y <- labelKeep(y)
#y <- y[mcols(y)$keep,]
if (keep){
	y <- labelKeep(y)
	y <- y[mcols(y)$keep,]
	}

cat("scaleInfReps-labelKeep OK!", "\n")
set.seed(seed)
cat("Start swish!")
y <- swish(y, x="condition", nperms=100) #x="condition"
y
head(mcols(y))
cat("swish OK!", "\n")



##04 Output data
Gstat <-data.frame(mcols(y))
annot <- Gstat %>% select(gene_id, gene_name) #select.() to select()

TPM<-data.frame((assays(gse)$abundance));TPM$gene_id<-rownames(TPM); TPM <- annot %>% inner_join(TPM, by = "gene_id")
count<-data.frame((assays(gse)$count));count$gene_id<-rownames(count); count <- annot %>% inner_join(count, by = "gene_id")

fwrite(count, paste0("Gcount_", outName, ".txt"), sep = "\t")
fwrite(TPM, paste0("Gtpm_", outName, ".txt"), sep = "\t")
fwrite(Gstat, paste0("Gstat_", outName, ".txt"), sep = "\t")

rm(se, gse, y, Gstat, annot, TPM, count) #free up memory
cat(outName, "done!", "\n")

###############################################
##3.4 output all
###############################################
#filter the samples
cg <- "Sham"
eg <- "TAC"
outName <- "all_tac112"

#not filter
coldata1 <- coldata ################ %>% drop_na(etiology) #Note coldata1 !!!!!!!!!!!!!!
se <- tximeta(coldata1) # delete skipMeta=TRUE
assayNames(se)
head(rownames(se))
se
class(se)
cat("tximeta OK!", "\n")
##03 Differential expression analysis at gene level
gse <- summarizeToGene(se)

#transcript
Tcount <- assays(se)$count
Ttpm <- assays(se)$abundance
#gene    
Gcount <- assays(gse)$count
Gtpm <- assays(gse)$abundance

fwrite(as_tidytable(Tcount, .keep_rownames="TXID"),paste0("Tcount_", outName, ".txt"), sep = "\t")
fwrite(as_tidytable(Ttpm, .keep_rownames="TXID"),paste0("Ttpm_", outName, ".txt"), sep = "\t")
fwrite(as_tidytable(Gcount, .keep_rownames="GENEID"),paste0("Gcount_", outName, ".txt"), sep = "\t")
fwrite(as_tidytable(Gtpm, .keep_rownames="GENEID"),paste0("Gtpm_", outName, ".txt"), sep = "\t")

cat("All done!", "\n")
sessionInfo()
