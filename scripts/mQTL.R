# Auxiliary functions
##################################################################################################
extract_features <- function(file_name, new_name){
  met_data <- read.table(file_name,header = TRUE,stringsAsFactors = FALSE)
  features <- as.vector(features$V1)
  features <- as.character(features)
  if(length(intersect(rownames(met_data),features))==0){
    rownames(met_data)<-met_data[,1]
  }
  else {
    features <- intersect(rownames(met_data),features)
  }
  met_data <- met_data[features,]
  write.table(met_data,new_name,sep="\t",quote = FALSE)
}
########################
extract_features_results <- function(file_name, new_name){
  met_data <- read.table(file_name,header = TRUE,stringsAsFactors = FALSE)
  features <- read.table("/home/natalja/Documents/EMIF/KCL/mQTL/mQTL_data/features.txt",stringsAsFactors = FALSE)
  features <- as.vector(features$V1)
  features <- as.character(features)
  features <- intersect(met_data$metabolite,features)
  if (length(features)>0){
    met_data <- met_data[met_data$metabolite %in% features,]
    write.table(met_data,new_name,sep="\t",quote = FALSE)
  }
  else{
    print("No such features!")
  }
}
########################
match_matrices <- function(matrix1_name, matrix2_name, matrix2_new_name){
  matrix1 <- read.table(matrix1_name, header = TRUE,stringsAsFactors = FALSE ,skip=0, nrows=1)
  matrix2 <- read.table(matrix2_name,header = TRUE,stringsAsFactors = FALSE)
  matrix2 <- matrix2[,intersect(colnames(matrix1),colnames(matrix2))]
  samples <- colnames(matrix2) 
  matrix2 <- cbind(rownames(matrix2),matrix2)
  colnames(matrix2) <- c("id",samples)
  write.table(matrix2,matrix2_new_name,sep="\t",quote = FALSE,row.names = FALSE)
}
########################
check_matrices <- function(dna_matrix_name,met_matrix_name,cov_name){
  matrix1 <- read.table(dna_matrix_name, header = TRUE,stringsAsFactors = FALSE ,skip=0, nrows=1)
  matrix2 <- read.table(met_matrix_name, header = TRUE,stringsAsFactors = FALSE ,skip=0, nrows=1)
  matrix3 <- read.table(cov_name, header = TRUE,stringsAsFactors = FALSE ,skip=0, nrows=1)
  a <- all.equal(colnames(matrix1)[-1],colnames(matrix2)[-1])
  b <- all.equal(colnames(matrix2)[-1],colnames(matrix3)[-1])
  c <- all.equal(colnames(matrix1)[-1],colnames(matrix3)[-1])
  return(c(a,b,c))
}
########################
# Bonferroni correction
bonferroni_correction <- function(file_name){
  mqtl <- read.table(file_name,header = TRUE,stringsAsFactors = FALSE)
  dim(mqtl)
  mqtl$bonf <- p.adjust(mqtl$p.value,"bonferroni")
  mqtl <- mqtl[mqtl$bonf<0.01,]
  dim(mqtl)
  names(mqtl) <- c("SNP","metabolite","beta", "t.stat","p.value" ,"FDR","bonf") 
  mqtl <- mqtl[-c(6)]
  write.table(mqtl,file_name,sep="\t",row.names = FALSE)
}
########################
# Permutation
permute_matrix <- function(file_input_name,file_output_name){
  met_data <- read.table(file_input_name,header = TRUE,stringsAsFactors = FALSE)
  met_data2 <- met_data[,c(1,sample(2:ncol(met_data)))]
  colnames(met_data2) <- colnames(met_data)
  write.table(met_data2,file_output_name,sep="\t",row.names = FALSE)
}
##################################################################################################

##################################################################################################
# NORMALIZATION
##################################################################################################
# EigenMS normalization
##################################################################################################
EigenMS_normalize <- function(file_name, output_file_name){
  base.dir = "/home/natalja/Documents/EMIF/KCL/"
  setwd(base.dir)
  # cov1 <- read.table("./mQTL/Covariates_ADC_vs_CTR.tsv",header = TRUE,stringsAsFactors = FALSE)
  # cov2 <- read.table("./mQTL/Covariates_cMCI_vs_sMCI.tsv",header = TRUE,stringsAsFactors = FALSE)
  # cov1 <- t(cov1)
  # cov2 <- t(cov2)
  # colnames(cov1) <- cov1[1,]
  # cov1 <- cov1[-1,]
  # colnames(cov2) <- cov2[1,]
  # cov2 <- cov2[-1,]
  # cov2 <- as.data.frame(cov2)
  # cov2$group <- gsub("0","2",cov2$group) # sMCI - 2
  # cov2$group <- gsub("1","3",cov2$group) # cMCI - 3
  # cov <- rbind(cov1,cov2)
  
  cov <- read.table("/home/natalja/Documents/EMIF/KCL/mQTL/mQTL_data/Covariates_EigenMS.txt",header = T,sep="\t")
  rownames(cov) <- cov$Sample
  
  ddata = read.table(file_name, header=TRUE, row.names=1, sep=",") 
  #ddata <- ddata[,-(2:5)]
  ddata <- ddata[,intersect(colnames(ddata),rownames(cov))] 
  ddata <- ddata[,order(intersect(colnames(ddata),rownames(cov)))]
  
  #grps = read.table("./mQTL_data/Covariates_ADC_vs_CTR_SHPOS.tsv", header=TRUE, row.names=1) 
  #ddata <- ddata[,colnames(grps)]
  
  cov <- cov[intersect(colnames(ddata),rownames(cov)),]
  cov <- cov[order(intersect(colnames(ddata),rownames(cov))),]
  grps2 <- as.factor(cov$EigenMS_factor)
  #names(grps2) <- cov$name
  
  dim(ddata) # 2802 metabolites x 227 samples
  scaleShift=abs(min(ddata, na.rm = TRUE))+1  
  m_logInts = log(ddata+scaleShift)
  
  # if not on log scale and if has 0's for missing values, uncomment the following 3 lines:
  # m_logInts[m_logInts==0] = NA  #  4.3% missing values, remove 0's, replace with NA's
  # m_logInts = log2(m_logInts)
  m_nummiss = sum(is.na(m_logInts)) # 2000 total values, 700 missing values
  
  # plot boxplots for each sample
  #par(mar=c(10,3,3,3)) # allows to have nice vertical labels!!! this provides room for them
  #par(mfcol=c(1,1))
  #boxplot(m_logInts, las=2) 
  
  m_numtot = dim(m_logInts)[1] * dim(m_logInts)[2] # 8000 total observations
  m_percmiss = m_nummiss/m_numtot  # 8.75% percent missing observations
  # plot numbr of missing values for each sample
  #par(mfcol=c(1,1))
  #barplot(colSums(is.na(m_logInts)), main="Numbers of missing values in samples (grp order)") #570 x 600
  
  # source in the EigenMS and correlation heatmap functions
  source("/home/natalja/Documents/R_projects/EigenMS.R")
  
  # define parameter prot.info, 2 column data frame with IDs for metabolites or peptides
  # in case of matabolites the 2 columns are identical. For peptides 2nd column can contain protein IDs
  m_prot.info = cbind(rownames(ddata),rownames(ddata)) # all unique metabolite/peptide IDs, otherwise not possible to have as row names
  
  # set up treatent group information
  #grps <- t(grps)
  #grps <- as.data.frame(grps)
  #grps2 <- as.factor(grps$group)
  
  #grps2 <- as.factor(cov$group)
  
  # Example part 1 - single factor
  m_ints_eig1 = eig_norm1(m=m_logInts,treatment=grps2,prot.info=m_prot.info,write_to_file = output_file_name)
  m_ints_norm1 = eig_norm2(rv=m_ints_eig1) 
  #par(mfcol=c(1,1))
  #boxplot(m_ints_norm1$norm_m, las=2)
  
  write.table(m_ints_norm1$norm_m,output_file_name,sep="\t")
}
##################################################################################################
# Quantile Normalization of metabolomic data
##################################################################################################
QN_normalize <- function(file_name, output_file_name){
  met <-read.table(file_name,header=T,sep="\t")
  for (i in 1:nrow(met)){
    mat <- met[i,]
    mat = t(apply(mat, 1, rank, ties.method = "average"));
    mat = qnorm(mat / (ncol(met)+1))
    met[i,] <- mat 
  }
  write.table(met,output_file_name,sep="\t")
}
##################################################################################################


##################################################################################################
# MatrixEQTL
##################################################################################################
# load snps 
library(MatrixEQTL)
SNP_file_name <- "/home/natalja/Documents/EMIF/KCL/mQTL/mQTL_data/dna_matrix/dna_matrix_imputed_SHPOS.tsv"
snps = SlicedData$new();
snps$fileDelimiter = "\t"; # the TAB character
snps$fileOmitCharacters = 'NA' ;# denote missing values;
snps$fileSkipRows = 1; # one row of column labels
snps$fileSkipColumns = 1; # one column of row labels
snps$fileSliceSize = 10000; # read file in pieces of 10,000 rows
snps$LoadFile(SNP_file_name);

# Outliers in genotype. Minor allele frequency filtering (was 10 percent)
maf.list = vector('list', length(snps))
for(sl in 1:length(snps)) {
  slice = snps[[sl]];
  maf.list[[sl]] = rowMeans(slice,na.rm=TRUE)/2;
  maf.list[[sl]] = pmin(maf.list[[sl]],1-maf.list[[sl]]);
}
maf = unlist(maf.list)

## Look at the distribution of MAF
# hist(maf[maf<0.1],seq(0,0.1,length.out=100))

cat('SNPs before filtering:',nrow(snps))
# snps$RowReorderSimple(maf>0.05);
snps$RowReorder(maf>0.05);
cat('SNPs before filtering:',nrow(snps))

# Run analysis
##################################################################################################
mqtl_analysis <- function(met_matrix_name,cov_name,output_file_name){ 
  require(MatrixEQTL)
  #base.dir = "/home/natalja/Documents/EMIF/KCL/mQTL/"
  #setwd(base.dir)

  useModel = modelLINEAR
  #SNP_file_name = paste(base.dir, dna_matrix_name, sep="")
  #snps_location_file_name = paste(base.dir, "/dna_matrix_pos.tsv", sep="");
  metabolomics_file_name = met_matrix_name #paste(base.dir, met_matrix_name, sep="")
  #gene_location_file_name = paste(base.dir, "/metloc.txt", sep="");
  covariates_file_name = cov_name #paste(base.dir, cov_name, sep="")
  
  pvOutputThreshold = 1e-5;
  errorCovariance = numeric();
  min.pv.by.genesnp = TRUE 
  
  tic_load = proc.time()[3];
  
  #snps = SlicedData$new();
  #snps$fileDelimiter = "\t"; # the TAB character
  #snps$fileOmitCharacters = 'NA' ;# denote missing values;
  #snps$fileSkipRows = 1; # one row of column labels
  #snps$fileSkipColumns = 1; # one column of row labels
  #snps$fileSliceSize = 10000; # read file in pieces of 10,000 rows
  #snps$LoadFile(SNP_file_name);
  
  ## Load metabolomics data
  
  metabolites = SlicedData$new();
  metabolites$fileDelimiter = '\t'; # the TAB character
  metabolites$fileOmitCharacters = 'NA'; # denote missing values;
  metabolites$fileSkipRows = 1; # one row of column labels
  metabolites$fileSkipColumns = 1; # one column of row labels
  metabolites$fileSliceSize = 10000; # read file in pieces of 10,000 rows
  metabolites$LoadFile(metabolomics_file_name);
  
  ## Load covariates
  
  cvrt = SlicedData$new();
  cvrt$fileDelimiter = '\t'; # the TAB character
  cvrt$fileOmitCharacters = 'NA'; # denote missing values;
  cvrt$fileSkipRows = 1; # one row of column labels
  cvrt$fileSkipColumns = 1; # one column of row labels
  cvrt$fileSliceSize = snps$nCols()+1; # read file in one piece
  if(length(covariates_file_name)>0) {
    cvrt$LoadFile(covariates_file_name);
  }
  
  toc_load = proc.time()[3];
  #cat('eQTL time: ', toc_load-tic_load, ' sec\n');
  
  ## Run the analysis
  {
    tic_mqtl = proc.time()[3];
    me=Matrix_eQTL_engine(snps, metabolites,     cvrt,output_file_name,pvOutputThreshold,useModel, errorCovariance, verbose=TRUE,pvalue.hist = "qqplot");
    # pvalue.hist = 100
    toc_mqtl = proc.time()[3];
  }
  #plot(me, col="grey")
  plot(me, pch = 16, cex = 0.7)
  
}  
##################################################################################################

bases <- c("SLNEG","SLPOS","UHPOS","URNEG","URPOS") # "SHPOS","SLNEG","SLPOS","UHPOS","URNEG","URPOS"
# Normalize all metabolic data
for(base in bases){
  #met <- paste("/home/natalja/Documents/EMIF/KCL/NPC/NPC/original_data/",base,"_original.csv",sep="")
  met <- paste("/home/natalja/Documents/EMIF/KCL/NPC/NPC/normalized_data/MetNormalizer/",base,".csv",sep="")
  met_EigenMS <- paste("/home/natalja/Documents/EMIF/KCL/NPC/NPC/normalized_data/EigenMS/",base,"_EigenMSnormalized_Met.tsv",sep="")
  met_EigenMS_qn <- paste("/home/natalja/Documents/EMIF/KCL/NPC/NPC/normalized_data/EigenMS/",base,"_EigenMSnormalized_Met_qn.tsv",sep="")
  
  EigenMS_normalize(met,met_EigenMS)
  QN_normalize(met_EigenMS,met_EigenMS_qn)
}

# Match matrices (genomic, metabolomic and covariates)
genomic_matrix_all <- "dna_matrix_10percentMAF.tsv"
for(base in bases){
  # source filter_columns
  source("/home/natalja/Documents/R_projects/filter_columns.R")
  
}
################################################################################


match_matrices("/home/natalja/Documents/EMIF/KCL/mQTL/mQTL_data/dna_matrix/dna_matrix_imputed_SHPOS.tsv",
               "/home/natalja/Documents/EMIF/KCL/NPC/NPC/normalized_data/EigenMS/SHPOS_EigenMSnormalized_qn.tsv",
               "/home/natalja/Documents/EMIF/KCL/NPC/NPC/normalized_data/EigenMS_genotype_matched/SHPOS_EigenMSnormalized_qn_gm.tsv")

check_matrices("/home/natalja/Documents/EMIF/KCL/mQTL/mQTL_data/dna_matrix/dna_matrix_imputed_SHPOS.tsv",
               "/home/natalja/Documents/EMIF/KCL/NPC/NPC/normalized_data/EigenMS_genotype_matched/SHPOS_EigenMSnormalized_gm.tsv",
               "/home/natalja/Documents/EMIF/KCL/mQTL/mQTL_data/Covariates_SHPOS.tsv")

met_matrix_name <-"/home/natalja/Documents/EMIF/KCL/NPC/NPC/normalized_data/EigenMS_genotype_matched/SHPOS_EigenMSnormalized_qn_gm.tsv"
cov_name <- "/home/natalja/Documents/EMIF/KCL/mQTL/mQTL_data/Covariates/Covariates_all_SHPOS.txt"
output_file_name <- "/home/natalja/Documents/EMIF/KCL/mQTL/mQTL_data/results/mQTL_EigenMS_qn_SHPOS.txt"
mqtl_analysis(met_matrix_name,cov_name,output_file_name)


permuted_met_matrix_name <- paste(gsub(".tsv","",met_matrix_name),"_permuted.tsv",sep="")
permute_matrix(met_matrix_name,permuted_met_matrix_name)
permuted_output_file_name <- "/home/natalja/Documents/EMIF/KCL/mQTL/mQTL_data/results/mQTL_EigenMS_permuted_qn_SHPOS.txt"

mqtl_analysis(permuted_met_matrix_name,cov_name,permuted_output_file_name)

library(plyr)
cov <- read.table("/home/natalja/Documents/EMIF/KCL/sample_info.tsv",sep="\t",header = T, stringsAsFactors = F)
cov <- cov[,1:7]
rownames(cov) <- cov$Sample
cov <-cov[,-1]
cov$Gender <- revalue(cov$Gender, c("Male"="0", "Female"="1"))
#cov$Centre <- revalue(cov$Centre, c("London"="1", "Lodz"="2","Kuopio"="3","Perugia"="4","Thessaloniki"="5","Toulouse"="6"))
cov$Cohort <- revalue(cov$Cohort, c("Addneuromed"="1", "DCR"="2","ART"="3"))

dna_samples <- read.table("/home/natalja/Documents/EMIF/KCL/mQTL/mQTL_data/dna_matrix/samples.txt",sep="\t",header = T)
cov <- cov[rownames(cov) %in% dna_samples$Sample,]

cov$Cohort <- as.numeric(cov$Cohort) - 1
cov <- cbind(cov,with(cov, model.matrix(~ Centre + 0)))
cov <-cov[,-4]

cov$Genotype_batch <- as.character(cov$Genotype_batch)
cov <- cbind(cov,with(cov, model.matrix(~ Genotype_batch + 0)))
cov <-cov[,-5]

# check if all values are the same 
length(unique(as.list(cov))) == 1

write.table(t(cov),"/home/natalja/Documents/EMIF/KCL/mQTL/mQTL_data/Covariates_all.txt",sep="\t",quote = F)

match_matrices("/home/natalja/Documents/EMIF/KCL/mQTL/mQTL_data/dna_matrix/dna_matrix_imputed_SHPOS.tsv",
               "/home/natalja/Documents/EMIF/KCL/NPC/NPC/normalized_data/EigenMS/met_matrix.tsv",
               "/home/natalja/Documents/EMIF/KCL/NPC/NPC/normalized_data/EigenMS/met_matrix_gm.tsv")

match_matrices("/home/natalja/Documents/EMIF/KCL/NPC/NPC/normalized_data/EigenMS/met_matrix_gm.tsv",
               "/home/natalja/Documents/EMIF/KCL/mQTL/mQTL_data/Covariates/Covariates_all.txt",
               "/home/natalja/Documents/EMIF/KCL/mQTL/mQTL_data/Covariates_all_met.txt")

# Combine matrices
met_SHPOS <- read.table("/home/natalja/Documents/EMIF/KCL/NPC/NPC/normalized_data/EigenMS/SHPOS_EigenMSnormalized.tsv",sep="\t",header = T)
met_SLNEG <- read.table("/home/natalja/Documents/EMIF/KCL/NPC/NPC/normalized_data/EigenMS/SLNEG_EigenMSnormalized.tsv",sep="\t",header = T)
met_SLPOS <- read.table("/home/natalja/Documents/EMIF/KCL/NPC/NPC/normalized_data/EigenMS/SLPOS_EigenMSnormalized.tsv",sep="\t",header = T)
met_UHPOS <- read.table("/home/natalja/Documents/EMIF/KCL/NPC/NPC/normalized_data/EigenMS/UHPOS_EigenMSnormalized.tsv",sep="\t",header = T)
met_URNEG <- read.table("/home/natalja/Documents/EMIF/KCL/NPC/NPC/normalized_data/EigenMS/URNEG_EigenMSnormalized.tsv",sep="\t",header = T)
met_URPOS <- read.table("/home/natalja/Documents/EMIF/KCL/NPC/NPC/normalized_data/EigenMS/URPOS_EigenMSnormalized.tsv",sep="\t",header = T)

rownames(met_SHPOS) <- paste("SHPOS_",rownames(met_SHPOS),sep="")
rownames(met_SLNEG) <- paste("SLNEG_",rownames(met_SLNEG),sep="")
rownames(met_SLPOS) <- paste("SLPOS_",rownames(met_SLPOS),sep="")
rownames(met_UHPOS) <- paste("UHPOS_",rownames(met_UHPOS),sep="")
rownames(met_URNEG) <- paste("URNEG_",rownames(met_URNEG),sep="")
rownames(met_URPOS) <- paste("URPOS_",rownames(met_URPOS),sep="")

serum_samples <- intersect(intersect(colnames(met_SHPOS),colnames(met_SLPOS)),colnames(met_SLNEG)) # 578
urine_samples <- intersect(intersect(colnames(met_UHPOS),colnames(met_URPOS)),colnames(met_URNEG)) # 549
all_samples <- intersect(serum_samples,urine_samples) # 416

library(plyr)
met_matrix <- rbind.fill(met_SHPOS,met_SLNEG)
rownames(met_matrix) <- c(rownames(met_SHPOS),rownames(met_SLNEG))
met_matrix2 <- rbind.fill(met_matrix,met_SLPOS)
rownames(met_matrix2) <- c(rownames(met_matrix),rownames(met_SLPOS))
met_matrix <- rbind.fill(met_matrix2,met_UHPOS)
rownames(met_matrix) <- c(rownames(met_matrix2),rownames(met_UHPOS))
met_matrix2 <- rbind.fill(met_matrix,met_URNEG)
rownames(met_matrix2) <- c(rownames(met_matrix),rownames(met_URNEG))
met_matrix <- rbind.fill(met_matrix2,met_URPOS)
rownames(met_matrix) <- c(rownames(met_matrix2),rownames(met_URPOS))

met_matrix2 <- met_matrix[,colSums(is.na(met_matrix)) == 0] 
# 129 samples
met_matrix2 <- met_matrix[,colnames(met_matrix) %in% all_samples]

met_matrix <- met_matrix2[,colnames(met_matrix2) %in% serum_samples]
serum_matrix<- log2(met_matrix)
# Mean and SD
met_matrix3 <- transform(met_matrix, SD=apply(met_matrix2,1, sd, na.rm = TRUE),MEAN=apply(met_matrix2,1, mean, na.rm = TRUE)) 
serum_stat <- met_matrix3[,579:580]
serum_stat$genetic_coefficient_of_variation <- serum_stat$SD/serum_stat$MEAN

write.table(serum_matrix,"/home/natalja/Documents/EMIF/serum_matrix.tsv",sep="\t",quote = F,row.names = T)

write.table(serum_matrix,"/home/natalja/Documents/EMIF/serum_matrix.tsv",sep="\t",quote = F,row.names = T)

# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}
sm <- t(serum_matrix)
met_corr <- rcorr(sm, type="pearson")
flattenCorrMatrix(met_corr$r, met_corr$P)
write.table(met_corr$r,"/home/natalja/Documents/EMIF/serum_corr.tsv",sep="\t",quote = F,row.names = T)

################################################################################
# Plot how normalization worked
################################################################################
library(fitdistrplus)
library(logspline)
quantile_normalize <- function(x,n){
  mat <- x
  mat = t(apply(mat, 1, rank, ties.method = "average"));
  mat = qnorm(mat / (n+1))
  return (mat)
}


base <- "NMR_Lipidomics_Positive_Urine"
met_orig <-read.table(paste("/home/natalja/Documents/EMIF/KCL/NPC/NPC/NMR/original_data/",base,"_original.txt",sep=""),header=T,sep="\t")
met <-read.table(paste("/home/natalja/Documents/EMIF/KCL/NPC/NPC/NMR/original_data//",base,"_qn_gm.txt",sep=""),header=T,sep="\t")

i<-2
x_orig <- met_orig[i,]
x_eigenMS <- met[i,]
x_orig_qn <- quantile_normalize(x_orig,ncol(met_orig))
x_eigenMS_qn <- quantile_normalize(x_eigenMS,ncol(met))


i<-1
x_orig <- nmr_orig[i,]
x_norm <- nmr_vsn[i,]
x_norm2 <- nmr_pqn[i,]
descdist(as.numeric(x_orig[]), discrete = FALSE)
descdist(as.numeric(x_norm[]), discrete = FALSE)
descdist(as.numeric(x_norm2[]), discrete = FALSE)

png("/home/natalja/Documents/EMIF/images/met_data_orig_distribution.png")
descdist(as.numeric(x_orig[]), discrete = FALSE)
dev.off()
png("/home/natalja/Documents/EMIF/images/met_data_orig_qn_distribution.png")
descdist(as.numeric(x_orig_qn[]), discrete = FALSE)
dev.off()
png("/home/natalja/Documents/EMIF/images/met_data_eigenMS_distribution.png")
descdist(as.numeric(x_eigenMS[]), discrete = FALSE)
dev.off()
png("/home/natalja/Documents/EMIF/images/met_data_eigenMS_qn_distribution.png")
descdist(as.numeric(x_eigenMS_qn[]), discrete = FALSE)
dev.off()

png("/home/natalja/Documents/EMIF/images_mQTL/met_data_orig_fit_norm_NMR.png")
fit.norm <- fitdist(as.numeric(x_orig[]), "norm")
plot(fit.norm)
dev.off()
png("/home/natalja/Documents/EMIF/images_mQTL/met_data_orig_qn_fit_norm_NMR.png")
fit.norm <- fitdist(as.numeric(x_orig_qn[]), "norm")
plot(fit.norm)
dev.off()
png("/home/natalja/Documents/EMIF/images/met_data_eigenMS_fit_norm.png")
fit.norm <- fitdist(as.numeric(x_eigenMS[]), "norm")
plot(fit.norm)
dev.off()
png("/home/natalja/Documents/EMIF/images/met_data_eigenMS_qn_fit_norm.png")
fit.norm <- fitdist(as.numeric(x_eigenMS_qn[]), "norm")
plot(fit.norm)
dev.off()
################################################################################

# Circos
base <- "SHPOS"
threshold <- 0.01

res <- read.table(paste("/home/natalja/Documents/EMIF/KCL/mQTL/mQTL_data/results/mQTL_EigenMS_qn_",base,".txt",sep=""),header=T,sep="\t")
res_sig <- res[res$FDR<threshold,]
res_sig$SNP <- gsub('_2','',res_sig$SNP)
res_sig$SNP <- gsub('_1','',res_sig$SNP)



library(biomaRt)
snp.db <- useMart("ENSEMBL_MART_SNP", dataset="hsapiens_snp")

nt.biomart <- getBM(c("refsnp_id","allele","chr_name","chrom_start","ensembl_gene_stable_id"),
                    filters="snp_filter",
                    values=res_sig$SNP,
                    mart=snp.db)
genes <- unique(nt.biomart$ensembl_gene_stable_id)

ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
genes_annot <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol','chromosome_name','description'), 
                     filters = 'ensembl_gene_id', values = genes, mart = ensembl)

res_annot <- merge(res_sig,nt.biomart,by.x="SNP",by.y="refsnp_id",all.x=T,all.y=F)
res_annot <- merge(res_annot,genes_annot,by.x="ensembl_gene_stable_id",by.y="ensembl_gene_id",all.x=T,all.y=F)

res_annot$nr <- 0
res <- res_annot[0,]

for(metabolite in as.character(unique(res_annot$gene))){
  met_group <- res_annot[res_annot$gene==metabolite,]
  for(chr in as.character(unique(met_group$chr_name))){
    if(!(is.na(chr))){
      met_chr_group <-  met_group[met_group$chr_name==chr,]
      snp <- met_chr_group[which.min(met_chr_group$FDR),]
      snp$nr <- nrow(met_chr_group)
      res <- rbind(res,snp)
    }
  }
}
write.table(res,"/home/natalja/Documents/EMIF/KCL/mQTL/mQTL_data/results/SHPOS_result_annot.txt",sep="\t",row.names = FALSE)
write.table(res_annot,"/home/natalja/Documents/EMIF/KCL/mQTL/mQTL_data/results/SHPOS_result_annot_full.txt",sep="\t",row.names = FALSE)

library(RCircos);
data(UCSC.HG19.Human.CytoBandIdeogram);
cyto.info <- UCSC.HG19.Human.CytoBandIdeogram;
chr.exclude <- NULL;
num.inside <- 5;
RCircos.Set.Core.Components(cyto.info,chr.exclude, num.inside, num.outside);


# Covariates again
cov <-  read.table("/home/natalja/Documents/EMIF/KCL/mQTL/mQTL_data/Covariates/cov.dat",sep="\t",header=T)
cov_adc <- cov[!(cov$Diagnosis %in% c("CTL","sMCI")),]
rownames(cov_adc) <- cov_adc$id
cov_adc2 <- t(cov_adc)
write.table(cov_adc2,"/home/natalja/Documents/EMIF/KCL/mQTL/mQTL_data/Covariates_test_SHPOS_disease.txt",sep="\t",quote = F)


#k-means
# ? repeat 100 times with random starting points, average the result
cov <-  read.table("/home/natalja/Documents/EMIF/KCL/pheno.txt",sep="\t",header=T)
met_qn <-read.table(paste("/home/natalja/Documents/EMIF/KCL/NPC/NPC/normalized_data/EigenMS/",base,"_EigenMSnormalized_qn.tsv",sep=""),header=T,sep="\t")

rownames(cov) <- cov$id
cov <- cov[,3:5]
cov <-cov[!is.na(cov$Diagnosis),]
cov <- cov[cov$Diagnosis %in% c("ADC","CTL"),]
e <- met[,names(met) %in% rownames(cov)]

e <- met_qn[,names(met_qn) %in% rownames(cov)]
cov <- cov[rownames(cov) %in% colnames(e),]
x <- t(e)
wss <- (nrow(x)-1)*sum(apply(x,2,var))
for (i in 2:20) wss[i] <- sum(kmeans(x,centers=i)$withinss)
plot(1:20, wss, type="b", xlab="Number of Clusters",ylab="Within groups sum of squares")

(kc <- kmeans(x, 2)) 

table(cov$Diagnosis, kc$cluster)

write.table(kc$cluster,"/home/natalja/Documents/EMIF/KCL/mQTL/mQTL_data/results/cluster.txt",sep="\t")

# heatmap
library(RColorBrewer)
library(genefilter)
e <- as.matrix(e)
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
rv <- rowVars(e) #a genefilter function - extract 
head(e)
head(rv)
idx <- order(-rv)[1:1000]  #get 40 most variable genes
heatmap(t(e[idx,]),col=hmcol,labRow  = cov$Diagnosis)
cols <- palette(brewer.pal(8, "Dark2"))[as.numeric(cov$Diagnosis)]
heatmap.2(e,trace = 'none', dendrogram='none',ColSideColors=cols,labCol = cov$Diagnosis)

# heatmap with mqtl results
cov <-  read.table("/home/natalja/Documents/EMIF/KCL/pheno.txt",sep="\t",header=T)
met_qn <-read.table(paste("/home/natalja/Documents/EMIF/KCL/NPC/NPC/normalized_data/EigenMS/",base,"_EigenMSnormalized_qn.tsv",sep=""),
                    header=T,sep="\t")
res <- read.table(paste("/home/natalja/Documents/EMIF/KCL/mQTL/mQTL_data/results/result_",base,"_annot.txt",sep=""),
                  header=T,sep="\t",stringsAsFactors = F)

rownames(cov) <- cov$FID
cov <-cov[!is.na(cov$Diagnosis),]
cov <- cov[cov$Diagnosis %in% c("ADC","CTL"),]
e <- met_qn[rownames(met_qn) %in% res$gene,names(met_qn) %in% rownames(cov)]
cov <- cov[rownames(cov) %in% colnames(e),]
cov <- cov[colnames(e),]

library(RColorBrewer)
e <- as.matrix(e)
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
cols <- palette(brewer.pal(8, "Accent"))[as.numeric(as.factor(cov$Diagnosis))]
heatmap(t(e),col=hmcol,labRow  = cov$Diagnosis,RowSideColors=cols,
        Colv=FALSE,Rowv=FALSE,hclustfun=function(d) hclust(d, method="ward.D"))
library(gplots)
heatmap.2(e,col = hmcol,dendrogram="none",trace="none",scale="row",
          Rowv=FALSE,ColSideColors=cols,labCol = cov$Diagnosis)

legend("bottom",      
       legend = unique(as.factor(cov$Diagnosis)),
       col = unique(as.numeric(as.factor(cov$Diagnosis))), 
       lty= 1,             
       lwd = 5,           
       cex=.7
)

# Save table of metabolites that have genomic associations
#########################################################################################################################################
cov <-  read.table("/home/natalja/Documents/EMIF/KCL/pheno.txt",sep="\t",header=T)
met_qn <-read.table(paste("/home/natalja/Documents/EMIF/KCL/NPC/NPC/normalized_data/EigenMS/",base,"_EigenMSnormalized_qn.tsv",sep=""),
                    header=T,sep="\t")
res <- read.table(paste("/home/natalja/Documents/EMIF/KCL/mQTL/mQTL_data/results/result_",base,"_annot.txt",sep=""),
                  header=T,sep="\t",stringsAsFactors = F)

rownames(cov) <- cov$FID
cov <-cov[!is.na(cov$Diagnosis),]
cov <- cov[cov$Diagnosis %in% c("ADC","CTL"),]
e <- met_qn[rownames(met_qn) %in% res$gene,names(met_qn) %in% rownames(cov)]
cov <- cov[rownames(cov) %in% colnames(e),]
cov <- cov[colnames(e),]

e <- t(e)
e <- as.data.frame(e)
e$diagnosis <- cov$Diagnosis
write.table(e,paste("/home/natalja/Documents/EMIF/KCL/mQTL/mQTL_data/results/metabolites_",base,".txt",sep=""),sep="\t")
#########################################################################################################################################
res <- read.table(paste("/home/natalja/Documents/EMIF/KCL/mQTL/mQTL_data/results/result_",base,"_annot_full.txt",sep=""),
                  header=T,sep="\t",stringsAsFactors = F)
list <- read.table(paste("/home/natalja/Documents/EMIF/KCL/mQTL/mQTL_data/results/metabolites/",base,"_ADC_CTL_metabolites.txt",sep=""),
                                      header=F,sep="\t",stringsAsFactors = F)

res_list <- res[res$gene %in% list$V1,]
