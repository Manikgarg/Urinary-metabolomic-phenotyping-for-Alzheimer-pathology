################################################################################
# Feature selection by Linear Regression with diagnosis in a model 
################################################################################
threshold <- 0.01
base <- "UHPOS"

#res <- read.table(paste("/hps/nobackup/ma/natalja/data/mQTL_data/NMR_Results/previous_results/mqtl_",base,".txt",sep=""),header=T,sep="\t")
dna_matrices_dir<-"/hps/nobackup/ma/natalja/data/mQTL_data/dna_matrices/"
met_matrices_dir<-"/hps/nobackup/ma/natalja/data/mQTL_data/met_matrices/normalization_final/" #/NMR/
cov_matrices_dir<-"/hps/nobackup/ma/natalja/data/mQTL_data/covariates/"

res <- read.table(paste("mqtl_",base,".txt",sep=""),header=T,sep="\t")
res_sig <- res[res$FDR<threshold,]
#res_sig$SNP <- gsub('_2','',res_sig$SNP)
#res_sig$SNP <- gsub('_1','',res_sig$SNP)

# Preparation
dim(res_sig)
length(unique(res_sig$SNP))
length(unique(res_sig$gene))

write.table(unique(res_sig$SNP),paste(base,"_sig_SNP.txt",sep=""),sep="\t",row.names=F,quote=F)
write.table(unique(res_sig$gene),paste(base,"_sig_met.txt",sep=""),sep="\t",row.names=F,quote=F)

print(paste("grep -wFf ",base,"_sig_SNP.txt ",dna_matrices_dir,"dna_matrix_",base,".tsv > dna_",base,"_sig.txt",sep=""))
print(paste("grep -wFf ",base,"_sig_met.txt ",met_matrices_dir,base,".tsv > met_",base,"_sig.txt",sep=""))
print(paste("head -1 ",met_matrices_dir,base,".tsv > ",base,"_header.txt",sep=""))
print(paste("cat ",base,"_header.txt met_",base,"_sig.txt > met_",base,"_sig2.txt",sep=""))
print(paste("cat ",base,"_header.txt dna_",base,"_sig.txt > dna_",base,"_sig2.txt",sep=""))


# Linear Regression
expr <- read.table(paste("met_",base,"_sig2.txt",sep=""),sep="\t",header=T,row.names=1,stringsAsFactors = F)
snps <- read.table(paste("dna_",base,"_sig2.txt",sep=""),sep="\t",header=T,row.names=1,stringsAsFactors = F)
cvrt = read.table(paste(cov_matrices_dir,"Covariates_diagnosis.txt",sep=""),sep="\t",header=T,row.names=1)

cvrt <- cvrt[,colnames(cvrt) %in% colnames(expr),]
cvrt <- cvrt[,colnames(expr)]

res_sig$s2 <- 0
res_sig$Diagnosis_ADC <- 0
res_sig$Diagnosis_CTL <- 0
result <-res_sig[0,]
r <- 1
#common_snps <- read.table("urine_common_snps.txt",sep="\t",header=TRUE,stringsAsFactors = FALSE)
for(index in 1:nrow(res_sig)){
  i <- res_sig[index,]
  s1 <- as.numeric(snps[as.character(i$SNP),])
  e1 <- as.numeric(expr[as.character(i$gene),])
  
  s2<-s1[!is.na(s1)]
  e2<-e1[!is.na(s1)]
  cvrt2<-cvrt[,!is.na(s1)]
  age <- as.numeric(cvrt2[1,])
  age[is.na(age)] <- 75
  
  gender <- as.numeric(cvrt2[2,])
  CentreLodz <- as.numeric(cvrt2[3,])
  CentrePerugia <- as.numeric(cvrt2[4,])
  CentreThessaloniki <- as.numeric(cvrt2[5,])
  CentreToulouse <- as.numeric(cvrt2[6,])
  Genotype_batch1<- as.numeric(cvrt2[7,])
  Genotype_batch2 <- as.numeric(cvrt2[8,])
  Diagnosis_ADC <- as.numeric(cvrt2[9,])
  Diagnosis_cMCI <- as.numeric(cvrt2[10,])
  Diagnosis_CTL <- as.numeric(cvrt2[11,])
  
  lm2 = lm(e2~s2+age+gender+CentreLodz+CentrePerugia+CentreThessaloniki+CentreToulouse+Genotype_batch1+Genotype_batch2+Diagnosis_ADC+Diagnosis_cMCI+Diagnosis_CTL)
  
  res <- data.frame(summary(lm2)$coef)
  significant_coefficient_names <- paste(rownames(res[res[,4]<=0.01,]), collapse="_")
  if (grepl("Diagnosis_ADC",   significant_coefficient_names) || grepl("Diagnosis_CTL",   significant_coefficient_names)){
    #print(as.character(i$SNP))
    #print(as.character(i$gene))
    #print(summary(lm2))
    i$s2 <- res[2,4]
    i$Diagnosis_ADC <- res[11,4]
    i$Diagnosis_CTL <- res[13,4]
    result[r,] <- i
    r <- r + 1
    
     i$SNP <- gsub('_2','',i$SNP)
     i$SNP <- gsub('_1','',i$SNP)
     
     #if (i$SNP %in% common_snps$x){
       lm3 = lm(e2~s2)
       png(paste("./plots_individual_associations_with_diagnosis_URPOS/snp_",as.character(i$SNP),"_",gsub('/','',as.character(i$gene)),".png",sep=""), width=950, height=400)
         plot(e2 ~ jitter(s2),
              col=(s2+1),xaxt="n",xlab="Genotype",ylab="Expression")
         axis(1,at=c(0:2),labels=c("AA","Aa","aa"))
         lines(lm3$fitted ~ s2,type="b",pch=15,col="darkgrey")
       dev.off()
     #}
  }
}


result$SNP <- gsub('_2','',result$SNP)
result$SNP <- gsub('_1','',result$SNP)

# Annotation
library(biomaRt)
snp.db <- useMart("ENSEMBL_MART_SNP", dataset="hsapiens_snp")
nt.biomart <- getBM(c("refsnp_id","chr_name","ensembl_gene_stable_id"),
                    filters="snp_filter",
                    values=result$SNP,
                    mart=snp.db)
res_annot <- merge(result,nt.biomart,by.x="SNP",by.y="refsnp_id",all.x=T,all.y=F)
res_annot <- res_annot[-grep("CHR_", res_annot$chr_name),]
genes <- unique(nt.biomart$ensembl_gene_stable_id)
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
genes_annot <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol','chromosome_name','description'), 
                     filters = 'ensembl_gene_id', values = genes, mart = ensembl)


res_annot <- merge(res_annot,genes_annot,by.x="ensembl_gene_stable_id",by.y="ensembl_gene_id",all.x=T,all.y=F)
write.table(res_annot,paste("diagnosis_result_annot_",base,".txt",sep=""),sep="\t")

write.table(result,paste("diagnosis_result_",base,".txt",sep=""),sep="\t")



################################################################################
# Feature selection by correlation - WEKA 
################################################################################
threshold <- 0.01
drugs <- read.table("/hps/nobackup/ma/natalja/data/mQTL_data/drugs.txt",sep="\t",header = T,stringsAsFactors = F)

# LC-MS
base <- "UHPOS"
main_dir <- "/hps/nobackup/ma/natalja/data/mQTL_data/results/normalization_final/"
output_file_name <- paste(main_dir,"result_",base,"_annot_full.txt",sep="")
res_all <- read.table(output_file_name,header=T,sep="\t")
expr <- read.table(paste("/hps/nobackup/ma/natalja/data/mQTL_data/results/normalization_final/met_",base,"_sig2.txt",sep=""),sep="\t",header=T,row.names=1,stringsAsFactors = F)

res_urine_all <-res_all[0,]
res_urine_all$base <= ""
expr_urine_all <-expr[0,]
#expr_urine_all$base <= ""
for(base in c("UHPOS","URPOS","URNEG")){
  main_dir <- "/hps/nobackup/ma/natalja/data/mQTL_data/results/normalization_final/"
  output_file_name <- paste(main_dir,"mqtl_",base,".txt",sep="")
  res_all <- read.table(output_file_name,header=T,sep="\t",stringsAsFactors = F)
  res_all$base <- base
  expr <- read.table(paste("/hps/nobackup/ma/natalja/data/mQTL_data/results/normalization_final/met_",base,"_sig2.txt",sep=""),sep="\t",header=T,row.names=1,stringsAsFactors = F)
  #expr$base <- base
  res_urine_all <- rbind(res_urine_all,res_all)
  expr_urine_all <- rbind(expr_urine_all,expr)
}
res_sig_lcms <- res_urine_all[res_urine_all$FDR<=threshold,]
res_sig_lcms$SNP <- gsub('_2','',res_sig_lcms$SNP)
res_sig_lcms$SNP <- gsub('_1','',res_sig_lcms$SNP)

# NMR 
base <-"NMR_Nosey_Urine"
expr_nmr <- read.table(paste("/hps/nobackup/ma/natalja/data/mQTL_data/NMR_Results/normalization_final2/met_",base,"_sig2.txt",sep=""),sep="\t",header=T,row.names=1,stringsAsFactors = F)
res <- read.table(paste("/hps/nobackup/ma/natalja/data/mQTL_data/NMR_Results/normalization_final2/mqtl_",base,".txt",sep=""),header=T,sep="\t")
res$base <- base
res_sig_nmr <- res[res$FDR<=threshold,]
res_sig_nmr$SNP <- gsub('_2','',res_sig_nmr$SNP)
res_sig_nmr$SNP <- gsub('_1','',res_sig_nmr$SNP)

#write.table(unique(res_sig_lcms$SNP),"/hps/nobackup/ma/natalja/data/mQTL_data/final_results_urine/selected_matrices/snps_lcms.txt",sep="\t",row.names=F,quote=F)
#write.table(unique(res_sig_nmr$SNP),"/hps/nobackup/ma/natalja/data/mQTL_data/final_results_urine/selected_matrices/snps_nmr.txt",sep="\t",row.names=F,quote=F)

# All SNPs
res_all <- rbind(res_sig_lcms, res_sig_nmr)
length(unique(res_all$SNP))
write.table(unique(res_all$SNP),"all_sig_SNP.txt",sep="\t",row.names=F,quote=F)

print(paste("grep -wFf all_sig_SNP.txt ",dna_matrices_dir,"dna_matrix_NMR_Nosey_Urine.tsv > dna_all_sig.txt",sep=""))
print(paste("head -1 ",met_matrices_dir,"NMR_Nosey_Urine.tsv > all_header.txt",sep=""))
print(paste("cat all_header.txt dna_all_sig.txt > dna_all_sig2.txt",sep=""))


# Remove drugs
expr_lcms <- expr_urine_all[!rownames(expr_urine_all) %in% drugs$Metabolite,]
dim(expr_lcms)
# There are not know drugs in NMR
dim(expr_nmr)

# The common samples for NMR and LC-MS
expr_nmr <- expr_nmr[,colnames(expr_nmr) %in% colnames(expr_lcms),]
expr_lcms <- expr_lcms[,colnames(expr_lcms) %in% colnames(expr_nmr),]
expr_nmr <- expr_nmr[,colnames(expr_lcms)]

# Save all selected matrices
#write.table(expr_nmr,"/hps/nobackup/ma/natalja/data/mQTL_data/final_results_urine/selected_matrices/met_nmr.txt",sep="\t",row.names=F,quote=F)
#write.table(expr_lcms,"/hps/nobackup/ma/natalja/data/mQTL_data/final_results_urine/selected_matrices/met_lcms.txt",sep="\t",row.names=F,quote=F)
#write.table(res_sig_lcms,"/hps/nobackup/ma/natalja/data/mQTL_data/final_results_urine/selected_matrices/res_sig_lcms.txt",row.names=F,quote=F)
#write.table(res_sig_nmr,"/hps/nobackup/ma/natalja/data/mQTL_data/final_results_urine/selected_matrices/res_sig_nmr.txt",row.names=F,quote=F)
#write.table(unique(res_sig_lcms$SNP),"/hps/nobackup/ma/natalja/data/mQTL_data/final_results_urine/selected_matrices/snps_lcms.txt",sep="\t",row.names=F,quote=F)
#write.table(unique(res_sig_nmr$SNP),"/hps/nobackup/ma/natalja/data/mQTL_data/final_results_urine/selected_matrices/snps_nmr.txt",sep="\t",row.names=F,quote=F)

################################################################################
# WEKA - data export in CSV format
################################################################################
cvrt = read.table("/hps/nobackup/ma/natalja/data/mQTL_data/covariates/Covariates_all_vertical_text2.txt",sep="\t",header=T,row.names=1)
cvrt <- cvrt[rownames(cvrt) %in% colnames(expr_lcms),]
cvrt <- cvrt[colnames(expr_lcms),]
df <- as.data.frame(t(expr_lcms))
df$Diagnosis <- as.character(cvrt$Diagnosis)
df$Diagnosis2 <- as.character(cvrt$Diagnosis)
df$Diagnosis2[df$Diagnosis2=="cMCI"] <- "ADC"
df$Diagnosis2[df$Diagnosis2=="sMCI"] <- "CTL"
write.csv(df,"ML_data_lcms.csv",row.names = FALSE)

cvrt = read.table("/hps/nobackup/ma/natalja/data/mQTL_data/covariates/Covariates_all_vertical_text2.txt",sep="\t",header=T,row.names=1)
cvrt <- cvrt[rownames(cvrt) %in% colnames(expr_nmr),]
cvrt <- cvrt[colnames(expr_nmr),]
df <- as.data.frame(t(expr_nmr))
df$Diagnosis <- as.character(cvrt$Diagnosis)
write.csv(df,"ML_data_nmr.csv",row.names = FALSE)

expr <- rbind(expr_lcms,expr_nmr)
cvrt = read.table("/hps/nobackup/ma/natalja/data/mQTL_data/covariates/Covariates_all_vertical_text2.txt",sep="\t",header=T,row.names=1)
cvrt <- cvrt[rownames(cvrt) %in% colnames(expr),]
cvrt <- cvrt[colnames(expr),]
df <- as.data.frame(t(expr))
df$Diagnosis <- as.character(cvrt$Diagnosis)
df$Diagnosis2 <- as.character(cvrt$Diagnosis)
df$Diagnosis2[df$Diagnosis2=="cMCI"] <- "ADC"
df$Diagnosis2[df$Diagnosis2=="sMCI"] <- "CTL"
write.csv(df,"ML_data_lcms_nmr.csv",row.names = FALSE)


# filter by selection results
# Original expr size is 1541 metabolites, dna_all - 6923 SNPs
# It's too big matrix for ML procedures - filtering
method <- "LR"
selected_mets <- read.table(paste("/hps/nobackup/ma/natalja/data/mQTL_data/final_results_urine/feature_selection_",method,".txt",sep=""),sep="\t", header=FALSE, stringsAsFactors=F)
ADC_snps <- read.table(paste("/hps/nobackup/ma/natalja/data/mQTL_data/final_results_urine/ADC_snps.txt",sep=""),sep="\t", header=FALSE, stringsAsFactors=F)

expr_selected <- expr[rownames(expr) %in% selected_mets$V1,]
# 393 metabolites
# Now select SNPs for these metabolites
res_all$SNP2 <- res_all$SNP
res_all$SNP2 <- gsub('_2','',res_all$SNP2)
res_all$SNP2 <- gsub('_1','',res_all$SNP2)
res_all_selected <- res_all[res_all$gene %in% selected_mets$V1,]
length(unique(res_all_selected$SNP2))

dna_all <- read.table("/hps/nobackup/ma/natalja/data/mQTL_data/final_results_urine/selected_matrices/dna_all.txt",sep="\t",header=T,row.names=1,stringsAsFactors = F)
# The common samples for NMR and LC-MS
dna_all<- dna_all[,colnames(dna_all) %in% colnames(expr),]
dna_all <- dna_all[,colnames(expr)]
dna_all [is.na(dna_all )] <- 0

# a) Metabolites only
df <- as.data.frame(t(expr))
df$Diagnosis <- as.character(cvrt$Diagnosis)
df$Diagnosis2 <- as.character(cvrt$Diagnosis)
df$Diagnosis3 <- as.character(cvrt$Diagnosis)
df$Diagnosis2[df$Diagnosis2=="cMCI"] <- "ADC"
df$Diagnosis2[df$Diagnosis2=="sMCI"] <- "CTL"
df$Diagnosis3[df$Diagnosis3=="cMCI"] <- "ADC"
df$Diagnosis3[df$Diagnosis3=="sMCI"] <- "ADC"
write.csv(df,"ML_data_a.csv",row.names = FALSE)

df <- as.data.frame(t(expr))
df$Diagnosis <- as.character(cvrt$Diagnosis)
write.csv(df[cvrt$Diagnosis %in% c("CTL","ADC"),],"ML_data_a_ADC_CTL.csv",row.names = FALSE)

# b) Metabolites and SNPs
df <- as.data.frame(cbind(t(expr),t(dna_all)))
df$Diagnosis <- as.character(cvrt$Diagnosis)
df$Diagnosis2 <- as.character(cvrt$Diagnosis)
df$Diagnosis3 <- as.character(cvrt$Diagnosis)
df$Diagnosis2[df$Diagnosis2=="cMCI"] <- "ADC"
df$Diagnosis2[df$Diagnosis2=="sMCI"] <- "CTL"
df$Diagnosis3[df$Diagnosis3=="cMCI"] <- "ADC"
df$Diagnosis3[df$Diagnosis3=="sMCI"] <- "ADC"
write.csv(df,"ML_data_b.csv",row.names = FALSE)

df <- as.data.frame(cbind(t(expr),t(dna_all)))
df$Diagnosis <- as.character(cvrt$Diagnosis)
write.csv(df[cvrt$Diagnosis %in% c("CTL","ADC"),],"ML_data_b_ADC_CTL.csv",row.names = FALSE)

# c) Metabolites, SNPs and covariates
df <- as.data.frame(cbind(t(expr),t(dna_all)))
df$Diagnosis <- as.character(cvrt$Diagnosis)
df$Diagnosis2 <- as.character(cvrt$Diagnosis)
df$Diagnosis3 <- as.character(cvrt$Diagnosis)
df$Diagnosis2[df$Diagnosis2=="cMCI"] <- "ADC"
df$Diagnosis2[df$Diagnosis2=="sMCI"] <- "CTL"
df$Diagnosis3[df$Diagnosis3=="cMCI"] <- "ADC"
df$Diagnosis3[df$Diagnosis3=="sMCI"] <- "ADC"
df$Gender <- cvrt$Gender
df$Centre <- cvrt$Centre
df$Age <- cvrt$Age

df <- as.data.frame(cbind(t(expr),t(dna_all)))
df$Diagnosis <- as.character(cvrt$Diagnosis)
write.csv(df[cvrt$Diagnosis %in% c("CTL","ADC"),],"ML_data_c_ADC_CTL.csv",row.names = FALSE)
write.csv(df,"ML_data_c.csv",row.names = FALSE)

# d) Metabolites and covariates
df <- as.data.frame(t(expr))
df$Gender <- cvrt$Gender
df$Centre <- cvrt$Centre
df$Age <- cvrt$Age
df$Diagnosis <- as.character(cvrt$Diagnosis)
df$Diagnosis2 <- as.character(cvrt$Diagnosis)
df$Diagnosis3 <- as.character(cvrt$Diagnosis)
df$Diagnosis4 <- as.character(cvrt$Diagnosis)
df$Diagnosis2[df$Diagnosis2=="cMCI"] <- "ADC"
df$Diagnosis2[df$Diagnosis2=="sMCI"] <- "CTL"
df$Diagnosis3[df$Diagnosis3=="cMCI"] <- "ADC"
df$Diagnosis3[df$Diagnosis3=="sMCI"] <- "ADC"
df$Diagnosis4[df$Diagnosis4=="cMCI"] <- "MCI"
df$Diagnosis4[df$Diagnosis4=="sMCI"] <- "MCI"

write.csv(df,"ML_data_d.csv",row.names = FALSE)

write.csv(df[cvrt$Diagnosis %in% c("CTL","ADC"),],"ML_data_d_ADC_CTL.csv",row.names = FALSE)

#dna_all_selected <- dna_all[rownames(dna_all) %in% res_all_selected$SNP | rownames(dna_all) %in% ADC_snps$V1,]
# 4131 SNPs
#dna_all_selected <- dna_all_selected[grep('^rs',rownames(dna_all_selected)),]
# 3856 SNPs
#dna_all_selected <- dna_all_selected [complete.cases(dna_all_selected ), ]
# 1365
cvrt = read.table("/hps/nobackup/ma/natalja/data/mQTL_data/covariates/Covariates_all_vertical_text2.txt",sep="\t",header=T,row.names=1)
cvrt <- cvrt[rownames(cvrt) %in% colnames(expr),]
cvrt <- cvrt[colnames(expr),]
df <- as.data.frame(cbind(t(expr),t(dna_all)))
df$Gender <- cvrt$Gender
df$Centre <- cvrt$Centre
df$Age <- cvrt$Age
df$Diagnosis <- as.character(cvrt$Diagnosis)
df$Diagnosis2 <- as.character(cvrt$Diagnosis)
df$Diagnosis2[df$Diagnosis2=="cMCI"] <- "ADC"
df$Diagnosis2[df$Diagnosis2=="sMCI"] <- "CTL"
write.csv(df,"ML_data_lcms_nmr_genotype.csv",row.names = FALSE)
write.csv(df[cvrt$Diagnosis %in% c("CTL","ADC"),],"ML_data_lcms_nmr_genotypeADC_CTL.csv",row.names = FALSE)



method <- "IG"
selected_mets <- read.table(paste("/hps/nobackup/ma/natalja/data/mQTL_data/final_results_urine/feature_selection_",method,".txt",sep=""),sep="\t", header=FALSE, stringsAsFactors=F)

df_selected <- df[,colnames(df) %in% selected_mets$V1 | colnames(df) %in% c("Diagnosis2","Diagnosis") ]
write.csv(df_selected,"ML_data_lcms_nmr_genotype_IG.csv",row.names = FALSE)


# Original met data
met <-read.table("/hps/nobackup/ma/natalja/data/mQTL_data/final_results_urine/selected_matrices/original_met/NMR_Nosey_Urine_original.txt",header=T,sep="\t")
cvrt = read.table("/hps/nobackup/ma/natalja/data/mQTL_data/covariates/Covariates_all_vertical_text2.txt",sep="\t",header=T,row.names=1)
cvrt <- cvrt[rownames(cvrt) %in% colnames(met),]
cvrt <- cvrt[colnames(met),]

e <- met[,names(met) %in% rownames(cvrt)]
cvrt <- cvrt[rownames(cvrt) %in% colnames(e),]
cvrt <- cvrt[colnames(e),]

e_selected <- e[rownames(e) %in% unique(res_sig_nmr$gene),]

e_selected_ADC_CTL <- e_selected[,cvrt$Diagnosis %in% c("CTL","ADC")]
cvrt_ADC_CTL <- cvrt[colnames(e_selected_ADC_CTL),]
df <- as.data.frame(t(e_selected_ADC_CTL ))
df$Diagnosis <- as.character(cvrt_ADC_CTL$Diagnosis)
write.csv(df,"ML_data_nmr_original_ADC_CTL.csv",row.names = FALSE)


df <- as.data.frame(t(e_selected1))
df$Diagnosis <- as.character(cvrt$Diagnosis)


write.csv(df,"ML_data_nmr_original.csv",row.names = FALSE)
################################################################################



################################################################################
# Images for feature selection results
################################################################################
# feature selection by correlation (WEKA result)
library(gplots)
library(RColorBrewer) 
library(ggplot2)
# colors in ggplot2 style
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
cols <- gg_color_hue(4) 
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
hmcol <- brewer.pal(11,"RdBu")
as.fumeric <- function(x,levels=unique(x)) {
  if(!is.character(x)) stop("'x' must be a character")
  as.numeric(factor(x,levels=levels))
}

method <- "corr" # "LR"
selected_mets <- read.table(paste("/hps/nobackup/ma/natalja/data/mQTL_data/final_results_urine/feature_selection_",method,".txt",sep=""),sep="\t", header=FALSE, stringsAsFactors=F)

# Subset from normalized metabolite matrix 
expr_selected_lcms <- expr_lcms[rownames(expr_lcms) %in% selected_mets$V1,]
expr_selected_nmr <- expr_nmr[rownames(expr_nmr) %in% selected_mets$V1,]
expr_selected <- rbind(expr_selected_lcms,expr_selected_nmr)
cvrt = read.table("/hps/nobackup/ma/natalja/data/mQTL_data/covariates/Covariates_all_vertical_text2.txt",sep="\t",header=T,row.names=1)
cvrt <- cvrt[rownames(cvrt) %in% colnames(expr_selected),]
cvrt <- cvrt[colnames(expr_selected),]
df <- as.data.frame(t(expr_selected))
df$Diagnosis <- as.character(cvrt$Diagnosis)

###############################
# Boxplot of selected features 
###############################

df_ADC <- df[df$Diagnosis %in% c("ADC"),]
df_ADC <- df_ADC[,-(ncol(df_ADC))]
res <- data.frame(values = c(as.matrix(df_ADC)), metabolites = rep(colnames(df_ADC),each=nrow(df_ADC)),Diagnosis = "ADC")
df_CTL <- df[df$Diagnosis %in% c("CTL"),]
df_CTL <- df_CTL[,-(ncol(df_CTL))]
res_CTL <- data.frame(values = c(as.matrix(df_CTL)), metabolites = rep(colnames(df_CTL),each=nrow(df_CTL)),Diagnosis = "CTL")
df_sMCI <- df[df$Diagnosis %in% c("sMCI"),]
df_sMCI <- df_sMCI[,-(ncol(df_sMCI))]
res_sMCI <- data.frame(values = c(as.matrix(df_sMCI)), metabolites = rep(colnames(df_sMCI),each=nrow(df_sMCI)),Diagnosis = "sMCI")
df_cMCI <- df[df$Diagnosis %in% c("cMCI"),]
df_cMCI <- df_cMCI[,-(ncol(df_cMCI))]
res_cMCI <- data.frame(values = c(as.matrix(df_cMCI)), metabolites = rep(colnames(df_cMCI),each=nrow(df_cMCI)),Diagnosis = "cMCI")

res <- rbind(res,res_CTL)
res <- rbind(res,res_sMCI)
res <- rbind(res,res_cMCI)

pdf(paste(method,"_boxplot.pdf",sep=""))
p10 <- ggplot(res,aes(x=metabolites,y=values,fill=Diagnosis)) + geom_boxplot(alpha=0.7,outlier.size = 0, coef = 0) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p10
dev.off()
###############################


###############################
# Plot mean values of selected features 
###############################
dfm <- aggregate(x = as.data.frame(t(expr_selected)), by = list(df$Diagnosis), FUN = "mean")
rownames(dfm) <- dfm$Group.1
dfm <- dfm [,-1]
fit <- kmeans(as.data.frame(expr_selected), 5)
dfm_t <- as.data.frame(t(dfm))
dfm_t <- dfm_t[order(fit$cluster),]

cols <- gg_color_hue(4) 
#c("firebrick4","firebrick3","darkgreen","steelblue4")
diagnoses <- c("ADC","CTL","sMCI","cMCI")
ltys <- c(1,1,3,3) 

pdf(paste(method,"_mean_values.pdf",sep=""))
plot(dfm_t$ADC, type="l", col=cols[1],ylim=range(dfm_t), axes=F, ann=T,xlab="", ylab="Mean values",cex.lab=0.8, lwd=2)
grid(lty = 6, col = "cornsilk2") 
axis(1, at=seq(1,nrow(dfm_t),by=1),labels=rownames(dfm_t), las = 2,cex.axis=0.6)
# Plot y axis with smaller horizontal labels 
axis(2, las=1, cex.axis=0.8)
box()
lines(dfm_t$CTL, type="l", lwd=2, col=cols[2])

lines(dfm_t$sMCI, type="l", lty=3, lwd=2, col=cols[3])
lines(dfm_t$cMCI, type="l", lty=3, lwd=2, col=cols[4])

legend("topright", diagnoses, cex=0.8, col=cols, 
       lty=ltys, lwd=2, bty="n");
dev.off()
###############################


###############################
# Heatmap of selected features 
###############################
#re-order levels to correspond to ggplot2 colors and order
cvrt$Diagnosis <- factor(cvrt$Diagnosis, levels = diagnoses)

pdf("heatmap_4D_CF_COND.pdf")
t <- heatmap.2(as.matrix(expr_selected[order(cvrt$Diagnosis)]),col=brewer.pal(11,"RdBu"), 
               trace="none",Colv="none", dendrogram = "row",cexRow=0.5, key=F, main="Metabolites and diagnoses",
               colCol=cols[cvrt[order(cvrt$Diagnosis),]$Diagnosis])

par(lend = 1)           # square line ends for the color legend
legend("topleft",      # location of the legend on the heatmap plot
       legend = diagnoses, # category labels
       col = cols,  # color key
       lty= 1,             # line style
       lwd = 10            # line width
)
dev.off()
################################################################################
# Random forests for feature selection
################################################################################
# More data!
# Metabolomic matrices without sync with genomic data 
# 470 samples instead of 343
expr1 <- read.table("/hps/nobackup/ma/natalja/data/mQTL_data/final_results_urine/selected_matrices/UHPOS_met_sig2.txt",sep="\t",header=T,row.names=1,stringsAsFactors = F)
expr2 <- read.table("/hps/nobackup/ma/natalja/data/mQTL_data/final_results_urine/selected_matrices/URPOS_met_sig2.txt",sep="\t",header=T,row.names=1,stringsAsFactors = F)
expr3 <- read.table("/hps/nobackup/ma/natalja/data/mQTL_data/final_results_urine/selected_matrices/URNEG_met_sig2.txt",sep="\t",header=T,row.names=1,stringsAsFactors = F)
expr_nmr <- read.table("/hps/nobackup/ma/natalja/data/mQTL_data/final_results_urine/selected_matrices/NMR_met_sig2.txt",sep="\t",header=T,row.names=1,stringsAsFactors = F)

expr_nmr <- expr_nmr[,colnames(expr_nmr) %in% colnames(expr1)]
expr1 <- expr1[,colnames(expr1) %in% colnames(expr_nmr)]
expr2 <- expr2[,colnames(expr2) %in% colnames(expr_nmr)]
expr3 <- expr3[,colnames(expr3) %in% colnames(expr_nmr)]

expr_nmr <- expr_nmr[,colnames(expr1)]
expr2 <- expr2[,colnames(expr1)]
expr3 <- expr3[,colnames(expr1)]

expr <- rbind(expr_nmr,expr1)
expr <- rbind(expr,expr2)
expr <- rbind(expr,expr3)

cvrt = read.table("/hps/nobackup/ma/natalja/data/mQTL_data/covariates/Covariates_all_vertical_text2.txt",sep="\t",header=T,row.names=1)
cvrt <- cvrt[rownames(cvrt) %in% colnames(expr),]
cvrt <- cvrt[colnames(expr),]

# a) Metabolites only and d) Mteabolites and covariates
df <- as.data.frame(t(expr))
df$Gender <- cvrt$Gender
df$Centre <- cvrt$Centre
df$Age <- cvrt$Age
df$Diagnosis1 <- as.character(cvrt$Diagnosis)
df$Diagnosis2 <- as.character(cvrt$Diagnosis)
df$Diagnosis3 <- as.character(cvrt$Diagnosis)
df$Diagnosis2[df$Diagnosis2=="cMCI"] <- "ADC"
df$Diagnosis2[df$Diagnosis2=="sMCI"] <- "CTL"
df$Diagnosis3[df$Diagnosis3=="cMCI"] <- "ADC"
df$Diagnosis3[df$Diagnosis3=="sMCI"] <- "ADC"

# a) 
drops <- c("Gender","Centre","Age","Diagnosis2","Diagnosis3")
df1_a <- df[,!(names(df) %in% drops)]
names(df1_a)[names(df1_a) == 'Diagnosis1'] <- 'Diagnosis'

drops <- c("Gender","Centre","Age","Diagnosis1","Diagnosis3")
df2_a <- df[,!(names(df) %in% drops)]
names(df2_a)[names(df2_a) == 'Diagnosis2'] <- 'Diagnosis'

drops <- c("Gender","Centre","Age","Diagnosis1","Diagnosis2")
df3_a <- df[,!(names(df) %in% drops)]
names(df3_a)[names(df3_a) == 'Diagnosis3'] <- 'Diagnosis'

drops <- c("Gender","Centre","Age","Diagnosis2","Diagnosis3")
df4_a <- df[cvrt$Diagnosis %in% c("CTL","ADC"),!(names(df) %in% drops)]
names(df4_a)[names(df4_a) == 'Diagnosis1'] <- 'Diagnosis'

# or set d
drops <- c("Diagnosis2","Diagnosis3")
df1_d <- df[,!(names(df) %in% drops)]
names(df1_d)[names(df1_d) == 'Diagnosis1'] <- 'Diagnosis'

drops <- c("Diagnosis1","Diagnosis3")
df2_d <- df[,!(names(df) %in% drops)]
names(df2_d)[names(df2_d) == 'Diagnosis2'] <- 'Diagnosis'

drops <- c("Diagnosis1","Diagnosis2")
df3_d <- df[,!(names(df) %in% drops)]
names(df3_d)[names(df3_d) == 'Diagnosis3'] <- 'Diagnosis'

drops <- c("Diagnosis2","Diagnosis3")
df4_d <- df[cvrt$Diagnosis %in% c("CTL","ADC"),!(names(df) %in% drops)]
names(df4_d)[names(df4_d) == 'Diagnosis1'] <- 'Diagnosis'

# set b
dna_all <- read.table("/hps/nobackup/ma/natalja/data/mQTL_data/final_results_urine/selected_matrices/dna_all.txt",sep="\t",header=T,row.names=1,stringsAsFactors = F)
# The common samples for NMR and LC-MS
dna_all<- dna_all[,colnames(dna_all) %in% colnames(expr)]
expr_dna <- expr[,colnames(expr) %in% colnames(dna_all)]
dna_all <- dna_all[,colnames(expr_dna)]
dna_all [is.na(dna_all )] <- 0
cvrt = read.table("/hps/nobackup/ma/natalja/data/mQTL_data/covariates/Covariates_all_vertical_text2.txt",sep="\t",header=T,row.names=1)
cvrt <- cvrt[rownames(cvrt) %in% colnames(expr_dna),]
cvrt <- cvrt[colnames(expr_dna),]

df <- as.data.frame(cbind(t(expr_dna),t(dna_all)))
df$Gender <- cvrt$Gender
df$Centre <- cvrt$Centre
df$Age <- cvrt$Age
df$Diagnosis1 <- as.character(cvrt$Diagnosis)
df$Diagnosis2 <- as.character(cvrt$Diagnosis)
df$Diagnosis3 <- as.character(cvrt$Diagnosis)
df$Diagnosis2[df$Diagnosis2=="cMCI"] <- "ADC"
df$Diagnosis2[df$Diagnosis2=="sMCI"] <- "CTL"
df$Diagnosis3[df$Diagnosis3=="cMCI"] <- "ADC"
df$Diagnosis3[df$Diagnosis3=="sMCI"] <- "ADC"

# b) 
drops <- c("Gender","Centre","Age","Diagnosis2","Diagnosis3")
df1_b <- df[,!(names(df) %in% drops)]
names(df1_b)[names(df1_b) == 'Diagnosis1'] <- 'Diagnosis'

drops <- c("Gender","Centre","Age","Diagnosis1","Diagnosis3")
df2_b <- df[,!(names(df) %in% drops)]
names(df2_b)[names(df2_b) == 'Diagnosis2'] <- 'Diagnosis'

drops <- c("Gender","Centre","Age","Diagnosis1","Diagnosis2")
df3_b <- df[,!(names(df) %in% drops)]
names(df3_b)[names(df3_b) == 'Diagnosis3'] <- 'Diagnosis'

drops <- c("Gender","Centre","Age","Diagnosis2","Diagnosis3")
df4_b <- df[cvrt$Diagnosis %in% c("CTL","ADC"),!(names(df) %in% drops)]
names(df4_b)[names(df4_b) == 'Diagnosis1'] <- 'Diagnosis'

# or set c
drops <- c("Diagnosis2","Diagnosis3")
df1_c <- df[,!(names(df) %in% drops)]
names(df1_c)[names(df1_c) == 'Diagnosis1'] <- 'Diagnosis'

drops <- c("Diagnosis1","Diagnosis3")
df2_c <- df[,!(names(df) %in% drops)]
names(df2_c)[names(df2_c) == 'Diagnosis2'] <- 'Diagnosis'

drops <- c("Diagnosis1","Diagnosis2")
df3_c <- df[,!(names(df) %in% drops)]
names(df3_c)[names(df3_c) == 'Diagnosis3'] <- 'Diagnosis'

drops <- c("Diagnosis2","Diagnosis3")
df4_c <- df[cvrt$Diagnosis %in% c("CTL","ADC"),!(names(df) %in% drops)]
names(df4_c)[names(df4_c) == 'Diagnosis1'] <- 'Diagnosis'

df1_a$Diagnosis<- as.factor(df1_a$Diagnosis)
df1_b$Diagnosis<- as.factor(df1_b$Diagnosis)
df1_c$Diagnosis<- as.factor(df1_c$Diagnosis)
df1_d$Diagnosis<- as.factor(df1_d$Diagnosis)
df2_a$Diagnosis<- as.factor(df2_a$Diagnosis)
df2_b$Diagnosis<- as.factor(df2_b$Diagnosis)
df2_c$Diagnosis<- as.factor(df2_c$Diagnosis)
df2_d$Diagnosis<- as.factor(df2_d$Diagnosis)
df3_a$Diagnosis<- as.factor(df3_a$Diagnosis)
df3_b$Diagnosis<- as.factor(df3_b$Diagnosis)
df3_c$Diagnosis<- as.factor(df3_c$Diagnosis)
df3_d$Diagnosis<- as.factor(df3_d$Diagnosis)
df4_a$Diagnosis<- as.factor(df4_a$Diagnosis)
df4_b$Diagnosis<- as.factor(df4_b$Diagnosis)
df4_c$Diagnosis<- as.factor(df4_c$Diagnosis)
df4_d$Diagnosis<- as.factor(df4_d$Diagnosis)
names(df1_a) <- make.names(names(df1_a))
names(df1_b) <- make.names(names(df1_b))
names(df1_c) <- make.names(names(df1_c))
names(df1_d) <- make.names(names(df1_d))
names(df2_a) <- make.names(names(df2_a))
names(df2_b) <- make.names(names(df2_b))
names(df2_c) <- make.names(names(df2_c))
names(df2_d) <- make.names(names(df2_d))
names(df3_a) <- make.names(names(df3_a))
names(df3_b) <- make.names(names(df3_b))
names(df3_c) <- make.names(names(df3_c))
names(df3_d) <- make.names(names(df3_d))
names(df4_a) <- make.names(names(df4_a))
names(df4_b) <- make.names(names(df4_b))
names(df4_c) <- make.names(names(df4_c))
names(df4_d) <- make.names(names(df4_d)) 


nr <- 500
################################################################################
# Set A
############################
df_selected <- df1_a
names(df_selected) <- make.names(names(df_selected))
df_selected$Diagnosis <- as.factor(df_selected$Diagnosis)
samp_1_a <- sample(nrow(df_selected), 0.6 * nrow(df_selected))
test <- df_selected[samp_1_a,]
train <- df_selected[-samp_1_a,]
# dfeault number of trees 500, default mtry = 39
rf1_a=randomForest(Diagnosis ~ ., data = df_selected, ntree=nr)
df_selected <- df2_a
names(df_selected) <- make.names(names(df_selected))
df_selected$Diagnosis <- as.factor(df_selected$Diagnosis)
samp_2_a <- sample(nrow(df_selected), 0.6 * nrow(df_selected))
test <- df_selected[samp_2_a,]
train <- df_selected[-samp_2_a,]
# dfeault number of trees 500, default mtry = 39
rf2_a=randomForest(Diagnosis ~ ., data = df_selected, ntree=nr)
df_selected <- df3_a
names(df_selected) <- make.names(names(df_selected))
df_selected$Diagnosis <- as.factor(df_selected$Diagnosis)
samp_3_a <- sample(nrow(df_selected), 0.6 * nrow(df_selected))
test <- df_selected[samp_3_a,]
train <- df_selected[-samp_3_a,]
# dfeault number of trees 500, default mtry = 39
rf3_a=randomForest(Diagnosis ~ ., data = df_selected, ntree=nr)
df_selected <- df4_a
names(df_selected) <- make.names(names(df_selected))
df_selected$Diagnosis <- as.factor(df_selected$Diagnosis)
samp_4_a <- sample(nrow(df_selected), 0.6 * nrow(df_selected))
test <- df_selected[samp_4_a,]
train <- df_selected[-samp_4_a,]
# dfeault number of trees 500, default mtry = 39
rf4_a=randomForest(Diagnosis ~ ., data = df_selected, ntree=nr)

################################################################################
# Set B
############################
df_selected <- df1_b
names(df_selected) <- make.names(names(df_selected))
df_selected$Diagnosis <- as.factor(df_selected$Diagnosis)
samp_1_b <- sample(nrow(df_selected), 0.6 * nrow(df_selected))
test <- df_selected[samp_1_b,]
train <- df_selected[-samp_1_b,]
# dfeault number of trees 500, default mtry = 39
rf1_b=randomForest(Diagnosis ~ ., data = df_selected, ntree=nr)
df_selected <- df2_b
names(df_selected) <- make.names(names(df_selected))
df_selected$Diagnosis <- as.factor(df_selected$Diagnosis)
samp_2_b <- sample(nrow(df_selected), 0.6 * nrow(df_selected))
test <- df_selected[samp_2_b,]
train <- df_selected[-samp_2_b,]
# dfeault number of trees 500, default mtry = 39
rf2_b=randomForest(Diagnosis ~ ., data = df_selected, ntree=nr)
df_selected <- df3_b
names(df_selected) <- make.names(names(df_selected))
df_selected$Diagnosis <- as.factor(df_selected$Diagnosis)
samp_3_b <- sample(nrow(df_selected), 0.6 * nrow(df_selected))
test <- df_selected[samp_3_b,]
train <- df_selected[-samp_3_b,]
# dfeault number of trees 500, default mtry = 39
rf3_b=randomForest(Diagnosis ~ ., data = df_selected, ntree=nr)
df_selected <- df4_b
names(df_selected) <- make.names(names(df_selected))
df_selected$Diagnosis <- as.factor(df_selected$Diagnosis)
samp_4_b <- sample(nrow(df_selected), 0.6 * nrow(df_selected))
test <- df_selected[samp_4_b,]
train <- df_selected[-samp_4_b,]
# dfeault number of trees 500, default mtry = 39
rf4_b=randomForest(Diagnosis ~ ., data = df_selected, ntree=nr)

################################################################################
# Set C
############################
df_selected <- df1_c
names(df_selected) <- make.names(names(df_selected))
df_selected$Diagnosis <- as.factor(df_selected$Diagnosis)
samp_1_c <- sample(nrow(df_selected), 0.6 * nrow(df_selected))
test <- df_selected[samp_1_c,]
train <- df_selected[-samp_1_c,]
# dfeault number of trees 500, default mtry = 39
rf1_c=randomForest(Diagnosis ~ ., data = df_selected, ntree=nr)
df_selected <- df2_c
names(df_selected) <- make.names(names(df_selected))
df_selected$Diagnosis <- as.factor(df_selected$Diagnosis)
samp_2_c <- sample(nrow(df_selected), 0.6 * nrow(df_selected))
test <- df_selected[samp_2_c,]
train <- df_selected[-samp_2_c,]
# dfeault number of trees 500, default mtry = 39
rf2_c=randomForest(Diagnosis ~ ., data = df_selected, ntree=nr)
df_selected <- df3_c
names(df_selected) <- make.names(names(df_selected))
df_selected$Diagnosis <- as.factor(df_selected$Diagnosis)
samp_3_c <- sample(nrow(df_selected), 0.6 * nrow(df_selected))
test <- df_selected[samp_3_c,]
train <- df_selected[-samp_3_c,]
# dfeault number of trees 500, default mtry = 39
rf3_c=randomForest(Diagnosis ~ ., data = df_selected, ntree=nr)
df_selected <- df4_c
names(df_selected) <- make.names(names(df_selected))
df_selected$Diagnosis <- as.factor(df_selected$Diagnosis)
samp_4_c <- sample(nrow(df_selected), 0.6 * nrow(df_selected))
test <- df_selected[samp_4_c,]
train <- df_selected[-samp_4_c,]
# dfeault number of trees 500, default mtry = 39
rf4_c=randomForest(Diagnosis ~ ., data = df_selected, ntree=nr)

################################################################################
# Set D
############################
df_selected <- df1_d
names(df_selected) <- make.names(names(df_selected))
df_selected$Diagnosis <- as.factor(df_selected$Diagnosis)
samp_1_d <- sample(nrow(df_selected), 0.6 * nrow(df_selected))
test <- df_selected[samp_1_d,]
train <- df_selected[-samp_1_d,]
# dfeault number of trees 500, default mtry = 39
rf1_d=randomForest(Diagnosis ~ ., data = df_selected, ntree=nr)
df_selected <- df2_d
names(df_selected) <- make.names(names(df_selected))
df_selected$Diagnosis <- as.factor(df_selected$Diagnosis)
samp_2_d <- sample(nrow(df_selected), 0.6 * nrow(df_selected))
test <- df_selected[samp_2_d,]
train <- df_selected[-samp_2_d,]
# dfeault number of trees 500, default mtry = 39
rf2_d=randomForest(Diagnosis ~ ., data = df_selected, ntree=nr)
df_selected <- df3_d
names(df_selected) <- make.names(names(df_selected))
df_selected$Diagnosis <- as.factor(df_selected$Diagnosis)
samp_3_d <- sample(nrow(df_selected), 0.6 * nrow(df_selected))
test <- df_selected[samp_3_d,]
train <- df_selected[-samp_3_d,]
# dfeault number of trees 500, default mtry = 39
rf3_d=randomForest(Diagnosis ~ ., data = df_selected, ntree=nr)
df_selected <- df4_d
names(df_selected) <- make.names(names(df_selected))
df_selected$Diagnosis <- as.factor(df_selected$Diagnosis)
samp_4_d <- sample(nrow(df_selected), 0.6 * nrow(df_selected))
test <- df_selected[samp_4_d,]
train <- df_selected[-samp_4_d,]
# dfeault number of trees 500, default mtry = 39
rf4_d=randomForest(Diagnosis ~ ., data = df_selected, ntree=nr)

library(randomForest)
source("/hps/nobackup/ma/natalja/data/mQTL_data/scripts/load_data.R")

err_rates <- c(0)
for(i in 1:10){
  rf=randomForest(Diagnosis ~ ., data = df4_b)
  err_rates <- rbind(err_rates,rf$err.rate[,1])
}
err_rates <- err_rates[-1,]
mean_res4_b <- apply(err_rates,2,mean)

err_rates <- c(0)
for(i in 1:10){
  rf=randomForest(Diagnosis ~ ., data = df4_a)
  err_rates <- rbind(err_rates,rf$err.rate[,1])
}
err_rates <- err_rates[-1,]
mean_res4_a <- apply(err_rates,2,mean)

err_rates <- c(0)
for(i in 1:10){
  rf=randomForest(Diagnosis ~ ., data = df4_c)
  err_rates <- rbind(err_rates,rf$err.rate[,1])
}
err_rates <- err_rates[-1,]
mean_res4_c <- apply(err_rates,2,mean)

err_rates <- c(0)
for(i in 1:10){
  rf=randomForest(Diagnosis ~ ., data = df4_d)
  err_rates <- rbind(err_rates,rf$err.rate[,1])
}
err_rates <- err_rates[-1,]
mean_res4_d <- apply(err_rates,2,mean)


mean_B4 <- mean_res
mean_ADC <- mean_res
err_rates <- cbind(mean_ADC,mean_B4)

cols <- c("#f58231",
          "#000080",
          "#e6194b"
)

ltys <- c(1,1,1)
nr<-500
pdf("RF_ADC_B4_B.pdf")
matplot(1:nr , err_rates, col=cols,type="l",lty=ltys, ylab="Error",xlab="trees")
legend("topright",
       legend=c("Set A with 9 GWAS SNPs, Classifier IV","Set A with B4 SNPs, Classifier IV","Set B, Classifier IV")
       ,lty=ltys, 
       col=cols,cex=0.5)
dev.off()


rf1_a <- randomForest(Diagnosis ~ ., data = df1_a)
rf2_a <- randomForest(Diagnosis ~ ., data = df2_a)
rf3_a <- randomForest(Diagnosis ~ ., data = df3_a)

rf1_b <- randomForest(Diagnosis ~ ., data = df1_b)
rf2_b <- randomForest(Diagnosis ~ ., data = df2_b)
rf3_b <- randomForest(Diagnosis ~ ., data = df3_b)

rf1_c <- randomForest(Diagnosis ~ ., data = df1_c)
rf2_c <- randomForest(Diagnosis ~ ., data = df2_c)
rf3_c <- randomForest(Diagnosis ~ ., data = df3_c)

rf1_d <- randomForest(Diagnosis ~ ., data = df1_d)
rf2_d <- randomForest(Diagnosis ~ ., data = df2_d)
rf3_d <- randomForest(Diagnosis ~ ., data = df3_d)

rf4_a <- randomForest(Diagnosis ~ ., data = df4_a)
rf4_b <- randomForest(Diagnosis ~ ., data = df4_b)
rf4_c <- randomForest(Diagnosis ~ ., data = df4_c)
rf4_d <- randomForest(Diagnosis ~ ., data = df4_d)

err_rates <- cbind(rf1_a$err.rate[,1],rf2_a$err.rate[,1])
err_rates <- cbind(err_rates,rf3_a$err.rate[,1])
err_rates <- cbind(err_rates,mean_res4_a)
err_rates <- cbind(err_rates,rf1_b$err.rate[,1])
err_rates <- cbind(err_rates,rf2_b$err.rate[,1])
err_rates <- cbind(err_rates,rf3_b$err.rate[,1])
err_rates <- cbind(err_rates,mean_res4_b)
err_rates <- cbind(err_rates,rf1_c$err.rate[,1])
err_rates <- cbind(err_rates,rf2_c$err.rate[,1])
err_rates <- cbind(err_rates,rf3_c$err.rate[,1])
err_rates <- cbind(err_rates,mean_res4_c)
err_rates <- cbind(err_rates,rf1_d$err.rate[,1])
err_rates <- cbind(err_rates,rf2_d$err.rate[,1])
err_rates <- cbind(err_rates,rf3_d$err.rate[,1])
err_rates <- cbind(err_rates,mean_res4_d)

cols <- c("#f58231",
  "#f58231",
  "#f58231",
  "#f58231",
  "#000080",
  "#000080",
  "#000080",
  "#000080",
  "#808000",
  "#808000",
  "#808000",
  "#808000",
  "#e6194b",
  "#e6194b",
  "#e6194b",
  "#e6194b")

ltys <- c(2,3,4,1,2,3,4,1,2,3,4,1,2,3,4,1)

png("RF_all_sets.png")
matplot(1:nr , err_rates, col=cols,type="l",lty=ltys, ylab="OOB error",xlab="trees")
legend("topright",
       legend=c("Set A, Classifier I","Set A, Classifier II","Set A, Classifier III","Set A, Classifier IV",
                "Set B, Classifier I","Set B, Classifier II","Set B, Classifier III","Set B, Classifier IV",
                "Set C, Classifier I","Set C, Classifier II","Set C, Classifier III","Set C, Classifier IV",
                "Set D, Classifier I","Set D, Classifier II","Set D, Classifier III","Set D, Classifier IV")
       ,lty=ltys, 
       col=cols,cex=0.8)
dev.off()

min(err_rates)
#...
# Set B Diagnosis IV - smallest error = 0.03125
which(rf4_b$err.rate[,1]==min(rf4_b$err.rate[,1]))
#309



df_selected <- df4_a
names(df_selected) <- make.names(names(df_selected))
df_selected$Diagnosis <- as.factor(df_selected$Diagnosis)
test <- df_selected[samp_4_b,]
train <- df_selected[-samp_4_b,]

# default number of trees 500, default mtry = 39
rf=randomForest(Diagnosis ~ ., data = train)

which(Diagnosis.rf$err.rate[,1]==min(Diagnosis.rf$err.rate[,1]))
# 116 149 150

#Plotting the Error vs Number of Trees Graph.
pdf("RF_Errors_vs_NumberofTrees_B4.pdf")
plot(rf4_b, main="Set B, diagnosis IV")
legend("topright",
       legend=c("ADC class",
                "CTL class",
                "Out of Bag")
       ,lty=c(2,3,1), 
       col=c("red","green","black"))
dev.off()

nr <-400

oob.err=double(nr)
test.err_in=double(nr)
test.err_out=double(nr)

#mtry is no of Variables randomly chosen at each split
for(mtry in  1:nr) 
{
    rf=randomForest(Diagnosis ~ . , data = df_selected,mtry=mtry,ntree=309) 
    oob.err[mtry] = rf$err.rate[309,1] #Error of all Trees fitted
   
    pred<-predict(rf,train) #Predictions on Test Set for each Tree
    res <- table(pred,train$Diagnosis)
    test.err_in[mtry]= 1 - (res[1,1]+res[2,2])/sum(res)

  
    pred<-predict(rf,test) #Predictions on Test Set for each Tree
    res <- table(pred,test$Diagnosis)
    test.err_out[mtry]= 1 - (res[1,1]+res[2,2])/sum(res)
  
    cat(mtry," ") #printing the output to the console
}

res <- cbind(oob.err,test.err_in,test.err_out)
res <- res[which(rowSums(res) > 0),] 



which(oob.err==min(oob.err))
#Plotting the Error vs Number of Trees Graph.
pdf("RF_OutofBagError_B4.pdf")
matplot(1:nr , oob.err, type="l",lty=1, col=c("red"),ylab="Mean Squared Error",xlab="Number of Predictors Considered at each Split")
#legend("topright",legend=c("Out of Bag Error"),pch=19, col=c("red"))
dev.off()
#A4 mtry = 91, 93
#B4 mtry = 199, 322

# rf=randomForest(Diagnosis ~ . , data = df_selected,ntree=442, mtry=91)
#randomForest(formula = Diagnosis ~ ., data = df_selected, ntree = 500,      mtry = 91) 
#Type of random forest: classification
#Number of trees: 500
#No. of variables tried at each split: 91

#OOB estimate of  error rate: 2.37%
#Confusion matrix:
#  ADC CTL class.error
#ADC 160   2  0.01234568
#CTL   5 128  0.03759398



#Number of trees: 150 (ntree)
#No. of variables tried at each split: 159 (mtry)


rf=randomForest(Diagnosis ~ . , data = train ,mtry=159,ntree=150,importance=TRUE) 
# oob.err = 2.82%
pred <- predict(rf, newdata = test)

table(pred, test$Diagnosis)

importance    <- importance(rf,type=1,scale=FALSE)
varImportance <- data.frame(Variables = row.names(importance), MDA = round(importance[ ,'MeanDecreaseAccuracy'],5), MDG = round(importance[ ,'MeanDecreaseGini'],5))

write.csv(varImportance,"feature_selected_D4.txt", row.names = FALSE)

importance_cfa_conditional <- varimp(cfa,conditional=TRUE)

varImportance_cfa <- data.frame(Variables = names(importance_cfa), CF = round(as.numeric(importance_cfa),5), 
                                CF_COND=round(as.numeric(importance_cfa_conditional),5))
write.csv(varImportance_cfa,"feature_selected_D4_cf.txt", row.names = FALSE)

#Importance scores
df_selected <- df4_b
names(df_selected) <- make.names(names(df_selected))
df_selected$Diagnosis <- as.factor(df_selected$Diagnosis)
df_selected_b <- df_selected

df_selected <- df4_a
names(df_selected) <- make.names(names(df_selected))
df_selected$Diagnosis <- as.factor(df_selected$Diagnosis)
df_selected_a <- df_selected

# mtry and ntree - parameters tuning results
rf_a=randomForest(Diagnosis ~ . , data = df_selected_a ,mtry=91,ntree=443,importance=TRUE)
rf_b=randomForest(Diagnosis ~ . , data = df_selected_b ,mtry=199,ntree=309,importance=TRUE)

importance    <- importance(rf_b)
varImportance <- data.frame(Variables = row.names(importance), MDA_b = round(importance[ ,'MeanDecreaseAccuracy'],5), MDG_b = round(importance[ ,'MeanDecreaseGini'],5))
varImportance$MDA_a <- 0
varImportance$MDG_a <- 0
importance    <- importance(rf_a)
varImportance_a <- data.frame(Variables = row.names(importance), MDA_a = round(importance[ ,'MeanDecreaseAccuracy'],5), MDG_a = round(importance[ ,'MeanDecreaseGini'],5))

for(v in varImportance_a$Variables){
  varImportance[varImportance$Variables==v,]$MDA_a <- varImportance_a[varImportance_a$Variables==v,]$MDA_a
  varImportance[varImportance$Variables==v,]$MDG_a <- varImportance_a[varImportance_a$Variables==v,]$MDG_a
}

library(party)
cfa <- cforest(Diagnosis ~ . ,data=df_selected_a,control=cforest_unbiased(mtry=91,ntree=443))

library(party)
cfa <- cforest(Diagnosis ~ . ,data=df4_d,control=cforest_unbiased(mtry=70,ntree=680))

importance_cfa <- varimp(cfa)

varImportance$CF_a <- 0
varImportance$CFCOND_a <- 0
varImportance$CF_b <- 0
varImportance$CFCOND_b <- 0
importance_cfa_conditional <- varimp(cfa,conditional=TRUE)

varImportance_cfa <- data.frame(Variables = names(importance_cfa), CF = round(as.numeric(importance_cfa),5), 
                                CF_COND=round(as.numeric(importance_cfa_conditional),5))

for(v in names(importance_cfa)){
  varImportance[varImportance$Variables==v,]$CF_a <-importance_cfa[names(importance_cfa)==v,]
  varImportance[varImportance$Variables==v,]$MDG_a <- varImportance_a[varImportance_a$Variables==v,]$MDG_a
}

cfb <- cforest(Diagnosis ~ . ,data=df_selected_b,control=cforest_unbiased(mtry=199,ntree=309))
importance_cfb <- varimp(cfb)
importance_cfb_conditional <- varimp(cfb,conditional=TRUE)
varImportance_cfb <- data.frame(Variables = names(importance_cfb), CF = round(as.numeric(importance_cfb),5), 
                                CF_COND=round(as.numeric(importance_cfb_conditional),5))

for(v in names(importance_cfb)){
  varImportance[varImportance$Variables==v,]$CF_a <-varImportance_cfa[varImportance_cfa$Variables==v,]$CF
  varImportance[varImportance$Variables==v,]$CFCOND_a <- varImportance_cfa[varImportance_cfa$Variables==v,]$CF_COND
  varImportance[varImportance$Variables==v,]$CF_b <-varImportance_cfb[varImportance_cfb$Variables==v,]$CF
  varImportance[varImportance$Variables==v,]$CFCOND_b <- varImportance_cfb[varImportance_cfb$Variables==v,]$CF_COND
}
# Recursive Feature Elimination incorporating resampling
# https://topepo.github.io/caret/recursive-feature-elimination.html#rfe
library(caret)
# define the control using a random forest selection function for 10 fold cross-validation
control<-rfeControl(functions=rfFuncs, method="cv", number=10)
dim(df_selected)
# run the RFE algorithm
results <- rfe(df_selected[,1:1542],df_selected[,1543], sizes=c(1:1542),rfeControl=control)
# summarize the results
print(results)
# list the chosen features
predictors(results)
# plot the results
pdf("RF_RFE_A4.pdf")
plot(results, type=c("g", "o"))
dev.off()
write.csv(results,"feature_selected_RFE_A4.txt", row.names = FALSE)


df_selected <- df4_b
names(df_selected) <- make.names(names(df_selected))
df_selected$Diagnosis <- as.factor(df_selected$Diagnosis)
rf=randomForest(Diagnosis ~ . , data = df_selected) 
df_selected <- df1_a
names(df_selected) <- make.names(names(df_selected))
df_selected$Diagnosis <- as.factor(df_selected$Diagnosis)

test <- df_selected[df_selected$Diagnosis %in% c("sMCI","cMCI"),]
test$pred <- predict(rf,test)

table(test[,8475], test[,8476])

# B4 model, OOB 2.94% 
#      ADC CTL
#cMCI  11  12
#sMCI  12  68



# 48% of cMCI samples are classified as ADC
# 85% of sMCI samples are classified as CTL

df_selected <- df4_a
names(df_selected) <- make.names(names(df_selected))
df_selected$Diagnosis <- as.factor(df_selected$Diagnosis)
rf=randomForest(Diagnosis ~ . , data = df_selected) 

# A4 model, OOB 2.71%
#      ADC CTL
#cMCI  16   7
#sMCI  17  63

#     ADC CTL
#cMCI  26  10
#sMCI  34 106


# 72% of cMCI samples are classified as ADC
# 75% of sMCI samples are classified as CTL

library(Boruta)
df_selected <- df4_a
# df_selected <- df4_b
names(df_selected) <- make.names(names(df_selected))
df_selected$Diagnosis <- as.factor(df_selected$Diagnosis)
# for B4 - genotypes 
# to treat them as factor data type
convert <- c(1543:8474)
df_selected[,convert] <- data.frame(apply(df_selected[convert], 2, as.factor)) 


boruta.train <- Boruta(Diagnosis~., data = df_selected, doTrace = 2)

pdf("Boruta_B4.pdf")
plot(boruta.train, xlab = "", xaxt = "n")
dev.off()
# The tentative attributes will be classified as confirmed or rejected 
# by comparing the median Z score of the attributes with the median Z score of the best shadow attribute
final.boruta <- TentativeRoughFix(boruta.train)
print(final.boruta)
getSelectedAttributes(final.boruta, withTentative = F)

varImportance$Boruta1 <- 0
varImportance$Boruta2 <- 0
varImportance[varImportance$Variables %in% getSelectedAttributes(boruta.train, withTentative = F),]$Boruta1 <- 1
varImportance[varImportance$Variables %in% getSelectedAttributes(boruta.train, withTentative = T),]$Boruta2 <- 1

boruta.df <- attStats(final.boruta)
print(boruta.df)

write.csv(getSelectedAttributes(boruta.train, withTentative = T),"Boruta_feature_selected_A4.txt", row.names = FALSE)


library(mlr) # based on weka?
df_selected <- df4_a
names(df_selected) <- make.names(names(df_selected))
df_selected$Diagnosis <- as.factor(df_selected$Diagnosis)
df_task <- makeClassifTask(data = df_selected, target = "Diagnosis")
fv2 = generateFilterValuesData(df_task, method = c("information.gain", "chi.squared"))


pdf("FilterMethods_A4.pdf")
plotFilterValues(fv2)
dev.off()

library(caret)
# In order to remove redundant features the caret R package has been used. 
# It allows to analyze a correlation matrix of features (metabolites, genotypes 
# and covariates ???) and found the attributes that can be removed.  
# correlation of features except diagnoses 
drops <- c("Diagnosis","Diagnosis3","Diagnosis4")
df_selected <- df[,!(names(df) %in% drops)]
df_selected$Gender <- as.numeric(df_selected$Gender)
df_selected$Centre <- as.numeric(df_selected$Centre)

correlationMatrix <- cor(df_selected)
highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=0.75)
filteredMatrix <- df_selected[,-highlyCorrelated]

df_selected <- cbind(filteredMatrix,df[,8475])
colnames(df_selected)[ncol(df_selected)] <- "Diagnosis"

names(df_selected) <- make.names(names(df_selected))
df_selected$Diagnosis <- as.factor(df_selected$Diagnosis)
df_selected$Diagnosis2 <- as.factor(df_selected$Diagnosis2)





control <- trainControl(method="repeatedcv", number=10, repeats=3)
# train the model
model <- train(Diagnosis ~ . , data=df_selected, method="lvq", preProcess="scale", trControl=control)
# estimate variable importance
importance <- varImp(model, scale=FALSE)
# summarize importance
print(importance)
# plot importance
plot(importance)

control <- rfeControl(functions=rfFuncs, method="cv", number=10)
# Recursive Feature Elimination algorithm
# A Random Forest algorithm is used on each iteration to evaluate the model. 
# The algorithm is configured to explore all possible subsets of the attributes. 
results <- rfe(df_selected[,1:ncol(df_selected)], df_selected[,ncol(df_selected)], sizes=c(1:ncol(df_selected)), rfeControl=control)


library(randomForest)

df_selected <- df[,colnames(df) %in% selected_mets$V1 | colnames(df) %in% c("Diagnosis2","Diagnosis") ]

df_selected <- df

df_selected <- df[cvrt$Diagnosis %in% c("CTL","ADC"),]

names(df_selected) <- make.names(names(df_selected))
df_selected$Diagnosis <- as.factor(df_selected$Diagnosis)
df_selected$Diagnosis2 <- as.factor(df_selected$Diagnosis2)
df_selected$Gender <- as.factor(df_selected$Gender)
df_selected$Centre <- as.factor(df_selected$Centre)

samp <- sample(nrow(df_selected), 0.6 * nrow(df_selected))
test <- df_selected[samp,]
train <- df_selected[-samp,]



drops <- c("Gender","Centre")
df_selected <- df_selected[,!(names(df_selected) %in% drops)]

samp <- sample(nrow(df_selected), 0.6 * nrow(df_selected))
test <- df_selected[samp,]
train <- df_selected[-samp,]

model <- randomForest(Diagnosis2 ~ . - Diagnosis, data = train, importance=TRUE,ntree=500)
pred <- predict(model, newdata = test)

table(pred, test$Diagnosis2)
png("RandomForest.png")
varImpPlot(model)
dev.off()

model <- randomForest(Diagnosis2 ~ . - Diagnosis, data = train, importance=TRUE,ntree=2000)
importance    <- importance(model)
varImportance <- data.frame(Variables = row.names(importance), Importance = round(importance[ ,'MeanDecreaseGini'],2))

#Create a rank variable based on importance
library(dplyr)
rankImportance <- varImportance %>% mutate(Rank = paste0('#',dense_rank(desc(Importance))))

#Use ggplot2 to visualize the relative importance of variables
library(ggplot2)
library(ggthemes)
ggplot(rankImportance, aes(x = reorder(Variables, Importance), 
                           y = Importance, fill = Importance)) +
  geom_bar(stat='identity') + 
  geom_text(aes(x = Variables, y = 0.5, label = Rank),
            hjust=0, vjust=0.55, size = 4, colour = 'red') +
  labs(x = 'Variables') +
  coord_flip() + 
  theme_few()

library(party)
fit <- cforest(Diagnosis ~ . - Diagnosis2,data = train,controls=cforest_unbiased(ntree=2000, mtry=3))
fit2 <- cforest(Diagnosis2 ~ . - Diagnosis,data = train,controls=cforest_unbiased(ntree=2000, mtry=3))
pred <-  predict(fit, test, OOB=TRUE, type = "response")
table(pred,test$Diagnosis)
pred2 <-  predict(fit2, test, OOB=TRUE, type = "response")
table(pred2,test$Diagnosis2)

png("RandomForest2.png")
varImpPlot(fit)
dev.off()


#########
#clusters
#A fundamental question is how to determine the value of the parameter k. 
#If we looks at the percentage of variance explained as a function of the number of clusters: 
#One should choose a number of clusters so that adding another cluster doesn’t give much better modeling of the data. 
#More precisely, if one plots the percentage of variance explained by the clusters against the number of clusters, 
#the first clusters will add much information (explain a lot of variance), but at some point the marginal gain will drop, 
# giving an angle in the graph. The number of clusters is chosen at this point, hence the “elbow criterion”.
wssplot <- function(data, nc=500, seed=1234){
  wss <- (nrow(data)-1)*sum(apply(data,2,var))
  for (i in 2:nc){
    set.seed(seed)
    wss[i] <- sum(kmeans(data, centers=i)$withinss)}
  plot(1:nc, wss, type="b", xlab="Number of Clusters",
       ylab="Within groups sum of squares")}

png("WSSplot_NMR.png")
wssplot(expr_nmr)
dev.off()

expr_lcms <- rbind(expr1,expr2)
expr_lcms <- rbind(expr_lcms,expr3)

png("WSSplot_LCMS.png")
wssplot(expr_lcms)
dev.off()

cvrt = read.table("/hps/nobackup/ma/natalja/data/mQTL_data/covariates/Covariates_all_vertical_text2.txt",sep="\t",header=T,row.names=1)
cvrt <- cvrt[rownames(cvrt) %in% colnames(expr_nmr),]
cvrt <- cvrt[colnames(expr_nmr),]


cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}


mat_data <- data.matrix(expr_nmr)
rownames(mat_data) <- gsub("X\\.","",rownames(mat_data))
rownames(mat_data) <- gsub("X","",rownames(mat_data))
M <- cor(t(mat_data))
p.mat <- cor.mtest(t(mat_data))
# matrix of the p-value of the correlation

head(p.mat[, 1:5])



png("Corr_plot2.png")
corrplot(M, type="upper", order="hclust", 
         p.mat = p.mat, sig.level = 0.01, insig = "blank",addgrid.col = NA)

dev.off()


################################################################################
#First set the mtry to the default value (sqrt of total number of all predictors) 
#and search for the optimal ntree value. 
#To find the number of trees that correspond to a stable classifier, 
#we build random forest with different ntree values (100, 200, 300….,1,000). 
#We build 10 RF classifiers for each ntree value, record the OOB error rate and 
#see the number of trees where the out of bag error rate stabilizes and reach minimum.
library(caret)
library(party)
source("/hps/nobackup/ma/natalja/data/mQTL_data/scripts/load_data.R")
tuning_results <- data.frame(set="",nr=0, ntree=0, kappa=0,accuracy=0)
control <- trainControl(method="oob")
ntree_seq <- c(seq(0, 2500, by = 100)[-1])
sets <-c("df4_b") #c("df1_d","df2_d","df3_d","df4_d")
# c("df1_b","df2_b","df3_b","df4_b")
# "df1_c","df2_c","df3_c","df4_c",
# "df1_d","df2_d","df3_d","df4_d")
#"df1_a","df2_a","df3_a","df4_a",
for (set in sets){
  tuning_results <- data.frame(set="",nr=0, ntree=0, kappa=0,accuracy=0)
  for (ntree in ntree_seq) {
    data_set <- eval(parse(text = set))
    mtry_fixed <- sqrt(ncol(data_set))
    for(i in 1:10){
      mod <- train(Diagnosis ~ .,
                   data = data_set,
                   method = "cforest",
                   tuneGrid = data.frame(.mtry = mtry_fixed),
                   #trControl = control,
                   controls = cforest_unbiased(ntree = ntree))
      mod_res <- data.frame(set=set,nr=i,ntree=ntree,kappa=mod$results$Kappa, accuracy=mod$results$Accuracy)
      tuning_results <- rbind(tuning_results,mod_res)
    }

    print(ntree)
  }
  print(set)
  write.csv(tuning_results,paste(set,"_results.csv",sep=""),row.names = F)
} 

library(randomForest)
source("/hps/nobackup/ma/natalja/data/mQTL_data/scripts/load_data.R")
tuning_results <- data.frame(set="",nr=0, ntree=0, oob=0)
control <- trainControl(method="oob")
ntree_seq <- c(seq(0, 1000, by = 100)[-1])
sets <- #c("df1_d","df2_d","df3_d","df4_d")
 c("df4_d")
for (set in sets){
  tuning_results <- data.frame(set="",nr=0, ntree=0, oob=0)
  for (ntree in ntree_seq) {
    data_set <- eval(parse(text = set))
    mtry_fixed <- sqrt(ncol(data_set))
    for(i in 1:10){
      rf <- randomForest(Diagnosis ~ ., data = data_set, ntree=ntree)
      mod_res <- data.frame(set=set,nr=i,ntree=ntree,oob=rf$err.rate[nrow(rf$err.rate),1])
      tuning_results <- rbind(tuning_results,mod_res)
    }
    
    print(ntree)
  }
  print(set)
  write.csv(tuning_results,paste(set,"_results_rf.csv",sep=""),row.names = F)
} 


set <- "df4_d"
ntree_seq <- c(seq(200, 399, by = 10)[-1])
tuning_results <- data.frame(set="",nr=0, ntree=0, oob=0)
for (ntree in ntree_seq) {
  data_set <- eval(parse(text = set))
  mtry_fixed <- sqrt(ncol(data_set))
  for(i in 1:10){
    rf <- randomForest(Diagnosis ~ ., data = data_set, ntree=ntree)
    mod_res <- data.frame(set=set,nr=i,ntree=ntree,oob=rf$err.rate[nrow(rf$err.rate),1])
    tuning_results <- rbind(tuning_results,mod_res)
  }
  
  print(ntree)
}
print(set)
write.csv(tuning_results,paste(set,"_results_rf2.csv",sep=""),row.names = F)


set <- "df4_d"
ntree_seq <- c(seq(280, 299, by = 1)[-1])
tuning_results <- data.frame(set="",nr=0, ntree=0, oob=0)
for (ntree in ntree_seq) {
  data_set <- eval(parse(text = set))
  mtry_fixed <- sqrt(ncol(data_set))
  for(i in 1:10){
    rf <- randomForest(Diagnosis ~ ., data = data_set, ntree=ntree)
    mod_res <- data.frame(set=set,nr=i,ntree=ntree,oob=rf$err.rate[nrow(rf$err.rate),1])
    tuning_results <- rbind(tuning_results,mod_res)
  }
  
  print(ntree)
}
print(set)
write.csv(tuning_results,paste(set,"_results_rf3.csv",sep=""),row.names = F)


set <- "df4_d"
ntree_seq <- c(seq(310, 329, by = 1)[-1])
tuning_results <- data.frame(set="",nr=0, ntree=0, oob=0)
for (ntree in ntree_seq) {
  data_set <- eval(parse(text = set))
  mtry_fixed <- sqrt(ncol(data_set))
  for(i in 1:10){
    rf <- randomForest(Diagnosis ~ ., data = data_set, ntree=ntree)
    mod_res <- data.frame(set=set,nr=i,ntree=ntree,oob=rf$err.rate[nrow(rf$err.rate),1])
    tuning_results <- rbind(tuning_results,mod_res)
  }
  
  print(ntree)
}
print(set)
write.csv(tuning_results,paste(set,"_results_rf4.csv",sep=""),row.names = F)

set <- "df4_d"
data_set <- eval(parse(text = set))
mtry_fixed <- ncol(data_set)

mtry_seq <- c(seq(0, mtry_fixed , by = 10)[-1])
tuning_results <- data.frame(set="",nr=0, mtry=0, oob=0)
for (mtry in mtry_seq) {
  for(i in 1:10){
    rf <- randomForest(Diagnosis ~ ., data = data_set, ntree=680, mtry=mtry)
    mod_res <- data.frame(set=set,nr=i,mtry=mtry,oob=rf$err.rate[nrow(rf$err.rate),1])
    tuning_results <- rbind(tuning_results,mod_res)
  }
  
  print(mtry)
}
print(set)
write.csv(tuning_results,paste(set,"_results_rf_mtry.csv",sep=""),row.names = F)



#There are two ways to find the optimal mtry :
#Apply a similar procedure such that random forest is run 10 times. 
#The optimal number of predictors selected for split is selected for which out of bag error rate stabilizes and reach minimum.
#Experiment with including the (square root of total number of all predictors), 
#(half of this square root value), and (twice of the square root value). 
#And check which mtry returns maximum Area under curve. 
#Thus, for 1000 predictors the number of predictors to select for each node would be 16, 32, and 64 predictors.

###############################
#NTREE
res <- read.csv("df4_a_results_rf2.csv",header = T)
res <- res[-1,]
names <- unique(res$ntree)
# for randomForest
res <- matrix(res$oob,nrow=10)
# for cforest
#res <- matrix(res$kappa,nrow=10)

# Two or more results together
#res2 <- read.csv("df4_a_results_rf2.csv",header = T)
#res2 <- res2[-1,]
#res<- rbind(res,res2)
#res<- res[order(res$ntree,res$nr),]
#res_f <- aggregate(res$oob,by=list(ntree=res$ntree,nr=res$nr),mean)
#res_f<- res_f[order(res_f$ntree,res_f$nr),]
#names <- unique(res_f$ntree)
#res <- res_f
#colnames(res)<-c("ntree","nr","oob")
#res <- matrix(res$oob,nrow=10)

# SD and mean values for iterations
sd_res <- apply(res,2,sd)
mean_res <- apply(res,2,mean)
which(sd_res==min(sd_res))
which(mean_res==min(mean_res))
x_min <-names[which(mean_res==min(mean_res))]

# SD and mean values into dataframe
dat <- data.frame(ntree=names,OOB = mean_res, sd=sd_res)

#Plot version 1
plot(names,mean_res,pch=19,xlab="ntree",ylab="OOB error rate",main="Set A, classifier IV",
     ylim=c(min(mean_res-sd_res),max((mean_res+sd_res))))
lines(rbind(names,names,NA),rbind(mean_res-sd_res,mean_res+sd_res,NA))

#Plot version 2
library(ggplot2)
dat$highlight <- ifelse(dat$ntree %in% c(x_min), "red","black")
ticks <- data.frame (t = c(250, 500, 760, 1000))
png("ntree_tuning_4A.png")
ggplot(dat, aes(x = ntree, y = OOB)) +
  geom_line() +
  geom_ribbon(aes(ymin = OOB - sd,
                  ymax = OOB + sd), alpha = 0.2) +
  labs(y = "OOB error rate")+
  geom_point(aes(x=x_min,y=min(mean_res)), color="red") + 
  scale_x_continuous(breaks=c(ticks$t))
dev.off()

#MTRY
res <- read.csv("df4_a_results_rf_mtry.csv",header = T)
res <- res[-1,]
names <- unique(res$mtry)
res <- matrix(res$oob,nrow=10)
# SD and mean values for iterations
sd_res <- apply(res,2,sd)
mean_res <- apply(res,2,mean)
which(sd_res==min(sd_res))
x_min <-names[which(mean_res==min(mean_res))]
which(mean_res==min(mean_res))
# SD and mean values into dataframe
dat <- data.frame(mtry=names,OOB = mean_res, sd=sd_res)

#Plot version 1
plot(names,mean_res,pch=19,xlab="ntree",ylab="OOB error rate",main="Set A, classifier IV",
     ylim=c(min(mean_res-sd_res),max((mean_res+sd_res))))
lines(rbind(names,names,NA),rbind(mean_res-sd_res,mean_res+sd_res,NA))

#Plot version 2
library(ggplot2)
dat$highlight <- ifelse(dat$ntree %in% c(x_min), "red","black")
ticks <- data.frame (t = c(0, 70, 500,1000,1500))
png("mtry_tuning_4A.png")
ggplot(dat, aes(x = mtry, y = OOB)) +
  geom_line() +
  geom_ribbon(aes(ymin = OOB - sd,
                  ymax = OOB + sd), alpha = 0.2) +
  labs(y = "OOB error rate")+
  geom_point(aes(x=x_min,y=min(mean_res)), color="red") +
  scale_x_continuous(breaks=c(ticks$t))
#scale_x_continuous(breaks=seq(10, 1540, 100)) 
dev.off()