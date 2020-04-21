diagnosis <- ""
bases <- c("UHPOS","URPOS","URNEG")
bases <- c("SHPOS","SLPOS","SLNEG")
bases <- c("NMR_CPMG_Serum","NMR_Lipidomics_Negative_Serum","NMR_Lipidomics_Positive_Serum","NMR_Nosey_Serum","NMR_Lipidomics_Negative_Urine","NMR_Lipidomics_Positive_Urine","NMR_Nosey_Urine") 


base <- "UHPOS"
################################################################################
# Stats
################################################################################
for (base in c("NMR_CPMG_Serum","NMR_Lipidomics_Negative_Serum","NMR_Lipidomics_Positive_Serum","NMR_Nosey_Serum","NMR_Lipidomics_Negative_Urine","NMR_Lipidomics_Positive_Urine","NMR_Nosey_Urine") 
) {
  print(base)
  output_file_name <- paste("mqtl_",base,diagnosis,".txt",sep="") 
  permuted_output_file_name <- paste("mqtl_permuted_",base,diagnosis,".txt",sep="") 
  res <- read.table(output_file_name,header=T,sep="\t")
  res_permuted <- read.table(permuted_output_file_name,header=T,sep="\t")
  print(dim(res[res$p.value<0.01,]))
  print(dim(res[res$FDR<0.01,]))
  print(dim(res_permuted[res_permuted$p.value<0.01,]))
  print(dim(res_permuted[res_permuted$FDR<0.01,]))
}
################################################################################

################################################################################
# Historgrams and manhattan plots
################################################################################
source("/hps/nobackup/ma/natalja/data/mQTL_data/scripts/mqtl.R")
library(biomaRt)
library(qqman)
snp.db <- useMart("ENSEMBL_MART_SNP", dataset="hsapiens_snp")

for(base in c("SHPOS","SLPOS","SLNEG")){
  output_file_name <- paste("mqtl_",base,diagnosis,".txt",sep="") 
  res <- read.table(output_file_name,header=T,sep="\t")
  res$SNP <- gsub('_2','',res$SNP)
  res$SNP <- gsub('_1','',res$SNP)
  
  png(paste("histogram_p_values_",base,diagnosis,".png",sep=""))
  hist(res$p.value[!is.na(res$p.value)],breaks = 50,xlab = "P-values",col = "paleturquoise3",main = paste(base,diagnosis,sep=""))
  dev.off()
  
  res_split <- split(res, (seq(nrow(res))-1) %/% 20000) 
  res_all <- res[0,]
  for(i in 1:length(res_split)){
    res_part <- res_split[[i]]
    nt.biomart <- getBM(c("refsnp_id","chr_name","chrom_start"),
                        filters="snp_filter",
                        values=res_part$SNP,
                        mart=snp.db)
    res_all_annot <- merge(res_part,nt.biomart,by.x="SNP",by.y="refsnp_id",all.x=T,all.y=F)
    res_all_annot <- res_all_annot[-grep("CHR_", res_all_annot$chr_name),]
    res_all<- rbind(res_all,res_all_annot)
  }
  
  res_all <- res_all[!is.na(res_all$chr_name) & res_all$chr_name!="",]
  res_all <- res_all[!is.na(res_all$p.value) & res_all$p.value!="",]
  res_all <- res_all[!is.na(res_all$chrom_start) & res_all$chrom_start!="",]
  res_all <- res_all[!is.na(res_all$SNP) & res_all$SNP!="",]
  write.table(res_all,paste("mqtl_annot_",base,diagnosis,".txt",sep="")
              ,sep="\t",row.names = FALSE)
  
  # FDR - qqman manhattan plot 
  res_plot<-data.frame(SNP=res_all$SNP,CHR=as.character(res_all$chr_name),P=res_all$FDR,
                       BP=res_all$chrom_start,stringsAsFactors=FALSE)
  res_plot$CHR[res_plot$CHR == "X"] <- "23"
  res_plot$CHR <- as.numeric(res_plot$CHR)
  
  png(paste("manhattan_FDR_mqtl_",base,diagnosis,".png",sep=""), width=950, height=400)
  manhattan(res_plot, ylim=c(0,35), main = paste(base,diagnosis,sep=""),suggestiveline=-log10(1e-2),genomewideline=FALSE)
  dev.off()
  
  # raw p-values manhattan.plot 
  res_plot<-data.frame(SNP=res_all$SNP,CHR=as.character(res_all$chr_name),P=res_all$p.value,
                       BP=res_all$chrom_start,stringsAsFactors=FALSE)
  res_plot$CHR[res_plot$CHR == "X"] <- "23"
  res_plot$CHR <- as.numeric(res_plot$CHR)
  
  png(paste("manhattan_pvalues_mqtl_",base,diagnosis,".png",sep=""), width=950, height=400)
  manhattan.plot(res_plot$CHR,res_plot$BP,res_plot$P, ylim=c(0,35), main = paste(base,diagnosis,sep=""))
  dev.off()
}  


################################################################################
# Manhattan plots 2
################################################################################

source("/hps/nobackup/ma/natalja/data/mQTL_data/scripts/mqtl.R")
library(biomaRt)
library(qqman)
base <- "..."
output_file_name <- paste("mqtl_annot_",base,diagnosis,".txt",sep="") 
res_all <- read.table(output_file_name,header=T,sep="\t")
# FDR - qqman manhattan plot 
res_plot<-data.frame(SNP=res_all$SNP,CHR=as.character(res_all$chr_name),P=res_all$FDR,
                     BP=res_all$chrom_start,stringsAsFactors=FALSE)
res_plot$CHR[res_plot$CHR == "X"] <- "23"
res_plot$CHR <- as.numeric(res_plot$CHR)

png(paste("manhattan_FDR_mqtl_",base,diagnosis,".png",sep=""), width=950, height=400)
manhattan(res_plot, ylim=c(0,35), main = paste(base,diagnosis,sep=""),suggestiveline=-log10(1e-2),genomewideline=FALSE)
dev.off()

# raw p-values manhattan.plot 
res_plot<-data.frame(SNP=res_all$SNP,CHR=as.character(res_all$chr_name),P=res_all$p.value,
                     BP=res_all$chrom_start,stringsAsFactors=FALSE)
res_plot$CHR[res_plot$CHR == "X"] <- "23"
res_plot$CHR <- as.numeric(res_plot$CHR)

png(paste("manhattan_pvalues_mqtl_",base,diagnosis,".png",sep=""), width=950, height=400)
manhattan.plot(res_plot$CHR,res_plot$BP,res_plot$P, ylim=c(0,35), main = paste(base,diagnosis,sep=""))
dev.off()

################################################################################
# Previous studies 
################################################################################
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3446258/
ps <- read.table("/hps/nobackup/ma/natalja/data/mQTL_data/previous_studies.txt",sep="\t",header = T,stringsAsFactors = F)

for(base in c("UHPOS","URPOS","URNEG")){
  main_dir <- "/hps/nobackup/ma/natalja/data/mQTL_data/results/"
  
  output_file_name <- paste(main_dir,"result_",base,"_annot_full.txt",sep="") 
  res_all <- read.table(output_file_name,header=T,sep="\t")
  print(base)
  print(intersect(res_all$SNP, ps$SNP.ID))
  #write.table(res_all[res_all$SNP %in% intersect(res_all$SNP, ps$SNP.ID),],paste(main_dir,"previous_studies_",base,".txt",sep=""),sep="\t")
}

for(base in c("SHPOS","SLPOS","SLNEG")){
  main_dir <- "/hps/nobackup/ma/natalja/data/mQTL_data/results/normalization_final/"
  output_file_name <- paste(main_dir,"result_",base,"_annot_full.txt",sep="") 
  res_all <- read.table(output_file_name,header=T,sep="\t",stringsAsFactors = F)
  print(base)
  print(intersect(res_all$SNP, ps$SNP.ID))
  #write.table(res_all[res_all$SNP %in% intersect(res_all$SNP, ps$SNP.ID),],paste(main_dir,"previous_studies_",base,".txt",sep=""),sep="\t")

}
  
for (base in c("NMR_CPMG_Serum","NMR_Lipidomics_Negative_Serum","NMR_Lipidomics_Positive_Serum","NMR_Lipidomics_Negative_Urine","NMR_Lipidomics_Positive_Urine","NMR_Nosey_Urine"))
{
  main_dir <- "/hps/nobackup/ma/natalja/data/mQTL_data/NMR_Results/normalization_final/"
  output_file_name <- paste(main_dir,"result_",base,"_annot_full.txt",sep="")   
  res_all <- read.table(output_file_name,header=T,sep="\t",stringsAsFactors = F)
  print(base)
  print(intersect(res_all$SNP, ps$SNP.ID))
  
} 


################################################################################
# LC-MS and NMR intersection  
################################################################################

res_serum <-res_all[0,]
for(base in c("SHPOS","SLPOS","SLNEG")){
  main_dir <- "/hps/nobackup/ma/natalja/data/mQTL_data/results/EigenMS_QN_removed_cohort_centre/"
  output_file_name <- paste(main_dir,"result_",base,"_annot_full.txt",sep="")   
  res_all <- read.table(output_file_name,header=T,sep="\t",stringsAsFactors = F)
  res_serum <- rbind(res_serum,res_all)
}

res_serum_nmr <-res_all[0,]
for (base in c("NMR_CPMG_Serum","NMR_Lipidomics_Negative_Serum","NMR_Lipidomics_Positive_Serum"))
{
  main_dir <- "/hps/nobackup/ma/natalja/data/mQTL_data/NMR_Results/normalization_final/"
  output_file_name <- paste(main_dir,"result_",base,"_annot_full.txt",sep="") 
  res_all <- read.table(output_file_name,header=T,sep="\t",stringsAsFactors = F)
  res_serum_nmr <- rbind(res_serum_nmr,res_all)
} 

list_snps <- intersect(unique(res_serum$SNP),unique(res_serum_nmr$SNP))

unique(res_serum[res_serum$SNP %in% list_snps,]$gene)
unique(res_serum_nmr[res_serum_nmr$SNP %in% list_snps,]$gene)
# 94
# 68
res_urine <-res_all[0,]
for(base in c("UHPOS","URPOS","URNEG")){
  main_dir <- "/hps/nobackup/ma/natalja/data/mQTL_data/results/"
  output_file_name <- paste(main_dir,"result_",base,"_annot_full.txt",sep="")   
  res_all <- read.table(output_file_name,header=T,sep="\t",stringsAsFactors = F)
  res_urine <- rbind(res_urine,res_all)
}


res_urine_nmr <-res_all[0,]
for (base in c("NMR_Lipidomics_Negative_Urine","NMR_Lipidomics_Positive_Urine","NMR_Nosey_Urine"))
{
  main_dir <- "/hps/nobackup/ma/natalja/data/mQTL_data/NMR_Results/normalization_final/"
  output_file_name <- paste(main_dir,"result_",base,"_annot_full.txt",sep="") 
  res_all <- read.table(output_file_name,header=T,sep="\t",stringsAsFactors = F)
  res_urine_nmr <- rbind(res_urine_nmr,res_all)
} 
list_snps <- intersect(unique(res_urine$SNP),unique(res_urine_nmr$SNP))

unique(res_urine[res_urine$SNP %in% list_snps,]$gene)
unique(res_urine_nmr[res_urine_nmr$SNP %in% list_snps,]$gene)
# 379
# 261


################################################################################
# Linear Regression with diagnosis in a model 
################################################################################

threshold <- 0.01

res <- read.table(paste("/hps/nobackup/ma/natalja/data/mQTL_data/NMR_Results/mqtl_",base,".txt",sep=""),header=T,sep="\t")
res <- read.table(paste("mqtl_",base,".txt",sep=""),header=T,sep="\t")
res_sig <- res[res$FDR<threshold,]
res_sig$SNP <- gsub('_2','',res_sig$SNP)
res_sig$SNP <- gsub('_1','',res_sig$SNP)

dim(res_sig)
length(unique(res_sig$SNP))
length(unique(res_sig$gene))

write.table(unique(res_sig$SNP),paste(base,"_sig_SNP.txt",sep=""),sep="\t",row.names=F,quote=F)
write.table(unique(res_sig$gene),paste(base,"_sig_met.txt",sep=""),sep="\t",row.names=F,quote=F)

# grep -wFf sig_SNP.txt ../dna_matrices/dna_matrix_NMR_Lipidomics_Positive_Urine.tsv > dna_sig.txt
# grep -wFf sig_met.txt ../met_matrices/NMR/NMR_Lipidomics_Positive_Urine.tsv > met_sig.txt
# head -1  ../met_matrices/NMR/NMR_Lipidomics_Positive_Urine.tsv > header.txt
# cat header.txt met_sig.txt > met_sig2.txt 
# cat header.txt dna_sig.txt > dna_sig2.txt
expr <- read.table("met_sig2.txt",sep="\t",header=T,row.names=1,stringsAsFactors = F)
snps <- read.table("dna_sig2.txt",sep="\t",header=T,row.names=1,stringsAsFactors = F)
cvrt = read.table("/hps/nobackup/ma/natalja/data/mQTL_data/covariates/Covariates_NMR_Lipidomics_Positive_Urine.txt",sep="\t",header=T,row.names=1)



cov <-  read.table("/hps/nobackup/ma/natalja/data/mQTL_data/pheno.txt",sep="\t",header=T)
met_qn <-read.table(paste("/hps/nobackup/ma/natalja/data/mQTL_data/met_matrices/NMR/",base,".tsv",sep=""),
                    header=T,sep="\t")
res <- read.table(paste("/hps/nobackup/ma/natalja/data/mQTL_data/NMR_Results/result_",base,"_annot_full.txt",sep=""),
                  header=T,sep="\t",stringsAsFactors = F)

res <- read.table(paste("result_",base,"_annot_full.txt",sep=""),
                  header=T,sep="\t",stringsAsFactors = F)

rownames(cov) <- cov$FID
cov <-cov[!is.na(cov$Diagnosis),]
cov <- cov[cov$Diagnosis %in% c("ADC","CTL"),]

rownames(met_qn) <- met_qn$id
met_qn <- met_qn[,-1]

e <- met_qn[rownames(met_qn) %in% res$gene,names(met_qn) %in% rownames(cov)]
cov <- cov[rownames(cov) %in% colnames(e),]
cov <- cov[colnames(e),]

e <- t(e)
e <- as.data.frame(e)

w_t <- function(x){
  df <- structure(list(Diagnosis = cov$Diagnosis, 
                       Value = x), .Names = c("Diagnosis","Value"), class = "data.frame"); 
  t <- wilcox.test(formula = Value ~ Diagnosis, data = df);
  return(t$p.value)}

res <- apply(e, 2, w_t)
res_FDR <- p.adjust(res,"fdr")
e_imp <- e[,res_FDR<0.01] # important metabolites 

write.table(colnames(e_imp),paste("/hps/nobackup/ma/natalja/data/mQTL_data/NMR_Results/diagnosis_different_metabolites_",base,".txt",sep=""),sep="\t")

e$diagnosis <- cov$Diagnosis
write.table(e,paste("/home/natalja/Documents/EMIF/KCL/mQTL/mQTL_data/results/metabolites_",base,".txt",sep=""),sep="\t")