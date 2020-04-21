################################################################################
# Linear Regression with diagnosis in a model 
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
common_snps <- read.table("urine_common_snps.txt",sep="\t",header=TRUE,stringsAsFactors = FALSE)
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
     
     if (i$SNP %in% common_snps$x){
       lm3 = lm(e2~s2)
       png(paste("./plots_individual_associations_with_diagnosis_URPOS/snp_",as.character(i$SNP),"_",gsub('/','',as.character(i$gene)),".png",sep=""), width=950, height=400)
         plot(e2 ~ jitter(s2),
              col=(s2+1),xaxt="n",xlab="Genotype",ylab="Expression")
         axis(1,at=c(0:2),labels=c("AA","Aa","aa"))
         lines(lm3$fitted ~ s2,type="b",pch=15,col="darkgrey")
       dev.off()
     }
  }
}


result$SNP <- gsub('_2','',result$SNP)
result$SNP <- gsub('_1','',result$SNP)
length(unique(result$SNP))
length(unique(result[result$SNP %in% common_snps$x,]$SNP))

length(unique(result$gene))
length(unique(result[result$SNP %in% common_snps$x,]$gene))

result2 <- result[result$SNP %in% common_snps$x,]
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

library(biomaRt)
snp.db <- useMart("ENSEMBL_MART_SNP", dataset="hsapiens_snp")
nt.biomart <- getBM(c("refsnp_id","chr_name","ensembl_gene_stable_id"),
                    filters="snp_filter",
                    values=result2$SNP,
                    mart=snp.db)
res_annot <- merge(result,nt.biomart,by.x="SNP",by.y="refsnp_id",all.x=T,all.y=F)
res_annot <- res_annot[-grep("CHR_", res_annot$chr_name),]
genes <- unique(nt.biomart$ensembl_gene_stable_id)
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
genes_annot <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol','chromosome_name','description'), 
                     filters = 'ensembl_gene_id', values = genes, mart = ensembl)


res_annot <- merge(res_annot,genes_annot,by.x="ensembl_gene_stable_id",by.y="ensembl_gene_id",all.x=T,all.y=F)

write.table(res_annot,paste("diagnosis_result_annot_common_",base,".txt",sep=""),sep="\t")

write.table(result,paste("diagnosis_result_",base,".txt",sep=""),sep="\t")

###############################################
# ALMS1 and NAT2 SNPs  deeper investigation - heatmaps showing metabolite levels
###############################################
library(colorspace)
threshold <- 0.01

ALMS_snps <- read.table("/hps/nobackup/ma/natalja/data/mQTL_data/results/normalization_final/list_ALMS1.txt",sep="\t", header=FALSE, stringsAsFactors=F)
NAT2_snps <- read.table("/hps/nobackup/ma/natalja/data/mQTL_data/results/normalization_final/list_NAT2.txt",sep="\t", header=FALSE, stringsAsFactors=F)
ABO_snps <- read.table("/hps/nobackup/ma/natalja/data/mQTL_data/results/normalization_final/list_ABO.txt",sep="\t", header=FALSE, stringsAsFactors=F)
ALL_snps <- read.table("/hps/nobackup/ma/natalja/data/mQTL_data/results/normalization_final/list_ALL.txt",sep="\t", header=FALSE, stringsAsFactors=F)
selected_mets <- read.table("/hps/nobackup/ma/natalja/data/mQTL_data/results/normalization_final/list_selected_attributes.txt",sep="\t", header=FALSE, stringsAsFactors=F)
selected_mets2 <- read.table("/hps/nobackup/ma/natalja/data/mQTL_data/results/normalization_final/list_selected_attributes2.txt",sep="\t", header=FALSE, stringsAsFactors=F)


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

# The common samples for NMR and LC-MS
expr_nmr <- expr_nmr[,colnames(expr_nmr) %in% colnames(expr_urine_all),]
expr_urine_all <- expr_urine_all[,colnames(expr_urine_all) %in% colnames(expr_nmr),]
expr_nmr <- expr_nmr[,colnames(expr_urine_all)]

write.table(unique(res_sig_lcms$SNP),"/hps/nobackup/ma/natalja/data/mQTL_data/final_results_urine/selected_matrices/snps_lcms.txt",sep="\t",row.names=F,quote=F)
write.table(unique(res_sig_nmr$SNP),"/hps/nobackup/ma/natalja/data/mQTL_data/final_results_urine/selected_matrices/snps_nmr.txt",sep="\t",row.names=F,quote=F)

# ALMS1
res_sig_lcms_ALMS1 <- res_sig_lcms[res_sig_lcms$SNP %in% ALMS_snps$V1,]
expr_ALMS1_lcms <- expr_urine_all[rownames(expr_urine_all) %in% unique(res_sig_lcms_ALMS1$gene),]
res_sig_nmr_ALMS1 <- res_sig_nmr[res_sig_nmr$SNP %in% ALMS_snps$V1,]
expr_ALMS1_nmr <- expr_nmr[rownames(expr_nmr) %in% unique(res_sig_nmr_ALMS1$gene),]

# NAT2
res_sig_lcms_NAT2 <- res_sig_lcms[res_sig_lcms$SNP %in% NAT2_snps$V1,]
expr_NAT2_lcms <- expr_urine_all[rownames(expr_urine_all) %in% unique(res_sig_lcms_NAT2$gene),]
res_sig_nmr_NAT2 <- res_sig_nmr[res_sig_nmr$SNP %in% NAT2_snps$V1,]
expr_NAT2_nmr <- expr_nmr[rownames(expr_nmr) %in% unique(res_sig_nmr_NAT2$gene),]

# ABO
res_sig_lcms_ABO <- res_sig_lcms[res_sig_lcms$SNP %in% ABO_snps$V1,]
expr_ABO_lcms <- expr_urine_all[rownames(expr_urine_all) %in% unique(res_sig_lcms_ABO$gene),]
res_sig_nmr_ABO <- res_sig_nmr[res_sig_nmr$SNP %in% ABO_snps$V1,]
expr_ABO_nmr <- expr_nmr[rownames(expr_nmr) %in% unique(res_sig_nmr_ABO$gene),]


# SNPs
base <- "LCMS"
snps <- read.table(paste("/hps/nobackup/ma/natalja/data/mQTL_data/results/normalization_final/dna_",base,"_sig2.txt",sep=""),sep="\t",header=T,row.names=1,stringsAsFactors = F)
rownames(snps) <- gsub('_2','',rownames(snps))
rownames(snps) <- gsub('_1','',rownames(snps))

snps_ALMS1 <- snps[rownames(snps) %in% res_sig_lcms_ALMS1$SNP | rownames(snps) %in% res_sig_nmr_ALMS1$SNP,]
snps_ALMS1  <- snps_ALMS1 [,colnames(snps_ALMS1) %in% colnames(expr_ALMS1_nmr),]

snps_ABO <- snps[rownames(snps) %in% res_sig_lcms_ABO$SNP | rownames(snps) %in% res_sig_nmr_ABO$SNP,]
snps_ABO  <- snps_ABO [,colnames(snps_ABO) %in% colnames(expr_ABO_nmr),]

cvrt = read.table(paste(cov_matrices_dir,"Covariates_all_vertical_text2.txt",sep=""),sep="\t",header=T,row.names=1)
cvrt <- cvrt[rownames(cvrt) %in% colnames(expr_ALMS1_nmr),]
cvrt <- cvrt[colnames(expr_ALMS1_nmr),]

#################
## SNPs heatmap - ALMS1 genotypes patterns
# common samples betweeb NMR and LC-MS
# save test$colInd from clusering 
#################
expr_ALMS1 <- expr[rownames(expr) %in% unique(res_sig$gene),]
rownames(snps) <- gsub('_2','',rownames(snps))
rownames(snps) <- gsub('_1','',rownames(snps))
snps_ALMS1 <- snps[rownames(snps) %in% unique(res_sig$SNP),]

expr_nmr_ALMS1 <- expr_nmr[rownames(expr_nmr) %in% unique(res_sig$gene),]
rownames(snps) <- gsub('_2','',rownames(snps))
rownames(snps) <- gsub('_1','',rownames(snps))
snps_ALMS1 <- snps[rownames(snps) %in% unique(res_sig$SNP),]

colors_genotypes <- structure(c("grey36","red4","forestgreen"), names = c(0, 1, 2))
library(RColorBrewer)

png("ALMS1_SNPs.png")
test1 <- heatmap.2(as.matrix(snps_ALMS1),dendrogram="none", trace="none", scale="none", key=F, main="ALMS1 SNPs",col=colors_genotypes)

cvrt <- cvrt[test1$colInd,]

test <- heatmap.2(as.matrix(snps_ALMS1),dendrogram="none", trace="none", scale="none", key=F, main="ALMS1 SNPs",col=colors_genotypes,
                  labCol = c(""),ColSideColors=rainbow_hcl(4)[c(cvrt$Diagnosis)])

par(lend = 1)           # square line ends for the color legend
legend("topleft",      # location of the legend on the heatmap plot
       legend = c("AA", "Aa", "aa"), # category labels
       col = c("grey36","red4","forestgreen"),  # color key
       lty= 1,             # line style
       lwd = 10            # line width
)
dev.off()


cvrt = read.table(paste(cov_matrices_dir,"Covariates_all_vertical_text2.txt",sep=""),sep="\t",header=T,row.names=1)
cvrt <- cvrt[rownames(cvrt) %in% colnames(expr_ABO_lcms),]
cvrt <- cvrt[colnames(expr_ABO_lcms),]


png("ABO_SNPs.png")
test1 <- heatmap.2(as.matrix(snps_ABO),dendrogram="none", trace="none", scale="none", key=F, main="ABO SNPs",col=colors_genotypes)

cvrt <- cvrt[test1$colInd,]

test <- heatmap.2(as.matrix(snps_ABO),dendrogram="none", trace="none", scale="none", key=F, main="ABO SNPs",col=colors_genotypes,
                  labCol = c(""),ColSideColors=rainbow_hcl(4)[cvrt$Diagnosis])

par(lend = 1)           # square line ends for the color legend
legend("topleft",      # location of the legend on the heatmap plot
       legend = c("AA", "Aa", "aa"), # category labels
       col = c("grey36","red4","forestgreen"),  # color key
       lty= 1,             # line style
       lwd = 10            # line width
)
legend("topright", pch=16, cex=0.7, pt.cex = 1, col=rainbow_hcl(4),
      legend=unique(cvrt$Diagnosis))
dev.off()

  

png("ABO_LCMS_met.png")
heatmap.2(as.matrix(expr_ABO_lcms[,test1$colInd]),col=brewer.pal(11,"RdBu"),Colv=FALSE, dendrogram="row", trace="none")
dev.off()

png("ABO_NMR_met.png")
heatmap.2(as.matrix(expr_ABO_nmr[,test1$colInd]),col=brewer.pal(11,"RdBu"),Colv=FALSE, dendrogram="row", trace="none")
dev.off()


png("ABO_test.png")
t <- heatmap.2(as.matrix(expr_ABO_lcms),col=brewer.pal(11,"RdBu"), trace="none",colCol=c("firebrick3","plum4","springgreen4","plum4")[cvrt$Diagnosis])
cvrt <- cvrt[t$colInd,]
heatmap.2(as.matrix(expr_ABO_lcms),col=brewer.pal(11,"RdBu"), trace="none",ColSideColors=rainbow_hcl(4)[cvrt$Diagnosis])
dev.off()

png("ABO_test2.png")
t <- heatmap.2(as.matrix(expr_ABO_lcms[order(cvrt$Diagnosis)]),col=brewer.pal(11,"RdBu"), trace="none",Colv="none",colCol=c("firebrick3","plum4","springgreen4","plum4")[cvrt[order(cvrt$Diagnosis),]$Diagnosis])
dev.off()

png("ALMS1_LCMS_met.png")
heatmap.2(as.matrix(expr_ALMS1_lcms[,test$colInd]),col=brewer.pal(11,"RdBu"),Colv=FALSE, dendrogram="row", trace="none")
dev.off()
png("ALMS1_NMR_met.png")
heatmap.2(as.matrix(expr_ALMS1_nmr[,test$colInd]),col=brewer.pal(11,"RdBu"),Colv=FALSE, dendrogram="row", trace="none")
dev.off()

png("NAT2_LCMS_met.png")
heatmap.2(as.matrix(expr_NAT2_lcms[,test$colInd]),col=brewer.pal(11,"RdBu"),Colv=FALSE, dendrogram="row", trace="none")
dev.off()
png("NAT2_NMR_met.png")
heatmap.2(as.matrix(expr_NAT2_nmr[,test$colInd]),col=brewer.pal(11,"RdBu"),Colv=FALSE, dendrogram="row", trace="none")
dev.off()


#################### ALL SIGNIFICANT METABOLITES################################
drugs <- read.table("/hps/nobackup/ma/natalja/data/mQTL_data/drugs.txt",sep="\t",header = T,stringsAsFactors = F)
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
write.table(expr_nmr,"/hps/nobackup/ma/natalja/data/mQTL_data/final_results_urine/selected_matrices/met_nmr.txt",sep="\t",row.names=F,quote=F)
write.table(expr_lcms,"/hps/nobackup/ma/natalja/data/mQTL_data/final_results_urine/selected_matrices/met_lcms.txt",sep="\t",row.names=F,quote=F)
write.table(res_sig_lcms,"/hps/nobackup/ma/natalja/data/mQTL_data/final_results_urine/selected_matrices/res_sig_lcms.txt",row.names=F,quote=F)
write.table(res_sig_nmr,"/hps/nobackup/ma/natalja/data/mQTL_data/final_results_urine/selected_matrices/res_sig_nmr.txt",row.names=F,quote=F)
write.table(unique(res_sig_lcms$SNP),"/hps/nobackup/ma/natalja/data/mQTL_data/final_results_urine/selected_matrices/snps_lcms.txt",sep="\t",row.names=F,quote=F)
write.table(unique(res_sig_nmr$SNP),"/hps/nobackup/ma/natalja/data/mQTL_data/final_results_urine/selected_matrices/snps_nmr.txt",sep="\t",row.names=F,quote=F)


# SVM ????
library("e1071")
cvrt = read.table(paste(cov_matrices_dir,"Covariates_all_vertical_text2.txt",sep=""),sep="\t",header=T,row.names=1)
cvrt <- cvrt[rownames(cvrt) %in% colnames(expr_lcms),]
cvrt <- cvrt[colnames(expr_lcms),]
df <- as.data.frame(t(expr_lcms))
df$Diagnosis <- cvrt$Diagnosis

#Randomly shuffle the data
yourData<-yourData[sample(nrow(yourData)),]

#Create 10 equally size folds
folds <- cut(seq(1,nrow(yourData)),breaks=10,labels=FALSE)

#Perform 10 fold cross validation
for(i in 1:10){
  #Segement your data by fold using the which() function 
  testIndexes <- which(folds==i,arr.ind=TRUE)
  testData <- yourData[testIndexes, ]
  trainData <- yourData[-testIndexes, ]
  #Use the test and train data partitions however you desire...
}

attach(df)
x <- subset(df, select=-Diagnosis)
y <- Diagnosis
svm_model <- svm(Diagnosis ~ ., data=df)
summary(svm_model)
# or
svm_model1 <- svm(x,y)
summary(svm_model1)
pred <- predict(svm_model1,x)
table(pred,y)
svm_tune <- tune(svm, train.x=x, train.y=y, 
                 kernel="radial", ranges=list(cost=10^(-1:2), gamma=c(.5,1,2)))

tuned = tune.svm(churn~., data = trainset, gamma = 10^-2, cost = 10^2, tunecontrol=tune.control(cross=10))

print(svm_tune)
svm_model_after_tune <- svm(Diagnosis ~ ., data=df, kernel="radial", cost=0.1, gamma=0.5)
summary(svm_model_after_tune)
pred <- predict(svm_model_after_tune,x)

write.table(expr_lcms,paste("ML_data_lcms.txt",sep=""),sep="\t")

cvrt = read.table("/hps/nobackup/ma/natalja/data/mQTL_data/covariates/Covariates_all_vertical_text2.txt",sep="\t",header=T,row.names=1)
cvrt <- cvrt[rownames(cvrt) %in% colnames(expr_lcms),]
cvrt <- cvrt[colnames(expr_lcms),]
expr_lcms_ADC_CTL <- expr_lcms[,cvrt$Diagnosis %in% c("CTL","ADC")]
cvrt_ADC_CTL <- cvrt[colnames(expr_lcms_ADC_CTL),]
df <- as.data.frame(t(expr_lcms_ADC_CTL ))
df$Diagnosis <- as.character(cvrt_ADC_CTL$Diagnosis)
write.csv(df,"ML_data_lcms_ADC_CTL.csv",row.names = FALSE)

cvrt = read.table("/hps/nobackup/ma/natalja/data/mQTL_data/covariates/Covariates_all_vertical_text2.txt",sep="\t",header=T,row.names=1)
cvrt <- cvrt[rownames(cvrt) %in% colnames(expr_nmr),]
cvrt <- cvrt[colnames(expr_nmr),]
df <- as.data.frame(t(expr_nmr))
df$Diagnosis <- as.character(cvrt$Diagnosis)
write.csv(df,"ML_data_nmr.csv",row.names = FALSE)

expr_selected <- expr_lcms[rownames(expr_lcms) %in% selected_mets$V1,]
cvrt = read.table("/hps/nobackup/ma/natalja/data/mQTL_data/covariates/Covariates_all_vertical_text2.txt",sep="\t",header=T,row.names=1)
cvrt <- cvrt[rownames(cvrt) %in% colnames(expr_selected),]
cvrt <- cvrt[colnames(expr_selected),]
df <- as.data.frame(t(expr_selected))
df$Diagnosis <- as.character(cvrt$Diagnosis)

df_ADC <- df[df$Diagnosis %in% c("ADC"),]
df_ADC <- df_ADC[,-46]
res <- data.frame(values = c(as.matrix(df_ADC)), metabolites = rep(colnames(df_ADC),each=nrow(df_ADC)),Diagnosis = "ADC")
df_CTL <- df[df$Diagnosis %in% c("CTL"),]
df_CTL <- df_CTL[,-46]
res_CTL <- data.frame(values = c(as.matrix(df_CTL)), metabolites = rep(colnames(df_CTL),each=nrow(df_CTL)),Diagnosis = "CTL")
df_sMCI <- df[df$Diagnosis %in% c("sMCI"),]
df_sMCI <- df_sMCI[,-46]
res_sMCI <- data.frame(values = c(as.matrix(df_sMCI)), metabolites = rep(colnames(df_sMCI),each=nrow(df_sMCI)),Diagnosis = "sMCI")
df_cMCI <- df[df$Diagnosis %in% c("cMCI"),]
df_cMCI <- df_cMCI[,-46]
res_cMCI <- data.frame(values = c(as.matrix(df_cMCI)), metabolites = rep(colnames(df_cMCI),each=nrow(df_cMCI)),Diagnosis = "cMCI")

res <- rbind(res,res_CTL)
res <- rbind(res,res_sMCI)
res <- rbind(res,res_cMCI)
pdf("selected_metabolites_boxplot2.pdf")
p10 <- ggplot(res,aes(x=metabolites,y=values,fill=Diagnosis)) + geom_boxplot(alpha=0.7,outlier.size = 0, coef = 0) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p10
dev.off()


dfm <- aggregate(x = as.data.frame(t(expr_selected)), by = list(df$Diagnosis), FUN = "mean")
rownames(dfm) <- dfm$Group.1
dfm <- dfm [,-1]
fit <- kmeans(as.data.frame(expr_selected), 5)
dfm_t <- as.data.frame(t(dfm))
dfm_t <- dfm_t[order(fit$cluster),]

cols <- c("firebrick4","firebrick3","darkgreen","steelblue4")
ltys <- c(1,3,1,3) 
pdf("mean_values_selected_metabolites.pdf")
plot(dfm_t$ADC, type="l", col=cols[1],ylim=range(dfm_t), axes=F, ann=T,xlab="", ylab="Mean values",cex.lab=0.8, lwd=2)
grid(lty = 6, col = "cornsilk2") 
axis(1, at=seq(1,nrow(dfm_t),by=1),labels=rownames(dfm_t), las = 2,cex.axis=0.6)
# Plot y axis with smaller horizontal labels 
axis(2, las=1, cex.axis=0.8)
box()
lines(dfm_t$CTL, type="l", lwd=2, col=cols[3])

lines(dfm_t$sMCI, type="l", lty=3, lwd=2, col=cols[4])
lines(dfm_t$cMCI, type="l", lty=3, lwd=2, col=cols[2])

legend("topright", names(dfm_t), cex=0.8, col=cols, 
       lty=ltys, lwd=2, bty="n");
dev.off()

png("dfm.png")
plot(dfm_t$ADC, type="l", col="red",ylim=range(dfm_t), axes=F, ann=T,xlab="", ylab="Mean values",cex.lab=0.8, lwd=2)
axis(1, at=seq(1,nrow(dfm_t),by=1),labels=rownames(dfm_t), las = 2,cex.axis=0.6)
dev.off()
#1) UNKNOWN
# Remove known previously published SNP associations with metabolites
res_sig_lcms_unknown <- res_sig_lcms[!res_sig_lcms$SNP %in% ALL_snps$V1,]
expr_lcms_unknown <- expr_lcms[rownames(expr_lcms) %in% unique(res_sig_lcms_unknown$gene),]
res_sig_nmr_unknown <- res_sig_nmr[!res_sig_nmr$SNP %in% ALL_snps$V1,]
expr_nmr_unknown <- expr_nmr[rownames(expr_nmr) %in% unique(res_sig_nmr_ALMS1$gene),]

png("UNKNOWN_LCMS.png")
cvrt = read.table(paste(cov_matrices_dir,"Covariates_all_vertical_text2.txt",sep=""),sep="\t",header=T,row.names=1)
cvrt <- cvrt[rownames(cvrt) %in% colnames(expr_lcms_unknown),]
cvrt <- cvrt[colnames(expr_lcms_unknown),]

t <- heatmap.2(as.matrix(expr_lcms_unknown[order(cvrt$Diagnosis)]),col=brewer.pal(11,"RdBu"), 
               trace="none",Colv="none", dendrogram = "row",
               colCol=c("firebrick3","plum4","springgreen4","plum1")[cvrt[order(cvrt$Diagnosis),]$Diagnosis])
dev.off()


png("UNKNOWN_NMR.png")
cvrt = read.table(paste(cov_matrices_dir,"Covariates_all_vertical_text2.txt",sep=""),sep="\t",header=T,row.names=1)
cvrt <- cvrt[rownames(cvrt) %in% colnames(expr_nmr_unknown),]
cvrt <- cvrt[colnames(expr_nmr_unknown),]

t <- heatmap.2(as.matrix(expr_nmr_unknown[order(cvrt$Diagnosis)]),col=brewer.pal(11,"RdBu"), 
               trace="none",Colv="none", dendrogram = "row",
               colCol=c("firebrick3","plum4","springgreen4","plum1")[cvrt[order(cvrt$Diagnosis),]$Diagnosis])
dev.off()

png("UNKNOWN_NMR2.png")
cvrt = read.table(paste(cov_matrices_dir,"Covariates_all_vertical_text2.txt",sep=""),sep="\t",header=T,row.names=1)
cvrt <- cvrt[rownames(cvrt) %in% colnames(expr_nmr_unknown),]
cvrt <- cvrt[colnames(expr_nmr_unknown),]
expr_nmr_unknown2 <- expr_nmr_unknown[,cvrt$Diagnosis %in% c("CTL","ADC")]
cvrt <- cvrt[colnames(expr_nmr_unknown2),]
t <- heatmap.2(as.matrix(expr_nmr_unknown2),col=brewer.pal(11,"RdBu"), 
               trace="none",dendrogram = "both",
               colCol=c("firebrick3","springgreen4")[as.factor(as.character(cvrt$Diagnosis))])
dev.off()


#2) Alzheimer's related - diagnosis plays role in LR 
res_urine_all <-res_all[0,]
for(base in c("UHPOS","URPOS","URNEG")){
  main_dir <- "/hps/nobackup/ma/natalja/data/mQTL_data/results/normalization_final/"
  output_file_name <- paste(main_dir,"diagnosis_result_",base,".txt",sep="")
  res_all <- read.table(output_file_name,header=T,sep="\t",stringsAsFactors = F)
  res_urine_all <- rbind(res_urine_all,res_all)
}
diagnosis_SNPs_lcms <- unique(res_urine_all$SNP)

res_sig_diag_lcms <- res_sig_lcms[res_sig_lcms$SNP %in% unique(res_urine_all$SNP),]
expr_lcms_diag <- expr_lcms[rownames(expr_lcms) %in% unique(res_sig_diag_lcms$gene),]
expr_lcms_diag_unknown <- expr_lcms[rownames(expr_lcms) %in% unique(res_sig_lcms_unknown$gene) & rownames(expr_lcms) %in% unique(res_sig_diag_lcms$gene),]
expr_lcms_diag_known <- expr_lcms[!(rownames(expr_lcms) %in% unique(res_sig_lcms_unknown$gene)) & rownames(expr_lcms) %in% unique(res_sig_diag_lcms$gene),]

png("DIAG_LCMS.png")
cvrt = read.table(paste(cov_matrices_dir,"Covariates_all_vertical_text2.txt",sep=""),sep="\t",header=T,row.names=1)
cvrt <- cvrt[rownames(cvrt) %in% colnames(expr_lcms_diag),]
cvrt <- cvrt[colnames(expr_lcms_diag),]

t <- heatmap.2(as.matrix(expr_lcms_diag[order(cvrt$Diagnosis)]),col=brewer.pal(11,"RdBu"), 
               trace="none",Colv="none", dendrogram = "row",
               colCol=c("firebrick3","plum4","springgreen4","plum1")[cvrt[order(cvrt$Diagnosis),]$Diagnosis])
dev.off()
png("DIAG_LCMS2.png")
cvrt = read.table(paste(cov_matrices_dir,"Covariates_all_vertical_text2.txt",sep=""),sep="\t",header=T,row.names=1)
cvrt <- cvrt[rownames(cvrt) %in% colnames(expr_lcms_diag),]
cvrt <- cvrt[colnames(expr_lcms_diag),]

t <- heatmap.2(as.matrix(expr_lcms_diag),col=brewer.pal(11,"RdBu"), 
               trace="none",dendrogram = "row",
               colCol=c("firebrick3","plum4","springgreen4","plum1")[cvrt$Diagnosis])
dev.off()

png("DIAG_LCMS3.png")
cvrt = read.table(paste(cov_matrices_dir,"Covariates_all_vertical_text2.txt",sep=""),sep="\t",header=T,row.names=1)
cvrt <- cvrt[rownames(cvrt) %in% colnames(expr_lcms_diag),]
cvrt <- cvrt[colnames(expr_lcms_diag),]
expr_lcms_diag2 <- expr_lcms_diag[,cvrt$Diagnosis %in% c("CTL","ADC")]
cvrt <- cvrt[colnames(expr_lcms_diag2),]
t <- heatmap.2(as.matrix(expr_lcms_diag2),col=brewer.pal(11,"RdBu"), 
               trace="none",dendrogram = "both",
               colCol=c("firebrick3","springgreen4")[as.factor(as.character(cvrt$Diagnosis))])
dev.off()




pdf("DIAG_WEKA_SELECTED_LCMS.pdf")
expr_selected <- expr_lcms[rownames(expr_lcms) %in% selected_mets$V1,]
cvrt = read.table(paste(cov_matrices_dir,"Covariates_all_vertical_text2.txt",sep=""),sep="\t",header=T,row.names=1)
cvrt <- cvrt[rownames(cvrt) %in% colnames(expr_selected),]
cvrt <- cvrt[colnames(expr_selected),]

t <- heatmap.2(as.matrix(expr_selected[order(cvrt$Diagnosis)]),col=brewer.pal(11,"RdBu"), 
               trace="none",Colv="none", dendrogram = "row",cexRow=0.5, key=F, main="Metabolites and diagnoses",
               colCol=c("firebrick3","plum4","springgreen4","plum1")[cvrt[order(cvrt$Diagnosis),]$Diagnosis])

par(lend = 1)           # square line ends for the color legend
legend("topleft",      # location of the legend on the heatmap plot
       legend = c("ADC", "cMCI", "CTL", "sMCI"), # category labels
       col = c("firebrick3","plum4","springgreen4","plum1"),  # color key
       lty= 1,             # line style
       lwd = 10            # line width
)
dev.off()

library(RColorBrewer) 
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
hmcol <- brewer.pal(11,"RdBu")
as.fumeric <- function(x,levels=unique(x)) {
  if(!is.character(x)) stop("'x' must be a character")
  as.numeric(factor(x,levels=levels))
}

selected_mets2 <- read.table("/hps/nobackup/ma/natalja/data/mQTL_data/results/normalization_final/list_selected_attributes2.txt",sep="\t", header=FALSE, stringsAsFactors=F)

pdf("DIAG_WEKA_SELECTED_LCMS2.pdf")
expr_selected <- expr_lcms[rownames(expr_lcms) %in% selected_mets$V1,]

cvrt = read.table(paste(cov_matrices_dir,"Covariates_all_vertical_text2.txt",sep=""),sep="\t",header=T,row.names=1)
cvrt <- cvrt[rownames(cvrt) %in% colnames(expr_selected),]
cvrt <- cvrt[colnames(expr_selected),]

cols <- palette(brewer.pal(8, "Dark2"))[as.fumeric(as.character(cvrt$Diagnosis))]

cols <- c("firebrick3","plum4","springgreen4","plum1")

heatmap.2(as.matrix(expr_selected[order(cvrt$Diagnosis)]),col=hmcol, 
          trace="none",Colv="none", dendrogram = "row",cexRow=0.5, key=F, main="Metabolites and diagnoses",
          #trace="none", distfun=dist.pear, hclustfun=hclust.ave,scale="row",
          #cexRow=0.5, key=F, main="Metabolites and diagnoses",
          colCol=cols[cvrt[order(cvrt$Diagnosis),]$Diagnosis])

par(lend = 1)           # square line ends for the color legend
legend("topleft",      # location of the legend on the heatmap plot
       legend = c("ADC", "cMCI", "CTL", "sMCI"), # category labels
       col =cols,  # color key
       lty= 1,             # line style
       lwd = 10            # line width
)

t <- heatmap.2(as.matrix(expr_selected), 
                         labCol=cvrt$Diagnosis,
                         trace="none", 
                         ColSideColors=cols, 
                         col=hmcol,
                         cexRow=0.5, 
                         key=F, 
                         main="Metabolites and diagnoses")

par(lend = 1)           # square line ends for the color legend
legend("topleft",      # location of the legend on the heatmap plot
       legend = levels(cvrt$Diagnosis), # category labels
       col = palette(brewer.pal(8, "Dark2"))[c(1,2,3,4)],  # color key
       lty= 1,             # line style
       lwd = 10            # line width
)
dev.off()


res_sig_lcms_known <- res_sig_lcms[res_sig_lcms$SNP %in% ALL_snps$V1,]
expr_lcms_known <- expr_lcms[rownames(expr_lcms) %in% unique(res_sig_lcms_known$gene),]

png("ALL_test.png")
cvrt = read.table(paste(cov_matrices_dir,"Covariates_all_vertical_text2.txt",sep=""),sep="\t",header=T,row.names=1)
cvrt <- cvrt[rownames(cvrt) %in% colnames(expr_urine_all),]
cvrt <- cvrt[colnames(expr_urine_all),]

t <- heatmap.2(as.matrix(expr_urine_all[order(cvrt$Diagnosis)]),col=brewer.pal(11,"RdBu"), trace="none",Colv="none",colCol=c("firebrick3","plum4","springgreen4","plum1")[cvrt[order(cvrt$Diagnosis),]$Diagnosis])
dev.off()

png("ALL_test2.png")
cvrt = read.table(paste(cov_matrices_dir,"Covariates_all_vertical_text2.txt",sep=""),sep="\t",header=T,row.names=1)
cvrt <- cvrt[rownames(cvrt) %in% colnames(expr_nmr),]
cvrt <- cvrt[colnames(expr_nmr),]

t <- heatmap.2(as.matrix(expr_nmr[order(cvrt$Diagnosis)]),col=brewer.pal(11,"RdBu"), trace="none",colCol=c("firebrick3","plum4","springgreen4","plum1")[cvrt[order(cvrt$Diagnosis),]$Diagnosis])
dev.off()

library(ggplot2)
library(reshape2) # this is used for melt function
rnames<- rownames(snps_ABO)
df<-cbind(rnames,data.frame(snps_ABO))
df.m<-melt(df)
df.m$value <- factor(df.m$value)
png("ABO_SNPs.png")
rnames <- rownames(snps_ABO)
(p <- ggplot(df.m, aes(variable, rnames)) + geom_tile(aes(fill = value)) +  scale_fill_manual(
  values=c("grey36","red4","forestgreen"),
  breaks=c("0","1","2"),
  labels=c("AA","Aa","aa")))
dev.off()


