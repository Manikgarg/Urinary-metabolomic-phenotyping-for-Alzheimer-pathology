# NMR data
nmr_p1 <- read.table("NMR Nosey (SERUM)_part1.csv",header = T,sep=";",stringsAsFactors = FALSE,dec = ",")
nmr_p2 <- read.table("NMR Nosey (SERUM)_part2.csv",header = T,sep=";",stringsAsFactors = FALSE,dec= ",")


pheno_data1 <- nmr_p1[,1:12]
pheno_data2 <- nmr_p2[,1:12]
nmr_p1 <- nmr_p1[,-(2:12)]
nmr_p2 <- nmr_p2[,-(2:12)]
# For now remove QC
nmr_p1 <- nmr_p1[nmr_p1$Sample.ID!="QC",]
nmr_p2 <- nmr_p2[nmr_p2$Sample.ID!="QC",]
nmr_p1 <- nmr_p1[nmr_p1$Sample.ID!="",]
nmr_p2 <- nmr_p2[nmr_p2$Sample.ID!="",]
#637 samples
rownames(nmr_p1)<- nmr_p1$Sample.ID
rownames(nmr_p2)<- nmr_p2$Sample.ID
nmr <- merge(nmr_p1,nmr_p2,by="Sample.ID")
# 637 x 18647
nmr<-t(nmr)
colnames(nmr)<- nmr[1,]
nmr <- nmr[-1,]
class(nmr) <- "numeric"



# Lipidomics
nmr_p1 <- read.table("Lipidomics Positive (URINE)_part1.csv",header = T,sep=",",stringsAsFactors = FALSE,dec = ".")
nmr_p2 <- read.table("Lipidomics Positive (URINE)_part2.csv",header = T,sep=",",stringsAsFactors = FALSE,dec= ".")
nmr_p3 <- read.table("Lipidomics Positive (URINE)_part3.csv",header = T,sep=",",stringsAsFactors = FALSE,dec= ".")

nmr_p1 <- read.table("Lipidomics Negative (URINE)_part1.csv",header = T,sep=",",stringsAsFactors = FALSE,dec = ".")
nmr_p2 <- read.table("Lipidomics Negative (URINE)_part2.csv",header = T,sep=",",stringsAsFactors = FALSE,dec= ".")

nmr_p1 <- read.table("Lipidomics Positive (SERUM)_part1.csv",header = T,sep=",",stringsAsFactors = FALSE,dec = ".")
nmr_p2 <- read.table("Lipidomics Positive (SERUM)_part2.csv",header = T,sep=",",stringsAsFactors = FALSE,dec= ".")

rownames(nmr_p1)<- nmr_p1$Sample.ID
rownames(nmr_p2)<- nmr_p2$Sample.ID
rownames(nmr_p3)<- nmr_p3$Sample.ID
nmr12 <- merge(nmr_p1,nmr_p2,by="Sample.ID")
nmr <- merge(nmr12,nmr_p3,by="Sample.ID")
nmr<-t(nmr)
colnames(nmr)<- nmr[1,]
nmr <- nmr[-1,]
class(nmr) <- "numeric"



source(file="mswsd_resamp_publi.R")
unbal_reg(nmr)
resamp_mswsd(nmr)

norm.my.data<-norm_unbal(nmr,70,"VSN")
write.table(nmr,"NMR_Nosey_Serum_original.txt",sep="\t",quote = F,row.names = T)
write.table(norm.my.data,"NMR_Nosey_Urine_norm_vsn.txt",sep="\t",quote = F,row.names = T)
norm.my.data<-norm_unbal(nmr,60,"PQN")
write.table(norm.my.data,"NMR_Nosey_Urine_norm_pqn.txt",sep="\t",quote = F,row.names = T)

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
bases <- c("NMR_CPMG_Serum","NMR_Lipidomics_Negative_Urine","NMR_Lipidomics_Positive_Urine","NMR_Lipidomics_Negative_Serum","NMR_Lipidomics_Positive_Serum","NMR_Nosey_Serum","NMR_Nosey_Urine") 
# Normalize all metabolic data
for(base in bases){
  met <- paste("/home/natalja/Documents/EMIF/KCL/NPC/NPC/NMR/original_data/",base,"_original.txt",sep="")
  met_qn <- paste("/home/natalja/Documents/EMIF/KCL/NPC/NPC/NMR/normalized_data/",base,"_qn.tsv",sep="")
  QN_normalize(met,met_qn)
}


