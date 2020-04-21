# Annotate mqtl_<base>.txt files with SNPs info
# works in chunks - 20000 records
require(biomaRt)

args <- commandArgs(trailingOnly=TRUE)
main_dir <- args[1]
base <- args[2]

setwd(main_dir)

#reas mqtl results
output_file_name <- paste("mqtl_",base,".txt",sep="") 
res <- read.table(output_file_name,header=T,sep="\t")

#rs id format
res$SNP <- gsub('_2','',res$SNP)
res$SNP <- gsub('_1','',res$SNP)

#split results into chunks 
res_split <- split(res, (seq(nrow(res))-1) %/% 20000) 
res_all <- res[0,]

snp.db <- useMart("ENSEMBL_MART_SNP", dataset="hsapiens_snp")
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

#clean by removing NAs
res_all <- res_all[!is.na(res_all$chr_name) & res_all$chr_name!="",]
res_all <- res_all[!is.na(res_all$p.value) & res_all$p.value!="",]
res_all <- res_all[!is.na(res_all$chrom_start) & res_all$chrom_start!="",]
res_all <- res_all[!is.na(res_all$SNP) & res_all$SNP!="",]

write.table(res_all,paste("mqtl_",base,"_annot.txt",sep=""),sep="\t",quote = F,row.names = F)


