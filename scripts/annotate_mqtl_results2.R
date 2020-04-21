# Annotate significant results with gene info
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

#selects significant associations only
res_sig <- res[res$FDR<0.01,]
snp.db <- useMart("ENSEMBL_MART_SNP", dataset="hsapiens_snp")
#SNPs annotation
nt.biomart <- getBM(c("refsnp_id","allele","chr_name","chrom_start","chrom_end","ensembl_gene_stable_id"),
                    filters="snp_filter",
                    values=res_sig$SNP,
                    mart=snp.db)
res_annot <- merge(res_sig,nt.biomart,by.x="SNP",by.y="refsnp_id",all.x=T,all.y=F)
genes <- unique(nt.biomart$ensembl_gene_stable_id)

#Genes annotation
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
genes_annot <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol','chromosome_name','description'), 
                     filters = 'ensembl_gene_id', values = genes, mart = ensembl)

res_annot <- merge(res_annot,genes_annot,by.x="ensembl_gene_stable_id",by.y="ensembl_gene_id",all.x=T,all.y=F)
write.table(res_annot,paste("result_",base,"_annot_full.txt",sep=""),sep="\t",quote = F,row.names = F)


#A tag SNP is a representative single nucleotide polymorphism (SNP) in a region of the genome 
#with high linkage disequilibrium that represents a group of SNPs called a haplotype. 
#It is possible to identify genetic variation and association to phenotypes without genotyping every SNP in a chromosomal region.
#12105785 SNPs
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

write.table(res,paste("result_",base,"_annot.txt",sep=""),sep="\t",quote = F,row.names = F)
