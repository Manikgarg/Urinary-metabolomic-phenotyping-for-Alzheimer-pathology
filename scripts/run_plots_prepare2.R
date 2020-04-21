library(biomaRt)
library(qqman)

args <- commandArgs(trailingOnly=TRUE)
base <- args[1]
diagnosis <- args[2]


setwd("/hps/nobackup/ma/natalja/data/mQTL_data/results/")
source("/hps/nobackup/ma/natalja/data/mQTL_data/scripts/mqtl.R")

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



# CIRCOS TRACKS
res_sig <- res[res$FDR<0.01,]
nt.biomart <- getBM(c("refsnp_id","allele","chr_name","chrom_start","chrom_end","ensembl_gene_stable_id"),
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

res_genes <- res[!is.na(res$hgnc_symbol) & res$hgnc_symbol!="",]
res_genes<-res_genes[-grep("CHR_", res_genes$chr_name),]
#Text track with lead snps gene names
res_text<-data.frame(paste("hs",res_genes$chr_name,sep=""),
                     res_genes$chrom_start,
                     res_genes$chrom_end,res_genes$hgnc_symbol)
colnames(res_text)<-c("CHR","START","END","GENE")
res_text <- res_text[!duplicated(res_text$GENE),]
res_text$END<-res_text$END+200
write.table(res_text,paste("lead_snp_genes_",base,diagnosis,".txt",sep="")
            ,sep="\t",row.names = FALSE)

#P-values track with all significant snps
res_snps<-res_annot[-grep("CHR_", res_annot$chr_name),]
#res_snps$bonf <- p.adjust(res_snps$p.value,"bonferroni")
logp <- -log10(res_snps$p.value)
logp <- logp - min(logp)
res_pvalues<-data.frame(paste("hs",res_snps$chr_name,sep=""),
                        res_snps$chrom_start,
                        res_snps$chrom_end,logp)
colnames(res_pvalues)<-c("CHR","START","END","P_VALUE")
res_pvalues$END<-res_pvalues$END+200
write.table(res_pvalues,paste("log_p_values_",base,diagnosis,".txt",sep="")
            ,sep="\t",row.names = FALSE)


