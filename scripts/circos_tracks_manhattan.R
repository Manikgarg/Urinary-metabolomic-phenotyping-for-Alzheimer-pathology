# Circos tracks
args <- commandArgs(trailingOnly=TRUE)
main_dir <- args[1]
base <- args[2]

setwd(main_dir)

options("scipen"=100, "digits"=4) # to prevent scientific notation

# Manhattan plot track - all p-values, significant ADCvsCTL are coloured red
bases <- c( "SHPOS","SLPOS","UHPOS","URPOS","URNEG")
for(base in bases){
 output_file_name_sig <- paste("result_",base,"_annot_full.txt",sep="") 
 output_file_name <- paste("mqtl_annot_",base,".txt",sep="") # all p-values with annotations
 res_all <- read.table(output_file_name,header=T,sep="\t",stringsAsFactors=FALSE)
 res_sig <- read.table(output_file_name_sig,header=T,sep="\t",stringsAsFactors=FALSE)
 met_list <- read.table(paste("./metabolites/",base,"_ADC_CTL_metabolites.txt",sep=""),stringsAsFactors = FALSE,header = FALSE) # list of metabolites ADCvsCTL
   
 res_sig <- res_sig[res_sig$gene %in% met_list$V1,] # ADCvsCTL associations only
 
 res_all$col <- ""
 res_all[res_all$SNP %in% res_sig$SNP & res_all$gene %in% res_sig$gene,]$col <- "color=red"
 
 res_plot<-data.frame(SNP=res_all$SNP,
                      CHR=res_all$chr_name,
                      P=res_all$p.value, # FDR or original p-value?????????????????
                      FDR=res_all$FDR,
                      BP=res_all$chrom_start,
                      COL= res_all$col,stringsAsFactors=FALSE)
 res_plot<-res_plot[-grep("CHR_", res_plot$CHR),]
 res_plot$CHR[as.character(res_plot$CHR) == "X"] <- "23"
 res_plot$CHR <- as.numeric(res_plot$CHR)
 
 valmax <- max(res_plot$P)
 res_plot$P <- round(-log10(res_plot$P/valmax),4)
 
 valmaxFDR <- max(res_plot$FDR)
 res_plot$FDR <- round(-log10(res_plot$FDR/valmaxFDR),4)

 res_pvalues<-data.frame(paste("hs",res_plot$CHR,sep=""),
                         res_plot$BP,
                         res_plot$BP,
                         res_plot$P,
                         #res_plot$FDR,
                         res_plot$COL)
 colnames(res_pvalues)<-c("CHR","START","END","P_VALUE","COLOUR")
 res_pvalues$END<-res_pvalues$END+200
 
 write.table(res_pvalues,paste("log_p_values_",base,"_full.txt",sep="")
                                   ,sep="\t",col.names = F, row.names = FALSE,quote = FALSE)

}

# Gene names, significant ADCvsCTL genes are coloured red
for(base in bases){
  output_file_name_sig <- paste("result_",base,"_annot_full.txt",sep="") 
  res_sig <- read.table(output_file_name_sig,header=T,sep="\t",stringsAsFactors=FALSE)
  met_list <- read.table(paste("./metabolites/",base,"_ADC_CTL_metabolites.txt",sep=""),stringsAsFactors = FALSE,header = FALSE) # list of metabolites ADCvsCTL
  
  res_sig$col <- ""
  res_sig[res_sig$gene %in% met_list$V1,] $col <- "color=red"
  res_sig <- res_sig[!is.na(res_sig$hgnc_symbol),]

  
  res_genes<-data.frame(paste("hs",res_sig$chr_name,sep=""),
                          res_sig$chrom_start,
                          res_sig$chrom_start,
                          res_sig$hgnc_symbol,
                          #res_plot$FDR,
                          res_sig$col)
  colnames(res_genes)<-c("CHR","START","END","GENE","COLOUR")
  res_genes$END<-res_genes$END+200
  res_genes <- res_genes[!duplicated(res_genes$GENE),]
  write.table(res_genes,paste("text_",base,"_full.txt",sep="")
              ,sep="\t",col.names = F, row.names = FALSE,quote = FALSE)
  
}

base <- "SLNEG"
output_file_name_sig <- paste("result_",base,"_annot_full.txt",sep="") 
output_file_name <- paste("mqtl_annot_",base,".txt",sep="") # all p-values with annotations
res_all <- read.table(output_file_name,header=T,sep="\t",stringsAsFactors=FALSE)
res_sig <- read.table(output_file_name_sig,header=T,sep="\t",stringsAsFactors=FALSE)


res_plot<-data.frame(SNP=res_all$SNP,
                     CHR=res_all$chr_name,
                     P=res_all$p.value, # FDR or original p-value?????????????????
                     FDR=res_all$FDR,
                     BP=res_all$chrom_start,
                     stringsAsFactors=FALSE)
  
  res_plot<-res_plot[-grep("CHR_", res_plot$CHR),]
  res_plot$CHR[as.character(res_plot$CHR) == "X"] <- "23"
  res_plot$CHR <- as.numeric(res_plot$CHR)

valmax <- max(res_plot$P)
res_plot$P <- round(-log10(res_plot$P/valmax),4)

valmaxFDR <- max(res_plot$FDR)
res_plot$FDR <- round(-log10(res_plot$FDR/valmaxFDR),4)

res_pvalues<-data.frame(paste("hs",res_plot$CHR,sep=""),
                        res_plot$BP,
                        res_plot$BP,
                        res_plot$P)
                        #res_plot$FDR)
colnames(res_pvalues)<-c("CHR","START","END","P_VALUE")
res_pvalues$END<-res_pvalues$END+200

write.table(res_pvalues,paste("log_p_values_",base,"_full.txt",sep="")
            ,sep="\t",col.names = F, row.names = FALSE,quote = FALSE)

res_genes<-data.frame(paste("hs",res_sig$chr_name,sep=""),
                      res_sig$chrom_start,
                      res_sig$chrom_start,
                      res_sig$hgnc_symbol)
                      #res_plot$FDR)
colnames(res_genes)<-c("CHR","START","END","GENE")
res_genes$END<-res_genes$END+200
res_genes <- res_genes[!duplicated(res_genes$GENE),]
write.table(res_genes,paste("text_",base,"_full.txt",sep="")
            ,sep="\t",col.names = F, row.names = FALSE,quote = FALSE)