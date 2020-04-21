################################################################################
# Run metabolic QTL analysis for LC-MS data
# 1. mQTL with original matrices
# 2. mQTL with permuted matrices
# 3. Results are annotated with bioMart
################################################################################
source("/hps/nobackup/ma/natalja/data/mQTL_data/scripts/auxiliary_functions.R")
library(biomaRt)
library(MatrixEQTL)

args <- commandArgs(trailingOnly=TRUE)
base <- args[1]
diagnosis <- ""
################################################################################
# Function to run mQTL analysis
################################################################################
mqtl_analysis <- function(snps,met_matrix_name,cov_name,output_file_name){ 
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
  #plot(me, pch = 16, cex = 0.7)
}  
################################################################################
################################################################################

# load snps 
SNP_file_name <- paste("/hps/nobackup/ma/natalja/data/mQTL_data/dna_matrices/dna_matrix_",base,diagnosis,".tsv",sep="") 
snps = SlicedData$new();
snps$fileDelimiter = "\t"; # the TAB character
snps$fileOmitCharacters = 'NA' ;# denote missing values;
snps$fileSkipRows = 1; # one row of column labels
snps$fileSkipColumns = 1; # one column of row labels
snps$fileSliceSize = 10000; # read file in pieces of 10,000 rows
snps$LoadFile(SNP_file_name);

met_matrix_name <-paste("/hps/nobackup/ma/natalja/data/mQTL_data/met_matrices/normalization_final/",base,diagnosis,".tsv",sep="")
cov_name <- paste("/hps/nobackup/ma/natalja/data/mQTL_data/covariates/Covariates_",base,diagnosis,".txt",sep="") 
output_file_name <- paste("/hps/nobackup/ma/natalja/data/mQTL_data/results/normalization_final/mqtl_",base,diagnosis,".txt",sep="") 
  
#mQTL Analysis
mqtl_analysis(snps,met_matrix_name,cov_name,output_file_name)

permuted_met_matrix_name <- paste(gsub(".tsv","",met_matrix_name),"_permuted.tsv",sep="")
permute_matrix(met_matrix_name,permuted_met_matrix_name)
permuted_output_file_name <- paste("/hps/nobackup/ma/natalja/data/mQTL_data/results/normalization_final/mqtl_permuted_",base,diagnosis,".txt",sep="") 

#mQTL analysis after permutation 
mqtl_analysis(snps,permuted_met_matrix_name,cov_name,permuted_output_file_name)
  
res <- read.table(output_file_name,header=T,sep="\t")
res_permuted <- read.table(permuted_output_file_name,header=T,sep="\t")
res_sig <- res[res$FDR<0.01,]

if (nrow(res_permuted[res_permuted$FDR<0.01,])==0){
    res_sig$SNP <- gsub('_2','',res_sig$SNP)
    res_sig$SNP <- gsub('_1','',res_sig$SNP)
    
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
    write.table(res,paste("/hps/nobackup/ma/natalja/data/mQTL_data/results/normalization_final/result_",base,diagnosis,"_annot.txt",sep="") ,sep="\t",row.names = FALSE)
    write.table(res_annot,paste("/hps/nobackup/ma/natalja/data/mQTL_data/results/normalization_final/result_",base,diagnosis,"_annot_full.txt",sep=""),sep="\t",row.names = FALSE)
}


