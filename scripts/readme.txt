Data preparation for mqtl 

1) The sample samples in DNA. MET and Covariate matrices

Rscript /hps/nobackup/ma/natalja/data/imputed_genotyping_data/filter_columns.R /hps/nobackup/ma/natalja/data/imputed_genotyping_data/dna_matrix_10percentMAF.tsv /hps/nobackup/ma/natalja/data/mQTL_data/NMR/normalized_data/NMR_CPMG_Serum_qn.tsv /hps/nobackup/ma/natalja/data/mQTL_data/dna_matrices/dna_matrix_NMR_CPMG_Serum.tsv;

Rscript /hps/nobackup/ma/natalja/data/imputed_genotyping_data/filter_columns.R  /hps/nobackup/ma/natalja/data/mQTL_data/NMR/normalized_data/NMR_CPMG_Serum_qn.tsv /hps/nobackup/ma/natalja/data/mQTL_data/covariates/Covariates_all_gm_no_diagnosis.txt /hps/nobackup/ma/natalja/data/mQTL_data/NMR/normalized_data/NMR_CPMG_Serum_qn_gm.tsv

Rscript /hps/nobackup/ma/natalja/data/imputed_genotyping_data/filter_columns.R /hps/nobackup/ma/natalja/data/mQTL_data/covariates/Covariates_all_gm_no_diagnosis.txt /hps/nobackup/ma/natalja/data/mQTL_data/NMR/normalized_data/NMR_CPMG_Serum_qn.tsv /hps/nobackup/ma/natalja/data/mQTL_data/covariates/Covariates_NMR_CPMG_Serum_gm.txt

2) The same sample order in DNA. MET and Covariate matrices
match_matrices <- function(matrix1_name, matrix2_name, matrix2_new_name){
  matrix1 <- read.table(matrix1_name, header = TRUE,stringsAsFactors = FALSE ,skip=0, nrows=1)
  matrix2 <- read.table(matrix2_name,header = TRUE,stringsAsFactors = FALSE)
  matrix2 <- matrix2[,intersect(colnames(matrix1),colnames(matrix2))]
  samples <- colnames(matrix2) 
  matrix2 <- cbind(rownames(matrix2),matrix2)
  colnames(matrix2) <- c("id",samples)
  write.table(matrix2,matrix2_new_name,sep="\t",quote = FALSE,row.names = FALSE)
}

match_matrices("/hps/nobackup/ma/natalja/data/mQTL_data/dna_matrices/dna_matrix_NMR_CPMG_Serum.tsv",
               "/hps/nobackup/ma/natalja/data/mQTL_data/met_matrices/NMR/NMR_CPMG_Serum_qn_gm.tsv",
               "/hps/nobackup/ma/natalja/data/mQTL_data/met_matrices/NMR/NMR_CPMG_Serum.tsv")

match_matrices("/hps/nobackup/ma/natalja/data/mQTL_data/dna_matrices/dna_matrix_NMR_CPMG_Serum.tsv",
               "/hps/nobackup/ma/natalja/data/mQTL_data/covariates/Covariates_NMR_CPMG_Serum_gm.txt",
               "/hps/nobackup/ma/natalja/data/mQTL_data/covariates/Covariates_NMR_CPMG_Serum.txt")


Run mQTL analysis
3) 
bsub -M 50000 -R 'rusage[mem=50000]' 'Rscript /hps/nobackup/ma/natalja/data/mQTL_data/scripts/run_mqtl2.R "NMR_CPMG_Serum"'


