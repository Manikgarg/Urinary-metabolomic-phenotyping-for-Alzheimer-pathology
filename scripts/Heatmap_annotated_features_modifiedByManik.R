################################################################################
# Heatmap from selected metabolites
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
hmcol <- colorRampPalette(brewer.pal(9,"GnBu"))(100)
hmcol <- brewer.pal(11,"RdBu")
as.fumeric <- function(x,levels=unique(x)) {
  if(!is.character(x)) stop("'x' must be a character")
  as.numeric(factor(x,levels=levels))
}

source("./scripts/load_data.R")



#method <-"4D"#"corr" #"LR"

#selected_features <- c("5.01_535.1753m/z","5.51_551.1705m/z","2.58_151.0635n","2.31_288.0117n","2.59_318.0226n","3.07_283.1580n","6.13_685.2611m/z","6.13_600.2850m/z","6.13_600.4780m>","6.13_668.2724m/z","6.13_736.2562m/z","6.30_602.2977m/z","5.04_614.2623m/z","5.85_614.2613m/z","8.01_639.3614n","2.03_112.0501m/z","3.80_258.1087m/z","3.80_126.0658m/z","4.55_190.2> 5n","2.62_125.0837n","2.62_92.1133n","3.74_249.0860n","X1.2051","X1.1952","X2.873","X1.7312","1.21_229.2724m/z","1.72_229.2717m/z","X5.2831","X5.276","X5.4085","4.64_637.1565m/z","4.> _197.0803m/z","4.43_202.1581m/z","5.63_163.0859m/z","1.41_95.0125m/z","0.95_384.1467m/z","0.90_204.0860m/z")

# annotations <- c("Diltiazem glucuronide","Diltiazem-O-glucuronide","Paracetamol",
#"Paracetamol sulphate","3-methoxy-paracetamol-sulphate","Quinine",
#"Pregnenolone sulfate N-acetylglucosamine",
#"Pregnenolone sulfate N-acetylglucosamine","Pregnenolone sulfate N-acetylglucosamine",
#"Pregnenolone sulfate N-acetylglucosamine","Pregnenolone sulfate N-acetylglucosamine",
#"Pregnanediol sulfate N-acetylglucosamine",
#"hydroxylated pregnanolone sulfate N-acetylglucosamine",
#"hydroxylated pregnanolone sulfate N-acetylglucosamine",
#"Taurochenodeoxycholic acid N-acetylglucosaminide",
#"2-O-methylcytidine","5-methylcytidine","5-methylcytidine","Butyryl or isobutyryl carnitine","N-acetylisoputreanine-gamma-lactam","N-acetylisoputreanine-gamma-lactam","S-adenosyl-methionine","3-amino-isobutyrate","3-amino-isobutyrate","Trimethylamine","Lysine","N,N,N-trimethyl-l-alanyl-l-proline betaine","N,N,N-trimethyl-l-alanyl-l-proline betaine","Galactose","Galactose","Sucrose","Glutathione","3-hydroxyhippuric acid","Hippuric acid","5-Methoxyindoleacetate","Dimethyl sulfone","N-Acetyllactosamine","N-Acetyllactosamine")

selected_features <- c("5.01_535.1753m/z","5.51_551.1705m/z","2.58_151.0635n","2.31_288.0117n","2.59_318.0226n","3.07_283.1580n","6.13_685.2611m/z","6.13_600.2850m/z","6.13_600.4780m/z","6.13_668.2724m/z","6.13_736.2562m/z","6.30_602.2977m/z","5.04_614.2623m/z","5.85_614.2613m/z","8.01_639.3614n","2.03_112.0501m/z","3.80_258.1087m/z","3.80_126.0658m/z","4.55_190.2275n","2.62_125.0837n","2.62_92.1133n","3.74_249.0860n","X1.2051","X1.1952","X2.873","X1.7312","1.21_229.2724m/z","1.72_229.2717m/z","X5.2831","X5.276","X5.4085","4.64_637.1565m/z")
annotations <- c("Diltiazem glucuronide","Diltiazem-O-glucuronide","Paracetamol","Paracetamol sulphate","3-methoxy-paracetamol-sulphate","Quinine","Pregnenolone sulfate *","Pregnenolone sulfate *","Pregnenolone sulfate *","Pregnenolone sulfate *","Pregnenolone sulfate *","Pregnanediol sulfate *","Hydroxylated pregnanolone sulfate *","Hydroxylated pregnanolone sulfate *","Taurochenodeoxycholic or Taurodeoxycholic acid **","2-O-methylcytidine","5-methylcytidine","5-methylcytidine","Butyryl or isobutyryl carnitine","N-acetylisoputreanine-gamma-lactam","N-acetylisoputreanine-gamma-lactam","S-adenosyl-methionine","3-amino-isobutyrate","3-amino-isobutyrate","Trimethylamine","Lysine","N,N,N-trimethyl-l-alanyl-l-proline betaine","N,N,N-trimethyl-l-alanyl-l-proline betaine","Galactose","Galactose","Sucrose","Glutathione")

selected_mets <- data.frame("Feature" = selected_features,"Annotation" = annotations)
expr_selected <- expr[rownames(expr) %in% selected_mets$Feature,]

cvrt = read.table("./covariates/Covariates_all_vertical_text2.txt",sep="\t",header=T,row.names=1)
cvrt <- cvrt[rownames(cvrt) %in% colnames(expr_selected),]
cvrt <- cvrt[colnames(expr_selected),]
df <- as.data.frame(t(expr_selected))
df$Diagnosis <- as.character(cvrt$Diagnosis)

cvrt$Diagnosis<-as.character(cvrt$Diagnosis)
cvrt[cvrt$Diagnosis=="ADC",]$Diagnosis <-"AD"
cvrt$Diagnosis<-as.factor(cvrt$Diagnosis)
###############################
# Heatmap of selected features 
###############################
diagnoses <- c("CTL","sMCI","cMCI","AD")
cols <- c("#00AFBB","#E7B800","#FC4E07","#4F7942")
#re-order levels to correspond to ggplot2 colors and order
cvrt$Diagnosis <- factor(cvrt$Diagnosis, levels = diagnoses)

m <- as.matrix(expr_selected[order(cvrt$Diagnosis)])
rownames(selected_mets) <- selected_mets$Feature
selected_mets <- selected_mets[rownames(selected_mets) %in% rownames(m),]
selected_mets <- selected_mets[rownames(m),]
rownames(m) <- selected_mets$Annotation

# png("heatmap_annotated.png", width = 14, height = 9, units = 'in', res=300)
# #pdf("heatmap_annotated.pdf")
# t <- heatmap.2(m,
# 	colsep=c(133, 273, 309), sepwidth=2 , sepcolor=c("black","black","black"), 
# 	col=brewer.pal(11,"RdBu"), 
# 	trace="none",Colv="none", 
# 	#dendrogram ="both",
# 	#scale = "row",
# 	dendrogram ="row",
# 	cexRow=0.9,#key=F, 
# 	colCol=cols[cvrt[order(cvrt$Diagnosis),]$Diagnosis],
# 	margins=c(10,15),key=TRUE, keysize = 1, symkey=FALSE, density.info="none",srtRow=45)
# 
# # par(lend = 1)           # square line ends for the color legend
# # legend("topleft",      # location of the legend on the heatmap plot
# #        legend = diagnoses, # category labels
# #        col = cols,  # color key
# #        lty= 1,             # line style
# #        lwd = 10            # line width
# # )
# dev.off()
# 
# ##############Manik#################
# expr_selected = read.table("./Desktop/expr_selected.tsv", sep = "\t", quote = "", header = TRUE)
# cvrt = read.table("./Desktop/cvrt.tsv", sep = "\t", quote = "", header = TRUE)
# selected_mets = read.table("./Desktop/selected_mets.tsv", sep = "\t", quote = "", header = TRUE)
# 
# cvrt$Diagnosis <- factor(cvrt$Diagnosis, levels = diagnoses)
# 
# m <- as.matrix(expr_selected[order(cvrt$Diagnosis)])
# rownames(selected_mets) <- selected_mets$Feature
# selected_mets <- selected_mets[rownames(selected_mets) %in% rownames(m),]
# selected_mets <- selected_mets[rownames(m),]
# rownames(m) <- selected_mets$Annotation
# 
# library(pheatmap)
# ann_colors = list(
#   Gender = c(Male = "gray0", Female = "gray64"),
#   Diagnosis = c(CTL = "#00AFBB", sMCI = "#E7B800", cMCI = "#FC4E07", AD = "#4F7942")
#   )
# pheatmap(m, cluster_rows=TRUE,cluster_cols=FALSE, 
#          #clustering_distance_cols = "correlation",
#          annotation_col = cvrt[, c("Diagnosis", "Gender")], 
#          show_colnames = F, 
#          show_rownames = T, angle_col = "45",
#          annotation_colors = ann_colors,
#          color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
#          cellheight = 8,
#          #cellwidth = 2,
#          fontsize = 6
#          )

## The source code is taken from https://rdrr.io/bioc/Pigengene/man/pheatmap.type.html and modified
pheatmap.type <- function(Data, annRow, type=colnames(annRow)[1], doTranspose=FALSE, conditions="Auto", ...){
  ## annRow: A data frame with row names the same as row names of Data.
  ## type: The column name of annRow representing two or more conditions.
  ## This function first performs hierarchical clustering on samples
  ## (rows of Data) within each condition.
  ##^Then plots a heatmap without further clustering of rows. 
  ## ... are passed to pheatmap function.
  res <- list()
  
  ######## My modification to include gender information ###########
  annRowOriginal = annRow
  #################################################################
  
  annRow <- annRow[, type, drop=FALSE]
  ## QC:
  if(is.null(rownames(annRow)))
    stop("annrow must have row names!")
  if(any(! rownames(annRow) %in% rownames(Data)))
    stop("annRow has rows that are not present in Data!")
  ## Put all samples in the same condition together.
  annRow <- annRow[order(annRow[, 1]), , drop=FALSE]
  samplesOriginalOrder <- rownames(Data)
  ##Data <- Data[rownames(annRow), , drop=FALSE]
  if(conditions=="Auto")
    conditions <- unique(as.character(annRow[, 1]))
  if(any(!conditions %in% unique(as.character(annRow[, 1])))){
    warning("Some of the conditions are not in annRow.")
    conditions <- intersect(conditions, unique(as.character(annRow[, 1])))
  }
  pheatmapS <- list()
  dataPlot <- c()
  ann1 <- c()
  for(cond in conditions){
    condSamples <- rownames(annRow)[which(annRow==cond)]  
    pa <- pheatmap(Data[condSamples, , drop=FALSE], cluster_cols=FALSE, silent=TRUE)
    pheatmapS[[as.character(cond)]] <- pa
    o2 <- pa$tree_row$order
    dataPlot <- rbind(dataPlot, Data[condSamples[o2], , drop=FALSE])
    ann1 <- rbind(ann1, annRow[condSamples[o2], , drop=FALSE])
  }
  ######## My modification to include gender information ###########
  ann1 = data.frame(Diagnosis = ann1, Gender = annRowOriginal[rownames(ann1), "Gender"])
  ####################################
  if(!doTranspose){
    pAll <- pheatmap(dataPlot, annotation_row=ann1, cluster_rows=FALSE, ...)
  } else { ## Transpose
    pAll <- pheatmap(t(dataPlot), annotation_col=ann1, cluster_cols=FALSE, ...)
  }
  res[["pheatmapS"]] <- pheatmapS
  res[["pheat"]] <- pAll
  res[["ordering"]] <- match(rownames(dataPlot), samplesOriginalOrder)
  res[["annRowAll"]] <- ann1
  invisible(res)
  #return(ann1)
}

pheatmap.type(t(m), annRow = cvrt[, c("Diagnosis", "Gender")], 
              type = colnames(cvrt[, c("Diagnosis", "Gender")])[1],
              #doTranspose=FALSE, 
              conditions=c("CTL", "sMCI", "cMCI", "AD"),
              show_rownames = F, angle_col = "45",
              annotation_colors = ann_colors,
              color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
              #cellheight = 4,
              cellwidth = 12,
              fontsize = 10, 
              filename = "Heatmap_annotated_clusteredWithinGroupsWithGenderInfo.pdf")


