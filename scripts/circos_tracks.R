library(biomaRt)
snp.db <- useMart("ENSEMBL_MART_SNP", dataset="hsapiens_snp")

args <- commandArgs(trailingOnly=TRUE)
base <- args[1]
diagnosis <- args[2]

setwd("/home/natalja/Documents/EMIF/KCL/mQTL/mQTL_data/results/")
base <-"SLPOS"
diagnosis<-""
output_file_name <- paste("mqtl_",base,diagnosis,".txt",sep="") 
res <- read.table(output_file_name,header=T,sep="\t")

res$SNP <- gsub('_2','',res$SNP)
res$SNP <- gsub('_1','',res$SNP)

#hist(res$p.value,breaks = 50,xlab = "P-values",col = "paleturquoise3",main = paste(base,diagnosis,sep=""))

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

#Plot
library(qqman)
res_plot<-data.frame(SNP=res_all$SNP,CHR=res_all$chr_name,P=res_all$p.value,
                     BP=res_all$chrom_start,stringsAsFactors=FALSE)
res_plot$CHR[as.character(res_plot$CHR) == "X"] <- "23"
res_plot$CHR <- as.numeric(res_plot$CHR)


manhattan(res_plot, col = c("blue4", "orange3"), ymax = 12)





res_sig <- res[res$FDR<0.01,]

res_sig$SNP <- gsub('_2','',res_sig$SNP)
res_sig$SNP <- gsub('_1','',res_sig$SNP)
  
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

#Links track with all significant snps
res_links<-data.frame(paste("hs",res_snps$chr_name,sep=""),
                        res_snps$chrom_start,
                        res_snps$chrom_end+200,
                      paste("hs",res_snps$chr_name,sep=""),
                      res_snps$chrom_start,
                      res_snps$chrom_end+200)
colnames(res_links)<-c("CHR1","START1","END1","CHR2","START2","END2")

#GWAS catalogue
res_annot <- read.table(paste("result_",base,diagnosis,"_annot.txt",sep=""),header = T,sep="\t")
list <- read.table(paste("./metabolites/",base,"_ADC_CTL_metabolites.txt",sep=""),header = F,sep="\t",stringsAsFactors = F)
write.table(res_annot[res_annot$gene %in% list$V1,],paste("./metabolites/",base,"_selected.txt",sep=""),sep="\t",quote = F,row.names = F)
library(gwascat)
data(ebicat38)
#5671
res_genes <- res_annot[!is.na(res_annot$chr_name) & res_annot$chr_name!="",]
#5336
res_genes<-res_genes[-grep("CHR_", res_genes$chr_name),]
#5273
intersect(getRsids(ebicat38), res_genes$SNP)
#rs80088139
intr = ebicat38[ intersect(getRsids(ebicat38), res_genes$SNP) ]
sort(table(getTraits(intr)), decreasing=TRUE)[1:10]

#Traits by genomic location
gr6.0 = GRanges(seqnames=res_genes$chr_name,
                 IRanges(res_genes$chrom_start, width=1))
mcols(gr6.0)$rsid = as.character(res_genes$SNP)
#seqlevels(gr6.0) = paste("chr", seqlevels(gr6.0), sep="")
ag = function(x) as(x, "GRanges")
ovraw = suppressWarnings(subsetByOverlaps(ag(ebicat38), gr6.0))
length(ovraw)
ovaug = suppressWarnings(subsetByOverlaps(ag(ebicat38+1000), gr6.0))
length(ovaug)
rawrs = mcols(ovraw)$SNPS
augrs = mcols(ovaug)$SNPS
ebicat38[augrs]
getTraits(ebicat38[augrs])


# Manhattam plots
library(lattice)
manhattan.plot<-function(chr, pos, pvalue, 
                         sig.level=NA, annotate=NULL, ann.default=list(),
                         should.thin=T, thin.pos.places=2, thin.logp.places=2, 
                         xlab="Chromosome", ylab=expression(-log[10](p-value)),
                         col=c("gray","darkgray"), panel.extra=NULL, pch=20, cex=0.8,...) {
  
  if (length(chr)==0) stop("chromosome vector is empty")
  if (length(pos)==0) stop("position vector is empty")
  if (length(pvalue)==0) stop("pvalue vector is empty")
  
  #make sure we have an ordered factor
  if(!is.ordered(chr)) {
    chr <- ordered(chr)
  } else {
    chr <- chr[,drop=T]
  }
  
  #make sure positions are in kbp
  if (any(pos>1e6)) pos<-pos/1e6;
  
  #calculate absolute genomic position
  #from relative chromosomal positions
  posmin <- tapply(pos,chr, min);
  posmax <- tapply(pos,chr, max);
  posshift <- head(c(0,cumsum(posmax)),-1);
  names(posshift) <- levels(chr)
  genpos <- pos + posshift[chr];
  getGenPos<-function(cchr, cpos) {
    p<-posshift[as.character(cchr)]+cpos
    return(p)
  }
  
  #parse annotations
  grp <- NULL
  ann.settings <- list()
  label.default<-list(x="peak",y="peak",adj=NULL, pos=3, offset=0.5, 
                      col=NULL, fontface=NULL, fontsize=NULL, show=F)
  parse.label<-function(rawval, groupname) {
    r<-list(text=groupname)
    if(is.logical(rawval)) {
      if(!rawval) {r$show <- F}
    } else if (is.character(rawval) || is.expression(rawval)) {
      if(nchar(rawval)>=1) {
        r$text <- rawval
      }
    } else if (is.list(rawval)) {
      r <- modifyList(r, rawval)
    }
    return(r)
  }
  
  if(!is.null(annotate)) {
    if (is.list(annotate)) {
      grp <- annotate[[1]]
    } else {
      grp <- annotate
    } 
    if (!is.factor(grp)) {
      grp <- factor(grp)
    }
  } else {
    grp <- factor(rep(1, times=length(pvalue)))
  }
  
  ann.settings<-vector("list", length(levels(grp)))
  ann.settings[[1]]<-list(pch=pch, col=col, cex=cex, fill=col, label=label.default)
  
  if (length(ann.settings)>1) { 
    lcols<-trellis.par.get("superpose.symbol")$col 
    lfills<-trellis.par.get("superpose.symbol")$fill
    for(i in 2:length(levels(grp))) {
      ann.settings[[i]]<-list(pch=pch, 
                              col=lcols[(i-2) %% length(lcols) +1 ], 
                              fill=lfills[(i-2) %% length(lfills) +1 ], 
                              cex=cex, label=label.default);
      ann.settings[[i]]$label$show <- T
    }
    names(ann.settings)<-levels(grp)
  }
  for(i in 1:length(ann.settings)) {
    if (i>1) {ann.settings[[i]] <- modifyList(ann.settings[[i]], ann.default)}
    ann.settings[[i]]$label <- modifyList(ann.settings[[i]]$label, 
                                          parse.label(ann.settings[[i]]$label, levels(grp)[i]))
  }
  if(is.list(annotate) && length(annotate)>1) {
    user.cols <- 2:length(annotate)
    ann.cols <- c()
    if(!is.null(names(annotate[-1])) && all(names(annotate[-1])!="")) {
      ann.cols<-match(names(annotate)[-1], names(ann.settings))
    } else {
      ann.cols<-user.cols-1
    }
    for(i in seq_along(user.cols)) {
      if(!is.null(annotate[[user.cols[i]]]$label)) {
        annotate[[user.cols[i]]]$label<-parse.label(annotate[[user.cols[i]]]$label, 
                                                    levels(grp)[ann.cols[i]])
      }
      ann.settings[[ann.cols[i]]]<-modifyList(ann.settings[[ann.cols[i]]], 
                                              annotate[[user.cols[i]]])
    }
  }
  rm(annotate)
  
  #reduce number of points plotted
  if(should.thin) {
    logp <- round(-log10(pvalue),thin.logp.places)
    thinned <- unique(data.frame(
      logp=logp <- logp - min(logp), 
      pos=round(genpos,thin.pos.places), 
      chr=chr,
      grp=grp)
    )
    logp <- thinned$logp
    genpos <- thinned$pos
    chr <- thinned$chr
    grp <- thinned$grp
    rm(thinned)
  } else {
    logp <- -log10(pvalue)
    logp <- logp - min(logp)
  }
  rm(pos, pvalue)
  gc()
  
  #custom axis to print chromosome names
  axis.chr <- function(side,...) {
    if(side=="bottom") {
      panel.axis(side=side, outside=T,
                 at=((posmax+posmin)/2+posshift),
                 labels=levels(chr), 
                 ticks=F, rot=0,
                 check.overlap=F
      )
    } else if (side=="top" || side=="right") {
      panel.axis(side=side, draw.labels=F, ticks=F);
    }
    else {
      axis.default(side=side,...);
    }
  }
  
  #make sure the y-lim covers the range (plus a bit more to look nice)
  prepanel.chr<-function(x,y,...) { 
    A<-list();
    maxy<-ceiling(max(y, ifelse(!is.na(sig.level), -log10(sig.level), 0)))+.5;
    A$ylim=c(0,maxy);
    A;
  }
  
  xyplot(logp~genpos, chr=chr, groups=grp,
         axis=axis.chr, ann.settings=ann.settings, 
         prepanel=prepanel.chr, scales=list(axs="i"),
         panel=function(x, y, ..., getgenpos) {
           if(!is.na(sig.level)) {
             #add significance line (if requested)
             panel.abline(h=-log10(sig.level), lty=2);
           }
           panel.superpose(x, y, ..., getgenpos=getgenpos);
           if(!is.null(panel.extra)) {
             panel.extra(x,y, getgenpos, ...)
           }
         },
         panel.groups = function(x,y,..., subscripts, group.number) {
           A<-list(...)
           #allow for different annotation settings
           gs <- ann.settings[[group.number]]
           A$col.symbol <- gs$col[(as.numeric(chr[subscripts])-1) %% length(gs$col) + 1]    
           A$cex <- gs$cex[(as.numeric(chr[subscripts])-1) %% length(gs$cex) + 1]
           A$pch <- gs$pch[(as.numeric(chr[subscripts])-1) %% length(gs$pch) + 1]
           A$fill <- gs$fill[(as.numeric(chr[subscripts])-1) %% length(gs$fill) + 1]
           A$x <- x
           A$y <- y
           do.call("panel.xyplot", A)
           #draw labels (if requested)
           if(gs$label$show) {
             gt<-gs$label
             names(gt)[which(names(gt)=="text")]<-"labels"
             gt$show<-NULL
             if(is.character(gt$x) | is.character(gt$y)) {
               peak = which.max(y)
               center = mean(range(x))
               if (is.character(gt$x)) {
                 if(gt$x=="peak") {gt$x<-x[peak]}
                 if(gt$x=="center") {gt$x<-center}
               }
               if (is.character(gt$y)) {
                 if(gt$y=="peak") {gt$y<-y[peak]}
               }
             }
             if(is.list(gt$x)) {
               gt$x<-A$getgenpos(gt$x[[1]],gt$x[[2]])
             }
             do.call("panel.text", gt)
           }
         },
         xlab=xlab, ylab=ylab, 
         panel.extra=panel.extra, getgenpos=getGenPos, ...
  );
}

