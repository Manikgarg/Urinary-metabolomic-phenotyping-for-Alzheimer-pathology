#########################################
# Load metabolomic matrices 
#########################################
drugs <- read.table("/nfs/ma/home/natalja/EMIF/Multimodality_AD_Project/metabolomicQTL/mQTL_data/drugs.txt",sep="\t",header = T,stringsAsFactors = F)
threshold <- 0.01
source("/nfs/ma/home/natalja/EMIF/Multimodality_AD_Project/metabolomicQTL/mQTL_data/scripts/load_data.R")
cvrt = read.table("/nfs/ma/home/natalja/EMIF/Multimodality_AD_Project/metabolomicQTL/mQTL_data/covariates/Covariates_all_vertical_text2.txt",sep="\t",header=T,row.names=1)
cvrt <- cvrt[rownames(cvrt) %in% colnames(expr),]
cvrt <- cvrt[colnames(expr),]
cvrt$Diagnosis <- as.character(cvrt$Diagnosis)
cvrt$Diagnosis[cvrt$Diagnosis=="ADC"] <- "AD"
cvrt$Diagnosis <- as.factor(cvrt$Diagnosis)
#expr
#expr_orig
#expr_orig_qn
#########################################
# Drugs we are intersted in
#########################################
# Paracetamol 2.58_151.0635n
# Diltiazem 5.01_535.1753m/z, 5.51_551.1705m/z

#parecetamol and diltiazem
annotated_metabolites  <- c("2.31_288.0117n")
#diltiazem
#annotated_metabolites  <- c("5.01_535.1753m/z", "5.51_551.1705m/z")
#quinidine
#annotated_metabolites  <- c("3.07_283.1580n")
#acisoga
#annotated_metabolites  <- c("2.62_125.0837n","2.62_92.1133n")

#3-amino-isobutyrate
#annotated_metabolites  <- c("X1.2051","X1.1952")

#bile acids
#annotated_metabolites  <- c("6.13_600.2850m/z","6.13_600.4780m/z","6.13_668.2724m/z","6.13_685.2611m/z","6.13_736.2562m/z","6.30_602.2977m/z","8.01_639.3614n")

#others
#annotated_metabolites  <- c("3.80_258.1087m/z","3.80_126.0658m/z","2.03_112.0501m/z","4.55_190.2275n","3.74_249.0860n")
#methylcytidine 
#annotated_metabolites  <- c("3.80_258.1087m/z","3.80_126.0658m/z","2.03_112.0501m/z")

#annotated_metabolites <-  c("3.07_283.1580n","2.58_151.0635n","5.01_535.1753m/z", "5.51_551.1705m/z","6.30_602.2977m/z","6.13_685.2611m/z", "6.13_600.2850m/z", "6.13_600.4780m/z", "6.13_668.2724m/z", "3.80_258.1087m/z", "4.55_190.2275n", "3.80_126.0658m/z", "8.01_639.3614n", "2.62_125.0837n", "6.13_736.2562m/z", "2.62_92.1133n", "2.03_112.0501m/z", "X5.4085", "3.74_249.0860n", "X1.2051", "X2.873", "X1.7312", "X5.2831", "X5.276", "X3.0968", "X3.0962", "X1.1952", "X3.0979", "X3.0973")

met_norm<- expr[rownames(expr) %in% annotated_metabolites,]
# 2) original data
#met_orig<- expr_orig[rownames(expr_orig) %in% annotated_metabolites,]
# 3) original QN data
#met_orig_qn<- expr_orig_qn[rownames(expr_orig_qn) %in% annotated_metabolites,]

# 3-amino-isobutyrate/3-hydroxibutyrate 1,2


data <- met_norm
#data <- met_orig_qn
#for(i in 1:29){
  #rownames(data[12,])
  #my_data <- data.frame(value=as.numeric(met_norm[i,]),group=cvrt$Diagnosis)
  my_data <- data.frame(value=rowMeans(cbind(as.numeric(data[1,]))),group=cvrt$Diagnosis)
  #grouped <- group_by(my_data,group)
  #summarise(grouped,mean=mean(value),sd=sd(value), nr=n())
  # Compute the analysis of variance
  res.aov <- aov(value ~ group, data = my_data)
  # Summary of the analysis
 annotated_metabolites
  summary(res.aov)
  TukeyHSD(res.aov)
  #TukeyHSD(res.aov)$group[,4]
#}

annotated_metabolites  <- c("2.59_318.0226n")
met_norm<- expr[rownames(expr) %in% annotated_metabolites,]
data <- met_norm
my_data <- data.frame(value=rowMeans(cbind(as.numeric(data[1,]))),group=cvrt$Diagnosis)
res.aov <- aov(value ~ group, data = my_data)
annotated_metabolites
summary(res.aov)
TukeyHSD(res.aov)
######
annotated_metabolites  <- c("5.04_614.2623m/z")
met_norm<- expr[rownames(expr) %in% annotated_metabolites,]
data <- met_norm
my_data <- data.frame(value=rowMeans(cbind(as.numeric(data[1,]))),group=cvrt$Diagnosis)
res.aov <- aov(value ~ group, data = my_data)
annotated_metabolites
summary(res.aov)
TukeyHSD(res.aov)
######
annotated_metabolites  <- c("5.85_614.2613m/z")
met_norm<- expr[rownames(expr) %in% annotated_metabolites,]
data <- met_norm
my_data <- data.frame(value=rowMeans(cbind(as.numeric(data[1,]))),group=cvrt$Diagnosis)
res.aov <- aov(value ~ group, data = my_data)
annotated_metabolites
summary(res.aov)
TukeyHSD(res.aov)
#####
annotated_metabolites  <- c("X5.4085")
met_norm<- expr[rownames(expr) %in% annotated_metabolites,]
data <- met_norm
my_data <- data.frame(value=rowMeans(cbind(as.numeric(data[1,]))),group=cvrt$Diagnosis)
res.aov <- aov(value ~ group, data = my_data)
annotated_metabolites
summary(res.aov)
TukeyHSD(res.aov)
#####
annotated_metabolites  <- c("1.21_229.2724m/z","1.72_229.2717m/z")
met_norm<- expr[rownames(expr) %in% annotated_metabolites,]
data <- met_norm
my_data <- data.frame(value=rowMeans(cbind(as.numeric(data[1,]),as.numeric(data[2,]))),group=cvrt$Diagnosis)
res.aov <- aov(value ~ group, data = my_data)
annotated_metabolites
summary(res.aov)
TukeyHSD(res.aov)
#####
annotated_metabolites  <- c("X1.2051","X1.1952")
met_norm<- expr[rownames(expr) %in% annotated_metabolites,]
data <- met_norm
my_data <- data.frame(value=rowMeans(cbind(as.numeric(data[1,]),as.numeric(data[2,]))),group=cvrt$Diagnosis)
res.aov <- aov(value ~ group, data = my_data)
annotated_metabolites
summary(res.aov)
TukeyHSD(res.aov)
#####
annotated_metabolites  <- c("X2.873")
met_norm<- expr[rownames(expr) %in% annotated_metabolites,]
data <- met_norm
my_data <- data.frame(value=rowMeans(cbind(as.numeric(data[1,]))),group=cvrt$Diagnosis)
res.aov <- aov(value ~ group, data = my_data)
annotated_metabolites
summary(res.aov)
TukeyHSD(res.aov)
#####
annotated_metabolites  <- c("X1.7312")
met_norm<- expr[rownames(expr) %in% annotated_metabolites,]
data <- met_norm
my_data <- data.frame(value=rowMeans(cbind(as.numeric(data[1,]))),group=cvrt$Diagnosis)
res.aov <- aov(value ~ group, data = my_data)
annotated_metabolites
summary(res.aov)
TukeyHSD(res.aov)
#####
annotated_metabolites  <- c("X5.2831","X5.276")
met_norm<- expr[rownames(expr) %in% annotated_metabolites,]
data <- met_norm
my_data <- data.frame(value=rowMeans(cbind(as.numeric(data[1,]),as.numeric(data[2,]))),group=cvrt$Diagnosis)
res.aov <- aov(value ~ group, data = my_data)
annotated_metabolites
summary(res.aov)
TukeyHSD(res.aov)
#####
annotated_metabolites  <- c("4.84_197.0803m/z","1.08_218.1368m/z","2.19_157.0494m/z")
met_norm<- expr[rownames(expr) %in% annotated_metabolites,]
data <- met_norm
my_data <- data.frame(value=rowMeans(cbind(as.numeric(data[1,]),as.numeric(data[2,]),as.numeric(data[3,]))),group=cvrt$Diagnosis)
res.aov <- aov(value ~ group, data = my_data)
annotated_metabolites
summary(res.aov)
TukeyHSD(res.aov)
#####
annotated_metabolites  <- c("6.17_168.0663m/z","6.17_138.0552m/z","6.16_204.0867m/z","6.16_186.0759m/z")
met_norm<- expr[rownames(expr) %in% annotated_metabolites,]
data <- met_norm
my_data <- data.frame(value=rowMeans(cbind(as.numeric(data[1,]),as.numeric(data[2,]),as.numeric(data[3,]),as.numeric(data[4,]))),group=cvrt$Diagnosis)
res.aov <- aov(value ~ group, data = my_data)
annotated_metabolites
summary(res.aov)
TukeyHSD(res.aov)
#####
annotated_metabolites  <- c("4.43_202.1581m/z")
met_norm<- expr[rownames(expr) %in% annotated_metabolites,]
data <- met_norm
my_data <- data.frame(value=rowMeans(cbind(as.numeric(data[1,]))),group=cvrt$Diagnosis)
res.aov <- aov(value ~ group, data = my_data)
annotated_metabolites
summary(res.aov)
TukeyHSD(res.aov)
#####
annotated_metabolites  <- c("1.41_95.0125m/z")
met_norm<- expr[rownames(expr) %in% annotated_metabolites,]
data <- met_norm
my_data <- data.frame(value=rowMeans(cbind(as.numeric(data[1,]))),group=cvrt$Diagnosis)
res.aov <- aov(value ~ group, data = my_data)
annotated_metabolites
summary(res.aov)
TukeyHSD(res.aov)
#####
annotated_metabolites  <- c("2.05_281.1487m/z")
met_norm<- expr[rownames(expr) %in% annotated_metabolites,]
data <- met_norm
my_data <- data.frame(value=rowMeans(cbind(as.numeric(data[1,]))),group=cvrt$Diagnosis)
res.aov <- aov(value ~ group, data = my_data)
annotated_metabolites
summary(res.aov)
TukeyHSD(res.aov)
#####
annotated_metabolites  <- c("6.88_184.1322m/z")
met_norm<- expr[rownames(expr) %in% annotated_metabolites,]
data <- met_norm
my_data <- data.frame(value=rowMeans(cbind(as.numeric(data[1,]))),group=cvrt$Diagnosis)
res.aov <- aov(value ~ group, data = my_data)
annotated_metabolites
summary(res.aov)
TukeyHSD(res.aov)
#####
annotated_metabolites  <- c("1.88_181.0361m/z")
met_norm<- expr[rownames(expr) %in% annotated_metabolites,]
data <- met_norm
my_data <- data.frame(value=rowMeans(cbind(as.numeric(data[1,]))),group=cvrt$Diagnosis)
res.aov <- aov(value ~ group, data = my_data)
annotated_metabolites
summary(res.aov)
TukeyHSD(res.aov)
#####
annotated_metabolites  <- c("4.64_637.1565m/z")
met_norm<- expr[rownames(expr) %in% annotated_metabolites,]
data <- met_norm
my_data <- data.frame(value=rowMeans(cbind(as.numeric(data[1,]))),group=cvrt$Diagnosis)
res.aov <- aov(value ~ group, data = my_data)
annotated_metabolites
summary(res.aov)
TukeyHSD(res.aov)
ag <- aggregate(value ~ group, my_data, function(x) c(mean = mean(x), sd = sd(x)))
ag
#####
annotated_metabolites  <- c("8.01_639.3614n")
met_norm<- expr[rownames(expr) %in% annotated_metabolites,]
data <- met_norm
my_data <- data.frame(value=rowMeans(cbind(as.numeric(data[1,]))),group=cvrt$Diagnosis)
res.aov <- aov(value ~ group, data = my_data)
annotated_metabolites
summary(res.aov)
TukeyHSD(res.aov)
ag <- aggregate(value ~ group, my_data, function(x) c(mean = mean(x), sd = sd(x)))
ag
#####
