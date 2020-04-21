library(multtest)
# Warranty
# THE AUTHORS MAKE NO WARRANTIES, EXPRESSED OR IMPLIED, REGARDING THE FITNESS OF OUR SOFTWARE FOR ANY PARTICULAR PURPOSE. THE AUTHORS CLAIM NO
# LIABILITY FOR DATA LOSS OR OTHER PROBLEMS CAUSED DIRECTLY OR INDIRECTLY BY THE SOFTWARE. THE USER IS ASSUMING THE ENTIRE RISK AS TO THE SOFTWARES
# QUALITY AND ACCURACY. 
# Parameters
# Note features should be in rows, samples in columns

mswsd <- function(x.data, n)
{
	use_percent = n
	#cat("#### percentage of used features ######","\n")
	#cat("percentage of used features", use_percent,"\n")
	# Set input data
	data <- x.data
	# determination of invariant features
	rem<-floor(nrow(data)-(nrow(data)*use_percent/100))
	# calculate variance
	vari<-apply(data,1,var) #calculate variance for each feature across samples
	# remove features with highest variance
	sel<-order(vari)[1:(nrow(data)-rem)]
	data2<-data[sel,] #data matrix that contains only features with low variance
	rm(vari)
	# Now we want to compute whether enough features have been removed
	# For this we calculate for each spectrum feature wise scaling factors to a baseline spectrum
	# In case that only non-regulated features remain the standard deviation of these scaling factors should
	# be minimal. Note this is done for each spectrum separately since the scaling factors could be
	# drastically different for different spectra. To obtain a value for the complete data set we 
	# calculate the mean over all spectra
	linear.baseline <- apply(abs(data2),1,median) #compute baseline
	# compute for each feature in each spectrum a scaling factor to the baseline
	fac_feat<-matrix(nrow=nrow(data2), ncol=ncol(data2))
	for  (i in 1:ncol(data2)) #samples
	{
   		for (j in 1:nrow(data2)) # features
   		{
     		fac_feat[j,i]<-data2[j,i]/linear.baseline[j]
   	 	}
	}
	# standard deviation of factors witin sample shhould be ideally small
	sd_fac_feat<-numeric(ncol(data2))
	for (i in 1:ncol(data2))
	{
   		sd_fac_feat[i]<-sd(fac_feat[,i])
	}
	# mean standard deviation of factors across samples
	mean_sd_fac_feat=mean(sd_fac_feat)
	#cat("##### standard deviation across feature wise correction factors ####", "\n")
	#cat("mean standard deviation across correction factors", mean_sd_fac_feat,"\n")
	#cat("##### scaling factors ######","\n")
	
	return(mean_sd_fac_feat)
}

# Here comes the main function for resampling of mswsd values
# start with all features than in steps of 5% remove features until 10% of the original features remain
# start by removing the most variable ones first, in total we will have 19 different levels
# to test the stability of the mswsd values we employ a resampling approach with 100 iterations
resamp_mswsd<-function(metabo.data)
{
	s <- seq(10, 100, 5) #prepare steps for feature reduction
	vari <- matrix(ncol = 19, nrow = 100) #to store obtained mswsd values
        no_samp<-ncol(metabo.data) #number of samples
        part_samp<-floor(no_samp*0.66) # use two thirds of samples
        
	for(i in 1:length(s))
	{
		variables <- s[i]
		for(k in 1:100)	# 100 times resampling
		{
			sel <- sample(1:no_samp, part_samp) # randomly select two thirds out of all samples
			vari[k,i] <- mswsd(metabo.data[,sel], variables)
		}
                # normality check on total spectral areas of remaining features
                rem<-floor(nrow(metabo.data)-(nrow(metabo.data)*variables/100))
                # calculate variance
                vari2<-apply(metabo.data,1,var) #calculate variance for each feature across samples
                # remove features with highest variance
                sel2<-order(vari2)[1:(nrow(metabo.data)-rem)]
                metabo.data2<-metabo.data[sel2,] #data matrix that contains only features with low variance
                total <- apply(metabo.data2, 2, sum)
                pp<-shapiro.test(total)$p.value
                cat("pecentage of features used=",variables,"p-value of total areas=",pp,"\n")
                rm(vari2,rem,sel2,metabo.data2)
	}
	# prepare plot
	pdf("resamp_mswsd.pdf", onefile = FALSE)
	box.mswsd <- boxplot(vari, use.cols = TRUE, names = as.character(s), xlab = "percentage of features", ylab = "mswsd", main = "Resampling of mswsd")
	dev.off()
        boxplot(vari, use.cols = TRUE, names = as.character(s), xlab = "percentage of features", ylab = "mswsd", main = "Resampling of mswsd")
}

# From the resampling of mswsd values plot you may have selected the ideal percentage of features used for data normalization
# With this percentage we apply now the actual data normalization

norm_unbal<-function(metabo.data, use_percent,flag)
{       
        metabo.data<-abs(metabo.data)
	# Identify normalization reference features
	cat("#### percentage of used features ######","\n")
	cat("percentage of used features", use_percent,"\n")
	# determination of invariant features
	rem<-floor(nrow(metabo.data)-(nrow(metabo.data)*use_percent/100))
	# calculate variance
	varian<-apply(metabo.data,1,var) #calculate variance for each feature across samples
	# remove features with highest variance
	sel<-order(varian)[1:(nrow(metabo.data)-rem)]
	metabo.data2<-metabo.data[sel,] #data matrix that contains only features with low variance
        
        # perform actual normalization
        # linear baseline normalization based on mean values
        if(flag=="LBME")
        {
	   linear.baseline <- apply(metabo.data2,1,median) #compute baseline
	   baseline.mean <- mean(linear.baseline)
	   sample.means <- apply(metabo.data2,2,mean)
	   linear.scaling <- baseline.mean/sample.means
	   cat("mean lin sca=",mean(linear.scaling),"sd lin sca=",sd(linear.scaling),"rsd lin sca=",sd(linear.scaling)/mean(linear.scaling),"\n")
	   cat("min lin sca=",min(linear.scaling),"max lin sca=",max(linear.scaling),"\n")
	   norm.metabo.data <- t(t(metabo.data)*linear.scaling) 
 	}
        # linear baseline normalization based on median values
        if (flag=="LBMD")
        {
           linear.baseline <- apply(metabo.data2,1,median) #compute baseline
           baseline.median<-median(linear.baseline)
           sample.medians<-apply(metabo.data2,2,median)
           linear.scaling<-baseline.median/sample.medians
	   cat("mean lin sca=",mean(linear.scaling),"sd lin sca=",sd(linear.scaling),"rsd lin sca=",sd(linear.scaling)/mean(linear.scaling),"\n")
           cat("min lin sca=",min(linear.scaling),"max lin sca=",max(linear.scaling),"\n")
	   norm.metabo.data <- t(t(metabo.data)*linear.scaling)
        }
        #VSN
        if(flag=="VSN")
        {
           library(vsn)
           vsn.model<-vsn2(metabo.data2)
           norm.metabo.data<-predict(vsn.model,(metabo.data))
        }
        #PQN
        if(flag=="PQN")
        {
           reference<-apply(metabo.data2,1,median)
	   quotient<-metabo.data2/reference
	   quotient.median<-apply(quotient,2,median)
	   norm.metabo.data<-t(t(metabo.data)/quotient.median)
        }
        return(norm.metabo.data)
	rm(varian,sel,rem)
}

unbal_reg<-function(metabo.data)
{
   metabo.data<-abs(metabo.data)
   {
      cat("if p-value smaller than 0.05 total spectral areas are not normally distributed-> unbalanced regulation","\n")
      total <- apply(metabo.data, 2, sum)
      shapiro.test(total)
   }
}
