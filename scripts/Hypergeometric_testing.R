#List1 has 598 genes and List2 has 5500 genes and the total genes available in the pool from which these two are drawn is of size 23000 (say).
phyper(overlap-1,list1,PopSize-list1,list2,lower.tail = FALSE, log.p = FALSE)

#phyper(88,598,23000-598,5500,lower.tail = FALSE, log.p = FALSE)

m <- 10; n <- 7; k <- 8
x <- 0:(k+1)
rbind(phyper(x, m, n, k), dhyper(x, m, n, k))
all(phyper(x, m, n, k) == cumsum(dhyper(x, m, n, k)))  # FALSE
## but error is very small:
signif(phyper(x, m, n, k) - cumsum(dhyper(x, m, n, k)), digits = 3)

# SNPs 12,105,785
# List1 UPLC-MS 6,047 List2 NMR 876, intersection 838
overlap=838
list1=6047
list2=876
PopSize=12105785
format(phyper(overlap-1,list1,PopSize-list1,list2,lower.tail = FALSE, log.p = FALSE), scientific=TRUE)
phyper(overlap-1,list1,PopSize-list1,list2,lower.tail = FALSE, log.p = TRUE)
#P = exp(-6277.554)  which is smaller than the smallest 
#positive number in double precision, 
#which by the way is available in R as 

# .Machine$double.xmin 
#[1] 2.225074e-308 
