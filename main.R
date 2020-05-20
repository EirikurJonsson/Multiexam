library(st514)
library(tidyverse)
library(ggplot2)
library(GGally)
library(corrplot)
library(ggcorrplot)
library(ggpubr)
library(lattice)
data("T1-6")

#Renaming
colnames(tbl) <- c("Age", "V2", "V3", "V4", "V5","target")
tbl$target <- as.factor(tbl$target)

# Non-multiple-sclerosis (NMS) group data
nms <- subset(tbl, target == 0)
# Multiple-sclerosis (MS) group data
ms <- subset(tbl, target == 1)

###################DISCRIPTIVE STATISTICS#########################

descriptive <- function(x){
  "
  This function takes in one input parameter
  and that is a dataframe or matrix
  
  The output of this function will be
  
  1. Mean for each column in the dataset
  2. Covariance matrix
  3. Co/Variance Matrix
  4. Variance matrix
  5. Inverse of the covariance matrix
  6. Determinant of the covariance matrix
  "
  xmean <- colMeans(x)
  xcov <- cov(x)
  xcovVar <- cov(x) * (nrow(x) -1) / nrow(x)
  xvar <- var(x)
  xcovInv <- solve(xcov)
  xdet <- det(xcov)

  # Debug block
  
  # print("Mean Vector")
  # print(xmean)
  # print("-----------------------------------------------------")
  # print("Covariance")
  # print(xcov)
  # print("-----------------------------------------------------")
  # print("Xco/Variance")
  # print(xcovVar)
  # print("-----------------------------------------------------")
  # print("Variance")
  # print(xvar)
  # print("-----------------------------------------------------")
  # print("Covariance inverse")
  # print(xcovInv)
  # print("-----------------------------------------------------")
  # print("Determinant")
  # print(xdet)
  mylist = list("meanVector" = xmean, "Covariance" = xcov, "Co-Variance" = xcovVar, "Variance" = xvar, "CovInverse" = xcovInv, "Determinant" = xdet)
  return(mylist)
}


confRegionE <- function(x, f){
  "
  This function calculates the confidence region ellips.
  It takes in two imput parameters, dataframe/matrix and 
  a signifigance level. 
  
  For the output we get:
  1. qqnorm plot for each column
  2. qqline plot for each column
  3. A single qqplot using qchisq function which is the
     Mahalanobis distance.
  "
  cVar <- cov(x)
  mahl <- mahalanobis(x, colMeans(x), cVar)
  r <- eigen(cov(x))
  n <- nrow(x)
  p <- ncol(x)
  fstat <- qf(f, p, n-p)
  region <- rep(0, length(r$values))
  for (i in 1:p){
    region[i] <- 2*sqrt(r$values[i] * p * (n-1) * fstat / (n*(n-p)))
  }
  par(mfrow = c(3,2))
  for (i in 1:p) {
    qqnorm(x[,i], main = colnames(x)[i])
    qqline(x[,i])
  }
  qqplot(qchisq(ppoints(n, a = 0.5), df = p), mahl,
         ylab = "Mahalanobis distances",
         xlab = bquote("Quantiles of " ~ chi[.(p)]^2),
         main = bquote("Q-Q plot of Mahalanobis" * ~ D^2 * 
                         " vs. quantiles of" * ~ chi[.(p)]^2))

  
}
confRegionE(ms[1:5], 0.95)

################################################################################################################
#The Bonferroni method is an alternative method for multiple comparisons (developed from a probability inequality).
#It is used for type 1 errors (aka false positive etc. When doing multiple comparisons, the chance of committing type 1 errors increases).
#To get the Bonferroni adjusted p-value, we need to divide the original p-value with the number of analyses on the dependent variables (aka, the number n).
#N.B.: When doing the Bonferroni, you need ALL the decimals, to get the correct result!
#Example: 0.05/9 = 0.006 (where 0.05 was the original p-value and 9 was the number of n). 

bon <- function(x, y){
  # Our actual p value needs to be lower that the p-value returned from this function
  # We enter 1-x for our signifigance value and y is number of hypothesis
  
  # Example 0.05/2 where signifigance level of 0.05 or with 95% siginifigance.
  print("Bonfferonni adjusted p-value")
  return(x/y)
}

#################################################################################################################
"
Calculating the p-values of the F-statistics we get from the 
Tsquared function that is here down below.
"
print(pf(62.9, df1 = ncol(ms), df2 = nrow(ms)-ncol(ms), lower.tail = FALSE))
print(pf(7642.242, df1 = ncol(nms), df2 = nrow(nms) - ncol(nms), lower.tail = FALSE))

#################################################################################################################
#The two sample Hotellings T^2 is the multivariate extension of the common two group students t-test.
#In a t-test, differences in the mean response between two populations are studied. T^2 is used when the number of response variables are two or more. 
#If the observed statistical distance T^2 is too large (aka, if sample mean is "too far" from mu0), the H0:mu=mu0 is rejected.  

tsqared <- function(x,y,f){
  "
  Hotellings Tsquared function has 3 input parameters
  1. X = dataset
  2. Y = mean vector to test if exists as a mean vector in the dataset
  3. F = signifgance level for our F-statistics
  
  This function then returns:
  
  1. Results of the actual Hotellings Tsquared test
  2. F-statistic
  "
  n <- nrow(x)
  xmean <- colMeans(x)
  xcov <- cov(x)
  xcovin <- solve(xcov)
  xyvar <- xmean - y
  p1 <- xyvar %*% xcovin
  result <- round(p1%*%(n*xyvar),3)
  fstat <- qf(f, df1 = ncol(x), df2 = nrow(x)-ncol(x))
  if(result > fstat){
    print("We REJECT the Null Hypothesis that y is a possible mean vector in x")
  }else{
    print("We ACCEPT the Null Hypothesis and y is a possible mean vector in x")
  }
  return(c(result, fstat))
}


# Probability density function
pDensity <- function(x, y){
  # Probability density function
  # X = dataset
  # Y = mean vector
  p = ncol(x)
  xcov = cov(x)
  xdet = det(xcov)
  xcovInv = solve(xcov)
  mu = colMeans(x)
  result = (1/((2*pi)^(p/2)*sqrt(xdet)) %*% exp(-(y-mu)%*%(xcovInv)%*%(y-mu)/2))
  return(result)
  
}

##################### Plot section ##############################


# Non MS
boxplot(tbl$V2[tbl$target == 0])
boxplot(tbl$V3[tbl$target == 0])
boxplot(tbl$V4[tbl$target == 0])
boxplot(tbl$V5[tbl$target == 0])

# MS

boxplot(tbl$V2[tbl$target == 1])
boxplot(tbl$V3[tbl$target == 1])
boxplot(tbl$V4[tbl$target == 1])
boxplot(tbl$V5[tbl$target == 1])

# Colored correlation plot
ggpairs(data = tbl, columns = 1:5, mapping = aes(color=target))

marginal.dot.plot(nms$V2, nms$V4, pch = 1, col = "pink")

# QQ plots for norm and line
qqnorm(tbl[,1])
qqline(tbl[,1])
##QQ2 is not normally distributed 
qqnorm(tbl[,2])
qqline(tbl[,2])
##QQ3 not normally distributed
qqnorm(tbl[,3])
qqline(tbl[,3])
#QQ4 not normal
qqnorm(tbl[,4])
qqline(tbl[,4])
##QQ5 not normal
qqnorm(tbl[,5])
qqline(tbl[,5])
