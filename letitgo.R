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

#DISCRIPTIVE STATISTICS

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

####################################################################################################################
#DISCRIPTIVE STATISTICS FOR NMS GROUP:
descriptive(nms[1:5])

#DISCRIPTIVE STATISTICS FOR MS GROUP:
descriptive(ms[1:5])

############MEAN VECTORS
#Age
#According to the mean vectors for NMS and MS, respectively, it is evident that the mean age for the NMS group is lower than for the MS group.

#Mean V2 vs V4
#The NMS group generally has a lower response time to stimuli than the MS group.
#Mean V3 vs V5
#The difference between the left- and the right eye is suggestively larger in the MS group compared to the NMS group.

############COVARIANCE MATRIX
#If two quantities have a positive covariance, they increase/decrease together.
#If both variables tend to increase or decrease together, the coefficient is positive. If one variable tends to increase as the other decreases, the coefficient is negative.
#Finally, the differences between L and R are significantly higher in the MS group compared to the NMS group.

############VARIANCE MATRIX
#Variance measures how far a set of data is spread out. A variance of zero indicates that all the data values are identical. All non-zero variances are positive. A small variance indicates that the data points tend to be very close to the mean and to each other.

#For V2 and V4
#For the NMS group, the variance between V2 and V4 is low compared to the MS group. 

#For V3 and V5
#For the NMS group the density is high between V3 and V5 with a value of 0.50. In contrast the variance for V3 and V5 is high which indicates a widespread. 

####################################################################################################################

# qq norm calculates the univariate normlity of the different variables
# if the values lie along a straight line, the data follows a normal distribution with some mean and variance
#In these cases we see, that near the two extemes, the points deviate from the straight line, indicating that the assumption of normality may not be valid
#deviations from the straight line might indicate: heavier tails, skewness, outliers, and/or clustered data
#QQ1 is not normally distributes since the data deviates from the straight line
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

########################################################################################

#although univariate normality does not imply multivariate normality if any single variable fails to follow normality we cannot have joint multivariate normality
#For multivariate normality, both p-values of skewness and kurtosis statistics should be greater than 0.05

# WE NEED TO APPLY THIS FUNCTION TO OUR DATA IN ORDER TO DETERMINE, IF THE DATA IS MULTIVARIATE NORMALLY DISTRIBUTED
# THIS IS CHRISTIANS CODE FROM EXERCISE 4.26, I HAVE TRIED TO APPLY OUR DATA - IT COULD BE AWESOME THOUGH TO GET A QQLINE - I MANAGED TO GET ONE, BUT I WAS NOT ENTIRELY SURE, IF IT WAS CORRECT - SO I DELETED IT.
# Exercise 4.26
v <- as.matrix(tbl[,1:5])
v.mean <- as.vector(colMeans(v))
S <- cov(v)
n <- nrow(v)
D2 <- rep(0, n)
for (i in 1:n) {
  D2[i] <- t(v[i,] - v.mean) %*% solve(S) %*% (v[i,] - v.mean)
}
round(D2, 4)
which(D2 < qchisq(0.5, 2))
p <- ncol(v)
D2 <- mahalanobis(v, colMeans(v), S)
qqplot(qchisq(ppoints(n, a = 0.5), df = p), D2,
       ylab = "Mahalanobis distances",
       xlab = bquote("Quantiles of " ~ chi[.(p)]^2),
       main = bquote("Q-Q plot of Mahalanobis" * ~ D^2 * 
                       " vs. quantiles of" * ~ chi[.(p)]^2))

#####################################################################################################################

#A confidence region is a multi-dimensional generalization of a confidence interval (an interval for an individual variable).
#It is a set of point in a n-dimensional space, often represented as an ellipsoid around a point (which is an estimated solution to a problem). 
#To determine whether any mu0 lies within the confidence region (aka, is a possible value for mu),we need to compute the generalized squared distance and
#compare it with the F distribution. If the squared distance is larger than the F distribution, then mu0 is NOT in the confidence region.
#Since this is analogous to testing H0: mu = mu0 vs. H1: mu IS NOT EQUAL TO mu0, we see that the confidence region consists of all mu0 vectors for which the T^2 test
#would NOT reject H0 in favor of H1 at significance level, alpha. 


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

#########################################################################################################


# Plot section

## Boxplot section

### Non MS
boxplot(tbl$V2[tbl$target == 0])
boxplot(tbl$V3[tbl$target == 0])
boxplot(tbl$V4[tbl$target == 0])
boxplot(tbl$V5[tbl$target == 0])

### MS

boxplot(tbl$V2[tbl$target == 1])
boxplot(tbl$V3[tbl$target == 1])
boxplot(tbl$V4[tbl$target == 1])
boxplot(tbl$V5[tbl$target == 1])

## Colored correlation plot
ggpairs(data = tbl, columns = 1:5, mapping = aes(color=target))

marginal.dot.plot(nms$V2, nms$V4, pch = 1, col = "pink")


a = descriptive(nms[1:5])
b = descriptive(ms[1:5])


###############################################################################################################

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


