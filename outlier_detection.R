library(st514)
source("https://imada.sdu.dk/~chdj/ST514/F20/BoxM.r")
data("T1-6")
df <- tbl

#Make a dot plot for each variable
hist(df$V1)
hist(df$V2)
hist(df$V3)
hist(df$V4)
hist(df$V5)


#Outlier detection
colnames(df) <- c(paste("x", 1:5, sep=""), "d2")
dfz <- round(scale(df[1:5]), 1)
colnames(dfz) <- c(paste("z", 1:5, sep = ""))
outliers <- cbind(df[,1:5], 1:98, dfz, tbl[,5])
colnames(outliers)[c(6,12)] <- c("Obs", "d2")
max(outliers$d2)
which(outliers$d2 == max(outliers$d2))

#This function calculates the critical value
chi_crit_val <- qchisq(0.005, 5, lower.tail = F)
chi_crit_val

#This function prints out the outliers by observation number in relation to the treshold of 16.74
which(outliers$d2 > chi_crit_val)
which(outliers$d2 > 16.74)

