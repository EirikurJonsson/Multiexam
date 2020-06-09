library(st514)
library(MASS)
library(tidyverse)
library(ggplot2)
data("T1-6")
ggplot(tbl, mapping = aes(x = V3, y = V5, color = V6))+
  geom_point(size = 3)


#hvad christian har gjort i hans løsninger
# Exercise 6.37
#data("T6-9")
#source("https://imada.sdu.dk/~chdj/ST514/F20/BoxM.r")
#tbl
#BoxM(log.tbl[1:3], tbl[4], Ftest = FALSE)



#install.packages("biotools")
library(biotools)
?boxM
boxM(tbl[1:5], tbl[6])

