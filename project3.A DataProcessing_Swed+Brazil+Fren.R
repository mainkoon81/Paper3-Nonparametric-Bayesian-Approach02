
library(MCMCpack)
library(mvnfast)
library(tidyverse)

set.seed(1)


library(car) # Companion to Applied Regression": to Recodes a numeric, character vector as per specifications.
library(varhandle) # to use "unfactor()"

library(mvtnorm)
#library(mgcv)   # for GAM
library(splines) # for MARS
library(statmod) 
  


#---------------------------- DATA PROCESSING ---------------------------------#


##################
# with Brazil Data:+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##################

train_df <- read.csv("C:/Users/kimm4/Desktop/WORKSPACE/2nd 3rd paper/DATASET/Brazil_c.train.csv",
                      header=T, 
                      na.strings=c("."), 
                      stringsAsFactors=F)
test_df <- read.csv("C:/Users/kimm4/Desktop/WORKSPACE/2nd 3rd paper/DATASET/Brazil_c.test.csv",
                     header=T, 
                     na.strings=c("."), 
                     stringsAsFactors=F)

#> train set -------------------------------------------------------------------
head(train_df)
summary(train_df)
# subset(train_df1, subset=SumInsAvg==0) #.... almost..same as "ExposTotal==0", but it is a subset..
# sum(train_df1$SumInsAvg==0) # there are 2799 obv with the 0 inflation
# train_df1 <- subset(train_df1, subset=SumInsAvg>0) 
# train_df1$Affiliation = as.numeric(train_df1$Affiliation) # z (binary) should be numeric
# train_df1 <- drop_na(train_df1) # we r not interested in NA

train_df$Area <- factor(train_df$Area)
train_df$StateAb <- factor(train_df$StateAb)
summary(train_df)


# ---------------------------------------------------------------------------- #
cor(train_df$AggClaim, train_df$ExposTotal, method="pearson") # cor: 0.54

cor(train_df$AggClaim, train_df$SumInsAvg, method="pearson")  # cor: 0.10
cor(train_df$AggClaim, train_df$PremTotal, method="pearson")  # cor: 0.59
# ---------------------------------------------------------------------------- #


quantile(train_df$PremTotal, c(0.025, 0.975))
#      2.5%      97.5% 
#   226.022  77627.929 

quantile(train_df$SumInsAvg, c(0.025, 0.975))
#      2.5%       97.5% 
#   8081.121  100118.792


plot(train_df$ExposTotal, train_df$PremTotal) # seemingly, large exposure always gives large premium..
plot( log(train_df$ExposTotal), log(train_df$PremTotal) )
plot(train_df$ExposTotal, train_df$AggClaim) 


# premium can be the possible hidden variable??
sum( train_df$PremTotal<=225 ) # 1490 obv
sum( train_df$PremTotal<77600 & train_df$PremTotal>225 ) # 57069 obv
sum( train_df$PremTotal>77600 ) # 1504 obv

sum( test_df$PremTotal<=225 ) # 61 obv
sum( test_df$PremTotal<77600 & test_df$PremTotal>225 ) # 2325 obv
sum( test_df$PremTotal>77600 ) # 63 obv

# or

# SumInsAvg can be the possible hidden variable??
sum( train_df$SumInsAvg<=8000 ) # 1476 obv
sum( train_df$SumInsAvg<98000 & train_df$PremTotal>8000 ) # 20859 obv
sum( train_df$SumInsAvg>98000 ) # 1613 obv

sum( test_df$SumInsAvg<=8000 ) # 48 obv
sum( test_df$SumInsAvg<98000 & test_df$PremTotal>8000 ) # 837 obv
sum( test_df$SumInsAvg>98000 ) # 59 obv


################################################################################
### Now, let's add some errors

n.breaks = sqrt( nrow(train_df) ) #******Rule of thumb
hist(train_df$AggClaim, breaks=n.breaks)

hist(train_df$ExposTotal, breaks=n.breaks)# , xlim=c(0,100))
quantile(train_df$ExposTotal, c(0.025, 0.975))
#  2.5%  97.5% 
#0.240 78.759 

hist(train_df$PremTotal, breaks=n.breaks)# , xlim=c(0,100))
quantile(train_df$PremTotal, c(0.025, 0.975))
#  2.5%  97.5% 
#226.022 77627.929 






#-------------------------------------------------------------------------------
# 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1%
############################# Error For train ##################################
#table(train_df$Affiliation)
# #    0       1 
# # 2599   57464   
#subset(train_df, subset=Affiliation==0) # prone to error
# train_df[which(train_df$ExposTotal>100), ] # to face the truth, lots of outliers in Exposure | affiliation="1"

quantile(train_df$PremTotal, c(0.025, 0.975))
#  2.5%  97.5% 
#226.022 77627.929 

train.n <- length(train_df$AggClaim)
err <- numeric(train.n)

# latent unknown factor...possible? [PremTotal]
# uuu = 0.1
# uuu = 0.4

set.seed(1)

for (i in 1:train.n) {
  if(train_df$Affiliation[i]==0){
    # From unknown factor
    uuu <- ifelse( train_df$PremTotal[i]<77600 & train_df$PremTotal[i]>225, 0.1, 0.4 )

    # From clean covariate
    delta_sq <- ifelse( train_df$ExposTotal[i]<=0.5, 1, 
                        ifelse(train_df$ExposTotal[i]>35, 25, 3) )
    
    err[i] <- rnorm( n=1, mean=0, sd=sqrt( uuu*delta_sq ) )
  }
  else{
    err[i] <- 0
  }
}
hist(err,  breaks=n.breaks)

train_df$Error <- err

head(train_df, 10)

train_df$ln_Expos <- log(train_df$ExposTotal)
train_df$ln_Expos_err <- train_df$ln_Expos + train_df$Error 

hist(train_df$ln_Expos_err, breaks=n.breaks)

# check 
plot(train_df$ln_Expos, train_df$Error, cex=0.5)

plot(train_df$ln_Expos, train_df$ln_Expos_err, cex=0.5)


# 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1%
################################## For test ####################################
test.n <- length(test_df$SumInsAvg)
test.n <- length(test_df$ExposTotal)
err <- numeric(test.n)

# uuu = 0.1
# uuu = 0.4

set.seed(1)

for (i in 1:test.n) {
  if(test_df$Affiliation[i]==0){
    
    uuu <- ifelse( test_df$PremTotal[i]<77600 & test_df$PremTotal[i]>225, 0.1, 0.4 )

    delta_sq <- ifelse( test_df$ExposTotal[i]<=0.5, 1, 
                        ifelse(test_df$ExposTotal[i]>35, 100, 3) )
    
    err[i] <- rnorm( n=1, mean=0, sd=sqrt( uuu*delta_sq ) )
  }
  else{
    err[i] <- 0
  }
}
hist(err,  breaks=n.breaks)


test_df$Error <- err

head(test_df)
summary(test_df)

test_df$ln_Expos <- log(test_df$ExposTotal)
test_df$ln_Expos_err <- test_df$ln_Expos + test_df$Error 

hist(test_df$ln_Expos_err, breaks=n.breaks)

# check 
plot(test_df$ln_Expos, test_df$Error, cex=0.5)

plot(test_df$ln_Expos, test_df$ln_Expos_err, cex=0.5)





### Overall Error percentage ( 1%, 10%, 50% ? )

head(train_df)
train.overall_E.bra = sum(abs(train_df$Error)) /  sum(abs(train_df$ln_Expos)); train.overall_E.bra
# 1 %
test.overall_E.bra = sum(abs(test_df$Error)) /  sum(abs(test_df$ln_Expos)); test.overall_E.bra
# 0.8 % ------------------------------------------------------------------------
# Therefore.....................................................................
# SAVE....----------------------------------------------------------------------
# write.csv(train_df, "C:/Users/kimm4/Desktop/WORKSPACE/2nd 3rd paper/DATASET/Brazil.paper3-1.train.csv",
#           row.names=FALSE)
# write.csv(test_df, "C:/Users/kimm4/Desktop/WORKSPACE/2nd 3rd paper/DATASET/Brazil.paper3-1.test.csv",
#           row.names=FALSE)




#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# 10% 10% 10% 10% 10% 10% 10% 10% 10% 10% 10% 10% 10% 10% 10% 10% 10% 10% 10% 10%
############################# Error For train ##################################
quantile(train_df$PremTotal, c(0.025, 0.975))
#  2.5%  97.5% 
#226.022 77627.929 

table(test_df$Affiliation)
#   0*    1 
# 104  2345 



train.n <- length(train_df$AggClaim)
err <- numeric(train.n)

for (i in 1:train.n) {
  if(train_df$Affiliation[i]==0){
    # From unknown factor
    uuu <- ifelse( train_df$PremTotal[i]<77600 & train_df$PremTotal[i]>225, 5, 0.5 )
    
    # From clean covariate
    delta_sq <- ifelse( train_df$ExposTotal[i]<=0.5, 10, 
                        ifelse(train_df$ExposTotal[i]>35, 25, 30) )
    
    err[i] <- rnorm( n=1, mean=0, sd=sqrt( uuu*delta_sq ) )
  }
  else{
    err[i] <- 0
  }
}
hist(err,  breaks=n.breaks)

train_df$Error <- err

head(train_df, 10)

train_df$ln_Expos <- log(train_df$ExposTotal)
train_df$ln_Expos_err <- train_df$ln_Expos + train_df$Error 


# 10% 10% 10% 10% 10% 10% 10% 10% 10% 10% 10% 10% 10% 10% 10% 10% 10% 10% 10% 10%
################################## For test ####################################
test.n <- length(test_df$SumInsAvg)
test.n <- length(test_df$ExposTotal)
err <- numeric(test.n)

# uuu = 0.1
# uuu = 0.4

set.seed(1)

for (i in 1:test.n) {
  if(test_df$Affiliation[i]==0){
    
    uuu <- ifelse( test_df$PremTotal[i]<77600 & test_df$PremTotal[i]>225, 5, 0.5 )
    
    delta_sq <- ifelse( test_df$ExposTotal[i]<=0.5, 10, 
                        ifelse(test_df$ExposTotal[i]>35, 25, 30) )
    
    err[i] <- rnorm( n=1, mean=0, sd=sqrt( uuu*delta_sq ) )
  }
  else{
    err[i] <- 0
  }
}
hist(err,  breaks=n.breaks)


test_df$Error <- err

head(test_df)
summary(test_df)

test_df$ln_Expos <- log(test_df$ExposTotal)
test_df$ln_Expos_err <- test_df$ln_Expos + test_df$Error 


### Overall Error percentage ( 1%, 10%, 50% ? )

head(train_df)
train.overall_E.bra = sum(abs(train_df$Error)) /  sum(abs(train_df$ln_Expos)); train.overall_E.bra
# 17 %
test.overall_E.bra = sum(abs(test_df$Error)) /  sum(abs(test_df$ln_Expos)); test.overall_E.bra
# 15 % ------------------------------------------------------------------------
# Therefore.....................................................................
# SAVE....----------------------------------------------------------------------
# write.csv(train_df, "C:/Users/kimm4/Desktop/WORKSPACE/2nd 3rd paper/DATASET/Brazil.paper3-10.train.csv",
#           row.names=FALSE)
# write.csv(test_df, "C:/Users/kimm4/Desktop/WORKSPACE/2nd 3rd paper/DATASET/Brazil.paper3-10.test.csv",
#           row.names=FALSE)





#-------------------------------------------------------------------------------
# 25% 25% 25% 25% 25% 25% 25% 25% 25% 25% 25% 25% 25% 25% 25% 25% 25% 25% 25% 25%
############################# Error For train ##################################
#table(train_df$Affiliation)
# #    0       1 
# # 2599   57464   
#subset(train_df, subset=Affiliation==0) # prone to error
# train_df[which(train_df$ExposTotal>100), ] # to face the truth, lots of outliers in Exposure | affiliation="1"

quantile(train_df$PremTotal, c(0.025, 0.975))
#  2.5%  97.5% 
#226.022 77627.929 

train.n <- length(train_df$AggClaim)
err <- numeric(train.n)

# latent unknown factor...possible? [PremTotal]
# uuu = 0.1
# uuu = 0.4

set.seed(1)

for (i in 1:train.n) {
  if(train_df$Affiliation[i]==0){
    # From unknown factor
    uuu <- ifelse( train_df$PremTotal[i]<77600 & train_df$PremTotal[i]>225, 2, 10 )
    
    # From clean covariate
    delta_sq <- ifelse( train_df$ExposTotal[i]<=0.5, 80, 
                        ifelse(train_df$ExposTotal[i]>35, 200, 15) )
    
    err[i] <- rnorm( n=1, mean=0, sd=sqrt( uuu*delta_sq ) )
  }
  else{
    err[i] <- 0
  }
}


hist(err,  breaks=n.breaks)

train_df$Error <- err

head(train_df, 10)

train_df$ln_Expos <- log(train_df$ExposTotal)
train_df$ln_Expos_err <- train_df$ln_Expos + train_df$Error 

hist(train_df$ln_Expos_err, breaks=n.breaks)

# check 
plot(train_df$ln_Expos, train_df$Error, cex=0.5)

plot(train_df$ln_Expos, train_df$ln_Expos_err, cex=0.5)


# 25% 25% 25% 25% 25% 25% 25% 25% 25% 25% 25% 25% 25% 25% 25% 25% 25% 25% 25% 25%
################################## For test ####################################
test.n <- length(test_df$SumInsAvg)
test.n <- length(test_df$ExposTotal)
err <- numeric(test.n)

# uuu = 0.1
# uuu = 0.4

set.seed(1)

for (i in 1:test.n) {
  if(test_df$Affiliation[i]==0){
    
    uuu <- ifelse( test_df$PremTotal[i]<77600 & test_df$PremTotal[i]>225, 3, 15 )
    
    delta_sq <- ifelse( test_df$ExposTotal[i]<=0.5, 80, 
                        ifelse(test_df$ExposTotal[i]>35, 200, 15) )
    
    err[i] <- rnorm( n=1, mean=0, sd=sqrt( uuu*delta_sq ) )
  }
  else{
    err[i] <- 0
  }
}
hist(err,  breaks=n.breaks)


test_df$Error <- err

head(test_df)
summary(test_df)

test_df$ln_Expos <- log(test_df$ExposTotal)
test_df$ln_Expos_err <- test_df$ln_Expos + test_df$Error 

hist(test_df$ln_Expos_err, breaks=n.breaks)

# check 
plot(test_df$ln_Expos, test_df$Error, cex=0.5)

plot(test_df$ln_Expos, test_df$ln_Expos_err, cex=0.5)





### Overall Error percentage ( 1%, 10%, 50% ? )

head(train_df)
train.overall_E.bra = sum(abs(train_df$Error)) /  sum(abs(train_df$ln_Expos)); train.overall_E.bra
# 25 %
test.overall_E.bra = sum(abs(test_df$Error)) /  sum(abs(test_df$ln_Expos)); test.overall_E.bra
# 25 % ------------------------------------------------------------------------
# Therefore.....................................................................
# SAVE....----------------------------------------------------------------------
# write.csv(train_df, "C:/Users/kimm4/Desktop/WORKSPACE/2nd 3rd paper/DATASET/Brazil.paper3-25.train.csv",
#           row.names=FALSE)
# write.csv(test_df, "C:/Users/kimm4/Desktop/WORKSPACE/2nd 3rd paper/DATASET/Brazil.paper3-25.test.csv",
#           row.names=FALSE)






#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# 40% 40% 40% 40% 40% 40% 40% 40% 40% 40% 40% 40% 40% 40% 40% 40% 40% 40% 40% 40%
############################# Error For train ##################################
quantile(train_df$PremTotal, c(0.025, 0.975))
#  2.5%  97.5% 
#226.022 77627.929 

table(test_df$Affiliation)
#   0*    1 
# 104  2345 



train.n <- length(train_df$AggClaim)
err <- numeric(train.n)

for (i in 1:train.n) {
  if(train_df$Affiliation[i]==0){
    # From unknown factor
    uuu <- ifelse( train_df$PremTotal[i]<77600 & train_df$PremTotal[i]>225, 30, 8 )
    
    # From clean covariate
    delta_sq <- ifelse( train_df$ExposTotal[i]<=0.5, 10, 
                        ifelse(train_df$ExposTotal[i]>35, 25, 30) )
    
    err[i] <- rnorm( n=1, mean=0, sd=sqrt( uuu*delta_sq ) )
  }
  else{
    err[i] <- 0
  }
}
hist(err,  breaks=n.breaks)

train_df$Error <- err

head(train_df, 10)

train_df$ln_Expos <- log(train_df$ExposTotal)
train_df$ln_Expos_err <- train_df$ln_Expos + train_df$Error 


# 40% 40% 40% 40% 40% 40% 40% 40% 40% 40% 40% 40% 40% 40% 40% 40% 40% 40% 40% 40%
################################## For test ####################################
test.n <- length(test_df$SumInsAvg)
test.n <- length(test_df$ExposTotal)
err <- numeric(test.n)

# uuu = 0.1
# uuu = 0.4

set.seed(17)

for (i in 1:test.n) {
  if(test_df$Affiliation[i]==0){
    
    uuu <- ifelse( test_df$PremTotal[i]<77600 & test_df$PremTotal[i]>225, 30, 8 )
    
    delta_sq <- ifelse( test_df$ExposTotal[i]<=0.5, 10, 
                        ifelse(test_df$ExposTotal[i]>35, 25, 30) )
    
    err[i] <- rnorm( n=1, mean=0, sd=sqrt( uuu*delta_sq ) )
  }
  else{
    err[i] <- 0
  }
}


test_df$Error <- err

head(test_df)
summary(test_df)

test_df$ln_Expos <- log(test_df$ExposTotal)
test_df$ln_Expos_err <- test_df$ln_Expos + test_df$Error 


### Overall Error percentage ( 1%, 10%, 50% ? )

head(train_df)
train.overall_E.bra = sum(abs(train_df$Error)) /  sum(abs(train_df$ln_Expos)); train.overall_E.bra
# 44 %
test.overall_E.bra = sum(abs(test_df$Error)) /  sum(abs(test_df$ln_Expos)); test.overall_E.bra
# 43 % ------------------------------------------------------------------------
# Therefore.....................................................................
# SAVE....----------------------------------------------------------------------
# write.csv(train_df, "C:/Users/kimm4/Desktop/WORKSPACE/2nd 3rd paper/DATASET/Brazil.paper3-40.train.csv",
#           row.names=FALSE)
# write.csv(test_df, "C:/Users/kimm4/Desktop/WORKSPACE/2nd 3rd paper/DATASET/Brazil.paper3-40.test.csv",
#           row.names=FALSE)




#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

##################
# with Sweden Data: ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##################
train.swed <- read.csv("C:/Users/kimm4/Desktop/WORKSPACE/2nd 3rd paper/DATASET/BSwed.train.csv",
                      header=T, 
                      na.strings=c("."), 
                      stringsAsFactors=F)
test.swed <- read.csv("C:/Users/kimm4/Desktop/WORKSPACE/2nd 3rd paper/DATASET/BSwed.test.csv",
                     header=T, 
                     na.strings=c("."), 
                     stringsAsFactors=F)

head(train.swed)
summary(train.swed)
table(train.swed$Experience)
#  0*   1 
#658  894
table(train.swed$Risk)
#  1   2   3   4   5 (25,000 km)
#327 339 328 293 265 
cor(train.swed$TotalLoss, train.swed$SumInsAvg)

quantile(train.swed$SumInsAvg, c(0.025, 0.975))
#      2.5%      97.5% 
#   7.53875 9811.19775 


#### Add Error!

#-------------------------------------------------------------------------------
# 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1%
## > For train 
train.n <- length(train.swed$TotalLoss)
err <- numeric(train.n)

# latent unknown factor...possible? [Risk]
# uuu = 1
# uuu = 10
# uuu = 20

set.seed(17)
for (i in 1:train.n) {
  if(train.swed$Experience[i]==0){
    
    uuu <- ifelse( train.swed$Risk[i]>=4, 0.1, 
                   ifelse(train.swed$Risk[i]<=2, 1, 0.6) )
    delta_sq <- ifelse( train.swed$SumInsAvg[i]<9800 & train.swed$SumInsAvg[i]>7.5, 0.1, 0.2 )
    
    err[i] <- rnorm( n=1, mean=0, sd=sqrt( uuu*delta_sq ) )
  }
  else{
    err[i] <- 0
  }
}
n.breaks = sqrt( nrow(train.swed) ) #******Rule of thumb
hist(err,  breaks=n.breaks)

train.swed$Error <- err

head(train.swed)
summary(train.swed)

train.swed$ln_insured <- log(train.swed$SumInsAvg)
train.swed$ln_insured_err <- train.swed$ln_insured + train.swed$Error 

# 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1%
## > For test
test.n <- length(test.swed$TotalLoss)
err <- numeric(test.n)

set.seed(13)
for (i in 1:test.n) {
  if(train.swed$Experience[i]==0){
    
    uuu <- ifelse( test.swed$Risk[i]>=4, 0.1, 
                   ifelse(test.swed$Risk[i]<=2, 1, 0.6) )
    delta_sq <- ifelse( test.swed$SumInsAvg[i]<9800 & test.swed$SumInsAvg[i]>7.5, 0.1, 0.2 )
    
    err[i] <- rnorm( n=1, mean=0, sd=sqrt( uuu*delta_sq ) )
  }
  else{
    err[i] <- 0
  }
}

n.breaks = sqrt( nrow(test.swed) ) #******Rule of thumb
hist(err,  breaks=n.breaks)

test.swed$Error <- err

head(test.swed)
summary(test.swed)

test.swed$ln_insured <- log(test.swed$SumInsAvg)
test.swed$ln_insured_err <- test.swed$ln_insured + test.swed$Error 

head(train.swed)
train.overall_E.swe = sum(abs(train.swed$Error)) /  sum(abs(train.swed$ln_insured )); train.overall_E.swe
# 1.5 %
test.overall_E.swe = sum(abs(test.swed$Error)) /  sum(abs(test.swed$ln_insured )); test.overall_E.swe
# 1.7 % --------------------------------------------------------------------------
# Therefore.....................................................................
# SAVE....----------------------------------------------------------------------
# write.csv(train.swed, "C:/Users/kimm4/Desktop/WORKSPACE/2nd 3rd paper/DATASET/BSwed.paper3-1.train.csv",
#           row.names=FALSE)
# write.csv(test.swed, "C:/Users/kimm4/Desktop/WORKSPACE/2nd 3rd paper/DATASET/BSwed.paper3-1.test.csv",
#           row.names=FALSE)







#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# 10% 10% 10% 10% 10% 10% 10% 10% 10% 10% 10% 10% 10% 10% 10% 10% 10% 10% 10% 10%
## > For train 
train.n <- length(train.swed$TotalLoss)
err <- numeric(train.n)

# latent unknown factor...possible? [Risk]
# uuu = 1
# uuu = 10
# uuu = 20

set.seed(17)
for (i in 1:train.n) {
  if(train.swed$Experience[i]==0){
    
    uuu <- ifelse( train.swed$Risk[i]>=4, 1, 
                       ifelse(train.swed$Risk[i]<=2, 20, 10) )
    delta_sq <- ifelse( train.swed$SumInsAvg[i]<9800 & train.swed$SumInsAvg[i]>7.5, 0.1, 0.2 )

    err[i] <- rnorm( n=1, mean=0, sd=sqrt( uuu*delta_sq ) )
  }
  else{
    err[i] <- 0
  }
}
n.breaks = sqrt( nrow(train.swed) ) #******Rule of thumb
hist(err,  breaks=n.breaks)

train.swed$Error <- err

head(train.swed)
summary(train.swed)

train.swed$ln_insured <- log(train.swed$SumInsAvg)
train.swed$ln_insured_err <- train.swed$ln_insured + train.swed$Error 
hist(train.swed$ln_insured_err, breaks=n.breaks)

# check 
plot(train.swed$ln_insured, train.swed$Error, cex=0.5)
plot(train.swed$ln_insured, train.swed$ln_insured_err, cex=0.5)
plot(train.swed$Risk, train.swed$ln_insured, cex=0.5)



# 10% 10% 10% 10% 10% 10% 10% 10% 10% 10% 10% 10% 10% 10% 10% 10% 10% 10% 10% 10%
## > For test 

test.n <- length(test.swed$TotalLoss)
err <- numeric(test.n)

set.seed(13)
for (i in 1:test.n) {
  if(train.swed$Experience[i]==0){
    
    uuu <- ifelse( test.swed$Risk[i]>=4, 1, 
                       ifelse(test.swed$Risk[i]<=2, 20, 10) )
    delta_sq <- ifelse( test.swed$SumInsAvg[i]<9800 & test.swed$SumInsAvg[i]>7.5, 0.1, 0.2 )
    
    err[i] <- rnorm( n=1, mean=0, sd=sqrt( uuu*delta_sq ) )
  }
  else{
    err[i] <- 0
  }
}

n.breaks = sqrt( nrow(test.swed) ) #******Rule of thumb
hist(err,  breaks=n.breaks)

test.swed$Error <- err

head(test.swed)
summary(test.swed)

test.swed$ln_insured <- log(test.swed$SumInsAvg)
test.swed$ln_insured_err <- test.swed$ln_insured + test.swed$Error 
hist(test.swed$ln_insured_err, breaks=n.breaks)

# check 
plot(test.swed$ln_insured, test.swed$Error, cex=0.5)
plot(test.swed$ln_insured, test.swed$ln_insured_err, cex=0.5)


head(train.swe)
train.overall_E.swe = sum(abs(train.swed$Error)) /  sum(abs(train.swed$ln_insured )); train.overall_E.swe
# 6 %
test.overall_E.swe = sum(abs(test.swed$Error)) /  sum(abs(test.swed$ln_insured )); test.overall_E.swe
# 7 % --------------------------------------------------------------------------
# Therefore.....................................................................
# SAVE....----------------------------------------------------------------------
# write.csv(train.swed, "C:/Users/kimm4/Desktop/WORKSPACE/2nd 3rd paper/DATASET/BSwed.paper3-10.train.csv",
#           row.names=FALSE)
# write.csv(test.swed, "C:/Users/kimm4/Desktop/WORKSPACE/2nd 3rd paper/DATASET/BSwed.paper3-10.test.csv",
#           row.names=FALSE)



#-------------------------------------------------------------------------------
# 25% 25% 25% 25% 25% 25% 25% 25% 25% 25% 25% 25% 25% 25% 25% 25% 25% 25% 25% 25%
############################# Error For train ##################################
train.n <- length(train.swed$TotalLoss)
err <- numeric(train.n)

# latent unknown factor...possible? [Risk]
# uuu = 1
# uuu = 10
# uuu = 20

set.seed(17)
for (i in 1:train.n) {
  if(train.swed$Experience[i]==0){
    
    uuu <- ifelse( train.swed$Risk[i]>=4, 25, 
                   ifelse(train.swed$Risk[i]<=2, 50, 90) )
    delta_sq <- ifelse( train.swed$SumInsAvg[i]<9800 & train.swed$SumInsAvg[i]>7.5, 0.3, 0.5 )
    
    err[i] <- rnorm( n=1, mean=0, sd=sqrt( uuu*delta_sq ) )
  }
  else{
    err[i] <- 0
  }
}
n.breaks = sqrt( nrow(train.swed) ) #******Rule of thumb
hist(err,  breaks=n.breaks)

train.swed$Error <- err

head(train.swed)
summary(train.swed)

train.swed$ln_insured <- log(train.swed$SumInsAvg)
train.swed$ln_insured_err <- train.swed$ln_insured + train.swed$Error 
hist(train.swed$ln_insured_err, breaks=n.breaks)

# check 
plot(train.swed$ln_insured, train.swed$Error, cex=0.5)
plot(train.swed$ln_insured, train.swed$ln_insured_err, cex=0.5)
plot(train.swed$Risk, train.swed$ln_insured, cex=0.5)


# 25% 25% 25% 25% 25% 25% 25% 25% 25% 25% 25% 25% 25% 25% 25% 25% 25% 25% 25% 25%
################################## For test ####################################

test.n <- length(test.swed$TotalLoss)
err <- numeric(test.n)

set.seed(13)
for (i in 1:test.n) {
  if(train.swed$Experience[i]==0){
    
    uuu <- ifelse( test.swed$Risk[i]>=4, 25, 
                   ifelse(test.swed$Risk[i]<=2, 50, 90) )
    delta_sq <- ifelse( test.swed$SumInsAvg[i]<9800 & test.swed$SumInsAvg[i]>7.5, 0.3, 0.5 )
    
    err[i] <- rnorm( n=1, mean=0, sd=sqrt( uuu*delta_sq ) )
  }
  else{
    err[i] <- 0
  }
}

n.breaks = sqrt( nrow(test.swed) ) #******Rule of thumb
hist(err,  breaks=n.breaks)

test.swed$Error <- err

head(test.swed)
summary(test.swed)

test.swed$ln_insured <- log(test.swed$SumInsAvg)
test.swed$ln_insured_err <- test.swed$ln_insured + test.swed$Error 
hist(test.swed$ln_insured_err, breaks=n.breaks)

# check 
plot(test.swed$ln_insured, test.swed$Error, cex=0.5)
plot(test.swed$ln_insured, test.swed$ln_insured_err, cex=0.5)


head(train.swed)
train.overall_E.swe = sum(abs(train.swed$Error)) /  sum(abs(train.swed$ln_insured )); train.overall_E.swe
# 25 %
test.overall_E.swe = sum(abs(test.swed$Error)) /  sum(abs(test.swed$ln_insured )); test.overall_E.swe
# 25 % --------------------------------------------------------------------------
# Therefore.....................................................................
# SAVE....----------------------------------------------------------------------
# write.csv(train.swed, "C:/Users/kimm4/Desktop/WORKSPACE/2nd 3rd paper/DATASET/BSwed.paper3-25.train.csv",
#           row.names=FALSE)
# write.csv(test.swed, "C:/Users/kimm4/Desktop/WORKSPACE/2nd 3rd paper/DATASET/BSwed.paper3-25.test.csv",
#           row.names=FALSE)






















#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# 40% 40% 40% 40% 40% 40% 40% 40% 40% 40% 40% 40% 40% 40% 40% 40% 40% 40% 40% 40%
############################# Error For train ##################################
train.n <- length(train.swed$TotalLoss)
err <- numeric(train.n)

# latent unknown factor...possible? [Risk]
# uuu = 1
# uuu = 10
# uuu = 20

set.seed(17)
for (i in 1:train.n) {
  if(train.swed$Experience[i]==0){
    
    uuu <- ifelse( train.swed$Risk[i]>=4, 290, 
                   ifelse(train.swed$Risk[i]<=2, 390, 310) )
    delta_sq <- ifelse( train.swed$SumInsAvg[i]<9800 & train.swed$SumInsAvg[i]>7.5, 0.1, 0.2 )
    
    err[i] <- rnorm( n=1, mean=0, sd=sqrt( uuu*delta_sq ) )
  }
  else{
    err[i] <- 0
  }
}
n.breaks = sqrt( nrow(train.swed) ) #******Rule of thumb
hist(err,  breaks=n.breaks)

train.swed$Error <- err

head(train.swed)
summary(train.swed)

train.swed$ln_insured <- log(train.swed$SumInsAvg)
train.swed$ln_insured_err <- train.swed$ln_insured + train.swed$Error 
hist(train.swed$ln_insured_err, breaks=n.breaks)

# check 
plot(train.swed$ln_insured, train.swed$Error, cex=0.5)
plot(train.swed$ln_insured, train.swed$ln_insured_err, cex=0.5)
plot(train.swed$Risk, train.swed$ln_insured, cex=0.5)



# 40% 40% 40% 40% 40% 40% 40% 40% 40% 40% 40% 40% 40% 40% 40% 40% 40% 40% 40% 40%
## > For test 

test.n <- length(test.swed$TotalLoss)
err <- numeric(test.n)

set.seed(13)
for (i in 1:test.n) {
  if(train.swed$Experience[i]==0){
    
    uuu <- ifelse( test.swed$Risk[i]>=4, 290, 
                   ifelse(test.swed$Risk[i]<=2, 390, 310) )
    delta_sq <- ifelse( test.swed$SumInsAvg[i]<9800 & test.swed$SumInsAvg[i]>7.5, 0.1, 0.2 )
    
    err[i] <- rnorm( n=1, mean=0, sd=sqrt( uuu*delta_sq ) )
  }
  else{
    err[i] <- 0
  }
}

n.breaks = sqrt( nrow(test.swed) ) #******Rule of thumb
hist(err,  breaks=n.breaks)

test.swed$Error <- err

head(test.swed)
summary(test.swed)

test.swed$ln_insured <- log(test.swed$SumInsAvg)
test.swed$ln_insured_err <- test.swed$ln_insured + test.swed$Error 
hist(test.swed$ln_insured_err, breaks=n.breaks)

# check 
plot(test.swed$ln_insured, test.swed$Error, cex=0.5)
plot(test.swed$ln_insured, test.swed$ln_insured_err, cex=0.5)


head(train.swed)
train.overall_E.swe = sum(abs(train.swed$Error)) /  sum(abs(train.swed$ln_insured )); train.overall_E.swe
# 40 %
test.overall_E.swe = sum(abs(test.swed$Error)) /  sum(abs(test.swed$ln_insured )); test.overall_E.swe
# 44 % --------------------------------------------------------------------------
# Therefore.....................................................................
# SAVE....----------------------------------------------------------------------
# write.csv(train.swed, "C:/Users/kimm4/Desktop/WORKSPACE/2nd 3rd paper/DATASET/BSwed.paper3-40.train.csv",
#           row.names=FALSE)
# write.csv(test.swed, "C:/Users/kimm4/Desktop/WORKSPACE/2nd 3rd paper/DATASET/BSwed.paper3-40.test.csv",
#           row.names=FALSE)

























#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

###################
# with French Data: ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
###################

train.fren <- read.csv("C:/Users/kimm4/Desktop/WORKSPACE/2nd 3rd paper/DATASET/BFrench.train.csv",
                       header=T, 
                       na.strings=c("."), 
                       stringsAsFactors=F)
test.fren <- read.csv("C:/Users/kimm4/Desktop/WORKSPACE/2nd 3rd paper/DATASET/BFrench.test.csv",
                      header=T, 
                      na.strings=c("."), 
                      stringsAsFactors=F)

head(train.fren)
summary(train.fren)
table(train.fren$CarAge)
#    0     1* 
#26435 14078 
cor(train.fren$ClaimAmount, train.fren$Exposure) # cor: 0.02


quantile(train.fren$ClaimAmount, c(0.025, 0.975))
#      2.5%      97.5% 
#    68.000   7890.914  
quantile(train.fren$Density, c(0.025, 0.975))
#      2.5%      97.5% 
#       14      17140  
quantile(train.fren$Exposure, c(0.025, 0.975))
#      2.5%      97.5% 
#      0.07         1  


#### Add Error!
#-------------------------------------------------------------------------------
# 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1%
## > For train 
train.n <- length(train.fren$ClaimAmount)
err <- numeric(train.n)

set.seed(17)
for (i in 1:train.n) {
  if(train.fren$CarAge[i]==1){ #very old car..
    
    uuu <- ifelse( train.fren$Density[i]>=17000, 0.01, 
                   ifelse(train.fren$Density[i]<=20, 0.0012, 0.005) )
    delta_sq <- ifelse( train.fren$Exposure[i]<1 & train.fren$Exposure[i]>0.07, 0.1, 0.2 )
    
    err[i] <- rnorm( n=1, mean=0, sd=sqrt( uuu*delta_sq ) )
  }
  else{
    err[i] <- 0
  }
}
n.breaks = sqrt( nrow(train.fren) ) #******Rule of thumb
hist(err,  breaks=n.breaks)

train.fren$Error <- err

head(train.fren)
summary(train.fren)

train.fren$ln_expo <- log(train.fren$Exposure)
train.fren$ln_expo_err <- train.fren$ln_expo + train.fren$Error 


# 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1% 1%
## > For test 

test.n <- length(test.fren$ClaimAmount)
err <- numeric(test.n)

for (i in 1:test.n) {
  if(test.fren$CarAge[i]==1){ #very old car..
    
    uuu <- ifelse( test.fren$Density[i]>=17000, 0.01, 
                   ifelse(test.fren$Density[i]<=20, 0.0012, 0.05) )
    delta_sq <- ifelse( test.fren$Exposure[i]<1 & test.fren$Exposure[i]>0.07, 0.1, 0.2 )
    
    err[i] <- rnorm( n=1, mean=0, sd=sqrt( uuu*delta_sq ) )
  }
  else{
    err[i] <- 0
  }
}
n.breaks = sqrt( nrow(test.fren) ) #******Rule of thumb
hist(err,  breaks=n.breaks)

test.fren$Error <- err

head(test.fren)
summary(test.fren)

test.fren$ln_expo <- log(test.fren$Exposure)
test.fren$ln_expo_err <- test.fren$ln_expo + test.fren$Error

head(train.fren)
train.overall_E.fren = sum(abs(train.fren$Error)) /  sum(abs(train.fren$ln_expo )); train.overall_E.fren
# 1 %
test.overall_E.fren = sum(abs(test.fren$Error)) /  sum(abs(test.fren$ln_expo )); test.overall_E.fren
# 4 % -------------------------------------------------------------------------
#Therefore.....................................................................
#SAVE....----------------------------------------------------------------------
# write.csv(train.fren, "C:/Users/kimm4/Desktop/WORKSPACE/2nd 3rd paper/DATASET/BFren.paper3-1.train.csv",
#           row.names=FALSE)
# write.csv(test.fren, "C:/Users/kimm4/Desktop/WORKSPACE/2nd 3rd paper/DATASET/BFren.paper3-1.test.csv",
#           row.names=FALSE)




#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# 10% 10% 10% 10% 10% 10% 10% 10% 10% 10% 10% 10% 10% 10% 10% 10% 10% 10% 10% 10%
## > For train 
train.n <- length(train.fren$ClaimAmount)
err <- numeric(train.n)

set.seed(17)
for (i in 1:train.n) {
  if(train.fren$CarAge[i]==1){ #very old car..
    
    uuu <- ifelse( train.fren$Density[i]>=17000, 1, 
                   ifelse(train.fren$Density[i]<=20, 0.12, 0.5) )
    delta_sq <- ifelse( train.fren$Exposure[i]<1 & train.fren$Exposure[i]>0.07, 0.1, 0.2 )
    
    err[i] <- rnorm( n=1, mean=0, sd=sqrt( uuu*delta_sq ) )
  }
  else{
    err[i] <- 0
  }
}
n.breaks = sqrt( nrow(train.fren) ) #******Rule of thumb
hist(err,  breaks=n.breaks)

train.fren$Error <- err

head(train.fren)
summary(train.fren)

train.fren$ln_expo <- log(train.fren$Exposure)
train.fren$ln_expo_err <- train.fren$ln_expo + train.fren$Error 


# 10% 10% 10% 10% 10% 10% 10% 10% 10% 10% 10% 10% 10% 10% 10% 10% 10% 10% 10% 10%
## > For test 

test.n <- length(test.fren$ClaimAmount)
err <- numeric(test.n)

for (i in 1:test.n) {
  if(test.fren$CarAge[i]==1){ #very old car..
    
    uuu <- ifelse( test.fren$Density[i]>=17000, 1, 
                   ifelse(test.fren$Density[i]<=20, 0.12, 0.5) )
    delta_sq <- ifelse( test.fren$Exposure[i]<1 & test.fren$Exposure[i]>0.07, 0.1, 0.2 )
    
    err[i] <- rnorm( n=1, mean=0, sd=sqrt( uuu*delta_sq ) )
  }
  else{
    err[i] <- 0
  }
}
n.breaks = sqrt( nrow(test.fren) ) #******Rule of thumb
hist(err,  breaks=n.breaks)

test.fren$Error <- err

head(test.fren)
summary(test.fren)

test.fren$ln_expo <- log(test.fren$Exposure)
test.fren$ln_expo_err <- test.fren$ln_expo + test.fren$Error

head(train.fren)
train.overall_E.fren = sum(abs(train.fren$Error)) /  sum(abs(train.fren$ln_expo )); train.overall_E.fren
# 12.7 %
test.overall_E.fren = sum(abs(test.fren$Error)) /  sum(abs(test.fren$ln_expo )); test.overall_E.fren
# 13.2 % -------------------------------------------------------------------------
#Therefore.....................................................................
#SAVE....----------------------------------------------------------------------
# write.csv(train.fren, "C:/Users/kimm4/Desktop/WORKSPACE/2nd 3rd paper/DATASET/BFren.paper3-10.train.csv",
#           row.names=FALSE)
# write.csv(test.fren, "C:/Users/kimm4/Desktop/WORKSPACE/2nd 3rd paper/DATASET/BFren.paper3-10.test.csv",
#           row.names=FALSE)






#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# 40% 40% 40% 40% 40% 40% 40% 40% 40% 40% 40% 40% 40% 40% 40% 40% 40% 40% 40% 40%
############################# Error For train ##################################
train.n <- length(train.fren$ClaimAmount)
err <- numeric(train.n)

# latent unknown factor...possible? [Density]
# uuu = 2
# uuu = 10
# uuu = 20

set.seed(17)
for (i in 1:train.n) {
  if(train.fren$CarAge[i]==1){ #very old car..
    
    uuu <- ifelse( train.fren$Density[i]>=17000, 8, 
                       ifelse(train.fren$Density[i]<=20, 2, 5) )
    delta_sq <- ifelse( train.fren$Exposure[i]<1 & train.fren$Exposure[i]>0.07, 0.1, 0.2 )
    
    err[i] <- rnorm( n=1, mean=0, sd=sqrt( uuu*delta_sq ) )
  }
  else{
    err[i] <- 0
  }
}
n.breaks = sqrt( nrow(train.fren) ) #******Rule of thumb
hist(err,  breaks=n.breaks)

train.fren$Error <- err

head(train.fren)
summary(train.fren)

train.fren$ln_expo <- log(train.fren$Exposure)
train.fren$ln_expo_err <- train.fren$ln_expo + train.fren$Error 
hist(train.fren$ln_expo_err, breaks=n.breaks)

# check 
plot(train.fren$ln_expo, train.fren$Error, cex=0.5)
plot(train.fren$ln_expo, train.fren$ln_expo_err, cex=0.5)
plot(train.fren$Density, train.fren$ln_expo, cex=0.5)


# 40% 40% 40% 40% 40% 40% 40% 40% 40% 40% 40% 40% 40% 40% 40% 40% 40% 40% 40% 40%
## > For test 
test.n <- length(test.fren$ClaimAmount)
err <- numeric(test.n)

set.seed(17)
for (i in 1:test.n) {
  if(test.fren$CarAge[i]==1){ #very old car..
    
    uuu <- ifelse( test.fren$Density[i]>=17000, 8, 
                       ifelse(test.fren$Density[i]<=20, 2, 5) )
    delta_sq <- ifelse( test.fren$Exposure[i]<1 & test.fren$Exposure[i]>0.07, 0.1, 0.2 )
    
    err[i] <- rnorm( n=1, mean=0, sd=sqrt( uuu*delta_sq ) )
  }
  else{
    err[i] <- 0
  }
}
n.breaks = sqrt( nrow(test.fren) ) #******Rule of thumb
hist(err,  breaks=n.breaks)

test.fren$Error <- err

head(test.fren)
summary(test.fren)

test.fren$ln_expo <- log(test.fren$Exposure)
test.fren$ln_expo_err <- test.fren$ln_expo + test.fren$Error 
hist(test.fren$ln_expo_err, breaks=n.breaks)

# check 
plot(test.fren$ln_expo, test.fren$Error, cex=0.5)
plot(test.fren$ln_expo, test.fren$ln_expo_err, cex=0.5)
plot(test.fren$Density, test.fren$ln_expo, cex=0.5)


head(train.fren)
train.overall_E.fren = sum(abs(train.fren$Error)) /  sum(abs(train.fren$ln_expo )); train.overall_E.fren
# 40 %
test.overall_E.fren = sum(abs(test.fren$Error)) /  sum(abs(test.fren$ln_expo )); test.overall_E.fren
# 42 % -------------------------------------------------------------------------
#Therefore.....................................................................
#SAVE....----------------------------------------------------------------------
# write.csv(train.fren, "C:/Users/kimm4/Desktop/WORKSPACE/2nd 3rd paper/DATASET/BFren.paper3-40.train.csv",
#           row.names=FALSE)
# write.csv(test.fren, "C:/Users/kimm4/Desktop/WORKSPACE/2nd 3rd paper/DATASET/BFren.paper3-40.test.csv",
#           row.names=FALSE)



################################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++ Let's Begin ! +++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
################################## > DP.intro < ################################


train.bra <- read.csv("C:/Users/kimm4/Desktop/WORKSPACE/2nd 3rd paper/DATASET/Brazil.paper3-1.train.csv",
                      header=T, 
                      na.strings=c("."), 
                      stringsAsFactors=F)
test.bra <- read.csv("C:/Users/kimm4/Desktop/WORKSPACE/2nd 3rd paper/DATASET/Brazil.paper3-1.test.csv",
                     header=T, 
                     na.strings=c("."), 
                     stringsAsFactors=F)

train.swe <- read.csv("C:/Users/kimm4/Desktop/WORKSPACE/2nd 3rd paper/DATASET/BSwed.paper3-1.train.csv",
                      header=T, 
                      na.strings=c("."), 
                      stringsAsFactors=F)
test.swe <- read.csv("C:/Users/kimm4/Desktop/WORKSPACE/2nd 3rd paper/DATASET/BSwed.paper3-1.test.csv",
                     header=T, 
                     na.strings=c("."), 
                     stringsAsFactors=F)

train.fre <- read.csv("C:/Users/kimm4/Desktop/WORKSPACE/2nd 3rd paper/DATASET/BFren.paper3-1.train.csv",
                      header=T, 
                      na.strings=c("."), 
                      stringsAsFactors=F)
test.fre <- read.csv("C:/Users/kimm4/Desktop/WORKSPACE/2nd 3rd paper/DATASET/BFren.paper3-1.test.csv",
                     header=T, 
                     na.strings=c("."), 
                     stringsAsFactors=F)

# [note] Dont be confused. ----------------------------------------------------#
#sumInsAvg: the payment for loss will be in proportion to the value insured.
# Claim amount = (Actual loss × Insured amount) / Value of property
# pure_premium(loss cost) * exposure = total claim 
#
# The sum insured (planned loss?) is the amount of money that an insurer is obligated 
# to cover in the event of a covered damage. 
# The sum insured correlates directly to the amount of premium you pay, but not 
# always to the property's actual value or asset insured. If the sum insured is 
# less than the required amount to rebuild, it is often referred to as underinsurance.
#------------------------------------------------------------------------------#


