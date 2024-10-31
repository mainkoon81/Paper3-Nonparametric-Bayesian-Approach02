library(mvnfast)
library(tidyverse)
library(Rcpp)
library(RcppArmadillo)
#library(extraDistr)
library(MCMCpack)


#> C++ for Gold standard
sourceCpp("C:/Users/kimm4/Desktop/WORKSPACE/2nd 3rd paper/DATASET/Basic_training/project3.DPcpp.Basic.cpp")
#sourceCpp("project3.Benchmark.DPcpp.cpp")

#> C++ for Gustafson Correction
sourceCpp("C:/Users/kimm4/Desktop/WORKSPACE/2nd 3rd paper/DATASET/Basic_training/project3.DPcpp.NDB.cpp")
#sourceCpp("project3.DPcpp.NDB.cpp")





#-------------------------------------------------------------------------------
##### Swedish Data with overall error: [1%]
#-------------------------------------------------------------------------------
full.train.sweden <- read.csv("C:/Users/kimm4/Desktop/WORKSPACE/2nd 3rd paper/DATASET/BSwed.paper3-1.train.csv")
full.test.sweden <- read.csv("C:/Users/kimm4/Desktop/WORKSPACE/2nd 3rd paper/DATASET/BSwed.paper3-1.test.csv")
#full.train.sweden <- read.csv("BSwed.paper3-1.train.csv")
#full.test.sweden <- read.csv("BSwed.paper3-1.test.csv")
head(full.train.sweden)
str(full.train.sweden)
summary(full.train.sweden)


n.train <- nrow(full.train.sweden)
n.test <- nrow(full.test.sweden)

cor(x=full.train.sweden, method="pearson")
cor(x=full.train.sweden, method="spearman")


#%%%% train %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
##### Variable Definitions
S = full.train.sweden$TotalLoss
logS = log(S)
Z = full.train.sweden$Experience

# dirty data on log scale, divided by standard dev to be on same scale
# Xstar = full.train.sweden$ln_insured_err/sd(full.train.sweden$ln_insured_err) 
Xstar = full.train.sweden$ln_insured_err
X_tru = full.train.sweden$ln_insured
# or...X on log scale, divided by standard dev to be on same scale
#X = log(full.train.sweden$SumInsAvg) / sd( log(full.train.sweden$SumInsAvg) )  

matXstar = as.matrix(cbind(1, Xstar, Z))
matX_tru = as.matrix(cbind(1, X_tru, Z))
#matX = as.matrix(cbind(1, X, Z))

n.breaks.train = sqrt( nrow(full.train.sweden) ) #******Rule of thumb


#%%%% test %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
##### Variable Definitions
S.test = full.test.sweden$TotalLoss
logS.test = log(S.test)
Z.test = full.test.sweden$Experience

# dirty data on log scale, divided by standard dev to be on same scale
# Xstar = full.train.sweden$ln_insured_err/sd(full.train.sweden$ln_insured_err) 
Xstar.test = full.test.sweden$ln_insured_err
X_tru.test = full.test.sweden$ln_insured
#X.test = log(full.test.sweden$SumInsAvg) / sd( log(full.test.sweden$SumInsAvg) ) 

matXstar.test = as.matrix(cbind(1, Xstar.test, Z.test))
matX_tru.test = as.matrix(cbind(1, X_tru.test, Z.test))
#matX.test = as.matrix(cbind(1, X.test, Z.test))

n.breaks.test = sqrt( nrow(full.test.sweden) ) #******Rule of thumb


####] Now, define the outcome model!
# - Log-skewnormal density ? no no no ... our outcome is based on "LogS"

# - Skew Normal Density
dsknorm = function(ls, x, beta, sig2, xi) {
  mu = sum(x*beta) # x^T*B
  sig = sqrt(sig2)
  z = (ls-mu)/sig
  return(2/sig*dnorm(z)*pnorm(xi*z))
}

# - Log of SkewNormal...better.................???
dsknorm_log = function(ls, x, beta, sig2, xi) {
  mu = sum(x*beta) # x^T*B
  sig = sqrt(sig2)
  z = (ls-mu)/sig
  return(log(2) - log(sig) + dnorm(z, log=TRUE) + pnorm(xi*z, log=TRUE))
}

# - Logskew normal density outcome model
# .... LogSkewNormal.....too tiny....................???
dlogsknorm = function(s, x, beta, sig2, xi) {
  mu = sum(x*beta) 
  sig = sqrt(sig2)
  z = (log(s)-mu)/sig
  return(2/(s*sig)*dnorm(z)*pnorm(xi*z))
}
# .... Log of LogSkewNormal...better.................???
dlogsknorm_log = function(s, x, beta, sig2, xi) {
  mu = sum(x*beta)
  sig = sqrt(sig2)
  z = (log(s)-mu)/sig
  return(log(2)-log(s*sig) + dnorm(z, log = TRUE) + pnorm(xi*z, log = TRUE))
}

# # The delta: p(Y > 1|X): Sigmoid
# sigmoid = function(x) {
#   1/(1+exp(-x))
# } 
# curve(sigmoid(x), from=-5, to=5)
#-------------------------------------------------------------------------------




library(GGally)
library(ggplot2)
df = data.frame( S, X_tru, factor(Z) )

p <- ggpairs(df, diag = list(continuous = wrap("barDiag", fill = "white", color = "black")),
             lower = list(continuous = wrap("smooth", alpha = 0.3), 
                          combo = wrap("facethist", fill = "white", color = "black")),
             upper = list(continuous = wrap("cor", size = 10)))   # Adjust size as needed 
p + theme_minimal() + theme(panel.background = element_rect(fill = "white", color = "black"))




################################################################################
############################# S ~ DP LSN (clean) ###############################
################################################################################

################################################################################
### Step00> Define prior W/O consideration of clusters
################################################################################

# ----------- Outcome ----- # ::: for S ~ LSN( X*betaparam, sig2param, xiparam )

# PRIOR: "betaparam" ~ MVN(m0, SIG_b0)

#> Gaussian regression for initialize OUTCOME model parameter m0, SIG_b0
fit1 <- glm( S ~  X_tru + factor(Z), family=Gamma(link = "log") )
#fit.tes <- glm( S ~  X_tru + factor(Z), family=poisson(link = "log") )
# fit1 <- lm( logS ~ Xstar + factor(Z) )
options(scipen = 999)
summary(fit1)
#summary(fit.tes)
#vcov(fit.tes)


m0 = coef(fit1)    # 1x3 initial reg_coeff vector (sort of "means"): Regression result as "mean"
SIG_b0 = vcov(fit1)   # 3x3 initial cov matrix of reg_coeff
p = length(m0)
options(scipen = 999)
SIG_b0
SIG_b0inv = solve(a=SIG_b0) # inverse of cov matrix of reg_coeff for later use (for posterior on reg_coeff)!!

#****** for VAR Inflation factor
delta = 0.01             # equal importance on m0 and beta ????
varinf = n.train/100 # rule of thumb divide by 5

# #> SIG_b0
# #- it requires -------
# qu0 = p+2
# # LAMB ?????????????????????????????????????????????
# LAMB = SIG_b0


###> PRIOR: "sig2param" ~ IG( u0, v0 )
# Hyperprior               u0 ~ Fink( rho_u1, rho_u2 )
# Hyperprior                   v0 ~ Ga( rho_bv1, rho_v2 )

u0=1.7591
v0=7.113

# # Hyperparameter u0
# rho_u1 = 1/8
# rho_u2 = 3/2 
# 
# # Hyperparameter v0
# rho_v1 = 8 
# rho_v2 = 1 
# 
# # Based on the result? mean(unlist(list_u0param)), mean(unlist(list_v0param))
# # E[invGa] = scale/(shape-1) = 50
# u0_new = 1.7591 #3 #2 #1.1 # from 3/7  # from # summary(unlist(list_u0param)) # for u0_new                            
# v0_new = 7.113  #3 #5 #8.5 # from 3/7  # from # summary(unlist(list_v0param)) # for v0_new

# ## [CHECK later] ---------------------------------------------------------------
# # How can we design proposal for sig2 ~ InvGa(.) ?
# par(mfrow=c(1,1))
# plot( density(unlist(list_sig2param)) )                                   # sig2 we have...before
# curve( dinvgamma(x, shape=mean(unlist(list_u0param)), 
#                  scale=mean(unlist(list_v0param))), add=TRUE, col="red") # sig2 we obtain....after
# curve( dinvgamma(x, shape=u0_new, 
#                  scale=v0_new), add=TRUE, col="blue")                    # sig2 proposal adjustment
# 
# table(unlist(list_sig2param))
# sort(table(unlist(list_sig2param)),DESC=TRUE)                             # to see if they are appropriate...
# #-------------------------------------------------------------------------------

###> PRIOR: "xiparam" ~ t(loc, nu0, sca)
#loc = 0   # ? Location hyperparameter for T distribution
#nu0 = 1/2 # df for T distribution on pigtail parameter
#sca = 5   # ? Scale hyperparameter for T distribution 
nu0 = n.train-1
#- df refers to the number of independent observations (N-1):

# ------------------------------------------------------------------------------ #



# ::::::: # ------------------------ PRIOR ------------------------- # ::::::: #
# ------- Covariate ------ ::: for Z~Bin( 1, [piparam] ) where 
#                                             "piparam" ~ Beta(g0, h0) 
g0 = 0.5 
h0 = 0.5   

#] ------ Covariate ------ ::: for X ~ N( X_bar, [lambda2param] ) where 
# PRIOR                                            "lambda2param" ~ IG(c0, d0)
#------------------------------------------------------------------------------------------------------------------------
#X_bar = mean(X_tru) #????

#> lambda2
c0 = 0.5          
d0 = 0.5

# # ---- Covariate model ----# ::: for X|Z ~ N( K0+K1*Z, [lambda2param] ) where 
# #                                                        "lambda2param" ~ IG(c0, d0) 
# fit2 <- lm(X ~ Z)
# summary(fit2)
# 
# kappa_v0 <- coef(fit2)
# SIG_k0 <- vcov(fit2)*varinf
# SIG_inv_k0 <- solve(SIG_k0)
# 
# c0 = 0.5          
# d0 = 0.5    
# 
# KK1 = cbind(1, Z) # for the analytical solution for Parameter-Free covariate model development later
# # [note] 
# # - we need KK1 matrix for all obvs w/o concerning cluster
# # - KK1_j : cluster-wise 
# # - KK2_j         : cluster-wise


# ::::::: # ------------------------ PRIOR ------------------------- # ::::::: #
# ----------- precision -----------# :::  for alphaparam ~ Ga( gamma0, psi0 )
gamma0 = 1    #                                          
psi0 = 1    



################################################################################
### Step01> initialize cluster membership - Kmeans clustering
################################################################################
# J=3
# set.seed(1)
# cluster = kmeans(cbind(S,X_tru,Z), J)
# cl_membership = cluster$cluster 
# table(cl_membership)

cl_membership <- full.train.sweden$Risk
table(cl_membership)
cl_membership.test <- full.test.sweden$Risk
table(cl_membership.test)
J=5





################################################################################
### Step02> initialize major parametersssss by clusters using your "POSTERIOR"
#                                                 Based on kmean clustering: J=3
################################################################################
########################--- for Covariates + alpha ---##########################
################################################################################

##[A] for Z ~ Bin( 1, [piparamS_j] )
#                     "piparamS_j" ~ Beta(g0, h0)
#-------------------------------------------------------------------------------
# Covariate parameter from Posterior
piparam = numeric(J)

#> 1. Regular
for (j in 1:J) {
  Zj=Z[cl_membership==j]
  nj=length(Zj)
  piparam[j] = rbeta( n=1, shape1=g0+sum(Zj), shape2=h0+nj-sum(Zj) ) 
}
# #> 2. MAR? For Z ~ Bin( n, "prob_j" )
# Z.miss = is.na(Z)
# for (j in 1:J) {
#   Zj=Z[cl_membership==j & !Z.miss]
#   nj=length(Zj)
#   if(nj > 0) { #sample from posterior
#     piparam[j] = rbeta( n=1, shape1=g0+sum(Zj), shape2=h0+nj-sum(Zj) )
#   }
#   else { #sample from prior
#     piparam[j] = rbeta(n=1, shape1 = g0, shape2 = h0)
#   }
# }
piparam #% pi parameter sampled from posterior

#% [Check!] empirical pi ? using "aggregate(.)": investigation by splitting the data into subset
pi_empirical = aggregate(x=Z, by=list(cl_membership), FUN=mean, na.rm= T)$x; pi_empirical



##[B] for X ~ N( X_bar, [lambda2param] ) where 
# PRIOR                 "lambda2param" ~ IG(c0, d0)
#-------------------------------------------------------------------------------
# Covariate parameter from Posterior
lambda2param = numeric(J)

#> 1. Regular
for (j in 1:J) {
  Xj=X_tru[cl_membership==j]
  n_j=length(Xj)
  lambda2param[j]=rinvgamma( n=1, shape=(c0+n_j)/2, scale=(d0 + sum( (Xj - mean(Xj))^2 ) )/2 )  
}                                                                                                     
lambda2param #% lambda2 parameter (posterior sample)
#% empirical X_barF ? using "aggregate(.)": investigation by splitting the data into subset
X_bar = aggregate(x=X_tru, by=list(cl_membership), FUN=mean, na.rm= T)$x; X_bar
# ---------------------------------------------------------------- Perfect!!!!!!


# #> 2. NDB? for X|Z ~ N( "kappa0_j"+"kappa1_j"*Z, [lambda2param] ) 
# # Sample it from Normal(for kappa), ............ "lambda2param" ~ IG(c0, d0)
# kappaparam = matrix(nrow = J, ncol = length(kappa_v0)); kappaparam # vectors for each cluster
# lambda2param = numeric(J); lambda2param
# # [note] 
# # - kappa_v0: initial intercept, initial slope
# # - kappa0_j from kappaparam : intercepts for each cluster
# # - kappa1_j from kappaparam : slope for each cluster
# for (j in 1:J) {
#   Xj = X[cl_membership==j]
#   Zj = Z[cl_membership==j]
#   nj = length(Xj)
#   lambda2param[j]=rinvgamma( n=1, shape=(nj+c0)/2, scale=0.5*( sum((Xj - kappa_v0[1]-kappa_v0[2]*Zj)^2) + d0 ) )
#   
#   KK1j = cbind(1, Zj)
#   KK2j = matrix(c(sum(Xj), sum(Xj*Zj)), nrow =2, ncol = 1) # 2x1 vector
#   
#   SIG_inv_j = solve(SIG_inv_k0+(t(KK1j) %*% KK1j)) # solve( ..[matrix].. ) = [matrix]^(-1)
#   kappaparam[j,] <- rmvn(n=1, mu = SIG_inv_j %*% (SIG_inv_k0 %*% kappa_v0 + KK2j), sigma = lambda2param[j]*SIG_inv_j)      
# }
# kappaparam   #% kappa parameter sampled from posterior
# lambda2param #% lambda2 parameter sampled from posterior
# 
# #% [Check!] empirical kappa, lambda ? using "aggregate(.)": investigation by splitting the data into subset
# kappa_empirical = rbind(coef(lm(X[cl_membership==1] ~ Z[cl_membership==1])),
#                       coef(lm(X[cl_membership==2] ~ Z[cl_membership==2])),
#                       coef(lm(X[cl_membership==3] ~ Z[cl_membership==3]))); kappa_empirical
# lambda_empirical = aggregate(x=as.numeric(X), by=list(cl_membership), FUN=var, na.rm=T)$x; lambda_empirical



##[C] for alphaparam ~ Ga( shape=gamma0, rate=psi0 ) .. starting with J clusters for Mixed Gamma POSTERIOR
alpha0=2 # initialization: typically..1,2...
eta = rbeta(n=1, shape1=alpha0+1, shape2=n.train); eta
pi_eta = (gamma0+J-1)/(gamma0+J-1 + n.train*(psi0-log(eta))); pi_eta

# precision parameters and its done!
alphaparam = pi_eta*rgamma( n=1, shape=gamma0+J, 
                       rate=psi0-log(eta) ) + (1-pi_eta)*rgamma( n=1, shape=gamma0+J-1, rate=psi0-log(eta) ); alphaparam




################################################################################
#############################--- for outcome ---################################
################################################################################

##[D] ------ Outcome S ------ ::: [betaparam], [sig2param], [xiparam].... try MH-algorithm
# - No conjugacy...so prepare MM sampling

# This is where new parameters go.
betaparam = matrix(data=NA, nrow=J, ncol=length(m0)) 
sig2param = numeric(J) 
xiparam = numeric(J)

# Set in place the old parameters..sample from prior...now..it's temporary values
beta_prev = matrix( rep(rmvn(n=1, mu=m0, sigma=SIG_b0*varinf), J), nrow=J, byrow= TRUE )
sig2_prev = rep(rinvgamma(n=1, shape=u0, scale=v0), J)
xi_prev = rep(rt(n=1, df=nu0), J) # imagine we already have them...old days..


#### Metropolis Hastings to update betaparam, sig2param, xiparam
for(j in 1:J) {
  # Sample proposals from priors
  beta_prop = rmvn(n=1, mu=m0, sigma=SIG_b0*varinf)
  sig2_prop = rinvgamma(n=1, shape=u0, scale=v0)
  xi_prop = rt(n=1, df=nu0)
  
  # subsetting by cluster
  Sj = S[cl_membership==j]
  matXj = matX_tru[cl_membership==j, ]
  
  # #>>>>>>>>>>>>>>>>>>>> if you have MAR.....impute first! <<<<<<<<<<<<<<<<<<<#
  # Z.miss_j = Z.miss[cl_membership==j]
  # miss.indx = which(Z.miss_j)
  # # p0: joint where x1 = 0
  # # p1: joint where x1 = 1
  # for(i in miss.indx) {
  #   if (Sj[i]>1) {
  #     p0log = 
  #       dlogsknorm_log(y=Sj[i], 
  #                      x=t(as.matrix(c(matXj[i,1], 0, matXj[i,3]), nrow=1)), 
  #                      beta=beta_prev[j,], 
  #                      sig2=sig2_prev[j], 
  #                      xi=xi_prev[j]) + 
  #       dbinom(x=0, size=1, piparam[j], log=TRUE) + 0
  #       #log( 1-sigmoid( sum(c(matXj[i,1], 0, matXj[i,3])*betat_old[j, ]) ) )
  #     
  #     p1log = 
  #       dlogsknorm_log(y=Sj[i], 
  #                      x=t(as.matrix(c(matXj[i,1], 1, matXj[i,3]), nrow=1)), 
  #                      beta=beta_prev[j,], 
  #                      sig2=sig2_prev[j], 
  #                      xi=xi_prev[j]) + 
  #       dbinom(x=1, size=1, piparam[j], log=TRUE) + 0
  #       #log( 1-sigmoid( sum(c(matXj[i,1], 1, matXj[i,3])*betat_old[j, ]) ) )
  #   }
  #   else {
  #     p0log = 
  #       dbinom(x=0, size=1, piparam[j], log=TRUE) + 0 # (1 - pi_j)
  #       #log( sigmoid( sum(c(matXj[i,1], 0, matXj[i,3])*betat_old[j, ]) ) )
  #     
  #     p1log = 
  #       dbinom(x=1, size=1, piparam[j], log=TRUE) + 0 # (pi_j)
  #       #log( sigmoid( sum(c(matXj[i,1], 1, matXj[i,3])*betat_old[j, ]) ) )
  #   }
  #   
  #   # Let's impute!!!
  #   matXj[i,2] = rbinom(n=1, size=1, prob=1/(1+exp(p0log-p1log))) #imputing NA with using posterior Pi 
  # }
  
  ##> A. MH algorithm for: [betaparam]
  #-Note: "beta_prop" is from PRIOR proposal (new) while "beta_prev[j,]" is from the single POSTERIOR sample (old). 
  numerator=0
  denominator=0
  for (i in 1:length(Sj)){
    numerator = numerator + 
      dlogsknorm_log(Sj[i], matXj[i,], t(as.matrix(beta_prop)), sig2_prop, xi_prop)
    denominator = denominator + 
      dlogsknorm_log(Sj[i], matXj[i,], t(as.matrix(beta_prev[j,])), sig2_prev[j], xi_prev[j])
  }
  # compute the ratio
  ratio = min(exp(numerator-denominator), 1)
  
  U = runif(n = 1, min = 0, max = 1)
  if(U < ratio) {
    betaparam[j, ] = beta_prop
    sig2param[j] = sig2_prop
    xiparam[j] = xi_prop
  } 
  else {
    betaparam[j, ] = beta_prev[j, ]
    sig2param[j] = sig2_prev[j]
    xiparam[j] = xi_prev[j]
  }
}
betaparam   #(J,p) : cl.membership(J), predictors(p)
beta_prev #(J,p)
sig2param   #(J)
sig2_prev #(J)
xiparam     #(J)
xi_prev   #(J)



################################################################################
### Step03> Data model development #############################################
########## Using sampled parameters above, ---- discrete / continuous ##########
################################################################################

### [discrete Data model]
#
# covariate w/o NA: ----------> Z,X
#-> dbinom(x, size, prob=piparam)*dnorm(x, mean=X_bar, sd=lambda2param)
# covariate with NA in Z: ---> X only
#-> dnorm(x, mean=X_bar, sd=lambda2param)

# Outcome w/o NA: ------------> S|Z,X
#-> dlogsknorm_log(y, x, beta=betaparam, sig2=sig2param, xi=xiparam)
# Outcome with NA: -----------> Sh|Z,X then "MARGINALIZE" w.r.t "Z"
#-> dlogsknorm_log(y, x, beta=betaparam, sig2=sig2param, xi=xiparam) ???

################################################################################
######### [continuous (Parameter-free) Data model for training Data] ###########
################################################################################
# 1. param-free outcome density : f0S
# 2. param-free covariate density : f0X=f0x1(Z[i])*f0x2(X[i])
# 3. param-free E[outcome] value : E0S


##### For [Param-free covariate model]: 
# covariate model w/o NA: ----------> int Z,X... w.r.t "piparam","X_bar","lambda2param"
# ::: piparam[j] = rbeta( n=1, shape1=g0+sum(Zj), shape2=h0+nj-sum(Zj) )
# ::: lambda2param[j]=rinvgamma( n=1, shape=(c0+n_j)/2, scale=(d0 + sum( (Xj - mean(Xj))^2 ) )/2 )
# ::: X_bar?

#> Analytical form
f0x1 = function(Z) {
  beta(Z+g0, 1-Z+h0)/beta(g0,h0)
}
#> Analytical form
f0x2 = function(X, covar) {
  1/sqrt(pi)*gamma(c0+1)/gamma(c0)*d0^(c0/2)/( (X-mean(covar))^2 + d0 )^((c0+1)/2)
}
#> JOINT 
f0X = numeric(n.train) #% param-free covariate model f0(X,Z)



##### For [Param-free outcome model]: -------------------> MonteCarlo Integration
# Calculate Outcome and Covariate parameter free data model for each observation  
f0S = numeric(n.train) #% param-free outcome model f0(S|X,Z)
E0S = numeric(n.train) #% param-free E(S|X,Z) = Expected value of S|X,Z

M = 1000 # Number of Monte Carlo samples
sumS = numeric(M)
sumES = numeric(M)

for(h in 1:n.train) {
  
  ##> Analytical solution for f0(X, Z) for ### covariate model ###
  # : necessary parameters are supplied by pre-determined hyperparam values...
  f0X[h] = f0x1(Z[h])*f0x2(X=X_tru[h], covar=X_tru)  
  
  ##> Monte Carlo Integration for the ### outcome model ###
  for(j in 1:M) {
    xi_samplej = rt(n=1, df=nu0)                               # bullet on xi
    sig_samplej = rinvgamma(n=1, shape=u0, scale=v0)           # bullet on sig2
    beta_samplej = rmvn(n=1, mu=m0, sigma=SIG_b0*varinf)       # bullet on beta
    
    sumS[j] = dlogsknorm( s=S[h], x=matX_tru[h,], beta=beta_samplej, sig2=sig_samplej, xi=xi_samplej )* 
      dmvn(X=beta_samplej, mu=m0, sigma=SIG_b0*varinf)*        # to joint beta
      dinvgamma(x=sig_samplej, shape=u0, scale=v0)*            # to joint sig2
      dt(x=xi_samplej, df=nu0)                                 # to joint xi
    # 20*exp( betaparamclean[j,1] + betaparamclean[j,2]*cleanx + betaparamclean[j,3]*Z[i] - sig2paramclean[j]/2)*( 1-pnorm(-xiparamclean[j]*sqrt(sig2paramclean[j])/sqrt(xiparamclean[j]^2+1)))
    sumES[j] = 18*exp( sum(matX_tru[h,]*beta_samplej) - sig_samplej/2 )*(1-pnorm(-xi_samplej*sqrt(sig_samplej)/sqrt(xi_samplej^2+1)))
  }
  E0S[h] = sum(sumES)/M  # E[S] value of f0
  f0S[h] = sum(sumS)/M   # density of f0
  
  print(paste("h=",h))
}

par(mfrow=c(1,1))
plot( x=density(f0X) )            # param-free covariate model
plot( x=density(f0S) )            # param-free outcome model
# plot( x=density(E0S) ) # param-free E[outcome] model (wrong)
plot(x=sort(E0S), y=f0S[order(E0S)], type="l")

summary(S)    # original S
summary(E0S)  # E[S] based on f0(S|x) by train data
# - sort(.) : return a vector in ascending order.
# - order(.) : return the index of each element in a vector in sorted order
# - rank(.) : assign a rank to each element in a vector (smallest = 1)


#..............................................................................#
#..............................................................................#
#..............................................................................#
#................................... now ......................................# 
#..............................................................................#
#..............................................................................#
#..............................................................................#
################################################################################
### Step04> Gibbs Sampler --------- cl_membership and param update ---with J= ?
################################################################################
set.seed(1)

total_iter = 30000
# r_convergence = 1000 <- first run Gibbs sampler to determine r value for convergence
r_convergence = 28000

# loglikelihood = numeric(total_iter)                  # for monitor convergence
#loglikelihood = matrix(0, nrow = n.train, ncol = total_iter)
loglikelihood = numeric(total_iter)
# [note] why matrix? each col gives each iteration result

list_piparam = list()        #for Z
list_lambda2param = list()   #for X
list_alphaparam = list()     #for alpha

list_betaparam = list()    #for S
list_sig2param = list()    #for S
list_xiparam = list()      #for S

list_cl = list()

J_prev = J

piparam_prev = piparam
lambda2param_prev = lambda2param
alphaparam_prev = alphaparam

beta_prev = betaparam
sig2_prev = sig2param
xi_prev = xiparam

#### Let's run!!!
for (r in 1:total_iter) {
  ##############################################################################################################
  ##[1] cluster membership computation ----------------------------------------
  clusterout = clusterDPBasic(S, X_tru, Z, cl_membership,
                         piparam, lambda2param,
                         betaparam, sig2param, xiparam, alphaparam,
                         f0X, f0S, u0, v0, m0, SIG_b0, nu0, g0, h0, c0, d0, gamma0, psi0, varinf) #~~~~~Q.David
  
  
  cl_membership = clusterout$cl_membership
  J = length(unique(cl_membership))             #% in case, cluster scenario changes, reflecting them
  
  #piparam = clusterout$piparam
  #lambda2param = clusterout$lambda2param
  
  betaparam = clusterout$betaparam
  sig2param = clusterout$sig2param
  xiparam = clusterout$xiparam
  ##############################################################################################################
  
  #%%%%%%%%%%%%%%%%%%%%%%%% Heavy Computation ON/OFF %%%%%%%%%%%%%%%%%%%%%%%%%%#
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%# 
  # for (i in 1:n.train){
  #   cluster_si = cl_membership     #% current membership vector
  #   
  #   #%%% Before start, check if there is any cluster that has only a single value. 
  #   # if so...then remove the cluster and need to remove parameters as well and re-number clusters
  #   if( sum(cluster_si==cluster_si[i])==1 ) {
  #     j = cluster_si[i]
  #     cluster_si[cluster_si>j] = cluster_si[cluster_si>j] - 1 
  #     
  #     #% for Z
  #     piparam = piparam[-j]
  #     #% for X
  #     lambda2param = lambda2param[-j]
  #     #% for S
  #     sig2param = sig2param[-j]
  #     betaparam = betaparam[-j,]
  #     xiparam = xiparam[-j]
  #     
  #     J = J-1
  #   }
  #   
  #   # a)remove obv and initialize...?
  #   cluster_si[i] = 0                                  #% replace the membership value (the first obv) with "0"!
  #   nj = as.numeric( table(cluster_si[cluster_si>0]) ) #% number of observations in each cluster without observation i
  #   probs = numeric( length(nj)+1 )                    #% for P(s_i=j) ... it's c(31-1,41,28)?? so 3 + 1 ? total cluster number?
  # 
  #   # b)Iterate through each cluster and Calculate probability of staying the same: P(s_i=j)
  #   x_i = c(1, X_tru[i], Z[i])                                        #% c(1, Z, X)
  #   for(j in 1:length(nj)) {
  #     probs[j] = nj[j]/(n.train-1+alphaparam)*
  #       dlogsknorm(s=S[i], x=x_i, beta=betaparam[j,], sig2=sig2param[j], xi=xiparam[j])* 
  #       dbinom(x=x_i[3], size=1, prob = piparam[j])*                # Covariate Z
  #       dnorm(x=x_i[2], mean=mean(X_tru), sd=sqrt(lambda2param[j])) # Covariate X
  #   }
  #   
  #   # c)After iteration through each cluster, Adding the new probability of "Forming a new cluster": P(s_i=J+1) 
  #   probs[j+1] = alphaparam/(n.train-1+alphaparam)*f0S[i]*f0X[i]   #% so it gives...probs: c(prob, prob, prob, 0) -> c(prob, prob, prob, prob)
  # 
  #   # d)Finally, draw a new cluster for each datapoint from a multinomial distribution  #~~~~~~~~~~~~~~~~Q.David
  #   newclust = which( rmultinom(n=1, size=1, prob=probs)==1 ) #% "which(.)" gives the indices of (logic=True)
  #   
  #   #% reflect whats happening above
  #   cl_membership = cluster_si
  # 
  #   cl_membership[i] = newclust  #% assigning the cl_membership to each datapoint. WOW WOW WOW WOW ! 
  #   #% "Multinomial(.)" is the way the datapt accept/reject the new cluster????
  #   #% to face the truth, probabilities are: "probs/sum(probs)"
  #   #% but.."rmultinom(.)" automatically addresses "probs"
  #   
  #   # e)If new cluster is selected by a certain data point, 
  #   # then add new value (obtained from new cluster) to the existing parameter pool
  #   if( length(cl_membership[cl_membership==cl_membership[i]]) == 1 ) {
  #     #% for Z, append
  #     piparam = c( piparam, 
  #                  rbeta(n=1, shape1=g0, shape2=h0) )
  #     #% for X, append
  #     lambda2param = c( lambda2param, 
  #                       rinvgamma(n=1, shape=c0/2, 
  #                                 scale=d0/2) ) 
  # 
  #     beta_p = rmvn(n=1, mu=m0, sigma=SIG_b0)
  #     sig2_p = rinvgamma(n=1, shape=u0, scale=v0)
  #     xi_p = rt(n=1, df=nu0)
  #     betaparam = rbind(betaparam, beta_p)
  #     sig2param = c(sig2param, sig2_p)
  #     xiparam = c(xiparam, xi_p) # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Forget MH! 
  # 
  #     J = J+1
  #   }
  # }
  # J = length(unique(cl_membership))                  #% in case, cluster scenario changes, reflecting them
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%# 
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%# 
 

   

  ###[2] Updating Posterior in accordance with new cl_membership ------------------------------------------
  #% for Z - piparam from beta posterior
  piparam = numeric(J)
  for (j in 1:J) {
    Zj=Z[cl_membership==j]
    nj=length(Zj)
    piparam[j]=rbeta( n=1, shape1=g0+sum(Zj), shape2=h0+nj-sum(Zj) ) #posterior
  }
  #% for X - lambda2param update from MVN, InvGa posterior 
  lambda2param = numeric(J) #; lambda2param
  for (j in 1:J) {
    Xj=X_tru[cl_membership==j]
    Zj = Z[cl_membership==j]
    nj=length(Xj)
    
    lambda2param[j]=rinvgamma( n=1, shape=(nj+c0)/2, scale=( sum((Xj - mean(X_tru))^2) + d0 )/2 )
  }
  
  #% for alphaparam
  eta = rbeta(n=1, shape1=alpha0+1, shape2=n.train)
  pi_eta = (gamma0+J-1)/(gamma0+J-1 + n.train*(psi0-log(eta)))
  alphaparam = pi_eta*rgamma( n=1, shape=gamma0+J, 
                              rate=psi0-log(eta) ) + (1-pi_eta)*rgamma( n=1, shape=gamma0+J-1, rate=psi0-log(eta) )
  
  #% for S...based on the C++ results
  beta_old = betaparam
  sig2_old = sig2param
  xi_old = xiparam
  
  for(j in 1:J) {
    # Sample proposals from priors
    beta_p = rmvn(n=1, mu=m0, sigma=SIG_b0)
    sig2_p = rinvgamma(n=1, shape=u0, scale=v0)
    xi_p = rt(n=1, df=nu0)
    
    # subsetting by cluster
    Sj = S[cl_membership==j]
    if(length(Sj)==1) {
      matXj = matrix( matX_tru[cl_membership==j, ], nrow=1 ) # It is a vector added to the matrix
    } else {
      matXj = matX_tru[cl_membership==j, ]                   # It is a regular matrix
    } 
    
    #matXj = matX_tru[cl_membership==j, ] 
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # > MH algorithm (all together)
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # numerator=0
    # denominator=0
    # for (i in 1:length(Sj)){
    #   numerator = numerator + dlogsknorm_log(Sj[i], matXj[i,], 
    #                                          t(as.matrix(beta_p)), sig2_p, xi_p)
    #   denominator = denominator + dlogsknorm_log(Sj[i], matXj[i,], 
    #                                              t(as.matrix(beta_old[j,])), sig2_old[j], xi_old[j])
    #   
    #   ratio = min(exp(numerator-denominator), 1)
    #   U = runif(n = 1, min = 0, max = 1)
    #   if(U < ratio) {
    #     betaparam[j, ] = beta_p
    #     sig2param[j] = sig2_p
    #     xiparam[j] = xi_p
    #   } else {
    #     betaparam[j, ] = beta_old[j, ]
    #     sig2param[j] = sig2_old[j]
    #     xiparam[j] = xi_old[j]
    #   }
    # }
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # > MH algorithm (separate)
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    #>>> - Sifting [betaparam]_j #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    numerator=0
    denominator=0
    for (i in 1:length(Sj)){
      numerator = numerator +
        dlogsknorm_log(Sj[i], matXj[i,], t(as.matrix(beta_p)), sig2_old[j], xi_old[j]) + 0
      #dmvn(beta_p, mu=m0.st, sigma=SIG_b0.st*varinf) + 0
      #dmvn(beta_prev.st[j,], mu=m0.st, sigma=SIG_b0.st*varinf)
      
      denominator = denominator +
        dlogsknorm_log(Sj[i], matXj[i,], t(as.matrix(beta_old[j,])), sig2_old[j], xi_old[j]) + 0
      #dmvn(beta_prev.st[j,], mu=m0.st, sigma=SIG_b0.st*varinf) + 0
      #dmvn(beta_p, mu=m0.st, sigma=SIG_b0.st*varinf)
    }
    
    #> [CHECK.1] ratio is too small???
    ratio_beta = min(exp(numerator-denominator), 1)
    U = runif(n = 1, min = 0, max = 1)
    if(ratio_beta > U) {
      betaparam[j, ] = beta_p
    } else {
      betaparam[j,] = beta_old[j,]
    }
    
    #>>> - Sifting [sig2param]_j #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    numerator=0
    denominator=0
    for (i in 1:length(Sj)){
      numerator = numerator +
        dlogsknorm_log(Sj[i], matXj[i,], t(as.matrix(betaparam[j,])), sig2_p, xi_old[j]) + 0
      #log(dinvgamma(sig2_p, shape=u0.st, scale=v0.st)) +
      #log(dinvgamma(sig2_prev.st[j], shape=u0_new.st, scale=v0_new.st))
      
      denominator = denominator +
        dlogsknorm_log(Sj[i], matXj[i,], t(as.matrix(betaparam[j,])), sig2_old[j], xi_old[j]) + 0
      #log(dinvgamma(sig2_prev.st[j], shape=u0.st, scale=v0.st)) +
      #log(dinvgamma(sig2_p, shape=u0_new.st, scale=v0_new.st))
    }
    
    #> [CHECK.2] ratio is too small???
    ratio_sig2 = min(exp(numerator-denominator), 1)
    U = runif(n = 1, min = 0, max = 1)
    if(U < ratio_sig2) {
      sig2param[j] = sig2_p
    }
    else {
      sig2param[j] = sig2_old[j]
    }
    
    #>>> - Sifting [xiparam]_j #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    numerator=0
    denominator=0
    for (i in 1:length(Sj)){
      numerator = numerator +
        dlogsknorm_log(Sj[i], matXj[i,], t(as.matrix(betaparam[j,])), sig2param[j], xi_p) + 0
      #dt(xi_p, df=nu0.st, log = TRUE) + 0
      #dt(xi_prev.st[j], df=nu0.st, log = TRUE)
      
      denominator = denominator +
        dlogsknorm_log(Sj[i], matXj[i,], t(as.matrix(betaparam[j,])), sig2param[j], xi_old[j]) + 0
      #dt(xi_prev.st[j], df=nu0.st, log = TRUE) + 0
      #dt(xi_p, df=nu0.st, log = TRUE)
    }
    
    #> [CHECK.3] ratio is too small???
    ratio_xi = min(exp(numerator-denominator), 1)
    U = runif(n = 1, min = 0, max = 1)
    if(U < ratio_xi) {
      xiparam[j] = xi_p
    }
    else {
      xiparam[j] = xi_old[j]
    }
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
  } #%%% END of [j] %%% 
    

  
  ###[3] Calculating Loglikelihood ---------------------------------------------
  loglike_r = 0
  
  for(i in 1:n.train) {
    x_i = matX_tru[i, ]
    j = cl_membership[i]
    
    loglike_r = loglike_r +                                                     
      dlogsknorm_log(s=S[i], x=x_i, beta=betaparam[j,], sig2=sig2param[j], xi=xiparam[j]) +
      dbinom(x = x_i[3], size = 1, prob = piparam[j], log = TRUE) + 
      dnorm(x = x_i[2], mean = mean(X_tru), sd = sqrt(lambda2param[j]), log = TRUE)
  }
  loglikelihood[r] = loglike_r # cumulative sum
  
  if(r > r_convergence) {
    list_piparam[[r-r_convergence]] = piparam   #for Z
    
    list_lambda2param[[r-r_convergence]] = lambda2param  #for X
    
    list_alphaparam[[r-r_convergence]] = alphaparam     #for alpha
    
    list_betaparam[[r-r_convergence]] = betaparam    #for S
    list_sig2param[[r-r_convergence]] = sig2param    #for S
    list_xiparam[[r-r_convergence]] = xiparam      #for S
    
    list_cl[[r-r_convergence]] = cl_membership
  }
  print(paste("r=",r)) #% iteration progress
}

#plot(apply(loglikelihood, 2, sum), type="l")
plot(loglikelihood, type="l")
########################################################

table(cl_membership) # final number of cl...
lppd_EX_DP.S <- mean( loglikelihood ); lppd_EX_DP.S # -22207.38


#..............................................................................#
#..............................................................................#
############################### Predictions ####################################
#..............................................................................#
#..............................................................................#

################################################################################
###> Training data -------------------------------------------------------------
################################################################################
#% Ultimate liquid weight vector for each iteration
W_paramFree.vec = matrix(0, nrow = n.train, ncol = total_iter-r_convergence) 
density.train = matrix(0, nrow = n.train, ncol = total_iter-r_convergence)
expval.train = matrix(0, nrow = n.train, ncol = total_iter-r_convergence)

#% prediction for each iteration
for(r in 1:(total_iter-r_convergence)) {
  piparamclean = list_piparam[[r]]
  lambda2paramclean = list_lambda2param[[r]] #for Xstar
  alphaparamclean = list_alphaparam[[r]]      #for alpha
  
  betaparamclean = list_betaparam[[r]]      #for S
  sig2paramclean = list_sig2param[[r]]      #for S
  xiparamclean = list_xiparam[[r]]       #for S
  
  cl_membership = list_cl[[r]]
  
  J = nrow(betaparamclean)
  f_S = numeric(J)                  # density for S
  
  #% Ultimate solid weight matrix
  W_paramBase.matrix = matrix(0, nrow=n.train, ncol=J)

  for(h in 1:n.train) {
    
    fX_j = numeric(J)       # discrete covariate model for each j
    w_j.solid = numeric(J)  # discrete weight for each j  
    E_value = numeric(J)    # to save E[S|X] for each j
    
    for(j in 1:J) {
      # discrete covariate joint
      fX_j[j] = dbinom(x=Z[h], size=1, prob=piparamclean[j])*
        dnorm(x=X_tru[h], mean=mean(X_tru), sd=sqrt(lambda2paramclean[j]))
      
      #% solid weight ***** component
      w_j.solid[j] = length(S[cl_membership==j])/(alphaparamclean + n.train)*fX_j[j]
      
      #% compute pred.density
      f_S[j] = dlogsknorm(s=S[h], x=c(1, X_tru[h], Z[h]), 
                              beta=betaparamclean[j,], sig2=sig2paramclean[j], xi=xiparamclean[j])
      #% compute pred.value
      E_value[j] = 4*exp(sum( c(1, X_tru[h], Z[h])*betaparamclean[j,]) - sig2paramclean[j]/2)*
        (1-pnorm(-xiparamclean[j]*sqrt(sig2paramclean[j])/sqrt(xiparamclean[j]^2+1)))
    }
    
    #% liquid weight ***** component
    WJ1 = alphaparamclean/(alphaparamclean + n.train)*f0X[h] #; print(paste("liquid weight=", WJ1))
    
    #% Collecting values and constructing weights for each iteration
    W_paramFree.vec[h,r] = WJ1/(WJ1 + sum(w_j.solid))
    W_paramBase.matrix[h, ] = w_j.solid/(WJ1 + sum(w_j.solid)) 
    
    #% Collecting values for predictive density (via weighted AVG)
    density.train[h,r] = W_paramFree.vec[h,r]*f0S[h] + sum(W_paramBase.matrix[h,]*f_S)
    #% Collecting values for prediction (via weighted AVG)
    expval.train[h,r] = W_paramFree.vec[h,r]*E0S[h] + sum(W_paramBase.matrix[h,]*E_value)
  }
  print(paste("r=", r))
}

# AVG prediction
expS.DP.Avg.train <- apply(X=expval.train, MARGIN=1, FUN=mean)

summary(S)
summary(expS.DP.Avg.train)
summary(log(S))
summary(log(expS.DP.Avg.train))
#> summary(log(S))
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#4.277   9.282  10.690  10.713  11.955  16.719 
#> summary(log(expS.DP.train))
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#7.080   9.978  10.766  11.117  11.999  16.816 


format( aggregate(log(S), by=list(cl_membership), FUN=summary), scientific=FALSE)
format( aggregate(log(expS.DP.Avg.train), by=list(cl_membership), FUN=summary), scientific=FALSE)



################################################################################
###> Testing data -------------------------------------------------------------
################################################################################ 
###> Parameter-free component
#> Analytical form
f0x1 = function(Z) {
  beta(Z+g0, 1-Z+h0)/beta(g0,h0)
}
#> Analytical form
f0x2 = function(X, covar) {
  1/sqrt(pi)*gamma(c0+1)/gamma(c0)*d0^(c0/2)/( (X-mean(covar))^2 + d0 )^((c0+1)/2)
}
#> JOINT 
f0X.test = numeric(n.test) #% param-free covariate model f0(X,Z)

##### For [Param-free outcome model]: -------------------> MonteCarlo Integration
# Calculate Outcome and Covariate parameter free data model for each observation  
f0S.test = numeric(n.test) #% param-free outcome model f0(S|X,Z)
E0S.test = numeric(n.test) #% param-free E(S|X,Z) = Expected value of S|X,Z

M = 1000 # Number of Monte Carlo samples
sumS = numeric(M)
sumES = numeric(M)

for(h in 1:n.test) {
  
  ##> Analytical solution for f0(X, Z) for ### covariate model ###
  # : necessary parameters are supplied by pre-determined hyperparam values...
  f0X.test[h] = f0x1(Z.test[h])*f0x2(X_tru.test[h], X_tru.test)  
  
  ##> Monte Carlo Integration for the ### outcome model ###
  for(j in 1:M) {
    xi_samplej = rt(n=1, df=nu0)                               # bullet on xi
    sig_samplej = rinvgamma(n=1, shape=u0, scale=v0)           # bullet on sig2
    beta_samplej = rmvn(n=1, mu=m0, sigma=SIG_b0*varinf)       # bullet on beta
    
    sumS[j] = dlogsknorm( s=S.test[h], x=matX_tru.test[h,], beta=beta_samplej, sig2=sig_samplej, xi=xi_samplej )* 
      dmvn(X=beta_samplej, mu=m0, sigma=SIG_b0*varinf)*        # to joint beta
      dinvgamma(x=sig_samplej, shape=u0, scale=v0)*            # to joint sig2
      dt(x=xi_samplej, df=nu0)                                 # to joint xi
    # 20*exp( betaparamclean[j,1] + betaparamclean[j,2]*cleanx + betaparamclean[j,3]*Z[i] - sig2paramclean[j]/2)*( 1-pnorm(-xiparamclean[j]*sqrt(sig2paramclean[j])/sqrt(xiparamclean[j]^2+1)))
    sumES[j] = 18*exp( sum(matX_tru.test[h,]*beta_samplej) - sig_samplej/2 )*(1-pnorm(-xi_samplej*sqrt(sig_samplej)/sqrt(xi_samplej^2+1)))
  }
  E0S.test[h] = sum(sumES)/M  # E[S] value of f0
  f0S.test[h] = sum(sumS)/M   # density of f0
  
  print(paste("h=",h))
}
summary(S.test)    # original S
summary(E0S.test)  # E[S] based on f0(S|x) by train data





###> NOW....
# cl_membership.test <- full.test.sweden$Risk
# table(cl_membership.test)                    # this is wrong!!!!!!!!!!!!!!!!!

#% Ultimate liquid weight vector for each iteration
W_paramFree.vec.test = matrix(0, nrow = n.test, ncol = total_iter-r_convergence) 
density.test = matrix(0, nrow = n.test, ncol = total_iter-r_convergence)
expval.test = matrix(0, nrow = n.test, ncol = total_iter-r_convergence)


for(r in 1:(total_iter-r_convergence)) {
  piparamclean = list_piparam[[r]]
  lambda2paramclean = list_lambda2param[[r]] #for Xstar
  alphaparamclean = list_alphaparam[[r]]      #for alpha
  
  betaparamclean = list_betaparam[[r]]      #for S
  sig2paramclean = list_sig2param[[r]]      #for S
  xiparamclean = list_xiparam[[r]]       #for S
  
  cl_membership.test = list_cl[[r]]

  
  J = nrow(betaparamclean)
  
  
  f_S.test = numeric(J)                  # density for S
  
  #% Ultimate solid weight matrix
  W_paramBase.matrix.test = matrix(0, nrow=n.test, ncol=J)
  
  for(h in 1:n.test) {
    
    fX_j = numeric(J)       # discrete covariate model for each j
    w_j.solid = numeric(J)  # discrete weight for each j  
    E_value.test = numeric(J)    # to save E[S|X] for each j
    
    for(j in 1:J) {
      # discrete covariate joint
      fX_j[j] = dbinom(x=Z.test[h], size=1, prob=piparamclean[j])*
        dnorm(x=X_tru.test[h], mean=mean(X_tru.test), sd=sqrt(lambda2paramclean[j]))
      
      #% solid weight ***** component #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Q.David. How to determine "cl_membership" for test set?
      #w_j.solid[j] = length(S.test[cl_membership.test==j])/(alphaparamclean + n.test)*fX_j[j]
      w_j.solid[j] = sum(cl_membership==j) / (alphaparamclean + n.train)*fX_j[j] #~~~~~~~~~~~~~~~~~~~~~~~~~~~~QQQ
      #% compute pred.density
      f_S.test[j] = dlogsknorm(s=S.test[h], x=c(1, X_tru.test[h], Z.test[h]), 
                          beta=betaparamclean[j,], sig2=sig2paramclean[j], xi=xiparamclean[j])
      #% compute pred.value
      E_value.test[j] = 4*exp(sum( c(1, X_tru.test[h], Z.test[h])*betaparamclean[j,]) - sig2paramclean[j]/2)*
        (1-pnorm(-xiparamclean[j]*sqrt(sig2paramclean[j])/sqrt(xiparamclean[j]^2+1)))
    }
    
    #% liquid weight ***** component
    #WJ1 = alphaparamclean/(alphaparamclean + n.test)*f0X[h] #; print(paste("liquid weight=", WJ1))
    WJ1 = alphaparamclean/(alphaparamclean + n.train)*f0X[h] #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~QQQ
    #% Collecting values and constructing weights for each iteration
    W_paramFree.vec.test[h,r] = WJ1/(WJ1 + sum(w_j.solid))
    W_paramBase.matrix.test[h, ] = w_j.solid/(WJ1 + sum(w_j.solid)) 
    
    #% Collecting values for predictive density (via weighted AVG)
    density.test[h,r] = W_paramFree.vec.test[h,r]*f0S.test[h] + sum(W_paramBase.matrix.test[h,]*f_S.test)
    #% Collecting values for prediction (via weighted AVG)
    expval.test[h,r] = W_paramFree.vec.test[h,r]*E0S.test[h] + sum(W_paramBase.matrix.test[h,]*E_value.test)
  }
  print(paste("r=", r))
}




###----------------------------- Visual S -----------------------------------###

expval.train
density.train

par(mfrow=c(1,1))
plot(density(expval.train[,1])) # pointless
plot(density(expval.train[,2000])) # pointless

log_expS.DP.Avg.train = log( apply(X=expval.train, MARGIN=1, FUN=mean) )

### turn a matrix into a dataframe of 100 scenarios of E[Y|X] and their densities
den_df <- data.frame(density.train); head(den_df,2) 
val_df <- data.frame(expval.train); head(val_df,2) 

summary(val_df[, 1])    # 1st scenario 
summary(val_df[, 2000]) # 2000th scenario
summary(val_df[, 1950:2000])


hist(log(S), breaks=n.breaks.train, freq=F, xlab ="log(S)", main="Predictive loss density for a policy", 
     col="white", ylim=c(0,0.4))
for (i in 1:50) {
  lines(density( log(val_df[, i]) ), lwd=0.01, col="black" )
}
for (i in 100:150) {
  lines(density( log(val_df[, i]) ), lwd=0.01, col="grey" )
}
for (i in 200:250) {
  lines(density( log(val_df[, i]) ), lwd=0.01, col="grey" )
}
for (i in 1:2000) {
  lines(density( log(val_df[, i]) ), lwd=0.01, col="grey" )
}

lines(density(log_expS.DP.Avg.train), col="red", lwd=5)



# expval.test
# density.test
# 
# par(mfrow=c(1,1))
# plot(density(expval.test[,1])) # pointless
# plot(density(expval.test[,2000])) # pointless
# 
# log_expS.DP.Avg.test = log( apply(X=expval.test, MARGIN=1, FUN=mean) )
# 
# ### turn a matrix into a dataframe of 100 scenarios of E[Y|X] and their densities
# den_df.test <- data.frame(density.test); head(den_df.test,2) 
# val_df.test <- data.frame(expval.test); head(val_df.test,2) 
# 
# summary(val_df.test[, 1])    # 1st scenario 
# summary(val_df.test[, 2000]) # 2000th scenario
# summary(val_df.test[, 1950:2000])
# 
# 
# hist(log(S.test), breaks=n.breaks.test, freq=F, xlab ="log(S)", main="Predictive loss density for a policy", 
#      col="white", ylim=c(0,0.4))
# for (i in 1:50) {
#   lines(density( log(val_df.test[, i]) ), lwd=0.01, col="black" )
# }
# 
# for (i in 1:2000) {
#   lines(density( log(val_df.test[, i]) ), lwd=0.01, col="grey" )
# }
# 
# lines(density(log_expS.DP.Avg.test), col="red", lwd=5)














################################################################################
#- Out-of sample----------------------------------------------------------------
################################################################################
####> [AVG.vec Scenario] 
den.DP.Avg.test <- apply(X=density.test, MARGIN=1, FUN=mean)
expS.DP.Avg.test <- apply(X=expval.test, MARGIN=1, FUN=mean)

summary(S.test)
summary(expS.DP.Avg.test)
summary(log(S.test))
summary(log(expS.DP.Avg.test))

SSPE.DP.LSN <- sum( (log(expS.DP.Avg.test) - log(S.test))^2 ); SSPE.DP.LSN    #269.1253
SAPE.DP.LSN <- sum( abs(log(expS.DP.Avg.test) - log(S.test)) ); SAPE.DP.LSN   #196.2275
lppd_EX_DP.S                                                                  #-22207.38

par(mfrow=c(1,2))
hist(S.test, freq=F, xlab ="S_h", col="white", breaks=n.breaks.test*10) #, main="Predictive total loss density for a policy") 
lines(density(expS.DP.Avg.test), col="red", lwd=2)   # extremely skewed..so not meaningful
hist(log(S.test), freq=F, xlab ="log(S_h)", col="white", breaks=n.breaks.test*2) #, main="Predictive total loss density for a policy") 
lines(density(log(expS.DP.Avg.test)), col="red", lwd=2) 




##### 95% credible interval
compute_CI <- function(samples, alpha = 0.05) {
  quantile(samples, probs = c(alpha / 2, 1 - alpha / 2))
}


#> beta0
list_betaparam0 <- numeric(2000)
for (i in 1:2000){
  list_betaparam0[i] <- mean( list_betaparam[[i]][, 1] ) 
}
summary(list_betaparam0)# 7.076

compute_CI(list_betaparam0, alpha = 0.05)[1]
#    2.5% 
# 7.017344 
compute_CI(list_betaparam0, alpha = 0.05)[2]
#    97.5% 
# 7.141527


#> beta1
list_betaparam1 <- numeric(2000)
for (i in 1:2000){
  list_betaparam1[i] <- mean( list_betaparam[[i]][, 2] ) 
}
summary(list_betaparam1)# 0.8355  

compute_CI(list_betaparam1, alpha = 0.05)[1]
#    2.5% 
# 0.8242637  
compute_CI(list_betaparam1, alpha = 0.05)[2]
#    97.5% 
# 0.8464996 


#> beta2
list_betaparam2 <- numeric(2000)
for (i in 1:2000){
  list_betaparam2[i] <- mean( list_betaparam[[i]][, 3] ) 
}
summary(list_betaparam2) #-0.4604 

compute_CI(list_betaparam2, alpha = 0.05)[1]
#    2.5% 
# -0.5007453  
compute_CI(list_betaparam2, alpha = 0.05)[2]
#    97.5% 
# -0.4182603  


#> sig2
list_sig2test <- numeric(2000)
for (i in 1:2000) {
  list_sig2test[i] <- mean( list_sig2param[[i]] )
}
summary(list_sig2test ) #4.0459   

compute_CI(list_sig2test, alpha = 0.05)[1]
#    2.5% 
# 0.8985093   
compute_CI(list_sig2test, alpha = 0.05)[2]
#    97.5% 
# 12.33876  


#> xi
list_xitest <- numeric(2000)
for (i in 1:2000) {
  list_xitest[i] <- mean( list_xiparam[[i]] )
}
summary(list_xitest ) #-0.06749 

compute_CI(list_xitest, alpha = 0.05)[1]
#    2.5% 
# -0.6932471   
compute_CI(list_xitest, alpha = 0.05)[2]
#    97.5% 
# 0.6332571  



# --------------------------------------------------------------------------------- [END of the regular models]
#-------------------------- Experiments start here! ---------------------------#
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
########################### with Dirty data [1%] ###############################
################################################################################
### Step00> Define prior, hyperprior model and hyperparameter values
################################################################################
# ----------- Outcome ----- # ::: for S ~ LSN( X*betaparam, sig2param, xiparam )

# PRIOR: "betaparam" ~ MVN(m0, SIG_b0)

#> Gaussian regression for initialize OUTCOME model parameter m0, SIG_b0
fit1.xx <- glm( S ~  Xstar + factor(Z), family=Gamma(link = "log") )

options(scipen = 999)
summary(fit1.xx)


m0.xx = coef(fit1.xx)    # 1x3 initial reg_coeff vector (sort of "means"): Regression result as "mean"
SIG_b0.xx = vcov(fit1.xx); SIG_b0.xx   # 3x3 initial cov matrix of reg_coeff
p = length(m0.xx)
options(scipen = 999)
SIG_b0inv.xx = solve(a=SIG_b0.xx) # inverse of cov matrix of reg_coeff for later use (for posterior on reg_coeff)!!

#****** for VAR Inflation factor
delta = 0.01             # equal importance on m0 and beta ????
varinf = n.train/100 # rule of thumb divide by 5



###> PRIOR: "sig2param" ~ IG( u0, v0 )
# Hyperprior               u0 ~ Fink( rho_u1, rho_u2 )
# Hyperprior                   v0 ~ Ga( rho_bv1, rho_v2 )

u0.xx=1.7591
v0.xx=7.113


# ## [CHECK later] ---------------------------------------------------------------
# # How can we design proposal for sig2 ~ InvGa(.) ?
# par(mfrow=c(1,1))
# plot( density(unlist(list_sig2param)) )                                   # sig2 we have...before
# curve( dinvgamma(x, shape=mean(unlist(list_u0param)), 
#                  scale=mean(unlist(list_v0param))), add=TRUE, col="red") # sig2 we obtain....after
# curve( dinvgamma(x, shape=u0_new, 
#                  scale=v0_new), add=TRUE, col="blue")                    # sig2 proposal adjustment
# 
# table(unlist(list_sig2param))
# sort(table(unlist(list_sig2param)),DESC=TRUE)                             # to see if they are appropriate...
# #-------------------------------------------------------------------------------

###> PRIOR: "xiparam" ~ t(loc, nu0, sca)
#loc = 0   # ? Location hyperparameter for T distribution
#nu0 = 1/2 # df for T distribution on pigtail parameter
#sca = 5   # ? Scale hyperparameter for T distribution 
nu0.xx = n.train-1
#- df refers to the number of independent observations (N-1):

# ------------------------------------------------------------------------------ #



# ::::::: # ------------------------ PRIOR ------------------------- # ::::::: #
# ------- Covariate ------ ::: for Z~Bin( 1, [piparam] ) where 
#                                             "piparam" ~ Beta(g0, h0) 
g0.xx = 0.5 
h0.xx = 0.5   

#] ------ Covariate ------ ::: for X ~ N( X_bar, [lambda2param] ) where 
# PRIOR                                            "lambda2param" ~ IG(c0, d0)
#------------------------------------------------------------------------------------------------------------------------
#X_bar = mean(X_tru) #????

#> lambda2
c0.xx = 0.5          
d0.xx = 0.5

# # ---- Covariate model ----# ::: for X|Z ~ N( K0+K1*Z, [lambda2param] ) where 
# #                                                        "lambda2param" ~ IG(c0, d0) 
# fit2 <- lm(X ~ Z)
# summary(fit2)
# 
# kappa_v0 <- coef(fit2)
# SIG_k0 <- vcov(fit2)*varinf
# SIG_inv_k0 <- solve(SIG_k0)
# 
# c0 = 0.5          
# d0 = 0.5    
# 
# KK1 = cbind(1, Z) # for the analytical solution for Parameter-Free covariate model development later
# # [note] 
# # - we need KK1 matrix for all obvs w/o concerning cluster
# # - KK1_j : cluster-wise 
# # - KK2_j         : cluster-wise


# ::::::: # ------------------------ PRIOR ------------------------- # ::::::: #
# ----------- precision -----------# :::  for alphaparam ~ Ga( gamma0, psi0 )
gamma0.xx = 1    #                                          
psi0.xx = 1    



################################################################################
### Step01> initialize cluster membership - Kmeans clustering
################################################################################
# J=3
# set.seed(1)
# cluster = kmeans(cbind(S,X_tru,Z), J)
# cl_membership = cluster$cluster 
# table(cl_membership)

cl_membership.xx <- full.train.sweden$Risk
table(cl_membership.xx)

J=5





################################################################################
### Step02> initialize major parametersssss by clusters using your "POSTERIOR"
#                                                 Based on kmean clustering: J=3
################################################################################
########################--- for Covariates + alpha ---##########################
################################################################################

##[A] for Z ~ Bin( 1, [piparamS_j] )
#                     "piparamS_j" ~ Beta(g0, h0)
#-------------------------------------------------------------------------------
# Covariate parameter from Posterior
piparam.xx = numeric(J)

#> 1. Regular
for (j in 1:J) {
  Zj=Z[cl_membership.xx==j]
  nj=length(Zj)
  piparam.xx[j] = rbeta( n=1, shape1=g0.xx+sum(Zj), shape2=h0.xx+nj-sum(Zj) ) 
}
# #> 2. MAR? For Z ~ Bin( n, "prob_j" )
# Z.miss = is.na(Z)
# for (j in 1:J) {
#   Zj=Z[cl_membership==j & !Z.miss]
#   nj=length(Zj)
#   if(nj > 0) { #sample from posterior
#     piparam[j] = rbeta( n=1, shape1=g0+sum(Zj), shape2=h0+nj-sum(Zj) )
#   }
#   else { #sample from prior
#     piparam[j] = rbeta(n=1, shape1 = g0, shape2 = h0)
#   }
# }
piparam.xx #% pi parameter sampled from posterior

#% [Check!] empirical pi ? using "aggregate(.)": investigation by splitting the data into subset
pi_empirical = aggregate(x=Z, by=list(cl_membership.xx), FUN=mean, na.rm= T)$x; pi_empirical



##[B] for X ~ N( X_bar, [lambda2param] ) where 
# PRIOR                 "lambda2param" ~ IG(c0, d0)
#-------------------------------------------------------------------------------
# Covariate parameter from Posterior
lambda2param.xx = numeric(J)

#> 1. Regular
for (j in 1:J) {
  Xj=Xstar[cl_membership.xx==j]
  n_j=length(Xj)
  lambda2param.xx[j]=rinvgamma( n=1, shape=(c0.xx+n_j)/2, scale=(d0.xx + sum( (Xj - mean(Xj))^2 ) )/2 )  
}                                                                                                     
lambda2param.xx #% lambda2 parameter (posterior sample)
#% empirical X_barF ? using "aggregate(.)": investigation by splitting the data into subset
X_bar.xx = aggregate(x=Xstar, by=list(cl_membership.xx), FUN=mean, na.rm= T)$x; X_bar.xx
# ---------------------------------------------------------------- Perfect!!!!!!


# #> 2. NDB? for X|Z ~ N( "kappa0_j"+"kappa1_j"*Z, [lambda2param] ) 
# # Sample it from Normal(for kappa), ............ "lambda2param" ~ IG(c0, d0)
# kappaparam = matrix(nrow = J, ncol = length(kappa_v0)); kappaparam # vectors for each cluster
# lambda2param = numeric(J); lambda2param
# # [note] 
# # - kappa_v0: initial intercept, initial slope
# # - kappa0_j from kappaparam : intercepts for each cluster
# # - kappa1_j from kappaparam : slope for each cluster
# for (j in 1:J) {
#   Xj = X[cl_membership==j]
#   Zj = Z[cl_membership==j]
#   nj = length(Xj)
#   lambda2param[j]=rinvgamma( n=1, shape=(nj+c0)/2, scale=0.5*( sum((Xj - kappa_v0[1]-kappa_v0[2]*Zj)^2) + d0 ) )
#   
#   KK1j = cbind(1, Zj)
#   KK2j = matrix(c(sum(Xj), sum(Xj*Zj)), nrow =2, ncol = 1) # 2x1 vector
#   
#   SIG_inv_j = solve(SIG_inv_k0+(t(KK1j) %*% KK1j)) # solve( ..[matrix].. ) = [matrix]^(-1)
#   kappaparam[j,] <- rmvn(n=1, mu = SIG_inv_j %*% (SIG_inv_k0 %*% kappa_v0 + KK2j), sigma = lambda2param[j]*SIG_inv_j)      
# }
# kappaparam   #% kappa parameter sampled from posterior
# lambda2param #% lambda2 parameter sampled from posterior
# 
# #% [Check!] empirical kappa, lambda ? using "aggregate(.)": investigation by splitting the data into subset
# kappa_empirical = rbind(coef(lm(X[cl_membership==1] ~ Z[cl_membership==1])),
#                       coef(lm(X[cl_membership==2] ~ Z[cl_membership==2])),
#                       coef(lm(X[cl_membership==3] ~ Z[cl_membership==3]))); kappa_empirical
# lambda_empirical = aggregate(x=as.numeric(X), by=list(cl_membership), FUN=var, na.rm=T)$x; lambda_empirical



##[C] for alphaparam ~ Ga( shape=gamma0, rate=psi0 ) .. starting with J clusters for Mixed Gamma POSTERIOR
alpha0.xx=2 # initialization: typically..1,2...
eta.xx = rbeta(n=1, shape1=alpha0.xx+1, shape2=n.train); eta.xx
pi_eta.xx = (gamma0.xx+J-1)/(gamma0.xx+J-1 + n.train*(psi0.xx-log(eta.xx))); pi_eta.xx

# precision parameters and its done!
alphaparam.xx = pi_eta.xx*rgamma( n=1, shape=gamma0.xx+J, 
                                  rate=psi0.xx-log(eta.xx) ) + (1-pi_eta.xx)*rgamma( n=1, shape=gamma0.xx+J-1, rate=psi0.xx-log(eta.xx) ); alphaparam.xx




################################################################################
#############################--- for outcome ---################################
################################################################################

##[D] ------ Outcome S ------ ::: [betaparam], [sig2param], [xiparam].... try MH-algorithm
# - No conjugacy...so prepare MM sampling

# This is where new parameters go.
betaparam.xx = matrix(data=NA, nrow=J, ncol=length(m0.xx)) 
sig2param.xx = numeric(J) 
xiparam.xx = numeric(J)

# Set in place the old parameters..sample from prior...now..it's temporary values
beta_prev.xx = matrix( rep(rmvn(n=1, mu=m0.xx, sigma=SIG_b0.xx*varinf), J), nrow=J, byrow= TRUE )
sig2_prev.xx = rep(rinvgamma(n=1, shape=u0.xx, scale=v0.xx), J)
xi_prev.xx = rep(rt(n=1, df=nu0.xx), J) # imagine we already have them...old days..


#### Metropolis Hastings to update betaparam, sig2param, xiparam
for(j in 1:J) {
  # Sample proposals from priors
  beta_prop.xx = rmvn(n=1, mu=m0.xx, sigma=SIG_b0.xx*varinf)
  sig2_prop.xx = rinvgamma(n=1, shape=u0.xx, scale=v0.xx)
  xi_prop.xx = rt(n=1, df=nu0.xx)
  
  # subsetting by cluster
  Sj = S[cl_membership.xx==j]
  matXj = matXstar[cl_membership.xx==j, ]
  
  # #>>>>>>>>>>>>>>>>>>>> if you have MAR.....impute first! <<<<<<<<<<<<<<<<<<<#
  # Z.miss_j = Z.miss[cl_membership==j]
  # miss.indx = which(Z.miss_j)
  # # p0: joint where x1 = 0
  # # p1: joint where x1 = 1
  # for(i in miss.indx) {
  #   if (Sj[i]>1) {
  #     p0log = 
  #       dlogsknorm_log(y=Sj[i], 
  #                      x=t(as.matrix(c(matXj[i,1], 0, matXj[i,3]), nrow=1)), 
  #                      beta=beta_prev[j,], 
  #                      sig2=sig2_prev[j], 
  #                      xi=xi_prev[j]) + 
  #       dbinom(x=0, size=1, piparam[j], log=TRUE) + 0
  #       #log( 1-sigmoid( sum(c(matXj[i,1], 0, matXj[i,3])*betat_old[j, ]) ) )
  #     
  #     p1log = 
  #       dlogsknorm_log(y=Sj[i], 
  #                      x=t(as.matrix(c(matXj[i,1], 1, matXj[i,3]), nrow=1)), 
  #                      beta=beta_prev[j,], 
  #                      sig2=sig2_prev[j], 
  #                      xi=xi_prev[j]) + 
  #       dbinom(x=1, size=1, piparam[j], log=TRUE) + 0
  #       #log( 1-sigmoid( sum(c(matXj[i,1], 1, matXj[i,3])*betat_old[j, ]) ) )
  #   }
  #   else {
  #     p0log = 
  #       dbinom(x=0, size=1, piparam[j], log=TRUE) + 0 # (1 - pi_j)
  #       #log( sigmoid( sum(c(matXj[i,1], 0, matXj[i,3])*betat_old[j, ]) ) )
  #     
  #     p1log = 
  #       dbinom(x=1, size=1, piparam[j], log=TRUE) + 0 # (pi_j)
  #       #log( sigmoid( sum(c(matXj[i,1], 1, matXj[i,3])*betat_old[j, ]) ) )
  #   }
  #   
  #   # Let's impute!!!
  #   matXj[i,2] = rbinom(n=1, size=1, prob=1/(1+exp(p0log-p1log))) #imputing NA with using posterior Pi 
  # }
  
  ##> A. MH algorithm for: [betaparam]
  #-Note: "beta_prop" is from PRIOR proposal (new) while "beta_prev[j,]" is from the single POSTERIOR sample (old). 
  numerator=0
  denominator=0
  for (i in 1:length(Sj)){
    numerator = numerator + 
      dlogsknorm_log(Sj[i], matXj[i,], t(as.matrix(beta_prop.xx)), sig2_prop.xx, xi_prop.xx)
    denominator = denominator + 
      dlogsknorm_log(Sj[i], matXj[i,], t(as.matrix(beta_prev.xx[j,])), sig2_prev.xx[j], xi_prev.xx[j])
  }
  # compute the ratio
  ratio = min(exp(numerator-denominator), 1)
  
  U = runif(n = 1, min = 0, max = 1)
  if(U < ratio) {
    betaparam.xx[j, ] = beta_prop.xx
    sig2param.xx[j] = sig2_prop.xx
    xiparam.xx[j] = xi_prop.xx
  } 
  else {
    betaparam.xx[j, ] = beta_prev.xx[j, ]
    sig2param.xx[j] = sig2_prev.xx[j]
    xiparam.xx[j] = xi_prev.xx[j]
  }
}
betaparam.xx   #(J,p) : cl.membership(J), predictors(p)
beta_prev.xx #(J,p)
sig2param.xx   #(J)
sig2_prev.xx #(J)
xiparam.xx     #(J)
xi_prev.xx   #(J)



################################################################################
### Step03> Data model development #############################################
########## Using sampled parameters above, ---- discrete / continuous ##########
################################################################################

### [discrete Data model]
#
# covariate w/o NA: ----------> Z,X
#-> dbinom(x, size, prob=piparam)*dnorm(x, mean=X_bar, sd=lambda2param)
# covariate with NA in Z: ---> X only
#-> dnorm(x, mean=X_bar, sd=lambda2param)

# Outcome w/o NA: ------------> S|Z,X
#-> dlogsknorm_log(y, x, beta=betaparam, sig2=sig2param, xi=xiparam)
# Outcome with NA: -----------> Sh|Z,X then "MARGINALIZE" w.r.t "Z"
#-> dlogsknorm_log(y, x, beta=betaparam, sig2=sig2param, xi=xiparam) ???

################################################################################
######### [continuous (Parameter-free) Data model for training Data] ###########
################################################################################
# 1. param-free outcome density : f0S
# 2. param-free covariate density : f0X=f0x1(Z[i])*f0x2(X[i])
# 3. param-free E[outcome] value : E0S


##### For [Param-free covariate model]: 
# covariate model w/o NA: ----------> int Z,X... w.r.t "piparam","X_bar","lambda2param"
# ::: piparam[j] = rbeta( n=1, shape1=g0+sum(Zj), shape2=h0+nj-sum(Zj) )
# ::: lambda2param[j]=rinvgamma( n=1, shape=(c0+n_j)/2, scale=(d0 + sum( (Xj - mean(Xj))^2 ) )/2 )
# ::: X_bar?

#> Analytical form
f0x1 = function(Z) {
  beta(Z+g0.xx, 1-Z+h0.xx)/beta(g0.xx,h0.xx)
}
#> Analytical form
f0x2 = function(X) {
  1/sqrt(pi)*gamma(c0.xx+1)/gamma(c0.xx)*d0.xx^(c0.xx/2)/( (X-mean(X))^2 + d0.xx )^((c0.xx+1)/2)
}
#> JOINT 
f0X.xx = numeric(n.train) #% param-free covariate model f0(X,Z)



##### For [Param-free outcome model]: -------------------> MonteCarlo Integration
# Calculate Outcome and Covariate parameter free data model for each observation  
f0S.xx = numeric(n.train) #% param-free outcome model f0(S|X,Z)
E0S.xx = numeric(n.train) #% param-free E(S|X,Z) = Expected value of S|X,Z

M = 1000 # Number of Monte Carlo samples
sumS.xx = numeric(M)
sumES.xx = numeric(M)

for(h in 1:n.train) {
  
  ##> Analytical solution for f0(X, Z) for ### covariate model ###
  # : necessary parameters are supplied by pre-determined hyperparam values...
  f0X.xx[h] = f0x1(Z[h])*f0x2(Xstar[h])  
  
  ##> Monte Carlo Integration for the ### outcome model ###
  for(j in 1:M) {
    xi_samplej = rt(n=1, df=nu0.xx)                               # bullet on xi
    sig_samplej = rinvgamma(n=1, shape=u0.xx, scale=v0.xx)           # bullet on sig2
    beta_samplej = rmvn(n=1, mu=m0.xx, sigma=SIG_b0.xx*varinf)       # bullet on beta
    
    sumS.xx[j] = dlogsknorm( s=S[h], x=matXstar[h,], beta=beta_samplej, sig2=sig_samplej, xi=xi_samplej )* 
      dmvn(X=beta_samplej, mu=m0.xx, sigma=SIG_b0.xx*varinf)*        # to joint beta
      dinvgamma(x=sig_samplej, shape=u0.xx, scale=v0.xx)*            # to joint sig2
      dt(x=xi_samplej, df=nu0.xx)                                 # to joint xi
    # 20*exp( betaparamclean[j,1] + betaparamclean[j,2]*cleanx + betaparamclean[j,3]*Z[i] - sig2paramclean[j]/2)*( 1-pnorm(-xiparamclean[j]*sqrt(sig2paramclean[j])/sqrt(xiparamclean[j]^2+1)))
    sumES.xx[j] = 18*exp( sum(matXstar[h,]*beta_samplej) - sig_samplej/2 )*(1-pnorm(-xi_samplej*sqrt(sig_samplej)/sqrt(xi_samplej^2+1)))
  }
  E0S.xx[h] = sum(sumES.xx)/M  # E[S] value of f0
  f0S.xx[h] = sum(sumS.xx)/M   # density of f0
  
  print(paste("h=",h))
}

par(mfrow=c(1,1))
plot( x=density(f0X.xx) )            # param-free covariate model
plot( x=density(f0S.xx) )            # param-free outcome model
# plot( x=density(E0S) ) # param-free E[outcome] model (wrong)
plot(x=sort(E0S.xx), y=f0S.xx[order(E0S.xx)], type="l")

summary(S)    # original S
summary(E0S.xx)  # E[S] based on f0(S|x) by train data
# - sort(.) : return a vector in ascending order.
# - order(.) : return the index of each element in a vector in sorted order
# - rank(.) : assign a rank to each element in a vector (smallest = 1)


#..............................................................................#
#..............................................................................#
#..............................................................................#
#................................... now ......................................# 
#..............................................................................#
#..............................................................................#
#..............................................................................#
################################################################################
### Step04> Gibbs Sampler --------- cl_membership and param update ---with J= ?
################################################################################
set.seed(1)

total_iter = 30000
# r_convergence = 1000 <- first run Gibbs sampler to determine r value for convergence
r_convergence = 28000

#loglikelihood.xx = matrix(0, nrow = n.train, ncol = total_iter)
loglikelihood.xx = numeric(total_iter)

# [note] why matrix? each col gives each iteration result

list_piparam.xx = list()        #for Z
list_lambda2param.xx = list()   #for X
list_alphaparam.xx = list()     #for alpha

list_betaparam.xx = list()    #for S
list_sig2param.xx = list()    #for S
list_xiparam.xx = list()      #for S

list_cl.xx = list()

J_prev = J

piparam_prev.xx = piparam.xx
lambda2param_prev.xx = lambda2param.xx
alphaparam_prev.xx = alphaparam.xx

beta_prev.xx = betaparam.xx
sig2_prev.xx = sig2param.xx
xi_prev.xx = xiparam.xx

#### Let's run!!!
for (r in 1:total_iter) {
  ##############################################################################################################
  ##[1] cluster membership computation ----------------------------------------
  clusterout = clusterDPBasic(S, Xstar, Z, cl_membership.xx,
                         piparam.xx, lambda2param.xx,
                         betaparam.xx, sig2param.xx, xiparam.xx, alphaparam.xx,
                         f0X.xx, f0S.xx, u0.xx, v0.xx, m0.xx, SIG_b0.xx, 
                         nu0.xx, g0.xx, h0.xx, c0.xx, d0.xx, gamma0.xx, psi0.xx, varinf) 
  
  cl_membership.xx = clusterout$cl_membership
  J = length(unique(cl_membership.xx))             #% in case, cluster scenario changes, reflecting them
  
  #piparam.xx = clusterout$piparam
  #lambda2param.xx = clusterout$lambda2param
  
  betaparam.xx = clusterout$betaparam
  sig2param.xx = clusterout$sig2param
  xiparam.xx = clusterout$xiparam
  ##############################################################################################################
  
  
  #%%%%%%%%%%%%%%%%%%%%%%%% Heavy Computation ON/OFF %%%%%%%%%%%%%%%%%%%%%%%%%%#
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%# 
  # for (i in 1:n.train){
  #   cluster_si = cl_membership.xx     #% current membership vector
  #   
  #   #%%% Before start, check if there is any cluster that has only a single value. 
  #   # if so...then remove the cluster and need to remove parameters as well and re-number clusters
  #   if( sum(cluster_si==cluster_si[i])==1 ) {
  #     j = cluster_si[i]
  #     cluster_si[cluster_si>j] = cluster_si[cluster_si>j] - 1 
  #     
  #     #% for Z
  #     piparam.xx = piparam.xx[-j]
  #     #% for X
  #     lambda2param.xx = lambda2param.xx[-j]
  #     #% for S
  #     sig2param.xx = sig2param.xx[-j]
  #     betaparam.xx = betaparam.xx[-j,]
  #     xiparam.xx = xiparam.xx[-j]
  #     
  #     J = J-1
  #   }
  #   
  #   # a)remove obv and initialize...?
  #   cluster_si[i] = 0                                  #% replace the membership value (the first obv) with "0"!
  #   nj = as.numeric( table(cluster_si[cluster_si>0]) ) #% number of observations in each cluster without observation i
  #   probs = numeric( length(nj)+1 )                    #% for P(s_i=j) ... it's c(31-1,41,28)?? so 3 + 1 ? total cluster number?
  #   
  #   # b)Iterate through each cluster and Calculate probability of staying the same: P(s_i=j)
  #   x_i = c(1, Xstar[i], Z[i])                                        #% c(1, Z, X)
  #   for(j in 1:length(nj)) {
  #     probs[j] = nj[j]/(n.train-1+alphaparam.xx)*
  #       dlogsknorm(s=S[i], x=x_i, beta=betaparam.xx[j,], sig2=sig2param.xx[j], xi=xiparam.xx[j])* 
  #       dbinom(x=x_i[3], size=1, prob = piparam.xx[j])*                # Covariate Z
  #       dnorm(x=x_i[2], mean=mean(Xstar), sd=sqrt(lambda2param.xx[j])) # Covariate X
  #   }
  #   
  #   # c)After iteration through each cluster, Adding the new probability of "Forming a new cluster": P(s_i=J+1) 
  #   probs[j+1] = alphaparam.xx/(n.train-1+alphaparam.xx)*f0S.xx[i]*f0X.xx[i]   #% so it gives...probs: c(prob, prob, prob, 0) -> c(prob, prob, prob, prob)
  #   
  #   # d)Finally, draw a new cluster for each datapoint from a multinomial distribution  #~~~~~~~~~~~~~~~~Q.David
  #   newclust = which( rmultinom(n=1, size=1, prob=probs)==1 ) #% "which(.)" gives the indices of (logic=True)
  #   
  #   #% reflect whats happening above
  #   cl_membership.xx = cluster_si
  #   
  #   cl_membership.xx[i] = newclust  #% assigning the cl_membership to each datapoint. WOW WOW WOW WOW ! 
  #   #% "Multinomial(.)" is the way the datapt accept/reject the new cluster????
  #   #% to face the truth, probabilities are: "probs/sum(probs)"
  #   #% but.."rmultinom(.)" automatically addresses "probs"
  #   
  #   # e)If new cluster is selected by a certain data point, 
  #   # then add new value (obtained from new cluster) to the existing parameter pool
  #   if( length(cl_membership.xx[cl_membership.xx==cl_membership.xx[i]]) == 1 ) {
  #     #% for Z, append
  #     piparam.xx = c( piparam.xx, 
  #                  rbeta(n=1, shape1=g0.xx, shape2=h0.xx) )
  #     #% for X, append
  #     lambda2param.xx = c( lambda2param.xx, 
  #                       rinvgamma(n=1, shape=c0.xx/2, 
  #                                 scale=d0.xx/2) ) 
  #     #% for S, append
  #     # beta_old_j = rmvn(n=1, mu=m0.xx, sigma=SIG_b0.xx)
  #     # sig2_old_j = rinvgamma(n=1, shape=u0.xx, scale=v0.xx)
  #     # xi_old_j = rt(n=1, df=nu0.xx) # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Q.David
  #     
  #     beta_p = rmvn(n=1, mu=m0.xx, sigma=SIG_b0.xx)
  #     sig2_p = rinvgamma(n=1, shape=u0.xx, scale=v0.xx)
  #     xi_p = rt(n=1, df=nu0.xx)
  #     betaparam.xx = rbind(betaparam.xx, beta_p)
  #     sig2param.xx = c(sig2param.xx, sig2_p)
  #     xiparam.xx = c(xiparam.xx, xi_p) # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Forget MH! 
  #     
  #     # numerator = dlogsknorm_log(S[i], x_i, t(as.matrix(beta_p)), sig2_p, xi_p)
  #     # denominator = dlogsknorm_log(S[i], x_i, t(as.matrix(beta_old_j)), sig2_old_j, xi_old_j)
  #     # ratio = min(exp(numerator-denominator), 1)
  #     # U = runif(n = 1, min = 0, max = 1)
  #     # if(U < ratio) {
  #     #   betaparam.xx = rbind(betaparam.xx, beta_p)
  #     #   sig2param.xx = c(sig2param.xx, sig2_p)
  #     #   xiparam.xx = c(xiparam.xx, xi_p)
  #     # } else {
  #     #   betaparam.xx = rbind(betaparam.xx, beta_old_j)
  #     #   sig2param.xx = c(sig2param.xx, sig2_old_j)
  #     #   xiparam.xx = c(xiparam.xx, xi_old_j)
  #     # }
  #     J = J+1
  #   }
  # }
  # J = length(unique(cl_membership.xx))    #% in case, cluster scenario changes, reflecting them
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%# 
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#  
  
  
  ###[2] Updating Posterior in accordance with new cl_membership ------------------------------------------
  #% for Z - piparam from beta posterior
  piparam.xx = numeric(J)
  for (j in 1:J) {
    Zj=Z[cl_membership.xx==j]
    nj=length(Zj)
    piparam.xx[j]=rbeta( n=1, shape1=g0.xx+sum(Zj), shape2=h0.xx+nj-sum(Zj) ) #posterior
  }
  #% for X - lambda2param update from MVN, InvGa posterior 
  lambda2param.xx = numeric(J) #; lambda2param
  for (j in 1:J) {
    Xj=Xstar[cl_membership.xx==j]
    Zj = Z[cl_membership.xx==j]
    nj=length(Xj)
    
    lambda2param.xx[j]=rinvgamma( n=1, shape=(nj+c0.xx)/2, scale=( sum((Xj - mean(Xstar))^2) + d0.xx )/2 )
  }
  
  #% for alphaparam
  eta.xx = rbeta(n=1, shape1=alpha0.xx+1, shape2=n.train)
  pi_eta.xx = (gamma0.xx+J-1)/(gamma0.xx+J-1 + n.train*(psi0.xx-log(eta.xx)))
  alphaparam.xx = pi_eta.xx*rgamma( n=1, shape=gamma0.xx+J, 
                              rate=psi0.xx-log(eta.xx) ) + (1-pi_eta.xx)*rgamma( n=1, shape=gamma0.xx+J-1, rate=psi0.xx-log(eta.xx) )
  
  #% for S
  beta_old = betaparam.xx
  sig2_old = sig2param.xx
  xi_old = xiparam.xx
  
  
  for(j in 1:J) {
    # Sample proposals from priors
    beta_p = rmvn(n=1, mu=m0.xx, sigma=SIG_b0.xx)
    sig2_p = rinvgamma(n=1, shape=u0.xx, scale=v0.xx)
    xi_p = rt(n=1, df=nu0.xx)
    
    # subsetting by cluster
    Sj = S[cl_membership.xx==j]
    if(length(Sj)==1) {
      matXj = matrix( matXstar[cl_membership.xx==j, ], nrow=1 ) # It is a vector added to the matrix
    } else {
      matXj = matXstar[cl_membership.xx==j, ]                   # It is a regular matrix
    } # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Q.David
    
    #matXj = matX_tru[cl_membership==j, ] 
    
    
    # > MH algorithm
    numerator=0
    denominator=0
    for (i in 1:length(Sj)){
      numerator = numerator + dlogsknorm_log(Sj[i], matXj[i,], 
                                             t(as.matrix(beta_p)), sig2_p, xi_p)
      denominator = denominator + dlogsknorm_log(Sj[i], matXj[i,], 
                                                 t(as.matrix(beta_old[j,])), sig2_old[j], xi_old[j])
      
      ratio = min(exp(numerator-denominator), 1)
      U = runif(n = 1, min = 0, max = 1)
      if(U < ratio) {
        betaparam.xx[j, ] = beta_p
        sig2param.xx[j] = sig2_p
        xiparam.xx[j] = xi_p
      } else {
        betaparam.xx[j, ] = beta_old[j, ]
        sig2param.xx[j] = sig2_old[j]
        xiparam.xx[j] = xi_old[j]
      }
    }
  }
  
  
  
  ###[3] Calculating Loglikelihood ---------------------------------------------
  loglike_r = 0
  
  for(i in 1:n.train) {
    x_i = matXstar[i, ]
    j = cl_membership.xx[i]
    
    loglike_r = loglike_r +                                                     
      dlogsknorm_log(s=S[i], x=x_i, beta=betaparam.xx[j,], sig2=sig2param.xx[j], xi=xiparam.xx[j]) +
      dbinom(x = x_i[3], size = 1, prob = piparam.xx[j], log = TRUE) + 
      dnorm(x = x_i[2], mean = mean(Xstar), sd = sqrt(lambda2param.xx[j]), log = TRUE)
  }
  #loglikelihood.xx[, r] = loglike_r
  loglikelihood.xx[r] = loglike_r
  
  if(r > r_convergence) {
    list_piparam.xx[[r-r_convergence]] = piparam.xx   #for Z
    
    list_lambda2param.xx[[r-r_convergence]] = lambda2param.xx  #for X
    
    list_alphaparam.xx[[r-r_convergence]] = alphaparam.xx     #for alpha
    
    list_betaparam.xx[[r-r_convergence]] = betaparam.xx    #for S
    list_sig2param.xx[[r-r_convergence]] = sig2param.xx    #for S
    list_xiparam.xx[[r-r_convergence]] = xiparam.xx      #for S
    
    list_cl.xx[[r-r_convergence]] = cl_membership.xx
  }
  print(paste("r=",r)) #% iteration progress
}
#plot(apply(loglikelihood.xx, 2, sum), type="l")
#lppd_EX_DP.S.xx <- mean( colSums(loglikelihood.xx) ); lppd_EX_DP.S.xx # -35,615,910

plot( loglikelihood.xx, type="l")
lppd_EX_DP.S.xx <- mean( loglikelihood.xx ); lppd_EX_DP.S.xx # -22067.68




########################################################
#### n=1552,Error: Mat::operator(): index out of bounds
#### n=1552,Error in if (U < ratio) { : missing value where TRUE/FALSE needed

#..............................................................................#
#..............................................................................#
############################### Predictions ####################################
#..............................................................................#
#..............................................................................#

################################################################################
###> Training data -------------------------------------------------------------
################################################################################
#% Ultimate liquid weight vector for each iteration
W_paramFree.vec.xx = matrix(0, nrow = n.train, ncol = total_iter-r_convergence) 
density.train.xx = matrix(0, nrow = n.train, ncol = total_iter-r_convergence)
expval.train.xx = matrix(0, nrow = n.train, ncol = total_iter-r_convergence)

for(r in 1:(total_iter-r_convergence)) {
  piparamclean.xx = list_piparam.xx[[r]]
  lambda2paramclean.xx = list_lambda2param.xx[[r]] #for Xstar
  alphaparamclean.xx = list_alphaparam.xx[[r]]      #for alpha
  
  betaparamclean.xx = list_betaparam.xx[[r]]      #for S
  sig2paramclean.xx = list_sig2param.xx[[r]]      #for S
  xiparamclean.xx = list_xiparam.xx[[r]]       #for S
  
  cl_membership.xx = list_cl.xx[[r]]
  
  J = nrow(betaparamclean.xx)
  f_S.xx = numeric(J)                  # density for S
  
  #% Ultimate solid weight matrix
  W_paramBase.matrix.xx = matrix(0, nrow=n.train, ncol=J)
  
  for(h in 1:n.train) {
    
    fX_j.xx = numeric(J)       # discrete covariate model for each j
    w_j.solid.xx = numeric(J)  # discrete weight for each j  
    E_value.xx = numeric(J)    # to save E[S|X] for each j
    
    for(j in 1:J) {
      # discrete covariate joint
      fX_j.xx[j] = dbinom(x=Z[h], size=1, prob=piparamclean.xx[j])*
        dnorm(x=Xstar[h], mean=mean(Xstar), sd=sqrt(lambda2paramclean.xx[j]))
      
      #% solid weight ***** component
      w_j.solid.xx[j] = length(S[cl_membership.xx==j])/(alphaparamclean.xx + n.train)*fX_j.xx[j]
      
      #% compute pred.density
      f_S.xx[j] = dlogsknorm(s=S[h], x=c(1, Xstar[h], Z[h]), 
                          beta=betaparamclean.xx[j,], sig2=sig2paramclean.xx[j], xi=xiparamclean.xx[j])
      #% compute pred.value
      E_value.xx[j] = 4*exp(sum( c(1, Xstar[h], Z[h])*betaparamclean.xx[j,]) - sig2paramclean.xx[j]/2)*
        (1-pnorm(-xiparamclean.xx[j]*sqrt(sig2paramclean.xx[j])/sqrt(xiparamclean.xx[j]^2+1)))
    }
    
    #% liquid weight ***** component
    WJ1.xx = alphaparamclean.xx/(alphaparamclean.xx + n.train)*f0X.xx[h] #; print(paste("liquid weight=", WJ1))
    
    #% Collecting values and constructing weights for each iteration
    W_paramFree.vec.xx[h,r] = WJ1.xx/(WJ1.xx + sum(w_j.solid.xx))
    W_paramBase.matrix.xx[h, ] = w_j.solid.xx/(WJ1.xx + sum(w_j.solid.xx)) 
    
    #% Collecting values for predictive density (via weighted AVG)
    density.train.xx[h,r] = W_paramFree.vec.xx[h,r]*f0S.xx[h] + sum(W_paramBase.matrix.xx[h,]*f_S.xx)
    #% Collecting values for prediction (via weighted AVG)
    expval.train.xx[h,r] = W_paramFree.vec.xx[h,r]*E0S.xx[h] + sum(W_paramBase.matrix.xx[h,]*E_value.xx)
  }
  print(paste("r=", r))
}

expS.DP.Avg.train.xx <- apply(X=expval.train.xx, MARGIN=1, FUN=mean)

summary(S)
summary(expS.DP.Avg.train.xx)
summary(log(S))
summary(log(expS.DP.Avg.train.xx))
#> summary(log(S))
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#4.277   9.282  10.690  10.713  11.955  16.719 
#> summary(log(expS.DP.train))
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#7.047   9.862  10.683  11.051  11.940  16.907


format( aggregate(log(S), by=list(cl_membership.xx), FUN=summary), scientific=FALSE)
format( aggregate(log(expS.DP.Avg.train.xx), by=list(cl_membership.xx), FUN=summary), scientific=FALSE)





################################################################################
###> Testing data -------------------------------------------------------------
################################################################################ 
###> Parameter-free component
#> Analytical form
f0x1 = function(Z) {
  beta(Z+g0.xx, 1-Z+h0.xx)/beta(g0.xx,h0.xx)
}
#> Analytical form
f0x2 = function(X) {
  1/sqrt(pi)*gamma(c0.xx+1)/gamma(c0.xx)*d0.xx^(c0.xx/2)/( (X-mean(X))^2 + d0.xx )^((c0.xx+1)/2)
}
#> JOINT 
f0X.test.xx = numeric(n.test) #% param-free covariate model f0(X,Z)

##### For [Param-free outcome model]: -------------------> MonteCarlo Integration
# Calculate Outcome and Covariate parameter free data model for each observation  
f0S.test.xx = numeric(n.test) #% param-free outcome model f0(S|X,Z)
E0S.test.xx = numeric(n.test) #% param-free E(S|X,Z) = Expected value of S|X,Z

M = 1000 # Number of Monte Carlo samples
sumS.xx = numeric(M)
sumES.xx = numeric(M)

for(h in 1:n.test) {
  
  ##> Analytical solution for f0(X, Z) for ### covariate model ###
  # : necessary parameters are supplied by pre-determined hyperparam values...
  f0X.test.xx[h] = f0x1(Z.test[h])*f0x2(Xstar.test[h])  
  
  ##> Monte Carlo Integration for the ### outcome model ###
  for(j in 1:M) {
    xi_samplej = rt(n=1, df=nu0.xx)                               # bullet on xi
    sig_samplej = rinvgamma(n=1, shape=u0.xx, scale=v0.xx)           # bullet on sig2
    beta_samplej = rmvn(n=1, mu=m0.xx, sigma=SIG_b0.xx*varinf)       # bullet on beta
    
    sumS.xx[j] = dlogsknorm( s=S.test[h], x=matXstar.test[h,], beta=beta_samplej, sig2=sig_samplej, xi=xi_samplej )* 
      dmvn(X=beta_samplej, mu=m0.xx, sigma=SIG_b0.xx*varinf)*        # to joint beta
      dinvgamma(x=sig_samplej, shape=u0.xx, scale=v0.xx)*            # to joint sig2
      dt(x=xi_samplej, df=nu0.xx)                                 # to joint xi
    # 20*exp( betaparamclean[j,1] + betaparamclean[j,2]*cleanx + betaparamclean[j,3]*Z[i] - sig2paramclean[j]/2)*( 1-pnorm(-xiparamclean[j]*sqrt(sig2paramclean[j])/sqrt(xiparamclean[j]^2+1)))
    sumES.xx[j] = 18*exp( sum(matXstar.test[h,]*beta_samplej) - sig_samplej/2 )*(1-pnorm(-xi_samplej*sqrt(sig_samplej)/sqrt(xi_samplej^2+1)))
  }
  E0S.test.xx[h] = sum(sumES.xx)/M  # E[S] value of f0
  f0S.test.xx[h] = sum(sumS.xx)/M   # density of f0
  
  print(paste("h=",h))
}
summary(S.test)    # original S
summary(E0S.test.xx)  # E[S] based on f0(S|x) by train data





###> NOW....
# cl_membership.test <- full.test.sweden$Risk
# table(cl_membership.test)                    # this is wrong!!!!!!!!!!!!!!!!!

#% Ultimate liquid weight vector for each iteration
W_paramFree.vec.test.xx = matrix(0, nrow = n.test, ncol = total_iter-r_convergence) 

density.test.xx = matrix(0, nrow = n.test, ncol = total_iter-r_convergence)
expval.test.xx = matrix(0, nrow = n.test, ncol = total_iter-r_convergence)


for(r in 1:(total_iter-r_convergence)) {
  piparamclean.xx = list_piparam.xx[[r]]
  lambda2paramclean.xx = list_lambda2param.xx[[r]] #for Xstar
  alphaparamclean.xx = list_alphaparam.xx[[r]]      #for alpha
  
  betaparamclean.xx = list_betaparam.xx[[r]]      #for S
  sig2paramclean.xx = list_sig2param.xx[[r]]      #for S
  xiparamclean.xx = list_xiparam.xx[[r]]       #for S
  
  cl_membership.test.xx = list_cl.xx[[r]]
  

  J = nrow(betaparamclean.xx)
  
  
  f_S.test.xx = numeric(J)                  # density for S
  
  #% Ultimate solid weight matrix
  W_paramBase.matrix.test.xx = matrix(0, nrow=n.test, ncol=J)
  
  for(h in 1:n.test) {
    
    fX_j.xx = numeric(J)       # discrete covariate model for each j
    w_j.solid.xx = numeric(J)  # discrete weight for each j  
    E_value.test.xx = numeric(J)    # to save E[S|X] for each j
    
    for(j in 1:J) {
      # discrete covariate joint
      fX_j.xx[j] = dbinom(x=Z.test[h], size=1, prob=piparamclean.xx[j])*
        dnorm(x=Xstar.test[h], mean=mean(Xstar.test), sd=sqrt(lambda2paramclean.xx[j]))
      
      #% solid weight ***** component #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Q.David. How to determine "cl_membership" for test set?
      w_j.solid.xx[j] = length(S.test[cl_membership.test.xx==j])/(alphaparamclean.xx + n.test)*fX_j.xx[j]
      
      #% compute pred.density
      f_S.test.xx[j] = dlogsknorm(s=S.test[h], x=c(1, Xstar.test[h], Z.test[h]), 
                               beta=betaparamclean.xx[j,], sig2=sig2paramclean.xx[j], xi=xiparamclean.xx[j])
      #% compute pred.value
      E_value.test.xx[j] = 4*exp(sum( c(1, Xstar.test[h], Z.test[h])*betaparamclean.xx[j,]) - sig2paramclean.xx[j]/2)*
        (1-pnorm(-xiparamclean.xx[j]*sqrt(sig2paramclean.xx[j])/sqrt(xiparamclean.xx[j]^2+1)))
    }
    
    #% liquid weight ***** component
    WJ1.xx = alphaparamclean.xx/(alphaparamclean.xx + n.test)*f0X.xx[h] #; print(paste("liquid weight=", WJ1))
    
    #% Collecting values and constructing weights for each iteration
    W_paramFree.vec.test.xx[h,r] = WJ1.xx/(WJ1.xx + sum(w_j.solid.xx))
    W_paramBase.matrix.test.xx[h, ] = w_j.solid.xx/(WJ1.xx + sum(w_j.solid.xx)) 
    
    #% Collecting values for predictive density (via weighted AVG)
    density.test.xx[h,r] = W_paramFree.vec.test.xx[h,r]*f0S.test.xx[h] + sum(W_paramBase.matrix.test.xx[h,]*f_S.test.xx)
    #% Collecting values for prediction (via weighted AVG)
    expval.test.xx[h,r] = W_paramFree.vec.test.xx[h,r]*E0S.test.xx[h] + sum(W_paramBase.matrix.test.xx[h,]*E_value.test.xx)
  }
  print(paste("r=", r))
}


####> [AVG.vec Scenario] 
den.DP.Avg.test.xx <- apply(X=density.test.xx, MARGIN=1, FUN=mean)
expS.DP.Avg.test.xx <- apply(X=expval.test.xx, MARGIN=1, FUN=mean)

summary(S.test)
summary(expS.DP.Avg.test.xx)
summary(log(S.test))
summary(log(expS.DP.Avg.test.xx))

SSPE.DP.LSN.xx <- sum( (log(expS.DP.Avg.test.xx) - log(S.test))^2 ); SSPE.DP.LSN.xx   #276.0281
SAPE.DP.LSN.xx <- sum( abs(log(expS.DP.Avg.test.xx) - log(S.test)) ); SAPE.DP.LSN.xx  #199.5169





#> beta0
list_betaparam0.xx <- numeric(2000)
for (i in 1:2000){
  list_betaparam0.xx[i] <- mean( list_betaparam.xx[[i]][, 1] ) 
}
summary(list_betaparam0.xx)# 7.094

#> beta1
list_betaparam1.xx <- numeric(2000)
for (i in 1:2000){
  list_betaparam1.xx[i] <- mean( list_betaparam.xx[[i]][, 2] ) 
}
summary(list_betaparam1.xx)# 0.8387    

#> beta2
list_betaparam2.xx <- numeric(2000)
for (i in 1:2000){
  list_betaparam2.xx[i] <- mean( list_betaparam.xx[[i]][, 3] ) 
}
summary(list_betaparam2.xx) #-0.5041 

#> sig2
list_sig2test.xx <- numeric(2000)
for (i in 1:2000) {
  list_sig2test.xx[i] <- mean( list_sig2param.xx[[i]] )
}
summary(list_sig2test.xx ) #4.4266   


#> xi
list_xitest.xx <- numeric(2000)
for (i in 1:2000) {
  list_xitest.xx[i] <- mean( list_xiparam.xx[[i]] )
}
summary(list_xitest.xx ) #-0.05656 









# NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW
# NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW
# NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW
# NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW
# NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW
# NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW
# NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW
# NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW
# NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW
# NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW
# NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW
# NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW
# NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW
# NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW
# NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW
# NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW
# NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW
# NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW
# NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW
# NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW
# NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW
# NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW
# NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW
# NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW
# NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW
# NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW
# NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW
# NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW
# NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW
# NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW
# NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW
# NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW
# NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW
# NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW
# NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW

################################################################################
########################## S ~ DP LSN with Gustafson ###########################
################################################################################
################################################################################
### Step00> Define prior
################################################################################

# ::::::: # ------------------------ PRIOR ------------------------- # ::::::: #
# ------------- Outcome ----- # ::: for S ~ LSN( X*betaparam, sig2param, xiparam )

###> PRIOR: "betaparam" ~ MVN(m0, 1/delta*SIG_b0)
#> Gaussian regression for initialize OUTCOME model parameter m0, SIG_b0
fit1.st <- glm( S ~ Xstar + factor(Z) , family=Gamma(link = "log") )
# fit1 <- lm( logS ~ Xstar + factor(Z) )
summary(fit1.st)

m0.st = coef(fit1.st)    # 1x3 initial reg_coeff vector (sort of "means"): Regression result as "mean"
SIG_b0.st = vcov(fit1.st)   # 3x3 initial cov matrix of reg_coeff
p = length(m0.st)
options(scipen = 999)
SIG_b0.st
SIG_b0inv.st = solve(a=SIG_b0.st) # inverse of cov matrix of reg_coeff for later use (for posterior on reg_coeff)!!

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#****** for VAR Inflation factor
delta = 0.01             # equal importance on m0 and beta ????
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#varinf = n.train/50; varinf ###> better sampling stability
varinf = n.train/100 

#varinf.ka = n.train/5; varinf.ka  ###> toooooooooooooo much clusters production 
#varinf.ka = n.train/150; varinf.ka  
varinf.ka = n.train/250; varinf.ka ###> moderate clusters production
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(loglikelihood.st, type="l")




###> PRIOR: "sig2param" ~ IG( u0, v0 )
u0.st = 1.7591                          
v0.st = 7.113 





# ## [CHECK later] ---------------------------------------------------------------
# # How can we design proposal for sig2 ~ InvGa(.) ?
# par(mfrow=c(1,1))
plot( density(unlist(list_sig2param)) ) 
# plot( density(unlist(list_sig2param)), xlim=c(0,50) )                        # sig2 we have...before
curve( dinvgamma(x, shape=u0.st,
                 scale=v0.st), add=TRUE, col="blue", lwd=3)                    # sig2 proposal adjustment

table(unlist(list_sig2param))
sort(table(unlist(list_sig2param)),DESC=TRUE)                             # to see if they are appropriate...
# #-------------------------------------------------------------------------------

###> PRIOR: "xiparam"~t(loc, nu0, sca)
#loc = 0   # ? Location hyperparameter for T distribution
#nu0 = 1/2 # df for T distribution on pigtail parameter
#sca = 5   # ? Scale hyperparameter for T distribution 
nu0.st = n.train-1
#- df refers to the number of independent observations (N-1):
#nu0.st = n.train-100

# ## [CHECK later] ---------------------------------------------------------------
plot( density(unlist(list_xiparam)) ) 
curve( dt(x, df=nu0.st), add=TRUE, col="blue", lwd=3) 
# #-------------------------------------------------------------------------------
 



#] ------ Covariate ------ ::: for Z ~ Bin( 1, [piparam.st] )  
# PRIOR                                         "piparam.st" ~ Beta(g0.st, h0.st)
#------------------------------------------------------------------------------------------------------------------------
g0.st = 0.5 
h0.st = 0.5 

#] ------ Covariate ------ ::: for X ~ N( X_bar, [lambda2param] ) where 
# PRIOR                                            "lambda2param" ~ IG(c0, d0)..... nope!
#------------------------------------------------------------------------------------------------------------------------
# -------------------------::: for X|Z ~ N( K0+K1*Z, [lambda2param.st] ) where 
#                                                  "lambda2param.st" ~ IG(c0.st, d0.st)
fit2.st <- lm(Xstar ~ factor(Z))
summary(fit2.st) # 

#> KAPPA ????
# dirty covariate x is NOT tied to z under the NDB assumption, by this logic, regression of [x|z] and [xstar|z] should return 
# almost the same regression parameter-results...other than parameter-variance...
# It is for the sake of investigation or computation or mathematical convenience ?!?! ...since we need model (x,z)=(x|z)(z)
# It is for the connection between [clean] - [dirty] exposure model ?!?!
# We d not know VAR(Xstar). but by conditioning VAR(Xstar|Z), 
# which is "\hat{lambda}^2", and by comparing VAR(Xstar|Z) vs VAR(X|Z), we can perceive.... tau^2   
kappa_v0 <- coef(fit2.st)
SIG_k0 <- vcov(fit2.st) #*varinf.ka
SIG_inv_k0 <- solve(SIG_k0*varinf.ka); SIG_inv_k0  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#> lambda2param.st
c0.st = 0.5          
d0.st = 0.5    

KK1 = cbind(1, Z) # for the analytical solution for Parameter-Free covariate model development later
# [note] 
# - we need KK1 matrix for all obvs w/o concerning cluster
# - KK1_j  
# - KK2_j         

# # ::::::: # ------------------------ PRIOR ------------------------- # ::::::: #
# # ---- Measurement model ----# ::: for X*|X ~ N( X, [tau2] ) where "tau^2" ~ (1-zeta)*IG(c0.st, d0.st)
# u0 = 1 # for "Unit Information Prior"!!!!!
# 
# #::::::::::::::::::::: NEED TO CHANGE THIS EACH TIME ::::::::::::::::::::::::: #
# # this is tilde: tau2_0
# # seq(from = 0.1, to = 2.0, by = 0.1)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#zeta_scale = 0.93
zeta_scale = 0.95 
#zeta_scale = 0.99
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## for Xstar|X ~ N( X, "tau^2") 
# tau2param <- rep(rinvgamma(n = 1, shape = u0, scale = u0*tau2_0),J)
# tau2param.st = (1-zeta_scale)*lambda2param.st


# ::::::: # ------------------------ PRIOR ------------------------- # ::::::: #
#] ------- precision -----------# :::  for alpha ~ Mix of Ga( gamma0, psi0 )
gamma0.st = 1    #                                          
psi0.st = 1    





################################################################################
### Step01> initialize cluster membership - Kmeans clustering
################################################################################
# J=3
# set.seed(1)
# cluster = kmeans(cbind(S,X_tru,Z), J)
# cl_membership = cluster$cluster 
# table(cl_membership)

cl_membership.st <- full.train.sweden$Risk
table(cl_membership.st)
J=5





################################################################################
### Step02> initialize major parametersssss by clusters using your "POSTERIOR"
#                                                 Based on kmean clustering: J=3
################################################################################

# ---------------- for Exposure model + Measurement Model + alpha --------------


# : Exposure Model (Z only)
##[A] for Z ~ Bin( n=1, "prob_j" ) .. starting with 3 clusters for beta POSTERIOR for Z
set.seed(1)

piparam.st = numeric(J)
for (j in 1:J) {
  Zj=Z[cl_membership.st==j]
  nj=length(Zj)
  piparam.st[j] = rbeta( n=1, shape1=g0.st+sum(Zj), shape2=h0.st+nj-sum(Zj) ) 
}
piparam.st #% pi parameter sampled from posterior


# : Exposure model Xstar|Z
##[B] for Xstar|Z ~ N( "kappa0_j"+"kappa1_j"*Z, "lambda_j" ) .. starting with 3 clusters for Normal/IG POSTERIOR
kappaparam.st = matrix(nrow = J, ncol = length(kappa_v0)); kappaparam.st
lambda2param.st = numeric(J); lambda2param.st
for (j in 1:J) {
  Xstar_j=Xstar[cl_membership.st==j]
  Zj = Z[cl_membership.st==j]
  nj=length(Xstar_j)
  lambda2param.st[j]=rinvgamma( n=1, shape=(nj+c0.st)/2, scale=0.5*( sum((Xstar_j - kappa_v0[1]-kappa_v0[2]*Zj)^2) + d0.st ) )
  KK1j = cbind(1, Zj)
  KK2j = matrix(c(sum(Xstar_j), sum(Xstar_j*Zj)), nrow =2, ncol = 1) # 2x1 vector
  
  SIG_inv_j = solve(SIG_inv_k0+(t(KK1j) %*% KK1j)) # solve( ..[matrix].. ) = [matrix]^(-1)
  kappaparam.st[j,] <- rmvn(n=1, mu = SIG_inv_j %*% (SIG_inv_k0 %*% kappa_v0 + KK2j), sigma = lambda2param.st[j]*SIG_inv_j)      
}
kappaparam.st
lambda2param.st

#% [Check!] empirical kappa, lambda ? using "aggregate(.)": investigation by splitting the data into subset
kappa_empirical=rbind(coef(lm(Xstar[cl_membership.st==1] ~ Z[cl_membership.st==1])),
                      coef(lm(Xstar[cl_membership.st==2] ~ Z[cl_membership.st==2])),
                      coef(lm(Xstar[cl_membership.st==3] ~ Z[cl_membership.st==3]))); kappa_empirical

lambda_empirical = aggregate(x=as.numeric(Xstar), by=list(cl_membership.st), FUN=var, na.rm=T)$x; lambda_empirical




##[C] for alphaparam.st ~ Ga( shape=gamma0.st, rate=psi0.st ) .. starting with 3 clusters for Mixed Gamma POSTERIOR

# alpha0=2 # initialization: typically..1,2...
# eta = rbeta(n=1, shape1=alpha0+1, shape2=n.train); eta
# pi_eta = (g0+J-1)/(g0+J-1+n.train*(h0-log(eta))); pi_eta
# 
# # precision parameters and its done!
# alpha = pi_eta*rgamma( n=1, shape=g0+J, 
#                        rate=h0-log(eta) ) + (1-pi_eta)*rgamma(n=1, shape=g0+J-1, rate=h0-log(eta) ); alpha

alpha0.st=2 # initialization: typically..1,2...
eta.st = rbeta(n=1, shape1=alpha0.st+1, shape2=n.train); eta.st
pi_eta.st = (gamma0.st+J-1)/(gamma0.st+J-1 + n.train*(psi0.st-log(eta.st))); pi_eta.st

# precision parameters and its done!
alphaparam.st = pi_eta.st*rgamma( n=1, shape=gamma0.st+J, 
                               rate=psi0.st-log(eta.st) ) + (1-pi_eta.st)*rgamma(n=1, shape=gamma0.st+J-1, rate=psi0.st-log(eta.st) )
alphaparam.st






################################################################################
#############################--- for outcome ---################################
################################################################################

##[D] ------ Outcome S ------ ::: [betaparam.st], [sig2param.st], [xiparam.st].... try MH-algorithm
# - No conjugacy...so prepare MM sampling

# This is where new parameters go.
betaparam.st = matrix(data=NA, nrow=J, ncol=length(m0.st)) 
sig2param.st = numeric(J) 
xiparam.st = numeric(J)

# Set in place the old parameters..sample from prior...
beta_prev.st = matrix( rep(rmvn(n=1, mu=m0.st, sigma=SIG_b0.st*varinf), J), nrow=J, byrow= TRUE )
sig2_prev.st = rep(rinvgamma(n=1, shape=u0.st, scale=v0.st), J)
xi_prev.st = rep(rt(n=1, df=nu0.st), J) # imagine we already have them...old days..


#### Metropolis Hastings to update betaparam, sig2param, xiparam
for(j in 1:J) {
  # Sample proposals from priors
  beta_p = rmvn(n=1, mu=m0.st, sigma=SIG_b0.st*varinf)
  sig2_p = rinvgamma(n=1, shape=u0.st, scale=v0.st)
  xi_p = rt(n=1, df=nu0.st)
  
  # subsetting by cluster
  Sj = S[cl_membership.st==j]
  matXj = matXstar[cl_membership.st==j, ]
  
  # > MH algorithm
  # prepare components
  numerator=0
  denominator=0
  for (i in 1:length(Sj)){
    numerator = numerator + 
      dlogsknorm_log(Sj[i], matXj[i,], t(as.matrix(beta_p)), sig2_p, xi_p)
    denominator = denominator + 
      dlogsknorm_log(Sj[i], matXj[i,], t(as.matrix(beta_prev.st[j,])), sig2_prev.st[j], xi_prev.st[j])
  }
  # compute the ratio
  ratio = min(exp(numerator-denominator), 1)
  
  U = runif(n = 1, min = 0, max = 1)
  if(U < ratio) {
    betaparam.st[j, ] = beta_p
    sig2param.st[j] = sig2_p
    xiparam.st[j] = xi_p
  } 
  else {
    betaparam.st[j, ] = beta_prev.st[j, ]
    sig2param.st[j] = sig2_prev.st[j]
    xiparam.st[j] = xi_prev.st[j]
  }
}

betaparam.st   #(J,p) : membership(J), predictors(p)
beta_prev.st #(J,p)
sig2param.st   #(J)
sig2_prev.st #(J)
xiparam.st     #(J)
xi_prev.st   #(J)





################################################################################
######### [continuous (Parameter-free) Data model with training Data] ##########
################################################################################

# ------------------------------> MonteCarlo Integration
# Calculate parameter free Outcome and Covariate data model for each observation  

##### For Param-free covariate model: 
# ---> int_{}^{} Xstar,z,piparam,kappaparam,lambda2param w.r.t "piparam","kappaparam","lambda2param"
#> JOINT 
f0X.st = numeric(n.train) #% to save param-free covariate model: JOINT: f0(Xstar, Z)
# pre-discovered...by the analytical soltution...
constX = (d0.st/2)^(c0.st/2)*gamma((c0.st+1)/2)/(sqrt(2*pi)*gamma(c0.st/2)*det(SIG_k0)^0.5)


##### For [Param-free outcome model]: -------------------> MonteCarlo Integration
# Calculate Outcome and Covariate parameter free data model for each observation 
f0S.st = numeric(n.train) #% param-free outcome model f0(S|x)
E0S.st = numeric(n.train) # E(S|x) = Expected value of S|x ~ f0(S|x)


##### Now let's solve the integralsssss
set.seed(1)
M = 1000 # iterations of the Monte Carlo integral
for(h in 1:n.train) {
  # Analytical solution for f0(Xstar, Z) for ### covariate model ###
  # : necessary parameters are supplied by pre-determined hyperparam values...so ~~~~~~~~~~~~~~~~~~~~~~#Q.David
  f0X.st[h] = constX*beta(g0.st+Z[h],h0.st+1-Z[h])/beta(g0.st,h0.st)*det(KK1[h,] %*% t(KK1[h,])+SIG_inv_k0)^(-0.5)/
    (1/2*(d0.st+(Xstar[h]-KK1[h,] %*% kappa_v0)^2/(1+t(KK1[h,]) %*% SIG_k0 %*% KK1[h,])))^((c0.st+1)/2) 
  
  # Monte Carlo integration for S for ### outcome model ###
  # : necessary parameters are supplied by sampling from prior density
  sumS.st = numeric(M)
  sumES.st = numeric(M)
  for(j in 1:M) {
    xi_samplej = rt(n = 1, df = nu0.st)                               # prior on xi
    sig_samplej = rinvgamma(n = 1, shape = u0.st, scale = v0.st)         # prior on sig2
    beta_samplej = rmvn(n = 1, mu = m0.st, sigma = SIG_b0.st*varinf*2)  # prior on beta
    
    sumS.st[j] = 
      dlogsknorm( s=S[h],
                  x=matXstar[h,],
                  beta=beta_samplej,
                  sig2=sig_samplej,
                  xi=xi_samplej )*                                 # outcome with complete
      dmvn(X = beta_samplej, mu = m0.st, sigma = SIG_b0.st*varinf*2)*   # to joint beta
      dinvgamma(x = sig_samplej, shape = u0.st, scale = v0.st)*          # to joint sig2
      dt(x = xi_samplej, df = nu0.st)                                 # to joint xi
    sumES.st[j] = 18*exp(sum(matXstar[h,]*beta_samplej) - sig_samplej/2)*(1-pnorm(-xi_samplej*sqrt(sig_samplej)/sqrt(xi_samplej^2+1)))
    
  }
  E0S.st[h] = sum(sumES.st)/M  # E[S] value of f0
  f0S.st[h] = sum(sumS.st)/M   # density of f0
  
  print(paste("h=",h))
}

plot( x=density(f0S.st) )
plot( x=density(f0X.st) )

summary(E0S.st)  # E[S] based on f0(S|x) by train data
summary(S)    # original S
summary(f0X.st)

plot(x=sort(E0S.st), y=f0S.st[order(E0S.st)], type="l")
# - sort(.) : return a vector in ascending order.
# - order(.) : return the index of each element in a vector in sorted order
# - rank(.) : assign a rank to each element in a vector (smallest = 1)




#..............................................................................#
#..............................................................................#
#..............................................................................#
#................................... now ......................................# 
#..............................................................................#
#..............................................................................#
#..............................................................................#
################################################################################
### Step04> Gibbs Sampler --------- cl_membership and param update ---with J= ?
################################################################################

##### [Gibbs Sampler main game] ::: with correction for X
set.seed(1)
total_iter = 30000
#total_iter = 200
# r_convergence <- first run Gibbs sampler to determine r value for convergence
r_convergence = 28000
#r_convergence = 0

loglikelihood.st = numeric(total_iter)                  # for monitor convergence
#loglikelihood.st = matrix(0, nrow = n.train, ncol = total_iter)
# [note] why matrix? each col gives each iteration result

list_piparam.st = list()   #for Z

list_kappaparam.st = list()     #for Xstar
list_lambda2param.st = list()    #for Xstar

list_tau2param.st = list()       #for measurement model

list_alphaparam.st = list()     #for alpha

list_betaparam.st = list()    #for S
list_sig2param.st = list()    #for S
list_xiparam.st = list()      #for S


list_cl.st = list()

# -------------------- after correction ---------------------- # ~~~~~~~~~~~~~~~~ Q.David
list_betaparam.clean.st = list()    
list_sig2param.clean.st = list()    
list_xiparam.clean.st = list()

# -------------------- after correction ---------------------- # ~~~~~~~~~~~~~~~~ Q.David
#list_piparam.clean = list()        #for Z
list_lambda2param.clean.st = list()
list_kappaparam.clean.st = list()     #for X|Z


#--------------------------------------------------------------------- ALWAYS KEEP and NO-update !!!! --------------------------
### What we have so far................................................................................ for {initialization}
###] Prepare for the current versions of single "cluster-wise" parameters as the "previous" +++++++++++++++++++++++++++++++++

J_prev = J

kappaparam_prev.st = kappaparam.st
lambda2param_prev.st = lambda2param.st

tau2param.st = (1-zeta_scale)*lambda2param.st
tau2_prev.st <- tau2param.st

beta_prev.st = betaparam.st
sig2_prev.st = sig2param.st
xi_prev.st = xiparam.st



###] Miscellaneous ???
max_iter = 10   # let's say.. investigation( sampling rejection ) limit...?
counts_mat <- matrix(0, nrow = J*total_iter, ncol = 3) # count..How many rejected???

n_j.st = table(cl_membership.st)

beta_p_accepted.st <- 0
sig2_p_accepted.st <- 0
xi_p_accepted.st <- 0

beta_p_accepted.clean <- 0
sig2_p_accepted.clean <- 0
xi_p_accepted.clean <- 0


#> Dealing with "Kappa" is so trickyyyyyyyyyyyyyyyyyy?!?!?!?!
#> update problem in beta, sig2, xi ???????????? too small MH ratios ??????????? Q. HOW TO SORT OUT????
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
################################## START #######################################
#### Let's run!!!
for (r in 1:total_iter) {

  #################################################################################################################
  #[1] cluster membership computation ----------------------------------------# ~~~ do we need tau2 ????
  clusterout = clusterDP(S, Xstar, Z, cl_membership.st,
                         piparam.st, kappaparam.st, lambda2param.st,
                         betaparam.st, sig2param.st, xiparam.st, alphaparam.st,
                         f0X.st, f0S.st, u0.st, v0.st, m0.st, SIG_b0.st, nu0.st, g0.st, h0.st,
                         kappa_v0, SIG_k0, SIG_inv_k0, c0.st, d0.st, gamma0.st, psi0.st, varinf)

  cl_membership.st = clusterout$cl_membership
  J = length(unique(cl_membership.st))                  #% in case, cluster scenario changes, reflecting them

  #piparam.st = clusterout$piparam
  #lambda2param.st = clusterout$lambda2param

  kappaparam_prev.st = clusterout$kappaparam
  beta_prev.st = clusterout$betaparam
  sig2_prev.st = clusterout$sig2param
  xi_prev.st = clusterout$xiparam
  #################################################################################################################  
  
  #%%%%%%%%%%%%%%%%%%%%%%%% Heavy Computation ON/OFF %%%%%%%%%%%%%%%%%%%%%%%%%%#
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%# 
  # for (i in 1:n.train){
  #   cluster_si = cl_membership.st     #% current membership vector
  # 
  #   #%%% Before start, check if there is any cluster that has only a single value.
  #   # if so...then remove the cluster and need to remove parameters as well and re-number clusters
  #   if( sum(cluster_si==cluster_si[i])==1 ) {
  #     j = cluster_si[i]
  #     cluster_si[cluster_si>j] = cluster_si[cluster_si>j] - 1
  # 
  #     #% for Z
  #     piparam.st = piparam.st[-j]
  #     #% for X
  #     lambda2param.st = lambda2param.st[-j]
  #     kappaparam.st = kappaparam.st[-j, ]
  #     #% for S
  #     sig2param.st = sig2param.st[-j]
  #     betaparam.st = betaparam.st[-j, ]
  #     xiparam.st = xiparam.st[-j]
  # 
  #     J = J-1
  #   }
  # 
  #   # a)remove obv and initialize...?
  #   cluster_si[i] = 0                                  #% replace the membership value (the first obv) with "0"!
  #   nj = as.numeric( table(cluster_si[cluster_si>0]) ) #% number of observations in each cluster without observation i
  #   probs = numeric( length(nj)+1 )                    #% for P(s_i=j) ... it's c(31-1,41,28)?? so 3 + 1 ? total cluster number?
  # 
  #   # b)Iterate through each cluster and Calculate probability of staying the same: P(s_i=j)
  #   x_i = c(1, Xstar[i], Z[i])                                        #% c(1, X, Z)
  #   for(j in 1:length(nj)) {
  #     probs[j] = nj[j]/(n.train-1+alphaparam.st)*
  #       dlogsknorm(s=S[i], x=x_i, beta=betaparam.st[j,], sig2=sig2param.st[j], xi=xiparam.st[j])*
  #       dbinom(x=x_i[3], size=1, prob = piparam.st[j])*                # Covariate Z
  #       dnorm(x=x_i[2], mean=kappaparam.st[j,1]+kappaparam.st[j,2]*x_i[3], sd=sqrt(lambda2param.st[j])) # Covariate X
  #   }
  # 
  #   # c)After iteration through each cluster, Adding the new probability of "Forming a new cluster": P(s_i=J+1)
  #   probs[j+1] = alphaparam.st/(n.train-1+alphaparam.st)*f0S.st[i]*f0X.st[i]   #% so it gives...probs: c(prob, prob, prob, 0) -> c(prob, prob, prob, prob)
  # 
  #   # d)Finally, draw a new cluster for each datapoint from a multinomial distribution
  #   newclust = which( rmultinom(n=1, size=1, prob=probs)==1 ) #% "which(.)" gives the indices of (logic=True)
  # 
  #   #% reflect whats happening above
  #   cl_membership.st = cluster_si
  # 
  #   cl_membership.st[i] = newclust  #% assigning the cl_membership to each datapoint. WOW WOW WOW WOW !
  #   #% "Multinomial(.)" is the way the datapt accept/reject the new cluster????
  #   #% to face the truth, probabilities are: "probs/sum(probs)"
  #   #% but.."rmultinom(.)" automatically addresses "probs"
  # 
  #   # e)If new cluster is selected by a certain data point,
  #   # then add new value (obtained from new cluster) to the existing parameter pool
  #   if( length(cl_membership.st[cl_membership.st==cl_membership.st[i]]) == 1 ) {
  #     #% for Z, append
  #     piparam.st = c( piparam.st,
  #                  rbeta(n=1, shape1=g0.st, shape2=h0.st) )
  #     
  #     #% for X, append
  #     lambda2param.st = c( lambda2param.st, rinvgamma(n=1, shape=c0.st/2, scale=d0.st/2) )
  # 
  #     kappaparam.st = rbind( kappaparam.st, rmvn(n=1, mu=SIG_inv_j %*% (SIG_inv_k0 %*% kappa_v0 + KK2j),
  #                                                sigma=lambda2param.st[j]*SIG_inv_j) ) #~~~~~~~~~~~~~~~ CHECK
  #     #~~~~~~~~~~~~~~~~~~~~~~~ Q.David: warning: chol(): given matrix is not symmetric: Error: chol(): decomposition failed
  #     #~~~~~~~~~~~~~~~~~~~~~~~ if you forgot the initializing cl_membership, then of course, matrix operation error.......
  # 
  #     #% for S, append
  #     beta_p = rmvn(n=1, mu=m0.st, sigma=SIG_b0.st*varinf) #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ CHECK
  #     sig2_p = rinvgamma(n=1, shape=u0.st, scale=v0.st)
  #     xi_p = rt(n=1, df=nu0.st)
  # 
  #     betaparam.st = rbind(betaparam.st, beta_p)
  #     sig2param.st = c(sig2param.st, sig2_p)
  #     xiparam.st = c(xiparam.st, xi_p) # ~~ Forget MH!
  # 
  #     J = J+1
  #   }
  # }
  # J = length(unique(cl_membership.st))    #% in case, cluster scenario changes, reflecting them
  # 
  # kappaparam_prev.st = kappaparam.st
  # beta_prev.st = betaparam.st  
  # sig2_prev.st = sig2param.st
  # xi_prev.st = xiparam.st
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%# 
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#    

  
  
  ###[2] Updating Posterior ----------------------------------------------------

  
  #> Exception: % for Z - update piparam from beta posterior here coz...it's always clean
  piparam.st = numeric(J)
  for (j in 1:J) {
    Zj=Z[cl_membership.st==j]
    nj=length(Zj)
    piparam.st[j]=rbeta( n=1, shape1=g0.st+sum(Zj), shape2=h0.st+nj-sum(Zj) ) #posterior
  }

  ####% for Xstar|Z - INITIALIZE covariate parameters 
  kappaparam.st = matrix(nrow = J, ncol = length(kappa_v0)) #; kappaparam
  lambda2param.st = numeric(J) #; lambda2param
  
  ####% for Xstar|X - INITIALIZE covariate parameters
  tau2param.st = numeric(J) 
  
  ####% for S - INITIALIZE outcome parameters
  betaparam.st = matrix(data=NA, nrow=J, ncol=length(m0.st)) 
  sig2param.st = numeric(J) 
  xiparam.st = numeric(J)
  
  
  
  #> Exception:#### Extra work!!!! ####
  # If we meet with more clusters generated...then so be it!!!!!!!!!!!!!!!!!!!!!
  if(J_prev < J) {
    # define them
    betaprior = matrix( rep(rmvn(n=1, mu=m0.st, sigma=SIG_b0.st*varinf), J-J_prev), nrow=J-J_prev, byrow=TRUE ) #~~~ CHECK
    sig2prior = rep(rinvgamma(n=1, shape=u0.st, scale=v0.st), J-J_prev)
    xiprior = rep(rt(n=1, df=nu0.st), J-J_prev) # imagine we already have them...like old days..
    
    kappaprior = matrix( rep(rmvn(n=1, mu=kappa_v0, sigma=SIG_k0*varinf.ka), J-J_prev), nrow=J-J_prev, byrow=TRUE )#~~ CHECK
    
    # then append them to the previous pool and match up the dimension 
    beta_prev.st = rbind(beta_prev.st, betaprior)
    sig2_prev.st = c(sig2_prev.st, sig2prior)
    xi_prev.st = c(xi_prev.st, xiprior)

    kappaparam_prev.st = rbind(kappaparam_prev.st, kappaprior)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Q. But why Kappaparam ONLY? what about lambda2param or tau2param ?
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~> 
  }
  # At least ...bring the new cluster into life...!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  
  
  ####% for Xstar - INITIALIZE covariate parameters 
  kappaparamclean.st = matrix(0, nrow = J, ncol = ncol(kappaparam.st))     #for X|Z vs Xstar|Z
  lambda2paramclean.st = numeric(J)                                        #for X|Z vs Xstar|Z
  
  ####% for S - INITIALIZE outcome parameters
  betaparamclean.st =  matrix(0, nrow = J, ncol = ncol(betaparam.st))      #for S
  sig2paramclean.st = numeric(J)                                           #for S
  xiparamclean.st = numeric(J)                                            #for S
  
  ################################ Real GAME ###################################
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  for (j in 1:J) {
    
    Xstar_j=Xstar[cl_membership.st==j]
    Zj = Z[cl_membership.st==j]
    nj=length(Xstar_j)
    
    
    # #> Exception:#### Extra work!!!! ####
    # # If we meet with more clusters generated...then so be it!!!!!!!!!!!!!!!!!!!
    # if(j <= J_prev) {
    #   kappaparam_prev_j = kappaparam_prev.st[j,]
    #   #tau2_prev_j = tau2_prev.st[j]
    # } else {
    #   kappaparam_prev_j = kappa_v0 # point
    #   #tau2_prev_j = tau2_0
    # }
    # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Q. But why Kappaparam ONLY?
    

    ##############################################################
    ##############################################################
    ##### > To meet the conditions for the SYSTEM EQUATION < #####
    difference1 = -1
    difference2 = -1
    difference3 = -1
    
    count_iterations = 0
    
    ###> Don't stop until all the three positive values are collected.
    while(difference1 < 0 | difference2 < 0 | difference3 < 0) {
      
      
      # #> - proposal Sampling [piparamF_j] using proposals OF POSTERIOR {with clusters} 
      # piparam.st[j] = rbeta( n=1, shape1=g0.st+sum(Zj), shape2=h0.st+n_j.st[j]-sum(Zj) )
      # # no need...previous value...? Unaffected by NDB ..............................
      

      #> - proposal Sampling [lambda2paramF_j] using proposals OF POSTERIOR {with clusters}
      lambda2param.st[j]=rinvgamma( n=1, 
                                    shape=(nj+c0.st)/2, 
                                    scale=( sum((Xstar_j - kappaparam_prev.st[j,1]-kappaparam_prev.st[j,2]*Zj)^2) + d0.st )/2 
                                    )
      #####clean_x_j <- rnorm(n=nj, mean=Xstar_j, sd=sqrt(tau2_prev_j))
      #####tau2param.st[j] = rinvgamma(n = 1, shape = (nj+u0)/2, scale = 1/2*(sum((Xstar_j-clean_x_j)^2)+u0*tau2_0))
      tau2param.st[j]=rinvgamma( n=1,
                                    shape=(nj+c0.st)/2,
                                    scale=(1-zeta_scale)*( sum((Xstar_j - kappaparam_prev.st[j,1]-kappaparam_prev.st[j,2]*Zj)^2) + d0.st )/2
                                 )
      #> - tau? tau? tau? tau? tau? tau? tau? tau? tau? tau? tau? tau? tau? tau? tau? tau? tau? tau? tau? tau? tau? 
      #tau2param.st[j] = (1-zeta_scale)*lambda2param.st[j] #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Q....
      
      #> - proposal Sampling [kappaparam]_j ????????????????????????????????????
      KK1j = cbind(1, Zj)
      KK2j = matrix(c(sum(Xstar_j), sum(Xstar_j*Zj)), nrow =2, ncol = 1) # 2x1 vector
      SIG_inv_j = solve(SIG_inv_k0+(t(KK1j) %*% KK1j)) # solve( ..[matrix].. ) = [matrix]^(-1)
      kappaparam.st[j,] <- rmvn( n=1, mu=SIG_inv_j %*% (SIG_inv_k0 %*% kappa_v0 + KK2j), 
                                 sigma=lambda2param.st[j]*SIG_inv_j ) #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ CHECK

      #---------------------------------------------------------------------------------------------------------------
      #---------------------------------------------------------------------------------------------------------------
      #>>> Next,... Sample ****[Main parameters_j]******************************
      
      #>>> subsetting by cluster
      S_j = S[cl_membership.st==j]                                        # outcome
      matX_j.st = matrix(matXstar[cl_membership.st==j, ], nrow=length(S_j))   # covariate
      
      #>>> Sample "new" outcome paramters from PROPOSALS (priors)
      beta_p = rmvn(n=1, mu=m0.st, sigma=SIG_b0.st*varinf) #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ CHECK
      sig2_p = rinvgamma(n=1, shape=u0.st, scale=v0.st)
      xi_p = rt(n=1, df=nu0.st)
      
      
      
      #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      # > MH algorithm.I (all together version)
      #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      # numerator=0
      # denominator=0
      # for (i in 1:length(S_j)){
      #   numerator = numerator +
      #     dlogsknorm_log(S_j[i], matX_j.st[i,], t(as.matrix(beta_p)), sig2_p, xi_p)
      # 
      #   denominator = denominator +
      #     dlogsknorm_log(S_j[i], matX_j.st[i,], t(as.matrix(beta_prev.st[j,])), sig2_prev.st[j], xi_prev.st[j])
      # }
      # # compute the ratio
      # ratio = min(exp(numerator-denominator), 1)
      # #print(ratio)
      # 
      # if(count_iterations < 10) {
      #   U = runif(n = 1, min = 0, max = 1)
      #   if(U < ratio) {
      #     betaparam.st[j, ] = beta_p
      #     sig2param.st[j] = sig2_p
      #     xiparam.st[j] = xi_p
      #   }
      #   else {
      #     betaparam.st[j, ] = beta_prev.st[j, ]
      #     sig2param.st[j] = sig2_prev.st[j]
      #     xiparam.st[j] = xi_prev.st[j]
      #   }
      # } else {
      #   betaparam.st[j, ] = beta_p
      #   sig2param.st[j] = sig2_p
      #   xiparam.st[j] = xi_p
      # }
      #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ##########################################################################
      #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      
      #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ################# > MH algorithm.II (separate version) ###################
      #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      #>>> - Sifting [betaparam]_j #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      numerator=0
      denominator=0
      for (i in 1:length(S_j)){
        numerator = numerator +
          dlogsknorm_log(S_j[i], matX_j.st[i,], t(as.matrix(beta_p)), sig2_prev.st[j], xi_prev.st[j]) + 0
          #dmvn(beta_p, mu=m0.st, sigma=SIG_b0.st*varinf) + 0
          #dmvn(beta_prev.st[j,], mu=m0.st, sigma=SIG_b0.st*varinf)

        denominator = denominator +
          dlogsknorm_log(S_j[i], matX_j.st[i,], t(as.matrix(beta_prev.st[j,])), sig2_prev.st[j], xi_prev.st[j]) + 0
          #dmvn(beta_prev.st[j,], mu=m0.st, sigma=SIG_b0.st*varinf) + 0
          #dmvn(beta_p, mu=m0.st, sigma=SIG_b0.st*varinf)
      }

      #> [CHECK.1] ratio is too small???
      ratio_beta.st = min(exp(numerator-denominator), 1)
      U = runif(n = 1, min = 0, max = 1)
      if(ratio_beta.st > U) {
        betaparam.st[j, ] = beta_p
        beta_p_accepted.st <- beta_p_accepted.st+1
      } else {
        betaparam.st[j,] = beta_prev.st[j,]
      }

      #>>> - Sifting [sig2param]_j #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      numerator=0
      denominator=0
      for (i in 1:length(S_j)){
        numerator = numerator +
          dlogsknorm_log(S_j[i], matX_j.st[i,], t(as.matrix(betaparam.st[j,])), sig2_p, xi_prev.st[j]) + 0
          #log(dinvgamma(sig2_p, shape=u0.st, scale=v0.st)) +
          #log(dinvgamma(sig2_prev.st[j], shape=u0_new.st, scale=v0_new.st))

        denominator = denominator +
          dlogsknorm_log(S_j[i], matX_j.st[i,], t(as.matrix(betaparam.st[j,])), sig2_prev.st[j], xi_prev.st[j]) + 0
          #log(dinvgamma(sig2_prev.st[j], shape=u0.st, scale=v0.st)) +
          #log(dinvgamma(sig2_p, shape=u0_new.st, scale=v0_new.st))
      }

      #> [CHECK.2] ratio is too small???
      ratio_sig2.st = min(exp(numerator-denominator), 1)
      U = runif(n = 1, min = 0, max = 1)
      if(U < ratio_sig2.st) {
        sig2param.st[j] = sig2_p
        sig2_p_accepted.st <- sig2_p_accepted.st+1
      }
      else {
        sig2param.st[j] = sig2_prev.st[j]
      }

      #>>> - Sifting [xiparam]_j #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      numerator=0
      denominator=0
      for (i in 1:length(S_j)){
        numerator = numerator +
          dlogsknorm_log(S_j[i], matX_j.st[i,], t(as.matrix(betaparam.st[j,])), sig2param.st[j], xi_p) + 0
          #dt(xi_p, df=nu0.st, log = TRUE) + 0
          #dt(xi_prev.st[j], df=nu0.st, log = TRUE)

        denominator = denominator +
          dlogsknorm_log(S_j[i], matX_j.st[i,], t(as.matrix(betaparam.st[j,])), sig2param.st[j], xi_prev.st[j]) + 0
          #dt(xi_prev.st[j], df=nu0.st, log = TRUE) + 0
          #dt(xi_p, df=nu0.st, log = TRUE)
      }

      #> [CHECK.3] ratio is too small???
      ratio_xi.st = min(exp(numerator-denominator), 1)
      U = runif(n = 1, min = 0, max = 1)
      if(U < ratio_xi.st) {
        xiparam.st[j] = xi_p
        xi_p_accepted.st <- xi_p_accepted.st+1
      }
      else {
        xiparam.st[j] = xi_prev.st[j]
      }
      #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      
      ### Apply the system Equations (I) for beta_j_1 and sig2_j.... Q. what about xi? kappa? lambda2?
      # ::: pick up the clean samples (beta_j_1, sig2_j) first during Metropolis Hasting
      
      #################################################################################################################
      #################################################################################################################
      #################################################################################################################
      #################################################################################################################
      ###> [2] Calculate differences to ensure VALIDITY of [System Equation]
      #################################################################################################################
      #>>> [main] I. For clean beta_j_1 (the problematic***slope on clean X)
      betaparamclean_j = betaparam.st[j,2]*lambda2param.st[j]/(lambda2param.st[j]-tau2param.st[j])
      
      #>>> [main] II. For clean sig2_j (the problematic***variance param for LSN on clean X) 
      sig2paramclean_j = sig2param.st[j] - betaparamclean_j^2/( 1/tau2param.st[j] + 1/(lambda2param.st[j]-tau2param.st[j]) )
      
      #>>> [main] III. For clean xi_j
      # ?
      
      #>>> [supplementary] How to ensure They are positive ?
      difference1 = lambda2param.st[j]-tau2param.st[j]
      difference2 = sig2paramclean_j
      difference3 = 1/tau2param.st[j]+1/(lambda2param.st[j]-tau2param.st[j]) - xiparam.st[j]^2*betaparamclean_j^2/sig2paramclean_j
      #print(difference)
      count_iterations = count_iterations+1
      #print(c(count_iterations, difference1, difference2, difference3))
    } #%%% END of WHILE %%%
    #print(paste("count=", count_iterations))


    
    # Values that meet the condition are obtained. Now....
    ###> Collect the params and compute the clean parameters using [System Equation].
    #################################################################################################################
    lambda2paramclean.st[j] = lambda2param.st[j] - tau2param.st[j]    
    
    kappaparamclean.st[j,] = kappaparam.st[j,]
    
    betaparamclean.st[j,2] = betaparam.st[j,2]*lambda2param.st[j]/(lambda2param.st[j]-tau2param.st[j]) #for Beta1
    
    betaparamclean.st[j,1] = betaparam.st[j,1] - 
      betaparamclean.st[j,2]*kappaparam.st[j,1]*tau2param.st[j]/(lambda2param.st[j]-tau2param.st[j]) #for Beta0
    
    betaparamclean.st[j,3] = betaparam.st[j,3] - 
      betaparamclean.st[j,2]*kappaparam.st[j,2]*tau2param.st[j]/(lambda2param.st[j]-tau2param.st[j]) #for Beta2
    
    sig2paramclean.st[j] = sig2param.st[j] - 
      betaparamclean.st[j,2]^2/( 1/tau2param.st[j] + 1/lambda2paramclean.st[j] )   #for sig2
    
    xiparamclean.st[j] = sign(xiparam.st[j])*sqrt( xiparam.st[j]^2*(betaparamclean.st[j,2]^2/sig2paramclean.st[j]+1/tau2param.st[j]+1/lambda2paramclean.st[j])/
                                                  (1/tau2param.st[j]+1/lambda2paramclean.st[j] - xiparam.st[j]^2*betaparamclean.st[j,2]^2/sig2paramclean.st[j]) )
    

    #################################################################################################################
  } #%%% END of [j] %%% 
  
  #### For NEXT!!!
  ###> 1. Store the clean data parameters
  ###> 2. Set previous set of parameter to the [[[current]]] for the next iteration
  ## - Manner! main param for next iteration
  
  

  tau2_prev.st = tau2param.st
  
  kappaparam_prev.st = kappaparam.st
  beta_prev.st = betaparam.st
  sig2_prev.st = sig2param.st
  xi_prev.st = xiparam.st
  
  J_prev = J
  

  
  ###############################################################################################################
  ###[3] Calculating Loglikelihood to monitor the chain convergence ---------------------------------------------
  ###############################################################################################################
  loglike_r = 0
  #loglike_r = numeric(n.train)
  
  for(i in 1:n.train) {
    x_i = matXstar[i, ]
    j = cl_membership.st[i]
    
    loglike_r = loglike_r +                                                     
    #loglike_r[i] =   
      # dlogsknorm_log(s = S[i], x = x_i, beta = betaparam.st[j,], sig2 = sig2param.st[j], xi = xiparam.st[j]) +
      # dnorm(x = x_i[2], mean = kappaparam.st[j,1] + kappaparam.st[j,2]*Z[i], sd = sqrt(lambda2param.st[j]), log = TRUE) +
      # dbinom(x = x_i[3], size = 1, prob = piparam.st[j], log = TRUE)
      dlogsknorm_log(s = S[i], x = x_i, beta = betaparamclean.st[j,], sig2 = sig2paramclean.st[j], xi = xiparamclean.st[j]) +
      dnorm(x = x_i[2], mean = kappaparamclean.st[j,1] + kappaparamclean.st[j,2]*Z[i], sd = sqrt(lambda2paramclean.st[j]), log = TRUE) +
      dbinom(x = x_i[3], size = 1, prob = piparam.st[j], log = TRUE)
  }  
  #print(loglike_r)
  loglikelihood.st[r] = loglike_r
  #loglikelihood.st[,r] = loglike_r
  
  if(r > r_convergence) {
    
    list_piparam.st[[r-r_convergence]] = piparam.st      #for Z
    list_kappaparam.st[[r-r_convergence]] = kappaparam.st     #for Xstar
    list_lambda2param.st[[r-r_convergence]] = lambda2param.st  #for Xstar
    list_tau2param.st[[r-r_convergence]] = tau2param.st        # For measurement model
    
    list_alphaparam.st[[r-r_convergence]] = alphaparam.st          #for alpha
    
    list_betaparam.st[[r-r_convergence]] = betaparam.st  #for S
    list_sig2param.st[[r-r_convergence]] = sig2param.st  #for S
    list_xiparam.st[[r-r_convergence]] = xiparam.st      #for S

    list_cl.st[[r-r_convergence]] = cl_membership.st
    
    ################# STORE the CLEAN parameters ##################
    list_kappaparam.clean.st[[r-r_convergence]] = kappaparamclean.st     #for Xstar
    list_lambda2param.clean.st[[r-r_convergence]] = lambda2paramclean.st    #for Xstar
    list_betaparam.clean.st[[r-r_convergence]] = betaparamclean.st    #for S
    list_sig2param.clean.st[[r-r_convergence]] = sig2paramclean.st    #for S
    list_xiparam.clean.st[[r-r_convergence]] = xiparamclean.st        #for S

  }
  print(paste("r=",r))
}
#plot(apply(loglikelihood.st, 2, sum), type="l")
plot(loglikelihood.st, type="l")

beta_p_accepted.st/ length(unlist(list_sig2param.st))
sig2_p_accepted.st / length(unlist(list_sig2param.st))
xi_p_accepted.st/ length(unlist(list_xiparam.st))







#list_cl.st
unlist( lapply( list_cl.st, max ) ) # how much??
list_betaparam.st
list_sig2param.st
list_xiparam.st

lppd_EX_DP.S.st = mean( loglikelihood.st ); lppd_EX_DP.S.st # ;-21604.54 (0.95); -21954 (0.99) 
# crucial 2000
#plot(loglikelihood.st[28000:30000], type="l")
#lppd_EX_DP.S.st = mean( loglikelihood.st[28000:30000] ); lppd_EX_DP.S.st  # -21597.88


#..............................................................................#
#..............................................................................#
############################### Predictions ####################################
#..............................................................................#
#..............................................................................#

################################################################################
###> Training data -------------------------------------------------------------
################################################################################

#% Ultimate liquid weight vector for each iteration
W_paramFree.vec.st = matrix(0, nrow = n.train, ncol = total_iter-r_convergence) 
density.train.st = matrix(0, nrow = n.train, ncol = total_iter-r_convergence)
expval.train.st = matrix(0, nrow = n.train, ncol = total_iter-r_convergence)

for(r in 1:(total_iter-r_convergence)) {
  piparamclean.st = list_piparam.st[[r]]
  lambda2paramclean.st = list_lambda2param.st[[r]] #for Xstar
  alphaparamclean.st = list_alphaparam.st[[r]]      #for alpha
  
  betaparamclean.st = list_betaparam.st[[r]]      #for S
  sig2paramclean.st = list_sig2param.st[[r]]      #for S
  xiparamclean.st = list_xiparam.st[[r]]       #for S
  
  cl_membership.st = list_cl.st[[r]]
  
  J = nrow(betaparamclean.st)
  f_S.st = numeric(J)                  # density for S
  
  #% Ultimate solid weight matrix
  W_paramBase.matrix.st = matrix(0, nrow=n.train, ncol=J)
  
  for(h in 1:n.train) {
    
    fX_j.st = numeric(J)       # discrete covariate model for each j
    w_j.solid.st = numeric(J)  # discrete weight for each j  
    E_value.st = numeric(J)    # to save E[S|X] for each j
    
    for(j in 1:J) {
      # discrete covariate joint
      fX_j.st[j] = dbinom(x=Z[h], size=1, prob=piparamclean.st[j])*
        dnorm(x=Xstar[h], mean=mean(Xstar), sd=sqrt(lambda2paramclean.st[j]))
      
      #% solid weight ***** component
      w_j.solid.st[j] = length(S[cl_membership.st==j])/(alphaparamclean.st + n.train)*fX_j.st[j]
      
      #% compute pred.density
      f_S.st[j] = dlogsknorm(s=S[h], x=c(1, Xstar[h], Z[h]), 
                             beta=betaparamclean.st[j,], sig2=sig2paramclean.st[j], xi=xiparamclean.st[j])
      #% compute pred.value
      E_value.st[j] = 4*exp(sum( c(1, Xstar[h], Z[h])*betaparamclean.st[j,])-sig2paramclean.st[j]/2)*
        (1-pnorm(-xiparamclean.st[j]*sqrt(sig2paramclean.st[j])/sqrt(xiparamclean.st[j]^2+1)))
    }
    
    #% liquid weight ***** component
    WJ1.st = alphaparamclean.st/(alphaparamclean.st + n.train)*f0X.st[h] #; print(paste("liquid weight=", WJ1))
    
    #% Collecting values and constructing weights for each iteration
    W_paramFree.vec.st[h,r] = WJ1.st/(WJ1.st + sum(w_j.solid.st))
    W_paramBase.matrix.st[h, ] = w_j.solid.st/(WJ1.st + sum(w_j.solid.st)) 
    
    #% Collecting values for predictive density (via weighted AVG)
    density.train.st[h,r] = W_paramFree.vec.st[h,r]*f0S.st[h] + sum(W_paramBase.matrix.st[h,]*f_S.st)
    #% Collecting values for prediction (via weighted AVG)
    expval.train.st[h,r] = W_paramFree.vec.st[h,r]*E0S.st[h] + sum(W_paramBase.matrix.st[h,]*E_value.st)
  }
  print(paste("r=", r))
}

expS.DP.Avg.train.st <- apply(X=expval.train.st, MARGIN=1, FUN=mean)

summary(S)
summary(expS.DP.Avg.train.st)
summary(log(S))
summary(log(expS.DP.Avg.train.st))
#> summary(log(S))
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#4.277   9.282  10.690  10.713  11.955  16.719 
#> summary(log(expS.DP.train))
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#6.636   9.956  10.576  11.028  11.811  16.840
format( aggregate(log(S), by=list(cl_membership.st), FUN=summary), scientific=FALSE)
format( aggregate(log(expS.DP.Avg.train.st), by=list(cl_membership.st), FUN=summary), scientific=FALSE)



# # Clusterwise predictive density plots
# par(mfrow = c(1,3))
# 
# for(j in 1:J) {
#   hist(log(S)[cl_membership.st==j], freq=FALSE, breaks=50, col="white", xlim=c(5,15), ylim=c(0,1.6),
#        xlab = paste("Cluster", j))
#   lines(density(log(expS.DP.train.st)[cl_membership.st==j]), col="red", ylim=c(0,1.6), lwd=2)
# }


################################################################################
###> Testing data -------------------------------------------------------------
################################################################################
##### For Param-free covariate model: 
# ---> int_{}^{} Xstar,z,piparam,kappaparam,lambda2param w.r.t "piparam","kappaparam","lambda2param"
#> JOINT 
f0X.test.st = numeric(n.test) #% to save param-free covariate model: JOINT: f0(Xstar, Z)
# pre-discovered...by the analytical soltution...
constX = (d0.st/2)^(c0.st/2)*gamma((c0.st+1)/2)/(sqrt(2*pi)*gamma(c0.st/2)*det(SIG_k0)^0.5)


##### For [Param-free outcome model]: -------------------> MonteCarlo Integration
# Calculate Outcome and Covariate parameter free data model for each observation 
f0S.test.st = numeric(n.test) #% param-free outcome model f0(S|x)
E0S.test.st = numeric(n.test) # E(S|x) = Expected value of S|x ~ f0(S|x)

##### Now let's solve the integralsssss
set.seed(1)
M = 1000 # iterations of the Monte Carlo integral
for(h in 1:n.test) {
  # Analytical solution for f0(Xstar, Z) for ### covariate model ###
  # : necessary parameters are supplied by pre-determined hyperparam values...so ~~~~~~~~~~~~~~~~~~~~~~#Q.David
  f0X.test.st[h] = constX*beta(g0.st+Z.test[h], h0.st+1-Z[h])/beta(g0.st, h0.st)*det(KK1[h,] %*% t(KK1[h,])+SIG_inv_k0)^(-0.5)/
    (1/2*(d0.st+(Xstar.test[h]-KK1[h,] %*% kappa_v0)^2/(1+t(KK1[h,]) %*% SIG_k0 %*% KK1[h,])))^((c0.st+1)/2) 
  
  # Monte Carlo integration for S for ### outcome model ###
  # : necessary parameters are supplied by sampling from prior density
  sumS.st = numeric(M)
  sumES.st = numeric(M)
  for(j in 1:M) {
    xi_samplej = rt(n = 1, df = nu0.st)                               # prior on xi
    sig_samplej = rinvgamma(n = 1, shape = u0.st, scale = v0.st)         # prior on sig2
    beta_samplej = rmvn(n = 1, mu = m0.st, sigma = SIG_b0.st*varinf)  # prior on beta
    
    sumS.st[j] = 
      dlogsknorm( s=S.test[h],
                  x=matXstar.test[h,],
                  beta=beta_samplej,
                  sig2=sig_samplej,
                  xi=xi_samplej )*                                 # outcome with complete
      dmvn(X = beta_samplej, mu = m0.st, sigma = SIG_b0.st*varinf)*   # to joint beta
      dinvgamma(x = sig_samplej, shape = u0.st, scale = v0.st)*          # to joint sig2
      dt(x = xi_samplej, df = nu0.st)                                 # to joint xi
    sumES.st[j] = 18*exp(sum(matXstar.test[h,]*beta_samplej) - sig_samplej/2)*(1-pnorm(-xi_samplej*sqrt(sig_samplej)/sqrt(xi_samplej^2+1)))
    
  }
  E0S.test.st[h] = sum(sumES.st)/M  # E[S] value of f0
  f0S.test.st[h] = sum(sumS.st)/M   # density of f0
  
  print(paste("h=",h))
}
plot( x=density(f0S.test.st) )
plot( x=density(f0X.test.st) )
summary(E0S.test.st)  # E[S] based on f0(S|x) by train data
summary(S)    # original S
summary(f0X.test.st)
plot(x=sort(E0S.test.st), y=f0S.st[order(E0S.test.st)], type="l")



############################### Out-of-sample ##################################
#% Ultimate liquid weight vector for each iteration
W_paramFree.vec.test.st = matrix(0, nrow = n.test, ncol = total_iter-r_convergence) 
density.test.st = matrix(0, nrow = n.test, ncol = total_iter-r_convergence)
expval.test.st = matrix(0, nrow = n.test, ncol = total_iter-r_convergence)

for(r in 1:(total_iter-r_convergence)) {
  piparamclean.st = list_piparam.st[[r]]
  lambda2paramclean.st = list_lambda2param.st[[r]] #for Xstar
  alphaparamclean.st = list_alphaparam.st[[r]]      #for alpha
  
  betaparamclean.st = list_betaparam.st[[r]]      #for S
  sig2paramclean.st = list_sig2param.st[[r]]      #for S
  xiparamclean.st = list_xiparam.st[[r]]       #for S
  
  cl_membership.test.st = list_cl.st[[r]]
  
  J = nrow(betaparamclean.st)
  f_S.test.st = numeric(J)                  # density for S
  
  #% Ultimate solid weight matrix
  W_paramBase.matrix.test.st = matrix(0, nrow=n.test, ncol=J)
  
  for(h in 1:n.test) {
    
    fX_j.test.st = numeric(J)       # discrete covariate model for each j
    w_j.solid.test.st = numeric(J)  # discrete weight for each j  
    E_value.test.st = numeric(J)    # to save E[S|X] for each j
    
    for(j in 1:J) {
      # discrete covariate joint
      fX_j.test.st[j] = dbinom(x=Z.test[h], size=1, prob=piparamclean.st[j])*
        dnorm(x=Xstar.test[h], mean=mean(Xstar.test), sd=sqrt(lambda2paramclean.st[j]))
      
      #% solid weight ***** component
      #w_j.solid.test.st[j] = length(S.test[cl_membership.test.st==j])/(alphaparamclean.st + n.test)*fX_j.test.st[j]
      w_j.solid.test.st[j] = sum(cl_membership.st==j) / (alphaparamclean.st + n.train)*fX_j.test.st[j] #~~~~~~~~QQQ
      #% compute pred.density
      f_S.test.st[j] = dlogsknorm(s=S.test[h], x=c(1, Xstar.test[h], Z.test[h]), 
                             beta=betaparamclean.st[j,], sig2=sig2paramclean.st[j], xi=xiparamclean.st[j])
      #% compute pred.value
      E_value.test.st[j] = 4*exp(sum( c(1, Xstar.test[h], Z.test[h])*betaparamclean.st[j,])-sig2paramclean.st[j]/2)*
        (1-pnorm(-xiparamclean.st[j]*sqrt(sig2paramclean.st[j])/sqrt(xiparamclean.st[j]^2+1)))
    }
    
    #% liquid weight ***** component
    #WJ1.test.st = alphaparamclean.st/(alphaparamclean.st + n.test)*f0X.test.st[h] #; print(paste("liquid weight=", WJ1))
    WJ1.test.st = alphaparamclean.st/(alphaparamclean.st + n.train)*f0X.test.st[h] #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~QQQ
    #% Collecting values and constructing weights for each iteration
    W_paramFree.vec.test.st[h,r] = WJ1.test.st/(WJ1.test.st + sum(w_j.solid.test.st))
    W_paramBase.matrix.test.st[h, ] = w_j.solid.test.st/(WJ1.test.st + sum(w_j.solid.test.st)) 
    
    #% Collecting values for predictive density (via weighted AVG)
    density.test.st[h,r] = W_paramFree.vec.test.st[h,r]*f0S.test.st[h] + sum(W_paramBase.matrix.test.st[h,]*f_S.test.st)
    #% Collecting values for prediction (via weighted AVG)
    expval.test.st[h,r] = W_paramFree.vec.test.st[h,r]*E0S.test.st[h] + sum(W_paramBase.matrix.test.st[h,]*E_value.test.st)
  }
  print(paste("r=", r))
}

expS.DP.Avg.test.st <- apply(X=expval.test.st, MARGIN=1, FUN=mean)

summary(S.test)
summary(expS.DP.Avg.test.st)
summary(log(S.test))
summary(log(expS.DP.Avg.test.st))

SSPE.DP.LSN.st <- sum( (log(expS.DP.Avg.test.st) - log(S.test))^2 ); SSPE.DP.LSN.st #274.619 (0.9); #285.6494 (0.95); 277.4901 (0.99)
SAPE.DP.LSN.st <- sum( abs(log(expS.DP.Avg.test.st) - log(S.test)) ); SAPE.DP.LSN.st #198.8607 (0.9); #204.0293 (0.95); 201.5277 (0.99)




#> beta0
list_betaparam0.st <- numeric(2000)
for (i in 1:2000){
  list_betaparam0.st[i] <- mean( list_betaparam.st[[i]][, 1] ) 
}
summary(list_betaparam0.st)# 7.128(zeta0.85);  7.040(zeta0.9);  7.030(zeta0.95)

list_betaparam0.stst <- numeric(2000)
for (i in 1:2000){
  list_betaparam0.stst[i] <- mean( list_betaparam.clean.st[[i]][, 1] ) 
}
summary(list_betaparam0.stst)# 6.822(zeta0.85); 7.040(zeta0.9);  7.027(zeta0.95)


#> beta1
list_betaparam1.st <- numeric(2000)
for (i in 1:2000){
  list_betaparam1.st[i] <- mean( list_betaparam.st[[i]][, 2] ) 
}
summary(list_betaparam1.st)# 0.8286(zeta0.85); 0.8559(zeta0.9); 0.8512(zeta0.95)

list_betaparam1.stst <- numeric(2000)
for (i in 1:2000){
  list_betaparam1.stst[i] <- mean( list_betaparam.clean.st[[i]][, 2] ) 
}
summary(list_betaparam1.stst)# 0.8910(zeta0.85); 0.8560(zeta0.9); 0.8619(zeta0.95)



#> beta2
list_betaparam2.st <- numeric(2000)
for (i in 1:2000){
  list_betaparam2.st[i] <- mean( list_betaparam.st[[i]][, 3] ) 
}
summary(list_betaparam2.st) #-0.4800(zeta0.85); -0.5227(zeta0.9); -0.5196 (zeta0.95)

list_betaparam2.stst <- numeric(2000)
for (i in 1:2000){
  list_betaparam2.stst[i] <- mean( list_betaparam.clean.st[[i]][, 3] ) 
}
summary(list_betaparam2.stst) #-0.5359(zeta0.85); -0.5228(zeta0.9); -0.5286 (zeta0.95)



#> sig2
list_sig2test.st <- numeric(2000)
for (i in 1:2000) {
  list_sig2test.st[i] <- mean( list_sig2param.st[[i]] )
}
summary(list_sig2test.st ) #1.6458(zeta0.85); 1.578(zeta0.9); 1.494(zeta0.95)

list_sig2test.stst <- numeric(2000)
for (i in 1:2000) {
  list_sig2test.stst[i] <- mean( list_sig2param.clean.st[[i]] )
}
summary(list_sig2test.stst ) #1.4899(zeta0.85); 1.577(zeta0.9); 1.4681(zeta0.95)


#> xi
list_xitest.st <- numeric(2000)
for (i in 1:2000) {
  list_xitest.st[i] <- mean( list_xiparam.st[[i]] )
}
summary(list_xitest.st ) #-0.3238(zeta0.85); -0.8971(zeta0.9); -0.6816 (zeta0.95)

list_xitest.stst <- numeric(2000)
for (i in 1:2000) {
  list_xitest.stst[i] <- mean( list_xiparam.clean.st[[i]] )
}
summary(list_xitest.stst ) #-1.9626(zeta0.85); -0.8982(zeta0.9); -1.3782(zeta0.95)



#> tau2
list_tau2test.st <- numeric(2000)
for (i in 1:2000) {
  list_tau2test.st[i] <- mean( list_tau2param.st[[i]] )
}
summary(list_tau2test.st ) #0.2121(zeta0.85); 0.0007512(zeta0.9); 0.03422(zeta0.95)






























# ----------------------------------- RIVAL ------------------------------------

#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************

########################## Rival = Hierarchy ###################################
################################################################################
##################### S ~ hierarchical LSN with Gustafson ######################
################################################################################

#> Choose one for clustering: Risk, Zone, or k-means
# 1) Clustering Based on Risk
cl_membership.HB <- full.train.sweden$Risk
table(cl_membership.HB)

cl_membership.HB.test <- full.test.sweden$Risk
table(cl_membership.HB.test)






################################################################################
### Step01> Define prior, hyperprior model and hyperparameter values
################################################################################

# ------------- Outcome ----- # ::: for S ~ LSN( X*betaparam, sig2param, xiparam )

# PRIOR: "betaparam" ~ MVN(beta0_j, SIG2_b0_j)
# Hyperprior                beta0_j ~ MVN(m0, 1/delta*SIG_b0_j)
# Hyperprior                          SIG2_b0_j ~ IW( qu0, LAMB )

# ::: for "beta0_j" ~ MVN( m0, 1/delta*SIG_b0_j )
#> Gaussian regression for initialize OUTCOME model parameter m0, SIG_b0
fit1.HB.st <- glm( S ~ Xstar + factor(Z) , family=Gamma(link = "log") )
# fit1 <- lm( logS ~ Xstar + factor(Z) )
summary(fit1.HB.st)

m0.HB.st = coef(fit1.HB.st)    # 1x3 initial reg_coeff vector (sort of "means"): Regression result as "mean"
SIG_b0.HB.st = vcov(fit1.HB.st)   # 3x3 initial cov matrix of reg_coeff
p = length(m0.HB.st)
options(scipen = 999)
SIG_b0.HB.st
SIG_b0inv.HB.st = solve(a=SIG_b0.HB.st) # inverse of cov matrix of reg_coeff for later use (for posterior on reg_coeff)!!

#****** for VAR Inflation factor
delta = 0.01             # equal importance on m0 and beta ????
varinf = n.train/100 # rule of thumb divide by 5


#> SIG_b0
#- it requires -------
qu0.HB.st = p+2
# LAMB ?????????????????????????????????????????????
LAMB.HB.st = SIG_b0.HB.st


###> PRIOR: "sig2param" ~ IG( u0, v0 )
# Hyperprior               u0 ~ Fink( rho_u1, rho_u2 )
# Hyperprior                   v0 ~ Ga( rho_bv1, rho_v2 )

# Hyperparameter u0
rho_u1.HB.st = 1/8
rho_u2.HB.st = 3/2 

# Hyperparameter v0
rho_v1.HB.st = 8 
rho_v2.HB.st = 1 

# Based on the result? mean(unlist(list_u0param)), mean(unlist(list_v0param))
# E[invGa] = scale/(shape-1) = 50
u0_new.HB.st = 1.7591 #3 #2 #1.1 # from 3/7  # from # summary(unlist(list_u0param)) # for u0_new                            
v0_new.HB.st = 7.113  #3 #5 #8.5 # from 3/7  # from # summary(unlist(list_v0param)) # for v0_new

## [CHECK later] ---------------------------------------------------------------
# How can we design proposal for sig2 ~ InvGa(.) ?
par(mfrow=c(1,1))
plot( density(unlist(list_sig2param.HB.st)) )                                   # sig2 we have...before
curve( dinvgamma(x, shape=mean(unlist(list_u0param.HB.st)), 
                 scale=mean(unlist(list_v0param.HB.st))), add=TRUE, col="red") # sig2 we obtain....after
curve( dinvgamma(x, shape=u0_new.HB.st, 
                 scale=v0_new.HB.st), add=TRUE, col="blue")                    # sig2 proposal adjustment

table(unlist(list_sig2param.HB.st))
sort(table(unlist(list_sig2param.HB.st)),DESC=TRUE)                             # to see if they are appropriate...
#-------------------------------------------------------------------------------

###> PRIOR: "xiparam"~t(loc, nu0, sca)
#loc = 0   # ? Location hyperparameter for T distribution
#nu0 = 1/2 # df for T distribution on pigtail parameter
#sca = 5   # ? Scale hyperparameter for T distribution 
nu0.HB.st = n.train-1
#- df refers to the number of independent observations (N-1):


# ------------------------------------------------------------------------------ #











#] ------ Covariate ------ ::: for Z ~ Bin( 1, [piparam.HB.st] )  
# PRIOR                                         "piparam.HB.st" ~ Beta(g0.HB.st, h0.HB.st)
#------------------------------------------------------------------------------------------------------------------------
g0.HB.st = 0.5 
h0.HB.st = 0.5 

#] ------ Covariate ------ ::: for X ~ N( X_bar, [lambda2param] ) where 
# PRIOR                                            "lambda2param" ~ IG(c0, d0)..... nope!
#------------------------------------------------------------------------------------------------------------------------
# -------------------------::: for X|Z ~ N( K0+K1*Z, [lambda2param.HB.st] ) where 
#                                                  "lambda2param.HB.st" ~ IG(c0.HB.st, d0.HB.st)
fit2.HB.st <- lm(Xstar ~ factor(Z))
summary(fit2.HB.st) # 

#> KAPPA ????
# dirty covariate x is NOT tied to z under the NDB assumption, by this logic, regression of [x|z] and [xstar|z] should return 
# almost the same regression parameter-results...other than parameter-variance...
# It is for the sake of investigation or computation or mathematical convenience ?!?! ...since we need model (x,z)=(x|z)(z)
# It is for the connection between [clean] - [dirty] exposure model ?!?!
# We d not know VAR(Xstar). but by conditioning VAR(Xstar|Z), 
# which is "\hat{lambda}^2", and by comparing VAR(Xstar|Z) vs VAR(X|Z), we can perceive.... tau^2   

kappa_v0.HB <- coef(fit2.HB.st)
SIG_k0.HB <- vcov(fit2.HB.st)*varinf
SIG_inv_k0.HB <- solve(SIG_k0.HB)


#> lambda2param.HB.st
c0.HB.st = 0.5          
d0.HB.st = 0.5    

KK1.HB = cbind(1, Z) # for the analytical solution for Parameter-Free covariate model development later
# [note] 
# - we need KK1 matrix for all obvs w/o concerning cluster
# - KK1_j  
# - KK2_j         

# # ::::::: # ------------------------ PRIOR ------------------------- # ::::::: #
# # ---- Measurement model ----# ::: for X*|X ~ N( X, [tau2] ) where "tau^2" ~ (1-zeta)*IG(c0.HB.st, d0.HB.st)
# u0 = 1 # for "Unit Information Prior"!!!!!
# 
# #::::::::::::::::::::: NEED TO CHANGE THIS EACH TIME ::::::::::::::::::::::::: #
# # this is tilde: tau2_0
# # seq(from = 0.1, to = 2.0, by = 0.1)

# tau2 = (1-zeta_scale)*lambda2dirty
# zeta_scale = c(0.01,0.05,0.1,0.2,0.5)
# zeta_scale = c(0.60,0.70,0.80,0.90,0.95)

#zeta_scale = 0.95



################################################################################
### Step02> initialize major parametersssss by clusters using your "POSTERIOR"
#                                                 Based on natural clustering
################################################################################
J.HB = length(table(cl_membership.HB))
set.seed(1)

# ---------------- for Exposure model + Measurement Model --------------
# : Exposure Model (Z only)
##[A] for Z ~ Bin( n=1, "prob_j" ) .. starting with 3 clusters for beta POSTERIOR for Z

piparam.HB.st = numeric(J.HB)
for (j in 1:J.HB) {
  Zj=Z[cl_membership.HB==j]
  nj=length(Zj)
  piparam.HB.st[j] = rbeta( n=1, shape1=g0.HB.st+sum(Zj), shape2=h0.HB.st+nj-sum(Zj) ) 
}

piparam.HB.st #% pi parameter sampled from posterior

#% [Check!] empirical pi ? using "aggregate(.)": investigation by splitting the data into subset
pi_empirical = aggregate(x=Z, by=list(cl_membership.HB), FUN=mean, na.rm= T)$x; pi_empirical
piparam.HB.st - pi_empirical 




# : Exposure model Xstar|Z
##[B] for Xstar|Z ~ N( "kappa0_j"+"kappa1_j"*Z, "lambda_j" ) .. starting with 3 clusters for Normal/IG POSTERIOR
kappaparam.HB.st = matrix(nrow = J.HB, ncol = length(kappa_v0.HB)); kappaparam.HB.st
lambda2param.HB.st = numeric(J.HB); lambda2param.HB.st
for (j in 1:J.HB) {
  Xstar_j=Xstar[cl_membership.HB==j]
  Zj = Z[cl_membership.HB==j]
  nj.HB.st=length(Xstar_j)
  lambda2param.HB.st[j]=rinvgamma( n=1, shape=(nj.HB.st+c0.HB.st)/2, scale=( sum((Xstar_j - kappa_v0.HB[1]-kappa_v0.HB[2]*Zj)^2) + d0.HB.st )/2 )
  KK1.HBj = cbind(1, Zj)
  KK2.HBj = matrix(c(sum(Xstar_j), sum(Xstar_j*Zj)), nrow =2, ncol = 1) # 2x1 vector
  
  SIG_inv_j = solve(SIG_inv_k0.HB+(t(KK1.HBj) %*% KK1.HBj)) # solve( ..[matrix].. ) = [matrix]^(-1)
  kappaparam.HB.st[j,] <- rmvn(n=1, mu = SIG_inv_j %*% (SIG_inv_k0.HB %*% kappa_v0.HB + KK2.HBj), sigma = lambda2param.HB.st[j]*SIG_inv_j)      
}


kappaparam.HB.st
lambda2param.HB.st

#% [Check!] empirical kappa, lambda ? using "aggregate(.)": investigation by splitting the data into subset
kappa_empirical=rbind(coef(lm(Xstar[cl_membership.HB==1] ~ Z[cl_membership.HB==1])),
                      coef(lm(Xstar[cl_membership.HB==2] ~ Z[cl_membership.HB==2])),
                      coef(lm(Xstar[cl_membership.HB==3] ~ Z[cl_membership.HB==3]))); kappa_empirical

lambda_empirical.HB.st = aggregate(x=as.numeric(Xstar), by=list(cl_membership.HB), FUN=var, na.rm=T)$x; lambda_empirical.HB.st




# : Measurement Model initial parameter tau^2
##[D] for Xstar|X ~ N( X, "tau^2") 
# tau2param.HB <- rep(rinvgamma(n = 1, shape = u0, scale = u0*tau2_0),J.HB)
tau2param.HB = (1-zeta_scale)*lambda2param.HB.st


# ---------------------------- for Outcome -------------------------------------

###> Outcome parameters:  betaparam, sig2param, xiparam

# No conjugacy...so prepare MM sampling

# This is where new parameters go.
betaparam.HB.st = matrix(data=NA, nrow=J.HB, ncol=p) 
sig2param.HB.st = numeric(J.HB)
sig2paramfull.HB.st = 0
xiparam.HB.st = numeric(J.HB)



### (A) beta0, SIG_b0 For [betaparam.HB.st]
# 1) beta0 ~ MVN(m0, 1/delta*SIG_b0)  
# 2) SIG_b0 ~ ~ IW( qu0, LAMB )
# initial hyperparameter values from hyperprior +++++++++++++++++++++++(no need to care cluster)
# Set in place the old parameters..sample from prior...
# SIG_b0_prev_J = riwish(v = qu0, S = SIG_b0*varinf) # take one sample for SIG_b0, this will be initial parameter for all J clusters
# SIG_b0_prev = array(rep(SIG_b0_prev_J,J),dim=c(p,p,J))
# beta0_prev_J = rmvn(n = 1, mu = m0, sigma = 1/delta*SIG_b0_prev_J)
# beta0_prev = matrix(rep(beta0_prev_J,J),nrow=J,ncol=p,byrow = TRUE)
set.seed(1)
SIG_b0_prev.HB.st = riwish(v=qu0.HB.st, S=LAMB.HB.st) 
beta0_prev.HB.st = rmvn(n=1, mu=m0.HB.st, sigma = 1/delta*SIG_b0_prev.HB.st)



### (B) u0.HB.st, v0.HB.st For [sig2param.HB.st]
# 1) u0.HB.st ~ Fink(rho_u1.HB.st, rho_u2.HB.st)
# 2) v0.HB.st ~ Ga(rho_v1.HB.st, rho_v2.HB.st)
# Fink function
fink = function(u0, A, B) {
  A^(u0/2-1)/gamma(u0/2)^B
}

u0_proposed_shape.HB.st = 3 # 3
u0_proposed_rate.HB.st = 2 # 1

###> initial hyperparameter values from hyperprior +++++++++++++++++++++++(no need to care cluster)
u0_prev.HB.st = u0_proposed_shape.HB.st/u0_proposed_rate.HB.st     # instead of sampling from Fink function???
v0_prev.HB.st = rgamma(n=1, shape=rho_v1.HB.st, rate=rho_v2.HB.st)

#jpeg(file="plot.test.jpeg", width=1000, height=800)
#-------------------------------------------------------------------------------
# See? Fink function is close to Ga(3.75,2) ?
curve(fink(x,rho_u1.HB.st,rho_u2.HB.st)/integrate(function(x) fink(x, rho_u1.HB.st, rho_u2.HB.st),0,Inf)$value, 
      from = 0.01, to = 10, ylim = c(0,0.5))
curve(dgamma(x, shape = u0_proposed_shape.HB.st, rate = u0_proposed_rate.HB.st), add=TRUE, col="blue")
#-------------------------------------------------------------------------------


### (C) nu0 For [xiparam] ~ t(loc, nu0, sca)
# nu0 ~ Exponetial(.) ??? or constant? N-1 : dt(x, df, ncp, log = FALSE) : the central t distribution
curve( dt(x, df=nu0.HB.st) )



#-------------------------------------------------------------------------------
##### Now....ready for Single SAMPLE FROM the PRIOR
beta_prev_full.HB.st = rmvn(n=1, mu=beta0_prev.HB.st, sigma=SIG_b0_prev.HB.st)
beta_prev.HB.st = matrix( rep(beta_prev_full.HB.st, J.HB), nrow=J.HB, byrow= TRUE )
sig2_prev_full.HB.st = rinvgamma(n=1, shape=u0_prev.HB.st, scale=v0_prev.HB.st)
sig2_prev.HB.st = rep(sig2_prev_full.HB.st, J.HB)
xi_prev_full.HB.st = rt(n=1, df=nu0.HB.st)
xi_prev.HB.st = rep(xi_prev_full.HB.st, J.HB) # imagine we already have them...old days..


#####> Second, Prepare for Metropolis Hastings for [betaparam], [sig2param], [xiparam] ~ LSN(.)*MVN(_)*InvGa(_)*T(_)
###########################################################################
#> It's time to use......... 
# - "beta0_prev", "SIG_b0_prev", "u0_prev", "v0_prev", "nu0_prev"             
#    for MVN(_), InvGa(_) :based on population
# - "beta_prev", "sig2_prev", "xi_prev", "beta_prev_full", "sig2_prev_full", "xi_prev_full"  
#    for LSN(.)         :based on cluster-wise

#> update hyperparameter first, then ..... update a single sampling: [main parameter]
# 1) beta0param, SIG_b0param for: [betaparam]
# 2) u0param, v0param for: [sig2param]
# 3) nu0param for: [xiparam] 
#---------------------------------------------------------------------------------------------------------------
######################## from posterior of hyperprior ##########################

##> Sample SIG_b0param from posterior inverse wishart
SIG_b0param.HB.st = riwish(v = qu0.HB.st+2,
                           S = t(beta0_prev.HB.st-beta_prev_full.HB.st) %*% (beta0_prev.HB.st-beta_prev_full.HB.st)+
                             delta*t(beta0_prev.HB.st - m0.HB.st) %*% (beta0_prev.HB.st - m0.HB.st) + LAMB.HB.st)

##> Sample beta0param from multivariate normal
beta0param.HB.st = rmvn(n = 1, mu = delta/(delta+1)*m0.HB.st + 1/(delta+1)*beta_prev_full.HB.st,
                        sigma = 1/(delta+1)*SIG_b0param.HB.st)

##> Sample u0 from Metropolis Hastings using Gamma(3,1) as proposal
u0_proposal.HB.st = rgamma(n = 1, shape = u0_proposed_shape.HB.st, rate = u0_proposed_rate.HB.st)
ratio_u0.HB.st = min(fink(u0_proposal.HB.st, A = sig2_prev_full.HB.st*rho_u1.HB.st/(v0_prev.HB.st/2), B = rho_u2.HB.st+1)/
                       fink(u0_prev.HB.st, A = sig2_prev_full.HB.st*rho_u1.HB.st/(v0_prev.HB.st/2), B = rho_u2.HB.st+1)*
                       dgamma(x = u0_prev.HB.st, shape = u0_proposed_shape.HB.st, rate = u0_proposed_rate.HB.st)/
                       dgamma(x = u0_proposal.HB.st, shape = u0_proposed_shape.HB.st, rate = u0_proposed_rate.HB.st),1)
if(runif(n = 1, min = 0, max = 1) < ratio_u0.HB.st) {
  u0param.HB.st = u0_proposal.HB.st
} else {
  u0param.HB.st = u0_prev.HB.st
}

## Sample v0 from gamma distribution
v0param.HB.st = rgamma(n=1, shape=rho_v1.HB.st + u0param.HB.st/2, rate=rho_v2.HB.st + 0.5/sig2_prev_full.HB.st)

SIG_b0param.HB.st; beta0param.HB.st; u0param.HB.st; v0param.HB.st  
# for sharing...#### Now we have single samples to start MH engine.


#####> From here!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
############################################################################### "Cluster-wise"!!!!!!!
##] Run your MH to obtain single posterior samples of serious main parameters: [betaparamS], [sig2param] 

for(j in 1:J.HB) {
  
  # Sample proposals from priors
  beta_prop.HB.st = rmvn(n=1, mu=beta0param.HB.st, sigma=SIG_b0param.HB.st) #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~(CHECK)
  sig2_prop.HB.st = rinvgamma(n=1, shape=u0param.HB.st, scale=v0param.HB.st)
  xi_prop.HB.st = rt(n=1, df=nu0.HB.st)
  #sig2_prop = rinvgamma(n=1, shape=u0param*3, scale=v0param*2)
  #xi_prop = rnorm(n = 1, mean = 0, sd = nu0.HB.st/10)
  
  # subsetting by cluster
  Sj = S[cl_membership.HB==j]
  #logSj = logS[cl_membership.HB==j]
  matXj = matX_tru[cl_membership.HB==j, ]
  
  
  
  ##> A. MH algorithm for: [betaparam]
  #-Note: "beta_prop" is from PRIOR proposal (new) while "beta_prev[j,]" is from the single POSTERIOR sample (old). 
  numerator=0
  denominator=0
  for (i in 1:length(Sj)){
    #for (i in 1:length(logSj)){
    numerator = numerator + 
      #dsknorm_log(logSj[i], matXj[i,], t(as.matrix(beta_p)), sig2_prev[j], xi_prev[j])
      #dsknorm_log(Sj[i], matXj[i,], t(as.matrix(beta_p)), sig2_p, xi_p)
      dlogsknorm_log(Sj[i], matXj[i,], t(as.matrix(beta_prop.HB.st)), sig2_prev.HB.st[j], xi_prev.HB.st[j])
    
    denominator = denominator + 
      #dsknorm_log(logSj[i], matXj[i,], t(as.matrix(beta_prev[j,])), sig2_prev[j], xi_prev[j])
      #dsknorm_log(Sj[i], matXj[i,], t(as.matrix(beta_prev[j,])), sig2_prev[j], xi_prev[j])
      dlogsknorm_log(Sj[i], matXj[i,], t(as.matrix(beta_prev.HB.st[j,])), sig2_prev.HB.st[j], xi_prev.HB.st[j])
  }
  # compute the ratio
  ratio_beta.HB.st = min(exp(numerator-denominator), 1)
  
  U = runif(n = 1, min = 0, max = 1)
  if(U < ratio_beta.HB.st) {
    betaparam.HB.st[j, ] = beta_prop.HB.st
  } 
  else {
    betaparam.HB.st[j, ] = beta_prev.HB.st[j, ]
  }
  
  
  ##> B. MH algorithm for: [sig2param]
  # prepare components
  numerator=0
  denominator=0
  for (i in 1:length(Sj)){
    #for (i in 1:length(logSj)){
    numerator = numerator + 
      dlogsknorm_log(Sj[i], matXj[i,], t(as.matrix(betaparam.HB.st[j,])), sig2_prop.HB.st, xi_prev.HB.st[j]) + 
      #dsknorm_log(logSj[i], matXj[i,], t(as.matrix(betaparam[j,])), sig2_p, xi_prev[j]) + 
      log(dinvgamma(sig2_prop.HB.st, shape = u0param.HB.st, scale = v0param.HB.st)) + 
      log(dinvgamma(sig2_prev.HB.st[j], shape = u0param.HB.st, scale = v0param.HB.st))
    #log(dinvgamma(sig2_prev[j], shape = 3*u0param, scale = 2*v0param)) #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~(CHECK)
    
    denominator = denominator + 
      dlogsknorm_log(Sj[i], matXj[i,], t(as.matrix(betaparam.HB.st[j,])), sig2_prev.HB.st[j], xi_prev.HB.st[j])+ 
      #dsknorm_log(logSj[i], matXj[i,], t(as.matrix(betaparam[j,])), sig2_prev[j], xi_prev[j])+ 
      log(dinvgamma(sig2_prev.HB.st[j], shape = u0param.HB.st, scale = v0param.HB.st)) + 
      log(dinvgamma(sig2_prop.HB.st, shape = u0param.HB.st, scale = v0param.HB.st))
    #log(dinvgamma(sig2_prop, shape = 3*u0param, scale = 2*v0param))   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~(CHECK)
  }
  # compute the ratio
  ratio_sig2.HB.st = min(exp(numerator-denominator), 1)
  
  U = runif(n = 1, min = 0, max = 1)
  if(U < ratio_sig2.HB.st) {
    sig2param.HB.st[j] = sig2_prop.HB.st
  } 
  else {
    sig2param.HB.st[j] = sig2_prev.HB.st[j]
  }
  
  
  
  ##> C. MH algorithm for: [xiparam]
  # prepare components
  numerator=0
  denominator=0
  for (i in 1:length(Sj)){
    #for (i in 1:length(logSj)){
    numerator = numerator + 
      dlogsknorm_log(Sj[i], matXj[i,], t(as.matrix(betaparam.HB.st[j,])), sig2param.HB.st[j], xi_prop.HB.st) +
      #dsknorm_log(logSj[i], matXj[i,], t(as.matrix(betaparam[j,])), sig2param[j], xi_p) +
      dt(xi_prop.HB.st, df=nu0.HB.st, log = TRUE) + dt(xi_prev.HB.st[j], df=nu0.HB.st, log = TRUE)
    #dt(xi_prop, df=nu0, log = TRUE) + dnorm(xi_prev[j], mean = 0, sd = nu0, log = TRUE) 
    #dt(xi_prop, df = nu0, log = TRUE) + dnorm(xi_prev[j], mean = 0, sd = nu0/10, log = TRUE) #~~~~~~(CHECK)
    denominator = denominator + 
      dlogsknorm_log(Sj[i], matXj[i,], t(as.matrix(betaparam.HB.st[j,])), sig2param.HB.st[j], xi_prev.HB.st[j]) +
      #dsknorm_log(logSj[i], matXj[i,], t(as.matrix(betaparam[j,])), sig2param[j], xi_prev[j]) +
      dt(xi_prev.HB.st[j], df=nu0.HB.st, log = TRUE) + dt(xi_prop.HB.st, df=nu0.HB.st, log = TRUE)
    #dt(xi_prev[j], df=nu0, log = TRUE) + dnorm(xi_prop, mean = 0, sd = nu0, log = TRUE)
    #dt(xi_prev[j], df=nu0, log = TRUE) + dnorm(xi_prop, mean = 0, sd = nu0/10, log = TRUE)  #~~~~~~(CHECK)
  }
  # compute the ratio
  ratio_xi.HB.st = min(exp(numerator-denominator), 1)
  
  U = runif(n = 1, min = 0, max = 1)
  if(U < ratio_xi.HB.st) {
    xiparam.HB.st[j] = xi_prop.HB.st
  } 
  else {
    xiparam.HB.st[j] = xi_prev.HB.st[j]
  }
}
betaparam.HB.st; sig2param.HB.st; xiparam.HB.st
beta_prev.HB.st; sig2_prev.HB.st; xi_prev.HB.st










#..............................................................................#
#..............................................................................#
#..............................................................................#
#................................... now ......................................# 
#..............................................................................#
#..............................................................................#
#..............................................................................#
################################################################################
### Step04> Gibbs Sampler --------- cl_membership and param update ---with J.HB= ?
################################################################################

set.seed(1)
total_iter = 30000
# r_convergence = 1000 <- first run Gibbs sampler to determine r value for convergence
r_convergence = 28000

# loglikelihood = numeric(total_iter)                  # for monitor convergence
loglikelihood.HB.st = matrix(0, nrow = n.train, ncol = total_iter)
# [note] why matrix? each col gives each iteration result



#########################] Prepare a "pool"(list) to store the ultimate hyperparameters, parameters [##########################

#> w/o cluster ++++++++++++++++++++++++++++++++ (hyperparameter)
list_beta0param.HB.st = list() # For beta0 hyperparameter
list_SIG_b0param.HB.st = list() # For SIG_b0 hyperparameter
list_u0param.HB.st = list() # for u0 hyperparameter
list_v0param.HB.st = list() # for v0 hyperparameter

#> w/o cluster +++++++++++++++++++++++++++++++++++++++ (outcome)
list_betaparamfull.HB.st = list() # for S - full data
list_sig2paramfull.HB.st = list()  # for S - full data
list_xiparamfull.HB.st = list()
#> w/o cluster ++++++++++++++++++++++++++++++++++++++ (covariate)
# N/A
# N/A
# N/A

#-----------------------------------------------------------
#-----------------------------------------------------------
#-----------------------------------------------------------
#-----------------------------------------------------------
#-----------------------------------------------------------

####> with cluster-wise +++++++++++++++++++++++++++++++ (outcome) 
list_betaparam.HB.st = list()    #for S
list_sig2param.HB.st = list()    #for S
list_xiparam.HB.st = list()

# -------------------- after correction ---------------------- #
list_betaparam.HB.clean = list()    
list_sig2param.HB.clean = list()    
list_xiparam.HB.clean = list()

####> with cluster-wise +++++++++++++++++++++++++++++ (covariate) 
list_piparam.HB.st = list()       #for ZS
list_lambda2param.HB.st = list()  #for XS
list_kappaparam.HB.st = list()             #for XS

list_tau2param.HB = list()         #for measurement model

# -------------------- after correction ---------------------- #
list_piparam.HB.clean = list()        #for Z
list_lambda2param.HB.clean = list()
list_kappaparam.HB.clean = list()     #for X|Z

###########################################################################################################################


#--------------------------------------------------------------------- ALWAYS KEEP and NO-update !!!! --------------------------
### What we have so far................................................................................ for {initialization}
###] Prepare for the current versions of single "cluster-wise" parameters as the "previous" +++++++++++++++++++++++++++++++++

#> hyperparameter
beta0_prev.HB.st = beta0param.HB.st
SIG_b0_prev.HB.st = SIG_b0param.HB.st
u0_prev.HB.st = u0param.HB.st
v0_prev.HB.st = v0param.HB.st

#> outcome
beta_prev.HB.st = betaparam.HB.st
sig2_prev.HB.st = sig2param.HB.st
xi_prev.HB.st = xiparam.HB.st
#beta_prev_fullS.HB.st = rmvn(n=1, mu=beta0_prevS.HB.st, sigma=SIG_b0_prevS.HB.st) # prebviously given!!!!
#sig2_prev_full.HB.st = rgamma(n=1, shape=u0_prevS.HB.st, rate=v0_prevS.HB.st)     # prebviously given!!!!


#> covariate ---------------------------------------------------------
piparam_prev.HB.st = piparam.HB.st           # we do not need. unaffected.
piparam_prev.HB.clean = piparam.HB.st        # we do not need. unaffected.

lambda2param_prev.HB.st = lambda2param.HB.st # do we need???? 
lambda2param_prev.HB.clean = lambda2param.HB.st #????

kappaparam_prev.HB.st = kappaparam.HB.st
kappaparam_prev.HB.clean = kappaparam.HB.st

tau2_prev = tau2param.HB  # tau2param:::(1-zeta_scale)*lambda2paramS.HB.st
#---------------------------------------------------------------------
#--------------------------------------------------------------------- ALWAYS KEEP and NO-update !!!! --------------------------


###] Miscellaneous ???

max_iter = 10   # let's say.. investigation( sampling rejection ) limit...?
counts_mat <- matrix(0, nrow = J.HB*total_iter, ncol = 3) # count..How many rejected???

n_j.HB.st = table(cl_membership.HB)

beta_p_accepted.HB.st <- 0
sig2_p_accepted.HB.st <- 0
xi_p_accepted.HB.st <- 0

beta_p_accepted.HB.clean <- 0
sig2_p_accepted.HB.clean <- 0
xi_p_accepted.HB.clean <- 0

#diff_lambda2.HB <- matrix(0, nrow = n.train, ncol = total_iter)
#diff_tau2.HB <- matrix(0, nrow = n.train, ncol = total_iter)

n_j.HB.st # just temporary...for cluster-wise sampling...

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
################################## START #######################################

for (r in 1:total_iter) {
  
  ####% for Xstar - INITIALIZE covariate parameters 
  kappaparamclean = matrix(0, nrow = J.HB, ncol = ncol(kappaparam.HB.st))     #for X|Z vs Xstar|Z
  lambda2paramclean = numeric(J.HB)                                        #for X|Z vs Xstar|Z
  
  ####% for S - INITIALIZE outcome parameters
  betaparamclean =  matrix(0, nrow = J.HB, ncol = ncol(betaparam.HB.st))      #for S
  sig2paramclean = numeric(J.HB)                                           #for S
  xiparamclean = numeric(J.HB)                                            #for S
  
  xclean = numeric(n.train)                                             #for X
  
  
  
  
  ################################ Real GAME ###################################
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  for(j in 1:J.HB) {
    
    Xstar_j=Xstar[cl_membership.HB==j]
    Z_j = Z[cl_membership.HB==j]
    
    ##############################################################
    ##############################################################
    ##### > To meet the conditions for the SYSTEM EQUATION < #####
    difference1 = -1
    difference2 = -1
    difference3 = -1
    
    count_iterations = 0
    
    ###> Don't stop until all the three positive values are collected.  
    while( difference1 < 0 | difference2 < 0 | difference3 < 0) {
      
      ###[1] Updating Posterior -------------------------------------------------------------------------------
      
      #>>> First, sample ***HYPERPARAMETERS*** from hyperPosterior based on full data {w/o clusters}
      
      ## hyperparam for [betaparamS]
      SIG_b0param.HB.st = riwish(v=qu0.HB.st+2, 
                                 S=t(beta0_prev.HB.st-beta_prev_full.HB.st) %*% (beta0_prev.HB.st-beta_prev_full.HB.st) +
                                   delta*t(beta0_prev.HB.st - m0.HB.st) %*% (beta0_prev.HB.st - m0.HB.st) + LAMB.HB.st) 
      
      beta0param.HB.st = rmvn(n=1, mu=delta/(delta+1)*m0.HB.st + 1/(delta+1)*beta_prev_full.HB.st, 
                              sigma=1/(delta+1)*SIG_b0param.HB.st) 
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~(CHECK)
      
      ## hyperparam for [sig2param]  
      #> - Sample u0 from MH using Gamma(3,1) as proposal 
      u0_prop.HB.st = rgamma(n=1, shape=u0_proposed_shape.HB.st, rate=u0_proposed_rate.HB.st)
      ratio_u0.HB.st = min(fink(u0_prop.HB.st, A = sig2_prev_full.HB.st*rho_u1.HB.st/(v0_prev.HB.st/2), B = rho_u2.HB.st+1)/
                             fink(u0_prev.HB.st, A = sig2_prev_full.HB.st*rho_u1.HB.st/(v0_prev.HB.st/2), B = rho_u2.HB.st+1)*
                             dgamma(x = u0_prev.HB.st, shape = u0_proposed_shape.HB.st, rate = u0_proposed_rate.HB.st)/
                             dgamma(x = u0_prop.HB.st, shape = u0_proposed_shape.HB.st, rate = u0_proposed_rate.HB.st), 1)
      if(runif(n = 1, min = 0, max = 1) < ratio_u0.HB.st) {
        u0param.HB.st = u0_prop.HB.st
      } else {
        u0param.HB.st = u0_prev.HB.st
      }
      
      #> - Sample v0 from Gamma distribution 
      v0param.HB.st = rgamma(n=1, shape = rho_v1.HB.st + u0param.HB.st/2, rate = rho_v2.HB.st+1/sig2_prev_full.HB.st)
      
      
      #------------------------------------------------------------------------------------------------------------------
      #------------------------------------------------------------------------------------------------------------------
      
      #>>> Next I, Sample ****[Main parameters]**** 
      
      #> - Proposal Sampling using proposals OF PRIORS {w/o clusters} 
      beta_prop_full.HB.st = rmvn(n=1, mu=beta0param.HB.st, sigma=SIG_b0param.HB.st*varinf) #~~~~~~~~~~~~~~~~~~~(CHECK)
      sig2_prop_full.HB.st = rinvgamma(n=1, shape=u0param.HB.st, scale=v0param.HB.st)
      xi_prop_full.HB.st = rt(n=1, df=nu0.HB.st)
      
      #> - Sifting a single [betaparam] based on FULL data 
      numerator=0
      denominator=0
      # plugging prior sample first...
      for (i in 1:n.train){
        numerator = numerator + 
          dlogsknorm_log(S[i], matXstar[i,], t(as.matrix(beta_prop_full.HB.st)), sig2_prev_full.HB.st, xi_prev_full.HB.st)
        
        denominator = denominator + 
          dlogsknorm_log(S[i], matXstar[i,], t(as.matrix(beta_prev_full.HB.st)), sig2_prev_full.HB.st, xi_prev_full.HB.st)
      }
      ratio_beta.HB.st = min(exp(numerator-denominator), 1)
      
      U = runif(n = 1, min = 0, max = 1)
      if(U < ratio_beta.HB.st) {
        betaparam_full.HB.st = beta_prop_full.HB.st
      } 
      else {
        betaparam_full.HB.st = beta_prev_full.HB.st
      }
      
      
      #> - Sifting a single [sig2param] based on FULL data 
      numerator=0
      denominator=0
      # plugging prior sample first...
      for (i in 1:n.train){
        numerator = numerator + 
          dlogsknorm_log(S[i], matXstar[i,], t(as.matrix(betaparam_full.HB.st)), sig2_prop_full.HB.st, xi_prev_full.HB.st)
        
        denominator = denominator + 
          dlogsknorm_log(S[i], matXstar[i,], t(as.matrix(betaparam_full.HB.st)), sig2_prev_full.HB.st, xi_prev_full.HB.st)
      }
      ratio_sig2.HB.st = min(exp(numerator-denominator), 1)
      
      U = runif(n = 1, min = 0, max = 1)
      if(U < ratio_sig2.HB.st) {
        sig2param_full.HB.st = sig2_prop_full.HB.st
      } 
      else {
        sig2param_full.HB.st = sig2_prev_full.HB.st
      }
      
      
      #> - Sifting a single [xiparam] based on FULL data 
      numerator=0
      denominator=0
      # plugging prior sample first...
      for (i in 1:n.train){
        numerator = numerator + 
          dlogsknorm_log(S[i], matXstar[i,], t(as.matrix(betaparam_full.HB.st)), sig2param_full.HB.st, xi_prop_full.HB.st)
        denominator = denominator + 
          dlogsknorm_log(S[i], matXstar[i,], t(as.matrix(betaparam_full.HB.st)), sig2param_full.HB.st, xi_prev_full.HB.st)
      }
      # compute the ratio
      ratio_xi.HB.st = min(exp(numerator-denominator), 1)
      
      U = runif(n = 1, min = 0, max = 1)
      if(U < ratio_xi.HB.st) {
        xiparam_full.HB.st = xi_prop_full.HB.st
      } 
      else {
        xiparam_full.HB.st = xi_prev_full.HB.st
      }
      
      
      
      
      #--------------------------------------------------------------------------------------------------------#
      #--------------------------------------------------------------------------------------------------------#
      #-----------------------------------Now it's time for cluster-wise---------------------------------------#
      #--------------------------------------------------------------------------------------------------------#
      #--------------------------------------------------------------------------------------------------------#
      
      
      #>>> First.., Sample **[covariate parameters_j]**
      
      #Xstar_j=Xstar[cl_membership.HB==j] #: already defined
      #Z_j = Z[cl_membership.HB==j]        #: already defined
      #n_j.st[j]                         #: already defined
      
      #> - proposal Sampling [piparamF_j] using proposals OF POSTERIOR {with clusters} 
      piparam.HB.st[j] = rbeta( n=1, shape1=g0.HB.st+sum(Z_j), shape2=h0.HB.st+n_j.HB.st[j]-sum(Z_j) )
      # no need...previous value...?
      
      #> - proposal Sampling [lambda2paramF_j] using proposals OF POSTERIOR {with clusters} 
      lambda2param.HB.st[j] = rinvgamma( n=1, 
                                         shape=(c0.HB.st+n_j.HB.st[j])/2, 
                                         scale=(d0.HB.st + 
                                                  sum( (Xstar_j-kappaparam_prev.HB.st[j,1]-kappaparam_prev.HB.st[j,2]*Z_j)^2 ))/2 )
      
      #> - proposal Sampling [tau2param_j] using proposals OF POSTERIOR {with clusters}  
      # + VERSION 01 +                              
      # tau2param[j] = rinvgamma( n=1, 
      #                           shape=(c0S.HB.st+nS_j.HB.st[j])/2, 
      #                           scale=(1-zeta_scale)*(d0S.HB.st + 
      #                                    sum( (XSstar_j-kappaparam_prev.HB.st[j,1]-kappaparam_prev.HB.st[j,2]*ZSj)^2 ))/2
      # )
      # + VERSION 02 + 
      tau2param.HB[j] = rinvgamma( n=1, 
                                   shape=(c0.HB.st+n_j.HB.st[j])/2, 
                                   scale=(1-zeta_scale)*(d0.HB.st + 
                                                           sum( (Xstar_j-kappaparam_prev.HB.st[j,1]-kappaparam_prev.HB.st[j,2]*Z_j)^2 ))/2)
      #> - proposal Sampling [kappaparam]_j
      KK1.HBj = cbind(1, Z_j)
      KK2.HBj = matrix(c(sum(Xstar_j), sum(Xstar_j*Z_j)), nrow=2, ncol=1) # 2x1 vector
      SIG_inv_j = solve(SIG_inv_k0.HB+(t(KK1.HBj) %*% KK1.HBj)) # solve( ..[matrix].. ) = [matrix]^(-1)
      kappaparam.HB.st[j,] <- rmvn(n=1, mu=SIG_inv_j %*% (SIG_inv_k0.HB %*% kappa_v0.HB + KK2.HBj), 
                                   sigma = lambda2param.HB.st[j]*SIG_inv_j) #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~(CHECK)
      
      
      #---------------------------------------------------------------------------------------------------------------
      #---------------------------------------------------------------------------------------------------------------
      #>>> Next,... Sample ****[Main parameters_j]******************************
      
      #>>> subsetting by cluster
      S_j = S[cl_membership.HB==j]                  #outcome
      matX_j.HB.st = matXstar[cl_membership.HB==j, ] #covariate
      
      #>>> Sample "new" outcome paramters from PROPOSALS (priors)
      beta_prop.HB.st = rmvn(n=1, mu=beta0param.HB.st, sigma=SIG_b0param.HB.st*varinf) #~~~~~~~~~~~~~~~~~~~~~~~~~(CHECK)
      sig2_prop.HB.st = rinvgamma(n=1, shape=u0_new.HB.st, scale=v0_new.HB.st)
      xi_prop.HB.st = rt(n=1, df=nu0.HB.st)
      
      
      #>>> - Sifting [betaparam]_j
      numerator=0
      denominator=0
      for (i in 1:length(S_j)){
        numerator = numerator + 
          dlogsknorm_log(S_j[i], matX_j.HB.st[i,], t(as.matrix(beta_prop.HB.st)), sig2_prev.HB.st[j], xi_prev.HB.st[j])
        
        denominator = denominator + 
          dlogsknorm_log(S_j[i], matX_j.HB.st[i,], t(as.matrix(beta_prev.HB.st[j,])), sig2_prev.HB.st[j], xi_prev.HB.st[j])
      }
      ratio_beta.HB.st = min(exp(numerator-denominator), 1)
      U = runif(n = 1, min = 0, max = 1)
      if(ratio_beta.HB.st > U) {
        betaparam.HB.st[j, ] = beta_prop.HB.st
        beta_p_accepted.HB.st <- beta_p_accepted.HB.st+1
      } else {
        betaparam.HB.st[j,] = beta_prev.HB.st[j,]
      }
      
      
      #>>> - Sifting [sig2param]_j
      numerator=0
      denominator=0
      for (i in 1:length(S_j)){
        numerator = numerator + 
          dlogsknorm_log(S_j[i], matX_j.HB.st[i,], t(as.matrix(betaparam.HB.st[j,])), sig2_prop.HB.st, xi_prev.HB.st[j]) +
          log(dinvgamma(sig2_prop.HB.st, shape=u0param.HB.st, scale=v0param.HB.st)) + 
          log(dinvgamma(sig2_prev.HB.st[j], shape=u0_new.HB.st, scale=v0_new.HB.st))
        
        denominator = denominator + 
          dlogsknorm_log(S_j[i], matX_j.HB.st[i,], t(as.matrix(betaparam.HB.st[j,])), sig2_prev.HB.st[j], xi_prev.HB.st[j]) +
          log(dinvgamma(sig2_prev.HB.st[j], shape=u0param.HB.st, scale=v0param.HB.st)) + 
          log(dinvgamma(sig2_prop.HB.st, shape=u0_new.HB.st, scale=v0_new.HB.st))
      }
      ratio_sig2.HB.st = min(exp(numerator-denominator), 1)
      U = runif(n = 1, min = 0, max = 1)
      if(U < ratio_sig2.HB.st) {
        sig2param.HB.st[j] = sig2_prop.HB.st
        sig2_p_accepted.HB.st <- sig2_p_accepted.HB.st+1
      }
      else {
        sig2param.HB.st[j] = sig2_prev.HB.st[j]
      }
      
      
      #>>> - Sifting [xiparam]_j
      numerator=0
      denominator=0
      for (i in 1:length(S_j)){
        numerator = numerator + 
          dlogsknorm_log(S_j[i], matX_j.HB.st[i,], t(as.matrix(betaparam.HB.st[j,])), sig2param.HB.st[j], xi_prop.HB.st)
        denominator = denominator + 
          dlogsknorm_log(S_j[i], matX_j.HB.st[i,], t(as.matrix(betaparam.HB.st[j,])), sig2param.HB.st[j], xi_prev.HB.st[j])
      }
      # compute the ratio
      ratio_xi.HB.st = min(exp(numerator-denominator), 1)
      
      U = runif(n = 1, min = 0, max = 1)
      if(U < ratio_xi.HB.st) {
        xiparam.HB.st[j] = xi_prop.HB.st
        xi_p_accepted.HB.st <- xi_p_accepted.HB.st+1
      } 
      else {
        xiparam.HB.st[j] = xi_prev.HB.st[j]
      }
      
      
      #################################################################################################################
      #################################################################################################################
      #################################################################################################################
      #################################################################################################################
      ###> [2] Calculate differences to ensure VALIDITY of [System Equation]
      #################################################################################################################
      #>>> [main] I. For clean beta_j_1 (the problematic***slope on clean X)
      betaparamclean_j = betaparam.HB.st[j,2]*lambda2param.HB.st[j]/(lambda2param.HB.st[j]-tau2param.HB[j])
      
      #>>> [main] II. For clean sig2_j (the problematic***variance param for LSN on clean X)
      sig2paramclean_j = sig2param.HB.st[j] - betaparamclean_j^2/( 1/tau2param.HB[j] + 1/(lambda2param.HB.st[j]-tau2param.HB[j]) )
      
      #>>> [main] III. For clean xi_j
      # ?
      
      
      #>>> [supplementary] How to ensure They are positive ? 
      difference1 = lambda2param.HB.st[j]-tau2param.HB[j] # ensure clean lambda2 is positive
      difference2 = sig2paramclean_j
      difference3 = 1/tau2param.HB[j]+1/(lambda2param.HB.st[j]-tau2param.HB[j]) - xiparam.HB.st[j]^2*betaparamclean_j^2/sig2paramclean_j
      
      
      #> -while looping count- iteration for each j
      print(c(difference1, difference2, difference3))
      count_iterations = count_iterations+1
    } #%%% END of WHILE %%% 
    print(paste("count=", count_iterations))
    
    # Values that meet the condition are obtained. Now....
    ###> Collect the params and compute the clean parameters using [System Equation].
    #################################################################################################################
    
    #kappaparamclean = matrix(0, nrow = JS, ncol = ncol(kappaparam.st))     # already defined
    #lambda2paramcleanS = numeric(JS)                                       # already defined
    #betaparamcleanS =  matrix(0, nrow = JS, ncol = ncol(betaparamS.st))    # already defined
    #sig2paramclean =numeric(JS)                                           # already defined
    
    #xSclean = numeric(n.train)                                             # already defined
    
    
    lambda2paramclean[j] = lambda2param.HB.st[j] - tau2param.HB[j]    
    
    kappaparamclean[j,] = kappaparam.HB.st[j,]
    
    betaparamclean[j,2] = betaparam.HB.st[j,2]*lambda2param.HB.st[j]/(lambda2param.HB.st[j]-tau2param.HB[j]) #for Beta1
    
    betaparamclean[j,1] = betaparam.HB.st[j,1] - 
      betaparamclean[j,2]*kappaparam.HB.st[j,1]*tau2param.HB[j]/(lambda2param.HB.st[j]-tau2param.HB[j]) #for Beta0
    
    betaparamclean[j,3] = betaparam.HB.st[j,3] - 
      betaparamclean[j,2]*kappaparam.HB.st[j,2]*tau2param.HB[j]/(lambda2param.HB.st[j]-tau2param.HB[j]) #for Beta2
    
    sig2paramclean[j] = sig2param.HB.st[j] - 
      betaparamclean[j,2]^2/( 1/tau2param.HB[j] + 1/lambda2paramclean[j] )   #for sig2
    
    xiparamclean[j] = sign(xiparam.HB.st[j])*sqrt( xiparam.HB.st[j]^2*(betaparamclean[j,2]^2/sig2paramclean[j]+1/tau2param.HB[j]+1/lambda2paramclean[j])/
                                                     (1/tau2param.HB[j]+1/lambda2paramclean[j] - xiparam.HB.st[j]^2*betaparamclean[j,2]^2/sig2paramclean[j]) )
    
    #################################################################################################################
  } #%%% END of [j] %%%
  
  
  #### For NEXT!!!
  ###> 1. Store the clean data parameters
  ###> 2. Set previous set of parameter to the [[[current]]] for the next iteration
  ## - Manner! main param for next iteration
  
  ##]]]]]] - Let's SEE..
  #-beta0_prevS.st = beta0paramS.st #-------------------------------------------------already defined
  #-SIG_b0_prevS.st = SIG_b0paramS.st #-----------------------------------------------already defined
  #-u0_prevS.st = u0paramS.st #-------------------------------------------------------already defined
  #-v0_prevS.st = v0paramS.st #-------------------------------------------------------already defined
  
  #-beta_prevS.st = betaparamS.st #---------------------------------------------------already defined
  #-sig2_prev.st = sig2param.st #-----------------------------------------------------already defined
  #-beta_prev_fullS.st = rmvn(n=1, mu=beta0_prevS.st, sigma=SIG_b0_prevS.st) # prebviously given!!!!
  #-sig2_prev_full.st = rgamma(n=1, shape=u0_prevS.st, rate=v0_prevS.st)     # prebviously given!!!!
  
  #piparam_prevS.st = piparamS.st           # we do not need. unaffected.
  #piparam_prevS.clean = piparamS.st        # we do not need. unaffected.
  
  #-lambda2param_prevS.st = lambda2paramS.st # do we need???? #-----------------------already defined
  #-lambda2param_prevS.clean = lambda2paramS.st #???? #-------------------------------already defined
  #-kappaparam_prev.st = kappaparam.st #----------------------------------------------already defined
  #-kappaparam_prev.clean = kappaparam.st #-------------------------------------------already defined
  
  #-tau2_prev = tau2param.HB  # tau2param.HB:::(1-zeta_scale)*lambda2paramS.st #------------already defined
  
  
  ##]]]]]] - The truth is.......                                     
  #>>>> WHY? ++++++++++++++++++++++++++++++++++++++++++ to sample [betaparamS] and [sig2param] next time....[O]: ONLY USEFUL
  beta0_prev.HB.st = beta0param.HB.st #:::::::::::::::: useful for proposal sample of "SIG_b0paramS.HB.st" 
  SIG_b0_prev.HB.st = SIG_b0param.HB.st
  u0_prev.HB.st = u0param.HB.st #::::::::::::::::::::::::::::::::::::::::: useful to sift "u0paramS.HB.st"
  v0_prev.HB.st = v0param.HB.st #::::::::::::::::::::::::::::::::::::::::: useful to sift "u0paramS.HB.st"
  
  #>>>> WHY? ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ to sample [Y] next time....but [X]
  beta_prev.HB.st = betaparam.HB.st # ::::::::::::::::::::::::::::::::::: useful to sift "betaparamS.HB.st"
  sig2_prev.HB.st = sig2param.HB.st # ::::::::::::::::::::::::::::::::::::: useful to sift "betaparamS.HB.st" and "sig2param.HB.st"
  xi_prev.HB.st = xiparam.HB.st
  #beta_prev_fullS.HB.st  # this is not updated!!!
  #sig2_prev_full.HB.st   # this is not updated!!!
  
  #>>>> WHY? +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ to sample [ZS] next time....but [X]
  #piparam_prevS.HB.st = piparamS.HB.st           
  #piparam_prevS.clean = piparamS.HB.st                            
  
  #>>>> WHY? ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ to sample [XS|ZS] next time....but [X]
  lambda2param_prev.HB.st = lambda2param.HB.st
  #lambda2param_prev.clean = lambda2paramclean #***********
  kappaparam_prev.HB.st = kappaparam.HB.st # ::::::::::::: useful for proposal sample of "lambda2paramS", and thus "tau2param.HB"     
  #kappaparam_prev.clean = kappaparamclean #*****************
  
  #>>>> WHY? ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ to sample [XSstar|XS] next time....but [X]                                       
  tau2_prev = tau2param.HB                                  
  
  
  
  
  ###############################################################################################################
  ###[3] Calculating Loglikelihood to monitor the chain convergence ---------------------------------------------
  ###############################################################################################################                                        
  #>>> First
  loglike_r = numeric(n.train)
  
  for(i in 1:n.train) {
    x_i = matXstar[i, ]      # first row..{1,x,z}
    j = cl_membership.HB[i] # membership ID
    
    loglike_r[i] = dlogsknorm_log(S[i], x_i, betaparam.HB.st[j,], sig2param.HB.st[j], xiparam.HB.st[j]) + 
      dnorm(x=x_i[2], mean=kappaparam.HB.st[j,1] + kappaparam.HB.st[j,2]*Z[i], sd=sqrt(lambda2param.HB.st[j]), log=TRUE) +
      dbinom(x=x_i[3], size=1, prob=piparam.HB.st[j], log=TRUE) 
  }  
  loglikelihood.HB.st[,r] = loglike_r
  
  if(r > r_convergence) {
    
    #>>> W/O CLUSTER
    list_beta0param.HB.st[[r-r_convergence]] = beta0param.HB.st
    list_SIG_b0param.HB.st[[r-r_convergence]] = SIG_b0param.HB.st
    list_u0param.HB.st[[r-r_convergence]] = u0param.HB.st
    list_v0param.HB.st[[r-r_convergence]] = v0param.HB.st
    
    list_betaparamfull.HB.st[[r-r_convergence]] = betaparam_full.HB.st
    list_sig2paramfull.HB.st[[r-r_convergence]] = sig2param_full.HB.st
    list_xiparamfull.HB.st[[r-r_convergence]] = xiparam_full.HB.st
    
    #>>> With Cluster-wise
    list_piparam.HB.st[[r-r_convergence]] = piparam.HB.st     
    list_lambda2param.HB.st[[r-r_convergence]] = lambda2param.HB.st
    list_kappaparam.HB.st[[r-r_convergence]] = kappaparam.HB.st
    list_tau2param.HB[[r-r_convergence]] = tau2param.HB
    
    list_betaparam.HB.st[[r-r_convergence]] = betaparam.HB.st  
    list_sig2param.HB.st[[r-r_convergence]] = sig2param.HB.st
    list_xiparam.HB.st[[r-r_convergence]] = xiparam.HB.st
    
    ################# STORE the CLEAN parameters ##################
    list_kappaparam.HB.clean[[r-r_convergence]] = kappaparamclean     #for Xstar
    list_lambda2param.HB.clean[[r-r_convergence]] = lambda2paramclean    #for Xstar
    list_betaparam.HB.clean[[r-r_convergence]] = betaparamclean    #for S
    list_sig2param.HB.clean[[r-r_convergence]] = sig2paramclean    #for S
    list_xiparam.HB.clean[[r-r_convergence]] = xiparamclean        #for S
    
    # list_xclean[[r-r_convergence]] = xclean #~xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx Q. (can we impute?)
    
  }
  print(paste("r=",r))
  trace(what = "print", where = getNamespace("base"), exit = flush.console, print = FALSE)
}

########################################
####### End of Gibbs Sampler ###########
########################################
beta_p_accepted.HB.st/(total_iter*J.HB)
sig2_p_accepted.HB.st/(total_iter*J.HB)
xi_p_accepted.HB.st/(total_iter*J.HB)

par(mfrow=c(1,1))
plot(colSums(loglikelihood.HB.st),type="l")

lppd_EX_HB.S.st = mean( colSums(loglikelihood.HB.st) ); lppd_EX_HB.S.st   # -23107.04(0.95), 






#-------------------------------------------------------------------------------

############################### Predictions ####################################

# Expected Value of LSN distribution (Wang 2019)
# 2*exp(sum(x_i*beta_j[j,]) + sig2_j[j]/2)*(1-pnorm(-xi_j[j]*sqrt(sig2_j[j])/sqrt(xi_j[j]^2+1)))


###> Training data
expval.HB.train.st = matrix(0, nrow = n.train, ncol = total_iter-r_convergence)

for(r in 1:(total_iter-r_convergence)) {
  betaparamclean = list_betaparam.HB.st[[r]]
  sig2paramclean = list_sig2param.HB.st[[r]]
  xiparamclean = list_xiparam.HB.st[[r]]
  
  for(i in 1:n.train) {
    j = cl_membership.HB[i]
    #cleanx = rnorm(n = 1, mean = kappaparamclean[j,1] + kappaparamclean[j,2]*Z[i], sd = sqrt(lambda2paramclean[j]))
    cleanx = Xstar[i]
    #expval.train[i,r] <- 2*exp( betaparamclean[j,1] + betaparamclean[j,2]*cleanx + betaparamclean[j,3]*Z[i] + 
    #                              sig2paramclean[j]/2)*( 1-pnorm(-xiparamclean[j]*
    #                                                               sqrt(sig2paramclean[j])/sqrt(xiparamclean[j]^2+1)))
    expval.HB.train.st[i,r] <- 20*exp( betaparamclean[j,1] + betaparamclean[j,2]*cleanx + betaparamclean[j,3]*Z[i] - 
                                         sig2paramclean[j]/2)*( 1-pnorm(-xiparamclean[j]*
                                                                          sqrt(sig2paramclean[j])/sqrt(xiparamclean[j]^2+1)))
  }
}

expS.HB.train.st <- apply(X=expval.HB.train.st, MARGIN=1, FUN=mean)

summary(S)
summary(expS.HB.train.st)

format( aggregate(log(S), by=list(cl_membership.HB), FUN=summary), scientific=FALSE)
format( aggregate(log(expS.HB.train.st), by=list(cl_membership.HB), FUN=summary), scientific=FALSE)

# Clusterwise predictive density plots
# par(mfrow = c(2,3))
# 
# for(j in 1:J.HB) {
#   hist(log(S)[cl_membership.HB==j], freq=FALSE, breaks=50, col="white", xlim=c(5,15), ylim=c(0,1.6),
#        xlab = paste("Cluster", j))
#   lines(density(log(expS.HB.train.st)[cl_membership.HB==j]), col="red", ylim=c(0,1.6), lwd=2)
# }


###> Testing data
expval.HB.test.st = matrix(0, nrow = n.test, ncol = total_iter-r_convergence)

for(r in 1:(total_iter-r_convergence)) {
  betaparamclean = list_betaparam.HB.st[[r]]
  sig2paramclean = list_sig2param.HB.st[[r]]
  xiparamclean = list_xiparam.HB.st[[r]]
  
  for(i in 1:n.test) {
    j = cl_membership.HB[i]
    #cleanx = rnorm(n = 1, mean = kappaparamclean[j,1] + kappaparamclean[j,2]*Z[i], sd = sqrt(lambda2paramclean[j]))
    cleanx = Xstar.test[i]
    expval.HB.test.st[i,r] <- 20*exp( betaparamclean[j,1] + betaparamclean[j,2]*cleanx + betaparamclean[j,3]*Z[i] - 
                                        sig2paramclean[j]/2)*( 1-pnorm(-xiparamclean[j]*
                                                                         sqrt(sig2paramclean[j])/sqrt(xiparamclean[j]^2+1)))
  }
}
expS.HB.test.st <- apply(X=expval.HB.test.st, MARGIN=1, FUN=mean)

summary(S.test)
summary(expS.HB.test.st)

format( aggregate(log(S.test), by=list(cl_membership.HB.test), FUN=summary), scientific=FALSE)
format( aggregate(log(expS.HB.test.st), by=list(cl_membership.HB.test), FUN=summary), scientific=FALSE)

# Clusterwise predictive density plots
# par(mfrow = c(2,3))
# 
# for(j in 1:J.HB) {
#   hist(log(S.test)[cl_membership.HB.test==j], freq=FALSE, breaks=50, col="white", xlim=c(5,15), ylim=c(0,1.6),
#        xlab = paste("Cluster", j))
#   lines(density(log(expS.HB.test.st)[cl_membership.HB.test==j]), col="red", ylim=c(0,1.6), lwd=2)
# }




################### predictive density plots #######################
par(mfrow = c(1,1))

###> with train
hist(log(S), breaks=n.breaks.train, freq=F, xlab ="log(S)", main="Predictive loss density for a policy", 
     col="white", ylim=c(0,0.7), cex.axis=2, cex.lab=2)

lines(density( log(expS.DP.Avg.train) ), col="red", lwd=8)
lines(density( log(expS.DP.Avg.train.xx) ), col="black", lty = "dotted", lwd=2) # dirty
lines(density( log(expS.DP.Avg.train.st) ), col="blue", lwd=7) # dirty with correction
lines(density( log(expS.HB.train.st) ), col="green", lwd=4) # dirty with correction

legend("topright", legend = c("Gold standard(DP)", "Before correction(DP)", "After correction(DP)", "After correction(HB)"), #, "Correction Medium", "Correction Small"),
       col = c("red","black","blue","green"), #, "darkmagenta", "darkgreen"), 
       lty = c(1, 2, 1, 1), 
       lwd = c(8, 2, 7, 4), 
       cex = 0.8
)


###> with test
hist(log(S.test), breaks=n.breaks.test, freq=F, xlab ="log(S)", main="Predictive loss density for a policy", 
     col="white", ylim=c(0,0.4))

lines(density( log(expS.DP.Avg.test) ), col="red", lwd=5)
lines(density( log(expS.DP.Avg.test.xx) ), col="black", lty = "dotted") # dirty
lines(density( log(expS.DP.Avg.test.st) ), col="blue", lwd=1) # dirty with correction
lines(density( log(expS.HB.test.st) ), col="green", lwd=1) # dirty with correction

legend("topright", legend = c("Gold standard(DP)", "Before correction(DP)", "After correction(DP)", "After correction(HB)"), #, "Correction Medium", "Correction Small"),
       col = c("red","black","blue","green"), #, "darkmagenta", "darkgreen"), 
       lty = c(1, 2, 1, 1), 
       lwd = c(3, 1, 1, 1), 
       cex = 0.8
)




################################################################################
lppd_EX_DP.S    #-22207.38 
#                - 
lppd_EX_DP.S.xx #-22667.67
#                - 
lppd_EX_DP.S.st #-33475.79 (0.9); -21604.54 (0.95); -21954.74 (0.99)

#lppd_EX_HB.S
#lppd_EX_HB.S.xx
lppd_EX_HB.S.st #-23107.04





###> SSSE ----------------------------------------------------------------------
##> The standard
expS.DP.Avg.test
SSPE.DP.LSN <- sum( (log(expS.DP.Avg.test) - log(S.test))^2 ); SSPE.DP.LSN   #269.1253
SAPE.DP.LSN <- sum( abs(log(expS.DP.Avg.test) - log(S.test)) ); SAPE.DP.LSN  #196.2275
#expS.HB.test
#SSPE.HB.LSN <- sum( (log(expS.HB.test) - log(S.test))^2 ); SSPE.HB.LSN   #
#SAPE.HB.LSN <- sum( abs(log(expS.HB.test) - log(S.test)) ); SAPE.HB.LSN  #


##> Modelrisk
expS.DP.Avg.test.xx
SSPE.DP.LSN.xx <- sum( (log(expS.DP.Avg.test.xx) - log(S.test))^2 ); SSPE.DP.LSN.xx   #276.0281
SAPE.DP.LSN.xx <- sum( abs(log(expS.DP.Avg.test.xx) - log(S.test)) ); SAPE.DP.LSN.xx  #199.5169
#expS.HB.test.xx
#SSPE.HB.LSN.xx <- sum( (log(expS.HB.test.xx) - log(S.test))^2 ); SSPE.HB.LSN.xx   #
#SAPE.HB.LSN.xx <- sum( abs(log(expS.HB.test.xx) - log(S.test)) ); SAPE.HB.LSN.xx  #



##> Gustafson Correction (0.95)
expS.DP.Avg.test.st
SSPE.DP.LSN.st <- sum( (log(expS.DP.Avg.test.st) - log(S.test))^2 ); SSPE.DP.LSN.st  #274.619 (0.9); 285.6494 (0.95); 277.4901 (0.99)
SAPE.DP.LSN.st <- sum( abs(log(expS.DP.Avg.test.st) - log(S.test)) ); SAPE.DP.LSN.st  #198.8607 (0.9); 204.0293 (0.95); 201.5277 (0.99)
#                                                                                                     
expS.HB.test.st
SSPE.HB.LSN.st <- sum( (log(expS.HB.test.st) - log(S.test))^2 ); SSPE.HB.LSN.st   # 246.0146
SAPE.HB.LSN.st <- sum( abs(log(expS.HB.test.st) - log(S.test)) ); SAPE.HB.LSN.st  # 189.5098
# small error -> HB is better that DP????






##> Gustafson Correction (0.8)














###> CTE (for S) ---------------------------------------------------------------
#summary( expS.HB.train )
summary( expS.HB.train.st )
summary( expS.DP.Avg.train )
summary( expS.DP.Avg.train.xx )
summary( expS.DP.Avg.train.st )

#summary( expS.HB.test )
summary( expS.HB.test.st )
summary( expS.DP.Avg.test )
summary( expS.DP.Avg.test.xx )
summary( expS.DP.Avg.test.st )

# - Calculate the quantile threshold
alalpha = 0.95
alalpha = 0.9
alalpha = 0.5
alalpha = 0.1

#quantile_threshold.HB <- quantile(expS.HB.train, alalpha); quantile_threshold.HB
quantile_threshold.HB.st <- quantile(expS.HB.train.st, alalpha); quantile_threshold.HB.st
#tail_data <- expS.HB.train[expS.HB.train > quantile_threshold]; tail_data # NA?????
tail_data.HB.st <- expS.HB.train.st[expS.HB.train.st > quantile_threshold.HB.st]; tail_data.HB.st # NA?????
#CTE.HB <- mean(tail_data, na.rm = TRUE); CTE.HB          # 4300316  
CTE.HB.st <- mean(tail_data.HB.st, na.rm = TRUE); CTE.HB.st # 3854798

quantile_threshold.DP <- quantile(expS.DP.Avg.train, alalpha); quantile_threshold.DP
quantile_threshold.DP.xx <- quantile(expS.DP.Avg.train.xx, alalpha); quantile_threshold.DP.xx
quantile_threshold.DP.st <- quantile(expS.DP.Avg.train.st, alalpha); quantile_threshold.DP.st

tail_data.DP <- expS.DP.Avg.train[expS.DP.Avg.train > quantile_threshold.DP]; tail_data.DP # NA?????
tail_data.DP.xx <- expS.DP.Avg.train.xx[expS.DP.Avg.train.xx > quantile_threshold.DP.xx]; tail_data.DP.xx
tail_data.DP.st <- expS.DP.Avg.train.st[expS.DP.Avg.train.st > quantile_threshold.DP.st]; tail_data.DP.st # NA?????

CTE.DP <- mean(tail_data.DP, na.rm = TRUE); CTE.DP          #   
CTE.DP.xx <- mean(tail_data.DP.xx, na.rm = TRUE); CTE.DP.xx #
CTE.DP.st <- mean(tail_data.DP.st, na.rm = TRUE); CTE.DP.st # 














###> Calculate the KL divergence (for S) ---------------------------------------
library(entropy)

pdf <- density(expS.DP.Avg.train)$y
#pdf.HB <- density(expS.HB.train)$y
pdf.DP.xx <- density(expS.DP.Avg.train.xx)$y
pdf.DP.st <- density(expS.DP.Avg.train.st)$y
pdf.HB.st <- density(expS.HB.train.st)$y

kl.divergence.DP.xx <- entropy::KL.empirical(pdf, pdf.DP.xx)
print(kl.divergence.DP.xx) 

kl.divergence.DP.st <- entropy::KL.empirical(pdf, pdf.DP.st)
print(kl.divergence.DP.st) 

kl.divergence.HB.st <- entropy::KL.empirical(pdf, pdf.HB.st)
print(kl.divergence.HB.st) 


