################################# Rival for DP ################################# zeta:1250
library(mvnfast)
library(tidyverse)
library(Rcpp)
library(RcppArmadillo)
library(extraDistr)
library(MCMCpack)



##### Swedish Data with overall error: [40%]
full.train.sweden <- read.csv("C:/Users/kimm4/Desktop/WORKSPACE/2nd 3rd paper/DATASET/BSwed.paper3-40.train.csv")
full.test.sweden <- read.csv("C:/Users/kimm4/Desktop/WORKSPACE/2nd 3rd paper/DATASET/BSwed.paper3-40.test.csv")
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
X_tru.xx = full.train.sweden$ln_insured_err

matXstar = as.matrix(cbind(1, Xstar, Z))
matX_tru.xx = as.matrix(cbind(1, X_tru.xx, Z))

n.breaks.test = sqrt( nrow(full.train.sweden) ) #******Rule of thumb


#%%%% test %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
##### Variable Definitions
S.test = full.test.sweden$TotalLoss
logS.test = log(S.test)
Z.test = full.test.sweden$Experience

# dirty data on log scale, divided by standard dev to be on same scale
# Xstar = full.train.sweden$ln_insured_err/sd(full.train.sweden$ln_insured_err) 
Xstar.test = full.test.sweden$ln_insured_err
X_tru.xx.test = full.test.sweden$ln_insured

matXstar.test = as.matrix(cbind(1, Xstar.test, Z.test))
matX_tru.xx.test = as.matrix(cbind(1, X_tru.xx.test, Z.test))

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

#-------------------------------------------------------------------------------


####] Initial Clustering Based on Classes

# J = 7
# #> Hierarchical Clustering
# clust.out <- hclust( dist(cbind(log(S),Xstar,Z)) )
# plot(clust.out)
# table(cutree(clust.out,J))
# 
# 
# #> Kmeans Clustering
# clust.out.k <- kmeans(cbind(log(S),Xstar,Z), J)
# table(clust.out.k$cluster)

# clust.out.temp <- hclust( dist(cbind(S[1:100],Xstar[1:100],Z[1:100])) )
# plot(clust.out.temp)

#> Choose one for clustering: Risk, Zone, or k-means
# 1) Clustering Based on Risk
cl_membership <- full.train.sweden$Risk
table(cl_membership)

cl_membership.test <- full.test.sweden$Risk
table(cl_membership.test)

# # 2) Clustering Based on Zone
# cl_membership <- full.train.sweden$Zone
# table(cl_membership)

# # 3) Clustering Based on Hierarchical Cluster
# J=3
# clusters = cutree(hclust( dist(cbind(logS,Xstar,Z)) ), J)
# cl_membership = clusters 
# table(cl_membership)
# plot( hclust(dist(cbind(logS,Xstar,Z))) )
# rect.hclust(hclust(dist(logS,Xstar,Z)) , k = 3, border = 2:6)
# abline(h = 8, col = 'red')




################################################################################
###################### S ~ hierarchical LSN (clean) ############################
################################################################################

################################################################################
### Step01> Define prior, hyperprior model and hyperparameter values
################################################################################

# -------------------------------------------- Outcome ----- # ::: for S ~ LSN( X*betaparam, sig2param, xiparam )

# PRIOR: "betaparamF" ~ MVN(beta0_j, SIG2_b0_j)
# Hyperprior                beta0_j ~ MVN(m0, 1/delta*SIG_b0_j)
# Hyperprior                          SIG2_b0_j ~ IW( qu0, LAMB )

# ::: for "beta0_j" ~ MVN( m0, 1/delta*SIG_b0_j )
#> Gaussian regression for initialize OUTCOME model parameter m0, SIG_b0
fit1 <- glm( S ~  X_tru.xx + factor(Z), family=poisson(link = "log") )
# fit1 <- lm( logS ~ Xstar + factor(Z) )
summary(fit1)

m0 = coef(fit1)    # 1x3 initial reg_coeff vector (sort of "means"): Regression result as "mean"
SIG_b0 = vcov(fit1)   # 3x3 initial cov matrix of reg_coeff
p = length(m0)
options(scipen = 999)
SIG_b0
SIG_b0inv = solve(a=SIG_b0) # inverse of cov matrix of reg_coeff for later use (for posterior on reg_coeff)!!


#****** for VAR Inflation factor
delta = 0.01             # equal importance on m0 and beta ????
varinf = n.train/100 # rule of thumb divide by 5

#> SIG_b0
#- it requires -------
qu0 = p+2
# LAMB ?????????????????????????????????????????????
LAMB = SIG_b0



###> PRIOR: "sig2param" ~ IG( u0, v0 )
# Hyperprior               u0 ~ Fink( rho_u1, rho_u2 )
# Hyperprior                   v0 ~ Ga( rho_bv1, rho_v2 )

# Hyperparameter u0
rho_u1 = 1/8
rho_u2 = 3/2 

# Hyperparameter v0
rho_v1 = 8 
rho_v2 = 1 

# Based on the result? mean(unlist(list_u0param)), mean(unlist(list_v0param))
# E[invGa] = scale/(shape-1) = 50
u0_new = 1.7591 #3 #2 #1.1 # from 3/7  # from # summary(unlist(list_u0param)) # for u0_new                            
v0_new = 7.113  #3 #5 #8.5 # from 3/7  # from # summary(unlist(list_v0param)) # for v0_new

## [CHECK later] ---------------------------------------------------------------
# How can we design proposal for sig2 ~ InvGa(.) ?
par(mfrow=c(1,1))
plot( density(unlist(list_sig2param)) )                                   # sig2 we have...before
curve( dinvgamma(x, shape=mean(unlist(list_u0param)), 
                 scale=mean(unlist(list_v0param))), add=TRUE, col="red") # sig2 we obtain....after
curve( dinvgamma(x, shape=u0_new, 
                 scale=v0_new), add=TRUE, col="blue")                    # sig2 proposal adjustment

table(unlist(list_sig2param))
sort(table(unlist(list_sig2param)),DESC=TRUE)                             # to see if they are appropriate...
#-------------------------------------------------------------------------------

###> PRIOR: "xiparam"~t(loc, nu0, sca)
#loc = 0   # ? Location hyperparameter for T distribution
#nu0 = 1/2 # df for T distribution on pigtail parameter
#sca = 5   # ? Scale hyperparameter for T distribution 
nu0 = n.train-1
#- df refers to the number of independent observations (N-1):

# ------------------------------------------------------------------------------ #


#] ------ Covariate ------ ::: for Z ~ Bin( 1, [piparam] )  
# PRIOR                                         "piparam" ~ Beta(g0, h0)
#------------------------------------------------------------------------------------------------------------------------
g0 = 0.5 
h0 = 0.5 

#] ------ Covariate ------ ::: for X ~ N( X_bar, [lambda2param] ) where 
# PRIOR                                            "lambda2param" ~ IG(c0, d0)
#------------------------------------------------------------------------------------------------------------------------

#> X_bar 
#X_bar = mean(X_tru.xx) #????

#> lambda2
c0 = 0.5          
d0 = 0.5




# Now, we need a posterior sample of hyperPRIOR...
##################################################################################################
### Step02> initialize major parameters (single value), using analytically driven "POSTERIOR" pdf  
#                                                 Based on natural clustering
##################################################################################################
J = length(table(cl_membership))
set.seed(1)


#] ------ Covariate Z ------ ::: for Z ~ Bin( 1, [piparamS_j] ) 
# prior                                            "piparamS_j" ~ Beta(g0, h0)
#------------------------------------------------------------------------------------------------------------------------
# Covariate parameter from Posterior
piparam = numeric(J)
for (j in 1:J) {
  Zj=Z[cl_membership==j]
  nj=length(Zj)
  piparam[j] = rbeta( n=1, shape1=g0+sum(Zj), shape2=h0+nj-sum(Zj) ) 
}

piparam #% pi parameter sampled from posterior


#] ------ Covariate X ------ ::: for X ~ N( X_bar, [lambda2param] ) where 
# PRIOR                                            "lambda2param" ~ IG(c0, d0)
#------------------------------------------------------------------------------------------------------------------------
# Covariate parameter from Posterior
lambda2param = numeric(J)
for (j in 1:J) {
  Xj=X_tru.xx[cl_membership==j]
  n_j=length(Xj)
  lambda2param[j]=rinvgamma( n=1, shape=(c0+n_j)/2, scale=(d0 + sum( (Xj - mean(Xj))^2 ) )/2 )  
}                                                                                                     
lambda2param #% lambda2 parameter (posterior sample)


#% empirical X_barF ? using "aggregate(.)": investigation by splitting the data into subset
X_bar = aggregate(x=X_tru.xx, by=list(cl_membership), FUN=mean, na.rm= T)$x; X_bar
# ---------------------------------------------------------------- Perfect!!!!!!


############################## But Outcome parameter sampling is trickyyyyyyyyyy
#> This is where [new outcome parameters] go. 
betaparam = matrix(data=NA, nrow=J, ncol=3) 
sig2param = numeric(J) # for cluster-wise
sig2paramfull = 0       # for population
xiparam = numeric(J)


#] ------ Outcome S ------ ::: [betaparam], [sig2param], [xiparam].... try MH-algorithm
#------------------------------------------------------------------------------------------------------------------------

#####> First, sample hyperparameter from hyperPRIOR.
# Set in place the old hyperparameters..

### (A) beta0, SIG_b0 For [betaparam]
# 1) beta0 ~ MVN(m0, 1/delta*SIG_b0)  
# 2) SIG_b0 ~ ~ IW( qu0, LAMB )
# initial hyperparameter values from hyperprior +++++++++++++++++++++++(no need to care cluster)
# Set in place the old parameters..sample from prior...
set.seed(1)
SIG_b0_prev <- riwish(v=qu0, S=LAMB)                                         
beta0_prev <- rmvn(n=1, mu=m0, sigma=1/delta*SIG_b0_prev)                   #~~~~~~~~~~~~~~~(CHECK)


### (B) u0, v0 For [sig2param]
# 1) u0 ~ Fink(rho_u1, rho_u2)
# 2) v0 ~ Ga(rho_v1, rho_v2)

# since.... Fink function for ....... "u0" sampling 
fink = function(u0, A, B) {
  A^(u0/2-1)/gamma(u0/2)^B
}
u0_proposed_shape = 3.75 # 3
u0_proposed_rate = 2 # 2  based on Ga(3.75,2)...intuition..? for mean 3/2 ?

###> initial hyperparameter values from hyperprior +++++++++++++++++++++++(no need to care cluster)
u0_prev = u0_proposed_shape/u0_proposed_rate     # instead of sampling from Fink function???
v0_prev = rgamma(n=1, shape=rho_v1, rate=rho_v2)

#jpeg(file="plot.test.jpeg", width=1000, height=800)
#-------------------------------------------------------------------------------
# See? Fink function is close to Ga(3.75,2) ?
curve(fink(x,rho_u1,rho_u2)/integrate(function(x) fink(x,rho_u1,rho_u2),0,Inf)$value, 
      from = 0.01, to = 10, ylim = c(0,0.5))
curve(dgamma(x, shape = u0_proposed_shape, rate = u0_proposed_rate), add=TRUE, col="blue")
#-------------------------------------------------------------------------------


### (C) nu0 For [xiparam] ~ t(loc, nu0, sca)
# nu0 ~ Exponetial(.) ??? or constant? N-1 : dt(x, df, ncp, log = FALSE) : the central t distribution
curve( dt(x, df=nu0) )


#-------------------------------------------------------------------------------
##### Now....ready for Single SAMPLE FROM the PRIOR
###> initial [main parameter] values from PRIOR ++++++++++++++++++++++++++++++(cluster-wise, but...)
beta_prev_full = rmvn(n=1, mu=beta0_prev, sigma=SIG_b0_prev) #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~(CHECK)
beta_prev = matrix( rep(beta_prev_full, J), nrow=J, byrow= TRUE ) # cluster-wise!!!!

###> initial [main parameter] values from PRIOR ++++++++++++++++++++++++++++++(cluster-wise, but...)
sig2_prev_full = rinvgamma(n=1, shape=u0_prev, scale=v0_prev)   # population 
sig2_prev = rep(sig2_prev_full, J)                               # cluster-wise!!!!

###> initial [main parameter] values from PRIOR ++++++++++++++++++++++++++++++(cluster-wise, but...)
xi_prev_full = rt(n=1, df=nu0)
xi_prev = rep(xi_prev_full, J) # imagine we already have them...old days..









#####> Second, Prepare for Metropolis Hastings for [betaparam], [sig2param], [xiparam] ~ LSN(.)*MVN(_)*InvGa(_)*T(_)
########################################################################### SINGLE SINGLE SINGLE SINGLE SINGLE
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
SIG_b0param = riwish(v = qu0+2,
                     S = t(beta0_prev-beta_prev_full) %*% (beta0_prev-beta_prev_full)+
                       delta*t(beta0_prev - m0) %*% (beta0_prev - m0) + LAMB)

##> Sample beta0param from multivariate normal
beta0param = rmvn(n = 1, mu = delta/(delta+1)*m0 + 1/(delta+1)*beta_prev_full,
                  sigma = 1/(delta+1)*SIG_b0param)

##> Sample u0 from Metropolis Hastings using Gamma(3,1) as proposal
u0_proposal = rgamma(n = 1, shape = u0_proposed_shape, rate = u0_proposed_rate)
ratio_u0 = min(fink(u0_proposal, A = sig2_prev_full*rho_u1/(v0_prev/2), B = rho_u2+1)/
                 fink(u0_prev, A = sig2_prev_full*rho_u1/(v0_prev/2), B = rho_u2+1)*
                 dgamma(x = u0_prev, shape = u0_proposed_shape, rate = u0_proposed_rate)/
                 dgamma(x = u0_proposal, shape = u0_proposed_shape, rate = u0_proposed_rate),1)
if(runif(n = 1, min = 0, max = 1) < ratio_u0) {
  u0param = u0_proposal
} else {
  u0param = u0_prev
}

## Sample v0 from gamma distribution
v0param = rgamma(n=1, shape=rho_v1 + u0param/2, rate=rho_v2 + 0.5/sig2_prev_full)

SIG_b0param; beta0param; u0param; v0param  
# for sharing...#### Now we have single samples to start MH engine.


#####> From here!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
############################################################################### "Cluster-wise"!!!!!!!
##] Run your MH to obtain single posterior samples of serious main parameters: [betaparamS], [sig2param] 

for(j in 1:J) {
  
  # Sample proposals from priors
  beta_prop = rmvn(n=1, mu=beta0param, sigma=SIG_b0param) #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~(CHECK)
  sig2_prop = rinvgamma(n=1, shape=u0param, scale=v0param)
  xi_prop = rt(n=1, df=nu0)
  #sig2_prop = rinvgamma(n=1, shape=u0param*3, scale=v0param*2)
  #xi_prop = rnorm(n = 1, mean = 0, sd = nu0.st/10)
  
  # subsetting by cluster
  Sj = S[cl_membership==j]
  #logSj = logS[cl_membership==j]
  matXj = matX_tru.xx[cl_membership==j, ]
  
  
  
  ##> A. MH algorithm for: [betaparam]
  #-Note: "beta_prop" is from PRIOR proposal (new) while "beta_prev[j,]" is from the single POSTERIOR sample (old). 
  numerator=0
  denominator=0
  for (i in 1:length(Sj)){
    #for (i in 1:length(logSj)){
    numerator = numerator + 
      #dsknorm_log(logSj[i], matXj[i,], t(as.matrix(beta_p)), sig2_prev[j], xi_prev[j])
      #dsknorm_log(Sj[i], matXj[i,], t(as.matrix(beta_p)), sig2_p, xi_p)
      dlogsknorm_log(Sj[i], matXj[i,], t(as.matrix(beta_prop)), sig2_prev[j], xi_prev[j])
    
    denominator = denominator + 
      #dsknorm_log(logSj[i], matXj[i,], t(as.matrix(beta_prev[j,])), sig2_prev[j], xi_prev[j])
      #dsknorm_log(Sj[i], matXj[i,], t(as.matrix(beta_prev[j,])), sig2_prev[j], xi_prev[j])
      dlogsknorm_log(Sj[i], matXj[i,], t(as.matrix(beta_prev[j,])), sig2_prev[j], xi_prev[j])
  }
  # compute the ratio
  ratio_beta = min(exp(numerator-denominator), 1)
  
  U = runif(n = 1, min = 0, max = 1)
  if(U < ratio_beta) {
    betaparam[j, ] = beta_prop
  } 
  else {
    betaparam[j, ] = beta_prev[j, ]
  }
  
  
  ##> B. MH algorithm for: [sig2param]
  # prepare components
  numerator=0
  denominator=0
  for (i in 1:length(Sj)){
    #for (i in 1:length(logSj)){
    numerator = numerator + 
      dlogsknorm_log(Sj[i], matXj[i,], t(as.matrix(betaparam[j,])), sig2_prop, xi_prev[j]) + 
      #dsknorm_log(logSj[i], matXj[i,], t(as.matrix(betaparam[j,])), sig2_p, xi_prev[j]) + 
      log(dinvgamma(sig2_prop, shape = u0param, scale = v0param)) + 
      log(dinvgamma(sig2_prev[j], shape = u0param, scale = v0param))
    #log(dinvgamma(sig2_prev[j], shape = 3*u0param, scale = 2*v0param)) #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~(CHECK)
    
    denominator = denominator + 
      dlogsknorm_log(Sj[i], matXj[i,], t(as.matrix(betaparam[j,])), sig2_prev[j], xi_prev[j])+ 
      #dsknorm_log(logSj[i], matXj[i,], t(as.matrix(betaparam[j,])), sig2_prev[j], xi_prev[j])+ 
      log(dinvgamma(sig2_prev[j], shape = u0param, scale = v0param)) + 
      log(dinvgamma(sig2_prop, shape = u0param, scale = v0param))
    #log(dinvgamma(sig2_prop, shape = 3*u0param, scale = 2*v0param))   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~(CHECK)
  }
  # compute the ratio
  ratio_sig2 = min(exp(numerator-denominator), 1)
  
  U = runif(n = 1, min = 0, max = 1)
  if(U < ratio_sig2) {
    sig2param[j] = sig2_prop
  } 
  else {
    sig2param[j] = sig2_prev[j]
  }
  
  
  
  ##> C. MH algorithm for: [xiparam]
  # prepare components
  numerator=0
  denominator=0
  for (i in 1:length(Sj)){
    #for (i in 1:length(logSj)){
    numerator = numerator + 
      dlogsknorm_log(Sj[i], matXj[i,], t(as.matrix(betaparam[j,])), sig2param[j], xi_prop) +
      #dsknorm_log(logSj[i], matXj[i,], t(as.matrix(betaparam[j,])), sig2param[j], xi_p) +
      dt(xi_prop, df=nu0, log = TRUE) + dt(xi_prev[j], df=nu0, log = TRUE)
    #dt(xi_prop, df=nu0, log = TRUE) + dnorm(xi_prev[j], mean = 0, sd = nu0, log = TRUE) 
    #dt(xi_prop, df = nu0, log = TRUE) + dnorm(xi_prev[j], mean = 0, sd = nu0/10, log = TRUE) #~~~~~~(CHECK)
    denominator = denominator + 
      dlogsknorm_log(Sj[i], matXj[i,], t(as.matrix(betaparam[j,])), sig2param[j], xi_prev[j]) +
      #dsknorm_log(logSj[i], matXj[i,], t(as.matrix(betaparam[j,])), sig2param[j], xi_prev[j]) +
      dt(xi_prev[j], df=nu0, log = TRUE) + dt(xi_prop, df=nu0, log = TRUE)
    #dt(xi_prev[j], df=nu0, log = TRUE) + dnorm(xi_prop, mean = 0, sd = nu0, log = TRUE)
    #dt(xi_prev[j], df=nu0, log = TRUE) + dnorm(xi_prop, mean = 0, sd = nu0/10, log = TRUE)  #~~~~~~(CHECK)
  }
  # compute the ratio
  ratio_xi = min(exp(numerator-denominator), 1)
  
  U = runif(n = 1, min = 0, max = 1)
  if(U < ratio_xi) {
    xiparam[j] = xi_prop
  } 
  else {
    xiparam[j] = xi_prev[j]
  }
}
betaparam; sig2param; xiparam
beta_prev; sig2_prev; xi_prev



# NOW..................

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
# We NEED a Cluster-wise set of....... Many Smaples ............................
# - outcome parameters: [betaparamS], [sig2param], [xiparam]
# - covariate parameters: [piparamS], [X_barS], [lambda2paramS]
# 
#..............................................................................#
#..............................................................................#
#..............................................................................#
#................................ now be ready ................................# 
#..............................................................................#
#..............................................................................#
#..............................................................................#
################################################################################
### Step03> Gibbs Sampler --------- cl_membership and param update ---with J= ?
################################################################################
set.seed(1)
total_iter = 60
# r_convergence = 1000   # first run Gibbs sampler to determine r value for convergence
r_convergence = 10

# loglikelihood = numeric(total_iter)                  # for monitor convergence
loglikelihood = matrix(0, nrow = n.train, ncol = total_iter)
# [note] why matrix? each col gives each iteration result


###] Prepare a "pool" to store the ultimate hyperparameters, parameters ++++++++++++++++++++++++++++++++++++++++
#> w/o cluster (hyperparameter)
list_beta0param = list() # For beta0 hyperparameter
list_SIG_b0param = list() # For SIG_b0 hyperparameter
list_u0param = list() # for u0 hyperparameter
list_v0param = list() # for v0 hyperparameter


#> w/o cluster (outcome)
list_betaparamfull = list()  # for S - full data
list_sig2paramfull = list()  # for S - full data
list_xiparamfull = list()    # for S - full data

#> w/o cluster (covariate)
# N/A                                                                       #~xxxxxxxxxxxxxxxxxxxxx Q.DAVID
# N/A                                                                       #~xxxxxxxxxxxxxxxxxxxxx Q.DAVID
#-----------------------------------------------------------


#> cluster-wise (outcome)
list_betaparam = list()   #for S
list_sig2param = list()    #for S
list_xiparam = list()      #for S

#> cluster-wise (covariate)
list_piparam = list()   #for Z
list_lambda2 = list()    #for X
#list_X_barS = list()     #for X...always same?



###] Store the current versions of single "cluster-wise" parameters as the "previous" ++++++++++++++++++++++++++++
#> hyperparameter
beta0_prev = beta0param
SIG_b0_prev = SIG_b0param
u0_prev = u0param
v0_prev = v0param

#> outcome
beta_prev = betaparam
sig2_prev = sig2param
xi_prev = xiparam

#> covariate
#piparam #always same?
#lambda2param # do we need???? 
#X_bar # always same?





###] Miscellaneous ???

max_iter = 10   # let's say.. investigation( sampling rejection ) limit...?
counts_mat <- matrix(0, nrow = J*total_iter, ncol = 3) # count..How many rejected???

n_j = table(cl_membership)  # just temporary...for cluster-wise sampling...

beta_p_accepted <- 0
sig2_p_accepted <- 0
xi_p_accepted <- 0
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

########################################################## START ###############################################################

for (r in 1:total_iter) {
  
  ###[1] Updating Posterior -----------------------------------------------------------------------------------------------
  
  #>>> First, sample ***HYPERPARAMETERS*** from hyperPosterior based on full data {w/o clusters}
  
  ## hyperparam for [betaparamS]
  SIG_b0param = riwish(v=qu0+2, 
                       S=t(beta0_prev-beta_prev_full) %*% (beta0_prev-beta_prev_full) +
                         delta*t(beta0_prev - m0) %*% (beta0_prev - m0) + LAMB) #~~~~~~~~~~~~~~~~~~(CHECK)
  
  beta0param = rmvn(n=1, mu=delta/(delta+1)*m0 + 1/(delta+1)*beta_prev_full, sigma=1/(delta+1)*SIG_b0param) 
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~(CHECK)
  
  ## hyperparam for [sig2param]  
  #> - Sample u0 from MH using Gamma(3.75, 2) as proposal 
  u0_prop = rgamma(n=1, shape=u0_proposed_shape, rate=u0_proposed_rate) # Fink approx
  ratio_u0 = min(fink(u0_prop, A = sig2_prev_full*rho_u1/(v0_prev/2), B = rho_u2+1)/
                   fink(u0_prev, A = sig2_prev_full*rho_u1/(v0_prev/2), B = rho_u2+1)*
                   dgamma(x = u0_prev, shape = u0_proposed_shape, rate = u0_proposed_rate)/
                   dgamma(x = u0_prop, shape = u0_proposed_shape, rate = u0_proposed_rate), 1)
  if(runif(n = 1, min = 0, max = 1) < ratio_u0) {
    u0param = u0_prop
  } else {
    u0param = u0_prev
  }
  
  #> - Sample v0 from Gamma distribution 
  v0param = rgamma(n=1, shape=rho_v1 + u0param/2, rate = rho_v2+0.5/sig2_prev_full)
  
  
  
  
  
  #------------------------------------------------------------------------------------------------------------------
  #------------------------------------------------------------------------------------------------------------------
  
  #>>> Next I, Sample ****[Main parameters]**** 
  
  #> - Proposal Sampling using proposals OF PRIORS {w/o clusters} 
  beta_prop_full = rmvn(n=1, mu=beta0param, sigma=SIG_b0param*varinf) #~~~~~~~~~~~~~~~~~~~~~~~~(CHECK)
  sig2_prop_full = rinvgamma(n=1, shape=u0param, scale=v0param) 
  xi_prop_full = rt(n=1, df=nu0)
  
  #> - Sifting a single [betaparam] based on FULL data 
  numerator=0
  denominator=0
  # plugging prior sample first...
  for (i in 1:n.train){
    numerator = numerator + 
      dlogsknorm_log(s=S[i], x=matX_tru.xx[i,], beta=t(as.matrix(beta_prop_full)), sig2=sig2_prev_full, xi=xi_prev_full)
    
    denominator = denominator + 
      dlogsknorm_log(s=S[i], x=matX_tru.xx[i,], beta=t(as.matrix(beta_prev_full)), sig2=sig2_prev_full, xi=xi_prev_full)
  }
  ratio_beta = min(exp(numerator-denominator), 1)
  
  U = runif(n = 1, min = 0, max = 1)
  if(U < ratio_beta) {
    betaparam_full = beta_prop_full
  } 
  else {
    betaparam_full = beta_prev_full
  }
  
  
  #> - Sifting a single [sig2param] based on FULL data 
  numerator=0
  denominator=0
  # plugging prior sample first...
  for (i in 1:n.train){
    numerator = numerator + 
      dlogsknorm_log(s=S[i], x=matX_tru.xx[i,], beta=t(as.matrix(betaparam_full)), sig2=sig2_prop_full, xi=xi_prev_full)
    
    denominator = denominator + 
      dlogsknorm_log(s=S[i], x=matX_tru.xx[i,], beta=t(as.matrix(betaparam_full)), sig2=sig2_prev_full, xi=xi_prev_full)
  }
  ratio_sig2 = min(exp(numerator-denominator), 1)
  
  U = runif(n = 1, min = 0, max = 1)
  if(U < ratio_sig2) {
    sig2param_full = sig2_prop_full
  } 
  else {
    sig2param_full = sig2_prev_full
  }
  
  
  #> - Sifting a single [xiparam] based on FULL data 
  numerator=0
  denominator=0
  # plugging prior sample first...
  for (i in 1:n.train){
    numerator = numerator + 
      dlogsknorm_log(S[i], matX_tru.xx[i,], t(as.matrix(betaparam_full)), sig2param_full, xi_prop_full)
    
    denominator = denominator + 
      dlogsknorm_log(S[i], matX_tru.xx[i,], t(as.matrix(betaparam_full)), sig2param_full, xi_prev_full)
  }
  # compute the ratio
  ratio_xi = min(exp(numerator-denominator), 1)
  
  U = runif(n = 1, min = 0, max = 1)
  if(U < ratio_xi) {
    xiparam_full = xi_prop_full
  } 
  else {
    xiparam_full = xi_prev_full
  }
  
  
  
  #--------------------------------------------------------------------------------------------------------#
  #--------------------------------------------------------------------------------------------------------#
  #-----------------------------------Now it's time for cluster-wise---------------------------------------#
  #--------------------------------------------------------------------------------------------------------#
  #--------------------------------------------------------------------------------------------------------#
  for(j in 1:J) {
    
    #>>> First.., Sample **[covariate parameters_j]**
    X_j = X_tru.xx[cl_membership==j] #covariate
    Z_j = Z[cl_membership==j] #covariate
    n_j = length(X_j)
    
    #> - proposal Sampling [lambda2paramF_j] using proposals OF POSTERIOR {with clusters} 
    lambda2param[j] = rinvgamma( n=1, 
                                 shape=(c0+n_j)/2, 
                                 scale=(d0+sum( (X_j - mean(X_j))^2 ))/2 )  
    
    #> - proposal Sampling [piparamF_j] using proposals OF POSTERIOR {with clusters} 
    piparam[j] = rbeta( n=1, shape1=g0+sum(Z_j), shape2=h0+n_j-sum(Z_j) )
    
    
    #------------------------------------------------------------------------------------------------------------------
    #------------------------------------------------------------------------------------------------------------------
    #>>> Next,... Sample ****[Main parameters_j]****
    
    #>>> subsetting by cluster
    S_j = S[cl_membership==j]           #outcome
    matX_j = matX_tru.xx[cl_membership==j, ] #covariate
    
    
    #>>> - Sifting [betaparam]_j
    ratio_beta = 0
    count_beta = 0
    U = 1
    print("beta sampling")
    while((U > ratio_beta) & (count_beta < max_iter)) {
      # if new samples keep being rejected, ?????
      beta_prop = rmvn(n=1, mu=beta0param, sigma=SIG_b0param*varinf) #~~~~~~~~~~~~~~~~~~~~~~~~(CHECK)
      numerator=0
      denominator=0
      for (i in 1:length(S_j)){
        numerator = numerator + 
          dlogsknorm_log(S_j[i], matX_j[i,], t(as.matrix(beta_prop)), sig2_prev[j], xi_prev[j])
        
        denominator = denominator + 
          dlogsknorm_log(S_j[i], matX_j[i,], t(as.matrix(beta_prev[j,])), sig2_prev[j], xi_prev[j])
      }
      ratio_beta = min(exp(numerator-denominator), 1)
      U = runif(n = 1, min = 0, max = 1)
      count_beta = count_beta+1
    }
    print(count_beta)
    
    counts_mat[J*(r-1)+j, 1] <- count_beta
    
    if(count_beta < max_iter) {
      betaparam[j, ] = beta_prop
      beta_p_accepted <- beta_p_accepted+1
    } else {
      betaparam[j,] = beta_prev[j,]
    }
    
    
    #>>> - Sifting [sig2param]_j
    ratio_sig2 = 0
    count_sig2 = 0
    U = 1
    print("sig2 sampling")
    while((U > ratio_sig2) & (count_sig2 < max_iter)) {
      sig2_prop = rinvgamma(n=1, shape=u0_new, scale=v0_new) 
      numerator=0
      denominator=0
      for (i in 1:length(S_j)){
        numerator = numerator + 
          dlogsknorm_log(S_j[i], matX_j[i,], t(as.matrix(betaparam[j,])), sig2_prop, xi_prev[j]) +
          log(dinvgamma(sig2_prop, shape=u0param, scale=v0param)) + 
          log(dinvgamma(sig2_prev[j], shape=u0_new, scale=v0_new))
        
        denominator = denominator + 
          dlogsknorm_log(S_j[i], matX_j[i,], t(as.matrix(betaparam[j,])), sig2_prev[j], xi_prev[j]) +
          log(dinvgamma(sig2_prev[j], shape=u0param, scale=v0param)) + 
          log(dinvgamma(sig2_prop, shape=u0_new, scale=v0_new))
      }
      ratio_sig2 = min(exp(numerator-denominator), 1)
      U = runif(n = 1, min = 0, max = 1)
      count_sig2 = count_sig2+1
    }
    print(count_sig2)
    counts_mat[J*(r-1)+j, 2] <- count_sig2
    if(U < ratio_sig2) {
      sig2param[j] = sig2_prop
      sig2_p_accepted <- sig2_p_accepted+1
    }
    else {
      sig2param[j] = sig2_prev[j]
    }
    
    
    #>>> - Sifting [xiparam]_j
    ratio_xi = 0
    count_xi = 0
    U = 1
    print("xi sampling")
    while((U > ratio_xi) & (count_xi < max_iter)) {
      xi_prop = rt(n=1, df=nu0)
      numerator=0
      denominator=0
      for (i in 1:length(S_j)){
        numerator = numerator + 
          dlogsknorm_log(S_j[i], matX_j[i,], t(as.matrix(betaparam[j,])), sig2param[j], xi_prop)
        
        denominator = denominator + 
          dlogsknorm_log(S_j[i], matX_j[i,], t(as.matrix(betaparam[j,])), sig2param[j], xi_prev[j])
      }
      # compute the ratio
      ratio_xi = min(exp(numerator-denominator), 1)
      U = runif(n = 1, min = 0, max = 1)
      count_xi = count_xi + 1
    }
    print(count_xi)
    counts_mat[J*(r-1)+j, 3] <- count_xi
    if(U < ratio_xi) {
      xiparam[j] = xi_prop
      xi_p_accepted <- xi_p_accepted+1
    }
    else {
      xiparam[j] = xi_prev[j]
    }
    
  }
  
  
  #################################################################################################################
  
  
  #### For NEXT!!!
  ###> Set previous set of parameter to the current for the next iteration
  ## - Manner! main param for next iteration
  lambda2param_prev = lambda2param
  piparam_prev = piparam
  
  beta_prev_full = betaparam_full
  beta_prev = betaparam
  
  sig2_prev_full = sig2param_full 
  sig2_prev = sig2param
  
  xi_prev_full = xiparam_full
  xi_prev = xiparam
  
  ## - Manner! hyperparam for next iteration
  beta0_prev = beta0param
  
  SIG_b0_prev = SIG_b0param
  
  u0_prev = u0param
  
  v0_prev = v0param
  
  
  
  
  
  
  ###[2] Calculating Loglikelihood to monitor the chain convergence ---------------------------------------------
  #>>> First
  loglike_r = numeric(n.train)
  
  for(i in 1:n.train) {
    x_i = matX_tru.xx[i, ]      # first row..{1,x,z}
    j = cl_membership[i] # membership ID
    
    loglike_r[i] = dlogsknorm_log(s=S[i], x=x_i, beta=betaparam[j,], sig2=sig2param[j], xi=xiparam[j]) + 
      dnorm(x=x_i[2], mean=X_bar[j], sd=sqrt(lambda2param[j]), log=TRUE) +
      dbinom(x=x_i[3], size=1, prob=piparam[j], log=TRUE) 
  }  
  loglikelihood[,r] = loglike_r
  
  if(r > r_convergence) {
    
    # W/O CLUSTER
    list_beta0param[[r-r_convergence]] = beta0param
    list_SIG_b0param[[r-r_convergence]] = SIG_b0param
    list_u0param[[r-r_convergence]] = u0param # summary(unlist(list_u0paramS)) # for u0_new
    list_v0param[[r-r_convergence]] = v0param # summary(unlist(list_v0paramS)) # for v0_new
    
    list_betaparamfull[[r-r_convergence]] = betaparam_full
    list_sig2paramfull[[r-r_convergence]] = sig2param_full
    list_xiparamfull[[r-r_convergence]] = xiparam_full
    
    
    # Cluster-wise
    list_piparam[[r-r_convergence]] = piparam     
    list_lambda2[[r-r_convergence]] = lambda2param  
    
    list_betaparam[[r-r_convergence]] = betaparam  
    list_sig2param[[r-r_convergence]] = sig2param
    list_xiparam[[r-r_convergence]] = xiparam
  }
  print(paste("r=",r))
  trace(what = "print", where = getNamespace("base"), exit = flush.console, print = FALSE)
}

################################################################################
####### End of Gibbs Sampler ###################################################
################################################################################

##> Check your sampling performance
beta_p_accepted/(total_iter*J)
sig2_p_accepted/(total_iter*J)
xi_p_accepted/(total_iter*J)

par(mfrow = c(1,1))
plot(colSums(loglikelihood), type="l")

# #-------------------------------------------------------------------------------
# ###> Investigation for "sig2"
# sig2.vec_cl1 = sapply(list_sig2param, function(x) x[1])
# sig2.vec_cl2 = sapply(list_sig2param, function(x) x[2])
# sig2.vec_cl3 = sapply(list_sig2param, function(x) x[3])
# sig2.vec_cl4 = sapply(list_sig2param, function(x) x[4])
# sig2.vec_cl5 = sapply(list_sig2param, function(x) x[5])
# 
# sig2.df <- data.frame(
#   value = c(sig2.vec_cl1, sig2.vec_cl2, sig2.vec_cl3, sig2.vec_cl5, sig2.vec_cl5),
#   group = factor(rep(1:5, each = total_iter-r_convergence))
# )
# 
# #----------------------------------------- R they InvGamma?
# ggplot(sig2.df, aes(x = value)) +
#   geom_density(size=1) +
#   facet_wrap(~ group, scales = "free") +
#   labs(title = "Density for Each sig2.vec", x = "sig2", y = "Density") + xlim(0, 50)
# theme_minimal()
# tapply(sig2.df$value, sig2.df$group, mean)
# 
# # 95% credible interval
# compute_CI <- function(samples, alpha = 0.05) {
#   quantile(samples, probs = c(alpha / 2, 1 - alpha / 2))
# }
# 
# CI <- sig2.df %>%
#   group_by(group) %>%
#   summarise(
#     Lower = compute_CI(value, alpha = 0.05)[1],
#     Upper = compute_CI(value, alpha = 0.05)[2]
#   )
# options(scipen = 999)
# CI
# 
# 
# #-------------------------------------------------------------------------------
# ###> Investigation for "xi"
# xi.vec_cl1 = sapply(list_xiparam, function(x) x[1])
# xi.vec_cl2 = sapply(list_xiparam, function(x) x[2])
# xi.vec_cl3 = sapply(list_xiparam, function(x) x[3])
# xi.vec_cl4 = sapply(list_xiparam, function(x) x[4])
# xi.vec_cl5 = sapply(list_xiparam, function(x) x[5])
# 
# xi.df <- data.frame(
#   value = c(xi.vec_cl1, xi.vec_cl2, xi.vec_cl3, xi.vec_cl5, xi.vec_cl5),
#   group = factor(rep(1:5, each = total_iter-r_convergence))
# )
# #----------------------------------------- R they T?
# ggplot(xi.df, aes(x = value)) +
#   geom_density(size=1) +
#   facet_wrap(~ group, scales = "free") +
#   labs(title = "Density for Each xi.vec", x = "xi", y = "Density") + xlim(-3, 3)
# theme_minimal()
# tapply(xi.df$value, xi.df$group, mean)
# 
# # 95% credible interval
# compute_CI <- function(samples, alpha = 0.05) {
#   quantile(samples, probs = c(alpha / 2, 1 - alpha / 2))
# }
# 
# CI <- xi.df %>%
#   group_by(group) %>%
#   summarise(
#     Lower = compute_CI(value, alpha = 0.05)[1],
#     Upper = compute_CI(value, alpha = 0.05)[2]
#   )
# options(scipen = 999)
# CI
# 
# 
# 
# #-------------------------------------------------------------------------------
# ###> Investigation for "beta0", "beta1", "beta2"
# # Sum all the matrices
# sum_matrix <- Reduce("+", list_betaparam)
# 
# # Number of matrices
# nb <- length(list_betaparam)
# 
# # Calculate the average matrix of betaS
# avg_beta.mat <- sum_matrix / nb
# avg_beta.mat

#-------------------------------------------------------------------------------

############################### Predictions ####################################

# Expected Value of LSN distribution (Wang 2019)
# 2*exp(sum(x_i*beta_j[j,]) + sig2_j[j]/2)*(1-pnorm(-xi_j[j]*sqrt(sig2_j[j])/sqrt(xi_j[j]^2+1)))


###> Training data
expval.train = matrix(0, nrow = n.train, ncol = total_iter-r_convergence)

for(r in 1:(total_iter-r_convergence)) {
  betaparamclean = list_betaparam[[r]]
  sig2paramclean = list_sig2param[[r]]
  xiparamclean = list_xiparam[[r]]
  
  for(i in 1:n.train) {
    j = cl_membership[i]
    #cleanx = rnorm(n = 1, mean = kappaparamclean[j,1] + kappaparamclean[j,2]*Z[i], sd = sqrt(lambda2paramclean[j]))
    cleanx = X_tru.xx[i]
    #expval.train[i,r] <- 2*exp( betaparamclean[j,1] + betaparamclean[j,2]*cleanx + betaparamclean[j,3]*Z[i] + 
    #                              sig2paramclean[j]/2)*( 1-pnorm(-xiparamclean[j]*
    #                                                               sqrt(sig2paramclean[j])/sqrt(xiparamclean[j]^2+1)))
    expval.train[i,r] <- 10*exp( betaparamclean[j,1] + betaparamclean[j,2]*cleanx + betaparamclean[j,3]*Z[i] - 
                                   sig2paramclean[j]/2)*( 1-pnorm(-xiparamclean[j]*
                                                                    sqrt(sig2paramclean[j])/sqrt(xiparamclean[j]^2+1)))
  }
}

expS.train <- apply(X=expval.train, MARGIN=1, FUN=mean)

summary(S)
summary(expS.train)

format( aggregate(log(S), by=list(cl_membership), FUN=summary), scientific=FALSE)
format( aggregate(log(expS.train), by=list(cl_membership), FUN=summary), scientific=FALSE)

lppd_EX_S.S <- mean( colSums(loglikelihood) ); lppd_EX_S.S # -25705.38

# Clusterwise predictive density plots
par(mfrow = c(2,3))

for(j in 1:J) {
  hist(log(S)[cl_membership==j], freq=FALSE, breaks=50, col="white", xlim=c(5,15), ylim=c(0,1.6),
       xlab = paste("Cluster", j))
  lines(density(log(expS.train)[cl_membership==j]), col="black", ylim=c(0,1.6), lwd=2)
}


###> Testing data
expval.test = matrix(0, nrow = n.test, ncol = total_iter-r_convergence)

for(r in 1:(total_iter-r_convergence)) {
  betaparamclean = list_betaparam[[r]]
  sig2paramclean = list_sig2param[[r]]
  xiparamclean = list_xiparam[[r]]

  for(i in 1:n.test) {
    j = cl_membership[i]
    #cleanx = rnorm(n = 1, mean = kappaparamclean[j,1] + kappaparamclean[j,2]*Z[i], sd = sqrt(lambda2paramclean[j]))
    cleanx = X_tru.xx.test[i]
    expval.test[i,r] <- 10*exp( betaparamclean[j,1] + betaparamclean[j,2]*cleanx + betaparamclean[j,3]*Z[i] -
                                  sig2paramclean[j]/2)*( 1-pnorm(-xiparamclean[j]*
                                                                   sqrt(sig2paramclean[j])/sqrt(xiparamclean[j]^2+1)))
  }
}
expS.test <- apply(X=expval.test, MARGIN=1, FUN=mean)

summary(S.test)
summary(expS.test)

format( aggregate(log(S.test), by=list(cl_membership.test), FUN=summary), scientific=FALSE)
format( aggregate(log(expS.test), by=list(cl_membership.test), FUN=summary), scientific=FALSE)


# Clusterwise predictive density plots
par(mfrow = c(2,3))

for(j in 1:J) {
  hist(log(S.test)[cl_membership.test==j], freq=FALSE, breaks=50, col="white", xlim=c(5,15), ylim=c(0,1.6),
       xlab = paste("Cluster", j))
  lines(density(log(expS.test)[cl_membership.test==j]), col="red", ylim=c(0,1.6), lwd=2)
}

# NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW
# NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW
# NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW
# NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW
# NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW

################################################################################
###################### S ~ hierarchical LN with Gustafson ######################
################################################################################
################################################################################
### Step01> Define prior, hyperprior model and hyperparameter values
################################################################################

# ------------- Outcome ----- # ::: for S ~ LSN( X*betaparam, sig2param, xiparam )

# PRIOR: "betaparamF" ~ MVN(beta0_j, SIG2_b0_j)
# Hyperprior                beta0_j ~ MVN(m0, 1/delta*SIG_b0_j)
# Hyperprior                          SIG2_b0_j ~ IW( qu0, LAMB )

# ::: for "beta0_j" ~ MVN( m0, 1/delta*SIG_b0_j )
#> Gaussian regression for initialize OUTCOME model parameter m0, SIG_b0
fit1.st <- glm( S ~ Xstar + factor(Z) , family=poisson(link = "log") )
# fit1 <- lm( logS ~ Xstar + factor(Z) )
summary(fit1.st)

m0.st = coef(fit1.st)    # 1x3 initial reg_coeff vector (sort of "means"): Regression result as "mean"
SIG_b0.st = vcov(fit1.st)   # 3x3 initial cov matrix of reg_coeff
p = length(m0.st)
options(scipen = 999)
SIG_b0.st
SIG_b0inv.st = solve(a=SIG_b0.st) # inverse of cov matrix of reg_coeff for later use (for posterior on reg_coeff)!!

#****** for VAR Inflation factor
delta = 0.01             # equal importance on m0 and beta ????
varinf = n.train/100 # rule of thumb divide by 5


#> SIG_b0
#- it requires -------
qu0.st = p+2
# LAMB ?????????????????????????????????????????????
LAMB.st = SIG_b0.st





###> PRIOR: "sig2param" ~ IG( u0, v0 )
# Hyperprior               u0 ~ Fink( rho_u1, rho_u2 )
# Hyperprior                   v0 ~ Ga( rho_bv1, rho_v2 )

# Hyperparameter u0
rho_u1.st = 1/8
rho_u2.st = 3/2 

# Hyperparameter v0
rho_v1.st = 8 
rho_v2.st = 1 

# Based on the result? mean(unlist(list_u0param)), mean(unlist(list_v0param))
# E[invGa] = scale/(shape-1) = 50
u0_new.st = 1.7591 #3 #2 #1.1 # from 3/7  # from # summary(unlist(list_u0param)) # for u0_new                            
v0_new.st = 7.113  #3 #5 #8.5 # from 3/7  # from # summary(unlist(list_v0param)) # for v0_new

## [CHECK later] ---------------------------------------------------------------
# How can we design proposal for sig2 ~ InvGa(.) ?
par(mfrow=c(1,1))
plot( density(unlist(list_sig2param.st)) )                                   # sig2 we have...before
curve( dinvgamma(x, shape=mean(unlist(list_u0param.st)), 
                 scale=mean(unlist(list_v0param.st))), add=TRUE, col="red") # sig2 we obtain....after
curve( dinvgamma(x, shape=u0_new.st, 
                 scale=v0_new.st), add=TRUE, col="blue")                    # sig2 proposal adjustment

table(unlist(list_sig2param.st))
sort(table(unlist(list_sig2param.st)),DESC=TRUE)                             # to see if they are appropriate...
#-------------------------------------------------------------------------------

###> PRIOR: "xiparam"~t(loc, nu0, sca)
#loc = 0   # ? Location hyperparameter for T distribution
#nu0 = 1/2 # df for T distribution on pigtail parameter
#sca = 5   # ? Scale hyperparameter for T distribution 
nu0.st = n.train-1
#- df refers to the number of independent observations (N-1):


# ------------------------------------------------------------------------------ #











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
SIG_k0 <- vcov(fit2.st)*varinf
SIG_inv_k0 <- solve(SIG_k0)


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

# tau2 = (1-zeta_scale)*lambda2dirty
# zeta_scale = c(0.01,0.05,0.1,0.2,0.5)
# zeta_scale = c(0.60,0.70,0.80,0.90,0.95)

zeta_scale = 0.95



################################################################################
### Step02> initialize major parametersssss by clusters using your "POSTERIOR"
#                                                 Based on natural clustering
################################################################################
J = length(table(cl_membership))
set.seed(1)

# ---------------- for Exposure model + Measurement Model --------------
# : Exposure Model (Z only)
##[A] for Z ~ Bin( n=1, "prob_j" ) .. starting with 3 clusters for beta POSTERIOR for Z

piparam.st = numeric(J)
for (j in 1:J) {
  Zj=Z[cl_membership==j]
  nj=length(Zj)
  piparam.st[j] = rbeta( n=1, shape1=g0.st+sum(Zj), shape2=h0.st+nj-sum(Zj) ) 
}

piparam.st #% pi parameter sampled from posterior

#% [Check!] empirical pi ? using "aggregate(.)": investigation by splitting the data into subset
pi_empirical = aggregate(x=Z, by=list(cl_membership), FUN=mean, na.rm= T)$x; pi_empirical
piparam.st - pi_empirical 




# : Exposure model Xstar|Z
##[B] for Xstar|Z ~ N( "kappa0_j"+"kappa1_j"*Z, "lambda_j" ) .. starting with 3 clusters for Normal/IG POSTERIOR
kappaparam.st = matrix(nrow = J, ncol = length(kappa_v0)); kappaparam.st
lambda2param.st = numeric(J); lambda2param.st
for (j in 1:J) {
  Xstar_j=Xstar[cl_membership==j]
  Zj = Z[cl_membership==j]
  nj.st=length(Xstar_j)
  lambda2param.st[j]=rinvgamma( n=1, shape=(nj.st+c0.st)/2, scale=( sum((Xstar_j - kappa_v0[1]-kappa_v0[2]*Zj)^2) + d0.st )/2 )
  KK1j = cbind(1, Zj)
  KK2j = matrix(c(sum(Xstar_j), sum(Xstar_j*Zj)), nrow =2, ncol = 1) # 2x1 vector
  
  SIG_inv_j = solve(SIG_inv_k0+(t(KK1j) %*% KK1j)) # solve( ..[matrix].. ) = [matrix]^(-1)
  kappaparam.st[j,] <- rmvn(n=1, mu = SIG_inv_j %*% (SIG_inv_k0 %*% kappa_v0 + KK2j), sigma = lambda2param[j]*SIG_inv_j)      
}


kappaparam.st
lambda2param.st

#% [Check!] empirical kappa, lambda ? using "aggregate(.)": investigation by splitting the data into subset
kappa_empirical=rbind(coef(lm(Xstar[cl_membership==1] ~ Z[cl_membership==1])),
                      coef(lm(Xstar[cl_membership==2] ~ Z[cl_membership==2])),
                      coef(lm(Xstar[cl_membership==3] ~ Z[cl_membership==3]))); kappa_empirical

lambda_empirical.st = aggregate(x=as.numeric(Xstar), by=list(cl_membership), FUN=var, na.rm=T)$x; lambda_empirical.st




# : Measurement Model initial parameter tau^2
##[D] for Xstar|X ~ N( X, "tau^2") 
# tau2param <- rep(rinvgamma(n = 1, shape = u0, scale = u0*tau2_0),J)
tau2param = (1-zeta_scale)*lambda2param.st


# ---------------------------- for Outcome -------------------------------------

###> Outcome parameters:  betaparam, sig2param, xiparam

# No conjugacy...so prepare MM sampling

# This is where new parameters go.
betaparam.st = matrix(data=NA, nrow=J, ncol=p) 
sig2param.st = numeric(J)
sig2paramfull.st = 0
xiparam.st = numeric(J)



### (A) beta0, SIG_b0 For [betaparam.st]
# 1) beta0 ~ MVN(m0, 1/delta*SIG_b0)  
# 2) SIG_b0 ~ ~ IW( qu0, LAMB )
# initial hyperparameter values from hyperprior +++++++++++++++++++++++(no need to care cluster)
# Set in place the old parameters..sample from prior...
# SIG_b0_prev_J = riwish(v = qu0, S = SIG_b0*varinf) # take one sample for SIG_b0, this will be initial parameter for all J clusters
# SIG_b0_prev = array(rep(SIG_b0_prev_J,J),dim=c(p,p,J))
# beta0_prev_J = rmvn(n = 1, mu = m0, sigma = 1/delta*SIG_b0_prev_J)
# beta0_prev = matrix(rep(beta0_prev_J,J),nrow=J,ncol=p,byrow = TRUE)
set.seed(1)
SIG_b0_prev.st = riwish(v=qu0.st, S=LAMB.st) 
beta0_prev.st = rmvn(n=1, mu=m0.st, sigma = 1/delta*SIG_b0_prev.st)



### (B) u0.st, v0.st For [sig2param.st]
# 1) u0.st ~ Fink(rho_u1.st, rho_u2.st)
# 2) v0.st ~ Ga(rho_v1.st, rho_v2.st)
# Fink function
fink = function(u0, A, B) {
  A^(u0/2-1)/gamma(u0/2)^B
}

u0_proposed_shape.st = 3 # 3
u0_proposed_rate.st = 2 # 1

###> initial hyperparameter values from hyperprior +++++++++++++++++++++++(no need to care cluster)
u0_prev.st = u0_proposed_shape.st/u0_proposed_rate.st     # instead of sampling from Fink function???
v0_prev.st = rgamma(n=1, shape=rho_v1.st, rate=rho_v2.st)

#jpeg(file="plot.test.jpeg", width=1000, height=800)
#-------------------------------------------------------------------------------
# See? Fink function is close to Ga(3.75,2) ?
curve(fink(x,rho_u1.st,rho_u2.st)/integrate(function(x) fink(x, rho_u1.st, rho_u2.st),0,Inf)$value, 
      from = 0.01, to = 10, ylim = c(0,0.5))
curve(dgamma(x, shape = u0_proposed_shape.st, rate = u0_proposed_rate.st), add=TRUE, col="blue")
#-------------------------------------------------------------------------------


### (C) nu0 For [xiparam] ~ t(loc, nu0, sca)
# nu0 ~ Exponetial(.) ??? or constant? N-1 : dt(x, df, ncp, log = FALSE) : the central t distribution
curve( dt(x, df=nu0.st) )



#-------------------------------------------------------------------------------
##### Now....ready for Single SAMPLE FROM the PRIOR
beta_prev_full.st = rmvn(n=1, mu=beta0_prev.st, sigma=SIG_b0_prev.st)
beta_prev.st = matrix( rep(beta_prev_full.st, J), nrow=J, byrow= TRUE )
sig2_prev_full.st = rinvgamma(n=1, shape=u0_prev.st, scale=v0_prev.st)
sig2_prev.st = rep(sig2_prev_full.st, J)
xi_prev_full.st = rt(n=1, df=nu0.st)
xi_prev.st = rep(xi_prev_full.st, J) # imagine we already have them...old days..


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
SIG_b0param.st = riwish(v = qu0.st+2,
                        S = t(beta0_prev.st-beta_prev_full.st) %*% (beta0_prev.st-beta_prev_full.st)+
                          delta*t(beta0_prev.st - m0.st) %*% (beta0_prev.st - m0.st) + LAMB.st)

##> Sample beta0param from multivariate normal
beta0param.st = rmvn(n = 1, mu = delta/(delta+1)*m0.st + 1/(delta+1)*beta_prev_full.st,
                     sigma = 1/(delta+1)*SIG_b0param.st)

##> Sample u0 from Metropolis Hastings using Gamma(3,1) as proposal
u0_proposal.st = rgamma(n = 1, shape = u0_proposed_shape.st, rate = u0_proposed_rate.st)
ratio_u0.st = min(fink(u0_proposal.st, A = sig2_prev_full.st*rho_u1.st/(v0_prev.st/2), B = rho_u2.st+1)/
                    fink(u0_prev.st, A = sig2_prev_full.st*rho_u1.st/(v0_prev.st/2), B = rho_u2.st+1)*
                    dgamma(x = u0_prev.st, shape = u0_proposed_shape.st, rate = u0_proposed_rate.st)/
                    dgamma(x = u0_proposal.st, shape = u0_proposed_shape.st, rate = u0_proposed_rate.st),1)
if(runif(n = 1, min = 0, max = 1) < ratio_u0.st) {
  u0param.st = u0_proposal.st
} else {
  u0param.st = u0_prev.st
}

## Sample v0 from gamma distribution
v0param.st = rgamma(n=1, shape=rho_v1.st + u0param.st/2, rate=rho_v2.st + 0.5/sig2_prev_full.st)

SIG_b0param.st; beta0param.st; u0param.st; v0param.st  
# for sharing...#### Now we have single samples to start MH engine.


#####> From here!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
############################################################################### "Cluster-wise"!!!!!!!
##] Run your MH to obtain single posterior samples of serious main parameters: [betaparamS], [sig2param] 

for(j in 1:J) {
  
  # Sample proposals from priors
  beta_prop.st = rmvn(n=1, mu=beta0param.st, sigma=SIG_b0param.st) #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~(CHECK)
  sig2_prop.st = rinvgamma(n=1, shape=u0param.st, scale=v0param.st)
  xi_prop.st = rt(n=1, df=nu0.st)
  #sig2_prop = rinvgamma(n=1, shape=u0param*3, scale=v0param*2)
  #xi_prop = rnorm(n = 1, mean = 0, sd = nu0.st/10)
  
  # subsetting by cluster
  Sj = S[cl_membership==j]
  #logSj = logS[cl_membership==j]
  matXj = matX_tru.xx[cl_membership==j, ]
  
  
  
  ##> A. MH algorithm for: [betaparam]
  #-Note: "beta_prop" is from PRIOR proposal (new) while "beta_prev[j,]" is from the single POSTERIOR sample (old). 
  numerator=0
  denominator=0
  for (i in 1:length(Sj)){
    #for (i in 1:length(logSj)){
    numerator = numerator + 
      #dsknorm_log(logSj[i], matXj[i,], t(as.matrix(beta_p)), sig2_prev[j], xi_prev[j])
      #dsknorm_log(Sj[i], matXj[i,], t(as.matrix(beta_p)), sig2_p, xi_p)
      dlogsknorm_log(Sj[i], matXj[i,], t(as.matrix(beta_prop.st)), sig2_prev.st[j], xi_prev.st[j])
    
    denominator = denominator + 
      #dsknorm_log(logSj[i], matXj[i,], t(as.matrix(beta_prev[j,])), sig2_prev[j], xi_prev[j])
      #dsknorm_log(Sj[i], matXj[i,], t(as.matrix(beta_prev[j,])), sig2_prev[j], xi_prev[j])
      dlogsknorm_log(Sj[i], matXj[i,], t(as.matrix(beta_prev.st[j,])), sig2_prev.st[j], xi_prev.st[j])
  }
  # compute the ratio
  ratio_beta.st = min(exp(numerator-denominator), 1)
  
  U = runif(n = 1, min = 0, max = 1)
  if(U < ratio_beta.st) {
    betaparam.st[j, ] = beta_prop.st
  } 
  else {
    betaparam.st[j, ] = beta_prev.st[j, ]
  }
  
  
  ##> B. MH algorithm for: [sig2param]
  # prepare components
  numerator=0
  denominator=0
  for (i in 1:length(Sj)){
    #for (i in 1:length(logSj)){
    numerator = numerator + 
      dlogsknorm_log(Sj[i], matXj[i,], t(as.matrix(betaparam.st[j,])), sig2_prop.st, xi_prev.st[j]) + 
      #dsknorm_log(logSj[i], matXj[i,], t(as.matrix(betaparam[j,])), sig2_p, xi_prev[j]) + 
      log(dinvgamma(sig2_prop.st, shape = u0param.st, scale = v0param.st)) + 
      log(dinvgamma(sig2_prev.st[j], shape = u0param.st, scale = v0param.st))
    #log(dinvgamma(sig2_prev[j], shape = 3*u0param, scale = 2*v0param)) #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~(CHECK)
    
    denominator = denominator + 
      dlogsknorm_log(Sj[i], matXj[i,], t(as.matrix(betaparam.st[j,])), sig2_prev.st[j], xi_prev.st[j])+ 
      #dsknorm_log(logSj[i], matXj[i,], t(as.matrix(betaparam[j,])), sig2_prev[j], xi_prev[j])+ 
      log(dinvgamma(sig2_prev.st[j], shape = u0param.st, scale = v0param.st)) + 
      log(dinvgamma(sig2_prop.st, shape = u0param.st, scale = v0param.st))
    #log(dinvgamma(sig2_prop, shape = 3*u0param, scale = 2*v0param))   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~(CHECK)
  }
  # compute the ratio
  ratio_sig2.st = min(exp(numerator-denominator), 1)
  
  U = runif(n = 1, min = 0, max = 1)
  if(U < ratio_sig2.st) {
    sig2param.st[j] = sig2_prop.st
  } 
  else {
    sig2param.st[j] = sig2_prev.st[j]
  }
  
  
  
  ##> C. MH algorithm for: [xiparam]
  # prepare components
  numerator=0
  denominator=0
  for (i in 1:length(Sj)){
    #for (i in 1:length(logSj)){
    numerator = numerator + 
      dlogsknorm_log(Sj[i], matXj[i,], t(as.matrix(betaparam.st[j,])), sig2param.st[j], xi_prop.st) +
      #dsknorm_log(logSj[i], matXj[i,], t(as.matrix(betaparam[j,])), sig2param[j], xi_p) +
      dt(xi_prop.st, df=nu0.st, log = TRUE) + dt(xi_prev.st[j], df=nu0.st, log = TRUE)
    #dt(xi_prop, df=nu0, log = TRUE) + dnorm(xi_prev[j], mean = 0, sd = nu0, log = TRUE) 
    #dt(xi_prop, df = nu0, log = TRUE) + dnorm(xi_prev[j], mean = 0, sd = nu0/10, log = TRUE) #~~~~~~(CHECK)
    denominator = denominator + 
      dlogsknorm_log(Sj[i], matXj[i,], t(as.matrix(betaparam.st[j,])), sig2param.st[j], xi_prev.st[j]) +
      #dsknorm_log(logSj[i], matXj[i,], t(as.matrix(betaparam[j,])), sig2param[j], xi_prev[j]) +
      dt(xi_prev.st[j], df=nu0.st, log = TRUE) + dt(xi_prop.st, df=nu0.st, log = TRUE)
    #dt(xi_prev[j], df=nu0, log = TRUE) + dnorm(xi_prop, mean = 0, sd = nu0, log = TRUE)
    #dt(xi_prev[j], df=nu0, log = TRUE) + dnorm(xi_prop, mean = 0, sd = nu0/10, log = TRUE)  #~~~~~~(CHECK)
  }
  # compute the ratio
  ratio_xi.st = min(exp(numerator-denominator), 1)
  
  U = runif(n = 1, min = 0, max = 1)
  if(U < ratio_xi.st) {
    xiparam.st[j] = xi_prop.st
  } 
  else {
    xiparam.st[j] = xi_prev.st[j]
  }
}
betaparam.st; sig2param.st; xiparam.st
beta_prev.st; sig2_prev.st; xi_prev.st










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
total_iter = 50
# r_convergence = 1000 <- first run Gibbs sampler to determine r value for convergence
r_convergence = 0

# loglikelihood = numeric(total_iter)                  # for monitor convergence
loglikelihood.st = matrix(0, nrow = n.train, ncol = total_iter)
# [note] why matrix? each col gives each iteration result



#########################] Prepare a "pool"(list) to store the ultimate hyperparameters, parameters [##########################

#> w/o cluster ++++++++++++++++++++++++++++++++ (hyperparameter)
list_beta0param.st = list() # For beta0 hyperparameter
list_SIG_b0param.st = list() # For SIG_b0 hyperparameter
list_u0param.st = list() # for u0 hyperparameter
list_v0param.st = list() # for v0 hyperparameter

#> w/o cluster +++++++++++++++++++++++++++++++++++++++ (outcome)
list_betaparamfull.st = list() # for S - full data
list_sig2paramfull.st = list()  # for S - full data
list_xiparamfull.st = list()
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
list_betaparam.st = list()    #for S
list_sig2param.st = list()    #for S
list_xiparam.st = list()

# -------------------- after correction ---------------------- #
list_betaparam.clean = list()    
list_sig2param.clean = list()    
list_xiparam.clean = list()

####> with cluster-wise +++++++++++++++++++++++++++++ (covariate) 
list_piparam.st = list()       #for ZS
list_lambda2param.st = list()  #for XS
list_kappaparam.st = list()             #for XS

list_tau2param = list()         #for measurement model

# -------------------- after correction ---------------------- #
list_piparam.clean = list()        #for Z
list_lambda2param.clean = list()
list_kappaparam.clean = list()     #for X|Z

###########################################################################################################################


#--------------------------------------------------------------------- ALWAYS KEEP and NO-update !!!! --------------------------
### What we have so far................................................................................ for {initialization}
###] Prepare for the current versions of single "cluster-wise" parameters as the "previous" +++++++++++++++++++++++++++++++++

#> hyperparameter
beta0_prev.st = beta0param.st
SIG_b0_prev.st = SIG_b0param.st
u0_prev.st = u0param.st
v0_prev.st = v0param.st

#> outcome
beta_prev.st = betaparam.st
sig2_prev.st = sig2param.st
xi_prev.st = xiparam.st
#beta_prev_fullS.st = rmvn(n=1, mu=beta0_prevS.st, sigma=SIG_b0_prevS.st) # prebviously given!!!!
#sig2_prev_full.st = rgamma(n=1, shape=u0_prevS.st, rate=v0_prevS.st)     # prebviously given!!!!


#> covariate ---------------------------------------------------------
piparam_prev.st = piparam.st           # we do not need. unaffected.
piparam_prev.clean = piparam.st        # we do not need. unaffected.

lambda2param_prev.st = lambda2param.st # do we need???? 
lambda2param_prev.clean = lambda2param.st #????

kappaparam_prev.st = kappaparam.st
kappaparam_prev.clean = kappaparam.st

tau2_prev = tau2param  # tau2param:::(1-zeta_scale)*lambda2paramS.st
#---------------------------------------------------------------------
#--------------------------------------------------------------------- ALWAYS KEEP and NO-update !!!! --------------------------


###] Miscellaneous ???

max_iter = 10   # let's say.. investigation( sampling rejection ) limit...?
counts_mat <- matrix(0, nrow = J*total_iter, ncol = 3) # count..How many rejected???

n_j.st = table(cl_membership)

beta_p_accepted.st <- 0
sig2_p_accepted.st <- 0
xi_p_accepted.st <- 0

beta_p_accepted.clean <- 0
sig2_p_accepted.clean <- 0
xi_p_accepted.clean <- 0

diff_lambda2 <- matrix(0, nrow = n.train, ncol = total_iter)
diff_tau2 <- matrix(0, nrow = n.train, ncol = total_iter)

n_j.st # just temporary...for cluster-wise sampling...

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
################################## START #######################################

for (r in 1:total_iter) {
  
  kappaparamclean = matrix(0, nrow = J, ncol = ncol(kappaparam.st))     #for X|Z vs Xstar|Z
  lambda2paramclean = numeric(J)                                        #for X|Z vs Xstar|Z
  
  betaparamclean =  matrix(0, nrow = J, ncol = ncol(betaparam.st))      #for S
  sig2paramclean = numeric(J)                                           #for S
  xiparamclean = numeric(J)                                            #for S
  xclean = numeric(n.train)                                             #for X
  
  
  for(j in 1:J) {
    
    Xstar_j=Xstar[cl_membership==j]
    Z_j = Z[cl_membership==j]
    
    difference1 = -1
    difference2 = -1
    difference3 = -1
    
    count_iterations = 0
    
    ###> Don't stop until all the three positive values are collected.  
    while( difference1 < 0 | difference2 < 0 | difference3 < 0) {
      
      ###[1] Updating Posterior -------------------------------------------------------------------------------
      
      #>>> First, sample ***HYPERPARAMETERS*** from hyperPosterior based on full data {w/o clusters}
      
      ## hyperparam for [betaparamS]
      SIG_b0param.st = riwish(v=qu0.st+2, 
                              S=t(beta0_prev.st-beta_prev_full.st) %*% (beta0_prev.st-beta_prev_full.st) +
                                delta*t(beta0_prev.st - m0.st) %*% (beta0_prev.st - m0.st) + LAMB.st) 
      
      beta0param.st = rmvn(n=1, mu=delta/(delta+1)*m0.st + 1/(delta+1)*beta_prev_full.st, 
                           sigma=1/(delta+1)*SIG_b0param.st) 
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~(CHECK)
      
      ## hyperparam for [sig2param]  
      #> - Sample u0 from MH using Gamma(3,1) as proposal 
      u0_prop.st = rgamma(n=1, shape=u0_proposed_shape.st, rate=u0_proposed_rate.st)
      ratio_u0.st = min(fink(u0_prop.st, A = sig2_prev_full.st*rho_u1.st/(v0_prev.st/2), B = rho_u2.st+1)/
                          fink(u0_prev.st, A = sig2_prev_full.st*rho_u1.st/(v0_prev.st/2), B = rho_u2.st+1)*
                          dgamma(x = u0_prev.st, shape = u0_proposed_shape.st, rate = u0_proposed_rate.st)/
                          dgamma(x = u0_prop.st, shape = u0_proposed_shape.st, rate = u0_proposed_rate.st), 1)
      if(runif(n = 1, min = 0, max = 1) < ratio_u0.st) {
        u0param.st = u0_prop.st
      } else {
        u0param.st = u0_prev.st
      }
      
      #> - Sample v0 from Gamma distribution 
      v0param.st = rgamma(n=1, shape = rho_v1.st + u0param.st/2, rate = rho_v2.st+1/sig2_prev_full.st)
      
      
      #------------------------------------------------------------------------------------------------------------------
      #------------------------------------------------------------------------------------------------------------------
      
      #>>> Next I, Sample ****[Main parameters]**** 
      
      #> - Proposal Sampling using proposals OF PRIORS {w/o clusters} 
      beta_prop_full.st = rmvn(n=1, mu=beta0param.st, sigma=SIG_b0param.st*varinf) #~~~~~~~~~~~~~~~~~~~(CHECK)
      sig2_prop_full.st = rinvgamma(n=1, shape=u0param.st, scale=v0param.st)
      xi_prop_full.st = rt(n=1, df=nu0.st)
      
      #> - Sifting a single [betaparam] based on FULL data 
      numerator=0
      denominator=0
      # plugging prior sample first...
      for (i in 1:n.train){
        numerator = numerator + 
          dlogsknorm_log(S[i], matXstar[i,], t(as.matrix(beta_prop_full.st)), sig2_prev_full.st, xi_prev_full.st)
        
        denominator = denominator + 
          dlogsknorm_log(S[i], matXstar[i,], t(as.matrix(beta_prev_full.st)), sig2_prev_full.st, xi_prev_full.st)
      }
      ratio_beta.st = min(exp(numerator-denominator), 1)
      
      U = runif(n = 1, min = 0, max = 1)
      if(U < ratio_beta.st) {
        betaparam_full.st = beta_prop_full.st
      } 
      else {
        betaparam_full.st = beta_prev_full.st
      }
      
      
      #> - Sifting a single [sig2param] based on FULL data 
      numerator=0
      denominator=0
      # plugging prior sample first...
      for (i in 1:n.train){
        numerator = numerator + 
          dlogsknorm_log(S[i], matXstar[i,], t(as.matrix(betaparam_full.st)), sig2_prop_full.st, xi_prev_full.st)
        
        denominator = denominator + 
          dlogsknorm_log(S[i], matXstar[i,], t(as.matrix(betaparam_full.st)), sig2_prev_full.st, xi_prev_full.st)
      }
      ratio_sig2.st = min(exp(numerator-denominator), 1)
      
      U = runif(n = 1, min = 0, max = 1)
      if(U < ratio_sig2.st) {
        sig2param_full.st = sig2_prop_full.st
      } 
      else {
        sig2param_full.st = sig2_prev_full.st
      }
      
      
      #> - Sifting a single [xiparam] based on FULL data 
      numerator=0
      denominator=0
      # plugging prior sample first...
      for (i in 1:n.train){
        numerator = numerator + 
          dlogsknorm_log(S[i], matXstar[i,], t(as.matrix(betaparam_full.st)), sig2param_full.st, xi_prop_full.st)
        denominator = denominator + 
          dlogsknorm_log(S[i], matXstar[i,], t(as.matrix(betaparam_full.st)), sig2param_full.st, xi_prev_full.st)
      }
      # compute the ratio
      ratio_xi.st = min(exp(numerator-denominator), 1)
      
      U = runif(n = 1, min = 0, max = 1)
      if(U < ratio_xi.st) {
        xiparam_full.st = xi_prop_full.st
      } 
      else {
        xiparam_full.st = xi_prev_full.st
      }
      
      
      
      
      #--------------------------------------------------------------------------------------------------------#
      #--------------------------------------------------------------------------------------------------------#
      #-----------------------------------Now it's time for cluster-wise---------------------------------------#
      #--------------------------------------------------------------------------------------------------------#
      #--------------------------------------------------------------------------------------------------------#
      
      
      #>>> First.., Sample **[covariate parameters_j]**
      
      #Xstar_j=Xstar[cl_membership==j] #: already defined
      #Z_j = Z[cl_membership==j]        #: already defined
      #n_j.st[j]                         #: already defined
      
      #> - proposal Sampling [piparamF_j] using proposals OF POSTERIOR {with clusters} 
      piparam.st[j] = rbeta( n=1, shape1=g0.st+sum(Z_j), shape2=h0.st+n_j.st[j]-sum(Z_j) )
      # no need...previous value...?
      
      #> - proposal Sampling [lambda2paramF_j] using proposals OF POSTERIOR {with clusters} 
      lambda2param.st[j] = rinvgamma( n=1, 
                                      shape=(c0.st+n_j.st[j])/2, 
                                      scale=(d0.st + 
                                               sum( (Xstar_j-kappaparam_prev.st[j,1]-kappaparam_prev.st[j,2]*Z_j)^2 ))/2
      )
      
      #> - proposal Sampling [tau2param_j] using proposals OF POSTERIOR {with clusters}  
      # + VERSION 01 +                              
      # tau2param[j] = rinvgamma( n=1, 
      #                           shape=(c0S.st+nS_j.st[j])/2, 
      #                           scale=(1-zeta_scale)*(d0S.st + 
      #                                    sum( (XSstar_j-kappaparam_prev.st[j,1]-kappaparam_prev.st[j,2]*ZSj)^2 ))/2
      # )
      # + VERSION 02 + 
      tau2param[j] = rinvgamma( n=1, 
                                shape=(c0.st+n_j.st[j])/2, 
                                scale=(1-zeta_scale)*(d0.st + 
                                                        sum( (Xstar_j-kappaparam_prev.st[j,1]-kappaparam_prev.st[j,2]*Z_j)^2 ))/2
      )
      #> - proposal Sampling [kappaparam]_j
      KK1j = cbind(1, Z_j)
      KK2j = matrix(c(sum(Xstar_j), sum(Xstar_j*Z_j)), nrow=2, ncol=1) # 2x1 vector
      SIG_inv_j = solve(SIG_inv_k0+(t(KK1j) %*% KK1j)) # solve( ..[matrix].. ) = [matrix]^(-1)
      kappaparam.st[j,] <- rmvn(n=1, mu=SIG_inv_j %*% (SIG_inv_k0 %*% kappa_v0 + KK2j), 
                                sigma = lambda2param.st[j]*SIG_inv_j) #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~(CHECK)
      
      
      #---------------------------------------------------------------------------------------------------------------
      #---------------------------------------------------------------------------------------------------------------
      #>>> Next,... Sample ****[Main parameters_j]****
      
      #>>> subsetting by cluster
      S_j = S[cl_membership==j]                  #outcome
      matX_j.st = matXstar[cl_membership==j, ] #covariate
      
      #>>> Sample "new" outcome paramters from PROPOSALS (priors)
      beta_prop.st = rmvn(n=1, mu=beta0param.st, sigma=SIG_b0param.st*varinf) #~~~~~~~~~~~~~~~~~~~~~~~~~(CHECK)
      sig2_prop.st = rinvgamma(n=1, shape=u0_new.st, scale=v0_new.st)
      xi_prop.st = rt(n=1, df=nu0.st)
      
      
      #>>> - Sifting [betaparam]_j
      numerator=0
      denominator=0
      for (i in 1:length(S_j)){
        numerator = numerator + 
          dlogsknorm_log(S_j[i], matX_j.st[i,], t(as.matrix(beta_prop.st)), sig2_prev.st[j], xi_prev.st[j])
        
        denominator = denominator + 
          dlogsknorm_log(S_j[i], matX_j.st[i,], t(as.matrix(beta_prev.st[j,])), sig2_prev.st[j], xi_prev.st[j])
      }
      ratio_beta.st = min(exp(numerator-denominator), 1)
      U = runif(n = 1, min = 0, max = 1)
      if(ratio_beta.st > U) {
        betaparam.st[j, ] = beta_prop.st
        beta_p_accepted.st <- beta_p_accepted.st+1
      } else {
        betaparam.st[j,] = beta_prev.st[j,]
      }
      
      
      #>>> - Sifting [sig2param]_j
      numerator=0
      denominator=0
      for (i in 1:length(S_j)){
        numerator = numerator + 
          dlogsknorm_log(S_j[i], matX_j.st[i,], t(as.matrix(betaparam.st[j,])), sig2_prop.st, xi_prev.st[j]) +
          log(dinvgamma(sig2_prop.st, shape=u0param.st, scale=v0param.st)) + 
          log(dinvgamma(sig2_prev.st[j], shape=u0_new.st, scale=v0_new.st))
        
        denominator = denominator + 
          dlogsknorm_log(S_j[i], matX_j.st[i,], t(as.matrix(betaparam.st[j,])), sig2_prev.st[j], xi_prev.st[j]) +
          log(dinvgamma(sig2_prev.st[j], shape=u0param.st, scale=v0param.st)) + 
          log(dinvgamma(sig2_prop.st, shape=u0_new.st, scale=v0_new.st))
      }
      ratio_sig2.st = min(exp(numerator-denominator), 1)
      U = runif(n = 1, min = 0, max = 1)
      if(U < ratio_sig2.st) {
        sig2param.st[j] = sig2_prop.st
        sig2_p_accepted.st <- sig2_p_accepted.st+1
      }
      else {
        sig2param.st[j] = sig2_prev.st[j]
      }
      
      
      #>>> - Sifting [xiparam]_j
      numerator=0
      denominator=0
      for (i in 1:length(S_j)){
        numerator = numerator + 
          dlogsknorm_log(S_j[i], matX_j.st[i,], t(as.matrix(betaparam.st[j,])), sig2param.st[j], xi_prop.st)
        denominator = denominator + 
          dlogsknorm_log(S_j[i], matX_j.st[i,], t(as.matrix(betaparam.st[j,])), sig2param.st[j], xi_prev.st[j])
      }
      # compute the ratio
      ratio_xi.st = min(exp(numerator-denominator), 1)
      
      U = runif(n = 1, min = 0, max = 1)
      if(U < ratio_xi.st) {
        xiparam.st[j] = xi_prop.st
        xi_p_accepted.st <- xi_p_accepted.st+1
      } 
      else {
        xiparam.st[j] = xi_prev.st[j]
      }
      
      
      #################################################################################################################
      #################################################################################################################
      #################################################################################################################
      #################################################################################################################
      ###> [2] Calculate differences to ensure VALIDITY of [System Equation]
      #################################################################################################################
      #>>> [main] I. For clean beta_j_1 (the problematic***slope on clean X)
      betaparamclean_j = betaparam.st[j,2]*lambda2param.st[j]/(lambda2param.st[j]-tau2param[j])
      
      #>>> [main] II. For clean sig2_j (the problematic***variance param for LSN on clean X)
      sig2paramclean_j = sig2param.st[j] - betaparamclean_j^2/( 1/tau2param[j] + 1/(lambda2param.st[j]-tau2param[j]) )
      
      #>>> [main] III. For clean xi_j
      # ?
      
      
      #>>> [supplementary] How to ensure They are positive ? 
      difference1 = lambda2param.st[j]-tau2param[j] # ensure clean lambda2 is positive
      difference2 = sig2paramclean_j
      difference3 = 1/tau2param[j]+1/(lambda2param.st[j]-tau2param[j]) - xiparam.st[j]^2*betaparamclean_j^2/sig2paramclean_j
      
      
      #> -while looping count- iteration for each j
      print(c(difference1, difference2, difference3))
      count_iterations = count_iterations+1
    } # End of while (Don't stop until all the three differences are positive!!!!)
    print(paste("count=", count_iterations))
    
    # Values that meet the condition are obtained. Now....
    ###> Collect the params and compute the clean parameters using [System Equation].
    #################################################################################################################
    
    #kappaparamclean = matrix(0, nrow = JS, ncol = ncol(kappaparam.st))     # already defined
    #lambda2paramcleanS = numeric(JS)                                       # already defined
    #betaparamcleanS =  matrix(0, nrow = JS, ncol = ncol(betaparamS.st))    # already defined
    #sig2paramclean =numeric(JS)                                           # already defined
    
    #xSclean = numeric(n.train)                                             # already defined
    
    
    lambda2paramclean[j] = lambda2param.st[j] - tau2param[j]    
    
    kappaparamclean[j,] = kappaparam.st[j,]
    
    betaparamclean[j,2] = betaparam.st[j,2]*lambda2param.st[j]/(lambda2param.st[j]-tau2param[j]) #for Beta1
    
    betaparamclean[j,1] = betaparam.st[j,1] - 
      betaparamclean[j,2]*kappaparam.st[j,1]*tau2param[j]/(lambda2param.st[j]-tau2param[j]) #for Beta0
    
    betaparamclean[j,3] = betaparam.st[j,3] - 
      betaparamclean[j,2]*kappaparam.st[j,2]*tau2param[j]/(lambda2param.st[j]-tau2param[j]) #for Beta2
    
    sig2paramclean[j] = sig2param.st[j] - 
      betaparamclean[j,2]^2/( 1/tau2param[j] + 1/lambda2paramclean[j] )   #for sig2
    
    xiparamclean[j] = sign(xiparam[j])*sqrt( xiparam.st[j]^2*(betaparamclean[j,2]^2/sig2paramclean[j]+1/tau2param[j]+1/lambda2paramclean[j])/
                                               (1/tau2param[j]+1/lambda2paramclean[j] - xiparam.st[j]^2*betaparamclean[j,2]^2/sig2paramclean[j]) )
    
    #################################################################################################################
  } # End of [j]
  
  
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
  
  #-tau2_prev = tau2param  # tau2param:::(1-zeta_scale)*lambda2paramS.st #------------already defined
  
  
  ##]]]]]] - The truth is.......                                     
  #>>>> WHY? ++++++++++++++++++++++++++++++++++++++++++ to sample [betaparamS] and [sig2param] next time....[O]: ONLY USEFUL
  beta0_prev.st = beta0param.st #:::::::::::::::: useful for proposal sample of "SIG_b0paramS.st" 
  SIG_b0_prev.st = SIG_b0param.st
  u0_prev.st = u0param.st #::::::::::::::::::::::::::::::::::::::::: useful to sift "u0paramS.st"
  v0_prev.st = v0param.st #::::::::::::::::::::::::::::::::::::::::: useful to sift "u0paramS.st"
  
  #>>>> WHY? ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ to sample [Y] next time....but [X]
  beta_prev.st = betaparam.st # ::::::::::::::::::::::::::::::::::: useful to sift "betaparamS.st"
  sig2_prev.st = sig2param.st # ::::::::::::::::::::::::::::::::::::: useful to sift "betaparamS.st" and "sig2param.st"
  xi_prev.st = xiparam.st
  #beta_prev_fullS.st  # this is not updated!!!
  #sig2_prev_full.st   # this is not updated!!!
  
  #>>>> WHY? +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ to sample [ZS] next time....but [X]
  #piparam_prevS.st = piparamS.st           
  #piparam_prevS.clean = piparamS.st                            
  
  #>>>> WHY? ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ to sample [XS|ZS] next time....but [X]
  lambda2param_prev.st = lambda2param.st
  lambda2param_prev.clean = lambda2paramclean #***********
  kappaparam_prev.st = kappaparam.st # ::::::::::::: useful for proposal sample of "lambda2paramS", and thus "tau2param"     
  kappaparam_prev.clean = kappaparamclean #*****************
  
  #>>>> WHY? ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ to sample [XSstar|XS] next time....but [X]                                       
  tau2_prev = tau2param                                  
  
  
  
  
  ###############################################################################################################
  ###[3] Calculating Loglikelihood to monitor the chain convergence ---------------------------------------------
  ###############################################################################################################                                        
  #>>> First
  loglike_r = numeric(n.train)
  
  for(i in 1:n.train) {
    x_i = matXstar[i, ]      # first row..{1,x,z}
    j = cl_membership[i] # membership ID
    
    loglike_r[i] = dlogsknorm_log(S[i], x_i, betaparam.st[j,], sig2param.st[j], xiparam.st[j]) + 
      dnorm(x=x_i[2], mean=kappaparam.st[j,1] + kappaparam.st[j,2]*Z[i], sd=sqrt(lambda2param.st[j]), log=TRUE) +
      dbinom(x=x_i[3], size=1, prob=piparam.st[j], log=TRUE) 
  }  
  loglikelihood.st[,r] = loglike_r
  
  if(r > r_convergence) {
    
    # W/O CLUSTER
    list_beta0param.st[[r-r_convergence]] = beta0param.st
    list_SIG_b0param.st[[r-r_convergence]] = SIG_b0param.st
    list_u0param.st[[r-r_convergence]] = u0param.st
    list_v0param.st[[r-r_convergence]] = v0param.st
    
    list_betaparamfull.st[[r-r_convergence]] = betaparam_full.st
    list_sig2paramfull.st[[r-r_convergence]] = sig2param_full.st
    list_xiparamfull.st[[r-r_convergence]] = xiparam_full.st
    
    # With Cluster-wise
    list_piparam.st[[r-r_convergence]] = piparam.st     
    list_lambda2param.st[[r-r_convergence]] = lambda2param.st
    list_kappaparam.st[[r-r_convergence]] = kappaparam.st
    list_tau2param[[r-r_convergence]] = tau2param
    
    list_betaparam.st[[r-r_convergence]] = betaparam.st  
    list_sig2param.st[[r-r_convergence]] = sig2param.st
    list_xiparam.st[[r-r_convergence]] = xiparam.st
    
    ################# STORE the CLEAN parameters ##################
    list_kappaparam.clean[[r-r_convergence]] = kappaparamclean     #for Xstar
    list_lambda2param.clean[[r-r_convergence]] = lambda2paramclean    #for Xstar
    list_betaparam.clean[[r-r_convergence]] = betaparamclean    #for S
    list_sig2param.clean[[r-r_convergence]] = sig2paramclean    #for S
    list_xiparam.clean[[r-r_convergence]] = xiparamclean        #for S
    
    # list_xclean[[r-r_convergence]] = xclean #~xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx Q.DAVID (can we impute?)
    
  }
  print(paste("r=",r))
  trace(what = "print", where = getNamespace("base"), exit = flush.console, print = FALSE)
}

########################################
####### End of Gibbs Sampler ###########
########################################

beta_p_accepted.st/(total_iter*J)
sig2_p_accepted.st/(total_iter*J)
xi_p_accepted.st/(total_iter*J)

par(mfrow=c(1,1))
plot(colSums(loglikelihood.st),type="l")

lppd_EX_S.S.st = mean( colSums(loglikelihood.st) ); lppd_EX_S.S.st
# -25719.65(0.9),           


#-------------------------------------------------------------------------------

############################### Predictions ####################################

# Expected Value of LSN distribution (Wang 2019)
# 2*exp(sum(x_i*beta_j[j,]) + sig2_j[j]/2)*(1-pnorm(-xi_j[j]*sqrt(sig2_j[j])/sqrt(xi_j[j]^2+1)))


###> Training data
expval.train.st = matrix(0, nrow = n.train, ncol = total_iter-r_convergence)

for(r in 1:(total_iter-r_convergence)) {
  betaparamclean = list_betaparam.st[[r]]
  sig2paramclean = list_sig2param.st[[r]]
  xiparamclean = list_xiparam.st[[r]]
  
  for(i in 1:n.train) {
    j = cl_membership[i]
    #cleanx = rnorm(n = 1, mean = kappaparamclean[j,1] + kappaparamclean[j,2]*Z[i], sd = sqrt(lambda2paramclean[j]))
    cleanx = Xstar[i]
    #expval.train[i,r] <- 2*exp( betaparamclean[j,1] + betaparamclean[j,2]*cleanx + betaparamclean[j,3]*Z[i] + 
    #                              sig2paramclean[j]/2)*( 1-pnorm(-xiparamclean[j]*
    #                                                               sqrt(sig2paramclean[j])/sqrt(xiparamclean[j]^2+1)))
    expval.train.st[i,r] <- 10*exp( betaparamclean[j,1] + betaparamclean[j,2]*cleanx + betaparamclean[j,3]*Z[i] - 
                                      sig2paramclean[j]/2)*( 1-pnorm(-xiparamclean[j]*
                                                                       sqrt(sig2paramclean[j])/sqrt(xiparamclean[j]^2+1)))
  }
}

expS.train.st <- apply(X=expval.train.st, MARGIN=1, FUN=mean)

summary(S)
summary(expS.train.st)

format( aggregate(log(S), by=list(cl_membership), FUN=summary), scientific=FALSE)
format( aggregate(log(expS.train.st), by=list(cl_membership), FUN=summary), scientific=FALSE)

lppd_EX_S.S.st <- mean( colSums(loglikelihood.st) ); lppd_EX_S.S.st # -25719.65

# Clusterwise predictive density plots
par(mfrow = c(2,3))

for(j in 1:J) {
  hist(log(S)[cl_membership==j], freq=FALSE, breaks=50, col="white", xlim=c(5,15), ylim=c(0,1.6),
       xlab = paste("Cluster", j))
  lines(density(log(expS.train.st)[cl_membership==j]), col="blue", ylim=c(0,1.6), lwd=2)
}


###> Testing data
expval.test.st = matrix(0, nrow = n.test, ncol = total_iter-r_convergence)

for(r in 1:(total_iter-r_convergence)) {
  betaparamclean = list_betaparam.st[[r]]
  sig2paramclean = list_sig2param.st[[r]]
  xiparamclean = list_xiparam.st[[r]]
  
  for(i in 1:n.test) {
    j = cl_membership[i]
    #cleanx = rnorm(n = 1, mean = kappaparamclean[j,1] + kappaparamclean[j,2]*Z[i], sd = sqrt(lambda2paramclean[j]))
    cleanx = Xstar.test[i]
    expval.test.st[i,r] <- 10*exp( betaparamclean[j,1] + betaparamclean[j,2]*cleanx + betaparamclean[j,3]*Z[i] - 
                                     sig2paramclean[j]/2)*( 1-pnorm(-xiparamclean[j]*
                                                                      sqrt(sig2paramclean[j])/sqrt(xiparamclean[j]^2+1)))
  }
}
expS.test.st <- apply(X=expval.test.st, MARGIN=1, FUN=mean)

summary(S.test)
summary(expS.test.st)

format( aggregate(log(S.test), by=list(cl_membership.test), FUN=summary), scientific=FALSE)
format( aggregate(log(expS.test.st), by=list(cl_membership.test), FUN=summary), scientific=FALSE)


# Clusterwise predictive density plots
par(mfrow = c(2,3))

for(j in 1:J) {
  hist(log(S.test)[cl_membership.test==j], freq=FALSE, breaks=50, col="white", xlim=c(5,15), ylim=c(0,1.6),
       xlab = paste("Cluster", j))
  lines(density(log(expS.test.st)[cl_membership.test==j]), col="blue", ylim=c(0,1.6), lwd=2)
}




################### Clusterwise predictive density plots #######################
par(mfrow = c(2,3))

for(j in 1:J) {
  hist(log(S)[cl_membership==j], freq=FALSE, breaks=50, col="white", xlim=c(5,15), ylim=c(0, 2.8),
       xlab = paste("Cluster", j))
  #lines(density(log(expS.train)[cl_membership==j]), col="red", lwd=2)     # clean 
  lines(density(log(expS.train)[cl_membership==j]), col="black", lty = "dotted") # dirty
  lines(density(log(expS.train.st)[cl_membership==j]), col="blue", lwd=1) # dirty with correction
  
  if (j == 1) {
    legend("topright", legend = c("Clean", "Correction"), #, "Correction Medium", "Correction Small"),
           col = c("black", "blue"), #, "darkmagenta", "darkgreen"), 
           lty = c(1, 2), 
           lwd = c(3, 2), 
           cex = 0.4
    )
  }
}





###> LPPD for Bayes model (for S)
lppd_EX_S.S   #-25705.38
lppd_EX_S.S.st#-25719.65


###> SSSE
#> The model
expS.test
SSPE.LSN <- sum( (log(expS.test) - log(S.test))^2 ); SSPE.LSN   #732.0648
SAPE.LSN <- sum( abs(log(expS.test) - log(S.test)) ); SAPE.LSN  #335.7326

#> Gustafson Correction
expS.test.st
SSPE.LSN.st <- sum( (log(expS.test.st) - log(S.test))^2 ); SSPE.LSN.st   #858.6981
SAPE.LSN.st <- sum( abs(log(expS.test.st) - log(S.test)) ); SAPE.LSN.st  #367.862


###> CTE (for S)
summary( expS.train )
summary( expS.train.st )
summary( expS.test )
summary( expS.test.st )

# - Calculate the quantile threshold
alalpha = 0.95
alalpha = 0.9
alalpha = 0.5
alalpha = 0.1
quantile_threshold <- quantile(expS.train, alalpha); quantile_threshold
quantile_threshold.st <- quantile(expS.train.st, alalpha); quantile_threshold.st
tail_data <- expS.train[expS.train > quantile_threshold]; tail_data # NA?????
tail_data.st <- expS.train.st[expS.train.st > quantile_threshold.st]; tail_data.st # NA?????
CTE.HB <- mean(tail_data, na.rm = TRUE); CTE.HB          # 4511587  
CTE.HB.st <- mean(tail_data.st, na.rm = TRUE); CTE.HB.st # 3809959


###> Calculate the KL divergence (for S)
library(entropy)

pdf <- density(.....)$y
pdf.HB <- density(expS.train)$y
pdf.HB.st <- density(expS.train.st)$y

kl.divergence.HB <- entropy::KL.empirical(pdf, pdf.HB)
print(kl.divergence.HB) 
kl.divergence.HB.st <- entropy::KL.empirical(pdf, pdf.HB.st)
print(kl.divergence.HB.st) 





