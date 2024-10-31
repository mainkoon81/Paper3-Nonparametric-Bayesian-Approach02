############################### shard test #####################################
library(mvnfast)
library(tidyverse)
library(Rcpp)
library(RcppArmadillo)
#library(extraDistr)
library(MCMCpack)

library(statmod)
library(tweedie)

#> C++ for Gold standard
sourceCpp("C:/Users/kimm4/Desktop/WORKSPACE/2nd 3rd paper/DATASET/Basic_training/project3.DPcpp.Basic.cpp")
#sourceCpp("project3.Benchmark.DPcpp.cpp")

#> C++ for Gustafson Correction
sourceCpp("C:/Users/kimm4/Desktop/WORKSPACE/2nd 3rd paper/DATASET/Basic_training/project3.DPcpp.NDB.cpp")
#sourceCpp("project3.DPcpp.NDB.cpp")







####> SHARD set-up with "train set"
####> 
####> 
####> 
#------------------------------------------------------------------------------#
N <- n.train # Sample size for training data
ind <- 1:N # index for each observation in training

# Want the shard size to be (5006)
# 5006/2 -> subset size would be 2503
# 60063/2503 = 24 -- want 24 subsets and 23 shards -- 23 parallel CPUS

# Want the shard size to be (10412)
# 10412/2 -> subset size would be 5206
# 60063/5206 = 12 -- want 12 subsets and 11 shards -- 11 parallel CPUS ~~~~~~~~~~ Q. where is anchor????


###> Assign random subset memberships to each observation in the training set.
set.seed(1)
Nsh <- 11 # number of shards
subset.idx <- sample(x = 1:(Nsh+1), size = N, replace = TRUE); head(subset.idx)
table(subset.idx)

###> Collect the indices of anchor points
anchor <- ind[subset.idx==Nsh+1]; head(anchor)

###> Collect obv + anchor indices for each shards
shard.idx.list <- list()
for(i in 1:Nsh) {
  shard.idx.list[[i]] <- c(ind[subset.idx==i|subset.idx==Nsh+1])
}; head(shard.idx.list) 

head(shard.idx.list[[1]],20)
head(shard.idx.list[[2]],20)
unlist(lapply(shard.idx.list, length)) # Check sample size of each shard

#------------------------------------------------------------------------------#

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# choose a shard
sh <- 11 # first shard (1,2,3,4,5,6,7,8,9,10,11)
shard.sh.idx <- shard.idx.list[[sh]]
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# brazil1 <- read.csv("Brazil.paper3-1.train.csv")
# S <- brazil1$AggClaim
# Z <- brazil1$Affiliation
# X_tru <- brazil1$ln_Expos
# Xstar <- brazil1$ln_Expos_err

# Subset full sample by individual shard  
S.sh <- S[shard.sh.idx]
Z.sh <- Z[shard.sh.idx]
X_tru.sh <- X_tru[shard.sh.idx]
Xstar.sh <- Xstar[shard.sh.idx]

f0S.sh <- f0S[shard.sh.idx]
f0X.sh <- f0X[shard.sh.idx]

f0S.xx.sh <- f0S.xx[shard.sh.idx]
f0X.xx.sh <- f0X.xx[shard.sh.idx]

f0S.st.sh <- f0S.st[shard.sh.idx]
f0X.st.sh <- f0X.st[shard.sh.idx]

cl_membership.train.sh <- cl_membership.train[shard.sh.idx]
cl_membership.train.xx.sh <- cl_membership.train[shard.sh.idx]
cl_membership.train.st.sh <- cl_membership.train[shard.sh.idx]

cl_membership.HB.train.sh <- cl_membership.HB.train[shard.sh.idx]


matX_tru.sh <- as.matrix(cbind(1, X_tru.sh, Z.sh)) 
matXstar.sh = as.matrix(cbind(1, Xstar.sh, Z.sh))

n.train.sh <- length(S.sh)
varinf.sh <- n.train.sh/100
varinf.ka.sh <- n.train.sh/250 ###> moderate clusters production
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
table(cl_membership.train.sh)
























#####> Gibbs 01.Gold Standard
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Initialize Parameters based on cl_membership.train.sh
# Run Gibbs sampler to get the parameters for each cluster on shard sh
# Parameters for beta, sig2, xi, kappa, lambda2, tau2
J=4
set.seed(1)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
total_iter = 2#5000
# r_convergence = 1000 <- first run Gibbs sampler to determine r value for convergence
r_convergence = 0#4000
save_every = 1#10 # 1000/100
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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
  clusterout = clusterDPBasic(S = S.sh, X = X_tru.sh, Z = Z.sh, cl_membership = cl_membership.train.sh, 
                              piparam = piparam, lambda2param = lambda2param, 
                              betaparam = betaparam, sig2param = sig2param, xiparam = xiparam, 
                              alphaparam = alphaparam, 
                              f0X = f0X.sh, f0S = f0S.sh, u0 = u0, v0 = v0, m0 = m0, SIG_b0 = SIG_b0, 
                              nu0 = nu0, g0 = g0, h0 = h0, c0 = c0, d0 = d0, gamma0 = gamma0, psi0 = psi0, 
                              varinf = varinf.sh) #~~~~~Q.David
  
  
  cl_membership.train.sh = clusterout$cl_membership
  J = length(unique(cl_membership.train.sh))             #% in case, cluster scenario changes, reflecting them
  
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
    Zj=Z.sh[cl_membership.train.sh==j]
    nj=length(Zj)
    piparam[j]=rbeta( n=1, shape1=g0+sum(Zj), shape2=h0+nj-sum(Zj) ) #posterior
  }
  #% for X - lambda2param update from MVN, InvGa posterior 
  lambda2param = numeric(J) #; lambda2param
  for (j in 1:J) {
    Xj=X_tru.sh[cl_membership.train.sh==j]
    Zj = Z.sh[cl_membership.train.sh==j]
    nj=length(Xj)
    
    lambda2param[j]=rinvgamma( n=1, shape=(nj+c0)/2, scale=( sum((Xj - mean(X_tru.sh))^2) + d0 )/2 )
  }
  
  #% for alphaparam
  eta = rbeta(n=1, shape1=alpha0+1, shape2=n.train.sh)
  pi_eta = (gamma0+J-1)/(gamma0+J-1 + n.train.sh*(psi0-log(eta)))
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
    Sj = S.sh[cl_membership.train.sh==j]
    if(length(Sj)==1) {
      matXj = matrix( matX_tru.sh[cl_membership.train.sh==j, ], nrow=1 ) # It is a vector added to the matrix
    } else {
      matXj = matX_tru.sh[cl_membership.train.sh==j, ]                   # It is a regular matrix
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
  
  for(i in 1:n.train.sh) {
    x_i = matX_tru.sh[i, ]
    j = cl_membership.train.sh[i]
    
    loglike_r = loglike_r +                                                     
      dlogsknorm_log(s=S.sh[i], x=x_i, beta=betaparam[j,], sig2=sig2param[j], xi=xiparam[j]) +
      dbinom(x = x_i[3], size = 1, prob = piparam[j], log = TRUE) + 
      dnorm(x = x_i[2], mean = mean(X_tru.sh), sd = sqrt(lambda2param[j]), log = TRUE)
  }
  loglikelihood[r] = loglike_r # cumulative sum
  
  if( (r > r_convergence) & (r %% save_every==0) ) {
    list_piparam[[(r-r_convergence)/save_every]] = piparam   #for Z
    
    list_lambda2param[[(r-r_convergence)/save_every]] = lambda2param  #for X
    
    list_alphaparam[[(r-r_convergence)/save_every]] = alphaparam     #for alpha
    
    list_betaparam[[(r-r_convergence)/save_every]] = betaparam    #for S
    list_sig2param[[(r-r_convergence)/save_every]] = sig2param    #for S
    list_xiparam[[(r-r_convergence)/save_every]] = xiparam      #for S
    
    list_cl[[(r-r_convergence)/save_every]] = cl_membership.train.sh
  }
  print(paste("r=",r)) #% iteration progress
}

#plot(apply(loglikelihood, 2, sum), type="l")
options(scipen = 999)
plot(loglikelihood, type="l")
#plot(loglikelihood, type="l", ylim=c(-113942,-108000))
lppd_EX_DP.S <- mean( loglikelihood ); lppd_EX_DP.S

# From Gibbs sampler we need cluster memberships and parameters - [ model(A): Gold standard ] - just once
filename <- paste("testtt/GS.shard(",sh,")data.rds",sep="")
saveRDS(object = list(shard.sh.idx, anchor, list_cl, list_piparam, list_lambda2param, list_alphaparam, list_betaparam, list_sig2param, list_xiparam, N), file = filename)













#####> Gibbs 02.xx - just three[1%, 10%, 25%]
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
J=4
set.seed(1)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
total_iter = 2#5000
# r_convergence = 1000 <- first run Gibbs sampler to determine r value for convergence
r_convergence = 0#4000
save_every = 1#10 # 1000/100
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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
  clusterout = clusterDPBasic(S.sh, Xstar.sh, Z.sh, cl_membership.train.xx.sh,
                              piparam.xx, lambda2param.xx,
                              betaparam.xx, sig2param.xx, xiparam.xx, alphaparam.xx,
                              f0X.xx.sh, f0S.xx.sh, u0.xx, v0.xx, m0.xx, SIG_b0.xx, 
                              nu0.xx, g0.xx, h0.xx, c0.xx, d0.xx, gamma0.xx, psi0.xx, varinf.sh) 
  
  cl_membership.train.xx.sh = clusterout$cl_membership
  J = length(unique(cl_membership.train.xx.sh))             #% in case, cluster scenario changes, reflecting them
  
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
    Zj=Z.sh[cl_membership.train.xx.sh==j]
    nj=length(Zj)
    piparam.xx[j]=rbeta( n=1, shape1=g0.xx+sum(Zj), shape2=h0.xx+nj-sum(Zj) ) #posterior
  }
  #% for X - lambda2param update from MVN, InvGa posterior 
  lambda2param.xx = numeric(J) #; lambda2param
  for (j in 1:J) {
    Xj=Xstar.sh[cl_membership.train.xx.sh==j]
    Zj = Z.sh[cl_membership.train.xx.sh==j]
    nj=length(Xj)
    
    lambda2param.xx[j]=rinvgamma( n=1, shape=(nj+c0.xx)/2, scale=( sum((Xj - mean(Xstar.sh))^2) + d0.xx )/2 )
  }
  
  #% for alphaparam
  eta.xx = rbeta(n=1, shape1=alpha0.xx+1, shape2=n.train.sh)
  pi_eta.xx = (gamma0.xx+J-1)/(gamma0.xx+J-1 + n.train.sh*(psi0.xx-log(eta.xx)))
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
    Sj = S.sh[cl_membership.train.xx.sh==j]
    if(length(Sj)==1) {
      matXj = matrix( matXstar.sh[cl_membership.train.xx.sh==j, ], nrow=1 ) # It is a vector added to the matrix
    } else {
      matXj = matXstar.sh[cl_membership.train.xx.sh==j, ]                   # It is a regular matrix
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
  
  for(i in 1:n.train.sh) {
    x_i = matXstar.sh[i, ]
    j = cl_membership.train.xx.sh[i]
    
    loglike_r = loglike_r +                                                     
      dlogsknorm_log(s=S.sh[i], x=x_i, beta=betaparam.xx[j,], sig2=sig2param.xx[j], xi=xiparam.xx[j]) +
      dbinom(x = x_i[3], size = 1, prob = piparam.xx[j], log = TRUE) + 
      dnorm(x = x_i[2], mean = mean(Xstar.sh), sd = sqrt(lambda2param.xx[j]), log = TRUE)
  }
  #loglikelihood.xx[, r] = loglike_r
  loglikelihood.xx[r] = loglike_r
  
  if( (r > r_convergence) & (r %% save_every==0) ) {
    list_piparam.xx[[(r-r_convergence)/save_every]] = piparam.xx   #for Z
    
    list_lambda2param.xx[[(r-r_convergence)/save_every]] = lambda2param.xx  #for X
    
    list_alphaparam.xx[[(r-r_convergence)/save_every]] = alphaparam.xx     #for alpha
    
    list_betaparam.xx[[(r-r_convergence)/save_every]] = betaparam.xx    #for S
    list_sig2param.xx[[(r-r_convergence)/save_every]] = sig2param.xx    #for S
    list_xiparam.xx[[(r-r_convergence)/save_every]] = xiparam.xx      #for S
    
    list_cl.xx[[(r-r_convergence)/save_every]] = cl_membership.train.xx.sh
  }
  print(paste("r=",r)) #% iteration progress
}
plot( loglikelihood.xx, type="l")
lppd_EX_DP.S.xx <- mean( loglikelihood.xx ); lppd_EX_DP.S.xx # -?

# From Gibbs sampler we need cluster memberships and parameters - [ model(B): Model risk ] - just three[1%, 10%, 25%]
filename <- paste("testtt/xx.shard_0.01(",sh,")data.rds",sep="")
saveRDS(object = list(shard.sh.idx, anchor, list_cl.xx, list_piparam.xx, list_lambda2param.xx, list_alphaparam.xx, 
                      list_betaparam.xx, list_sig2param.xx, list_xiparam.xx, N), file = filename)


















#####> Gibbs 03.st - multiple[1%X3, 10%X3, 25%X3]
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zeta_scale = 0.99
#zeta_scale = 0.95
#zeta_scale = 0.91 
#zeta_scale = 0.8
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
J=4
set.seed(1)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
total_iter = 2#5000
# r_convergence = 1000 <- first run Gibbs sampler to determine r value for convergence
r_convergence = 0#4000
save_every = 1#10 # 1000/100
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

loglikelihood.st = numeric(total_iter)                  # for monitor convergence
#loglikelihood.st = matrix(0, nrow = n.train.sh, ncol = total_iter)
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

n_j.st = table(cl_membership.train.st.sh)

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
  clusterout = clusterDP(S.sh, Xstar.sh, Z.sh, cl_membership.train.st.sh,
                         piparam.st, kappaparam.st, lambda2param.st,
                         betaparam.st, sig2param.st, xiparam.st, alphaparam.st,
                         f0X.st.sh, f0S.st.sh, u0.st, v0.st, m0.st, SIG_b0.st, nu0.st, g0.st, h0.st,
                         kappa_v0, SIG_k0, SIG_inv_k0, c0.st, d0.st, gamma0.st, psi0.st, varinf.sh)
  
  cl_membership.train.st.sh = clusterout$cl_membership
  J = length(unique(cl_membership.train.st.sh))                  #% in case, cluster scenario changes, reflecting them
  
  #piparam.st = clusterout$piparam
  #lambda2param.st = clusterout$lambda2param
  
  kappaparam_prev.st = clusterout$kappaparam
  beta_prev.st = clusterout$betaparam
  sig2_prev.st = clusterout$sig2param
  xi_prev.st = clusterout$xiparam
  #################################################################################################################  
  
  #%%%%%%%%%%%%%%%%%%%%%%%% Heavy Computation ON/OFF %%%%%%%%%%%%%%%%%%%%%%%%%%#
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%# 
  # for (i in 1:n.train.sh){
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
  #     probs[j] = nj[j]/(n.train.sh-1+alphaparam.st)*
  #       dlogsknorm(s=S[i], x=x_i, beta=betaparam.st[j,], sig2=sig2param.st[j], xi=xiparam.st[j])*
  #       dbinom(x=x_i[3], size=1, prob = piparam.st[j])*                # Covariate Z
  #       dnorm(x=x_i[2], mean=kappaparam.st[j,1]+kappaparam.st[j,2]*x_i[3], sd=sqrt(lambda2param.st[j])) # Covariate X
  #   }
  # 
  #   # c)After iteration through each cluster, Adding the new probability of "Forming a new cluster": P(s_i=J+1)
  #   probs[j+1] = alphaparam.st/(n.train.sh-1+alphaparam.st)*f0S.st[i]*f0X.st[i]   #% so it gives...probs: c(prob, prob, prob, 0) -> c(prob, prob, prob, prob)
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
    Zj=Z.sh[cl_membership.train.st.sh==j]
    nj=length(Zj)
    piparam.st[j]=rbeta( n=1, shape1=g0.st+sum(Zj), shape2=h0.st+nj-sum(Zj) ) #posterior
  }
  
  #% for alphaparam
  eta.st = rbeta(n=1, shape1=alpha0.st+1, shape2=n.train.sh)
  pi_eta.st = (gamma0.st+J-1)/(gamma0.st+J-1 + n.train.sh*(psi0.st-log(eta.st)))
  alphaparam.st = pi_eta.st*rgamma( n=1, shape=gamma0.st+J, 
                                    rate=psi0.st-log(eta.st) ) + (1-pi_eta.st)*rgamma( n=1, shape=gamma0.st+J-1, rate=psi0.st-log(eta.st) )
  
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
    betaprior = matrix( rep(rmvn(n=1, mu=m0.st, sigma=SIG_b0.st*varinf.sh), J-J_prev), nrow=J-J_prev, byrow=TRUE ) #~~~ CHECK
    sig2prior = rep(rinvgamma(n=1, shape=u0.st, scale=v0.st), J-J_prev)
    xiprior = rep(rt(n=1, df=nu0.st), J-J_prev) # imagine we already have them...like old days..
    
    kappaprior = matrix( rep(rmvn(n=1, mu=kappa_v0, sigma=SIG_k0*varinf.ka.sh), J-J_prev), nrow=J-J_prev, byrow=TRUE )#~~ CHECK
    
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
    
    Xstar_j=Xstar.sh[cl_membership.train.st.sh==j]
    Zj = Z.sh[cl_membership.train.st.sh==j]
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
      #clean_x_j <- rnorm(n=nj, mean=Xstar_j, sd=sqrt(tau2_prev_j))
      #tau2param.st[j] = rinvgamma(n = 1, shape = (nj+u0)/2, scale = 1/2*(sum((Xstar_j-clean_x_j)^2)+u0*tau2_0))
      tau2param.st[j]=rinvgamma( n=1, 
                                 shape=(nj+c0.st)/2, 
                                 scale=(1-zeta_scale)*( sum((Xstar_j - kappaparam_prev.st[j,1]-kappaparam_prev.st[j,2]*Zj)^2) + d0.st )/2 
      )
      
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
      S_j = S.sh[cl_membership.train.st.sh==j]                                        # outcome
      matX_j.st = matrix(matXstar.sh[cl_membership.train.st.sh==j, ], nrow=length(S_j))   # covariate
      
      #>>> Sample "new" outcome paramters from PROPOSALS (priors)
      beta_p = rmvn(n=1, mu=m0.st, sigma=SIG_b0.st*varinf.sh) #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ CHECK
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
  #loglike_r = numeric(n.train.sh)
  
  for(i in 1:n.train.sh) {
    x_i = matXstar.sh[i, ]
    j = cl_membership.train.st.sh[i]
    
    loglike_r = loglike_r +                                                     
      #loglike_r[i] =   
      # dlogsknorm_log(s = S[i], x = x_i, beta = betaparam.st[j,], sig2 = sig2param.st[j], xi = xiparam.st[j]) +
      # dnorm(x = x_i[2], mean = kappaparam.st[j,1] + kappaparam.st[j,2]*Z[i], sd = sqrt(lambda2param.st[j]), log = TRUE) +
      # dbinom(x = x_i[3], size = 1, prob = piparam.st[j], log = TRUE)
      dlogsknorm_log(s = S.sh[i], x = x_i, beta = betaparamclean.st[j,], sig2 = sig2paramclean.st[j], xi = xiparamclean.st[j]) +
      dnorm(x = x_i[2], mean = kappaparamclean.st[j,1] + kappaparamclean.st[j,2]*Z.sh[i], sd = sqrt(lambda2paramclean.st[j]), log = TRUE) +
      dbinom(x = x_i[3], size = 1, prob = piparam.st[j], log = TRUE)
  }  
  #print(loglike_r)
  loglikelihood.st[r] = loglike_r
  #loglikelihood.st[,r] = loglike_r
  
  if( (r > r_convergence) & (r %% save_every==0) ) {
    
    list_piparam.st[[(r-r_convergence)/save_every]] = piparam.st      #for Z
    list_kappaparam.st[[(r-r_convergence)/save_every]] = kappaparam.st     #for Xstar
    list_lambda2param.st[[(r-r_convergence)/save_every]] = lambda2param.st  #for Xstar
    list_tau2param.st[[(r-r_convergence)/save_every]] = tau2param.st        # For measurement model
    
    list_alphaparam.st[[(r-r_convergence)/save_every]] = alphaparam.st          #for alpha
    
    list_betaparam.st[[(r-r_convergence)/save_every]] = betaparam.st  #for S
    list_sig2param.st[[(r-r_convergence)/save_every]] = sig2param.st  #for S
    list_xiparam.st[[(r-r_convergence)/save_every]] = xiparam.st      #for S
    
    list_cl.st[[(r-r_convergence)/save_every]] = cl_membership.train.st.sh
    
    ################# STORE the CLEAN parameters ##################
    list_kappaparam.clean.st[[(r-r_convergence)/save_every]] = kappaparamclean.st     #for Xstar
    list_lambda2param.clean.st[[(r-r_convergence)/save_every]] = lambda2paramclean.st    #for Xstar
    list_betaparam.clean.st[[(r-r_convergence)/save_every]] = betaparamclean.st    #for S
    list_sig2param.clean.st[[(r-r_convergence)/save_every]] = sig2paramclean.st    #for S
    list_xiparam.clean.st[[(r-r_convergence)/save_every]] = xiparamclean.st        #for S
    
  }
  print(paste("r=",r))
}
plot(loglikelihood.st, type="l")
lppd_EX_DP.S.st = mean( loglikelihood.st ); lppd_EX_DP.S.st

# From Gibbs sampler we need cluster memberships and parameters - [ model(C): Gustafson DPM ] - multiple[1%X3, 10%X3, 25%X3]
filename <- paste("testtt/st.shard_0.01_0.99(",sh,")data.rds",sep="")
saveRDS(object = list(shard.sh.idx, anchor, list_cl.st, list_piparam.st, list_lambda2param.st, list_kappaparam.st, list_tau2param.st,
                      list_alphaparam.st, list_betaparam.st, list_sig2param.st, list_xiparam.st, 
                      list_kappaparam.clean.st, list_lambda2param.clean.st,
                      list_betaparam.clean.st, list_sig2param.clean.st, list_xiparam.clean.st, N), file = filename)

























#####> Gibbs 04.HB.st - multiple[1%X3, 10%X3, 25%X3]
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zeta_scale = 0.99
#zeta_scale = 0.95
#zeta_scale = 0.91 
#zeta_scale = 0.8
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
total_iter = 2#5000
# r_convergence = 1000 <- first run Gibbs sampler to determine r value for convergence
r_convergence = 0#4000
save_every = 1#10 # 1000/10 = 100
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# loglikelihood = numeric(total_iter)                  # for monitor convergence
loglikelihood.HB.st = matrix(0, nrow = n.train.sh, ncol = total_iter)
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

n_j.HB.st = table(cl_membership.HB.train.sh)

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
  
  xclean = numeric(n.train.sh)                                             #for X
  
  
  
  
  ################################ Real GAME ###################################
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  for(j in 1:J.HB) {
    
    Xstar_j=Xstar.sh[cl_membership.HB.train.sh==j]
    Z_j = Z.sh[cl_membership.HB.train.sh==j]
    
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
      beta_prop_full.HB.st = rmvn(n=1, mu=beta0param.HB.st, sigma=SIG_b0param.HB.st*varinf.sh) #~~~~~~~~~~~~~~~~~~~(CHECK)
      sig2_prop_full.HB.st = rinvgamma(n=1, shape=u0param.HB.st, scale=v0param.HB.st)
      xi_prop_full.HB.st = rt(n=1, df=nu0.HB.st)
      
      #> - Sifting a single [betaparam] based on FULL data 
      numerator=0
      denominator=0
      # plugging prior sample first...
      for (i in 1:n.train.sh){
        numerator = numerator + 
          dlogsknorm_log(S.sh[i], matXstar.sh[i,], t(as.matrix(beta_prop_full.HB.st)), sig2_prev_full.HB.st, xi_prev_full.HB.st)
        
        denominator = denominator + 
          dlogsknorm_log(S.sh[i], matXstar.sh[i,], t(as.matrix(beta_prev_full.HB.st)), sig2_prev_full.HB.st, xi_prev_full.HB.st)
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
      for (i in 1:n.train.sh){
        numerator = numerator + 
          dlogsknorm_log(S.sh[i], matXstar.sh[i,], t(as.matrix(betaparam_full.HB.st)), sig2_prop_full.HB.st, xi_prev_full.HB.st)
        
        denominator = denominator + 
          dlogsknorm_log(S.sh[i], matXstar.sh[i,], t(as.matrix(betaparam_full.HB.st)), sig2_prev_full.HB.st, xi_prev_full.HB.st)
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
      for (i in 1:n.train.sh){
        numerator = numerator + 
          dlogsknorm_log(S.sh[i], matXstar.sh[i,], t(as.matrix(betaparam_full.HB.st)), sig2param_full.HB.st, xi_prop_full.HB.st)
        denominator = denominator + 
          dlogsknorm_log(S.sh[i], matXstar.sh[i,], t(as.matrix(betaparam_full.HB.st)), sig2param_full.HB.st, xi_prev_full.HB.st)
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
      S_j = S.sh[cl_membership.HB.train.sh==j]                  #outcome
      matX_j.HB.st = matXstar.sh[cl_membership.HB.train.sh==j, ] #covariate
      
      #>>> Sample "new" outcome paramters from PROPOSALS (priors)
      beta_prop.HB.st = rmvn(n=1, mu=beta0param.HB.st, sigma=SIG_b0param.HB.st*varinf.sh) #~~~~~~~~~~~~~~~~~~~~~~~~~(CHECK)
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
      #print(c(difference1, difference2, difference3))
      count_iterations = count_iterations+1
    } #%%% END of WHILE %%% 
    #print(paste("count=", count_iterations))
    
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
  loglike_r = numeric(n.train.sh)
  
  for(i in 1:n.train.sh) {
    x_i = matXstar.sh[i, ]      # first row..{1,x,z}
    j = cl_membership.HB.train.sh[i] # membership ID
    
    loglike_r[i] = dlogsknorm_log(S.sh[i], x_i, betaparam.HB.st[j,], sig2param.HB.st[j], xiparam.HB.st[j]) + 
      dnorm(x=x_i[2], mean=kappaparam.HB.st[j,1] + kappaparam.HB.st[j,2]*Z.sh[i], sd=sqrt(lambda2param.HB.st[j]), log=TRUE) +
      dbinom(x=x_i[3], size=1, prob=piparam.HB.st[j], log=TRUE) 
  }  
  loglikelihood.HB.st[,r] = loglike_r
  
  if( (r > r_convergence) & (r %% save_every==0) ) {
    #>>> W/O CLUSTER
    list_beta0param.HB.st[[(r-r_convergence)/save_every]] = beta0param.HB.st
    list_SIG_b0param.HB.st[[(r-r_convergence)/save_every]] = SIG_b0param.HB.st
    list_u0param.HB.st[[(r-r_convergence)/save_every]] = u0param.HB.st
    list_v0param.HB.st[[(r-r_convergence)/save_every]] = v0param.HB.st
    
    list_betaparamfull.HB.st[[(r-r_convergence)/save_every]] = betaparam_full.HB.st
    list_sig2paramfull.HB.st[[(r-r_convergence)/save_every]] = sig2param_full.HB.st
    list_xiparamfull.HB.st[[(r-r_convergence)/save_every]] = xiparam_full.HB.st
    
    #>>> With Cluster-wise
    list_piparam.HB.st[[(r-r_convergence)/save_every]] = piparam.HB.st     
    list_lambda2param.HB.st[[(r-r_convergence)/save_every]] = lambda2param.HB.st
    list_kappaparam.HB.st[[(r-r_convergence)/save_every]] = kappaparam.HB.st
    list_tau2param.HB[[(r-r_convergence)/save_every]] = tau2param.HB
    
    list_betaparam.HB.st[[(r-r_convergence)/save_every]] = betaparam.HB.st  
    list_sig2param.HB.st[[(r-r_convergence)/save_every]] = sig2param.HB.st
    list_xiparam.HB.st[[(r-r_convergence)/save_every]] = xiparam.HB.st
    
    ################# STORE the CLEAN parameters ##################
    list_kappaparam.HB.clean[[(r-r_convergence)/save_every]] = kappaparamclean     #for Xstar
    list_lambda2param.HB.clean[[(r-r_convergence)/save_every]] = lambda2paramclean    #for Xstar
    list_betaparam.HB.clean[[(r-r_convergence)/save_every]] = betaparamclean    #for S
    list_sig2param.HB.clean[[(r-r_convergence)/save_every]] = sig2paramclean    #for S
    list_xiparam.HB.clean[[(r-r_convergence)/save_every]] = xiparamclean        #for S
    
    # list_xclean[[r-r_convergence]] = xclean #~xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx Q. (can we impute?)
    
  }
  print(paste("r=",r))
  #trace(what = "print", where = getNamespace("base"), exit = flush.console, print = FALSE)
}

plot(colSums(loglikelihood.HB.st),type="l")
lppd_EX_HB.S.st = mean( colSums(loglikelihood.HB.st) ); lppd_EX_HB.S.st

### ?
list_cl.HB.st <- list()  
for (i in 1:( (total_iter-r_convergence)/save_every ) ) {
  list_cl.HB.st[[i]] <- cl_membership.HB.train  # Assign the vector to each list element
}


# From Gibbs sampler we need cluster memberships and parameters - [ model(D): Gustafson HB ] - multiple[1%X3, 10%X3, 25%X3]
filename <- paste("testtt/HB.st.shard_0.01_0.99(",sh,")data.rds",sep="")
saveRDS(object = list(shard.sh.idx, anchor, list_cl.HB.st, list_piparam.HB.st, list_lambda2param.HB.st, 
                      list_kappaparam.HB.st, list_tau2param.HB,
                      list_betaparam.HB.st, list_sig2param.HB.st, list_xiparam.HB.st, 
                      list_kappaparam.HB.clean, list_lambda2param.HB.clean,
                      list_betaparam.HB.clean, list_sig2param.HB.clean, list_xiparam.HB.clean, N), file = filename)




