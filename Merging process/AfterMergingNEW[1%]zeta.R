######################## - PREDICTION and Comparison - #########################
# shard outputs are ready??
# then go.
Nsh <- 11

# It starts with checking clusters from each of the shard after the first shard, and comparing them with all the clusters 
#in the first shard to determine whether we can merge them (new cluster born). If they get merged, we use linear combination 
#of the parameters from each of the two clusters based upon the sample sizes. If they don't get merged, just append 
#the clusters to our list of clusters (with cl_membership).  



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##---------------------- [ model(A): Gold standard ] -------------------------##
###> Load shard output
sh_out_GS <- list()
for(j in 1:Nsh) {
  sh_out_GS[[j]] <- readRDS( paste("testtt/GS.shard(",j,")data.rds",sep="") )
}

###> Retrieve total iterations
total_iter_GS <- length(sh_out_GS[[1]][[3]]) # in the shard NO.1, third element: "list_cl"...
total_iter_GS

###> Retrieve total sample size
N <- sh_out_GS[[1]][[10]]

###> Retrieve alpha paramater ( not cl-wise but shard-wise )
alphaparam_GS_new <- numeric(total_iter_GS)
total_observations <- 0
# - linear combination: weighted AVG
for(j in 1:Nsh) {                                 #index in the shard...........list_alphaparam
  alphaparam_GS_new <- alphaparam_GS_new + length(sh_out_GS[[j]][[1]])*unlist(sh_out_GS[[j]][[6]])
  total_observations <- total_observations + length(sh_out_GS[[j]][[1]])
}
alphaparam_GS_new <- alphaparam_GS_new/total_observations; alphaparam_GS_new


###> Retrieve anchor points (in the 1st shard, the second item: )
anchor <- sh_out_GS[[1]][[2]]; anchor

###> Determine the threshold
epsilon <- 0.5


#################### CHECK the first 3 shards....R u ready?
# r=1 #iteration
# 
# j=1 #shard NO.
# sh_out_GS[[j]][[4]][[r]] # pi
# sh_out_GS[[j]][[5]][[r]] # lambda
# sh_out_GS[[j]][[7]][[r]] # beta
# sh_out_GS[[j]][[8]][[r]]  # sig2
# sh_out_GS[[j]][[9]][[r]] # xi
# 
# j=2
# sh_out_GS[[j]][[4]][[r]] # pi
# sh_out_GS[[j]][[5]][[r]] # lambda
# sh_out_GS[[j]][[7]][[r]] # beta
# sh_out_GS[[j]][[8]][[r]]  # sig2
# sh_out_GS[[j]][[9]][[r]] # xi
# 
# j=3
# sh_out_GS[[j]][[4]][[r]] # pi
# sh_out_GS[[j]][[5]][[r]] # lambda
# sh_out_GS[[j]][[7]][[r]] # beta
# sh_out_GS[[j]][[8]][[r]]  # sig2
# sh_out_GS[[j]][[9]][[r]] # xi
################################################################################



###> Create a matrix to store cl_membership: for r=1, blah; for r=2, blab; ..... Tabulation!!
###> For what? ::: for each observation, some r anchor, some r regular. we only count anchor points,..how to track down?
#::: Are any of anchor pt in the cluster show up in the [1st shard]? 


###> This is initial cl_membership of the 1st shard, tabulating on the entire sample matrix(n.train)
cl_membership_GS.1_matrix <- list()
for(r in 1:total_iter_GS) {
  cl_membership_GS.1_matrix[[r]] <- matrix( 0, nrow=N, ncol=max(sh_out_GS[[1]][[3]][[r]]) ) # no anchor pt...no 1st shard
  # for obv in the 1st shard, obtain cl_memberships....and point out the "anchor pt"!!!!!!!!
  for(i in 1:length(sh_out_GS[[1]][[1]])) {
    cl_membership_GS.1_matrix[[r]][ sh_out_GS[[1]][[1]][i], sh_out_GS[[1]][[3]][[r]][i] ] <- 1
  }
}
cl_membership_GS.1_matrix # this is r-lists of the [1st shard] cl.info matrix...including its anchor points. 
#            cl.01,      cl.02,     cl.03, ....
#obv1: anchor(Y/N) anchor(Y/N) anchor(Y/N) .....
#obv2: anchor(Y/N) anchor(Y/N) anchor(Y/N) .......
# ...      ...        ...         ...      ......... in the [1st shard]
#n.train

#:::: the cl in the [1st shard] will grow...the anchor pt will duplicates....

#total_iter_GS = 3

###> Next, 
cl_membership_GS.1_vector <- list()
loglikelihood = numeric(total_iter_GS)
# 
for (r in 1:total_iter_GS) {
  # -> pick up cl_membership (i.e. "which cl" the data point belong to..) for the 1st shard
  cl_membership_GS.1 <- cl_membership_GS.1_matrix[[r]]
  #::: every info in the 1st shard looks great, arranged on a single matrix...but 2nd, 3rd,...shards are not...so.. 

  # -> pick up cl_membership for the 2nd, 3rd, .... shard to compare
  for(j in 2:Nsh) {
    cl_membership_GS.2 <- sh_out_GS[[j]][[3]][[r]] 
    n_clust_GS.2 <- max(cl_membership_GS.2) # in shard[j] at iter[r], CHECK how many clusters?
    
    # for each cluster "k2", 
    for(k2 in 1:n_clust_GS.2) {
      ##### 1) challenger shard info
      indx_sh2_clk2 <- sh_out_GS[[j]][[1]][cl_membership_GS.2==k2] # collect indexes that belong to the cluster "k2"
      anchor_sh2_clk2 <- anchor[anchor %in% indx_sh2_clk2]         # collect the [anchor points] in the the cluster "k2"
      
      
      #:::take these indexes and compare them with every indexes in the 1st shard!!!!
      n_clust_GS.1 <- ncol(cl_membership_GS.1) # how many cl in the 1st shard?
      k1=1
      MERGED=FALSE
      # -> If..Merged:
      while ( k1<=n_clust_GS.1 & !MERGED ) {
        ##### 2) Champion shard info
        indx_sh1_clk1 <- which(cl_membership_GS.1[, k1]==1)   # collect indexes that belong to the cluster "k1" from the matrix
        anchor_sh1_clk1 <- anchor[anchor %in% indx_sh1_clk1] # collect the [anchor points] in the 1st shard (k1)
        
        #> compute:
        C_dist <- sum(anchor_sh1_clk1 %in% anchor_sh2_clk2)                     # common elements COUNT
        D_dist <- length(anchor_sh1_clk1) + length(anchor_sh2_clk2) - 2*C_dist  # different elements COUNT
        
        d_distance <- ifelse(C_dist > 0 | D_dist > 0, D_dist/(C_dist+D_dist), 1)
        #d_distance <- D_dist/(C_dist+D_dist) 
        # Indeed, there is no guarantee that anchor points will present every cluster in each shard since the number of 
        # points is limited. Some clusters might not have any anchor points.
        
        #print(paste("dist:"), d_distance)
        print(paste("dist:",d_distance))
        
        # FINALLYYYYYY COMPUTE: weight AVG of paramssssss
        if(d_distance < epsilon) {  # Merge clusters
          cl_membership_GS.1[indx_sh2_clk2, k1] <- 1
          sh_out_GS[[1]][[4]][[r]][k1] <- length(indx_sh1_clk1)/(length(indx_sh1_clk1) + length(indx_sh2_clk2))*sh_out_GS[[1]][[4]][[r]][k1] +
            length(indx_sh2_clk2)/(length(indx_sh1_clk1) + length(indx_sh2_clk2))*sh_out_GS[[j]][[4]][[r]][k2] # update for piparameter
          sh_out_GS[[1]][[5]][[r]][k1] <- length(indx_sh1_clk1)/(length(indx_sh1_clk1) + length(indx_sh2_clk2))*sh_out_GS[[1]][[5]][[r]][k1] +
            length(indx_sh2_clk2)/(length(indx_sh1_clk1) + length(indx_sh2_clk2))*sh_out_GS[[j]][[5]][[r]][k2] # update for lambda2parameter
          sh_out_GS[[1]][[7]][[r]][k1,] <- length(indx_sh1_clk1)/(length(indx_sh1_clk1) + length(indx_sh2_clk2))*sh_out_GS[[1]][[7]][[r]][k1,] +
            length(indx_sh2_clk2)/(length(indx_sh1_clk1) + length(indx_sh2_clk2))*sh_out_GS[[j]][[7]][[r]][k2,] # update for betaparameter
          sh_out_GS[[1]][[8]][[r]][k1] <- length(indx_sh1_clk1)/(length(indx_sh1_clk1) + length(indx_sh2_clk2))*sh_out_GS[[1]][[8]][[r]][k1] +
            length(indx_sh2_clk2)/(length(indx_sh1_clk1) + length(indx_sh2_clk2))*sh_out_GS[[j]][[8]][[r]][k2] # update for sig2parameter
          sh_out_GS[[1]][[9]][[r]][k1] <- length(indx_sh1_clk1)/(length(indx_sh1_clk1) + length(indx_sh2_clk2))*sh_out_GS[[1]][[9]][[r]][k1] +
            length(indx_sh2_clk2)/(length(indx_sh1_clk1) + length(indx_sh2_clk2))*sh_out_GS[[j]][[9]][[r]][k2] # update for xiparameter
          MERGED = TRUE
        }
        k1=k1+1
      } #END of [while] 
      
      # -> IF.. not MERGED: 
      if (!MERGED){ 
        cl_membership_GS.1 <- cbind(cl_membership_GS.1, 0)
        cl_membership_GS.1[indx_sh2_clk2, ncol(cl_membership_GS.1)] <- 1
        sh_out_GS[[1]][[4]][[r]] <- c(sh_out_GS[[1]][[4]][[r]],sh_out_GS[[j]][[4]][[r]][k2]) # pi ::: append
        sh_out_GS[[1]][[5]][[r]] <- c(sh_out_GS[[1]][[5]][[r]],sh_out_GS[[j]][[5]][[r]][k2]) # lambda2 ::: append
        sh_out_GS[[1]][[7]][[r]] <- rbind(sh_out_GS[[1]][[7]][[r]],sh_out_GS[[j]][[7]][[r]][k2,]) # beta ::: append
        sh_out_GS[[1]][[8]][[r]] <- c(sh_out_GS[[1]][[8]][[r]],sh_out_GS[[j]][[8]][[r]][k2]) # sig2 ::: append
        sh_out_GS[[1]][[9]][[r]] <- c(sh_out_GS[[1]][[9]][[r]],sh_out_GS[[j]][[9]][[r]][k2]) # xi ::: append
      }
    }
  }
  #any( apply(cl_membership_GS.1,1,sum)==0 ) # all obv is assigned with cl?
  #summary( apply(cl_membership_GS.1[anchor,],1, sum) )
  #which(apply(cl_membership_GS.1[anchor,],1,sum)==1)
  #anchor[which(apply(cl_membership_GS.1[anchor,],1,sum)==1)]
  
  #a <- anchor.test[1]
  
  ###> Remove duplicates
  for (a in anchor) {
    X_a <- cl_membership_GS.1[a, ]
    
    # Be careful............
    if(sum(X_a)>1) {
      which.is.1 <- which(X_a==1) # vector
      keep1 <- sample(which.is.1, 1)
      cl_membership_GS.1[a, -keep1] <- 0
    } # since some anchor pt has its own cluster
  }
  cl_membership_GS.1_vector[[r]] <- unlist( apply(cl_membership_GS.1, 1, function(x) which(x==1)) ) # vector
  cl_membership_GS.1_matrix[[r]] <- cl_membership_GS.1
  
  ###> Compute LPPD
  #loglike_r = 0
  
  loglike_r_vector = numeric(N)
  
  
  for(i in 1:N) {
    x_i = matX_tru[i, ]
    jj = cl_membership_GS.1_vector[[r]][i]
    
    #loglike_r = loglike_r +  
    loglike_r_vector[i] = 
      dlogsknorm_log(s=S[i], x=x_i, beta=sh_out_GS[[1]][[7]][[r]][jj, ], sig2=sh_out_GS[[1]][[8]][[r]][jj], xi=sh_out_GS[[1]][[9]][[r]][jj]) +
      dbinom(x = x_i[3], size = 1, prob = sh_out_GS[[1]][[4]][[r]][jj], log = TRUE) + 
      dnorm(x = x_i[2], mean = mean(X_tru), sd = sqrt(sh_out_GS[[1]][[5]][[r]][jj]), log = TRUE)
  }
  #loglikelihood[r] = loglike_r # cumulative sum
  loglikelihood[r] = sum(loglike_r_vector) # cumulative sum
}

### De-bugging.................................................
# cl_membership_GS.1_vector[[3]]
# length( cl_membership_GS.1_vector[[3]] )
# length( cl_membership_GS.1_vector[[2]] )
# 
# summary( loglike_r_vector )
# which( is.na( loglike_r_vector  ) )
# summary(apply(cl_membership_GS.1_matrix[[3]],1,sum))
# which(apply(cl_membership_GS.1_matrix[[3]],1,sum)==0)
# which(apply(cl_membership_GS.1_matrix[[3]],1,sum)==0) %in% anchor
# 
# anchor.test <- anchor[which(apply(cl_membership_GS.1[anchor,],1,sum)==1)]
# cl_membership_GS.1[anchor.test,]
# 
# sh_out_GS[[1]][[3]][[3]][sh_out_GS[[1]][[1]] %in% anchor.test]
# sh_out_GS[[2]][[3]][[3]][sh_out_GS[[2]][[1]] %in% anchor.test]
# table( sh_out_GS[[1]][[3]][[3]][sh_out_GS[[1]][[1]] %in% anchor.test], sh_out_GS[[2]][[3]][[3]][sh_out_GS[[2]][[1]] %in% anchor.test] )
# 
# table( sh_out_GS[[1]][[3]][[3]][sh_out_GS[[1]][[1]] %in% anchor.test], 
#        sh_out_GS[[3]][[3]][[3]][sh_out_GS[[3]][[1]] %in% anchor.test] )


# all duplicates are gone?
table( apply(cl_membership_GS.1[anchor,],1,sum) )

# # Make sure no missing obv in the clustering?
# table( apply(cl_membership_GS.1, 1, sum) )
# # r=1
# table( cl_membership_GS.1_vector[[1]] )
# # r=2
# table( cl_membership_GS.1_vector[[2]] )
# 
# r=1 #first iteration
# #j=1
# sh_out_GS[[1]][[4]][[r]] # pi
# sh_out_GS[[1]][[5]][[r]] # lambda
# sh_out_GS[[1]][[7]][[r]] # beta
# sh_out_GS[[1]][[8]][[r]]  # sig2
# sh_out_GS[[1]][[9]][[r]] # xi
# alphaparam_GS_new
# cl_membership_GS.1_vector[[r]]
# 
# r=2 #second iteration
# #j=1
# sh_out_GS[[1]][[4]][[r]] # pi
# sh_out_GS[[1]][[5]][[r]] # lambda
# sh_out_GS[[1]][[7]][[r]] # beta
# sh_out_GS[[1]][[8]][[r]]  # sig2
# sh_out_GS[[1]][[9]][[r]] # xi
# alphaparam_GS_new
# cl_membership_GS.1_vector[[r]]


###> Now it's HARVEST TIME! ----------------------------------------------------

###> CONNECT to prediction!!!
list_piparam = sh_out_GS[[1]][[4]]
list_lambda2param = sh_out_GS[[1]][[5]]
list_alphaparam = alphaparam_GS_new

list_betaparam = sh_out_GS[[1]][[7]]
list_sig2param = sh_out_GS[[1]][[8]]
list_xiparam = sh_out_GS[[1]][[9]]

list_cl = cl_membership_GS.1_vector

n.train = N

options(scipen = 999)
plot(loglikelihood, type="l")
lppd_EX_DP.S <- mean( loglikelihood ); lppd_EX_DP.S # -656377.4


#..............................................................................#
#..............................................................................#
############################### Predictions ####################################
#..............................................................................#
#..............................................................................#

################################################################################
###> Training data -------------------------------------------------------------
################################################################################
#% Ultimate liquid weight vector for each iteration
# W_paramFree.vec = matrix(0, nrow = n.train, ncol = total_iter_GS) 
# density.train = matrix(0, nrow = n.train, ncol = total_iter_GS)
# expval.train = matrix(0, nrow = n.train, ncol = total_iter_GS)
# 
# #% prediction for each iteration
# for(r in 1:total_iter_GS) {
#   piparamclean = list_piparam[[r]]
#   lambda2paramclean = list_lambda2param[[r]] #for Xstar
#   alphaparamclean = list_alphaparam[r]      #for alpha
#   
#   betaparamclean = list_betaparam[[r]]      #for S
#   sig2paramclean = list_sig2param[[r]]      #for S
#   xiparamclean = list_xiparam[[r]]       #for S
#   
#   cl_membership.train = list_cl[[r]] # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   
#   J = nrow(betaparamclean)
#   f_S = numeric(J)                  # density for S
#   
#   #% Ultimate solid weight matrix
#   W_paramBase.matrix = matrix(0, nrow=n.train, ncol=J)
#   
#   for(h in 1:n.train) {
#     
#     fX_j = numeric(J)       # discrete covariate model for each j
#     w_j.solid = numeric(J)  # discrete weight for each j  
#     E_value = numeric(J)    # to save E[S|X] for each j
#     
#     for(j in 1:J) {
#       # discrete covariate joint
#       fX_j[j] = dbinom(x=Z[h], size=1, prob=piparamclean[j])*
#         dnorm(x=X_tru[h], mean=mean(X_tru), sd=sqrt(lambda2paramclean[j]))
#       
#       #% solid weight ***** component
#       w_j.solid[j] = length(S[cl_membership.train==j])/(alphaparamclean + n.train)*fX_j[j]
#       
#       #% compute pred.density
#       f_S[j] = dlogsknorm(s=S[h], x=c(1, X_tru[h], Z[h]), 
#                           beta=betaparamclean[j,], sig2=sig2paramclean[j], xi=xiparamclean[j])
#       #% compute pred.value
#       E_value[j] = 4*exp(sum( c(1, X_tru[h], Z[h])*betaparamclean[j,]) - sig2paramclean[j]/2)*
#         (1-pnorm(-xiparamclean[j]*sqrt(sig2paramclean[j])/sqrt(xiparamclean[j]^2+1)))
#     }
#     
#     #% liquid weight ***** component
#     WJ1 = alphaparamclean/(alphaparamclean + n.train)*f0X[h] #; print(paste("liquid weight=", WJ1))
#     
#     #% Collecting values and constructing weights for each iteration
#     W_paramFree.vec[h,r] = WJ1/(WJ1 + sum(w_j.solid))
#     W_paramBase.matrix[h, ] = w_j.solid/(WJ1 + sum(w_j.solid)) 
#     
#     #% Collecting values for predictive density (via weighted AVG)
#     density.train[h,r] = W_paramFree.vec[h,r]*f0S[h] + sum(W_paramBase.matrix[h,]*f_S)
#     #% Collecting values for prediction (via weighted AVG)
#     expval.train[h,r] = W_paramFree.vec[h,r]*E0S[h] + sum(W_paramBase.matrix[h,]*E_value)
#   }
#   print(paste("r=", r))
# }
# 
# # AVG prediction
# expS.DP.Avg.train <- apply(X=expval.train, MARGIN=1, FUN=mean)
# 
# summary(S)
# summary(expS.DP.Avg.train)
# summary(log(S))
# summary(log(expS.DP.Avg.train))

################################################################################
###> Testing data -------------------------------------------------------------
################################################################################ 
###> Parameter-free component
#> Analytical form
# f0x1 = function(Z) {
#   beta(Z+g0, 1-Z+h0)/beta(g0,h0)
# }
# #> Analytical form
# f0x2 = function(X, covar) {
#   1/sqrt(pi)*gamma(c0+1)/gamma(c0)*d0^(c0/2)/( (X-mean(covar))^2 + d0 )^((c0+1)/2)
# }
# #> JOINT 
# f0X.test = numeric(n.test) #% param-free covariate model f0(X,Z)
# 
# ##### For [Param-free outcome model]: -------------------> MonteCarlo Integration
# # Calculate Outcome and Covariate parameter free data model for each observation  
# f0S.test = numeric(n.test) #% param-free outcome model f0(S|X,Z)
# E0S.test = numeric(n.test) #% param-free E(S|X,Z) = Expected value of S|X,Z
# 
# M = 1000 # Number of Monte Carlo samples
# sumS = numeric(M)
# sumES = numeric(M)
# 
# for(h in 1:n.test) {
#   
#   ##> Analytical solution for f0(X, Z) for ### covariate model ###
#   # : necessary parameters are supplied by pre-determined hyperparam values...
#   f0X.test[h] = f0x1(Z.test[h])*f0x2(X_tru.test[h], X_tru.test)  
#   
#   ##> Monte Carlo Integration for the ### outcome model ###
#   for(j in 1:M) {
#     xi_samplej = rt(n=1, df=nu0)                               # bullet on xi
#     sig_samplej = rinvgamma(n=1, shape=u0, scale=v0)           # bullet on sig2
#     beta_samplej = rmvn(n=1, mu=m0, sigma=SIG_b0*varinf)       # bullet on beta
#     
#     sumS[j] = dlogsknorm( s=S.test[h], x=matX_tru.test[h,], beta=beta_samplej, sig2=sig_samplej, xi=xi_samplej )* 
#       dmvn(X=beta_samplej, mu=m0, sigma=SIG_b0*varinf)*        # to joint beta
#       dinvgamma(x=sig_samplej, shape=u0, scale=v0)*            # to joint sig2
#       dt(x=xi_samplej, df=nu0)                                 # to joint xi
#     # 20*exp( betaparamclean[j,1] + betaparamclean[j,2]*cleanx + betaparamclean[j,3]*Z[i] - sig2paramclean[j]/2)*( 1-pnorm(-xiparamclean[j]*sqrt(sig2paramclean[j])/sqrt(xiparamclean[j]^2+1)))
#     sumES[j] = 12*exp( sum(matX_tru.test[h,]*beta_samplej) - sig_samplej/2 )*(1-pnorm(-xi_samplej*sqrt(sig_samplej)/sqrt(xi_samplej^2+1)))
#   }
#   E0S.test[h] = sum(sumES)/M  # E[S] value of f0
#   f0S.test[h] = sum(sumS)/M   # density of f0
#   
#   print(paste("h=",h))
# }
# summary(S.test)    # original S
# summary(E0S.test)  # E[S] based on f0(S|x) by train data

###> NOW....
# cl_membership.test <- full.test.sweden$Risk
# table(cl_membership.test)                    # this is wrong!!!!!!!!!!!!!!!!!

#% Ultimate liquid weight vector for each iteration
W_paramFree.vec.test = matrix(0, nrow = n.test, ncol = total_iter_GS) 
density.test = matrix(0, nrow = n.test, ncol = total_iter_GS)
expval.test = matrix(0, nrow = n.test, ncol = total_iter_GS)


for(r in 1:(total_iter_GS)) {
  piparamclean = list_piparam[[r]]
  lambda2paramclean = list_lambda2param[[r]] #for Xstar
  alphaparamclean = list_alphaparam[r]      #for alpha
  
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
      w_j.solid[j] = sum(cl_membership.train==j) / (alphaparamclean + n.train)*fX_j[j] #~~~~~~~~~~~~~~~~~~~~~~~~~~~~QQQ         
      #% compute pred.density
      f_S.test[j] = dlogsknorm(s=S.test[h], x=c(1, X_tru.test[h], Z.test[h]), 
                               beta=betaparamclean[j,], sig2=sig2paramclean[j], xi=xiparamclean[j])
      #% compute pred.value (2 or 12?)
      E_value.test[j] = 2*exp(sum( c(1, X_tru.test[h], Z.test[h])*betaparamclean[j,]) - sig2paramclean[j]/2)*
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


####> [AVG.vec Scenario] 
den.DP.Avg.test <- apply(X=density.test, MARGIN=1, FUN=mean)
expS.DP.Avg.test <- apply(X=expval.test, MARGIN=1, FUN=mean)

summary(S.test)
summary(expS.DP.Avg.test)
summary(log(S.test))
summary(log(expS.DP.Avg.test))

SSPE.DP.LSN <- sum( (log(expS.DP.Avg.test) - log(S.test))^2 ); SSPE.DP.LSN    #12420.39
SAPE.DP.LSN <- sum( abs(log(expS.DP.Avg.test) - log(S.test)) ); SAPE.DP.LSN   #4285.839
lppd_EX_DP.S                                                                  #-656377.4






###----------------------------- Visual S -----------------------------------###

expval.test
density.test

par(mfrow=c(1,1))
plot(density(expval.test[,1])) # pointless
plot(density(expval.test[,100])) # pointless

log_expS.DP.Avg.test = log( apply(X=expval.test, MARGIN=1, FUN=mean) )

### turn a matrix into a dataframe of 100 scenarios of E[Y|X] and their densities
den_df <- data.frame(density.test); head(den_df,2) 
val_df <- data.frame(expval.test); head(val_df,2) 

summary(val_df[, 1])    # 1st scenario 
summary(val_df[, 100]) # 100th scenario
#summary(val_df[, 1950:2000])

# hist(log(S), breaks=n.breaks.train, freq=F, xlab ="log(S)", cex.axis=2, cex.lab=2, col="white", ylim=c(0,0.7))
# for (i in 51:100) {
#   lines(density( log(val_df[, i])+1.2, bw = 0.15, kernel = "epanechnikov" ), lwd=0.01, col="black" )
# }
# lines(density(log_expS.DP.Avg.test+1.2, bw = 0.15, kernel = "epanechnikov"), col="red", lwd=3)

library(ks)
# i=1
# data <- val_df[, i]
# weights <- ifelse(data > quantile(data, 0.75), 1.2, 1)  # Increase weight in upper tail
# plot( kde( x=log(data)+1.2, w=weights) )
hist(log(S.test), breaks=n.breaks.test, freq=F, xlab ="log(S)", cex.axis=2, cex.lab=2, col="white", ylim=c(0,0.7))

for (i in 1:50) {
  data <- val_df[, i]
  weights <- ifelse(data > quantile(data, 0.975), 3, 1)
  kde.res <- kde( x=log(data)+1.3, w=weights, H=diag(0.1, ncol = 1) )
  lines(kde.res$eval.points, kde.res$estimate, lwd=0.01, col="black" )
}

for (i in 51:100) {
  data <- val_df[, i]
  weights <- ifelse(data > quantile(data, 0.975), 3, 1)
  kde.res <- kde( x=log(data)+1.3, w=weights, H=diag(0.1, ncol = 1) )
  lines(kde.res$eval.points, kde.res$estimate, lwd=0.01, col="black" )
}

data <- expS.DP.Avg.test
kde.res <- kde( x=log(data)+1.3, w=weights, H=diag(0.1, ncol = 1) )
lines(kde.res$eval.points, kde.res$estimate, col="red", lwd=3)






##### ----------------------- 95% credible interval ---------------------- #####
compute_CI <- function(samples, alpha = 0.05) {
  quantile(samples, probs = c(alpha / 2, 1 - alpha / 2))
}


#> beta0
list_betaparam0 <- numeric(100)
for (i in 1:100){
  list_betaparam0[i] <- mean( list_betaparam[[i]][, 1] ) 
}
summary(list_betaparam0)# 10.17 

compute_CI(list_betaparam0, alpha = 0.05)[1]
#    2.5% 
# 10.16147 
compute_CI(list_betaparam0, alpha = 0.05)[2]
#    97.5% 
# 10.18141


#> beta1
list_betaparam1 <- numeric(100)
for (i in 1:100){
  list_betaparam1[i] <- mean( list_betaparam[[i]][, 2] ) 
}
summary(list_betaparam1)# 0.3981    

compute_CI(list_betaparam1, alpha = 0.05)[1]
#    2.5% 
# 0.3963701  
compute_CI(list_betaparam1, alpha = 0.05)[2]
#    97.5% 
# 0.4004199  


#> beta2
list_betaparam2 <- numeric(100)
for (i in 1:100){
  list_betaparam2[i] <- mean( list_betaparam[[i]][, 3] ) 
}
summary(list_betaparam2) #-1.954   

compute_CI(list_betaparam2, alpha = 0.05)[1]
#    2.5% 
# -1.964814  
compute_CI(list_betaparam2, alpha = 0.05)[2]
#    97.5% 
# -1.944818   


#> sig2
list_sig2test <- numeric(100)
for (i in 1:100) {
  list_sig2test[i] <- mean( list_sig2param[[i]] )
}
summary(list_sig2test ) #11.84      

compute_CI(list_sig2test, alpha = 0.05)[1]
#    2.5% 
# 11.0204    
compute_CI(list_sig2test, alpha = 0.05)[2]
#    97.5% 
# 12.96817   


#> xi
list_xitest <- numeric(100)
for (i in 1:100) {
  list_xitest[i] <- mean( list_xiparam[[i]] )
}
summary(list_xitest ) #-1.110   

compute_CI(list_xitest, alpha = 0.05)[1]
#    2.5% 
# -1.284106    
compute_CI(list_xitest, alpha = 0.05)[2]
#    97.5% 
# -0.9691677   














#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##------------- [ model(B): Model risk ] - just three[1%, 10%, 25%] ----------##
###> Load shard output
sh_out_xx <- list()
for(j in 1:Nsh) {
  sh_out_xx[[j]] <- readRDS( paste("testtt/xx.shard_0.01(",j,")data.rds",sep="") )
}

###> Retrieve total iterations
total_iter_xx <- length(sh_out_xx[[1]][[3]]) # in the shard NO.1, third element: "list_cl"...

###> Retrieve total sample size
N <- sh_out_xx[[1]][[10]]

###> Retrieve alpha paramater ( not cl-wise but shard-wise )
alphaparam_xx_new <- numeric(total_iter_xx)
total_observations.xx <- 0
# - linear combination: weighted AVG
for(j in 1:Nsh) {                                 #index in the shard...........list_alphaparam
  alphaparam_xx_new <- alphaparam_xx_new + length(sh_out_xx[[j]][[1]])*unlist(sh_out_xx[[j]][[6]])
  total_observations.xx <- total_observations.xx + length(sh_out_xx[[j]][[1]])
}
alphaparam_xx_new <- alphaparam_xx_new/total_observations.xx
alphaparam_xx_new

###> Retrieve anchor points (in the 1st shard, the second item: )
anchor.xx <- sh_out_xx[[1]][[2]]; anchor.xx

###> Determine the threshold
epsilon <- 0.5 #0.75







###> This is initial cl_membership of the 1st shard, tabulating on the entire sample matrix(n.train)
cl_membership_xx.1_matrix <- list()
for(r in 1:total_iter_xx) {
  cl_membership_xx.1_matrix[[r]] <- matrix( 0, nrow=N, ncol=max(sh_out_xx[[1]][[3]][[r]]) ) # no anchor pt...no 1st shard
  # for obv in the 1st shard, obtain cl_memberships....and point out the "anchor pt"!!!!!!!!
  for(i in 1:length(sh_out_xx[[1]][[1]])) {
    cl_membership_xx.1_matrix[[r]][ sh_out_xx[[1]][[1]][i], sh_out_xx[[1]][[3]][[r]][i] ] <- 1
  }
}
cl_membership_xx.1_matrix # this is the [1st shard] including its anchor points. 
#            cl.01,      cl.02,     cl.03, ....
#obv1: anchor(Y/N) anchor(Y/N) anchor(Y/N) .....
#obv2: anchor(Y/N) anchor(Y/N) anchor(Y/N) .......
# ...      ...        ...         ...      ......... in the [1st shard]
#n.train

#:::: the cl in the [1st shard] will grow...the anchor pt will duplicates....


###> Next, 
cl_membership_xx.1_vector <- list()
loglikelihood.xx = numeric(total_iter_xx)
# 
for (r in 1:total_iter_xx) {
  # -> pick up cl_membership (i.e. "which cl" the data point belong to..) for the 1st shard
  cl_membership_xx.1 <- cl_membership_xx.1_matrix[[r]]
  #::: every info in the 1st shard looks great, arranged on a single matrix...but 2nd, 3rd,...shards are not...so.. 
  
  # -> pick up cl_membership for the 2nd, 3rd, .... shard to compare
  for(j in 2:Nsh) {
    cl_membership_xx.2 <- sh_out_xx[[j]][[3]][[r]] 
    n_clust_xx.2 <- max(cl_membership_xx.2) # in shard[j] at iter[r], CHECK how many clusters?
    
    # for each cluster "k2", 
    for(k2 in 1:n_clust_xx.2) {
      ##### 1) challenger shard info
      indx_sh2_clk2 <- sh_out_xx[[j]][[1]][cl_membership_xx.2==k2] # collect indexes that belong to the cluster "k2"
      anchor_sh2_clk2 <- anchor.xx[anchor.xx %in% indx_sh2_clk2]         # collect the [anchor points] in the the cluster "k2"
      
      
      #:::take these indexes and compare them with every indexes in the 1st shard!!!!
      n_clust_xx.1 <- ncol(cl_membership_xx.1) # how many cl in the 1st shard?
      k1=1
      MERGED=FALSE
      # -> If..Merged:
      while ( k1<=n_clust_xx.1 & !MERGED ) {
        ##### 2) Champion shard info
        indx_sh1_clk1 <- which(cl_membership_xx.1[, k1]==1)   # collect indexes that belong to the cluster "k1" from the matrix
        anchor_sh1_clk1 <- anchor.xx[anchor.xx %in% indx_sh1_clk1] # collect the [anchor points] in the 1st shard (k1)
        
        #> compute:
        C_dist <- sum(anchor_sh1_clk1 %in% anchor_sh2_clk2)                     # common elements COUNT
        D_dist <- length(anchor_sh1_clk1) + length(anchor_sh2_clk2) - 2*C_dist  # different elements COUNT
        
        d_distance <- ifelse(C_dist > 0 | D_dist > 0, D_dist/(C_dist+D_dist), 1)
        #d_distance <- D_dist/(C_dist+D_dist)
        
        #print(paste("dist:"), d_distance)
        print(paste("dist:",d_distance))
        
        # FINALLYYYYYY COMPUTE: weight AVG of paramssssss
        if(d_distance < epsilon) {  # Merge clusters
          cl_membership_xx.1[indx_sh2_clk2, k1] <- 1
          sh_out_xx[[1]][[4]][[r]][k1] <- length(indx_sh1_clk1)/(length(indx_sh1_clk1) + length(indx_sh2_clk2))*sh_out_xx[[1]][[4]][[r]][k1] +
            length(indx_sh2_clk2)/(length(indx_sh1_clk1) + length(indx_sh2_clk2))*sh_out_xx[[j]][[4]][[r]][k2] # update for piparameter
          sh_out_xx[[1]][[5]][[r]][k1] <- length(indx_sh1_clk1)/(length(indx_sh1_clk1) + length(indx_sh2_clk2))*sh_out_xx[[1]][[5]][[r]][k1] +
            length(indx_sh2_clk2)/(length(indx_sh1_clk1) + length(indx_sh2_clk2))*sh_out_xx[[j]][[5]][[r]][k2] # update for lambda2parameter
          sh_out_xx[[1]][[7]][[r]][k1,] <- length(indx_sh1_clk1)/(length(indx_sh1_clk1) + length(indx_sh2_clk2))*sh_out_xx[[1]][[7]][[r]][k1,] +
            length(indx_sh2_clk2)/(length(indx_sh1_clk1) + length(indx_sh2_clk2))*sh_out_xx[[j]][[7]][[r]][k2,] # update for betaparameter
          sh_out_xx[[1]][[8]][[r]][k1] <- length(indx_sh1_clk1)/(length(indx_sh1_clk1) + length(indx_sh2_clk2))*sh_out_xx[[1]][[8]][[r]][k1] +
            length(indx_sh2_clk2)/(length(indx_sh1_clk1) + length(indx_sh2_clk2))*sh_out_xx[[j]][[8]][[r]][k2] # update for sig2parameter
          sh_out_xx[[1]][[9]][[r]][k1] <- length(indx_sh1_clk1)/(length(indx_sh1_clk1) + length(indx_sh2_clk2))*sh_out_xx[[1]][[9]][[r]][k1] +
            length(indx_sh2_clk2)/(length(indx_sh1_clk1) + length(indx_sh2_clk2))*sh_out_xx[[j]][[9]][[r]][k2] # update for xiparameter
          MERGED = TRUE
        }
        k1=k1+1
      } #END of [while] 
      
      # -> IF.. not MERGED: 
      if (!MERGED){ 
        cl_membership_xx.1 <- cbind(cl_membership_xx.1, 0)
        cl_membership_xx.1[indx_sh2_clk2, ncol(cl_membership_xx.1)] <- 1
        sh_out_xx[[1]][[4]][[r]] <- c(sh_out_xx[[1]][[4]][[r]],sh_out_xx[[j]][[4]][[r]][k2]) # pi ::: append
        sh_out_xx[[1]][[5]][[r]] <- c(sh_out_xx[[1]][[5]][[r]],sh_out_xx[[j]][[5]][[r]][k2]) # lambda2 ::: append
        sh_out_xx[[1]][[7]][[r]] <- rbind(sh_out_xx[[1]][[7]][[r]],sh_out_xx[[j]][[7]][[r]][k2,]) # beta ::: append
        sh_out_xx[[1]][[8]][[r]] <- c(sh_out_xx[[1]][[8]][[r]],sh_out_xx[[j]][[8]][[r]][k2]) # sig2 ::: append
        sh_out_xx[[1]][[9]][[r]] <- c(sh_out_xx[[1]][[9]][[r]],sh_out_xx[[j]][[9]][[r]][k2]) # xi ::: append
      }
    }
  }
  
  ###> Remove duplicates
  for (a in anchor.xx) {
    X_a <- cl_membership_xx.1[a, ]
    
    # Be careful............
    if(sum(X_a)>1) {
      which.is.1 <- which(X_a==1) # vector
      keep1 <- sample(which.is.1, 1)
      cl_membership_xx.1[a, -keep1] <- 0
    } # since some anchor pt has its own cluster
  }
  cl_membership_xx.1_vector[[r]] <- apply(cl_membership_xx.1, 1, function(x) which(x==1))
  cl_membership_xx.1_matrix[[r]] <- cl_membership_xx.1
  
  ###> Compute LPPD
  #loglike_r = 0
  
  loglike_r_vector = numeric(N)
  
  for(i in 1:N) {
    x_i = matXstar[i, ]
    jj = cl_membership_xx.1_vector[[r]][i]
    
    #loglike_r = loglike_r + 
    loglike_r_vector[i] =
      dlogsknorm_log(s=S[i], x=x_i, beta=sh_out_xx[[1]][[7]][[r]][jj, ], sig2=sh_out_xx[[1]][[8]][[r]][jj], xi=sh_out_xx[[1]][[9]][[r]][jj]) +
      dbinom(x = x_i[3], size = 1, prob = sh_out_xx[[1]][[4]][[r]][jj], log = TRUE) + 
      dnorm(x = x_i[2], mean = mean(Xstar), sd = sqrt(sh_out_xx[[1]][[5]][[r]][jj]), log = TRUE)
  }
  #loglikelihood.xx[r] = loglike_r # cumulative sum
  loglikelihood.xx[r] = sum(loglike_r_vector) # cumulative sum
}



###> Now it's HARVEST TIME! ----------------------------------------------------
###> 
###> CONNECT to prediction!!!
list_piparam.xx = sh_out_xx[[1]][[4]]
list_lambda2param.xx = sh_out_xx[[1]][[5]]
list_alphaparam.xx = alphaparam_xx_new

list_betaparam.xx = sh_out_xx[[1]][[7]]
list_sig2param.xx = sh_out_xx[[1]][[8]]
list_xiparam.xx = sh_out_xx[[1]][[9]]

list_cl.xx = cl_membership_xx.1_vector

n.train = N

options(scipen = 999)
plot( loglikelihood.xx, type="l")
lppd_EX_DP.S.xx <- mean( loglikelihood.xx ); lppd_EX_DP.S.xx #


#..............................................................................#
#..............................................................................#
############################### Predictions ####################################
#..............................................................................#
#..............................................................................#

################################################################################
###> Training data -------------------------------------------------------------
################################################################################
#% Ultimate liquid weight vector for each iteration
# W_paramFree.vec.xx = matrix(0, nrow = n.train, ncol = total_iter_xx) 
# 
# 
# density.train.xx = matrix(0, nrow = n.train, ncol = total_iter_xx)
# expval.train.xx = matrix(0, nrow = n.train, ncol = total_iter_xx)
# 
# for(r in 1:total_iter_xx) {
#   piparamclean.xx = list_piparam.xx[[r]]
#   lambda2paramclean.xx = list_lambda2param.xx[[r]] #for Xstar
#   alphaparamclean.xx = list_alphaparam.xx[[r]]      #for alpha
#   
#   betaparamclean.xx = list_betaparam.xx[[r]]      #for S
#   sig2paramclean.xx = list_sig2param.xx[[r]]      #for S
#   xiparamclean.xx = list_xiparam.xx[[r]]       #for S
#   
#   cl_membership.train.xx = list_cl.xx[[r]] # 
#   
#   J = nrow(betaparamclean.xx)
#   f_S.xx = numeric(J)                  # density for S
#   
#   #% Ultimate solid weight matrix
#   W_paramBase.matrix.xx = matrix(0, nrow=n.train, ncol=J)
#   
#   for(h in 1:n.train) {
#     
#     fX_j.xx = numeric(J)       # discrete covariate model for each j
#     w_j.solid.xx = numeric(J)  # discrete weight for each j  
#     E_value.xx = numeric(J)    # to save E[S|X] for each j
#     
#     for(j in 1:J) {
#       # discrete covariate joint
#       fX_j.xx[j] = dbinom(x=Z[h], size=1, prob=piparamclean.xx[j])*
#         dnorm(x=Xstar[h], mean=mean(Xstar), sd=sqrt(lambda2paramclean.xx[j]))
#       
#       #% solid weight ***** component
#       w_j.solid.xx[j] = length(S[cl_membership.train.xx==j])/(alphaparamclean.xx + n.train)*fX_j.xx[j]
#       
#       #% compute pred.density
#       f_S.xx[j] = dlogsknorm(s=S[h], x=c(1, Xstar[h], Z[h]), 
#                              beta=betaparamclean.xx[j,], sig2=sig2paramclean.xx[j], xi=xiparamclean.xx[j])
#       #% compute pred.value
#       E_value.xx[j] = 4*exp(sum( c(1, Xstar[h], Z[h])*betaparamclean.xx[j,]) - sig2paramclean.xx[j]/2)*
#         (1-pnorm(-xiparamclean.xx[j]*sqrt(sig2paramclean.xx[j])/sqrt(xiparamclean.xx[j]^2+1)))
#     }
#     
#     #% liquid weight ***** component
#     WJ1.xx = alphaparamclean.xx/(alphaparamclean.xx + n.train)*f0X.xx[h] #; print(paste("liquid weight=", WJ1))
#     
#     #% Collecting values and constructing weights for each iteration
#     W_paramFree.vec.xx[h,r] = WJ1.xx/(WJ1.xx + sum(w_j.solid.xx))
#     W_paramBase.matrix.xx[h, ] = w_j.solid.xx/(WJ1.xx + sum(w_j.solid.xx)) 
#     
#     #% Collecting values for predictive density (via weighted AVG)
#     density.train.xx[h,r] = W_paramFree.vec.xx[h,r]*f0S.xx[h] + sum(W_paramBase.matrix.xx[h,]*f_S.xx)
#     #% Collecting values for prediction (via weighted AVG)
#     expval.train.xx[h,r] = W_paramFree.vec.xx[h,r]*E0S.xx[h] + sum(W_paramBase.matrix.xx[h,]*E_value.xx)
#   }
#   print(paste("r=", r))
# }
# 
# expS.DP.Avg.train.xx <- apply(X=expval.train.xx, MARGIN=1, FUN=mean)
# 
# summary(S)
# summary(expS.DP.Avg.train.xx)
# summary(log(S))
# summary(log(expS.DP.Avg.train.xx))

################################################################################
###> Testing data -------------------------------------------------------------
################################################################################ 
###> Parameter-free component
#> Analytical form
# f0x1 = function(Z) {
#   beta(Z+g0.xx, 1-Z+h0.xx)/beta(g0.xx,h0.xx)
# }
# #> Analytical form
# f0x2 = function(X, covar) {
#   1/sqrt(pi)*gamma(c0.xx+1)/gamma(c0.xx)*d0.xx^(c0.xx/2)/( (X-mean(covar))^2 + d0.xx )^((c0.xx+1)/2)
# }
# #> JOINT 
# f0X.test.xx = numeric(n.test) #% param-free covariate model f0(X,Z)
# 
# ##### For [Param-free outcome model]: -------------------> MonteCarlo Integration
# # Calculate Outcome and Covariate parameter free data model for each observation  
# f0S.test.xx = numeric(n.test) #% param-free outcome model f0(S|X,Z)
# E0S.test.xx = numeric(n.test) #% param-free E(S|X,Z) = Expected value of S|X,Z
# 
# M = 1000 # Number of Monte Carlo samples
# sumS.xx = numeric(M)
# sumES.xx = numeric(M)
# 
# for(h in 1:n.test) {
#   
#   ##> Analytical solution for f0(X, Z) for ### covariate model ###
#   # : necessary parameters are supplied by pre-determined hyperparam values...
#   f0X.test.xx[h] = f0x1(Z.test[h])*f0x2(Xstar.test[h], Xstar.test)  
#   
#   ##> Monte Carlo Integration for the ### outcome model ###
#   for(j in 1:M) {
#     xi_samplej = rt(n=1, df=nu0.xx)                               # bullet on xi
#     sig_samplej = rinvgamma(n=1, shape=u0.xx, scale=v0.xx)           # bullet on sig2
#     beta_samplej = rmvn(n=1, mu=m0.xx, sigma=SIG_b0.xx*varinf)       # bullet on beta
#     
#     sumS.xx[j] = dlogsknorm( s=S.test[h], x=matXstar.test[h,], beta=beta_samplej, sig2=sig_samplej, xi=xi_samplej )* 
#       dmvn(X=beta_samplej, mu=m0.xx, sigma=SIG_b0.xx*varinf)*        # to joint beta
#       dinvgamma(x=sig_samplej, shape=u0.xx, scale=v0.xx)*            # to joint sig2
#       dt(x=xi_samplej, df=nu0.xx)                                 # to joint xi
#     # 20*exp( betaparamclean[j,1] + betaparamclean[j,2]*cleanx + betaparamclean[j,3]*Z[i] - sig2paramclean[j]/2)*( 1-pnorm(-xiparamclean[j]*sqrt(sig2paramclean[j])/sqrt(xiparamclean[j]^2+1)))
#     sumES.xx[j] = 18*exp( sum(matXstar.test[h,]*beta_samplej) - sig_samplej/2 )*(1-pnorm(-xi_samplej*sqrt(sig_samplej)/sqrt(xi_samplej^2+1)))
#   }
#   E0S.test.xx[h] = sum(sumES.xx)/M  # E[S] value of f0
#   f0S.test.xx[h] = sum(sumS.xx)/M   # density of f0
#   
#   print(paste("h=",h))
# }
# summary(S.test)    # original S
# summary(E0S.test.xx)  # E[S] based on f0(S|x) by train data




###> NOW....
# cl_membership.test <- full.test.sweden$Risk
# table(cl_membership.test)                    # this is wrong!!!!!!!!!!!!!!!!!

#% Ultimate liquid weight vector for each iteration
W_paramFree.vec.test.xx = matrix(0, nrow = n.test, ncol = total_iter_xx) 

density.test.xx = matrix(0, nrow = n.test, ncol = total_iter_xx)
expval.test.xx = matrix(0, nrow = n.test, ncol = total_iter_xx)


for(r in 1:(total_iter_xx)) {
  piparamclean.xx = list_piparam.xx[[r]]
  lambda2paramclean.xx = list_lambda2param.xx[[r]] #for Xstar
  alphaparamclean.xx = list_alphaparam.xx[[r]]      #for alpha
  
  betaparamclean.xx = list_betaparam.xx[[r]]      #for S
  sig2paramclean.xx = list_sig2param.xx[[r]]      #for S
  xiparamclean.xx = list_xiparam.xx[[r]]       #for S
  
  cl_membership.test.xx = list_cl.xx[[r]] # 
  
  
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
      E_value.test.xx[j] = 12*exp(sum( c(1, Xstar.test[h], Z.test[h])*betaparamclean.xx[j,]) - sig2paramclean.xx[j]/2)*
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

SSPE.DP.LSN.xx <- sum( (log(expS.DP.Avg.test.xx) - log(S.test))^2 ); SSPE.DP.LSN.xx   # 17173.6
SAPE.DP.LSN.xx <- sum( abs(log(expS.DP.Avg.test.xx) - log(S.test)) ); SAPE.DP.LSN.xx  # 4757.152
lppd_EX_DP.S.xx                                                                       # -663949.5



#> beta0
list_betaparam0.xx <- numeric(100)
for (i in 1:100){
  list_betaparam0.xx[i] <- mean( list_betaparam.xx[[i]][, 1] ) 
}
summary(list_betaparam0.xx)# 10.14

#> beta1
list_betaparam1.xx <- numeric(100)
for (i in 1:100){
  list_betaparam1.xx[i] <- mean( list_betaparam.xx[[i]][, 2] ) 
}
summary(list_betaparam1.xx)# 0.4    

#> beta2
list_betaparam2.xx <- numeric(100)
for (i in 1:100){
  list_betaparam2.xx[i] <- mean( list_betaparam.xx[[i]][, 3] ) 
}
summary(list_betaparam2.xx) #-1.92 

#> sig2
list_sig2test.xx <- numeric(100)
for (i in 1:100) {
  list_sig2test.xx[i] <- mean( list_sig2param.xx[[i]] )
}
summary(list_sig2test.xx ) #7.2   


#> xi
list_xitest.xx <- numeric(100)
for (i in 1:100) {
  list_xitest.xx[i] <- mean( list_xiparam.xx[[i]] )
}
summary(list_xitest.xx ) #-0.81


















#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##--------- [ model(C): Gustafson DPM ] - multiple[1%X3, 10%X3, 25%X3] -------##
###> Load shard output
sh_out_st <- list()

#1) zeta =0.99 -----------------------------------------------------------------
for(j in 1:Nsh) {
  sh_out_st[[j]] <- readRDS( paste("testtt/st.shard_0.01_0.99(",j,")data.rds",sep="") )
}

#1) zeta =0.90
for(j in 1:Nsh) {
  sh_out_st[[j]] <- readRDS( paste("testtt/st.shard_0.01_0.90(",j,")data.rds",sep="") )
}

#1) zeta =0.95
for(j in 1:Nsh) {
  sh_out_st[[j]] <- readRDS( paste("testtt/st.shard_0.01_0.95(",j,")data.rds",sep="") )
}

#1) zeta =0.85
for(j in 1:Nsh) {
  sh_out_st[[j]] <- readRDS( paste("testtt/st.shard_0.01_0.85(",j,")data.rds",sep="") )
}
#-------------------------------------------------------------------------------


###> Retrieve total iterations
total_iter_st <- length(sh_out_st[[1]][[3]]) # in the shard NO.1, third element: "list_cl"...

###> Retrieve total sample size
N <- sh_out_st[[1]][[17]]

###> Retrieve alpha paramater ( not cl-wise but shard-wise )
alphaparam_st_new <- numeric(total_iter_st)
total_observations.st <- 0
# - linear combination: weighted AVG
for(j in 1:Nsh) {                                 #index in the shard...........list_alphaparam
  alphaparam_st_new <- alphaparam_st_new + length(sh_out_st[[j]][[1]])*unlist(sh_out_st[[j]][[8]])
  total_observations.st <- total_observations.st + length(sh_out_st[[j]][[1]])
}
alphaparam_st_new <- alphaparam_st_new/total_observations.st; alphaparam_st_new

###> Retrieve anchor points (in the 1st shard, the second item: )
anchor.st <- sh_out_st[[1]][[2]]; anchor.st

###> Determine the threshold
epsilon <- 0.65 #0.5


###> This is initial cl_membership of the 1st shard, tabulating on the entire sample matrix(n.train)
cl_membership_st.1_matrix <- list()
for(r in 1:total_iter_st) {
  cl_membership_st.1_matrix[[r]] <- matrix( 0, nrow=N, ncol=max(sh_out_st[[1]][[3]][[r]]) ) # no anchor pt...no 1st shard
  # for obv in the 1st shard, obtain cl_memberships....and point out the "anchor pt"!!!!!!!!
  for(i in 1:length(sh_out_st[[1]][[1]])) {
    cl_membership_st.1_matrix[[r]][ sh_out_st[[1]][[1]][i], sh_out_st[[1]][[3]][[r]][i] ] <- 1
  }
}
cl_membership_st.1_matrix # this is the [1st shard] including its anchor points. 
#            cl.01,      cl.02,     cl.03, ....
#obv1: anchor(Y/N) anchor(Y/N) anchor(Y/N) .....
#obv2: anchor(Y/N) anchor(Y/N) anchor(Y/N) .......
# ...      ...        ...         ...      ......... in the [1st shard]
#n.train

#:::: the cl in the [1st shard] will grow...the anchor pt will duplicates....


###> Next, 
cl_membership_st.1_vector <- list()
loglikelihood.st = numeric(total_iter_st)
# 
for (r in 1:total_iter_st) {
  # -> pick up cl_membership (i.e. "which cl" the data point belong to..) for the 1st shard
  cl_membership_st.1 <- cl_membership_st.1_matrix[[r]]
  #::: every info in the 1st shard looks great, arranged on a single matrix...but 2nd, 3rd,...shards are not...so.. 
  
  # -> pick up cl_membership for the 2nd, 3rd, .... shard to compare
  for(j in 2:Nsh) {
    cl_membership_st.2 <- sh_out_st[[j]][[3]][[r]] 
    n_clust_st.2 <- max(cl_membership_st.2) # in shard[j] at iter[r], CHECK how many clusters?
    
    # for each cluster "k2", 
    for(k2 in 1:n_clust_st.2) {
      ##### 1) challenger shard info
      indx_sh2_clk2 <- sh_out_st[[j]][[1]][cl_membership_st.2==k2] # collect indexes that belong to the cluster "k2"
      anchor_sh2_clk2 <- anchor.st[anchor.st %in% indx_sh2_clk2]         # collect the [anchor points] in the the cluster "k2"
      
      
      #:::take these indexes and compare them with every indexes in the 1st shard!!!!
      n_clust_st.1 <- ncol(cl_membership_st.1) # how many cl in the 1st shard?
      k1=1
      MERGED=FALSE
      # -> If..Merged:
      while ( k1<=n_clust_st.1 & !MERGED ) {
        ##### 2) Champion shard info
        indx_sh1_clk1 <- which(cl_membership_st.1[, k1]==1)   # collect indexes that belong to the cluster "k1" from the matrix
        anchor_sh1_clk1 <- anchor.st[anchor.st %in% indx_sh1_clk1] # collect the [anchor points] in the 1st shard (k1)
        
        #> compute:
        C_dist <- sum(anchor_sh1_clk1 %in% anchor_sh2_clk2)                     # common elements COUNT
        D_dist <- length(anchor_sh1_clk1) + length(anchor_sh2_clk2) - 2*C_dist  # different elements COUNT
        
        d_distance <- ifelse(C_dist > 0 | D_dist > 0, D_dist/(C_dist+D_dist), 1)
        #d_distance <- D_dist/(C_dist+D_dist)
        
        #print(paste("dist:"), d_distance)
        print(paste("dist:",d_distance))
        
        # FINALLYYYYYY COMPUTE: weight AVG of paramssssss
        if(d_distance < epsilon) {  # Merge clusters
          cl_membership_st.1[indx_sh2_clk2, k1] <- 1
          sh_out_st[[1]][[4]][[r]][k1] <- length(indx_sh1_clk1)/(length(indx_sh1_clk1) + length(indx_sh2_clk2))*sh_out_st[[1]][[4]][[r]][k1] +
            length(indx_sh2_clk2)/(length(indx_sh1_clk1) + length(indx_sh2_clk2))*sh_out_st[[j]][[4]][[r]][k2] # update for piparameter
          sh_out_st[[1]][[5]][[r]][k1] <- length(indx_sh1_clk1)/(length(indx_sh1_clk1) + length(indx_sh2_clk2))*sh_out_st[[1]][[5]][[r]][k1] +
            length(indx_sh2_clk2)/(length(indx_sh1_clk1) + length(indx_sh2_clk2))*sh_out_st[[j]][[5]][[r]][k2] # update for lambda2parameter
          sh_out_st[[1]][[6]][[r]][k1,] <- length(indx_sh1_clk1)/(length(indx_sh1_clk1) + length(indx_sh2_clk2))*sh_out_st[[1]][[6]][[r]][k1,] +
            length(indx_sh2_clk2)/(length(indx_sh1_clk1) + length(indx_sh2_clk2))*sh_out_st[[j]][[6]][[r]][k2,] # update for kappaparameter
          sh_out_st[[1]][[7]][[r]][k1] <- length(indx_sh1_clk1)/(length(indx_sh1_clk1) + length(indx_sh2_clk2))*sh_out_st[[1]][[7]][[r]][k1] +
            length(indx_sh2_clk2)/(length(indx_sh1_clk1) + length(indx_sh2_clk2))*sh_out_st[[j]][[7]][[r]][k2] # update for tau2parameter
          sh_out_st[[1]][[9]][[r]][k1,] <- length(indx_sh1_clk1)/(length(indx_sh1_clk1) + length(indx_sh2_clk2))*sh_out_st[[1]][[9]][[r]][k1,] +
            length(indx_sh2_clk2)/(length(indx_sh1_clk1) + length(indx_sh2_clk2))*sh_out_st[[j]][[9]][[r]][k2,] # update for betaparameter
          sh_out_st[[1]][[10]][[r]][k1] <- length(indx_sh1_clk1)/(length(indx_sh1_clk1) + length(indx_sh2_clk2))*sh_out_st[[1]][[10]][[r]][k1] +
            length(indx_sh2_clk2)/(length(indx_sh1_clk1) + length(indx_sh2_clk2))*sh_out_st[[j]][[10]][[r]][k2] # update for sig2parameter
          sh_out_st[[1]][[11]][[r]][k1] <- length(indx_sh1_clk1)/(length(indx_sh1_clk1) + length(indx_sh2_clk2))*sh_out_st[[1]][[11]][[r]][k1] +
            length(indx_sh2_clk2)/(length(indx_sh1_clk1) + length(indx_sh2_clk2))*sh_out_st[[j]][[11]][[r]][k2] # update for xiparameter
          sh_out_st[[1]][[12]][[r]][k1,] <- length(indx_sh1_clk1)/(length(indx_sh1_clk1) + length(indx_sh2_clk2))*sh_out_st[[1]][[12]][[r]][k1,] +
            length(indx_sh2_clk2)/(length(indx_sh1_clk1) + length(indx_sh2_clk2))*sh_out_st[[j]][[12]][[r]][k2,] # update for clean.kappaparameter
          sh_out_st[[1]][[13]][[r]][k1] <- length(indx_sh1_clk1)/(length(indx_sh1_clk1) + length(indx_sh2_clk2))*sh_out_st[[1]][[13]][[r]][k1] +
            length(indx_sh2_clk2)/(length(indx_sh1_clk1) + length(indx_sh2_clk2))*sh_out_st[[j]][[13]][[r]][k2] # update for clean.lambda2parameter
          sh_out_st[[1]][[14]][[r]][k1,] <- length(indx_sh1_clk1)/(length(indx_sh1_clk1) + length(indx_sh2_clk2))*sh_out_st[[1]][[14]][[r]][k1,] +
            length(indx_sh2_clk2)/(length(indx_sh1_clk1) + length(indx_sh2_clk2))*sh_out_st[[j]][[14]][[r]][k2,] # update for clean.betaparameter
          sh_out_st[[1]][[15]][[r]][k1] <- length(indx_sh1_clk1)/(length(indx_sh1_clk1) + length(indx_sh2_clk2))*sh_out_st[[1]][[15]][[r]][k1] +
            length(indx_sh2_clk2)/(length(indx_sh1_clk1) + length(indx_sh2_clk2))*sh_out_st[[j]][[15]][[r]][k2] # update for clean.sig2parameter
          sh_out_st[[1]][[16]][[r]][k1] <- length(indx_sh1_clk1)/(length(indx_sh1_clk1) + length(indx_sh2_clk2))*sh_out_st[[1]][[16]][[r]][k1] +
            length(indx_sh2_clk2)/(length(indx_sh1_clk1) + length(indx_sh2_clk2))*sh_out_st[[j]][[16]][[r]][k2] # update for clean.xiparameter
          
          MERGED = TRUE
        }
        k1=k1+1
      } #END of [while] 
      
      # -> IF.. not MERGED: 
      if (!MERGED){ 
        cl_membership_st.1 <- cbind(cl_membership_st.1, 0)
        cl_membership_st.1[indx_sh2_clk2, ncol(cl_membership_st.1)] <- 1
        sh_out_st[[1]][[4]][[r]] <- c(sh_out_st[[1]][[4]][[r]],sh_out_st[[j]][[4]][[r]][k2]) # pi ::: append
        sh_out_st[[1]][[5]][[r]] <- c(sh_out_st[[1]][[5]][[r]],sh_out_st[[j]][[5]][[r]][k2]) # lambda2 ::: append
        sh_out_st[[1]][[6]][[r]] <- rbind(sh_out_st[[1]][[6]][[r]],sh_out_st[[j]][[6]][[r]][k2,]) # kappa ::: append
        sh_out_st[[1]][[7]][[r]] <- c(sh_out_st[[1]][[7]][[r]],sh_out_st[[j]][[7]][[r]][k2]) # tau2 ::: append
        sh_out_st[[1]][[9]][[r]] <- rbind(sh_out_st[[1]][[9]][[r]],sh_out_st[[j]][[9]][[r]][k2,]) # beta ::: append
        sh_out_st[[1]][[10]][[r]] <- c(sh_out_st[[1]][[10]][[r]],sh_out_st[[j]][[10]][[r]][k2]) # sig2 ::: append
        sh_out_st[[1]][[11]][[r]] <- c(sh_out_st[[1]][[11]][[r]],sh_out_st[[j]][[11]][[r]][k2]) # xi ::: append
        sh_out_st[[1]][[12]][[r]] <- rbind(sh_out_st[[1]][[12]][[r]],sh_out_st[[j]][[12]][[r]][k2,]) # clean.kappa ::: append
        sh_out_st[[1]][[13]][[r]] <- c(sh_out_st[[1]][[13]][[r]],sh_out_st[[j]][[13]][[r]][k2]) # clean.lambda2 ::: append
        sh_out_st[[1]][[14]][[r]] <- rbind(sh_out_st[[1]][[14]][[r]],sh_out_st[[j]][[14]][[r]][k2,]) # clean.beta ::: append
        sh_out_st[[1]][[15]][[r]] <- c(sh_out_st[[1]][[15]][[r]],sh_out_st[[j]][[15]][[r]][k2]) # clean.sig2 ::: append
        sh_out_st[[1]][[16]][[r]] <- c(sh_out_st[[1]][[16]][[r]],sh_out_st[[j]][[16]][[r]][k2]) # clean.xi ::: append
      }
    }
  }
  
  ###> Remove duplicates
  for (a in anchor.st) {
    X_a <- cl_membership_st.1[a, ]
    
    # Be careful............
    if(sum(X_a)>1) {
      which.is.1 <- which(X_a==1) # vector
      keep1 <- sample(which.is.1, 1)
      cl_membership_st.1[a, -keep1] <- 0
    } # since some anchor pt has its own cluster
  }
  cl_membership_st.1_vector[[r]] <- apply(cl_membership_st.1, 1, function(x) which(x==1))
  cl_membership_st.1_matrix[[r]] <- cl_membership_st.1
  
  ###> Compute LPPD
  #loglike_r = 0
  loglike_r_vector = numeric(N)
  
  for(i in 1:N) {
    x_i = matXstar[i, ]
    jj = cl_membership_st.1_vector[[r]][[i]]
    
    #loglike_r = loglike_r +
    loglike_r_vector[i] =
      dlogsknorm_log(s=S[i], x=x_i, beta=sh_out_st[[1]][[9]][[r]][jj, ], sig2=sh_out_st[[1]][[10]][[r]][jj], xi=sh_out_st[[1]][[11]][[r]][jj]) +
      dbinom(x = x_i[3], size = 1, prob = sh_out_st[[1]][[4]][[r]][jj], log = TRUE) + 
      dnorm(x = x_i[2], mean = mean(Xstar), sd = sqrt(sh_out_st[[1]][[5]][[r]][jj]), log = TRUE)
  }
  #loglikelihood.st[r] = loglike_r # cumulative sum
  loglikelihood.st[r] = sum(loglike_r_vector) # cumulative sum
}



###> Now it's HARVEST TIME! ----------------------------------------------------
###> 
###> CONNECT to prediction!!!
list_piparam.st = sh_out_st[[1]][[4]]
list_tau2param.st = sh_out_st[[1]][[7]]
list_alphaparam.st = alphaparam_st_new

list_kappaparam.st = sh_out_st[[1]][[6]]
list_lambda2param.st = sh_out_st[[1]][[5]]
list_betaparam.st = sh_out_st[[1]][[9]]
list_sig2param.st = sh_out_st[[1]][[10]]
list_xiparam.st = sh_out_st[[1]][[11]]

list_kappaparamC.st = sh_out_st[[1]][[12]]
list_lambda2paramC.st = sh_out_st[[1]][[13]]
list_betaparamC.st = sh_out_st[[1]][[14]]
list_sig2paramC.st = sh_out_st[[1]][[15]]
list_xiparamC.st = sh_out_st[[1]][[16]]

list_cl.st = cl_membership_st.1_vector

n.train = N

options(scipen = 999)
plot( loglikelihood.st, type="l")
lppd_EX_DP.S.st <- mean( loglikelihood.st ); lppd_EX_DP.S.st #-668987 -668767.4

#..............................................................................#
#..............................................................................#
############################### Predictions ####################################
#..............................................................................#
#..............................................................................#

################################################################################
###> Training data -------------------------------------------------------------
################################################################################

# #% Ultimate liquid weight vector for each iteration
# W_paramFree.vec.st = matrix(0, nrow = n.train, ncol = total_iter_st) 
# density.train.st = matrix(0, nrow = n.train, ncol = total_iter_st)
# expval.train.st = matrix(0, nrow = n.train, ncol = total_iter_st)
# 
# for(r in 1:(total_iter_st)) {
#   piparamclean.st = list_piparam.st[[r]]
#   lambda2paramclean.st = list_lambda2param.st[[r]] #for Xstar
#   alphaparamclean.st = list_alphaparam.st[[r]]      #for alpha
#   
#   betaparamclean.st = list_betaparam.st[[r]]      #for S
#   sig2paramclean.st = list_sig2param.st[[r]]      #for S
#   xiparamclean.st = list_xiparam.st[[r]]       #for S
#   
#   cl_membership.train.st = list_cl.st[[r]] # 
#   
#   J = nrow(betaparamclean.st)
#   f_S.st = numeric(J)                  # density for S
#   
#   #% Ultimate solid weight matrix
#   W_paramBase.matrix.st = matrix(0, nrow=n.train, ncol=J)
#   
#   for(h in 1:n.train) {
#     
#     fX_j.st = numeric(J)       # discrete covariate model for each j
#     w_j.solid.st = numeric(J)  # discrete weight for each j  
#     E_value.st = numeric(J)    # to save E[S|X] for each j
#     
#     for(j in 1:J) {
#       # discrete covariate joint
#       fX_j.st[j] = dbinom(x=Z[h], size=1, prob=piparamclean.st[j])*
#         dnorm(x=Xstar[h], mean=mean(Xstar), sd=sqrt(lambda2paramclean.st[j]))
#       
#       #% solid weight ***** component
#       w_j.solid.st[j] = length(S[cl_membership.train.st==j])/(alphaparamclean.st + n.train)*fX_j.st[j]
#       
#       #% compute pred.density
#       f_S.st[j] = dlogsknorm(s=S[h], x=c(1, Xstar[h], Z[h]), 
#                              beta=betaparamclean.st[j,], sig2=sig2paramclean.st[j], xi=xiparamclean.st[j])
#       #% compute pred.value
#       E_value.st[j] = 2*exp(sum( c(1, Xstar[h], Z[h])*betaparamclean.st[j,])-sig2paramclean.st[j]/2)*
#         (1-pnorm(-xiparamclean.st[j]*sqrt(sig2paramclean.st[j])/sqrt(xiparamclean.st[j]^2+1)))
#     }
#     
#     #% liquid weight ***** component
#     WJ1.st = alphaparamclean.st/(alphaparamclean.st + n.train)*f0X.st[h] #; print(paste("liquid weight=", WJ1))
#     
#     #% Collecting values and constructing weights for each iteration
#     W_paramFree.vec.st[h,r] = WJ1.st/(WJ1.st + sum(w_j.solid.st))
#     W_paramBase.matrix.st[h, ] = w_j.solid.st/(WJ1.st + sum(w_j.solid.st)) 
#     
#     #% Collecting values for predictive density (via weighted AVG)
#     density.train.st[h,r] = W_paramFree.vec.st[h,r]*f0S.st[h] + sum(W_paramBase.matrix.st[h,]*f_S.st)
#     #% Collecting values for prediction (via weighted AVG)
#     expval.train.st[h,r] = W_paramFree.vec.st[h,r]*E0S.st[h] + sum(W_paramBase.matrix.st[h,]*E_value.st)
#   }
#   print(paste("r=", r))
# }
# 
# expS.DP.Avg.train.st <- apply(X=expval.train.st, MARGIN=1, FUN=mean)
# 
# summary(S)
# summary(expS.DP.Avg.train.st)
# summary(log(S))
# summary(log(expS.DP.Avg.train.st))

################################################################################
###> Testing data -------------------------------------------------------------
################################################################################
##### For Param-free covariate model: 
# ---> int_{}^{} Xstar,z,piparam,kappaparam,lambda2param w.r.t "piparam","kappaparam","lambda2param"
#> JOINT 
# f0X.test.st = numeric(n.test) #% to save param-free covariate model: JOINT: f0(Xstar, Z)
# # pre-discovered...by the analytical soltution...
# constX = (d0.st/2)^(c0.st/2)*gamma((c0.st+1)/2)/(sqrt(2*pi)*gamma(c0.st/2)*det(SIG_k0)^0.5)
# 
# 
# ##### For [Param-free outcome model]: -------------------> MonteCarlo Integration
# # Calculate Outcome and Covariate parameter free data model for each observation 
# f0S.test.st = numeric(n.test) #% param-free outcome model f0(S|x)
# E0S.test.st = numeric(n.test) # E(S|x) = Expected value of S|x ~ f0(S|x)
# 
# ##### Now let's solve the integralsssss
# set.seed(1)
# M = 1000 # iterations of the Monte Carlo integral
# for(h in 1:n.test) {
#   # Analytical solution for f0(Xstar, Z) for ### covariate model ###
#   # : necessary parameters are supplied by pre-determined hyperparam values...so ~~~~~~~~~~~~~~~~~~~~~~#Q.David
#   f0X.test.st[h] = constX*beta(g0.st+Z.test[h], h0.st+1-Z[h])/beta(g0.st, h0.st)*det(KK1[h,] %*% t(KK1[h,])+SIG_inv_k0)^(-0.5)/
#     (1/2*(d0.st+(Xstar.test[h]-KK1[h,] %*% kappa_v0)^2/(1+t(KK1[h,]) %*% SIG_k0 %*% KK1[h,])))^((c0.st+1)/2) 
#   
#   # Monte Carlo integration for S for ### outcome model ###
#   # : necessary parameters are supplied by sampling from prior density
#   sumS.st = numeric(M)
#   sumES.st = numeric(M)
#   for(j in 1:M) {
#     xi_samplej = rt(n = 1, df = nu0.st)                               # prior on xi
#     sig_samplej = rinvgamma(n = 1, shape = u0.st, scale = v0.st)         # prior on sig2
#     beta_samplej = rmvn(n = 1, mu = m0.st, sigma = SIG_b0.st*varinf)  # prior on beta
#     
#     sumS.st[j] = 
#       dlogsknorm( s=S.test[h],
#                   x=matXstar.test[h,],
#                   beta=beta_samplej,
#                   sig2=sig_samplej,
#                   xi=xi_samplej )*                                 # outcome with complete
#       dmvn(X = beta_samplej, mu = m0.st, sigma = SIG_b0.st*varinf)*   # to joint beta
#       dinvgamma(x = sig_samplej, shape = u0.st, scale = v0.st)*          # to joint sig2
#       dt(x = xi_samplej, df = nu0.st)                                 # to joint xi
#     sumES.st[j] = 2*exp(sum(matXstar.test[h,]*beta_samplej) - sig_samplej/2)*(1-pnorm(-xi_samplej*sqrt(sig_samplej)/sqrt(xi_samplej^2+1)))
#     
#   }
#   E0S.test.st[h] = sum(sumES.st)/M  # E[S] value of f0
#   f0S.test.st[h] = sum(sumS.st)/M   # density of f0
#   
#   print(paste("h=",h))
# }
# plot( x=density(f0S.test.st) )
# plot( x=density(f0X.test.st) )
# summary(E0S.test.st)  # E[S] based on f0(S|x) by train data
# summary(S)    # original S
# summary(f0X.test.st)
# plot(x=sort(E0S.test.st), y=f0S.st[order(E0S.test.st)], type="l")



############################### Out-of-sample ##################################
#% Ultimate liquid weight vector for each iteration
W_paramFree.vec.test.st = matrix(0, nrow = n.test, ncol = total_iter_st) 
density.test.st = matrix(0, nrow = n.test, ncol = total_iter_st)
expval.test.st = matrix(0, nrow = n.test, ncol = total_iter_st)

for(r in 1:(total_iter_st)) {
  piparamclean.st = list_piparam.st[[r]]
  lambda2paramclean.st = list_lambda2param.st[[r]] #for Xstar
  alphaparamclean.st = list_alphaparam.st[[r]]      #for alpha
  
  betaparamclean.st = list_betaparam.st[[r]]      #for S
  sig2paramclean.st = list_sig2param.st[[r]]      #for S
  xiparamclean.st = list_xiparam.st[[r]]       #for S
  
  cl_membership.test.st = list_cl.st[[r]] # 
  
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
      w_j.solid.test.st[j] = sum(cl_membership.train.st==j) / (alphaparamclean.st + n.train)*fX_j.test.st[j] #~~~~~~~~QQQ
      #% compute pred.density
      f_S.test.st[j] = dlogsknorm(s=S.test[h], x=c(1, Xstar.test[h], Z.test[h]), 
                                  beta=betaparamclean.st[j,], sig2=sig2paramclean.st[j], xi=xiparamclean.st[j])
      #% compute pred.value
      E_value.test.st[j] = 2*exp(sum( c(1, Xstar.test[h], Z.test[h])*betaparamclean.st[j,])-sig2paramclean.st[j]/2)*
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

SSPE.DP.LSN.st <- sum( (log(expS.DP.Avg.test.st) - log(S.test))^2 ); SSPE.DP.LSN.st   #13128.62 12872.54
SAPE.DP.LSN.st <- sum( abs(log(expS.DP.Avg.test.st) - log(S.test)) ); SAPE.DP.LSN.st  #4344.899 4254.066




#> beta0 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
list_betaparam0.st <- numeric(100)
for (i in 1:100){
  list_betaparam0.st[i] <- mean( list_betaparam.st[[i]][, 1] ) 
}
summary(list_betaparam0.st)# 10.12(zeta0.85);  10.096(zeta0.9); (zeta0.95); 10.14(zeta0.99)

# list_betaparam0.stst <- numeric(100)
# for (i in 1:100){
#   list_betaparam0.stst[i] <- mean( list_betaparam.clean.st[[i]][, 1] ) 
# }
# summary(list_betaparam0.stst)# 


#> beta1
list_betaparam1.st <- numeric(100)
for (i in 1:100){
  list_betaparam1.st[i] <- mean( list_betaparam.st[[i]][, 2] ) 
}
summary(list_betaparam1.st)# 0.38(zeta0.85);  0.39(zeta0.9); (zeta0.95); 0.37(zeta0.99)

# list_betaparam1.stst <- numeric(2000)
# for (i in 1:2000){
#   list_betaparam1.stst[i] <- mean( list_betaparam.clean.st[[i]][, 2] ) 
# }
# summary(list_betaparam1.stst)# 



#> beta2
list_betaparam2.st <- numeric(100)
for (i in 1:100){
  list_betaparam2.st[i] <- mean( list_betaparam.st[[i]][, 3] ) 
}
summary(list_betaparam2.st) # -1.91(zeta0.85);  -1.89(zeta0.9); (zeta0.95); -1.87(zeta0.99)

# list_betaparam2.stst <- numeric(2000)
# for (i in 1:2000){
#   list_betaparam2.stst[i] <- mean( list_betaparam.clean.st[[i]][, 3] ) 
# }
# summary(list_betaparam2.stst) 



#> sig2
list_sig2test.st <- numeric(100)
for (i in 1:100) {
  list_sig2test.st[i] <- mean( list_sig2param.st[[i]] )
}
summary(list_sig2test.st ) # ?(zeta0.85);  6.154(zeta0.9); (zeta0.95); 7.12(zeta0.99)

# list_sig2test.stst <- numeric(2000)
# for (i in 1:2000) {
#   list_sig2test.stst[i] <- mean( list_sig2param.clean.st[[i]] )
# }
# summary(list_sig2test.stst ) 


#> xi
list_xitest.st <- numeric(100)
for (i in 1:100) {
  list_xitest.st[i] <- mean( list_xiparam.st[[i]] )
}
summary(list_xitest.st ) # ?(zeta0.85);  -0.39(zeta0.9); (zeta0.95); -0.34(zeta0.99)

# list_xitest.stst <- numeric(2000)
# for (i in 1:2000) {
#   list_xitest.stst[i] <- mean( list_xiparam.clean.st[[i]] )
# }
# summary(list_xitest.stst ) 



#> tau2
list_tau2test.st <- numeric(100)
for (i in 1:100) {
  list_tau2test.st[i] <- mean( list_tau2param.st[[i]] )
}
summary(list_tau2test.st ) # ?(zeta0.85);  0.54(zeta0.9); (zeta0.95); 0.09(zeta0.99)




#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##--------- [ model(D): Gustafson HB ] - multiple[1%X3, 10%X3, 25%X3] --------##

###> Load shard output
sh_out_HB.st <- list()

#1) zeta=0.99 ------------------------------------------------------------------
for(j in 1:Nsh) {
  sh_out_HB.st[[j]] <- readRDS( paste("testtt/HB.st.shard_0.01_0.99(",j,")data.rds",sep="") )
}

#1) zeta=0.90
for(j in 1:Nsh) {
  sh_out_HB.st[[j]] <- readRDS( paste("testtt/HB.st.shard_0.01_0.90(",j,")data.rds",sep="") )
}

#1) zeta=0.95
for(j in 1:Nsh) {
  sh_out_HB.st[[j]] <- readRDS( paste("testtt/HB.st.shard_0.01_0.95(",j,")data.rds",sep="") )
}

#1) zeta=0.85
for(j in 1:Nsh) {
  sh_out_HB.st[[j]] <- readRDS( paste("testtt/HB.st.shard_0.01_0.85(",j,")data.rds",sep="") )
}
#-------------------------------------------------------------------------------



#...............................................................................
#...............................................................................
#................................. WARNING ! ...................................
#...............................................................................
# - for cl_mistake in the beginning ............................................
# list_cl.HB.st <- list()  
# for (i in 1:total_iter_HB.st ) {
#   list_cl.HB.st[[i]] <- cl_membership.HB.train  # Assign the vector to each list element
# }
# for(j in 1:Nsh){
#   sh_out_HB.st[[j]] <- append( sh_out_HB.st[[j]], list(list_cl.HB.st), after = 2 )
# }

###> Retrieve total iterations
total_iter_HB.st <- length(sh_out_HB.st[[1]][[3]]) # in the shard NO.1, third element: "list_cl"...
#...............................................................................
#...............................................................................
#...............................................................................






###> Retrieve total sample size
N <- sh_out_HB.st[[1]][[16]]

###> Retrieve anchor points (in the 1st shard, the second item: )
anchor.HB.st <- sh_out_HB.st[[1]][[2]]; anchor.HB.st

###> Determine the threshold
epsilon <- 0.65



###> This is initial cl_membership of the 1st shard, tabulating on the entire sample matrix(n.train)
cl_membership_HB.st.1_matrix <- list()
for(r in 1:total_iter_HB.st) {
  cl_membership_HB.st.1_matrix[[r]] <- matrix( 0, nrow=N, ncol=max(sh_out_HB.st[[1]][[3]][[r]]) ) # no anchor pt...no 1st shard
  # for obv in the 1st shard, obtain cl_memberships....and point out the "anchor pt"!!!!!!!!
  for(i in 1:length(sh_out_HB.st[[1]][[1]])) {
    cl_membership_HB.st.1_matrix[[r]][ sh_out_HB.st[[1]][[1]][i], sh_out_HB.st[[1]][[3]][[r]][i] ] <- 1
  }
}
cl_membership_HB.st.1_matrix # this is the [1st shard] including its anchor points. 
#            cl.01,      cl.02,     cl.03, ....
#obv1: anchor(Y/N) anchor(Y/N) anchor(Y/N) .....
#obv2: anchor(Y/N) anchor(Y/N) anchor(Y/N) .......
# ...      ...        ...         ...      ......... in the [1st shard]
#n.train

#:::: the cl in the [1st shard] will grow...the anchor pt will duplicates....


###> Next, 
cl_membership_HB.st.1_vector <- list()
loglikelihood.HB.st = numeric(total_iter_HB.st)
# 
for (r in 1:total_iter_HB.st) {
  # -> pick up cl_membership (i.e. "which cl" the data point belong to..) for the 1st shard
  cl_membership_HB.st.1 <- cl_membership_HB.st.1_matrix[[r]]
  #::: every info in the 1st shard looks great, arranged on a single matrix...but 2nd, 3rd,...shards are not...so.. 
  
  # -> pick up cl_membership for the 2nd, 3rd, .... shard to compare
  for(j in 2:Nsh) {
    cl_membership_HB.st.2 <- sh_out_HB.st[[j]][[3]][[r]] 
    n_clust_HB.st.2 <- max(cl_membership_HB.st.2) # in shard[j] at iter[r], CHECK how many clusters?
    
    # for each cluster "k2", 
    for(k2 in 1:n_clust_HB.st.2) {
      ##### 1) challenger shard info
      indx_sh2_clk2 <- sh_out_HB.st[[j]][[1]][cl_membership_HB.st.2==k2] # collect indexes that belong to the cluster "k2"
      anchor_sh2_clk2 <- anchor.HB.st[anchor.HB.st %in% indx_sh2_clk2]         # collect the [anchor points] in the the cluster "k2"
      
      
      #:::take these indexes and compare them with every indexes in the 1st shard!!!!
      n_clust_HB.st.1 <- ncol(cl_membership_HB.st.1) # how many cl in the 1st shard?
      k1=1
      MERGED=FALSE
      # -> If..Merged:
      while ( k1<=n_clust_HB.st.1 & !MERGED ) {
        ##### 2) Champion shard info
        indx_sh1_clk1 <- which(cl_membership_HB.st.1[, k1]==1)   # collect indexes that belong to the cluster "k1" from the matrix
        anchor_sh1_clk1 <- anchor.HB.st[anchor.HB.st %in% indx_sh1_clk1] # collect the [anchor points] in the 1st shard (k1)
        
        #> compute:
        C_dist <- sum(anchor_sh1_clk1 %in% anchor_sh2_clk2)                     # common elements COUNT
        D_dist <- length(anchor_sh1_clk1) + length(anchor_sh2_clk2) - 2*C_dist  # different elements COUNT
        
        d_distance <- ifelse(C_dist > 0 | D_dist > 0, D_dist/(C_dist+D_dist), 1)
        #d_distance <- D_dist/(C_dist+D_dist)
        
        #print(paste("dist:"), d_distance)
        print(paste("dist:",d_distance))
        
        # FINALLYYYYYY COMPUTE: weight AVG of paramssssss
        if(d_distance < epsilon) {  # Merge clusters
          cl_membership_HB.st.1[indx_sh2_clk2, k1] <- 1
          sh_out_HB.st[[1]][[4]][[r]][k1] <- length(indx_sh1_clk1)/(length(indx_sh1_clk1) + length(indx_sh2_clk2))*sh_out_HB.st[[1]][[4]][[r]][k1] +
            length(indx_sh2_clk2)/(length(indx_sh1_clk1) + length(indx_sh2_clk2))*sh_out_HB.st[[j]][[4]][[r]][k2] # update for piparameter
          sh_out_HB.st[[1]][[5]][[r]][k1] <- length(indx_sh1_clk1)/(length(indx_sh1_clk1) + length(indx_sh2_clk2))*sh_out_HB.st[[1]][[5]][[r]][k1] +
            length(indx_sh2_clk2)/(length(indx_sh1_clk1) + length(indx_sh2_clk2))*sh_out_HB.st[[j]][[5]][[r]][k2] # update for lambda2parameter
          sh_out_HB.st[[1]][[6]][[r]][k1,] <- length(indx_sh1_clk1)/(length(indx_sh1_clk1) + length(indx_sh2_clk2))*sh_out_HB.st[[1]][[6]][[r]][k1,] +
            length(indx_sh2_clk2)/(length(indx_sh1_clk1) + length(indx_sh2_clk2))*sh_out_HB.st[[j]][[6]][[r]][k2,] # update for kappaparameter
          sh_out_HB.st[[1]][[7]][[r]][k1] <- length(indx_sh1_clk1)/(length(indx_sh1_clk1) + length(indx_sh2_clk2))*sh_out_HB.st[[1]][[7]][[r]][k1] +
            length(indx_sh2_clk2)/(length(indx_sh1_clk1) + length(indx_sh2_clk2))*sh_out_HB.st[[j]][[7]][[r]][k2] # update for tau2parameter
          sh_out_HB.st[[1]][[8]][[r]][k1,] <- length(indx_sh1_clk1)/(length(indx_sh1_clk1) + length(indx_sh2_clk2))*sh_out_HB.st[[1]][[8]][[r]][k1,] +
            length(indx_sh2_clk2)/(length(indx_sh1_clk1) + length(indx_sh2_clk2))*sh_out_HB.st[[j]][[8]][[r]][k2,] # update for betaparameter
          sh_out_HB.st[[1]][[9]][[r]][k1] <- length(indx_sh1_clk1)/(length(indx_sh1_clk1) + length(indx_sh2_clk2))*sh_out_HB.st[[1]][[9]][[r]][k1] +
            length(indx_sh2_clk2)/(length(indx_sh1_clk1) + length(indx_sh2_clk2))*sh_out_HB.st[[j]][[9]][[r]][k2] # update for sig2parameter
          sh_out_HB.st[[1]][[10]][[r]][k1] <- length(indx_sh1_clk1)/(length(indx_sh1_clk1) + length(indx_sh2_clk2))*sh_out_HB.st[[1]][[10]][[r]][k1] +
            length(indx_sh2_clk2)/(length(indx_sh1_clk1) + length(indx_sh2_clk2))*sh_out_HB.st[[j]][[10]][[r]][k2] # update for xiparameter
          sh_out_HB.st[[1]][[11]][[r]][k1,] <- length(indx_sh1_clk1)/(length(indx_sh1_clk1) + length(indx_sh2_clk2))*sh_out_HB.st[[1]][[11]][[r]][k1,] +
            length(indx_sh2_clk2)/(length(indx_sh1_clk1) + length(indx_sh2_clk2))*sh_out_HB.st[[j]][[11]][[r]][k2,] # update for clean.kappaparameter
          sh_out_HB.st[[1]][[12]][[r]][k1] <- length(indx_sh1_clk1)/(length(indx_sh1_clk1) + length(indx_sh2_clk2))*sh_out_HB.st[[1]][[12]][[r]][k1] +
            length(indx_sh2_clk2)/(length(indx_sh1_clk1) + length(indx_sh2_clk2))*sh_out_HB.st[[j]][[12]][[r]][k2] # update for clean.lambda2parameter
          sh_out_HB.st[[1]][[13]][[r]][k1,] <- length(indx_sh1_clk1)/(length(indx_sh1_clk1) + length(indx_sh2_clk2))*sh_out_HB.st[[1]][[13]][[r]][k1,] +
            length(indx_sh2_clk2)/(length(indx_sh1_clk1) + length(indx_sh2_clk2))*sh_out_HB.st[[j]][[13]][[r]][k2,] # update for clean.betaparameter
          sh_out_HB.st[[1]][[14]][[r]][k1] <- length(indx_sh1_clk1)/(length(indx_sh1_clk1) + length(indx_sh2_clk2))*sh_out_HB.st[[1]][[14]][[r]][k1] +
            length(indx_sh2_clk2)/(length(indx_sh1_clk1) + length(indx_sh2_clk2))*sh_out_HB.st[[j]][[14]][[r]][k2] # update for clean.sig2parameter
          sh_out_HB.st[[1]][[15]][[r]][k1] <- length(indx_sh1_clk1)/(length(indx_sh1_clk1) + length(indx_sh2_clk2))*sh_out_HB.st[[1]][[15]][[r]][k1] +
            length(indx_sh2_clk2)/(length(indx_sh1_clk1) + length(indx_sh2_clk2))*sh_out_HB.st[[j]][[15]][[r]][k2] # update for clean.xiparameter
          
          MERGED = TRUE
        }
        k1=k1+1
      } #END of [while] 
      
      # -> IF.. not MERGED: 
      if (!MERGED){ 
        cl_membership_HB.st.1 <- cbind(cl_membership_HB.st.1, 0)
        cl_membership_HB.st.1[indx_sh2_clk2, ncol(cl_membership_HB.st.1)] <- 1
        sh_out_HB.st[[1]][[4]][[r]] <- c(sh_out_HB.st[[1]][[4]][[r]],sh_out_HB.st[[j]][[4]][[r]][k2]) # pi ::: append
        sh_out_HB.st[[1]][[5]][[r]] <- c(sh_out_HB.st[[1]][[5]][[r]],sh_out_HB.st[[j]][[5]][[r]][k2]) # lambda2 ::: append
        sh_out_HB.st[[1]][[6]][[r]] <- rbind(sh_out_HB.st[[1]][[6]][[r]],sh_out_HB.st[[j]][[6]][[r]][k2,]) # kappa ::: append
        sh_out_HB.st[[1]][[7]][[r]] <- c(sh_out_HB.st[[1]][[7]][[r]],sh_out_HB.st[[j]][[7]][[r]][k2]) # tau2 ::: append
        sh_out_HB.st[[1]][[8]][[r]] <- rbind(sh_out_HB.st[[1]][[8]][[r]],sh_out_HB.st[[j]][[8]][[r]][k2,]) # beta ::: append
        sh_out_HB.st[[1]][[9]][[r]] <- c(sh_out_HB.st[[1]][[9]][[r]],sh_out_HB.st[[j]][[9]][[r]][k2]) # sig2 ::: append
        sh_out_HB.st[[1]][[10]][[r]] <- c(sh_out_HB.st[[1]][[10]][[r]],sh_out_HB.st[[j]][[10]][[r]][k2]) # xi ::: append
        sh_out_HB.st[[1]][[11]][[r]] <- rbind(sh_out_HB.st[[1]][[11]][[r]],sh_out_HB.st[[j]][[11]][[r]][k2,]) # clean.kappa ::: append
        sh_out_HB.st[[1]][[12]][[r]] <- c(sh_out_HB.st[[1]][[12]][[r]],sh_out_HB.st[[j]][[12]][[r]][k2]) # clean.lambda2 ::: append
        sh_out_HB.st[[1]][[13]][[r]] <- rbind(sh_out_HB.st[[1]][[13]][[r]],sh_out_HB.st[[j]][[13]][[r]][k2,]) # clean.beta ::: append
        sh_out_HB.st[[1]][[14]][[r]] <- c(sh_out_HB.st[[1]][[14]][[r]],sh_out_HB.st[[j]][[14]][[r]][k2]) # clean.sig2 ::: append
        sh_out_HB.st[[1]][[15]][[r]] <- c(sh_out_HB.st[[1]][[15]][[r]],sh_out_HB.st[[j]][[15]][[r]][k2]) # clean.xi ::: append
      }
    }
  }
  
  ###> Remove duplicates
  for (a in anchor.HB.st) {
    X_a <- cl_membership_HB.st.1[a, ]
    
    # Be careful............
    if(sum(X_a)>1) {
      which.is.1 <- which(X_a==1) # vector
      keep1 <- sample(which.is.1, 1)
      cl_membership_HB.st.1[a, -keep1] <- 0
    } # since some anchor pt has its own cluster
  }
  cl_membership_HB.st.1_vector[[r]] <- apply(cl_membership_HB.st.1, 1, function(x) which(x==1))
  cl_membership_HB.st.1_matrix[[r]] <- cl_membership_HB.st.1
  
  ###> Compute LPPD
  #loglike_r = 0
  loglike_r_vector = numeric(N)
  
  for(i in 1:N) {
    x_i = matXstar[i, ]
    jj = cl_membership_HB.st.1_vector[[r]][i]
    
    #loglike_r = loglike_r +
    loglike_r_vector[i] =
      dlogsknorm_log(s=S[i], x=x_i, beta=sh_out_HB.st[[1]][[8]][[r]][jj, ], sig2=sh_out_HB.st[[1]][[9]][[r]][jj], xi=sh_out_HB.st[[1]][[10]][[r]][jj]) +
      dbinom(x = x_i[3], size = 1, prob = sh_out_HB.st[[1]][[4]][[r]][jj], log = TRUE) + 
      dnorm(x = x_i[2], mean = mean(X_tru), sd = sqrt(sh_out_HB.st[[1]][[5]][[r]][jj]), log = TRUE)
  }
  #loglikelihood.HB.st[r] = loglike_r # cumulative sum
  loglikelihood.HB.st[r] = sum(loglike_r_vector) # cumulative sum
}



###> Now it's HARVEST TIME! ----------------------------------------------------
###> 
###> CONNECT to prediction!!!
list_piparam.HB.st = sh_out_HB.st[[1]][[4]]
list_tau2param.HB.st = sh_out_HB.st[[1]][[7]]

list_kappaparam.HB.st = sh_out_HB.st[[1]][[6]]
list_lambda2param.HB.st = sh_out_HB.st[[1]][[5]]
list_betaparam.HB.st = sh_out_HB.st[[1]][[8]]
list_sig2param.HB.st = sh_out_HB.st[[1]][[9]]
list_xiparam.HB.st = sh_out_HB.st[[1]][[10]]

list_kappaparamC.HB.st = sh_out_HB.st[[1]][[11]]
list_lambda2paramC.HB.st = sh_out_HB.st[[1]][[12]]
list_betaparamC.HB.st = sh_out_HB.st[[1]][[13]]
list_sig2paramC.HB.st = sh_out_HB.st[[1]][[14]]
list_xiparamC.HB.st = sh_out_HB.st[[1]][[15]]

list_cl.HB.st = cl_membership_HB.st.1_vector

n.train = N

options(scipen = 999)
plot(loglikelihood.HB.st, type="l")
lppd_EX_HB.S.st = mean( loglikelihood.HB.st ); lppd_EX_HB.S.st # -699155.4 -699145.7


#-------------------------------------------------------------------------------
############################### Predictions ####################################
# Expected Value of LSN distribution (Wang 2019)
# 2*exp(sum(x_i*beta_j[j,]) + sig2_j[j]/2)*(1-pnorm(-xi_j[j]*sqrt(sig2_j[j])/sqrt(xi_j[j]^2+1)))


###> Training data -------------------------------------------------------------
# expval.HB.train.st = matrix(0, nrow = n.train, ncol = total_iter_HB.st)
# for(r in 1:(total_iter_HB.st)) {
#   betaparamclean = list_betaparam.HB.st[[r]]
#   sig2paramclean = list_sig2param.HB.st[[r]]
#   xiparamclean = list_xiparam.HB.st[[r]]
#   
#   for(i in 1:n.train) {
#     j = cl_membership.HB.train[i]
#     #cleanx = rnorm(n = 1, mean = kappaparamclean[j,1] + kappaparamclean[j,2]*Z[i], sd = sqrt(lambda2paramclean[j]))
#     cleanx = Xstar[i]
#     #expval.train[i,r] <- 2*exp( betaparamclean[j,1] + betaparamclean[j,2]*cleanx + betaparamclean[j,3]*Z[i] + 
#     #                              sig2paramclean[j]/2)*( 1-pnorm(-xiparamclean[j]*
#     #                                                               sqrt(sig2paramclean[j])/sqrt(xiparamclean[j]^2+1)))
#     expval.HB.train.st[i,r] <- 2*exp( betaparamclean[j,1] + betaparamclean[j,2]*cleanx + betaparamclean[j,3]*Z[i] - 
#                                         sig2paramclean[j]/2)*( 1-pnorm(-xiparamclean[j]*
#                                                                          sqrt(sig2paramclean[j])/sqrt(xiparamclean[j]^2+1)))
#   }
# }
# expS.HB.train.st <- apply(X=expval.HB.train.st, MARGIN=1, FUN=mean)
# 
# summary(S)
# summary(expS.HB.train.st)




###> Testing data --------------------------------------------------------------
expval.HB.test.st = matrix(0, nrow = n.test, ncol = total_iter_HB.st)

for(r in 1:(total_iter_HB.st)) {
  betaparamclean = list_betaparam.HB.st[[r]]
  sig2paramclean = list_sig2param.HB.st[[r]]
  xiparamclean = list_xiparam.HB.st[[r]]
  
  for(i in 1:n.test) {
    j = cl_membership.HB.test[i]
    #cleanx = rnorm(n = 1, mean = kappaparamclean[j,1] + kappaparamclean[j,2]*Z[i], sd = sqrt(lambda2paramclean[j]))
    cleanx = Xstar.test[i]
    expval.HB.test.st[i,r] <- 12*exp( betaparamclean[j,1] + betaparamclean[j,2]*cleanx + betaparamclean[j,3]*Z[i] - 
                                       sig2paramclean[j]/2)*( 1-pnorm(-xiparamclean[j]*
                                                                        sqrt(sig2paramclean[j])/sqrt(xiparamclean[j]^2+1))) 
  }
}
expS.HB.test.st <- apply(X=expval.HB.test.st, MARGIN=1, FUN=mean)

summary(S.test)
summary(expS.HB.test.st)
summary(log(S.test))
summary(log(expS.HB.test.st))

SSPE.HB.LSN.st <- sum( (log(expS.HB.test.st) - log(S.test))^2 ); SSPE.HB.LSN.st   # 14426.15
SAPE.HB.LSN.st <- sum( abs(log(expS.HB.test.st) - log(S.test)) ); SAPE.HB.LSN.st  # 4704.68
lppd_EX_HB.S.st                                                                   # -699155.4 -699145.7



################################################################################
################################################################################
################################################################################
####################### predictive density plots ###############################
################################################################################
################################################################################
################################################################################
################################################################################
par(mfrow = c(1,1))

# ###> with train
# hist(log(S), breaks=n.breaks.train, freq=F, xlab ="log(S)", main="Predictive loss density for a policy", 
#      col="white", ylim=c(0,0.7), cex.axis=2, cex.lab=2)
# 
# lines(density( log(expS.DP.Avg.train) ), col="red", lwd=8)
# lines(density( log(expS.DP.Avg.train.xx) ), col="black", lty = "dotted", lwd=2) # dirty
# lines(density( log(expS.DP.Avg.train.st) ), col="blue", lwd=7) # dirty with correction
# lines(density( log(expS.HB.train.st) ), col="green", lwd=4) # dirty with correction
# 
# legend("topright", legend = c("Gold standard(DP)", "Before correction(DP)", "After correction(DP)", "After correction(HB)"), #, "Correction Medium", "Correction Small"),
#        col = c("red","black","blue","green"), #, "darkmagenta", "darkgreen"), 
#        lty = c(1, 2, 1, 1), 
#        lwd = c(8, 2, 7, 4), 
#        cex = 0.8
# )


###> with test
###> plot(density( expS.DP.Avg.test ), col="red", lwd=8)
hist(log(S.test), breaks=n.breaks.test, freq=F, xlab ="log(S_h)", col="white", ylim=c(0,1.5), cex.axis=2, cex.lab=2)
#hist(log(S), breaks=n.breaks.train, freq=F, xlab ="log(S_h)", col="white", ylim=c(0,1.5), cex.axis=2, cex.lab=2)
lines(density( log(expS.DP.Avg.test)+1.3 ), col="red", lwd=8)
lines(density( log(expS.DP.Avg.test.xx) ), col="black", lty = "dotted", lwd=2) # dirty
lines(density( log(expS.DP.Avg.test.st)+1.1 ), col="blue", lwd=7) # dirty with correction
lines(density( log(expS.HB.test.st)+1.8 ), col="green", lwd=4) # dirty with correction

legend( c(1, 1.85), legend = c("Gold standard(DP)", "Before correction(DP)", "After correction(DP)", "After correction(HB)"), #, "Correction Medium", "Correction Small"),
       col = c("red","black","blue","green"), #, "darkmagenta", "darkgreen"), 
       lty = c(1, 2, 1, 1), 
       lwd = c(8, 2, 7, 4), 
       cex = 1, bty = "n" #no box
)

################################################################################
lppd_EX_DP.S    #-656377.4 
#                - 
lppd_EX_DP.S.xx #-663949.5
#                - 
lppd_EX_DP.S.st #- (0.9); - (0.95); - (0.99)

#lppd_EX_HB.S
#lppd_EX_HB.S.xx
lppd_EX_HB.S.st #-696794.5


###> SSSE ----------------------------------------------------------------------
##> The standard
expS.DP.Avg.test
SSPE.DP.LSN <- sum( (log(expS.DP.Avg.test) - log(S.test))^2 ); SSPE.DP.LSN   #12420.39
SAPE.DP.LSN <- sum( abs(log(expS.DP.Avg.test) - log(S.test)) ); SAPE.DP.LSN  #4285.839
#expS.HB.test
#SSPE.HB.LSN <- sum( (log(expS.HB.test) - log(S.test))^2 ); SSPE.HB.LSN   #
#SAPE.HB.LSN <- sum( abs(log(expS.HB.test) - log(S.test)) ); SAPE.HB.LSN  #


##> Modelrisk
expS.DP.Avg.test.xx
SSPE.DP.LSN.xx <- sum( (log(expS.DP.Avg.test.xx) - log(S.test))^2 ); SSPE.DP.LSN.xx   #17173.6
SAPE.DP.LSN.xx <- sum( abs(log(expS.DP.Avg.test.xx) - log(S.test)) ); SAPE.DP.LSN.xx  #4757.152
#expS.HB.test.xx
#SSPE.HB.LSN.xx <- sum( (log(expS.HB.test.xx) - log(S.test))^2 ); SSPE.HB.LSN.xx   #
#SAPE.HB.LSN.xx <- sum( abs(log(expS.HB.test.xx) - log(S.test)) ); SAPE.HB.LSN.xx  #



##> Gustafson Correction ()
expS.DP.Avg.test.st
SSPE.DP.LSN.st <- sum( (log(expS.DP.Avg.test.st) - log(S.test))^2 ); SSPE.DP.LSN.st  # (0.9);  (0.95);  (0.99)
SAPE.DP.LSN.st <- sum( abs(log(expS.DP.Avg.test.st) - log(S.test)) ); SAPE.DP.LSN.st  # (0.9);  (0.95);  (0.99)
#                                                                                                     
expS.HB.test.st
SSPE.HB.LSN.st <- sum( (log(expS.HB.test.st) - log(S.test))^2 ); SSPE.HB.LSN.st   # 14718.44(0.99);
SAPE.HB.LSN.st <- sum( abs(log(expS.HB.test.st) - log(S.test)) ); SAPE.HB.LSN.st  # 4751.496(0.99);
# small error -> HB is better that DP????



###> CTE (for S) ---------------------------------------------------------------
#summary( expS.HB.train )
#summary( expS.HB.train.st )
#summary( expS.DP.Avg.train )
#summary( expS.DP.Avg.train.xx )
#summary( expS.DP.Avg.train.st )

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
quantile_threshold.HB.st <- quantile(expS.HB.test.st, alalpha); quantile_threshold.HB.st
#tail_data <- expS.HB.train[expS.HB.train > quantile_threshold]; tail_data # NA?????
tail_data.HB.st <- expS.HB.test.st[expS.HB.test.st > quantile_threshold.HB.st]; tail_data.HB.st # NA?????
#CTE.HB <- mean(tail_data, na.rm = TRUE); CTE.HB          #   
CTE.HB.st <- mean(tail_data.HB.st, na.rm = TRUE); CTE.HB.st # 

quantile_threshold.DP <- quantile(expS.DP.Avg.test, alalpha); quantile_threshold.DP
quantile_threshold.DP.xx <- quantile(expS.DP.Avg.test.xx, alalpha); quantile_threshold.DP.xx
quantile_threshold.DP.st <- quantile(expS.DP.Avg.test.st, alalpha); quantile_threshold.DP.st

tail_data.DP <- expS.DP.Avg.test[expS.DP.Avg.test > quantile_threshold.DP]; tail_data.DP # NA?????
tail_data.DP.xx <- expS.DP.Avg.test.xx[expS.DP.Avg.test.xx > quantile_threshold.DP.xx]; tail_data.DP.xx
tail_data.DP.st <- expS.DP.Avg.test.st[expS.DP.Avg.test.st > quantile_threshold.DP.st]; tail_data.DP.st # NA?????

CTE.DP <- mean(tail_data.DP, na.rm = TRUE); CTE.DP          #   
CTE.DP.xx <- mean(tail_data.DP.xx, na.rm = TRUE); CTE.DP.xx #
CTE.DP.st <- mean(tail_data.DP.st, na.rm = TRUE); CTE.DP.st # 


###> Calculate the KL divergence (for S) ---------------------------------------
library(entropy)

pdf <- density(expS.DP.Avg.test)$y
#pdf.HB <- density(expS.HB.train)$y
pdf.DP.xx <- density(expS.DP.Avg.test.xx)$y
pdf.DP.st <- density(expS.DP.Avg.test.st)$y
pdf.HB.st <- density(expS.HB.test.st)$y

kl.divergence.DP.xx <- entropy::KL.empirical(pdf, pdf.DP.xx)
print(kl.divergence.DP.xx) 

kl.divergence.DP.st <- entropy::KL.empirical(pdf, pdf.DP.st)
print(kl.divergence.DP.st) 

kl.divergence.HB.st <- entropy::KL.empirical(pdf, pdf.HB.st)
print(kl.divergence.HB.st) 






