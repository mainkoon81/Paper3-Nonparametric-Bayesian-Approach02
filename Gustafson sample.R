# P.Gustafson


################################################################################
# fits logistic regression without measurement error.
################################################################################



###################################################
### logistic regression with no meas. err       ###
###################################################
logreg0 <- function(d, w, x.oth, 
                    burn=500, size=5000, 
                    s=function(z) {z} ) {
  
  n <- length(d); p <- 2+dim(as.matrix(x.oth))[2]
  
  ##############################
  ### initial value for beta ###
  ##############################
  des.s <- cbind(1,s(w), x.oth)
  tmp <- glm(d~des.s[,-1], family=binomial())
  mn <- tmp$coef; vr <- summary(tmp)$cov.unscaled; jmp.beta <- sqrt(diag(vr)/2)
  vr.in <- solve(vr); vr.ch <- chol(vr)
  beta <- mn+t(vr.ch)%*%rnorm(p)
  
  #############################
  ### other pre-MCMC set-up ###
  #############################
  res <- matrix(0,size,p) 
  acc <- rep(0,n); acc1.beta <- acc2.beta <- acc1.den <- acc2.den <- 0
  
  ##################
  ### start MCMC ###
  ##################
  for (i in (-burn):size) {print(i)  
    
    ################################
    ### update beta,lik. approx, ###
    ################################
    acc1.den <- acc1.den+1
    beta.cand <- mn+t(vr.ch)%*%rnorm(3)
    tmp1 <- des.s %*% beta.cand; tmp2 <- d*tmp1
    logacc1 <- sum(tmp2)-sum(log(1+exp(tmp1))) +
      (0.5)* t(beta.cand-mn) %*% vr.in %*% (beta.cand-mn)
    tmp1 <- des.s %*% beta; tmp2 <- d*tmp1
    logacc2 <- sum(tmp2)-sum(log(1+exp(tmp1))) +
      (0.5)* t(beta-mn) %*% vr.in %*% (beta-mn) 
    if (log(runif(1))<(logacc1-logacc2)) {
      beta <- beta.cand; acc1.beta <- acc1.beta+1 }  
    
    ####################
    ### record stuff ###
    ####################
    if (i>0) {
      res[i,] <- beta } }
  
  ##############
  ### output ###
  ##############
  list(beta=res,
       acc1.beta=acc1.beta/acc1.den) 

}










################################################################################
# ---------- MCMC for logistic regression with measurement error ------------- #  
# supports both normal exposure model and flexible `reverse exposure' model (see Secs. 4.4, 4.5, 4.7 of MEMSE).
################################################################################


se <- function(x) { sqrt(var(x)/length(x)) }

##########################################################

rdirichlet <- function(param) {
  tmp <- rgamma(length(param),param)
  tmp/sum(tmp) }

###########################################################


###################################################
### main function                               ###
### method a, assume x|x.oth is normal          ###
### method b, assume x.oth|x normal, x flexible ###
### note currently only scalar x.oth supported  ###
###################################################
logreg1 <- function(d, w, x.oth, tau, 
                    method="a",
                    burn=500, size=5000, 
                    s=function(z) {z} ) {
  
  n <- length(d); p<-2+dim(as.matrix(x.oth))[2]
  
  #################################################
  ### set up grid, theta, x, gr.pr for method b ###
  #################################################
  if (method=="b") {
    gr.lft <- min(w)-2.5*tau; gr.rht <- max(w)+2.5*tau
    gr.val <- seq(from=gr.lft, to=gr.rht, by=tau/4)
    gr.len <- length(gr.val) 
    
    theta <- rep(0,n)
    for (i in 1:n) {
      theta[i] <- round( gr.len*(w[i]-gr.lft)/(gr.rht-gr.lft) ) }
    theta <- theta+sample( (-6:6), size=n, replace=T)
    ind <- (theta<1);      theta[ind] <- 1-theta[ind]
    ind <- (theta>gr.len); theta[ind] <- gr.len+1-(theta[ind]-gr.len) 
    x <- gr.val[theta]
    
    gr.pr <- rdirichlet(rep(1,gr.len)) }
  
  #############################
  ### set up x for method a ###
  #############################
  if (method=="a") {x <- rnorm(n,w,tau); gr.val <- NULL}
  
  ##############################
  ### initial value for beta ###
  ##############################
  des.s <- cbind(1,s(w), x.oth)
  tmp <- glm(d~des.s[,-1], family=binomial())
  mn <- tmp$coef; vr <- summary(tmp)$cov.unscaled; jmp.beta <- sqrt(diag(vr)/2)
  vr.in <- solve(vr); vr.ch <- chol(vr)
  beta <- mn+t(vr.ch)%*%rnorm(p)
  
  #########################################
  ### initial values for alph, lam2     ###                        
  ### different interps for methods a,b ###
  #########################################
  alph <- c(0,0); lam2 <- 1 
  
  #############################
  ### other pre-MCMC set-up ###
  #############################
  res <- matrix(0,size,p); res.cond <- matrix(0,size,3); res.x <- matrix(0,size,10) 
  res.pr <- NULL; if (method=="b") {res.pr <- matrix(0,size,gr.len)}
  acc <- rep(0,n); acc1.beta <- acc2.beta <- acc1.den <- acc2.den <- 0
  
  ##################
  ### start MCMC ###
  ##################
  for (i in (-burn):size) { print(i) 
    
    ##############################################
    ### update beta, either lik. approx, or RW ###
    ##############################################
    des.s <- cbind(1,s(x),x.oth)
    if (runif(1)<.075) {  
      acc1.den <- acc1.den+1
      tmp <- glm(d~des.s[,-1], family=binomial)
      mn <- tmp$coef; vr <- summary(tmp)$cov.unscaled
      vr.in <- solve(vr); vr.ch <- chol(vr)    
      beta.cand <- mn+t(vr.ch)%*%rnorm(3)
      tmp1 <- des.s %*% beta.cand; tmp2 <- d*tmp1
      logacc1 <- sum(tmp2)-sum(log(1+exp(tmp1))) +
        (0.5)* t(beta.cand-mn) %*% vr.in %*% (beta.cand-mn)
      tmp1 <- des.s %*% beta; tmp2 <- d*tmp1
      logacc2 <- sum(tmp2)-sum(log(1+exp(tmp1))) +
        (0.5)* t(beta-mn) %*% vr.in %*% (beta-mn) 
      if (log(runif(1))<(logacc1-logacc2)) {
        beta <- beta.cand; acc1.beta <- acc1.beta+1 } } 
    else {  
      acc2.den <- acc2.den+1
      beta.cand <- rnorm(3, beta, jmp.beta)
      tmp1 <- des.s %*% beta.cand; tmp2 <- d*tmp1
      logacc1 <- sum(tmp2)-sum(log(1+exp(tmp1)))
      tmp1 <- des.s %*% beta; tmp2 <- d*tmp1
      logacc2 <- sum(tmp2)-sum(log(1+exp(tmp1)))
      if (log(runif(1))<(logacc1-logacc2)) {
        beta <- beta.cand; acc2.beta <- acc2.beta+1 } } 
    
    #############################
    ### update x by method a  ###
    #############################
    if (method=="a") {
      x.cand <- x+rnorm(n,0,.75*tau); des.s.cand <- des.s; des.s.cand[,2] <- s(x.cand)
      tmp.cand <- des.s.cand%*%beta; tmp <- des.s%*%beta
      logacc1 <- (-0.5)*(w-x.cand)^2 / tau^2 +
        d*tmp.cand -  log(1+exp(tmp.cand)) + 
        (-0.5)*(x.cand-alph[1]-alph[2]*x.oth)^2 / lam2 
      logacc2 <- (-0.5)*(w-x)^2 / tau^2 +
        d*tmp -  log(1+exp(tmp)) + 
        (-0.5)*(x-alph[1]-alph[2]*x.oth)^2 / lam2 
      ind <- ( log(runif(n)) < (logacc1-logacc2) )
      x[ind] <- x.cand[ind]
      acc[ind] <- acc[ind]+1 }
    
    ############################
    ### update x by method b ###
    ############################
    if (method=="b") {
      theta.cand <- theta+sample((-6):6 ,size=n, replace=T)
      ind <- (theta.cand<1); theta.cand[ind] <- 1-theta.cand[ind]
      ind <- (theta.cand>gr.len); theta.cand[ind] <- gr.len+1-(theta.cand[ind]-gr.len)
      tmp <- des.s %*% beta
      x.cand <- gr.val[theta.cand]; des.s.cand <- des.s; des.s.cand[,2] <- s(x.cand); tmp.cand <- des.s.cand%*%beta
      logacc1 <- (-0.5)*(w-x.cand)^2 / tau^2 +
        (-0.5)*(x.oth-alph[1]-alph[2]*x.cand)^2 / lam2 +
        d*tmp.cand -  log(1+exp(tmp.cand)) + log(gr.pr[theta.cand])
      logacc2 <- (-0.5)*(w-x)^2 / tau^2 +
        (-0.5)*(x.oth-alph[1]-alph[2]*x)^2 / lam2 +
        d*tmp -  log(1+exp(tmp)) + log(gr.pr[theta])
      ind <- ( log(runif(n)) < (logacc1-logacc2) )
      theta[ind] <- theta.cand[ind]
      acc[ind] <- acc[ind]+1 
      x <- gr.val[theta]; des.s[,2] <- s(x) }
    
    ####################
    ### update gr.pr ###
    ####################
    if (method=="b") {
      hyp <- rep(1,gr.len)
      if (tau>0) {
        tmp <- table(theta)
        nonzero <- as.numeric(names(tmp))
        hyp[nonzero] <- hyp[nonzero]+tmp }
      gr.pr <- rdirichlet(hyp) } 
    
    ############################
    ### update alph and lam2 ###
    ############################
    if (method=="a") {
      desmat <- cbind(1,x.oth)
      tmp1 <- solve(t(desmat)%*%desmat); tmp2 <- chol(tmp1)
      alph <- tmp1%*%t(desmat)%*%x + sqrt(lam2)*(t(tmp2) %*% rnorm(2))
      lam2 <- (.5*sum((x-desmat%*%alph)^2)+.5)/rgamma(1, n/2+.5) }
    if (method=="b") {
      desmat <- cbind(1,x)
      tmp1 <- solve(t(desmat)%*%desmat); tmp2 <- chol(tmp1)
      alph <- tmp1%*%t(desmat)%*%x.oth + sqrt(lam2)*(t(tmp2) %*% rnorm(2))
      lam2 <- (.5*sum((x.oth-desmat%*%alph)^2)+.5)/rgamma(1, n/2+.5) }
    
    ####################
    ### record stuff ###
    ####################
    if (i>0) {
      res[i,] <- beta
      res.x[i,] <- x[1:10]
      if (method=="b") {res.pr[i,] <- gr.pr}
      res.cond[i,] <- c(alph,sqrt(lam2)) } }  
  
  ##############
  ### output ###
  ##############
  
  if (method=="b") {res.pr <- apply(res.pr,2,mean) }
  
  list(beta=res,
       x=res.x,
       pr=res.pr,
       gr=gr.val,
       cond=res.cond,
       acc.x=acc/(burn+size), 
       acc1.beta=acc1.beta/acc1.den,
       acc2.beta=acc2.beta/acc2.den ) 
}








