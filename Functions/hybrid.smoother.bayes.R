hybrid.smoother.bayes <- function(x,y,n.chains = 4, chain.length=1200, burn.in.percent = 1/6, lambda2.mn0 = 10,
                       lambda2.var0 = 10, phi2.mn0 = .00005,phi2.var0 = .001, s2.0=.00001, buffer.sd=0, at.least=0,construct.from="center"){
  
  init.time <- proc.time()
  
  #Needed Libraries
  library(fields)
  library(statmod)
  library(sns)
  
  X <- cbind(1,x)
  n <- length(x)
  y.mean <- mean(y)
  y <- y - y.mean
  
  # Covariance of GP
  cardinalX<- cbind(quantile( x, probs=c(.25,.75)) )
  Sigma <- Tps.cov(x,x, cardinalX = cardinalX)
  Sigma.inv <- solve(Sigma)
  
  # Basis Functions for Rough Section
  if(construct.from=="left"){
    #(Left to Right)
    D <- matrix(0,nrow=(n),ncol=(n))
    D[lower.tri(D,diag=TRUE)] <- 1
    D <- D[,-1]
  }
  
  if(construct.from=="center"){
    #(from center)
    m <- round(n/2)
    D <- matrix(0,nrow=(n),ncol=(n))
    D[1:m,1:m][upper.tri(D[1:m,1:m],diag=TRUE)] <- -1
    D1 <- matrix(0,nrow = n-m, ncol= n-m)
    D1[1:(n-m),1:(n-m)][lower.tri(D[1:(n-m),1:(n-m)],diag=TRUE)] <- 1
    D[(m+1):n,(m+1):n] <- D1
    D <- D[,-m]
  }
  gamma.length <- dim(D)[2]
  
  ###
  ### Speeding Up Computation
  ###
  tXX <- t(X)%*%X
  X.block <- solve(tXX)%*%t(X)
  PX <- X%*%X.block
  IPX <- diag(1,n)-PX
  G <- IPX%*%D
  tGG <- t(G)%*%G
  C <- diag(.001,n)[,-n]
  Gplus <- G + C
  Gplus.block <- solve(t(Gplus)%*%Gplus)%*%t(Gplus)
  PGplus <- G%*%Gplus.block
  PGplus.2 <- C%*%Gplus.block
  IPGplus <- diag(1,n)-PGplus
  K <- IPGplus %*% IPX + PGplus.2%*%(IPX)
  tKK <- t(K)%*%K
  J <- Gplus.block%*%IPX
  
  ###
  ### Subroutines
  ###
  
  invgammastrt <- function(igmn,igvar){
    q <- 2+(igmn^2)/igvar
    r <- 1/(igmn*(q-1))
    return(list(r=r,q=q))
  }
  
  gammastrt <- function(gmn,gvar){
    beta <- gmn/gvar
    alpha <- gmn*beta
    return(list(alpha=alpha,beta=beta))
  }
  
  sampleMVG <- function(A,b,n){
    tmp.chol <- chol(A)
    mvg.sample <- backsolve(tmp.chol,backsolve(tmp.chol,b,transpose=TRUE)+rnorm(n))
    return(mvg.sample)
  }
  
  ###
  ### Priors
  ###
  
  ### s2 ~ uniform
  
  ### beta ~ uniform
  
  ### lambda ~ gamma
  lambda2.prior <- gammastrt(lambda2.mn0,lambda2.var0)
  lambda2.alpha <- lambda2.prior$alpha
  lambda2.beta <- lambda2.prior$beta
  
  ### phi2 ~ inv-gamma
  phi2.prior <- invgammastrt(phi2.mn0,phi2.var0)
  phi2.r0 <- phi2.prior$r
  phi2.q0 <- phi2.prior$q
  
  
  ### gamma ~ multivariate normal, mean 0, covariance s2 * D.tau
  ### tau2's ~ independent normal distributions mean 0, var = 1/lambda^2
  
  ###
  ### MCMC
  ###
  
  chains <- vector(mode = "list", length = n.chains)
  parameters.to.check <- c("s2","lambda","phi2","beta.star0","beta.star1","gamma.star20","tau2.20","GP20")
  eff.ratio.results <- as.data.frame(matrix(NA,nrow=length(parameters.to.check),ncol=n.chains+1))
  names(eff.ratio.results) <- c("Parameter",1:n.chains)
  eff.ratio.results[,1] <- parameters.to.check
  eff.sample.size.results <- eff.ratio.results
  rhat.results <- eff.ratio.results
  
  for ( l in 1:n.chains){
    
    if(l==1){
      n.mcmc <- chain.length
      burn.in <- round(burn.in.percent*n.mcmc)
      N <- n.mcmc - burn.in
      
      
      ###
      ### Initial values
      ###
      
      ### lambda
      lambda <- sqrt(lambda2.mn0)
      
      ### GP
      f <- rnorm(n)
      Kf <- K%*%f
      Jf <- J%*%f
      
      ### gamma:  zeros except of top five percent of differences
      gamma.star <- diff(y) #c(diff(y))
      for (j in 1:(gamma.length)){
        gamma.star[j] <- max(c(0,(abs(gamma.star[j])-quantile(abs(gamma.star),.99))))*sign(gamma.star[j])
      }
      H.star <- G%*%gamma.star
      gamma <- gamma.star-Jf
      
      ### s2 
      s2 <- s2.0
      
      ### phi2 
      phi2 <- phi2.mn0
      
      ### F:  matrix of tau^2's
      tau2s <- rep(1,gamma.length)
      tau2s.inv <- c(1/rep(1,gamma.length))
      F.inv <- diag(tau2s.inv)
      F <- diag(tau2s)
      
      ### beta
      B.star <- lm(y-H.star ~ X-1)$coefficients
      trend.star <- X%*%B.star
    }
    
    
    
    if(l>1){
      n.mcmc <- chain.length
      burn.in <- round(burn.in.percent*n.mcmc)
      n.mcmc <- n.mcmc - burn.in
      burn.in <- 0
      N <- n.mcmc - burn.in
      
      
      ###
      ### Initial values
      ###
      
      ### lambda
      lambda <- mean(lambda.results)
      
      ### gamma:  zeros except of top five percent of differences
      gamma.star <- apply(gamma.star.results,1,mean)
      gamma <- apply(gamma.results,1,mean)
      
      ### s2 
      s2 <- mean(s2.results)
      
      ### phi2 
      phi2 <- mean(phi2.results)
      
      
      ### beta
      B.star <- apply(beta.star.results,1,mean)
      trend.star <- X%*%B.star
      
    }
    
    
    ### Initiate Matrices to store results
    gamma.star.results <- matrix(NA,gamma.length,n.mcmc-burn.in)
    gamma.results <- matrix(NA,gamma.length,n.mcmc-burn.in)
    s2.results <- rep(NA,n.mcmc-burn.in)
    phi2.results <- rep(NA,n.mcmc-burn.in)
    lambda.results <- rep(NA,n.mcmc-burn.in)
    tau2.results <- matrix(NA,gamma.length,n.mcmc-burn.in)
    beta.star.results <- matrix(NA,2,n.mcmc-burn.in)
    beta.results <- matrix(NA,2,n.mcmc-burn.in)
    GP.results <- matrix(NA,n,n.mcmc-burn.in)
    pred.H.results <- matrix(NA,n,n.mcmc-burn.in)
    pred.trend.results <- matrix(NA,n,n.mcmc-burn.in)
    pred.y.results <- matrix(NA,n,n.mcmc-burn.in)
    
    
    ### Speed Computation
    temp.q.s2 <- n/2 + (gamma.length)/2
    temp.q.phi2 <- n/2 + phi2.q0
    temp.alpha.lambda <- lambda2.alpha + gamma.length
    
    for (k in 1:n.mcmc){
      
      ###
      ### Update s2
      ###
      
      #temp.q <- see above
      temp.r <- t(y-trend.star-H.star-Kf)%*%(y-trend.star-H.star-Kf)/2 + (t(gamma)*tau2s.inv)%*%(gamma)/2
      s2 <- 1/rgamma(1,temp.q.s2,temp.r)
      

      
      ###
      ### Update phi2
      ###
      
      #temp.q <- see above
      temp.r <- t(f)%*%Sigma.inv%*%f/2 + 1/phi2.r0
      phi2 <- 1/rgamma(1,temp.q.phi2,temp.r)
      
      
      
      ### 
      ### Update f
      ###
      
      #tJF.inv <- t(J)%*%F.inv
      tJF.inv <- t(tau2s.inv*J)
      
      temp.A <- Sigma.inv/phi2 + tKK/s2 + tJF.inv%*%J/s2
      temp.b <- t(K)%*%(y-trend.star-H.star)/s2 + tJF.inv%*%gamma.star/s2
      f <- sampleMVG(temp.A,temp.b,n)
      Kf <- K%*%f
      Jf <- J%*%f
      
      
      
      ###
      ### Update lambda2 and lambda 
      ###
      
      #temp.alpha <- see above
      temp.beta <- sum(tau2s)/2 + lambda2.beta
      lambda2 <- rgamma(1,temp.alpha.lambda,temp.beta)
      lambda <- sqrt(lambda2)
      
      
      
      ### 
      ### Update gamma.star
      ###
      
      temp.A <- tGG/s2 + F.inv/s2
      temp.b <- t(G)%*%(y-trend.star-Kf)/s2 + tau2s.inv*Jf/s2
      gamma.star <- sampleMVG(temp.A,temp.b,gamma.length)
      
      H.star <- G%*%gamma.star 
      gamma <- gamma.star - Jf
      
      
      
      ### 
      ### Update beta.star
      ###
      
      temp.A <- tXX/s2
      temp.b <- t(X)%*%(y-H.star-Kf)/s2
      B.star <- sampleMVG(temp.A,temp.b,2)
      
      trend.star <- X%*%B.star
      
      
      
      ###
      ### Update tau2's in F
      ###
      
      temp.mean <- sqrt(s2*lambda^2/(gamma)^2)
      temp.shape <- lambda^2
      #tau2s.inv <- rinvgauss(gamma.length,temp.mean,,temp.shape)
      tau2s.inv <- rinvgauss(gamma.length,temp.mean,temp.shape)
      tau2s <- 1/tau2s.inv
      F.inv <- diag(tau2s.inv)
      F <- diag(tau2s)

      
      
      
      ###
      ### Predicted y
      ###
      
      pred.y <- trend.star + H.star + Kf + y.mean
      
      
      ###
      ### Undo Reparametrization
      ###
      
      H <- D%*%gamma
      B <- B.star - X.block%*%(H+f)
      trend <- X%*%B + y.mean
      
      
      ###
      ### Record Results
      ###
      
      if (k > burn.in){
        s2.results[k-burn.in] <- s2
        phi2.results[k-burn.in] <- phi2
        lambda.results[k-burn.in] <- lambda
        gamma.star.results[,k-burn.in] <- gamma.star
        beta.star.results[,k-burn.in] <- B.star
        tau2.results[,k-burn.in] <- diag(F)
        GP.results[,k-burn.in] <- f 
        gamma.results[,k-burn.in] <- gamma
        beta.results[,k-burn.in] <- B
        pred.H.results[,k-burn.in] <- H
        pred.trend.results[,k-burn.in] <- trend
        pred.y.results[,k-burn.in] <- pred.y
      }
      
      #End of each set of updates
    }
    
    # ###
    # ### Effectiveness of Mixing
    # ###
    # 
    # ### Effect Sample Size and Ratio
    # eff.sample.size.results[1,l+1] <- ess(s2.results)
    # eff.sample.size.results[3,l+1] <- ess(phi2.results)
    # eff.sample.size.results[2,l+1] <- ess(lambda.results)
    # eff.sample.size.results[4,l+1] <- ess(beta.star.results[1,])
    # eff.sample.size.results[5,l+1] <- ess(beta.star.results[2,])
    # eff.sample.size.results[6,l+1] <- ess(gamma.star.results[20,])
    # eff.sample.size.results[7,l+1] <- ess(tau2.results[20,])
    # eff.sample.size.results[8,l+1] <- ess(GP.results[20,])
    # 
    # eff.ratio.results[1:length(parameters.to.check),2:(n.chains+1)] <- eff.sample.size.results[1:length(parameters.to.check),2:(n.chains+1)]/N
    # 
    # 
    
    ###
    ### Combine Results
    ###
    chain <- list(gamma.star.results = gamma.star.results,
                  s2.results = s2.results,
                  phi2.results = phi2.results,
                  lambda.results = lambda.results,
                  beta.star.0.results = beta.star.results[1,],
                  beta.star.1.results = beta.star.results[2,],
                  pred.H.results = pred.H.results,
                  tau2.results = tau2.results,
                  pred.y.results = pred.y.results,
                  GP.results = GP.results,
                  beta0.results = beta.results[1,],
                  beta1.results = beta.results[2,],
                  gamma.results = gamma.results,
                  pred.trend.results = pred.trend.results)
    
    chains[[l]] <- chain
    
    #End of all chains
  }
  
  ###
  ### Combined chains
  ###
  
  all.gamma.star.results <- c()
  all.pred.H.results <- c()
  all.s2.results <- c()
  all.phi2.results <- c()
  all.beta.star.0.results <- c()
  all.beta.star.1.results <- c()
  all.lambda.results <- c()
  all.tau2.results <- c()
  all.pred.y.results <- c()
  all.GP.results <- c()
  all.gamma.results <- c()
  all.beta0.results <- c()
  all.beta1.results <- c()
  all.pred.trend.results <- c()
  
  for (l in 1:n.chains){
    all.gamma.star.results <- cbind(all.gamma.star.results,(chains[[l]])[[1]])
    all.s2.results <- c(all.s2.results,(chains[[l]])[[2]])
    all.phi2.results <- c(all.phi2.results,(chains[[l]])[[3]])
    all.lambda.results <- c(all.lambda.results,(chains[[l]])[[4]])
    all.beta.star.0.results <- c(all.beta.star.0.results,(chains[[l]])[[5]])
    all.beta.star.1.results <- c(all.beta.star.1.results,(chains[[l]])[[6]])
    all.pred.H.results <- cbind(all.pred.H.results,(chains[[l]])[[7]])
    all.tau2.results <- cbind(all.tau2.results,(chains[[l]])[[8]])
    all.pred.y.results <- cbind(all.pred.y.results,(chains[[l]])[[9]])
    all.GP.results <- cbind(all.GP.results,(chains[[l]])[[10]])
    all.beta0.results <- c(all.beta0.results,(chains[[l]])[[11]])
    all.beta1.results <- c(all.beta1.results,(chains[[l]])[[12]])
    all.gamma.results <- cbind(all.gamma.results,(chains[[l]])[[13]])
    all.pred.trend.results <- cbind(all.pred.trend.results,(chains[[l]])[[14]])
  }
  
  # 
  # 
  # ###
  # ### Chain Comparisons
  # ###
  # 
  # ###
  # ### Mixing 
  # ###
  # 
  # color = c("red", "orange","yellow","blue")
  # plot(NULL,xlim=c(1,N),ylim=c(min(all.s2.results),max(all.s2.results)),ylab="s2")
  # for (i in 1:n.chains){
  #   indices <- (1:N) + (i-1)*N
  #   points(all.s2.results[indices],type="l",col=color[i])
  # }
  # 
  # plot(NULL,xlim=c(1,N),ylim=c(min(all.phi2.results),max(all.phi2.results)),ylab="phi2")
  # for (i in 1:n.chains){
  #   indices <- (1:N) + (i-1)*N
  #   points(all.phi2.results[indices],type="l",col=color[i])
  # }
  # 
  # plot(NULL,xlim=c(1,N),ylim=c(min(all.lambda.results),max(all.lambda.results)),ylab="lambda")
  # for (i in 1:n.chains){
  #   indices <- (1:N) + (i-1)*N
  #   points(all.lambda.results[indices],type="l",col=color[i])
  # }
  # 
  # plot(NULL,xlim=c(1,N),ylim=c(min(all.beta.star.0.results),max(all.beta.star.0.results)),ylab="beta.star.0")
  # for (i in 1:n.chains){
  #   indices <- (1:N) + (i-1)*N
  #   points(all.beta.star.0.results[indices],type="l",col=color[i])
  # }
  # 
  # plot(NULL,xlim=c(1,N),ylim=c(min(all.beta.star.1.results),max(all.beta.star.1.results)),ylab="beta.star.1")
  # for (i in 1:n.chains){
  #   indices <- (1:N) + (i-1)*N
  #   points(all.beta.star.1.results[indices],type="l",col=color[i])
  # }
  # 
  # 
  # 
  # ###
  # ### Effectiveness of Mixing:  Rhat
  # ###
  # 
  # ### s2
  # combined.var <- var(all.s2.results)
  # for (i in 1:n.chains){
  #   indices <- (1:N) + (i-1)*N
  #   rhat <- sqrt(combined.var/var(all.s2.results[indices]))
  #   rhat.results[1,1+i] <- rhat
  # }
  # 
  # ### lambda
  # combined.var <- var(all.lambda.results)
  # for (i in 1:n.chains){
  #   indices <- (1:N) + (i-1)*N
  #   rhat <- sqrt(combined.var/var(all.lambda.results[indices]))
  #   rhat.results[2,1+i] <- rhat
  # }
  # 
  # ### phi2
  # combined.var <- var(all.phi2.results)
  # for (i in 1:n.chains){
  #   indices <- (1:N) + (i-1)*N
  #   rhat <- sqrt(combined.var/var(all.phi2.results[indices]))
  #   rhat.results[3,1+i] <- rhat
  # }
  # 
  # ### beta.star.0
  # combined.var <- var(all.beta.star.0.results)
  # for (i in 1:n.chains){
  #   indices <- (1:N) + (i-1)*N
  #   rhat <- sqrt(combined.var/var(all.beta.star.0.results[indices]))
  #   rhat.results[4,1+i] <- rhat
  # }
  # 
  # ### beta.star.1
  # combined.var <- var(all.beta.star.1.results)
  # for (i in 1:n.chains){
  #   indices <- (1:N) + (i-1)*N
  #   rhat <- sqrt(combined.var/var(all.beta.star.1.results[indices]))
  #   rhat.results[5,1+i] <- rhat
  # }
  # 
  # ### gamma20
  # combined.var <- var(all.gamma.star.results[20,])
  # for (i in 1:n.chains){
  #   indices <- (1:N) + (i-1)*N
  #   rhat <- sqrt(combined.var/var(all.gamma.star.results[20,indices]))
  #   rhat.results[6,1+i] <- rhat
  # }
  # 
  # ### tau2.20
  # combined.var <- var(all.tau2.results[20,])
  # for (i in 1:n.chains){
  #   indices <- (1:N) + (i-1)*N
  #   rhat <- sqrt(combined.var/var(all.tau2.results[20,indices]))
  #   rhat.results[7,1+i] <- rhat
  # }
  # 
  # ### GP.20
  # combined.var <- var(all.GP.results[20,])
  # for (i in 1:n.chains){
  #   indices <- (1:N) + (i-1)*N
  #   rhat <- sqrt(combined.var/var(all.GP.results[20,indices]))
  #   rhat.results[8,1+i] <- rhat
  # }
  # 
  # 
  # 
  # ###
  # ### Chain Densities
  # ###
  # 
  # ### s2
  # hist(all.s2.results,prob=TRUE,main = "s2 Distribution & Individual Chains")
  # for (i in 1:n.chains){
  #   s2.density <- density(chains[[i]]$s2.results)
  #   dx <- s2.density$x[2]-s2.density$x[1]
  #   k <- sum(s2.density$y*dx)
  #   points(s2.density$x,s2.density$y/k,type="l",col=color[i])
  # }
  # 
  # ### lambda
  # hist(all.lambda.results,prob=TRUE,main = "Lambda Distribution & Individual Chains")
  # for (i in 1:n.chains){
  #   lambda.density <- density(chains[[i]]$lambda.results)
  #   dx <- lambda.density$x[2]-lambda.density$x[1]
  #   k <- sum(lambda.density$y*dx)
  #   points(lambda.density$x,lambda.density$y/k,type="l",col=color[i])
  # }
  # 
  # ### phi2
  # hist(all.phi2.results,prob=TRUE,main = "phi2 Distribution & Individual Chains")
  # for (i in 1:n.chains){
  #   phi2.density <- density(chains[[i]]$phi2.results)
  #   dx <- phi2.density$x[2]-phi2.density$x[1]
  #   k <- sum(phi2.density$y*dx)
  #   points(phi2.density$x,phi2.density$y/k,type="l",col=color[i])
  # }
  # 
  # ### beta.star.0
  # hist(all.beta.star.0.results,prob=TRUE,main = "beta.star.0 Distribution & Individual Chains")
  # for (i in 1:n.chains){
  #   beta.star.0.density <- density(chains[[i]]$beta.star.0.results)
  #   dx <- beta.star.0.density$x[2]-beta.star.0.density$x[1]
  #   k <- sum(beta.star.0.density$y*dx)
  #   points(beta.star.0.density$x,beta.star.0.density$y/k,type="l",col=color[i])
  # }
  # 
  # ### beta.star.1
  # hist(all.beta.star.1.results,prob=TRUE,main = "beta.star.1 Distribution & Individual Chains")
  # for (i in 1:n.chains){
  #   beta.star.1.density <- density(chains[[i]]$beta.star.1.results)
  #   dx <- beta.star.1.density$x[2]-beta.star.1.density$x[1]
  #   k <- sum(beta.star.1.density$y*dx)
  #   points(beta.star.1.density$x,beta.star.1.density$y/k,type="l",col=color[i])
  # }
  
  
  
  ###
  ### Change Points
  ### 
  
  # quantiles on gamma 
  quantiles <- t(apply(all.gamma.results,1,quantile,c(.025,.975)))
  quantiles <- cbind(1:(gamma.length),quantiles)
  gamma.with.indices <- cbind(1:(gamma.length),all.gamma.results)
 
  
  # change points  
  change.points <- quantiles[!(quantiles[,2] < 0 & quantiles[,3] > 0),1]
  change.points.trimmed <- quantiles[ ((quantiles[,2] > buffer.sd*sqrt(mean(all.s2.results)) & quantiles[,3]> buffer.sd*sqrt(mean(all.s2.results))) | (quantiles[,2]< -buffer.sd*sqrt(mean(all.s2.results)) & quantiles[,3]< -buffer.sd*sqrt(mean(all.s2.results)))) & abs(rowMeans(all.gamma.results))>at.least,1 ]
  
  # ###
  # ### Separating out model parts
  # ###
  # 
  # plot(x,trend,type="l")
  # for (i in 1:50){
  #   index <- sample(1:(n.chains*(n.mcmc-burn.in)),1)
  #   points(x,all.pred.trend.results[,index],type="l")
  # }
  # plot(x,f,type="l",ylim = c(-2,2))
  # for (i in 1:50){
  #   index <- sample(1:(n.chains*(n.mcmc-burn.in)),1)
  #   points(x,all.GP.results[,index],type="l")
  # }
  # 
  # plot(x,trend+f,type="l")
  # for (i in 1:50){
  #   index <- sample(1:(n.chains*(n.mcmc-burn.in)),1)
  #   points(x,all.GP.results[,index]+all.pred.trend.results[,index],type="l")
  # }
  # 
  # plot(x,H,type="l",ylim = c(-3,3))
  # for (i in 1:50){
  #   index <- sample(1:(n.chains*(n.mcmc-burn.in)),1)
  #   points(x,all.pred.H.results[,index],type="l")
  # }
  # plot(x,pred.y,type="l")
  # for (i in 1:50){
  #   index <- sample(1:(n.chains*(n.mcmc-burn.in)),1)
  #   points(x,all.pred.y.results[,index],type="l")
  # }
  # points(x,y,cex=.05,col="pink")
  # abline(v=change.points)
  # 
  # plot(1:(gamma.length),quantiles[,2],cex=.2,ylim=c(min(quantiles[,2]),max(quantiles[,3])),ylab="Gamma",xlab="Time Index",main = "Identifying Change Points - 95% Credible Intervals")
  #   points(1:(gamma.length), quantiles[,3],cex=.2)
  #   apply(X=quantiles,1,function(X) lines(c(X[1],X[1]),c(X[2],X[3])))
  #   apply(X=gamma.with.indices,1,function(X) points(X[1],median(X[2:length(X)]),pch=4,col="blue"))
  #   abline(a=0,b=0,col="red")
  #   abline(a=-qnorm(.025/n)*sqrt(quantile(all.s2.results,.95)),b=0,col="red",lty=2)
  #   abline(a=qnorm(.025/n)*sqrt(quantile(all.s2.results,.95)),b=0,col="red",lty=2)

  
  time <- (proc.time()-init.time)[3]
  
  # return(list(chains=chains,

  #             all.gamma.results=all.gamma.results,
  #             change.points = change.points,
  #             change.points.trimmed = change.points.trimmed,
  #             rhat.results = rhat.results,
  #             eff.ratio.results = eff.ratio.results,
  #             eff.sample.size.results = eff.sample.size.results,
  #             time = time))
  
  
  return(list(s = sqrt(mean(all.s2.results)),
              pred.y.results = apply(pred.y.results,1,mean),
              pred.trend.results = apply(all.pred.trend.results,1,mean)+apply(all.GP.results,1,mean),
              pred.trend.chain = all.pred.trend.results+all.GP.results,
              gammas = apply(all.gamma.results,1,mean),
              gammas.chain = all.gamma.results,
              rough.basis = D,
              lambda.results=all.lambda.results,
              phi2.results=all.phi2.results,
              change.points = change.points,
              change.points.trimmed = change.points.trimmed,
              s = mean(all.s2.results),
              time = time))

}
