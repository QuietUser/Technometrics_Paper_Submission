hybrid.smoother.freq <- function(x,y,lambda.grid,omega.grid,tolerance = 1e-3,initial.gamma.guess=rep(0,length(x)-1),method="AICc",construct.from="center",at.least=0){
  init.time <- proc.time()
  
  library(splines)
  library(RSpectra)

  naturalSplineBasis <- function(sGrid,
                                 sKnots,
                                 degree = 3,
                                 derivative = 0) {
    boundaryKnots<- c( min(sKnots),max(sKnots))
    sKnots0<- c( rep( boundaryKnots[1],degree),sort(sKnots),
                 rep( boundaryKnots[2],degree) )
    testRight<- sGrid < min(sKnots)
    testLeft <- sGrid > max(sKnots)
    if( any(testRight |testLeft) )
    {stop("some points for evaluation outside knot range.")}

    basis <- splineDesign(sKnots0, sGrid,
                          ord= degree+1, outer.ok=TRUE,
                          derivs=derivative)
    # set up constraints to enforce natural BCs.
    const <- splineDesign(sKnots0, boundaryKnots, ord = degree+1,
                          derivs = c(2,2))
    qr.const <- qr(t(const))
    QBasis<- t(qr.qty( qr.const, t(basis) ))
    basis <- QBasis[,-(1:2)]
    basis

    return( basis )
  }

  naturalCubicSplineR<- function( sKnots){
    PhiKnots<- naturalSplineBasis( sKnots,sKnots,
                                   degree=3,
                                   derivative=2)
    K<- length( sKnots)
    Y<- PhiKnots[1:(K-1),]
    DELTA<- (PhiKnots[2:K,] - PhiKnots[1:(K-1),])
    h<- diff( sKnots)

    R<- (t(Y)%*%(h*Y) +
           (t(DELTA)%*%(h*Y) + t(Y)%*%(h*DELTA))/2 +
           t(DELTA)%*%(h*DELTA)/3)
    return(R)
  }

  gamma.FISTA <- function(gamma.coefs,tDV,tDVD,y,tau,lambda,tolerance){
    #For stopping threshold comparison
    gamma.coefs.old <- gamma.coefs

    #Initialize values for extra FISTA step
    t1 <- 1
    gamma.coefs.old.x <- gamma.coefs

    tau.tDVD <- tau * tDVD
    tau.tDVy <- tau * tDV%*%y


    for (k in 1:5000){

      #Gradient descent
      gamma.coefs <- gamma.coefs - tau.tDVD%*%gamma.coefs + tau.tDVy

      #Soft thresholding
      for (j in 1:p){
        gamma.coefs[j] <- max(c(0,(abs(gamma.coefs[j])-lambda*tau/2)))*sign(gamma.coefs[j])
      }

      #Momentum Jump (Extra FISTA step)
      gamma.coefs.new.x <- gamma.coefs
      t2 <- (1+sqrt(1+4*t1^2))/2
      gamma.coefs <- gamma.coefs.new.x + (t1-1)/t2 * (gamma.coefs.new.x - gamma.coefs.old.x)
      t1 <- t2
      gamma.coefs.old.x <- gamma.coefs.new.x

      #Break when tolerance is met
      # if (norm(gamma.coefs-gamma.coefs.old,type="2") < tolerance){
      #   break
      # }
      if (prod(gamma.coefs-gamma.coefs.old < tolerance)==1){
        break
      }

      gamma.coefs.old <- gamma.coefs
    }
    gamma.coefs
  }
  
  
  #Create Basis and R
  knots <- x
  basis <- naturalSplineBasis(knots, knots)
  der.basis <- naturalSplineBasis(knots,
                        knots,
                        degree = 3,
                        derivative = 1)
  tbb <- t(basis)%*%basis
  R <- naturalCubicSplineR(knots)
  n <- length(x)
  
  
  #Initialize Storage
  sparsity.0norm <- matrix(NA, nrow = length(lambda.grid), ncol = length(omega.grid))
  RMSE <- matrix(NA, nrow = length(lambda.grid), ncol = length(omega.grid))
  edf <- matrix(NA, nrow = length(lambda.grid), ncol = length(omega.grid))
  AIC <- matrix(NA, nrow = length(lambda.grid), ncol = length(omega.grid))
  AICc <- matrix(NA, nrow = length(lambda.grid), ncol = length(omega.grid))
  BIC <- matrix(NA, nrow = length(lambda.grid), ncol = length(omega.grid))

  if(construct.from=="left"){
    #Basis functions for rough portion (Left to Right)
    D <- matrix(0,nrow=(n),ncol=(n))
    D[lower.tri(D,diag=TRUE)] <- 1
    D <- D[,-1]
  }
  
  if(construct.from=="center"){
    #Basis functions for rough portion (from center)
    m <- round(n/2)
    D <- matrix(0,nrow=(n),ncol=(n))
    D[1:m,1:m][upper.tri(D[1:m,1:m],diag=TRUE)] <- -1
    D1 <- matrix(0,nrow = n-m, ncol= n-m)
    D1[1:(n-m),1:(n-m)][lower.tri(D[1:(n-m),1:(n-m)],diag=TRUE)] <- 1
    D[(m+1):n,(m+1):n] <- D1
    D <- D[,-m]
  }

  if(construct.from=="slope"){
    D <- matrix(0,nrow=(n),ncol=(n))
    D[lower.tri(D,diag=TRUE)] <- 1
    D <- D%*%D
    D <- D[,-1]
  }
  
  p <- dim(D)[2]
  
  #Prep for finding S - Simultaneous Diagonalization
  eig.tbb <- eigen(tbb, symmetric=TRUE)
  inv.sym.sqrt.tbb <-
    eig.tbb$vectors%*%diag(1/sqrt(eig.tbb$values))%*%
    t(eig.tbb$vectors)
  B <- t(inv.sym.sqrt.tbb)%*%R%*%inv.sym.sqrt.tbb
  eig.B <- eigen(B, symmetric=TRUE)
  Q <- eig.B$vectors
  D.eig <- eig.B$values
  left <- inv.sym.sqrt.tbb%*%Q
  
  #initial guess from gamma
  initial.gamma.guess <- rep(0,p)
  y.diff <- diff(y)
  top.1.percent <- quantile(abs(y.diff),.995)
  initial.gamma.guess[abs(y.diff)>top.1.percent] <- abs(y.diff)[abs(y.diff)>top.1.percent]
  
  #Find best values of lambda and omega
  
  for (w.j in 1: length(omega.grid)){
    w <- omega.grid[w.j]
    H<-  t(left)*sqrt(1/(1+D.eig*w)) 
    S <-  t(H) %*% H 
    # V <- (diag(rep(1,n)) - basis%*%S%*%t(basis))
    left.half <-  basis%*%t(H)
    V <- (diag(rep(1,n)) - left.half%*%t(left.half))
    tDV <- t(D)%*%V
    tDVD <- tDV%*%D
    
    #Step.size
    tau <- 1/eigs_sym(tDVD,1)$values
    
    
    for (l.i in 1:length(lambda.grid)){
      lambda <- lambda.grid[l.i]
      
      #Initialize coefficient estimates
      if (l.i == 1 & w.j == 1){
        gamma.coefs <- initial.gamma.guess
        }
      
      if (l.i == 1 & w.j != 1){gamma.coefs <- gamma.coefs.save}
      
      #Find gamma coefs
      gamma.coefs <- gamma.FISTA(gamma.coefs,tDV,tDVD,y,tau,lambda,tolerance)
      
      #Record info on each combination
      rough.function <- D %*% gamma.coefs
      smoother <- t(H)%*%t(left.half)
      basis.coefs <- smoother%*%(y-rough.function)
      SSR <- sum((y-rough.function-basis%*%basis.coefs)^2)
      RMSE[l.i,w.j]  <- sqrt(SSR/n)
      sparsity.0norm[l.i,w.j] <- sum(gamma.coefs != 0)
      hat.matrix <-  basis %*% smoother
      edf[l.i,w.j]  <- min(n , sum(diag( hat.matrix )) + sparsity.0norm[l.i,w.j] )
      s.est <- SSR/(n)
      AIC[l.i,w.j] <- log(s.est) + (n+2*edf[l.i,w.j])/n
      AICc[l.i,w.j] <- log(s.est) + (n+edf[l.i,w.j])/(n-edf[l.i,w.j])
      BIC[l.i,w.j] <- log(s.est) + (log(n)*edf[l.i,w.j])/n
      
      if(l.i == 1){gamma.coefs.save <- gamma.coefs}
    }
  }
  
  
  #Finding elbow point on plot of EDF, sqrt(SSE), and sparsity
  #There are clear upper and lower bounds on each of these quantities.
  #Create a plane through upper boundary points and compute distances to plane from each point.
  
  X <- cbind(1,x)
  ls <- lm(y~X-1)
  max.RMSE <- sqrt(sum((y-ls$fitted.values)^2)/n)
  elbow.distance <- c(1 - RMSE/max.RMSE - sparsity.0norm/p - edf/(n))
  
  lambda.vec <- rep(lambda.grid,length(omega.grid))
  omega.vec <- c(matrix(rep(omega.grid,length(lambda.grid)),byrow=TRUE,ncol=length(omega.grid)))
  edf.vec <- c(edf)
  elbow.info <- data.frame(lambda=lambda.vec, omega=omega.vec, complexity.of.gamma = c(sparsity.0norm)/p, edf = c(edf)/n, RMSE = c(RMSE)/max.RMSE, elbow.distance=elbow.distance)
  
  
  if(method=="elbow"){best.index <- which.max(elbow.distance)}
  if(method=="AIC"){best.index <- which.min(AIC)}
  if(method=="AICc"){best.index <- which.min(AICc)}
  if(method=="BIC"){best.index <- which.min(BIC)}
  
  
  best.lambda <- lambda.vec[best.index]
  best.omega <- omega.vec[best.index]
  best.edf <- edf[best.index]
  
  #Find model on best combination
  lambda <- best.lambda
  w <- best.omega
  gamma.coefs <- rep(0,p)
  
  S <- solve(tbb+w*R)
  V <- (diag(rep(1,n)) - basis%*%S%*%t(basis))
  tDV <- t(D)%*%V
  tDVD <- tDV%*%D
  
  #Step.size
  tau <- 1/eigs_sym(tDVD,1)$values
  
  
  #Find gamma coefs
  gamma.coefs <- gamma.FISTA(gamma.coefs,tDV,tDVD,y,tau,lambda,tolerance)
  
  #Results from best lambda and omega
  best.gamma.coefs <- gamma.coefs
  rough.function <- D %*% gamma.coefs
  smoother <- S%*%t(basis)
  basis.coefs <- smoother%*%(y-rough.function)
  smooth.trend <- basis%*%basis.coefs
  fitted.values <-  smooth.trend + rough.function
  residuals <- y - fitted.values
  
  
  #Identify change points
  change.points <- (2:n)[gamma.coefs != 0]
  change.points.trimmed <- (2:n)[abs(gamma.coefs) > at.least]
  
  time <- (proc.time()-init.time)[3]
  
  #Write Results
  return(list(
              fitted.values = fitted.values,
              basis.functions = basis,
              basis.coefs = basis.coefs,
              gamma.coefs = gamma.coefs,
              smooth.trend = smooth.trend,
              rough.function = rough.function,
              rough.function.basis = D,
              residuals = residuals,
              lambda.grid = lambda.grid,
              omega.grid = omega.grid,
              elbow.info = elbow.info,
              AIC.info = AIC,
              AICc.info = AICc,
              BIC.info = BIC,
              best.lambda = best.lambda,
              best.omega = best.omega,
              best.edf = best.edf,
              parameter.info = data.frame(lambda = lambda.vec,
                                          omega = omega.vec,
                                          edf = edf.vec,
                                          elbow.depth = elbow.distance),
              best.index = best.index,
              change.points = change.points,
              change.points.trimmed = change.points.trimmed,
              computation.time = time))
}


