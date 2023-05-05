#This is not a function I am planning on including in the package. I am including this to help simulate data for function testing.
sim.data <- function(n, p.con, p.mediators, rho.con.exp, rho.mediators,
                     r2.mediator, r2.outcome, total.effect.internal,
                     n.external.ratio, is.same.external, total.effect.external){
  
  #Generate internal data
  mediator.model.regressor.cov <- matrix(rho.con.exp, nrow = 1+p.con, ncol = 1+p.con)
  diag(mediator.model.regressor.cov) <- 1
  mediator.model.regressors <- mvrnorm(n = n, mu = rep(0,1+p.con), Sigma = mediator.model.regressor.cov)
  colnames(mediator.model.regressors) <- c("A","C1.m","C2.m","C3.m","C4.m","C5.m")
  alpha.c <- matrix(0.1, nrow = p.con, ncol = p.mediators)
  alpha.a.single <- 0.6
  alpha.a <- c(rep(alpha.a.single, 10), rep(0, 40))
  mediator.cor <- matrix(rho.mediators, nrow = p.mediators, ncol = p.mediators)
  mediator.cor[1:5,1:5] <- 0.3
  mediator.cor[6:10,6:10] <- 0.3
  mediator.cor[11:15,11:15] <- 0.3
  mediator.cor[16:50,16:50] <- 0.3
  diag(mediator.cor) <- 1
  sigma.sq.mediator <- ((1-r2.mediator)/r2.mediator)*(alpha.a.single*mediator.model.regressor.cov[1,1]*alpha.a.single)
  mediator.cov <- diag(sqrt(sigma.sq.mediator), p.mediators)%*%mediator.cor%*%diag(sqrt(sigma.sq.mediator), p.mediators)
  M.internal <- matrix(mediator.model.regressors[, colnames(mediator.model.regressors) == "A"], ncol = 1)%*%matrix(alpha.a, nrow = 1) +
    mediator.model.regressors[, colnames(mediator.model.regressors) != "A"]%*%alpha.c +
    mvrnorm(n = n, rep(0,p.mediators), Sigma = mediator.cov)
  A.internal <- as.vector(mediator.model.regressors[, colnames(mediator.model.regressors) == "A"])
  C.internal <- mediator.model.regressors[, colnames(mediator.model.regressors) != "A"]
  
  beta.c <- rep(0.1, p.con)
  beta.m.single <- 0.1
  beta.m <- c(rep(beta.m.single, 5), rep(0, 5), rep(beta.m.single, 5), rep(0, 35))
  beta.a <- total.effect.internal - sum(alpha.a*beta.m)
  if(beta.a < 0){
    stop('Total effect is negative!')
  }
  outcome.cov <- matrix(0, nrow = p.mediators + 1 + p.con, ncol = p.mediators + 1 + p.con)
  outcome.cov[1:p.mediators, 1:p.mediators] <- mediator.cov + t(as.matrix(rbind(alpha.a, alpha.c)))%*%mediator.model.regressor.cov%*%as.matrix(rbind(alpha.a, alpha.c))
  outcome.cov[(1+p.mediators):(p.mediators+1+p.con), (1+p.mediators):(p.mediators+1+p.con)] <- mediator.model.regressor.cov
  outcome.cov[1:p.mediators, (1+p.mediators):(p.mediators+1+p.con)] <- t(as.matrix(rbind(alpha.a, alpha.c)))%*%mediator.model.regressor.cov
  outcome.cov[(1+p.mediators):(p.mediators+1+p.con), 1:p.mediators] <- t(outcome.cov[1:p.mediators, (1+p.mediators):(p.mediators+1+p.con)])
  #Can show this using law of total variance and reconstructing joint multivariate normal distribution from conditional and marginal
  sigma.e.sq <- ((1-r2.outcome)/r2.outcome)*(matrix(beta.m, nrow = 1)%*%outcome.cov[1:p.mediators, 1:p.mediators]%*%matrix(beta.m, ncol = 1))
  Y.internal <- as.vector(cbind(M.internal,mediator.model.regressors)%*%matrix(c(beta.m,beta.a,beta.c), ncol = 1)) + rnorm(n = n, mean = 0, sd = sqrt(sigma.e.sq))
  
  int.mod <- lm(Y.internal ~ A.internal + C.internal)
  T.hat.internal <- coef(int.mod)[names(coef(int.mod)) == "A.internal"]
  var.T.hat.internal <- vcov(int.mod)[row.names(vcov(int.mod)) == "A.internal", colnames(vcov(int.mod)) == "A.internal"]
  
  #Generate external data
  if(is.same.external == TRUE){
    n.external <- n.external.ratio*n
    mediator.model.regressors.ext <- mvrnorm(n = n.external, mu = rep(0,1+p.con), Sigma = mediator.model.regressor.cov)
    colnames(mediator.model.regressors.ext) <- c("A","C1.m","C2.m","C3.m","C4.m","C5.m")
    sigma.ext.sq <- sigma.e.sq + as.numeric(matrix(beta.m, nrow = 1)%*%mediator.cov%*%matrix(beta.m, ncol = 1))
    total.effect.external <- total.effect.internal
    external.con.coeffs <- (alpha.c%*%matrix(beta.m, ncol = 1))
    A.external <- mediator.model.regressors.ext[, colnames(mediator.model.regressors.ext) == "A"]
    C.external <- mediator.model.regressors.ext[, colnames(mediator.model.regressors.ext) != "A"]
    Y.external <- A.external*total.effect.external + as.vector(C.external%*%external.con.coeffs) + rnorm(n = n.external, 0, sqrt(sigma.ext.sq))
    ext.mod <- lm(Y.external ~ A.external + C.external)
    T.hat.external <- coef(ext.mod)[names(coef(ext.mod)) == "A.external"]
    var.T.hat.external <- vcov(ext.mod)[row.names(vcov(ext.mod)) == "A.external", colnames(vcov(ext.mod)) == "A.external"]
  } else {
    n.external <- n.external.ratio*n
    mediator.model.regressors.ext <- mvrnorm(n = n.external, mu = rep(0,1+p.con), Sigma = mediator.model.regressor.cov)
    colnames(mediator.model.regressors.ext) <- c("A","C1.m","C2.m","C3.m","C4.m","C5.m")
    sigma.ext.sq <- sigma.e.sq + as.numeric(matrix(beta.m, nrow = 1)%*%mediator.cov%*%matrix(beta.m, ncol = 1))
    external.con.coeffs <- (alpha.c%*%matrix(beta.m, ncol = 1))
    A.external <- mediator.model.regressors.ext[, colnames(mediator.model.regressors.ext) == "A"]
    C.external <- mediator.model.regressors.ext[, colnames(mediator.model.regressors.ext) != "A"]
    Y.external <- A.external*total.effect.external + as.vector(C.external%*%external.con.coeffs) + rnorm(n = n.external, 0, sqrt(sigma.ext.sq))
    ext.mod <- lm(Y.external ~ A.external + C.external)
    T.hat.external <- coef(ext.mod)[names(coef(ext.mod)) == "A.external"]
    var.T.hat.external <- vcov(ext.mod)[row.names(vcov(ext.mod)) == "A.external", colnames(vcov(ext.mod)) == "A.external"]
  }
  
  return(list("Y" = Y.internal, "M" = M.internal, "A" = A.internal, "C" = C.internal,
              "R2.m" = summary(lm(M.internal[,1] ~ A.internal + C.internal))$r.squared,
              "R2.o" = summary(lm(Y.internal ~ M.internal + A.internal + C.internal))$r.squared,
              "alpha.a" = alpha.a, "alpha.c" = alpha.c,
              "beta.a" = beta.a, "beta.m" = beta.m, "beta.c" = beta.c,
              "sigma.e.sq" = sigma.e.sq, "Sigma.m" = mediator.cov,
              "n.external" = n.external, "T.hat.external" = T.hat.external,
              "var.T.hat.external" = var.T.hat.external,
              "T.hat.internal" = T.hat.internal,
              "var.T.hat.internal" = var.T.hat.internal,
              "R2.ext" = summary(ext.mod)$r.squared))
}

#' Estimate unconstrained model parameters.
#'
#' @param Y A (n x 1) continuous outcome vector.
#' @param M A (n x p_m) matrix of mediators.
#' @param A A (n x 1) vector of exposures.
#' @param C A (n x p_c) matrix of confounders and adjustment covariates. If there are no confounders or adjustment covariates set C = NULL.
#' @return A list containing point estimates of the unconstrained model parameters.
unconstrained.unpenalized <- function(Y, M, A, C = NULL){
  n <- nrow(M)
  p <- ncol(M)
  if(is.null(C)){
    k <- 0
    C <- matrix(1, nrow = n)
  } else {
    k <- ncol(C)
    C <- cbind(1,C)
  }
  X <- cbind(A,C)
  tX <- t(X)
  
  alpha.hat <- t(solve(tX%*%X)%*%tX%*%M)
  alpha.a.hat <- alpha.hat[,1]
  alpha.c.hat <- alpha.hat[,2:(k+2)]
  Sigma.m.hat <- matrix(0, nrow = p, ncol = p)
  for(i in 1:n){
    Sigma.m.hat <- Sigma.m.hat + matrix(M[i,] - alpha.hat%*%X[i,], ncol = 1)%*%matrix(M[i,] - alpha.hat%*%X[i,], nrow = 1)
  }
  Sigma.m.hat <- Sigma.m.hat/n
  
  X <- cbind(X, M)
  tX <- t(X)
  beta.hat <- as.vector(solve(tX%*%X)%*%tX%*%Y)
  beta.a.hat <- beta.hat[1]
  beta.c.hat <- beta.hat[2:(k+2)]
  beta.m.hat <- beta.hat[(k+3):length(beta.hat)]
  sigma.e.sq.hat <- sum((Y - X%*%matrix(beta.hat, ncol = 1))^2)/n
  
  return(list("alpha.a.hat" = alpha.a.hat, "alpha.c.hat" = alpha.c.hat, "Sigma.m.hat" = Sigma.m.hat,
              "beta.a.hat" = beta.a.hat, "beta.c.hat" = beta.c.hat, "beta.m.hat" = beta.m.hat,
              "sigma.e.sq.hat" = sigma.e.sq.hat))
}

#' Estimate hard constraint model parameters using cyclical coordinate descent.
#'
#' @param Y A (n x 1) continuous outcome vector.
#' @param M A (n x p_m) matrix of mediators.
#' @param A A (n x 1) vector of exposures.
#' @param C A (n x p_c) matrix of confounders and adjustment covariates. If there are no confounders or adjustment covariates set C = NULL.
#' @param T.hat.external External estimate of the total effect.
#' @param err.tol.out Termination condition for cyclical coordinate descent algorithm with respect to the outcome model parameters.
#' @param err.tol.med Termination condition for cyclical coordinate descent algorithm with respect to the mediator model parameters.
#' @param max.itr Maximum number of iterations for cyclical coordinate descent algorithm.
#' @return A list containing point estimates of the hard constraint model parameters and an indicator of whether the algorithm converges.
constrained.unpenalized <- function(Y, M, A, C = NULL, T.hat.external, err.tol.out = 1e-08, err.tol.med = 1e-08, max.itr = 10000){
  n <- nrow(M)
  p <- ncol(M)
  if(is.null(C)){
    k <- 0
    C <- matrix(1, nrow = n)
  } else {
    k <- ncol(C)
    C <- cbind(1,C)
  }
  X <- cbind(A,C)
  tX <- t(X)
  
  tC <- t(C)
  CtC.inv <- solve(tC%*%C)
  CtC.inv.tC <- CtC.inv%*%tC
  
  sum.sq.exp <- sum(A^2)
  
  #Calculate Initial Values Using Unconstrained Method
  alpha.hat.init <- t(solve(tX%*%X)%*%tX%*%M)
  alpha.a.hat.init <- alpha.hat.init[,1]
  alpha.c.hat.init <- alpha.hat.init[,2:(k+2)]
  Sigma.m.hat.init <- matrix(0, nrow = p, ncol = p)
  for(i in 1:n){
    Sigma.m.hat.init <- Sigma.m.hat.init + matrix(M[i,] - alpha.hat.init%*%X[i,], ncol = 1)%*%matrix(M[i,] - alpha.hat.init%*%X[i,], nrow = 1)
  }
  Sigma.m.hat.init <- Sigma.m.hat.init/n
  
  X <- cbind(X, M)
  tX <- t(X)
  beta.hat.init <- as.vector(solve(tX%*%X)%*%tX%*%Y)
  beta.a.hat.init <- beta.hat.init[1]
  beta.c.hat.init <- beta.hat.init[2:(k+2)]
  beta.m.hat.init <- beta.hat.init[(k+3):length(beta.hat.init)]
  sigma.e.sq.hat.init <- sum((Y - X%*%matrix(beta.hat.init, ncol = 1))^2)/n
  
  #Constrained coordinate descent algorithm
  beta.m.prev <- beta.m.hat.init
  beta.c.prev <- beta.c.hat.init
  alpha.c.prev <- alpha.c.hat.init
  alpha.a.prev <- alpha.a.hat.init
  sigma.e.sq.prev <- sigma.e.sq.hat.init
  Sigma.m.prev <- Sigma.m.hat.init
  beta.m.new <- rep(NA,p)
  beta.c.new <- rep(NA,k+1)
  alpha.c.new <- matrix(NA, nrow = p, ncol = k+1)
  alpha.a.new <- rep(NA,p)
  sigma.e.sq.new <- NA
  Sigma.m.new <- matrix(NA, nrow = p, ncol = p)
  delta.alpha <- 1
  delta.Sigma.m <- 1
  delta.beta <- 1
  delta.sigma.e.sq <- 1
  cnt <- 1
  while (((((delta.alpha > err.tol.med) | (delta.Sigma.m > err.tol.med)) | (delta.beta > err.tol.out)) | (delta.sigma.e.sq > err.tol.out)) & (cnt <= max.itr)) {
    alpha.c.new <- t(CtC.inv.tC%*%(M - matrix(A, ncol = 1)%*%matrix(alpha.a.prev, nrow = 1)))
    alpha.tmp.vec <- rep(0, p)
    for(i in 1:n){
      alpha.tmp.vec <- alpha.tmp.vec + A[i]*(solve(Sigma.m.prev)%*%matrix(M[i,] - as.vector(alpha.c.new%*%matrix(C[i,], ncol = 1)), ncol = 1) - (1/sigma.e.sq.prev)*(Y[i] - T.hat.external*A[i] - sum(M[i,]*beta.m.prev) - sum(C[i,]*beta.c.prev))*matrix(beta.m.prev, ncol = 1))
    }
    alpha.a.new <- as.vector((1/sum.sq.exp)*solve(solve(Sigma.m.prev) + (1/sigma.e.sq.prev)*matrix(beta.m.prev, ncol = 1)%*%matrix(beta.m.prev, nrow = 1))%*%matrix(alpha.tmp.vec, ncol = 1))
    Sigma.m.tmp <- matrix(0, nrow = p, ncol = p)
    for(i in 1:n){
      Sigma.m.tmp <- Sigma.m.tmp + matrix(M[i,] - A[i]*alpha.a.new - as.vector(alpha.c.new%*%matrix(C[i,], ncol = 1)), ncol = 1)%*%matrix(M[i,] - A[i]*alpha.a.new - as.vector(alpha.c.new%*%matrix(C[i,], ncol = 1)), nrow = 1)
    }
    Sigma.m.new <- Sigma.m.tmp/n
    
    beta.c.new <- as.vector(CtC.inv.tC%*%(Y - (T.hat.external - sum(alpha.a.new*beta.m.prev))*A - M%*%matrix(beta.m.prev, ncol = 1)))
    beta.tmp.vec <- rep(0, p)
    resid.med <- matrix(0, nrow = p, ncol = p)
    cor.mod <- rep(0, p)
    for(i in 1:n){
      resid.med <- resid.med + matrix(M[i,] - A[i]*alpha.a.new, ncol = 1)%*%matrix(M[i,] - A[i]*alpha.a.new, nrow = 1)
      cor.mod <- cor.mod + (Y[i] - A[i]*T.hat.external - sum(C[i,]*beta.c.new))*matrix(M[i,] - A[i]*alpha.a.new, ncol = 1)
    }
    beta.m.new <- as.vector(solve(resid.med)%*%cor.mod)
    sigma.e.sq.new <- sum((Y - (T.hat.external - sum(alpha.a.new*beta.m.new))*A - M%*%matrix(beta.m.new, ncol = 1) - C%*%matrix(beta.c.new, ncol = 1))^2)/n
    
    #Check convergence (average error tolerance per regression coefficient)
    delta.beta <- (sum((beta.m.new - beta.m.prev)^2) + sum((beta.c.new - beta.c.prev)^2))/(p+k+1)
    delta.alpha <- (sum((alpha.a.new - alpha.a.prev)^2) + sum((alpha.c.new - alpha.c.prev)^2))/(p*(k+1)+p)
    delta.sigma.e.sq <- abs(sigma.e.sq.new - sigma.e.sq.prev)
    delta.Sigma.m <- sum(abs(Sigma.m.new - Sigma.m.prev))/p^2
    
    #Update
    beta.m.prev <- beta.m.new
    beta.c.prev <- beta.c.new
    sigma.e.sq.prev <- sigma.e.sq.new
    alpha.c.prev <- alpha.c.new
    alpha.a.prev <- alpha.a.new
    Sigma.m.prev <- Sigma.m.new
    
    cnt <- cnt + 1
  }
  
  if(cnt > max.itr){
    converged <- 0
  } else {
    converged <- 1
  }
  
  return(list("alpha.a.hat" = alpha.a.new, "alpha.c.hat" = alpha.c.new, "Sigma.m.hat" = Sigma.m.new,
              "beta.m.hat" = beta.m.new, "beta.c.hat" = beta.c.new, "sigma.e.sq.hat" = sigma.e.sq.new,
              "converged" = converged))
}

#' Cyclical coordinate descent algorithm for the M-step in the EM Algorithm for the maximizing the soft constraint model likelihood.
#'
#' @param Y A (n x 1) continuous outcome vector.
#' @param M A (n x p_m) matrix of mediators.
#' @param A A (n x 1) vector of exposures.
#' @param C A (n x p_c) matrix of confounders and adjustment covariates. If there are no confounders or adjustment covariates set C = NULL.
#' @param first.moment Posterior expectation of the total effect parameter.
#' @param second.moment Posterior expection of the squared total effect parameter.
#' @param err.tol.out Termination condition for cyclical coordinate descent algorithm with respect to the outcome model parameters.
#' @param err.tol.med Termination condition for cyclical coordinate descent algorithm with respect to the mediator model parameters.
#' @param max.itr Maximum number of iterations for cyclical coordinate descent algorithm.
#' @return A list containing point estimates of the soft constraint model parameters and an indicator of whether the algorithm converges.
rand.eff.coord.desc.unpenalized <- function(Y, M, A, C = NULL, first.moment, second.moment, err.tol.out = 1e-08, err.tol.med = 1e-08, max.itr = 10000){
  init.vals <- unconstrained.unpenalized(Y = Y, M = M, A = A, C = C)
  alpha.a.prev <- init.vals$alpha.a.hat
  alpha.c.prev <- init.vals$alpha.c.hat
  Sigma.m.prev <- init.vals$Sigma.m.hat
  Sigma.m.inv.prev <- solve(Sigma.m.prev)
  beta.c.prev <- init.vals$beta.c.hat
  beta.m.prev <- init.vals$beta.m.hat
  sigma.e.sq.prev <- init.vals$sigma.e.sq.hat
  
  n <- nrow(M)
  p <- ncol(M)
  if(is.null(C)){
    k <- 0
    C <- matrix(1, nrow = n)
  } else {
    k <- ncol(C)
    C <- cbind(1,C)
  }
  X <- cbind(A,C)
  tX <- t(X)
  
  AtA <- sum(A^2)
  CtC.inv.tC <- solve(t(C)%*%C)%*%t(C)
  
  beta.m.new <- rep(NA,p)
  beta.c.new <- rep(NA,k+1)
  alpha.c.new <- matrix(NA, nrow = p, ncol = k+1)
  alpha.a.new <- rep(NA,p)
  sigma.e.sq.new <- NA
  Sigma.m.new <- matrix(NA, nrow = p, ncol = p)
  Sigma.m.inv.new <- matrix(NA, nrow = p, ncol = p)
  delta.alpha <- 1
  delta.Sigma.m <- 1
  delta.beta <- 1
  delta.sigma.e.sq <- 1
  cnt <- 1
  
  #Coordinate Descent Algorithm
  while (((((delta.alpha > err.tol.med) | (delta.Sigma.m > err.tol.med)) | (delta.beta > err.tol.out)) | (delta.sigma.e.sq > err.tol.out)) & (cnt <= max.itr)){
    
    #Update model parameters in mediator model
    alpha.c.new <- t(CtC.inv.tC%*%(M - matrix(A, ncol = 1)%*%matrix(alpha.a.prev, nrow = 1)))
    
    alpha.tmp.vec <- rep(0, p)
    
    for(i in 1:n){
      alpha.tmp.vec <- alpha.tmp.vec + A[i]*(Sigma.m.inv.prev%*%matrix(M[i,] - as.vector(alpha.c.new%*%matrix(C[i,], ncol = 1)), ncol = 1))
    }
    alpha.tmp.vec <- alpha.tmp.vec - (1/sigma.e.sq.prev)*sum(A*(Y - M%*%matrix(beta.m.prev, ncol = 1) - C%*%matrix(beta.c.prev, ncol = 1) - first.moment*A))*matrix(beta.m.prev, ncol = 1)
    alpha.a.new <- as.vector((1/AtA)*solve(Sigma.m.inv.prev + (1/sigma.e.sq.prev)*matrix(beta.m.prev, ncol = 1)%*%matrix(beta.m.prev, nrow = 1))%*%matrix(alpha.tmp.vec, ncol = 1))
    
    Sigma.m.tmp <- matrix(0, nrow = p, ncol = p)
    for(i in 1:n){
      Sigma.m.tmp <- Sigma.m.tmp + matrix(M[i,] - A[i]*alpha.a.new - as.vector(alpha.c.new%*%matrix(C[i,], ncol = 1)), ncol = 1)%*%matrix(M[i,] - A[i]*alpha.a.new - as.vector(alpha.c.new%*%matrix(C[i,], ncol = 1)), nrow = 1)
    }
    Sigma.m.new <- Sigma.m.tmp/n
    Sigma.m.inv.new <- solve(Sigma.m.new)
    
    #Update model parameters in outcome model
    beta.c.new <- CtC.inv.tC%*%(Y - M%*%matrix(beta.m.prev, ncol = 1) - (first.moment - sum(alpha.a.new*beta.m.prev))*A)
    beta.m.new <- solve(t(M - A%*%matrix(alpha.a.new, nrow = 1))%*%(M - A%*%matrix(alpha.a.new, nrow = 1)))%*%t(M - A%*%matrix(alpha.a.new, nrow = 1))%*%(Y - C%*%matrix(beta.c.new, ncol = 1) - first.moment*A)
    sigma.e.sq.tmp.vec <- Y - C%*%matrix(beta.c.new, ncol = 1) - M%*%matrix(beta.m.new, ncol = 1) + sum(alpha.a.new*beta.m.new)*A
    sigma.e.sq.new <- (sum(sigma.e.sq.tmp.vec^2) - 2*first.moment*sum(A*sigma.e.sq.tmp.vec) + AtA*second.moment)/n
    
    #Check convergence (average error tolerance)
    delta.beta <- (sum((beta.m.new - beta.m.prev)^2) + sum((beta.c.new - beta.c.prev)^2))/(p+k+1)
    delta.alpha <- (sum((alpha.a.new - alpha.a.prev)^2) + sum((alpha.c.new - alpha.c.prev)^2))/(p*(k+1)+p)
    delta.sigma.e.sq <- abs(sigma.e.sq.new - sigma.e.sq.prev)
    delta.Sigma.m <- sum(abs(Sigma.m.new - Sigma.m.prev))/p^2
    
    #Update
    beta.m.prev <- beta.m.new
    beta.c.prev <- beta.c.new
    sigma.e.sq.prev <- sigma.e.sq.new
    alpha.c.prev <- alpha.c.new
    alpha.a.prev <- alpha.a.new
    Sigma.m.prev <- Sigma.m.new
    Sigma.m.inv.prev <- Sigma.m.inv.new
    
    cnt <- cnt + 1
  }
  
  if(cnt > max.itr){
    converged <- 0
  } else {
    converged <- 1
  }
  
  return(list("alpha.a.hat" = alpha.a.new, "alpha.c.hat" = alpha.c.new, "Sigma.m.hat" = Sigma.m.new,
              "beta.m.hat" = beta.m.new, "beta.c.hat" = beta.c.new, "sigma.e.sq.hat" = sigma.e.sq.new,
              "converged" = converged))
}

#' Estimate soft constraint model parameters using the EM algorithm.
#'
#' @param Y A (n x 1) continuous outcome vector.
#' @param M A (n x p_m) matrix of mediators.
#' @param A A (n x 1) vector of exposures.
#' @param C A (n x p_c) matrix of confounders and adjustment covariates. If there are no confounders or adjustment covariates set C = NULL.
#' @param rand.eff.mean Mean of the random effects distribution for the internal total effect parameter.
#' @param rand.eff.var Variance of the random effects distribution for the internal total effect parameter.
#' @param T.hat.external External estimate of the total effect.
#' @param var.T.hat.external Estimated variance of the external total effect estimator.
#' @param err.tol.out Termination condition for cyclical coordinate descent algorithm with respect to the outcome model parameters.
#' @param err.tol.med Termination condition for cyclical coordinate descent algorithm with respect to the mediator model parameters.
#' @param max.itr Maximum number of iterations for cyclical coordinate descent algorithm.
#' @return A list containing point estimates of the soft constraint model parameters and an indicator of whether the algorithm converges.
rand.eff.unpenalized <- function(Y, M, A, C = NULL, rand.eff.mean, rand.eff.var, T.hat.external = T.hat.external, var.T.hat.external = var.T.hat.external, err.tol.out = 1e-08, err.tol.med = 1e-08, max.itr = 10000){
  init.vals <- unconstrained.unpenalized(Y = Y, M = M, A = A, C = C)
  alpha.a.prev <- init.vals$alpha.a.hat
  alpha.c.prev <- init.vals$alpha.c.hat
  Sigma.m.prev <- init.vals$Sigma.m.hat
  Sigma.m.inv.prev <- solve(Sigma.m.prev)
  beta.c.prev <- init.vals$beta.c.hat
  beta.m.prev <- init.vals$beta.m.hat
  sigma.e.sq.prev <- init.vals$sigma.e.sq.hat
  
  n <- nrow(M)
  p <- ncol(M)
  if(is.null(C)){
    k <- 0
    C <- matrix(1, nrow = n)
  } else {
    k <- ncol(C)
    C <- cbind(1,C)
  }
  
  AtA <- sum(A^2)
  
  beta.m.new <- rep(NA,p)
  beta.c.new <- rep(NA,k+1)
  alpha.c.new <- matrix(NA, nrow = p, ncol = k+1)
  alpha.a.new <- rep(NA,p)
  sigma.e.sq.new <- NA
  Sigma.m.new <- matrix(NA, nrow = p, ncol = p)
  Sigma.m.inv.new <- matrix(NA, nrow = p, ncol = p)
  delta.alpha <- 1
  delta.Sigma.m <- 1
  delta.beta <- 1
  delta.sigma.e.sq <- 1
  cnt <- 1
  
  #EM Algorithm
  while (((((delta.alpha > err.tol.med) | (delta.Sigma.m > err.tol.med)) | (delta.beta > err.tol.out)) | (delta.sigma.e.sq > err.tol.out)) & (cnt <= max.itr)){
    #E-Step
    first.moment <- (AtA/sigma.e.sq.prev + 1/rand.eff.var)^(-1)*((sum(A*(Y - C%*%matrix(beta.c.prev, ncol = 1) - M%*%matrix(beta.m.prev, ncol = 1) + A*sum(alpha.a.prev*beta.m.prev)))/sigma.e.sq.prev) + (rand.eff.mean/rand.eff.var))
    second.moment <- (AtA/sigma.e.sq.prev + 1/rand.eff.var)^(-1) + (first.moment)^2
    
    #M-Step (Coordinate Descent)
    if(k == 0){
      cd.out <- rand.eff.coord.desc.unpenalized(Y = Y, M = M, A = A, C = NULL, first.moment = first.moment, second.moment = second.moment)
    } else if(k == 1) {
      cd.out <- rand.eff.coord.desc.unpenalized(Y = Y, M = M, A = A, C = matrix(C[,-1], ncol = 1), first.moment = first.moment, second.moment = second.moment)
    } else if(k > 1){
      cd.out <- rand.eff.coord.desc.unpenalized(Y = Y, M = M, A = A, C = C[,-1], first.moment = first.moment, second.moment = second.moment)
    }
    alpha.c.new <- cd.out$alpha.c.hat
    alpha.a.new <- cd.out$alpha.a.hat
    Sigma.m.new <- cd.out$Sigma.m.hat
    Sigma.m.inv.new <- solve(Sigma.m.new)
    beta.c.new <- cd.out$beta.c.hat
    beta.m.new <- cd.out$beta.m.hat
    sigma.e.sq.new <- cd.out$sigma.e.sq.hat
    
    #Check convergence (average error tolerance)
    delta.beta <- (sum((beta.m.new - beta.m.prev)^2) + sum((beta.c.new - beta.c.prev)^2))/(p+k+1)
    delta.alpha <- (sum((alpha.a.new - alpha.a.prev)^2) + sum((alpha.c.new - alpha.c.prev)^2))/(p*(k+1)+p)
    delta.sigma.e.sq <- abs(sigma.e.sq.new - sigma.e.sq.prev)
    delta.Sigma.m <- sum(abs(Sigma.m.new - Sigma.m.prev))/p^2
    
    #Update
    beta.m.prev <- beta.m.new
    beta.c.prev <- beta.c.new
    sigma.e.sq.prev <- sigma.e.sq.new
    alpha.c.prev <- alpha.c.new
    alpha.a.prev <- alpha.a.new
    Sigma.m.prev <- Sigma.m.new
    
    cnt <- cnt + 1
  }
  
  first.moment <- (AtA/sigma.e.sq.new + 1/rand.eff.var)^(-1)*((sum(A*(Y - C%*%matrix(beta.c.new, ncol = 1) - M%*%matrix(beta.m.new, ncol = 1) + A*sum(alpha.a.new*beta.m.new)))/sigma.e.sq.new) + (rand.eff.mean/rand.eff.var))
  beta.a.new <- first.moment - sum(alpha.a.new*beta.m.new)
  
  if(cnt > max.itr){
    converged <- 0
  } else {
    converged <- 1
  }
  
  return(list("alpha.a.hat" = alpha.a.new, "alpha.c.hat" = alpha.c.new, "Sigma.m.hat" = Sigma.m.new,
              "beta.m.hat" = beta.m.new, "beta.c.hat" = beta.c.new, "sigma.e.sq.hat" = sigma.e.sq.new,
              "beta.a.hat" = beta.a.new, "post.var.beta.a.hat" = (AtA/sigma.e.sq.new + 1/rand.eff.var)^(-1),
              "converged" = converged))
}
