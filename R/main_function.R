# library(MASS)
# library(ggplot2)
# library(patchwork)

utils::globalVariables(c('confint.95', 'lbl', 'lcl.95', 'p.string ucl.95', 'ucl.95', 
                         'p.string'))

prog_bar_endcap = function() {
  cat('|')
}

prog_bar_progress = function(i, n) {
  if(i %% floor(0.02*n) == 0) cat('+')
}

#' Implementation of Mediation with External Summary Statistics Information (MESSI) from Boss et al. (2023).
#' 
#' @details
#' The Soft EB method should be the default method if the user is not sure which 
#' method to use.
#'
#' @param Y A (n x 1) continuous outcome vector.
#' @param M A (n x p_m) matrix of mediators.
#' @param A A (n x 1) vector of exposures.
#' @param C A (n x p_c) matrix of confounders and adjustment covariates. If there are no confounders or adjustment covariates set C = NULL.
#' @param method A string specifying which method to use. Options include 'Unconstrained', 'Hard', 'Soft EB', and 'Soft Fixed'. Default is 'Soft EB'.
#' @param T.hat.external External estimate of the total effect. Set to NULL if method = 'Unconstrained'.
#' @param var.T.hat.external Estimated variance of the external estimator of the total effect. Set to NULL if method = 'Unconstrained' or method = 'Hard'.
#' @param n.boot Number of parametric bootstrap draws for obtaining quantile-based confidence intervals for the TE and NDE. Relevant for method = 'Soft EB' and method = 'Soft Fixed'. Can set to NULL for method = 'Unconstrained' and method = 'Hard'.
#' @param s2.fixed Option to specify the tuning parameter s^2 in the soft constraint model. Only use if method = 'Soft Fixed'.
#' @return A list containing the (1) point estimates and confidence intervals for the natural direct effect, the natural indirect effect, and the total effect (2) point estimates for all mediation model parameters (3) the asymptotic variance covariance matrix corresponding to alpha_a and beta_m.
#' @examples
#' data(Med)
#' 
#' Y = Med$Y
#' M = Med$M
#' A = Med$A
#' C = Med$C
#' T.hat.external = Med$T.hat.external
#' var.T.hat.external = Med$var.T.hat.external
#' 
#' test <- messi(Y = Y, M = M, A = A, C = C, method = 'Unconstrained', T.hat.external = T.hat.external,
#'               var.T.hat.external = var.T.hat.external, s2.fixed = NULL)
#' 
#' n = Med$n
#' p = Med$p
#' 
#' plot_messi(n = n, alpha.a.hat = test$alpha.a.hat, beta.m.hat = test$beta.m.hat, 
#'            labels = paste0("M",1:p), asym.var.mat = test$asym.var.mat)
#'               
#' test <- messi(Y = Y, M = M, A = A, C = C, method = 'Hard', T.hat.external = T.hat.external,
#'               var.T.hat.external = var.T.hat.external, s2.fixed = NULL)
#' 
#' @importFrom stats coef lm pchisq pnorm quantile rchisq rnorm var vcov
#' @importFrom MASS mvrnorm
#' @importFrom ggplot2 ggplot aes geom_pointrange geom_vline ggtitle theme_classic scale_colour_identity scale_y_discrete theme element_blank element_text geom_text theme_void
#' @export
messi <- function(Y, M, A, C = NULL, method = "Soft EB", T.hat.external, var.T.hat.external, n.boot = 200, s2.fixed = NULL){
  #Conditional variance of A after regressing out C
  n <- length(Y)
  
  if(is.null(C)){
    sigma.A.sq <- var(A)
    te.internal.mod <- lm(Y ~ A)
  } else {
    sigma.A.sq <- var(lm(A ~ C)$residuals)
    te.internal.mod <- lm(Y ~ A + C)
  }
  
  T.hat.internal <- coef(te.internal.mod)[names(coef(te.internal.mod)) == "A"]
  var.T.hat.internal <- vcov(te.internal.mod)[names(coef(te.internal.mod)) == "A", names(coef(te.internal.mod)) == "A"]
  
  if(method == "Unconstrained"){
    #Obtain Point Estimates of Model Parameters
    final.fit <- unconstrained.unpenalized(Y = Y, M = M, A = A, C = C)
    
    alpha.a.hat <- final.fit$alpha.a.hat
    alpha.c.hat <- final.fit$alpha.c.hat
    beta.m.hat <- final.fit$beta.m.hat
    beta.a.hat <- final.fit$beta.a.hat
    beta.c.hat <- final.fit$beta.c.hat
    
    Sigma.m.hat <- final.fit$Sigma.m.hat
    Sigma.m.inv.hat <- solve(Sigma.m.hat)
    sigma.e.sq.hat <- final.fit$sigma.e.sq.hat
    
    #Obtain Point Estimates of Mediation Parameters
    nie.est <- sum(alpha.a.hat*beta.m.hat)
    nde.est <- beta.a.hat
    te.est <- nde.est + nie.est
    
    #Construct asymptotic confidence intervals for NDE and TE (Theorem 1 in paper)
    quad.form.beta.hat <- as.numeric(matrix(beta.m.hat, nrow = 1)%*%Sigma.m.hat%*%matrix(beta.m.hat, ncol = 1))
    quad.form.alpha.hat <- as.numeric(matrix(alpha.a.hat, nrow = 1)%*%Sigma.m.inv.hat%*%matrix(alpha.a.hat, ncol = 1))
    nde.var <- (sigma.e.sq.hat/sigma.A.sq) + sigma.e.sq.hat*quad.form.alpha.hat
    te.var <- (sigma.e.sq.hat + quad.form.beta.hat)/sigma.A.sq
    nde.ci95 <- c(nde.est - 1.96*sqrt(nde.var)/sqrt(n), nde.est + 1.96*sqrt(nde.var)/sqrt(n))
    te.ci95 <- c(te.est - 1.96*sqrt(te.var)/sqrt(n), te.est + 1.96*sqrt(te.var)/sqrt(n))
    
    #Wald test for testing if alpha_a = beta_m = 0
    asym.var.mat <- rbind(cbind(Sigma.m.hat/sigma.A.sq, matrix(0, nrow = nrow(Sigma.m.hat), ncol = ncol(Sigma.m.hat))),
                          cbind(matrix(0, nrow = nrow(Sigma.m.hat), ncol = ncol(Sigma.m.hat)), sigma.e.sq.hat*Sigma.m.inv.hat))
    wald.test.stat <- n*as.numeric(matrix(c(alpha.a.hat, beta.m.hat), nrow = 1)%*%solve(asym.var.mat)%*%matrix(c(alpha.a.hat, beta.m.hat), ncol = 1))
    p.wald.test <- pchisq(wald.test.stat, df = 2*ncol(M), ncp = 0, lower.tail = FALSE)
    
    if(p.wald.test <= 0.05){
      #If p-value is < 0.05 then construct asymptotic confidence intervals for NIE following Theorem 1
      nie.var <- (quad.form.beta.hat/sigma.A.sq) + sigma.e.sq.hat*quad.form.alpha.hat
      nie.ci95 <- c(nie.est - 1.96*sqrt(nie.var)/sqrt(n), nie.est + 1.96*sqrt(nie.var)/sqrt(n))
    } else if(p.wald.test > 0.05){
      #If p-value is > 0.05 then construct asymptotic confidence intervals for NIE following Theorem 3
      n.sim <- 10000
      simulate.ref.dist <- (1/2)*sqrt(sigma.e.sq.hat/sigma.A.sq)*(rchisq(n = n.sim, df = ncol(M)) - rchisq(n = n.sim, df = ncol(M)))
      nie.ci95 <- nie.est + quantile(simulate.ref.dist, probs = c(0.025, 0.975))/length(Y)
    }
    
    #Store Results
    med.summary <- data.frame("param" = c("NIE","NDE","TE"),
                              "est" = c(nie.est,nde.est,te.est),
                              "lcl95" = c(nie.ci95[1],nde.ci95[1],te.ci95[1]),
                              "ucl95" = c(nie.ci95[2],nde.ci95[2],te.ci95[2]))
    
  } else if(method == "Hard"){
    #Obtain Point Estimates of Model Parameters
    final.fit <- constrained.unpenalized(Y = Y, M = M, A = A, C = C, T.hat.external = T.hat.external)
    
    alpha.a.hat <- final.fit$alpha.a.hat
    alpha.c.hat <- final.fit$alpha.c.hat
    beta.m.hat <- final.fit$beta.m.hat
    beta.a.hat <- NULL
    beta.c.hat <- final.fit$beta.c.hat
    
    Sigma.m.hat <- final.fit$Sigma.m.hat
    Sigma.m.inv.hat <- solve(Sigma.m.hat)
    sigma.e.sq.hat <- final.fit$sigma.e.sq.hat
    
    #Obtain Point Estimates of Mediation Parameters
    nie.est <- sum(alpha.a.hat*beta.m.hat)
    nde.est <- T.hat.external - sum(alpha.a.hat*beta.m.hat)
    te.est <- nde.est + nie.est
    
    #Wald test for testing if alpha_a = beta_m = 0
    asym.var.mat <- rbind(cbind(solve(Sigma.m.inv.hat + matrix(beta.m.hat, ncol = 1)%*%matrix(beta.m.hat, nrow = 1)/sigma.e.sq.hat)/sigma.A.sq, matrix(0, nrow = nrow(Sigma.m.hat), ncol = ncol(Sigma.m.hat))),
                          cbind(matrix(0, nrow = nrow(Sigma.m.hat), ncol = ncol(Sigma.m.hat)), sigma.e.sq.hat*Sigma.m.inv.hat))
    wald.test.stat <- n*as.numeric(matrix(c(alpha.a.hat, beta.m.hat), nrow = 1)%*%solve(asym.var.mat)%*%matrix(c(alpha.a.hat, beta.m.hat), ncol = 1))
    p.wald.test <- pchisq(wald.test.stat, df = 2*ncol(M), ncp = 0, lower.tail = FALSE)
    
    #Construct asymptotic confidence intervals for NDE and NIE
    if(p.wald.test <= 0.05){
      #If p-value is < 0.05 then construct asymptotic confidence intervals for NDE and NIE following Theorem 2
      quad.form.beta.hat <- as.numeric(matrix(beta.m.hat, nrow = 1)%*%Sigma.m.hat%*%matrix(beta.m.hat, ncol = 1))
      quad.form.alpha.hat <- as.numeric(matrix(alpha.a.hat, nrow = 1)%*%Sigma.m.inv.hat%*%matrix(alpha.a.hat, ncol = 1))
      nie.var <- (quad.form.beta.hat/sigma.A.sq)*(sigma.e.sq.hat/(sigma.e.sq.hat+quad.form.beta.hat)) + sigma.e.sq.hat*quad.form.alpha.hat
      nde.var <- (sigma.e.sq.hat/sigma.A.sq)*(quad.form.beta.hat/(sigma.e.sq.hat+quad.form.beta.hat)) + sigma.e.sq.hat*quad.form.alpha.hat
      nie.ci95 <- c(nie.est - 1.96*sqrt(nie.var)/sqrt(n), nie.est + 1.96*sqrt(nie.var)/sqrt(n))
      nde.ci95 <- c(nde.est - 1.96*sqrt(nde.var)/sqrt(n), nde.est + 1.96*sqrt(nde.var)/sqrt(n))
    } else if(p.wald.test > 0.05){
      #If p-value is > 0.05 then construct asymptotic confidence intervals for NIE following Theorem 3
      n.sim <- 10000
      simulate.ref.dist <- (1/2)*sqrt(sigma.e.sq.hat/sigma.A.sq)*(rchisq(n = n.sim, df = ncol(M)) - rchisq(n = n.sim, df = ncol(M)))
      nie.ci95 <- nie.est + quantile(simulate.ref.dist, probs = c(0.025, 0.975))/length(Y)
      nde.ci95 <- nde.est + quantile(simulate.ref.dist, probs = c(0.025, 0.975))/length(Y)
    }
    
    #Hard constraint forces estimated TE to be equal to the external estimate
    te.ci95 <- c(te.est, te.est)
    
    #Store Results
    med.summary <- data.frame("param" = c("NIE","NDE","TE"),
                              "est" = c(nie.est,nde.est,te.est),
                              "lcl95" = c(nie.ci95[1],nde.ci95[1],te.ci95[1]),
                              "ucl95" = c(nie.ci95[2],nde.ci95[2],te.ci95[2]))
    
  } else if(method == "Soft EB"){
    #Obtain Point Estimates of Model Parameters (but also requires estimation of s^2 tuning parameter)
    s2.hat <- max(0,(T.hat.internal - T.hat.external)^2 - var.T.hat.internal)/var.T.hat.external
    if(s2.hat != 0){
      final.fit <- rand.eff.unpenalized(Y = Y, M = M, A = A, C = C, rand.eff.mean = T.hat.external, rand.eff.var = s2.hat*var.T.hat.external,
                                        T.hat.external = T.hat.external, var.T.hat.external = var.T.hat.external)
    } else if(s2.hat == 0){
      #If s^2 is estimated to be zero, then use soft constraint method for s^2 very close to zero
      s2.hat <- 0.000001
      final.fit <- rand.eff.unpenalized(Y = Y, M = M, A = A, C = C, rand.eff.mean = T.hat.external, rand.eff.var = s2.hat*var.T.hat.external,
                                        T.hat.external = T.hat.external, var.T.hat.external = var.T.hat.external)
    }
    
    alpha.a.hat <- final.fit$alpha.a.hat
    alpha.c.hat <- final.fit$alpha.c.hat
    beta.m.hat <- final.fit$beta.m.hat
    beta.a.hat <- final.fit$beta.a.hat
    beta.c.hat <- final.fit$beta.c.hat
    
    Sigma.m.hat <- final.fit$Sigma.m.hat
    Sigma.m.inv.hat <- solve(Sigma.m.hat)
    sigma.e.sq.hat <- final.fit$sigma.e.sq.hat
    
    #Obtain Point Estimates of Mediation Parameters
    nie.est <- sum(alpha.a.hat*beta.m.hat)
    nde.est <- beta.a.hat
    te.est <- nde.est + nie.est
    
    #Wald test for testing if alpha_a = beta_m = 0
    soft.const <- (1/(n*s2.hat*var.T.hat.external))*((sigma.A.sq/sigma.e.sq.hat)+(1/(n*s2.hat*var.T.hat.external)))^(-1)
    asym.var.mat <- rbind(cbind(solve(Sigma.m.inv.hat + soft.const*matrix(beta.m.hat, ncol = 1)%*%matrix(beta.m.hat, nrow = 1)/sigma.e.sq.hat)/sigma.A.sq, matrix(0, nrow = nrow(Sigma.m.hat), ncol = ncol(Sigma.m.hat))),
                          cbind(matrix(0, nrow = nrow(Sigma.m.hat), ncol = ncol(Sigma.m.hat)), sigma.e.sq.hat*Sigma.m.inv.hat))
    wald.test.stat <- n*as.numeric(matrix(c(alpha.a.hat, beta.m.hat), nrow = 1)%*%solve(asym.var.mat)%*%matrix(c(alpha.a.hat, beta.m.hat), ncol = 1))
    p.wald.test <- pchisq(wald.test.stat, df = 2*ncol(M), ncp = 0, lower.tail = FALSE)
    
    if(p.wald.test <= 0.05){
      #If p-value is < 0.05 then construct asymptotic confidence intervals for NIE following Theorem 4
      quad.form.beta.hat <- as.numeric(matrix(beta.m.hat, nrow = 1)%*%Sigma.m.hat%*%matrix(beta.m.hat, ncol = 1))
      quad.form.alpha.hat <- as.numeric(matrix(alpha.a.hat, nrow = 1)%*%Sigma.m.inv.hat%*%matrix(alpha.a.hat, ncol = 1))
      soft.const <- (1/(n*s2.hat*var.T.hat.external))*((sigma.A.sq/sigma.e.sq.hat)+(1/(n*s2.hat*var.T.hat.external)))^(-1)
      nie.var.asym <- (quad.form.beta.hat/sigma.A.sq)*(1+(soft.const*quad.form.beta.hat/sigma.e.sq.hat))^(-1) + sigma.e.sq.hat*quad.form.alpha.hat
      nie.ci95 <- c(nie.est - 1.96*sqrt(nie.var.asym/n), nie.est + 1.96*sqrt(nie.var.asym/n))
    } else if(p.wald.test > 0.05){
      #If p-value is > 0.05 then construct asymptotic confidence intervals for NIE following Theorem 3
      n.sim <- 10000
      simulate.ref.dist <- (1/2)*sqrt(sigma.e.sq.hat/sigma.A.sq)*(rchisq(n = n.sim, df = ncol(M)) - rchisq(n = n.sim, df = ncol(M)))
      nie.ci95 <- nie.est + quantile(simulate.ref.dist, probs = c(0.025, 0.975))/length(Y)
    }
    
    #Parametric Bootstrap to get Confidence Intervals for the NDE and TE
    n.boot <- n.boot
    s2.boot <- rep(NA, n.boot)
    nie.boot <- rep(NA, n.boot)
    nde.boot <- rep(NA, n.boot)
    te.boot <- rep(NA, n.boot)
    prog_bar_endcap()
    for(r in 1:n.boot){
      prog_bar_progress(r, n.boot)
      epsilon.y <- rnorm(n = n, mean = 0, sd = sqrt(sigma.e.sq.hat))
      if(is.null(C)){
        gen.M <- matrix(1, nrow = n, ncol = 1)%*%t(alpha.c.hat) + matrix(A, ncol = 1)%*%matrix(alpha.a.hat, nrow = 1) + mvrnorm(n = n, mu = rep(0, ncol(M)), Sigma = Sigma.m.hat)
        gen.Y <- A*beta.a.hat + matrix(1, nrow = n, ncol = 1)*as.vector(beta.c.hat) + as.vector(gen.M%*%beta.m.hat) + epsilon.y
        boot.mod <- lm(gen.Y ~ A)
      } else {
        gen.M <- cbind(1,C)%*%t(alpha.c.hat) + matrix(A, ncol = 1)%*%matrix(alpha.a.hat, nrow = 1) + mvrnorm(n = n, mu = rep(0, ncol(M)), Sigma = Sigma.m.hat)
        gen.Y <- A*beta.a.hat + as.vector(cbind(1,C)%*%beta.c.hat) + as.vector(gen.M%*%beta.m.hat) + epsilon.y
        boot.mod <- lm(gen.Y ~ A + C)
      }
      T.hat.internal.boot <- coef(boot.mod)[names(coef(boot.mod)) == "A"]
      var.T.hat.internal.boot <- vcov(boot.mod)[names(coef(boot.mod)) == "A", names(coef(boot.mod)) == "A"]
      s2.hat.boot <- max(0,(T.hat.internal.boot - T.hat.external)^2 - var.T.hat.internal.boot)/var.T.hat.external
      if(s2.hat.boot != 0){
        final.fit <- rand.eff.unpenalized(Y = gen.Y, M = gen.M, A = A, C = C, rand.eff.mean = T.hat.external, rand.eff.var = s2.hat.boot*var.T.hat.external,
                                          T.hat.external = T.hat.external, var.T.hat.external = var.T.hat.external)
      } else if(s2.hat.boot == 0){
        s2.hat.boot <- 0.000001
        final.fit <- rand.eff.unpenalized(Y = gen.Y, M = gen.M, A = A, C = C, rand.eff.mean = T.hat.external, rand.eff.var = s2.hat.boot*var.T.hat.external,
                                          T.hat.external = T.hat.external, var.T.hat.external = var.T.hat.external)
      }
      
      s2.boot[r] <- s2.hat.boot
      nde.boot[r] <- final.fit$beta.a.hat
      nie.boot[r] <- sum(final.fit$alpha.a.hat*final.fit$beta.m.hat)
      te.boot[r] <- nde.boot[r] + nie.boot[r]
    }
    prog_bar_endcap()
    delta.nde <- quantile(nde.boot - nde.est, probs = c(0.025,0.975))
    nde.ci95 <- c(nde.est - delta.nde[2], nde.est - delta.nde[1])
    
    delta.te <- quantile(te.boot - te.est, probs = c(0.025,0.975))
    te.ci95 <- c(te.est - delta.te[2], te.est - delta.te[1])
    
    #Store Results
    med.summary <- data.frame("param" = c("NIE","NDE","TE"),
                              "est" = c(nie.est,nde.est,te.est),
                              "lcl95" = c(nie.ci95[1],nde.ci95[1],te.ci95[1]),
                              "ucl95" = c(nie.ci95[2],nde.ci95[2],te.ci95[2]))
    
  } else if(method == "Soft Fixed"){
    #Obtain Point Estimates of Model Parameters
    final.fit <- rand.eff.unpenalized(Y = Y, M = M, A = A, C = C, rand.eff.mean = T.hat.external, rand.eff.var = s2.fixed*var.T.hat.external,
                                      T.hat.external = T.hat.external, var.T.hat.external = var.T.hat.external)
    
    alpha.a.hat <- final.fit$alpha.a.hat
    alpha.c.hat <- final.fit$alpha.c.hat
    beta.m.hat <- final.fit$beta.m.hat
    beta.a.hat <- final.fit$beta.a.hat
    beta.c.hat <- final.fit$beta.c.hat
    
    Sigma.m.hat <- final.fit$Sigma.m.hat
    Sigma.m.inv.hat <- solve(Sigma.m.hat)
    sigma.e.sq.hat <- final.fit$sigma.e.sq.hat
    
    #Obtain Point Estimates of Mediation Parameters
    nie.est <- sum(alpha.a.hat*beta.m.hat)
    nde.est <- beta.a.hat
    te.est <- nde.est + nie.est
    
    #Wald test for testing if alpha_a = beta_m = 0
    soft.const <- (1/(n*s2.fixed*var.T.hat.external))*((sigma.A.sq/sigma.e.sq.hat)+(1/(n*s2.fixed*var.T.hat.external)))^(-1)
    asym.var.mat <- rbind(cbind(solve(Sigma.m.inv.hat + soft.const*matrix(beta.m.hat, ncol = 1)%*%matrix(beta.m.hat, nrow = 1)/sigma.e.sq.hat)/sigma.A.sq,
                                matrix(0, nrow = nrow(Sigma.m.hat), ncol = ncol(Sigma.m.hat))), cbind(matrix(0, nrow = nrow(Sigma.m.hat), ncol = ncol(Sigma.m.hat)), sigma.e.sq.hat*Sigma.m.inv.hat))
    wald.test.stat <- n*as.numeric(matrix(c(alpha.a.hat, beta.m.hat), nrow = 1)%*%solve(asym.var.mat)%*%matrix(c(alpha.a.hat, beta.m.hat), ncol = 1))
    p.wald.test <- pchisq(wald.test.stat, df = 2*ncol(M), ncp = 0, lower.tail = FALSE)
    
    if(p.wald.test <= 0.05){
      #If p-value is < 0.05 then construct asymptotic confidence intervals for NIE following Theorem 4
      quad.form.beta.hat <- as.numeric(matrix(beta.m.hat, nrow = 1)%*%Sigma.m.hat%*%matrix(beta.m.hat, ncol = 1))
      quad.form.alpha.hat <- as.numeric(matrix(alpha.a.hat, nrow = 1)%*%Sigma.m.inv.hat%*%matrix(alpha.a.hat, ncol = 1))
      soft.const <- (1/(n*s2.fixed*var.T.hat.external))*((sigma.A.sq/sigma.e.sq.hat)+(1/(n*s2.fixed*var.T.hat.external)))^(-1)
      nie.var.asym <- (quad.form.beta.hat/sigma.A.sq)*(1+(soft.const*quad.form.beta.hat/sigma.e.sq.hat))^(-1) + sigma.e.sq.hat*quad.form.alpha.hat
      nie.ci95 <- c(nie.est - 1.96*sqrt(nie.var.asym/n), nie.est + 1.96*sqrt(nie.var.asym/n))
    } else if(p.wald.test > 0.05){
      #If p-value is > 0.05 then construct asymptotic confidence intervals for NIE following Theorem 3
      n.sim <- 10000
      simulate.ref.dist <- (1/2)*sqrt(sigma.e.sq.hat/sigma.A.sq)*(rchisq(n = n.sim, df = ncol(M)) - rchisq(n = n.sim, df = ncol(M)))
      nie.ci95 <- nie.est + quantile(simulate.ref.dist, probs = c(0.025, 0.975))/length(Y)
    }
    
    #Parametric Bootstrap to get Confidence Intervals for the NDE and TE
    n.boot <- n.boot
    nie.boot <- rep(NA, n.boot)
    nde.boot <- rep(NA, n.boot)
    te.boot <- rep(NA, n.boot)
    prog_bar_endcap()
    for(r in 1:n.boot){
      prog_bar_progress(r, n.boot)
      epsilon.y <- rnorm(n = n, mean = 0, sd = sqrt(sigma.e.sq.hat))
      if(is.null(C)){
        gen.M <- matrix(1, nrow = n, ncol = 1)%*%t(alpha.c.hat) + matrix(A, ncol = 1)%*%matrix(alpha.a.hat, nrow = 1) + mvrnorm(n = n, mu = rep(0, ncol(M)), Sigma = Sigma.m.hat)
        gen.Y <- A*beta.a.hat + matrix(1, nrow = n, ncol = 1)*as.vector(beta.c.hat) + as.vector(gen.M%*%beta.m.hat) + epsilon.y
      } else {
        gen.M <- cbind(1,C)%*%t(alpha.c.hat) + matrix(A, ncol = 1)%*%matrix(alpha.a.hat, nrow = 1) + mvrnorm(n = n, mu = rep(0, ncol(M)), Sigma = Sigma.m.hat)
        gen.Y <- A*beta.a.hat + as.vector(cbind(1, C)%*%beta.c.hat) + as.vector(gen.M%*%beta.m.hat) + epsilon.y
      }
      final.fit <- rand.eff.unpenalized(Y = gen.Y, M = gen.M, A = A, C = C, rand.eff.mean = T.hat.external, rand.eff.var = s2.fixed*var.T.hat.external,
                                        T.hat.external = T.hat.external, var.T.hat.external = var.T.hat.external)
      
      nde.boot[r] <- final.fit$beta.a.hat
      nie.boot[r] <- sum(final.fit$alpha.a.hat*final.fit$beta.m.hat)
      te.boot[r] <- nde.boot[r] + nie.boot[r]
    }
    prog_bar_endcap()
    delta.nde <- quantile(nde.boot - nde.est, probs = c(0.025,0.975))
    nde.ci95 <- c(nde.est - delta.nde[2], nde.est - delta.nde[1])
    
    delta.te <- quantile(te.boot - te.est, probs = c(0.025,0.975))
    te.ci95 <- c(te.est - delta.te[2], te.est - delta.te[1])
    
    #Store Results
    med.summary <- data.frame("param" = c("NIE","NDE","TE"),
                              "est" = c(nie.est,nde.est,te.est),
                              "lcl95" = c(nie.ci95[1],nde.ci95[1],te.ci95[1]),
                              "ucl95" = c(nie.ci95[2],nde.ci95[2],te.ci95[2]))
  }
  
  return(list("med.summary" = med.summary, "n" = n, "alpha.a.hat" = alpha.a.hat, "alpha.c.hat" = alpha.c.hat, "Sigma.m.hat" = Sigma.m.hat,
              "beta.m.hat" = beta.m.hat, "beta.c.hat" = beta.c.hat, "sigma.e.sq.hat" = sigma.e.sq.hat, "asym.var.mat" = asym.var.mat))
}

#' Forestplot to Summarize Estimation and Inference on alpha_a and beta_m.
#'
#' @param n Sample size of the analysis
#' @param alpha.a.hat Estimate of alpha_a, a (p_m x 1) vector.
#' @param beta.m.hat Estimate of beta_m, a (p_m x 1) vector.
#' @param labels A (p_m x 1) vector of mediator names. Make sure that the labels are in the same order as the mediators appear in the design matrix.
#' @param asym.var.mat Joint asymptotic variance-covariance matrix of alpha_a and beta_m, a (2p_m x 2p_m) matrix.
#' @examples
#' data(Med)
#' 
#' Y = Med$Y
#' M = Med$M
#' A = Med$A
#' C = Med$C
#' T.hat.external = Med$T.hat.external
#' var.T.hat.external = Med$var.T.hat.external
#' 
#' test <- messi(Y = Y, M = M, A = A, C = C, method = 'Unconstrained', T.hat.external = T.hat.external,
#'               var.T.hat.external = var.T.hat.external, s2.fixed = NULL)
#' 
#' n = Med$n
#' p = Med$p
#' 
#' plot_messi(n = n, alpha.a.hat = test$alpha.a.hat, beta.m.hat = test$beta.m.hat, 
#'            labels = paste0("M",1:p), asym.var.mat = test$asym.var.mat)
#'
#' 
#' @return Data frames and forestplots summarizing alpha_a and beta_m estimation.
#' @import patchwork
#' @export
plot_messi <- function(n, alpha.a.hat, beta.m.hat, labels, asym.var.mat){
  
  est <- c(alpha.a.hat, beta.m.hat)
  ci.95.low.bound <- est - 1.96*sqrt(diag(asym.var.mat)/n)
  ci.95.high.bound <- est + 1.96*sqrt(diag(asym.var.mat)/n)
  test.stat <- est/sqrt(diag(asym.var.mat)/n)
  pvals <- 2*pnorm(abs(test.stat), mean = 0, sd = 1, lower.tail = FALSE)
  
  plot.alpha <- data.frame("lbl" = labels, "est" = est[1:length(alpha.a.hat)],
                           "lcl.95" = ci.95.low.bound[1:length(alpha.a.hat)],
                           "ucl.95" = ci.95.high.bound[1:length(alpha.a.hat)],
                           "p" = pvals[1:length(alpha.a.hat)])
  
  plot.beta <- data.frame("lbl" = labels, "est" = est[(length(beta.m.hat)+1):(2*length(beta.m.hat))],
                          "lcl.95" = ci.95.low.bound[(length(beta.m.hat)+1):(2*length(beta.m.hat))],
                          "ucl.95" = ci.95.high.bound[(length(beta.m.hat)+1):(2*length(beta.m.hat))],
                          "p" = pvals[(length(beta.m.hat)+1):(2*length(beta.m.hat))])
  
  #Make forestplot for alpha.a
  plot.alpha$confint.95 <- paste0(format(round(plot.alpha$est, digits = 2), digits = 2), " (", format(round(plot.alpha$lcl.95, digits = 2), digits = 2), ",", format(round(plot.alpha$ucl.95, digits = 2), digits = 2), ")")
  plot.alpha$p.string <- ifelse(as.numeric(format(round(plot.alpha$p, digits = 3), digits = 3)) >= 0.001, format(round(plot.alpha$p, digits = 3), digits = 3), "<0.001")
  
  p_alpha <- ggplot(plot.alpha, aes(x = est, y = lbl, xmin = lcl.95, xmax = ucl.95)) +
    geom_pointrange(shape = 16, fill = "black") +
    geom_vline(xintercept = 0, linetype = 3) +
    ggtitle(expression(paste(hat(alpha)[a], " (95% CI)"))) +
    theme_classic() +
    scale_colour_identity() +
    scale_y_discrete(limits = rev(plot.alpha$lbl)) +
    theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.line.y = element_blank(),
          axis.ticks.y = element_blank(), axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5))
  
  data_table1 <- ggplot(data = plot.alpha, aes(y = lbl)) +
    geom_text(aes(x = 0, label = confint.95), hjust = 0.5) +
    scale_y_discrete(limits = rev(plot.alpha$lbl)) +
    ggtitle(expression(paste(hat(alpha)[a], " (95% CI)"))) +
    scale_colour_identity() +
    theme_void() + theme(plot.title = element_text(hjust = 0.5))
  
  data_table2 <- ggplot(data = plot.alpha, aes(y = lbl)) +
    geom_text(aes(x = 0, label = lbl), hjust = 0) +
    scale_y_discrete(limits = rev(plot.alpha$lbl)) +
    ggtitle("Mediators") +
    scale_colour_identity() +
    theme_void() + theme(plot.title = element_text(hjust = 0.5))
  
  data_table3 <- ggplot(data = plot.alpha, aes(y = lbl)) +
    geom_text(aes(x = 0, label = p.string), hjust = 0.5) +
    scale_y_discrete(limits = rev(plot.alpha$lbl)) +
    ggtitle("P-Value") +
    scale_colour_identity() +
    theme_void() + theme(plot.title = element_text(hjust = 0.5))
  
  p_alpha_combined <- data_table2 + p_alpha + (data_table1 + data_table3)
  
  #Make forestplot for beta.m
  plot.beta$confint.95 <- paste0(format(round(plot.beta$est, digits = 2), digits = 2), " (", format(round(plot.beta$lcl.95, digits = 2), digits = 2), ",", format(round(plot.beta$ucl.95, digits = 2), digits = 2), ")")
  plot.beta$p.string <- ifelse(as.numeric(format(round(plot.beta$p, digits = 3), digits = 3)) >= 0.001, format(round(plot.beta$p, digits = 3), digits = 3), "<0.001")
  
  p_beta <- ggplot(plot.beta, aes(x = est, y = lbl, xmin = lcl.95, xmax = ucl.95)) +
    geom_pointrange(shape = 16, fill = "black") +
    geom_vline(xintercept = 0, linetype = 3) +
    ggtitle(expression(paste(hat(beta)[m], " (95% CI)"))) +
    theme_classic() +
    scale_colour_identity() +
    scale_y_discrete(limits = rev(plot.beta$lbl)) +
    theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.line.y = element_blank(),
          axis.ticks.y = element_blank(), axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5))
  
  data_table1 <- ggplot(data = plot.beta, aes(y = lbl)) +
    geom_text(aes(x = 0, label = confint.95), hjust = 0.5) +
    scale_y_discrete(limits = rev(plot.beta$lbl)) +
    ggtitle(expression(paste(hat(beta)[m], " (95% CI)"))) +
    scale_colour_identity() +
    theme_void() + theme(plot.title = element_text(hjust = 0.5))
  
  data_table2 <- ggplot(data = plot.beta, aes(y = lbl)) +
    geom_text(aes(x = 0, label = lbl), hjust = 0) +
    scale_y_discrete(limits = rev(plot.beta$lbl)) +
    ggtitle("Mediators") +
    scale_colour_identity() +
    theme_void() + theme(plot.title = element_text(hjust = 0.5))
  
  data_table3 <- ggplot(data = plot.beta, aes(y = lbl)) +
    geom_text(aes(x = 0, label = p.string), hjust = 0.5) +
    scale_y_discrete(limits = rev(plot.beta$lbl)) +
    ggtitle("P-Value") +
    scale_colour_identity() +
    theme_void() + theme(plot.title = element_text(hjust = 0.5))
  
  p_beta_combined <- data_table2 + p_beta + (data_table1 + data_table3)
  
  #Invisible return matrix, but always plot
  
  #p_alpha_combined #Plot of alpha.a
  #p_beta_combined #Plot of beta.m
  return(list(data_table_alpha  = plot.alpha,
              data_table_beta  = plot.beta,
              plot_alpha_a = p_alpha_combined,
              plot_beta_m  = p_beta_combined))
}

# 
# #########################################################################################################
# ## Test case with null mediation effect
# #########################################################################################################
# 
# set.seed(20230419)
# n <- 500
# p <- 20
# Y <- rnorm(n = n, mean = 0, sd = 1)
# M <- mvrnorm(n = n, mu = rep(0, p), Sigma = diag(1, p))
# A <- rnorm(n = n, mean = 0, sd = 1)
# C <- NULL
# method <- "Unconstrained"
# s2.fixed <- NULL
# T.hat.external <- NULL
# var.T.hat.external <- NULL
# 
# test <- messi(Y = Y, M = M, A = A, C = C, method = method, T.hat.external = T.hat.external,
#               var.T.hat.external = var.T.hat.external, s2.fixed = s2.fixed)
# 
# plot_messi(n = n, alpha.a.hat = test$alpha.a.hat, beta.m.hat = test$beta.m.hat, labels = paste0("M",1:p), asym.var.mat = test$asym.var.mat)
# 
# set.seed(20230419)
# n <- 500
# p <- 20
# Y <- rnorm(n = n, mean = 0, sd = 1)
# M <- mvrnorm(n = n, mu = rep(0, p), Sigma = diag(1, p))
# A <- rnorm(n = n, mean = 0, sd = 1)
# C <- NULL
# method <- "Hard"
# s2.fixed <- NULL
# T.hat.external <- 0
# var.T.hat.external <- NULL
# 
# test <- messi(Y = Y, M = M, A = A, C = C, method = method, T.hat.external = T.hat.external,
#               var.T.hat.external = var.T.hat.external, s2.fixed = s2.fixed)
# 
# plot_messi(n = n, alpha.a.hat = test$alpha.a.hat, beta.m.hat = test$beta.m.hat, labels = paste0("M",1:p), asym.var.mat = test$asym.var.mat)
# 
# set.seed(20230419)
# n <- 500
# p <- 20
# Y <- rnorm(n = n, mean = 0, sd = 1)
# M <- mvrnorm(n = n, mu = rep(0, p), Sigma = diag(1, p))
# A <- rnorm(n = n, mean = 0, sd = 1)
# C <- NULL
# method <- "Soft EB"
# s2.fixed <- NULL
# T.hat.external <- 0
# var.T.hat.external <- 0.2
# 
# test <- messi(Y = Y, M = M, A = A, C = C, method = method, T.hat.external = T.hat.external,
#               var.T.hat.external = var.T.hat.external, s2.fixed = s2.fixed)
# 
# plot_messi(n = n, alpha.a.hat = test$alpha.a.hat, beta.m.hat = test$beta.m.hat, labels = paste0("M",1:p), asym.var.mat = test$asym.var.mat)
# 
# set.seed(20230419)
# n <- 500
# p <- 20
# Y <- rnorm(n = n, mean = 0, sd = 1)
# M <- mvrnorm(n = n, mu = rep(0, p), Sigma = diag(1, p))
# A <- rnorm(n = n, mean = 0, sd = 1)
# C <- NULL
# method <- "Soft Fixed"
# s2.fixed <- 1
# T.hat.external <- 0
# var.T.hat.external <- 0.2
# 
# test <- messi(Y = Y, M = M, A = A, C = C, method = method, T.hat.external = T.hat.external,
#               var.T.hat.external = var.T.hat.external, s2.fixed = s2.fixed)
# 
# plot_messi(n = n, alpha.a.hat = test$alpha.a.hat, beta.m.hat = test$beta.m.hat, labels = paste0("M",1:p), asym.var.mat = test$asym.var.mat)
# 
# #########################################################################################################
# ## Test case with non-null mediation effect + confounders
# #########################################################################################################
# 
# #Generate Data from Mediation Model
# set.seed(20220222)
# n <- 200
# p.con <- 5
# p.mediators <- 50
# rho.con.exp <- 0.2
# rho.mediators <- 0.2
# r2.mediator <- 0.05
# total.effect.internal <- 1
# r2.outcome <- 0.2
# n.external.ratio <- 100
# is.same.external <- TRUE
# total.effect.external <- total.effect.internal
# sim.dat <- sim.data(n = n, p.con = p.con, p.mediators = p.mediators, rho.con.exp = rho.con.exp,
#                     rho.mediators = rho.mediators, r2.mediator = r2.mediator, r2.outcome = r2.outcome,
#                     total.effect.internal = total.effect.internal,
#                     n.external.ratio = n.external.ratio, is.same.external = is.same.external,
#                     total.effect.external = total.effect.external)
# 
# #True parameter values that generated the data
# alpha.a.true <- sim.dat$alpha.a
# alpha.c.true <- sim.dat$alpha.c
# Sigma.m.true <- sim.dat$Sigma.m
# beta.m.true <- sim.dat$beta.m
# beta.a.true <- sim.dat$beta.a
# beta.c.true <- sim.dat$beta.c
# sigma.e.sq.true <- sim.dat$sigma.e.sq
# 
# #Estimated Internal and External Total Effects
# T.hat.internal <- sim.dat$T.hat.internal
# var.T.hat.internal <- sim.dat$var.T.hat.internal
# T.hat.external <- sim.dat$T.hat.external
# var.T.hat.external <- sim.dat$var.T.hat.external
# 
# set.seed(20230419)
# Y <- sim.dat$Y
# M <- sim.dat$M
# A <- sim.dat$A
# C <- sim.dat$C
# method <- "Unconstrained"
# s2.fixed <- NULL
# 
# test <- messi(Y = Y, M = M, A = A, C = C, method = method, T.hat.external = T.hat.external,
#               var.T.hat.external = var.T.hat.external, s2.fixed = s2.fixed)
# 
# plot_messi(n = n, alpha.a.hat = test$alpha.a.hat, beta.m.hat = test$beta.m.hat, labels = paste0("M",1:p), asym.var.mat = test$asym.var.mat)
# 
# set.seed(20230419)
# Y <- sim.dat$Y
# M <- sim.dat$M
# A <- sim.dat$A
# C <- sim.dat$C
# method <- "Hard"
# s2.fixed <- NULL
# 
# test <- messi(Y = Y, M = M, A = A, C = C, method = method, T.hat.external = T.hat.external,
#               var.T.hat.external = var.T.hat.external, s2.fixed = s2.fixed)
# 
# plot_messi(n = n, alpha.a.hat = test$alpha.a.hat, beta.m.hat = test$beta.m.hat, labels = paste0("M",1:p), asym.var.mat = test$asym.var.mat)
# 
# set.seed(20230419)
# Y <- sim.dat$Y
# M <- sim.dat$M
# A <- sim.dat$A
# C <- sim.dat$C
# method <- "Soft EB"
# s2.fixed <- NULL
# 
# test <- messi(Y = Y, M = M, A = A, C = C, method = method, T.hat.external = T.hat.external,
#               var.T.hat.external = var.T.hat.external, s2.fixed = s2.fixed)
# 
# plot_messi(n = n, alpha.a.hat = test$alpha.a.hat, beta.m.hat = test$beta.m.hat, labels = paste0("M",1:p), asym.var.mat = test$asym.var.mat)
# 
# set.seed(20230419)
# Y <- sim.dat$Y
# M <- sim.dat$M
# A <- sim.dat$A
# C <- sim.dat$C
# method <- "Soft Fixed"
# s2.fixed <- 1
# 
# test <- messi(Y = Y, M = M, A = A, C = C, method = method, T.hat.external = T.hat.external,
#               var.T.hat.external = var.T.hat.external, s2.fixed = s2.fixed)
# 
# plot_messi(n = n, alpha.a.hat = test$alpha.a.hat, beta.m.hat = test$beta.m.hat, labels = paste0("M",1:p), asym.var.mat = test$asym.var.mat)
