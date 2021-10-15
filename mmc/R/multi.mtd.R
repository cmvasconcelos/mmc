require(matrixcalc)

#Calculate Transition Matrices
ProbMatrixF <- function(s) {
  n = nrow(s)
  m1 = max(s)
  m0 = min(s)
  f = array(1, dim = c(m1, m1, (ncol(s) * ncol(s))))
  d = 0
  
  for (i in 1:ncol(s)) {
    for (j in 1:ncol(s)) {
      d = d + 1
      for (k1 in m0:m1) {
        for (k2 in m0:m1) {
          c = 0
          for (t in 2:n) {
            if (s[t - 1, j] == k1 && s[t, i] == k2) {
              c = c + 1 #count number of sequences of type k1 -> k2
            }
          }
          f[k2, k1, d] = c
        }
      }
    }
  }
  return(f)
}
#Calculate Probability Transition Matrices
ProbMatrixQ <- function(f, s) {
  n = nrow(s)
  m1 = max(s)
  m0 = min(s)
  q = array(0, dim = c(m1, m1, (ncol(s) * ncol(s))))
  d = 0
  
  for (i in 1:ncol(s)) {
    for (j in 1:ncol(s)) {
      d = d + 1
      for (k1 in m0:m1) {
        for (k2 in m0:m1) {
          if (sum(f[, k1, d]) > 0) {
            q[k2, k1, d] = f[k2, k1, d] / sum(f[, k1, d])
          }
        }
      }
    }
  }
  
  return(q)
}

#Array of number of sequences: n_i0_il
NumberOfSequences <- function(f, s) {
  m1 <- max(s)
  
  n = m1 * m1
  
  n_i0_il <- array(0, dim = c(n, (ncol(s) * ncol(s))))
  
  for (g in 1:(ncol(s) * ncol(s))) {
    n_i0_il[, g] = vec(f[, , g])
  }
  return(n_i0_il)
}
#Array of the probabilities values: q_i0_il
ArrayQ <- function(q, s) {
  m1 <- max(s)
  n = m1 * m1
  q_i0_il <- array(0, dim = c(n, (ncol(s) * ncol(s))))
  
  for (g in 1:(ncol(s) * ncol(s))) {
    q_i0_il[, g] = vec(q[, , g])
  }
  return(q_i0_il)
}

#Initial Values p. 8 de Berchtold 2001
CalculateInitialValues <- function(f, s) {
  #Calculate initial values for all equations
  n = nrow(s)
  m1 = max(s)
  u <- array(data = 0, dim = c(1, (ncol(s) * ncol(s))))
  sumu <- array(data = 0, dim = c(1, ncol(s)))
  lambda <- u
  
  for (g in 1:(ncol(s) * (ncol(s)))) {
    cg <- f[, , g]
    tc <- sum(cg)
    sr <- rowSums(cg)
    sc <- colSums(cg)
    
    num <- 0
    for (i in 1:m1) {
      for (j in 1:m1) {
        if (cg[i, j] != 0) {
          num <-
            num + cg[i, j] * (log2(sc[i]) + log2(sr[j]) - log2(cg[i, j]) - log2(tc))
        }
      }
    }
    
    den <- 0
    for (j in 1:m1) {
      if (sc[j] != 0 & tc != 0) {
        den <- den + sc[j] * (log2(sc[j]) - log2(tc))
      }
    }
    
    if (den != 0) {
      u[g] = num / den
    } else{
      u[g] = 0
    }
  }
  
  j <- 0
  for (i in 1:ncol(s)) {
    sumu[i] <- sum(u[(i + j):(ncol(s) * i)])
    j = j + ncol(s) - 1
  }
  
  sumu <- rep(sumu, each = ncol(s))
  
  for (j in 1:(ncol(s) * ncol(s))) {
    lambda[j] <- u[j] / sumu[j]
  }
  return(lambda = lambda)
}

#Calculate first partial derivatives for numerical maximization
PartialDerivatives <- function(n_i0_il, q_i0_il, lambda) {
  l <- length(lambda)
  pd_lambda <- rep(0, l)
  
  for (j in 1:l) {
    for (k in 1:length(n_i0_il[, j])) {
      if ((q_i0_il[k, ] %*% lambda) != 0) {
        pd_lambda[j] <-
          pd_lambda[j] + n_i0_il[k, j] * (q_i0_il[k, j] / (q_i0_il[k, ] %*% lambda))
      }
    }
  }
  return(pd_lambda)
}

#Log-likelihood
LogLikelihood <- function(lambda, n_i0_il, q_i0_il) {
  ll <- 0
  
  for (j in 1:ncol(q_i0_il)) {
    for (i in 1:length(n_i0_il[, j])) {
      if (q_i0_il[i, ] %*% lambda > 0) {
        ll <- ll + n_i0_il[i, j] * log(q_i0_il[i, ] %*% lambda)
      }
    }
  }
  ll
}

#Numerical maximization: Berchtold (2001)
OptimizeLambda <-
  function(lambda,
           delta,
           n_i0_il,
           q_i0_il,
           is_constrained,
           delta_stop) {
    delta_it <- delta
    
    ll <- 0
    ll <- LogLikelihood(lambda, n_i0_il, q_i0_il)
    pd_lambda <- 0
    pd_lambda <- PartialDerivatives(n_i0_il, q_i0_il, lambda)
    
    i_inc <- which.max(pd_lambda)
    i_dec <- which.min(pd_lambda)
    lambda_r <- lambda
    par_inc <- lambda_r[i_inc]
    par_dec <- lambda_r[i_dec]
    
    if (is_constrained) {
      if (par_inc == 1) {
        return(list(
          lambda = lambda,
          ll = ll,
          delta = delta
        ))
      }
      if (par_inc + delta_it > 1) {
        delta_it <- 1 - par_inc
      }
      if (par_dec == 0) {
        pd_lambda_sorted <- sort(pd_lambda, index.return = TRUE)
        i_dec <-
          pd_lambda_sorted$ix[min(which(lambda[pd_lambda_sorted$ix] > 0))]
        par_dec <- lambda_r[i_dec]
      }
      if (par_dec - delta_it < 0) {
        delta_it <- 1 - par_dec
      }
    }
    
    #Otimization loop
    while (TRUE) {
      if (is_constrained) {
        delta_it <- min(c(delta_it, 1 - par_inc, par_dec))
      }
      
      new_lambda <- lambda_r
      new_lambda[i_inc] <- par_inc + delta_it
      new_lambda[i_dec] <- par_dec - delta_it
      new_ll <- LogLikelihood(new_lambda, n_i0_il, q_i0_il)
      if (new_ll > ll) {
        return(list(
          lambda = new_lambda,
          ll = new_ll,
          delta = delta
        ))
      } else {
        if (delta_it <= delta_stop) {
          delta <- 2 * delta
          return(list(
            lambda = new_lambda,
            ll = new_ll,
            delta = delta
          ))
        }
        delta_it <- delta_it / 2
      }
    }
    
  }

#Second partial derivatives computation
Inference <- function(n_i0_il, q_i0_il, lambda) {
  l <- length(lambda)
  pd <- rep(0, l)
  hess <- matrix(0, l, l)
  var <- rep(0, l)
  se <- rep(0, l)
  zstat <- rep(0, l)
  pvalue <- rep(0, l)
  n <- length(n_i0_il)
  
  
  for (i in 1:l) {
    pd <- rep(0, l)
    for (j in 1:l) {
      for (k in 1:length(n_i0_il[, j])) {
        if ((q_i0_il[k, ] %*% lambda) != 0) {
          pd[j] <-
            pd[j] + n_i0_il[k, j] * (-(q_i0_il[k, i] * q_i0_il[k, j]) / (q_i0_il[k, ] %*% lambda) ^
                                       2)
        }
      }
      hess[i, j] <- pd[j]
    }
  }
  
  
  hessinv <- solve(-hess)
  
  var <- diag(hessinv)
  
  for (j in 1:l) {
    se[j] <- sqrt(var[j])
    zstat[j] <- lambda[j] / se[j]
    pvalue[j] <- 2 * (1 - pnorm(abs(zstat[j])))
  }
  
  return(l = list(
    se = se,
    zstat = zstat,
    pvalue = pvalue
  ))
}

#Output table
output.table <- function(estimates, se, zstat, pvalue) {
  stars <- rep("", length(pvalue))
  
  if (!is.character(se))
  {
    stars[pvalue <= 0.01] <- "***"
    stars[pvalue > 0.01 & pvalue <= 0.05] <- "**"
    stars[pvalue > 0.05 & pvalue <= 0.1] <- "*"
    
    se <- formatC(se, digits = 6, format = "f")
    zstat <- formatC(zstat, digits = 3, format = "f")
    pvalue <- formatC(pvalue, digits = 3, format = "f")
    stars <- format(stars, justify = "left")
  }
  else
  {
    se <- rep(".", length(pvalue))
    zstat <- rep(".", length(pvalue))
    pvalue <- rep(".", length(pvalue))
  }
  
  estimates <- formatC(estimates, digits = 6, format = "f")
  results <-
    data.frame(cbind(estimates, se, zstat, pvalue, stars), row.names = NULL)
  
  namcol <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)", "")
  colnames(results) <- namcol
  
  return(results)
  
}

#MMC MTD estimation
multi.mtd <-
  function(y,
           deltaStop = 0.0001,
           is_constrained = TRUE,
           delta = 0.1) {
    results <- list(NA)
    
    nr.eq <- ncol(y)
    m1 <- max(y)
    
    f <- ProbMatrixF(y)
    q <- ProbMatrixQ(f, y)
    
    n_i0_il <- NumberOfSequences(f, y)
    q_i0_il <- ArrayQ(q, y)
    
    lambda <- CalculateInitialValues(f, y)
    
    j = 0
    k <- 0
    for (i in 1:nr.eq) {
      ll <- 0
      
      n <- n_i0_il[, (i + j):(nr.eq * i)]
      q <- q_i0_il[, (i + j):(nr.eq * i)]
      l <- lambda[(i + j):(nr.eq * i)]
      
      opt <- OptimizeLambda(l, delta, n, q, is_constrained, deltaStop)
      
      lambdaOptim <- opt$lambda
      ll <- opt$ll
      
      inf <- Inference(n, q, lambdaOptim)
      
      results[[i + k]] <-
        output.table(lambdaOptim, inf$se, inf$zstat, inf$pvalue)
      results[[(i + k + 1)]] <- ll
      names(results)[[i + k]] <- paste("Equation", i, sep = " ")
      names(results)[[(i + k + 1)]] <- paste("LogLik", i, sep = " ")
      
      j = j + nr.eq - 1
      k = k + 1
    }
    
    return(results)
    
  }
