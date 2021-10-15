require(matrixcalc)
require(alabama)
require(fastDummies)
require(nnet)
require(Hmisc)

#Calculate conditional probabilities
ProbValuesXDependent <- function(s, x) {
  m1 <- max(s)
  t <- nrow(s)
  probs = matrix(NA, ncol = ncol(s) * ncol(s), nrow = t)

  if (m1 == 2) {
    g = 1
    for (i in 1:ncol(s)) {
      for (j in 1:ncol(s)) {
        px = list()
        px1 = list()

        s1 = ifelse(s[, i] == 1, 1, 0)
        l = 1

        sl = dummy_cols(s[, j])
        sl = sl[, -1]

        for (k in 1:m1) {
          s_l = Lag(sl[, k])
          res = glm(s1[s_l == 1] ~ x[s_l == 1],  family = binomial(link = 'logit'))
          px1[[l]] = cbind(fitted(res), 1 - fitted(res))
          px[[l]] = cbind(rep(0, length(s_l[-1])), rep(0, length(s_l[-1])))
          px[[l]][s_l[-1] == 1,] = px1[[l]]
          px[[l]] = rbind(c(NA, NA), px[[l]])
          l = l + 1
        }

        for (n in 2:(t - 1)) {
          probs[n - 1, g] = px[[s[n - 1, j]]][n, s[n, i]]
        }
        g = g + 1

      }
    }
    return(na.omit(probs))

  } else if (m1 > 2) {
    g = 1
    for (i in 1:ncol(s)) {
      for (j in 1:ncol(s)) {
        px = list()
        px1 = list()
        l = 1

        sl = dummy_cols(s[, j])
        sl = sl[, -1]

        for (k in 1:m1) {
          s_l = Lag(sl[, k])
          res = multinom(s[s_l == 1, i] ~ x[s_l == 1])
          px1[[l]] = matrix(fitted(res),
                            ncol = ncol(fitted(res)),
                            nrow = nrow(fitted(res)))
          px[[l]] = matrix(rep(0, m1 * length(s_l[-1])),
                           nrow = length(s_l[-1]),
                           ncol = m1)
          px[[l]][s_l[-1] == 1,] = px1[[l]]
          px[[l]] = rbind(rep(NA, m1), px[[l]])
          l = l + 1
        }

        for (n in 2:(t - 1)) {
          probs[n - 1, g] = px[[s[n - 1, j]]][n, s[n, i]]
        }
        g = g + 1
      }

    }
    return(na.omit(probs))
  }
}

LogLikelihood <- function(lambda, qi = q[, a]) {
  ll <- 0

  for (i in 1:nrow(qi)) {
    if (qi[i, ] %*% lambda > 0) {
      ll <- ll - log(qi[i, ] %*% lambda)
    }
  }
  ll
}


Inference <- function(hess, lambda) {
  l <- length(lambda)
  var <- rep(0, l)
  se <- rep(0, l)
  zstat <- rep(0, l)
  pvalue <- rep(0, l)

  if (is.singular.matrix(hess, tol = 1e-05) == FALSE) {
    hessinv <- solve(-hess)

    var <- diag(hessinv)

    if (any(var < 0)) {
      return(l = list(
        warning = 1,
        se = '.',
        zstat = '.',
        pvalue = '.'
      ))
    } else{
      for (j in 1:l) {
        se[j] <- sqrt(var[j])
        zstat[j] <- lambda[j] / se[j]
        pvalue[j] <- 2 * (1 - pnorm(abs(zstat[j])))
      }
      return(l = list(
        warning = 0,
        se = se,
        zstat = zstat,
        pvalue = pvalue
      ))
    }
  } else{
    return(l = list(
      warning = 1,
      se = '.',
      zstat = '.',
      pvalue = '.'
    ))
  }

}


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


mmcx <- function(y, x, initial) {
  results <- list(NA)
  nr.eq <- ncol(y)
  m1 <- max(y)

  q <- ProbValuesXDependent(s = y, x = x)

  j = 0
  k <- 0
  for (i in 1:nr.eq) {
    a = (i + j):(nr.eq * i)

    LogLikelihood <- function(lambda, qi = q[, a]) {
      ll <- 0

      for (i in 1:nrow(qi)) {
        if (qi[i, ] %*% lambda > 0) {
          ll <- ll - log(qi[i, ] %*% lambda)
        }
      }
      ll
    }

    he <- function(lambda) {
      h <- rep(NA, 1)
      h[1] <- sum(lambda) - 1
      h
    }

    hi <- function(lambda) {
      w <- rep(NA, 1)
      for (i in 1:length(lambda)) {
        w[i] <- lambda[i]
      }
      w
    }

    he.jac <- function(lambda) {
      j <- matrix(NA, 1, length(lambda))
      j[1,] <- rep(1, length(lambda))
      j
    }

    hi.jac <- function(lambda) {
      j <- diag(1, length(lambda), length(lambda))
      j
    }

    opt <-
      auglag(
        par = initial,
        fn = LogLikelihood,
        hin = hi,
        heq = he,
        heq.jac = he.jac,
        hin.jac = hi.jac,
        control.outer = list(trace = FALSE)
      )

    if (any(opt$ineq < 0) | any(round(opt$par, 1) == 1.0)) {
      ll <- '.'
      hessian <- '.'
      lambdahat <- '.'
      inf <- list(warning = 2)
    } else{
      hessian <- -opt$hessian
      ll <- -opt$value
      lambdahat <- opt$par

      inf <- Inference(hessian, lambdahat)
    }

    if (inf$warning == 2) {
      results[[i + k]] <- '.'
      results[[(i + k + 1)]] <- '.'
      names(results)[[i + k]] <-
        paste("Algorithm did not reach a solution with the constraints imposed.",
              sep = " ")
    } else if (inf$warning == 1) {
      results[[i + k]] <-
        output.table(lambdahat, inf$se, inf$zstat, inf$pvalue)
      results[[(i + k + 1)]] <- ll
      names(results)[[i + k]] <-
        paste("Equation",
              i,
              "-",
              "Hessian is singular, cannot compute standard errors.",
              sep = " ")
      names(results)[[(i + k + 1)]] <- paste("LogLik", i, sep = " ")
    } else{
      results[[i + k]] <-
        output.table(lambdahat, inf$se, inf$zstat, inf$pvalue)
      results[[(i + k + 1)]] <- ll
      names(results)[[i + k]] <- paste("Equation", i, sep = " ")
      names(results)[[(i + k + 1)]] <- paste("LogLik", i, sep = " ")
    }

    j = j + nr.eq - 1
    k = k + 1
  }

  return(results)

}
