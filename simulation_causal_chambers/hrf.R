#' Quadratic-inverse shrinkage
#'
#' Nonlinear shrinkage derived under Frobenius loss and its two cousins,
#' Inverse Stein’s loss and Minimum Variance loss, called quadratic-inverse
#' shrinkage (QIS). See Ledoit and Wolf (2022, Section 4.5).
#'
#' @param data (n*p): raw data matrix of n iid observations on p random
#' variables
#' @param k If k < 0, then the algorithm demeans the data by default, and
#' adjusts the effective sample size accordingly. If the user inputs k = 0,
#' then no demeaning takes place; if user inputs k = 1, then it signifies that
#' the data data have already been demeaned.
#'
#' @return sigmahat (p*p): the QIS covariance matrix estimate. An object of
#' class \code{matrix}.
#'
get.cov.qis <- function(data, k = -1) {
  dim.data <- dim(data)
  n <- dim.data[1]
  p <- dim.data[2]
  if (k < 0) {
    # demean the data and set k = 1
    data <- scale(data, scale = FALSE)
    k <- 1
  }
  n.free <- n - k
  c <- p / n.free
  sample <- (t(data) %*% data) / n.free
  sample <- (t(sample) + sample) / 2
  spectral <- eigen(sample, symmetric = TRUE)
  lambda <- spectral$values[p:1]
  u <- spectral$vectors[, p:1]
  h <- min(c^2, 1 / c^2)^0.35 / p^0.35
  inv.lambda <- 1 / lambda[max(1, p - n.free + 1):p]
  l.j <- matrix(rep(inv.lambda, each = min(p, n.free)), nrow = min(p, n.free))
  l.j.i <- l.j - t(l.j)
  theta <- rowMeans(l.j * l.j.i / (l.j.i^2 + h^2 * l.j^2))
  h.theta <- rowMeans(l.j * (h * l.j) / (l.j.i^2 + h^2 * l.j^2))
  a.theta2 <- theta^2 + h.theta^2
  if (p <= n.free) {
    delta <- 1 / ((1 - c)^2 * inv.lambda + 2 * c * (1 - c) * inv.lambda *
                    theta + c^2 * inv.lambda * a.theta2)
  } else {
    delta0 <- 1 / ((c - 1) * mean(inv.lambda))
    delta <- c(rep(delta0, p - n.free), 1 / (inv.lambda * a.theta2))
  }
  delta.qis <- delta * (sum(lambda) / sum(delta))
  sigma.hat <- u %*% diag(delta.qis) %*% t(u)
  return(sigma.hat)
}





hrf.orig <- function(x, y, rf.fit, num.iter = NULL, kappa = 2) {
  
  # Random forest predictions on all trees
  rf.predictions.all <- predict(
    rf.fit,
    x,
    predict.all = TRUE
  )$predictions[,"1",]
  
  # Calculate residuals
  y.numeric <- as.numeric(y)-1
  rf.residuals <- y.numeric - rf.predictions.all
  
  # Calculate covariance matrix and mean vector of residuals
  mean.vector <- colMeans(rf.residuals)
  cov.matrix <- get.cov.qis(rf.residuals)
  
  # Get optimized tree weights
  w <- CVXR::Variable(ncol(cov.matrix))
  objective <- CVXR::square((t(w) %*% mean.vector)) +
    CVXR::quad_form(w, cov.matrix)
  constraints <- list(
    sum(w) == 1,
    sum(abs(w)) <= kappa
  )
  prob <- CVXR::Problem(CVXR::Minimize(objective), constraints)
  solution <- CVXR::solve(prob, num_iter = num.iter)
  
  hedgedrf <- structure(
    list(
      tree.weights = solution$getValue(w),
      rf.fit = rf.fit
    ),
    class = "hedgedrf"
  )
  return(hedgedrf)
}


hrf.ood <- function(x, y, num.iter = NULL, kappa = 2) {
  
  num.env <- length(levels(x$Env))
  num.build.trees <- ceiling(num.env/2)
  build.trees.envs <- base::sample(levels(x$Env), size = num.build.trees, replace = FALSE)
  
  idx.trees <- which(x$Env %in% build.trees.envs)
  idx.weights <- which(!(x$Env %in% build.trees.envs))
  
  rf.fit <- ranger::ranger(x = x[idx.trees, 1:d], y = y[idx.trees], probability = TRUE, num.threads = 0, num.trees = 1000)
  
  # Random forest predictions on all trees
  rf.predictions.all <- predict(
    rf.fit,
    x[idx.weights, 1:d],
    predict.all = TRUE
  )$predictions[,"1",]
  
  # Calculate residuals
  y.numeric <- as.numeric(y[idx.weights])-1
  rf.residuals <- y.numeric - rf.predictions.all
  
  # Calculate covariance matrix and mean vector of residuals
  mean.vector <- colMeans(rf.residuals)
  cov.matrix <- get.cov.qis(rf.residuals)
  
  # Get optimized tree weights
  w <- CVXR::Variable(ncol(cov.matrix))
  objective <- CVXR::square((t(w) %*% mean.vector)) +
    CVXR::quad_form(w, cov.matrix)
  constraints <- list(
    sum(w) == 1,
    sum(abs(w)) <= kappa
  )
  prob <- CVXR::Problem(CVXR::Minimize(objective), constraints)
  solution <- CVXR::solve(prob, num_iter = num.iter)
  
  hedgedrf <- structure(
    list(
      tree.weights = solution$getValue(w),
      rf.fit = rf.fit
    ),
    class = "hedgedrf"
  )
  return(hedgedrf)
}




#' hedgedrf prediction
#'
#' @param object hedgedrf \code{hedgedrf} object.
#' @param data data New test data of class \code{data.frame} or
#' \code{gwaa.data} (GenABEL).
#' @param ... Additional arguments to pass to the \code{predict.ranger}
#' function.
#' @return The hedged random forest predictions. An object of class \code{matrix}.
#' @export
predict.hedgedrf <- function(object, data, ...) {
  # Random forest predictions on test data
  rf.predictions <- predict(
    object$rf.fit,
    data,
    predict.all = TRUE,
    probability = TRUE,
  )$predictions[, "1", ]
  hedgedrf.predictions <- as.matrix(rf.predictions) %*% object$tree.weights
  return(hedgedrf.predictions)
}





