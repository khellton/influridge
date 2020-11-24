#' Influence of single observations on the choice of the penalty parameter in ridge regression
#'
#' This function produces curve plot showing the optimal tuning paramer when up- or downweighting each observation
#' @param X numerical input matrix, of dimension nobs x nvars; each row is an observation vector.
#' @param y numerical output vector; response variable.
#' @param nw scalar, number of grid points for curve.
#' @param max.weight scalar, the maximum weight for the curve.
#' @export
#' plot_influridge

plot_influridge <- function(X, y, nw = 100, max.weight = 4, control.list = list(factr = 1e-4)){
  n <- dim(X)[1]
  weights <- seq(0, max.weight, length.out = nw) # Weight grid
  lambdaMatrix <- matrix(, n, nw)

  startLambda <- stats::optim(
    par = 1, tuning.cv.svd, w = rep(1,n), svd.int = svd(X), y.int = y,
    lower = -Inf, upper = Inf,
    method = "L-BFGS-B", control = control.list
  )$par

  for (i in 1:nw) {
    for (j in 1:n) {
      w <- rep(1, n)
      w[j] <- weights[i] # Weight to give observation j

            # find optimal lambda
      lambdaMatrix[j, i] <- stats::optim(
        par = startLambda, tuning.cv.svd, w = w / sum(w), svd.int = svd(X), y.int = y,
        lower = -Inf, upper = Inf,
        method = "L-BFGS-B", control = control.list
      )$par
    }
  }

  graphics::matplot(weights, t(lambdaMatrix),
          type = "l", ylim = c(min(lambdaMatrix), max(lambdaMatrix)), xlim = c(0, max(weights) + 0.1),
          xlab = "Weight of observation", ylab = "Tuning parameter", lty = 1, xaxs = "i", yaxs = "i",
          cex.lab = 1.7, mgp = c(2.8, 1, 0), cex.axis = 1.5, col = "grey"
  )
}

#' Weighted leave-one-out crossvalidation error
#'
#' This function produces the leave-one-out squared cross-validation error for a given tuning parameter and weight.
#' @param lambda scalar, ridge tuning parameter.
#' @param w vector of length n, weigths for all observations.
#' @param svd.int svd object, singular value decomposition of the input matrix.
#' @param y.int numerical output vector; response variable.
#' @export
#' tuning_cv_svd

tuning_cv_svd <- function(lambda, w, svd.int,y.int) {
  H <- svd.int$u %*% diag(svd.int$d^2 / (svd.int$d^2 + lambda)) %*% t(svd.int$u)
  e <- (diag(length(y)) - H) %*% y.int
  return(mean(w * (e / (1 - diag(H)))^2))
}
