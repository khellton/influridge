#' Influence of single observations on the choice of the penalty parameter in ridge regression
#'
#' This function produces a curve plot showing the optimal tuning paramer when up- or downweighting each observation.
#' @param X numerical input matrix, of dimension nobs x nvars; each row is an observation vector.
#' @param y numerical output vector; response variable.
#' @param nw scalar, number of grid points for curve.
#' @param max.weight scalar, the maximum weight for the curve.
#' @param noExpand scalar, number of expanders to be highlighted (less than 5).
#' @param noShrink scalar, number of shrinkers to be highlighted (less than 5).
#' @param degreeFreedom logical, should the degrees of freedom be plotted instead of the tuning parameter (default = FALSE).
#' @param control.list list, a list of control parameters for the optim function. See 'Details' under graphics::optim.
#' @export
#' influridge
#' @examples
#' p <- 5
#' n <- 20
#' sigma <- 1
#' beta <- rep(1, p)
#'
#' ## Simulating design matrix, X
#' set.seed(556)
#' X <- matrix(rnorm(n * p), n, p)
#'
#' ## Simulate outcome vector, Y
#' y <- X %*% beta + rnorm(n, 0, sigma)
#'
#' ## Plot curves (no highlighted shrinkers/expanders)
#' influridge(X, y)
#'
#' ## Adding a large positive residual to observation 10 creates an influential shrinker
#' y[10] <- y[10] + 3
#' influridge(X, y, noShrink = 1, nw = 20)
#'
#' ## Plot degrees of freedom
#' influridge(X, y, noShrink = 1, nw = 20, degreeFreedom = TRUE)
#'
#' \dontrun{
#' ## Make plot for Body Fat dataset
#' require(mfp)
#' data(bodyfat)
#' X <- bodyfat[, 6:17] # Omit non-continous age variable
#' y <- bodyfat$siri
#' n <- dim(X)[1]
#'
#' X <- scale(X, center = FALSE) # Scale data
#' X <- cbind(rep(1, n), X) # Add intercept to design matrix
#'
#' influridge(X, y, noShrink = 1, noExpand = 1, degreeFreedom = TRUE)
#' }
influridge <- function(X, y, nw = 40, max.weight = 4, noExpand = 0, noShrink = 0, degreeFreedom = FALSE, control.list = list(factr = 1e-4)) {
  if (noShrink > 5) {
    print("Number of highlighted shrinkers must be 5 or less")
  }
  if (noExpand > 5) {
    print("Number of highlighted expanders must be 5 or less")
  }

  n <- dim(X)[1]
  weights <- seq(0, max.weight, length.out = nw) # Weight grid
  lambdaMatrix <- matrix(, n, nw)

  startLambda <- stats::optim(
    par = 1, tuning_cv_svd, w = rep(1, n), svd.int = svd(X), y.int = y,
    lower = -Inf, upper = Inf,
    method = "L-BFGS-B", control = control.list
  )$par

  for (i in 1:nw) {
    for (j in 1:n) {
      w <- rep(1, n)
      w[j] <- weights[i] # Weight to give observation j

      # find optimal lambda
      lambdaMatrix[j, i] <- stats::optim(
        par = startLambda, tuning_cv_svd, w = w / sum(w), svd.int = svd(X), y.int = y,
        lower = 0, upper = Inf,
        method = "L-BFGS-B", control = control.list
      )$par
    }
  }

  # Set
  lwt <- rep(1, n)
  col <- rep("grey", n)
  plotIndex <- rep("", n)

  sortIndexEnd <- sort(lambdaMatrix[, dim(lambdaMatrix)[2]],
    decreasing = FALSE, index.return = TRUE
  )$ix

  if (noExpand > 0) {
    lwt[sortIndexEnd[1:(noExpand)]] <- 3
    col[sortIndexEnd[1:(noExpand)]] <- "blue"
    plotIndex[sortIndexEnd[1:(noExpand)]] <- as.character(sortIndexEnd[1:(noExpand)])
  }

  if (noShrink > 0) {
    lwt[sortIndexEnd[(n + 1 - noShrink):n]] <- 3
    col[sortIndexEnd[(n + 1 - noShrink):n]] <- "red"
    plotIndex[sortIndexEnd[(n + 1 - noShrink):n]] <- as.character(sortIndexEnd[(n + 1 - noShrink):n])
  }

  if (degreeFreedom == FALSE) {
    graphics::matplot(weights, t(lambdaMatrix),
      type = "l",
      ylim = c(min(lambdaMatrix), max(lambdaMatrix)),
      xlim = c(0, max(weights) + 0.1),
      xlab = "Weight of observation",
      ylab = "Tuning parameter",
      lty = 1,
      lwd = lwt,
      xaxs = "i", yaxs = "i",
      cex.lab = 1.7, mgp = c(2.8, 1, 0),
      cex.axis = 1.5, col = col
    )
    graphics::axis(4,
      at = lambdaMatrix[, dim(lambdaMatrix)[2]], labels = plotIndex,
      las = 2, gap.axis = -1, tick = FALSE, cex.axis = 1.5, hadj = 0.6
    )
    graphics::abline(v = 1, lty = 2, col = "gray")
  } else {
    svd.df <- svd(X)
    df <- apply(lambdaMatrix, c(1, 2), function(lam) sum(svd.df$d^2 / (svd.df$d^2 + lam)))
    graphics::matplot(weights, t(df),
      type = "l",
      ylim = c(min(df), max(df)),
      xlim = c(0, max(weights) + 0.1),
      xlab = "Weight of observation",
      ylab = "Tuning parameter",
      lty = 1,
      lwd = lwt,
      xaxs = "i", yaxs = "i",
      cex.lab = 1.7, mgp = c(2.8, 1, 0),
      cex.axis = 1.5, col = col
    )
    graphics::axis(4,
      at = df[, dim(df)[2]], labels = plotIndex,
      las = 2, gap.axis = -1, tick = FALSE, cex.axis = 1.5, hadj = 0.6
    )
    graphics::abline(v = 1, lty = 2, col = "gray")
  }
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

tuning_cv_svd <- function(lambda, w, svd.int, y.int) {
  H <- svd.int$u %*% diag(svd.int$d^2 / (svd.int$d^2 + lambda)) %*% t(svd.int$u)
  e <- (diag(length(y.int)) - H) %*% y.int
  return(mean(w * (e / (1 - diag(H)))^2))
}
