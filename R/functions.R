#' Influence of single observations on the choice of the penalty parameter in ridge regression
#'
#' This function produces curve plot showing the optimal tuning paramer when up- or downweighting each observation
#' @param X numerical input matrix, of dimension nobs x nvars; each row is an observation vector.
#' @param y numerical output vector; response variable.
#' @param nw scalar, number of grid points for curve.
#' @param max.weight scalar, the maximum weight for the curve.
#' @param noExpand scalar, number of expanders to be highlighted (less than 5). 
#' @param noShrink scalar, number of shrinkers to be highlighted (less than 5). 
#' @param degreeFreedom logical, should the degrees of freedom be plotted instead of the tuning parameter (default = FALSE).
#' @param control.list list, a list of control parameters for the optim function. See 'Details' under graphics::optim. 
#' 
#' @examples 
#' 
#' p <- 5 
#' n <- 50 
#' sigma <- 1 
#' beta <- rep(0.1,p)
#' 
#' Simulating design matrix, X
#' setseed(556)
#' X <- matrix(rnorm(n*p),n,p)
#'
#' Simulate outcome vector, Y
#' y <- X %*% beta + rnorm(n, 0, sigma)
#' 
#' Plot curves (no highlighted shrinkers/expanders)
#' influridge(X,y)
#' 
#' Simulate outcome vector Y, adding a large negative residual to observation 7.
#' y <- X %*% beta + c(rnorm(6, 0, sigma), -4, rnorm(n - 7, 0, sigma))
#' influridge(X,y,noShrink = 1)
#' 
#' Plot degrees of freedom
#' y <- X %*% beta + c(rnorm(6, 0, sigma), -4, rnorm(n - 7, 0, sigma))
#' influridge(X,y,noShrink = 1,degreeFreedom == TRUE)
#' 
#' @export
#' influridge

influridge <- function(X, y, nw = 100, max.weight = 4, 
                       noExpand = 0, noShrink = 0, 
                       degreeFreedom = FALSE,
                       control.list = list(factr = 1e-8)){
  
  if(noShrink > 5){print('Number of highlighted shrinkers must be 5 or less') }
  if(noExpand > 5){print('Number of highlighted expanders must be 5 or less') }
  
  n <- dim(X)[1]
  weights <- seq(0, max.weight, length.out = nw) # Weight grid
  lambdaMatrix <- matrix(, n, nw)

  startLambda <- stats::optim(
    par = 1, tuning_cv_svd, w = rep(1,n), svd.int = svd(X), y.int = y,
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
        lower = -Inf, upper = Inf,
        method = "L-BFGS-B", control = control.list
      )$par
    }
  }
  
  #Set 
  lwt <- rep(1,n)
  col <- rep("grey",n)
  
  sortIndexEnd <- sort(lambdaMatrix[,dim(lambdaMatrix)[2]],
                  decreasing = FALSE,index.return = TRUE)$ix
  
  if(noShrink > 0){
    lwt[sortIndexEnd[1:noExpand]] <- 3
    col[sortIndexEnd[1:noExpand]] <- 'blue'  }
  
  if(noShrink > 0){
    lwt[sortIndexEnd[(n-noShrink):n]] <- 3
    col[sortIndexEnd[(n-noShrink):n]] <- 'red'
  }
  
  if(degreeFreedom == FALSE){
    graphics::matplot(weights, t(lambdaMatrix), type = "l", 
                    ylim = c(min(lambdaMatrix), max(lambdaMatrix)), 
                    xlim = c(0, max(weights) + 0.1), 
                    xlab = "Weight of observation", 
                    ylab = "Tuning parameter", 
                    lty = 1, 
                    lwd = lwt,
                    xaxs = "i", yaxs = "i",
                    cex.lab = 1.7, mgp = c(2.8, 1, 0), 
                    cex.axis = 1.5, col = col
    )} else {
    svd.df <- svd(X)
    df <- apply(lambdaMatrix, c(1, 2), function(lam) sum(svd.df$d^2 / (svd.df$d^2 + lam)))
    graphics::matplot(weights, t(df), type = "l", 
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

tuning_cv_svd <- function(lambda, w, svd.int,y.int) {
  H <- svd.int$u %*% diag(svd.int$d^2 / (svd.int$d^2 + lambda)) %*% t(svd.int$u)
  e <- (diag(length(y)) - H) %*% y.int
  return(mean(w * (e / (1 - diag(H)))^2))
}
