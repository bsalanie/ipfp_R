################################################################################
####         the IPFP algorithm for the homoskedastic logit model            ###
####              (the marriage market of Choo and Siow  2006)               ###
####                                                                         ###
####           Alfred Galichon and Bernard Salanie                           ###
################################################################################


################################################################################
####     resmu = solveipfp(phi, p, q, maxiter=10000, epsconv=1e-6)          ####
####        phi: (nx, ny) matrix of joint surpluses                         ####
####           [assumes 0 joint surplus for singles]                        ####                      ####
####        p: nx-vector of margins for men                                 ####
####        q: ny-vector of margins for women                               ####
####                                                                        ####
####     returns resmu = list(status, muopt, niter, diter)                  ####
####      status = 0 if all went well                                       ####
####      muopt is a list (muxy, mux0, mu0y)                                ####  
####      niter the number of iterations                                    ####
####      diter how much the last two iterates differed                     ####
####                                                                        ####
####     iterations stop when they exceed maxiter (1000 by default)         ####
####      or when diter < epsconv (1e-6 by default)                         ####
################################################################################


## one IPFP step: updates (tx, ty) given z=exp(Phi/2) and margins
ipfpiter <- function(txy, z, p, q) {
  tz <- t(z)
  tx <- txy[[1]]
  lx <- crossprod(z, tx)
  ty1 <- (sqrt(lx*lx + 4*q)-lx)*0.5
  ly1 <- crossprod(tz, ty1)
  tx1 <- (sqrt(ly1*ly1 + 4*p)-ly1)*0.5
  ## we return
  list(tx1, ty1)
}


## weighted Euclidean distance between two iterations
##  (could be replaced with anything reasonable)
distxy <- function(txy1, txy, p, q) {
  tx1 <- txy1[[1]]
  tx <- txy[[1]]
  ty1 <- txy1[[2]]
  ty <- txy[[2]]
  distx2 <- sum(p*(tx-tx1)*(tx-tx1))/sum(p)
  disty2 <- sum(q*(ty-ty1)*(ty-ty1))/sum(q)
  ## we return
  sqrt(distx2 + disty2)
}


## computes the intermediate matrix z = exp(Phi/2)
compz <- function(phi) {
  z <- exp(phi*0.5)
  ## we return
  z
}

## solve IPFP by iterating over the square roots of mux0 and mu0y
solveipfp <- function(phi, p, q, epsconv=1e-6, maxiter=1000) {
  nx <- NROW(phi)
  ny <- NCOL(phi)
  z <- compz(phi)
  sumz <- sum(z)
  nIndivs <- sum(p)+sum(q)
  bigc <- sqrt(nIndivs/(2.0*sumz+nx+ny))
  ## initialize tx=sqrt(mux0) and ty=sqrt(mu0y)
  tx <- bigc + numeric(nx)
  ty <- bigc + numeric(ny)
  txy <- list(tx,ty)
  ## iterate
  diter <- +Inf
  niter <- 0
  while ((diter > epsconv) && (niter < maxiter)) {
    txy1 <- ipfpiter(txy, z, p, q)
    diter <- distxy(txy1, txy, p, q)
    niter <- niter + 1
    ##cat("Iter ", niter, ": diter=", diter, "\n")
    txy <- txy1
  }
  if (diter < epsconv) {
    status <- 0
  } else {
     status <- 1
  }
  tx <- as.vector(txy[[1]])
  ty <-  as.vector(txy[[2]])
  muXY <- (tx %o% ty)*z
  muX0 <- p - rowSums(muXY)
  mu0Y <- q - colSums(muXY)
  ## we return
  list(status, list(muXY, muX0, mu0Y), niter, diter)
}

############################################################
######                  Example                       ######
############################################################

runExample <- F

if (runExample) {
  nx <- 7
  ny <- 5
  p <- c(2, 4, 6, 3, 5, 7, 5)
  q <- c(3, 9, 7, 2, 4)
  x <- seq(nx)
  y <- seq(ny)
  dxy <-  outer(x, y, "-")
  Phi <- -dxy*dxy
  resmu <- solveipfp(Phi, p, q)
  status <- resmu[[1]]
  if (status == 0) {
    muopt <- resmu[[2]]
    niter <- resmu[[3]]
    diter <- resmu[[4]]
    cat("Done in ", niter, " iterations; up to ", diter, "\n")
    muxy <- muopt[[1]]
    print(muxy)
  } else {
    cat("Failed to converge in ", niter, "iterations; distance ", diter, "\n")
  }
}

## expected results:

## Done in  11  iterations; up to  2.630866e-07
##              [,1]         [,2]        [,3]       [,4]         [,5]
## [1,] 7.242520e-01 7.762220e-01 0.133029546 0.00362425 0.0001638729
## [2,] 6.121753e-01 1.783473e+00 0.830851061 0.06153010 0.0075626028
## [3,] 2.043798e-01 1.618539e+00 2.049624565 0.41260384 0.1378514666
## [4,] 1.233869e-02 2.656127e-01 0.914312141 0.50032006 0.4543813778
## [5,] 6.650319e-04 3.891497e-02 0.364130249 0.54163244 1.3371242327
## [6,] 1.062023e-05 1.689284e-03 0.042967216 0.17373191 1.1658459689
## [7,] 3.990231e-08 1.725286e-05 0.001192862 0.01311074 0.2391570180
