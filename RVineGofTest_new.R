##' RVineGofTest_new
##' 
##' Function to calculate goodness-of-fit (GOF) test statistics proposed in
##' 'Generalized Information Matrix Tests for Copulas'
##' The function is inspired by the RVineGofTest function of the VineCopula
##' package.
##' The package is required.
##' 
##' @author Dr. Ulf Schepsmeier
##' 
##' @param data matrix; copula data
##' @param RVM R-vine matrix object from the VineCopula package
##' @param method string; which GOF test statistic has to be calculated
##' "White2" = Determinant White Test
##' "White3" = Trace White Test
##' "IR2" = Log Determinant IR Test
##' "log_trace" = Log Trace IMT
##' "log_GAIC" = Log Generalized Akaike Information Criterion (GAIC) IMT
##' "eigen" = Log Eigenspectrum IMT 
##' "eigen2" = Eigenvalue Test 
##' @param log boolean; for the "IR2" test
##' 
##' @return test statistic
##' 
##' @note only a preliminary stage of the test statistics can be computed. 
##' Meaning that the necessary variance (-covariance matrix) has to be approximated
##' by bootstrapping and the test statistics have to be corrected.
##' see correctTestStatistic()
##' 
##' @seealso correctTestStatistic
##' 

RVineGofTest_new2 <- function(data,
                              RVM,
                              method = "eigen",
                              log = FALSE){
  
  if (any(!(RVM$family %in% c(0, 1:6, 13, 14, 16, 23, 24, 26, 33, 34, 36))))
    stop("Copula family not implemented.")
  
  if (is.vector(data)) {
    data <- t(as.matrix(data))
  } else{
    data <- as.matrix(data)
  }
  
  if (any(data > 1) || any(data < 0)) stop("Data has be in the interval [0,1].")
  N <- dim(data)[1]
  d <- dim(data)[2]
  
  out <- RVineHessian(data, RVM)
  C <- out$der
  H <- out$hessian
  
  Hbar <- H / N
  Cbar <- C / N
  p <- dim(Hbar)[1]
  
  if (method == "White2"){
    
    return(sum(diag(Hbar + Cbar)))
    
  }else if (method == "White3"){
    
    return(det(Hbar + Cbar))
    
  }else if (method == "IR2"){
    
    out <- RVineHessian(data, RVM)
    C <- out$der
    H <- out$hessian
    Z <- solve(-C, H)
    Test <- det(Z)
    
    if (log == TRUE){
      Test <- log(Test)
    }
    return(Test)		# it is s, and not the test statistic T
    
  }else if (method == "log_trace"){
    
    Z1 <- solve(-Hbar)
    Z2 <- solve(Cbar)
    Test <- log(sum(diag(Z1))) - log(sum(diag(Z2)))
    return(Test)
    
  }else if (method == "log_GAIC"){
    
    Z1 <- eigen(-Hbar, only.values = TRUE)$values
    Z2 <- eigen(Cbar, only.values = TRUE)$values
    Test <- log(1 / p * sum(1 / Z1 * Z2))
    return(Test)
    
  }else if (method == "eigen"){
    
    Z1 <- eigen(-Hbar, only.values = TRUE)$values
    Z2 <- eigen(Cbar, only.values = TRUE)$values
    Z3 <- log(Z2) - log(Z1)
    return(Z3)
    
  }else if (method == "eigen2"){
    
    out <- RVineHessian(data, RVM)
    C <- out$der
    H <- out$hessian
    Z <- solve(-H, C)
    Z1 <- eigen(Z, only.values = TRUE)$values - 1
    return(Z1)
    
  }
  
}
