##' correctTestStatistic
##' 
##' The test statistics we get from RVineGofTest_new are only a preliminary 
##' stage of the test statistics. They have to corrected with the variance
##' (covariance matrix).
##' 
##' @author Dr. Ulf Schepsmeier
##' 
##' @param I_b0 matrix or vector of test statistics; rows are samples of the d-dimensional
##' test statistic
##' @param N number of observations (not the number of boostraps)
##' @param method string; which GOF test statistic has to be calculated
##' "White2" = Determinant White Test
##' "White3" = Trace White Test
##' "IR2" = Log Determinant IR Test
##' "log_trace" = Log Trace IMT
##' "log_GAIC" = Log Generalized Akaike Information Criterion (GAIC) IMT
##' "eigen" = Log Eigenspectrum IMT 
##' "eigen2" = Eigenvalue Test 
##' 
##' @return I_b0_new vector of test statistics
##' 

correctTestStatistic <- function(I_b0, N, method){
  
  if(method == "eigen") return(correctTestStatisticEigen(I_b0 = I_b0, N = N))
  else if(method == "eigen2") return(correctTestStatisticEigen2(I_b0 = I_b0, N = N))
  else {
    sigma2_b0 <- var(I_b0)
    mu_b0 <- mean(I_b0)
    I_b0_new <- ((I_b0-mu_b0)/sqrt(sigma2_b0))^2
    return(I_b0_new)
  }
  
}



##' correctTestStatisticEigen
##' 
##' For the 'Log eigenspectrum test' the calculated test statistic we get from the 
##' RVineGofTest_new function is 
##' log(Z2)-log(Z1), where Z1=eigen(-Hbar, only.values=TRUE)$values and
##' Z2=eigen(Cbar, only.values=TRUE)$values
##' This has to be corrected by the variance-covariance matrix
##' It can only be approximated.
##' Thus as input we expect boostrapped 'Log eigenspectrum test' test statistics
##' Each test statistic we get from RVineGofTest_new2 is a vector!
##' 
##' @author Dr. Ulf Schepsmeier
##' 
##' @param I_b0 matrix of test statistics; rows are samples of the d-dimensional
##' test statistic
##' @param N number of observations (not the number of boostraps)
##' 
##' @return I_b0_new vector of test statistics
##' 

correctTestStatisticEigen <- function(I_b0, N){
  
  Boot <- nrow(I_b0)
  nn <- sqrt(N)
  sigma2_b0 <- cov(I_b0)
  
  I_b0_new <- rep(0, Boot)
  for(i in 1:Boot){
    I_b0_new[i] <- nn*(I_b0[i,])%*%solve(sigma2_b0)%*%(I_b0[i,])
  }
  
  return(I_b0_new)
}



##' correctTestStatisticEigen2
##' 
##' For the 'Eigenvalue test' the calculated test statistic we get from the 
##' RVineGofTest_new2 function is 
##' eigen(solve(-H,C), only.values=TRUE)$values-1
##' This has to be corrected.
##' It can only be approximated.
##' Thus as input we expect boostrapped 'Eigenvalue test' teststatistics
##' Each test statistic we get from RVineGofTest_new2 is a vector!
##' 
##' @author Dr. Ulf Schepsmeier
##' 
##' @param I_b0 matrix of test statistics; rows are samples of the d-dimensional
##' test statistic
##' @param N number of observations (not the number of boostraps)
##' 
##' @return I_b0_new vector of test statistics
##' 

correctTestStatisticEigen2 <- function(I_b0, N){
  
  Boot <- nrow(I_b0)
  nn <- sqrt(N)
  
  I_b0_new <- rep(0, Boot)
  for(i in 1:Boot){
    I_b0_new[i] <- nn*(I_b0[i,])%*%(I_b0[i,])
  }
  
  return(I_b0_new)
}