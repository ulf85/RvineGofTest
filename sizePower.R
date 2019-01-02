##' sizePower
##' 
##' Function to calculate/estimate the size and power of the GOF-tests
##' based on White and Golden et al.
##' 
##' @author Dr. ulf Schepsmeier
##' 
##' @param I_b0 vector of bootstrapped test statistics of the true model
##' @param I_bA1 vector of bootstrapped test statistics of the alternative model
##' 
##' @return list of size and power

sizePower <- function(I_b0, I_bA1)
{
  # p-values
  Boot <- length(I_b0)
  p_b0 <- rep(0, Boot)
  p_bA1 <- rep(0, Boot)
  
  for (j in 1:Boot)
  {
    p_b0[j] <- mean(I_b0 >= I_b0[j])
    p_bA1[j] <- mean(I_b0 >= I_bA1[j])
  }
  
  # steps
  xi_1 <- seq(0.001, 0.01, 0.001)
  xi_2 <- seq(0.015, 0.99, 0.005)
  xi_3 <- seq(0.991, 1, 0.001)
  xi <- c(xi_1, xi_2, xi_3)
  n <- length(xi)
  
  # distribution => size and power
  F_b0 <- rep(0, n)
  F_bA1 <- rep(0, n)
  
  for (j in 1:length(xi))
  {
    F_b0[j] <- mean(p_b0 <= xi[j])
    F_bA1[j] <- mean(p_bA1 <= xi[j])
  }
  
  
  out <- list(size = F_b0, power = F_bA1)
  return(out)
}

