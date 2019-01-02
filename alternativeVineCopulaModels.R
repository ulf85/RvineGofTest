# load the defined R-vine stored in a R-vine matrix objext (RVM)

library(VineCopula)

# Initial simulation run
simdata <- RVineSim(N = 1000, RVM = RVM)

# Estimate the three alternative vine copula models

# 1. Alternative: C-vine
RVM_CVine <- RVineStructureSelect(data  =simdata, 
                                  familyset = c(1,3:6,13,14,16,23,24,26,33,34,36), 
                                  type=1, # C-vine
                                  indeptest=TRUE, 
                                  trunclevel = 3) # to simplify we truncate at level three

# 2. Alternative: D-vine
# heuristic to get a D-vine: travelling salesman problem
library(TSP)
d <- dim(simdata)[2]
M <- 1 - abs(TauMatrix(simdata))
hamilton <- insert_dummy(TSP(M), label = "cut")
sol <- solve_TSP(hamilton, method = "repetitive_nn")
order <- cut_tour(sol, "cut")
DVM <- D2RVine(order, family = rep(0,d*(d-1)/2), par = rep(0,d*(d-1)/2))
RVM_DVine <- RVineCopSelect(simdata, familyset=c(1,3:6,13,14,16,23,24,26,33,34,36), 
                            DVM$Matrix, 
                            indeptest = TRUE, 
                            trunclevel=3)	

# 3. Alternative: R-vine with just Gaussian copulas
family2 <- matrix(0,d,d)
family2[lower.tri(family2)] <- 1
par2 <- matrix(0,d,d)
RVM_Gauss <- RVineMatrix(Matrix = RVM$Matrix, family = family2,
                         par = par2, par2 = par2)
RVM_Gauss <- RVineSeqEst(data = simdata, RVM = RVM_Gauss)$RVM			  

save.image(file = "Rvine_sample.RData")

