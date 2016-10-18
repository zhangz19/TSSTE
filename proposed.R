

# main function for the simulation

## convert mat file to R data file  
#library(R.matlab)
#pvec <- c(3,5,10)
#for(i in pvec){
#  N <- readMat(paste('unitsphere',i,'.mat', sep=''))$N
#  save(file=paste('unitsphere',i,sep=''), N)
#}

repl <- as.numeric(commandArgs(TRUE))
set.seed(repl)

require(methods)
library(mvtnorm)
library(cramer)

# ===========================  basic setup
# nsim <- 1e2  # number of replicates for calculting size/power
dyn.load('proposed.so') # load .so
M <- 5e3 # bootstrap size 
alpha <- 0.05 # significance level
ndesign <- length(designs <- c(1,2,3,4))
np <- length(pvec <- c(3,5,10))
# if(hpcrun) np <- length(pvec <- c(3,5,10)[case])
verbose <- FALSE

# for location and dispersion
dvec1 <- 0 #c(0, seq(.05,1.9, by=.3))
ndelta <- length(dvec1)

# for dependent case 
dvec2 <- 0

# ===========================  main function
ntest <- 2 # proposed, and cramer
pvec0 <- c(3,5)
designs0 <- c(4)

matAll <- numeric()
ptm <- Sys.time() # record the CPUtime        
for(p in pvec){
 mat <- matrix(NA, ndelta, ntest*ndesign)
 load(paste('unitsphere',p,sep='')) # N 
 storage.mode(N) <- 'double'
 
 for(design in designs){
   for(i in 1:ndelta){
    for(j in 1:ntest){ 
       
      if(verbose) cat(paste('p=',p, ', rep=',repl, ', design=',design, ', i=',i, ', j=',j,':\n',sep=''))
   
      # ===========================  1. data generation
      n <- 80; m <- 60; 
      
      if(design == 1){ #location
         x <- matrix(rnorm(n*p), nrow=n)
         y <- matrix(rnorm(m*p), nrow=m)
         y[,1] <- y[,1] + dvec1[i]
      }
      if(design == 2){ #dispersion
         x <- matrix(rnorm(n*p), nrow=n)
         if(i==1) y <- matrix(rnorm(m*p), nrow=m)
         if(i!=1) y <- dvec1[i]*matrix(rnorm(m*p), nrow=m)
      }
      if(design == 3){ #dependence alternatives: MVN
         x <- matrix(rnorm(n*p), nrow=n)
         A <- diag(p)
         A[1,1:2] <- sqrt(c(1-dvec2[i], dvec2[i]))
         A[2,1:2] <- sqrt(c(dvec2[i], 1-dvec2[i]))
         y <- matrix(rnorm(m*p), nrow=m)%*%t(A)
      }
      if(design == 4){ #dependence alternatives: MVT
         x <- rmvt(n, sigma = diag(p), df = 5)
         A <- diag(p)
         A[1,1:2] <- sqrt(c(1-dvec2[i], dvec2[i]))
         A[2,1:2] <- sqrt(c(dvec2[i], 1-dvec2[i]))
         y <- rmvt(m, sigma = diag(p), df = 5)%*%t(A)
      }
      
      # ===========================  2. run the test
      if(j==1){ # our proposed methods
      n <- nrow(x); m <- nrow(y); p <- ncol(x); len <- nrow(N)
      e <- matrix(rnorm(n*M), nrow=n)
      stats <- numeric(1);   T1 <- numeric(M)
      
      storage.mode(n) <- storage.mode(m) <- storage.mode(p) <- storage.mode(len) <- storage.mode(M) <- "integer"
      storage.mode(x) <- storage.mode(y) <- storage.mode(N) <- storage.mode(e) <- storage.mode(stats) <- storage.mode(T1) <- "double"
      
      junk <- .Call("proposed", x,y,n,m,p,N,len,e,M, stats, T1)
      CV <- quantile(T1, probs=1-alpha)
      sig <- as.integer(stats >= CV) # =1: --> H_0: F=G rejected
      }
      
      if(j==2) # cramer test
        sig <- cramer.test(x,y,sim="permutation",replicates=5e3)$result
      
        
      mat[i, (design-1)*ntest + j] <- sig
      
    }
   }
  }
  matAll <- cbind(matAll, mat)
  
}      
print(Sys.time()-ptm)    


mat <- matAll
save(file=paste('out',repl,sep=''), mat)

