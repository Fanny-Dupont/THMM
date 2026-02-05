###############################################################################
#                    Code to Simulate Data and Fit THMM for the paper:        #
#      Estimating the distance at which narwhal respond to disturbance:       #
#               a penalized threshold hidden Markov model                     #
#                                                                             #
#      Authors: Fanny Dupont, Marianne Marcoux, Nigel E. Hussey,              # 
#                       Jackie Dawson, Marie Auger-Méthé                      #
###############################################################################

# Load required libraries
library(momentuHMM)
library(boot)
library(tidyverse)
library(pracma)
library(parallel)
library(dirmult)
library(LaMa)

source("SourceFunctions.R")


#######################
#### GENERATE DATA ####
#######################


nbStates <- N <- 3            # number of states
M <- 1                        # number of individuals
nobs <- n <- 5e3              # number of observation per individuals: 1e3, 3e3, 5e3, 1e4
obsVect = rep(list(nobs),M)   # number of observations per animal (length of obsVect must be a factor of M)

# state-dependent distribution, only simulated gamma
dist <- list(step = "gamma")

# set parameters for simulated data
## state-dependent parameters 
mu0 <- c(10, 4,1)            # mean of step length
shape0 <- c(12, 10, 1.5)     # shape parameter of step length
sd0 <- mu0/sqrt(shape0)      # sd of step length


## hidden process parameters 
Gamma = list()
Gamma[[1]] = matrix(c(0.9, 0.05, 0.05, # Strong persistence for baseline behaviour 
                      0.05, 0.9, 0.05, # Strong persistence for baseline behaviour 
                      0.05, 0.05, 0.9), ncol = 3, nrow = 3) # Strong persistence for baseline behaviour 

Gamma[[2]] = t(matrix(c(0.9, 0.05, 0.05, # Strong persistence for disturbed behaviour 
                        0.15, 0.7, 0.15, # Weaker persistence for disturbed behaviour 
                        0.15, 0.15, 0.7), ncol = 3, nrow = 3))  # Weaker persistence for disturbed behaviour 


# Generate time series of covariates
time <- seq(1, n, length.out = n)

generate_time_series <- function(time) {
  # A combination of sine and cosine functions to have periodicity
   y <- 20 + 10*(sin(time / 150) + cos(time / 650))
  return(y)
}

# Generate the time series data
U <- generate_time_series(time)

#### If univariate covariate ####
# Define mixture probability:
nu = sapply(U, function(x)ifelse(x<21,0,1))

# Simulate null model: no threshold: unique(nu) = 0 # uncomment to simulate model under the null
# nu = as.numeric(sapply(U, function(x)ifelse(x<min(U),1,0)))


#### Uncomment for bivariate covariate 
# Cat <- simulate_sticky_bernoulli(n, p_same = 0.9, x0 = 1)

# first thresold is 21, other is 30
# threshold <- ifelse(Cat == 1, 21, 30)

# Both null
# threshold <- ifelse(Cat == 1, max(U)+1, max(U)+1)

# One threshold active
# threshold <- ifelse(Cat == 1, 21, max(U)+1)

# Define mixture porbability:
# nu <- ifelse(U < threshold, 0, 1)

# Obtain bivariate covariate
# UCat = ifelse(Cat,U,0)
# UnCat =  ifelse(Cat==0,U,0)
# U = cbind(UCat,UnCat)





# Initial state probability
delta=c(1/3,1/3,1/3)
# Define the resulting TPM: as a mixture of two component, for each time
GammaT = array(dim=c(n,3,3))
for(i in (1:n)){
  GammaT[i,,] <- (1-nu[i]) * (Gamma[[1]]) + Gamma[[2]] * nu[i]
}



# Generate data
# set.seed(id)  # id of the job on CC
data = sim_HMM_gamma(n,3,GammaT,delta,mu0,sd0)

# Standardize covariate between 0 and 1
truethreshold = 1/((21-min(U))/(max(U)-min(U)))
U = (U-min(U))/(max(U)-min(U))
obsVect = rep(list(n),M)   # number of observations per animal (length of obsVect must be a factor of M)
T_m = rep(obsVect,M/length(obsVect)) 
covIndex = c(1,n+1) # indices marking the start of each new individual (only one individual in simulations)

# Matrix of covariate in the hidden process (multinomial logit link)

X1 = X2 = NULL
# X1 is associated with TPM of first component. If no covariate, create only the intercept.
if (is.null(X1)) {
  X1 <- matrix(1, nrow = n, ncol = 1, dimnames = list(NULL, "Intercept"))
} else {
  if (length(X1[,1])!=n){stop(paste("Wrong dimension for X1"))}
  X1 <- cbind(Intercept = 1, as.matrix(X1))
}

# X1 is associated with TPM of second component. If no covariate, create only the intercept.
if (is.null(X2)) {
  X2 <- matrix(1, nrow = n, ncol = 1, dimnames = list(NULL, "Intercept"))
} else {
  if (length(X2[,1])!=n){stop(paste("Wrong dimension for X2"))}
  X2 <- cbind(Intercept = 1, as.matrix(X2))
}



#######################
###### Fit THMM #######
#######################

Nrep = 15       # Number of random initial values explored
mod = list()    # List of fitted models   


fit_model <- function(i) {
  
  ###########################################################################################################################################
  #
  # Below is a list of all the parameters that are estimated during the procedure and require initial value (NB: lasso penalty strength
  # parameter does not require initial values but is estimated). 
  # Log transformations are applied to ensure positivity and numerical stability when needed.
  #
  # 
  #   logmu      - Log of the gamma state-dependent means (mu) of the data streams, of dimension N (Number of states)
  #   logsigma   - Log of the gamma state-dependent standard deviations (sigma), of dimension N (Number of states)
  #   beta1      - Parameters for Gamma1 (the first transition probability matrix) on the logit scale, of dimension N*N-1 *nbcov
  #   beta2temp  - Parameters for Gamma2 (the second transition probability matrix), on the logit scale, of dimension N*N-1*nbcov - 1 
  #                We exclude the first element (intercept associated with the first state), which is computed from `beta1` and a mini-
  #                -mum separation constraint (to ensure minimum difference between the two components). Thus, `beta2temp` has one 
  #                fewer element than `beta1`.
  #   delta1     - Initial state probability for component 1
  #   delta2     - Initial state probability for component 2
  #   lbetaThr   - Log of the threshold parameter(s), referred to as beta0 in the main manuscript.
  #   logK       - Log of the separation parameter K used to enforce a minimum difference between normalizing constants for state 1 in both regimes. Specifically,
  #                the normalizing sum for row (state) 1 in regime D is forced to be increased by a fixed offset 0.2 + exp(logK) (i.e., minimum 0.2) relative to regime B. 
  #                The value 0.2 is the default used in the code and can be modified.
  #
  #                Note that this offset is fixed in the algorithm, while the resulting separation on the probability scale depends on the normalization and is therefore
  #                data-driven. In practice, for our applications, this constraint typically induces differences in persistence probabilities greater than 0.15, although the offset 
  #                value (0.2) can be adjusted.
  #
  ###########################################################################################################################################

  # INITIALIZE PARAMETERS 
  
  skip_to_next <- FALSE                                 # Initialize test for error in optimization                           
  j <- 3                                                # Number of states of the THMM
  beta01 <- beta02 <- array(NA, dim = c(j, j))          # Initialize coefficients of multinomial logit link in TPMs
  
  # Generate initial parameters for Gamma1 and Gamma2, from Pohle et al. (2017)
  for (l in 1:j) {
    alpha <- rep(0.2 / (j - 1), j)
    alpha[l] <- 0.9
    alpha <- alpha * 10
    beta01[l, ] <- dirmult::rdirichlet(1, alpha)
    beta02[l, ] <- dirmult::rdirichlet(1, alpha)
  }
  
  # Initial TPM parameters (on logit scale)
  beta01 <- logit(beta01)
  beta01 <- t(beta01)
  beta01 <- matrix(beta01[col(beta01) != row(beta01)], ncol = j * (j - 1), nrow = 1)
  
  # Define initial parameters for state-dependent parameter
  mu0.step <- abs(c(rnorm(1, mean = 1, sd = 0.5), rnorm(1, mean = 4, sd = 1), rnorm(1, mean = 10, sd = 1)))
  sigma0.step <- abs(c(rnorm(1, mean = 2.5, sd = 0.5), rnorm(1, mean = 1, sd = 0.5), rnorm(1, mean = 1.5, sd = 0.5)))

  betaThr0 = runif(1, 0.5, 5)                   # Initial threshold parameter
  lbetaThr0 = log(betaThr0)                     # Take the log-value
  K = abs(rnorm(1,0,0.5))                       # Initial value for K, offset 0.2 + exp(log(K)) is the forced separation between the two components
    
  Par0 = list(mu = mu0.step,
              sigma = sigma0.step,
              beta1 = c(beta01),
              beta2temp =  beta01[-1],
              lbetaThr= lbetaThr0,
              logK = log(K))
  
    # PROGRESSIVE SHARPNESS PROCEDURE 
  
    # Prepare for the progressive sharpness procedure to find a better initial value for betaThr

    # In dat: parameters that are fixed during the procedure but provide useful information
    dat = list(
      logmu = log(Par0$mu),             # Fix the state-dependent parameters
      logsigma = log(c(Par0$sigma)),    # Fix the state-dependent parameters
      step = data$step,                 # step-length data
      N = j,                            # Number of states
      X1 = X1,                          # Matrice of covariate associated with Gamma1
      X2 = X2,                          # Matrice of covariate associated with Gamma2
      U =  as.matrix(U, ncol = 1)       # Covariate associated with potential disturbance: in the step function
    )
    
    # In par: parameters that are estimated during the procedure
    par = list(
      beta1 = Par0$beta1,               # initial value for TPM parameters on logit scale  
      beta2temp = Par0$beta2temp,       # initial value for TPM parameters on logit scale  
      delta1 = c(0,0),                  # initial state probability 1/N for each state for component 1
      delta2 = c(0,0),                  # initial state probability 1/N for each state for component 2
      lbetaThr = Par0$lbetaThr,         # initial value for logarithm of threshold parameter
      logK = Par0$logK                  # initial value for log(K)
    )   
    
    # Get a better initial value for betaThr, based on the fixed state-dependent parameters
    result <- getbetaThr(par,dat,betaThr0) 
    
    
    
  # FIT LASSO THMM
    
    if(!is.null(result)){
    # Minimise lasso likelihood
      dat = list(
        step = data$step,   
        N = j,
        X1 = X1,
        X2 = X2,
        U =  as.matrix(U, ncol = 1),
        lambda_lasso = 10^(-1 + (6 * (i-1) / 49)) # approximation of what would be a grid search for lambda on log-scale. Divided by 49 because we explore 50 random initial values.
      )
      
      
      par = list(
        logmu = log(result$model$mu),     
        logsigma = log(c(result$model$sigma)),
        beta1 = result$model$beta1,         
        beta2temp =  result$model$beta2temp,
        delta1 = c(0,0),
        delta2 = c(0,0),
        lbetaThr = result$model$lbetaThr,
        logK = result$model$logK
      ) 
      
    tryCatch({
      result_lasso <- FitLassoTHMM(par,dat)   # Fit Lasso-THMM
      
    }, error = function(e) {
      skip_to_next <<- TRUE
    })

      # EXTRACT RESULTS
      
    if(!skip_to_next){return(list(model = result_lasso$model, npllk =  result_lasso$npllk, code = result_lasso$code,betaThr0 = betaThr0)        
    )}else{
      return(list(model = NULL, nllk = NA, code = NA , modH0 = NA,betaThr0=NA))
    }
  }else {
    return(list(model = NULL, nllk = NA, code = NA , modH0 = NA,betaThr0=NA))
  }
}

# Run the parallelized code
results = mclapply(1:Nrep, fit_model, mc.cores = 1)



##############################
###### Extract results #######
##############################

mod = lapply(results, function(res) res$model)
beta0 = lapply(results, function(res) res$betaThr0)
npllk = code = lbeta = array()
for(i in (1:Nrep)){
  if(is.null(mod[[i]])){ 
    npllk[i]= NA
    mod[[i]] = NA
    code[i] = NA
  }else{
    code[i] = results[[i]]$code
    if(is.na(code[i]) || code[i]!=0){ # Only keep models that had numerical convergence
      npllk[i] = NA
      lbeta[i] = NA
    }else{
      npllk[i] = results[[i]]$npllk
      lbeta[i] = results[[i]]$model$lbetaThr
    }
  }
}


idx = which.min(npllk)

m = list(model = mod[[idx]], nllk = npllk[idx],code = code[idx],conv=length(which(!is.na(npllk))))

# saveRDS(m, file=paste0("Results_2d_3s_",n, "_LASSO_pen_scenario3_",id, ".RData"))
# saveRDS(m, file=paste0("Results_2d_3s_",n, "_LASSO_pen_null_scenario3_",id, ".RData")), if null model is simulated



# References

# Pohle, J., Langrock, R., Van Beest, F. M., & Schmidt, N. M. (2017). Selecting the number of states in hidden Markov models: pragmatic solutions illustrated using animal movement. Journal of Agricultural, Biological and Environmental Statistics, 22(3), 270-293.
