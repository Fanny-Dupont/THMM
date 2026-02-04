###############################################################################
#                    Code to Simulate Data and Fit THMM for the paper:        #
#      Estimating the distance at which narwhals respond to disturbance:      #
#               a penalized threshold hidden Markov model                     #
#                                                                             #
#      Authors: Fanny Dupont, Marianne Marcoux, Nigel Hussey,                 # 
#                       Jackie Dawson, Marie Auger-Méthé                      #
###############################################################################
# Details and comments for all the functions are given in SourceFunctions.R
# only difference is in l.81, because we adapt the constraints in the TPM with 
# the fact that there are covariates in TPMs in the case study.
# Takes approx 4h to run with 50 replications, 8h with 100. 50 leads to a similar estimated threshold of 3.41 km.


library(momentuHMM)
library(boot)
library(tidyverse)
library(pracma)
library(parallel)
library(dirmult)
library(LaMa)
library(dplyr)
library(terra)
library(sf)
source("SourceFunctions.R")

# Load data 
data = read.csv("data_casestudy.csv")

# Distance above which we know narwhal can't be disturbed: 77km, which is reasonable based on previous studies about narwhals hearing/detection range.
indexQ = which(data$distance_km>=quantile(data$distance_km,0.6,na.rm=TRUE))


# Create matrix for both tpm's
X1 = as.matrix(cbind(Intercept=rep(1,length(data$dtoshore)),data$dtoshore))
X2 = as.matrix(cbind(Intercept=rep(1,length(data$dtoshore)),data$dtoshore))


# Create covariate u with potential disturbance effects
U = cbind(data$exposure*(1-data$noland),data$exposure*data$noland)
U = ifelse(is.na(U),0,U) # NAs means no vessel so no exposure


# Standardise the covariate
U <- apply(U, 2, function(x) {
  rng <- range(x, na.rm = TRUE)
  (x - rng[1]) / (rng[2] - rng[1])
})



getbetaThr_covinTPM <- function(par0,dat0){ 
  ii <- 1 
  nllk <- c(1e6, 0) 
  b <- 1
  par <- par0
  dat <- dat0
  while (abs(nllk[ii + 1] - nllk[ii]) > 1e-6) { 
    if (ii > 100){ return(list(NULL)) }
    if (b == 1e3){ break}
    if (b > 1e2) { b  = 1e3 }  
    fit_result =   tryCatch({
      nll_move = function(par) {
        getAll(par, dat) 
        betaThr = exp(lbetaThr) 
        eta <- U %*% betaThr - 1
        
        # Define mixture probability nu: force the first quantile of the data in baseline 
        nu <- c(1 / (1 + exp(-b * eta)))
        nu[indexQ] <- 0
        
        delta1 = c(1, exp(delta1))
        delta1 = delta1 / sum(delta1)
        delta2 = c(1, exp(delta2))
        delta2 = delta2 / sum(delta2)
        delta_comp = matrix(delta1 * (1-nu[1]) + delta2 * nu[1], ncol = N)
        
        # Force a minimum different between the two components. 
        K = 0.2+exp(logK)
        
        #### ONLY DIFFERENCE WITH CODE IN SourceFunctions.R, SINCE THERE ARE COVARIATES IN TPM 
        # Compute the difference between the sum of the first row of beta1 and beta2temp components,
        # then add the offset K. This ensures beta2[1] will be larger than beta1 by at least K.
        el = sum(exp(beta1[1,(1:(N-1))]))-sum(exp(beta2temp[(1:(N-2))]))+K
        
        # Use max(0, el) to ensure el is non-negative and numerically stable,
        el = max0_smooth(el)+1e-5  # avoid log(0) with 1e-5
        beta21 = log(el)
        
        # Combine
        beta2 = matrix(rbind(c(beta21,beta2temp),beta1[2,]),nrow=nrow(beta1))
        
        nu <- array(rep(nu, each = N * N), dim = c(N, N, length(nu)))
        
        # Calculate Gamma1
        if (ncol(X1) == 1) { 
          Gamma1 <- tpm(beta1, byrow = TRUE) # multinomial logit link
          Gamma1 <- array(Gamma1, dim = c(N, N, length(data$ID))) # Gamma2 is N by N by T (here n) matrix
          
        } else { 
          Gamma1 <- tpm_g(X1, t(beta1), byrow = TRUE) # multinomial logit link, effect of covariate X1 in Gamma1
        }
        
        # Calculate Gamma2
        if (ncol(X2) == 1) { 
          Gamma2 <- tpm(beta2, byrow = TRUE) # multinomial logit link
          Gamma2 <- array(Gamma2, dim = c(N, N, length(data$ID))) # Gamma2 is N by N by T matrix
          
        } else {
          Gamma2 <- tpm_g(X2, t(beta2), byrow = TRUE) # multinomial logit link, effect of covariate X1 in Gamma1
        }
        
        Gamma_comp <- (1-nu) * Gamma1 + nu * Gamma2
       
        # exponentiating because all parameters strictly positive
        mu = exp(logmu)
        sigma = exp(logsigma)
        kappa = exp(logkappa)
        mu.maxdepth = exp(logmu.maxdepth)
        sigma.maxdepth = exp(logsigma.maxdepth)
        
        # reporting statements for later use
        REPORT(mu); ADREPORT(mu)
        REPORT(sigma); ADREPORT(sigma)
        REPORT(mu.maxdepth); ADREPORT(mu.maxdepth)
        REPORT(sigma.maxdepth); ADREPORT(sigma.maxdepth)
        REPORT(kappa); ADREPORT(kappa)
        REPORT(lbetaThr); ADREPORT(lbetaThr)
        REPORT(beta1); ADREPORT(beta1)
        REPORT(beta2); ADREPORT(beta2)
        REPORT(beta2temp); ADREPORT(beta2temp)
        REPORT(logK); ADREPORT(logK)
        
        # calculating all state-dependent densities
        allprobs = matrix(1, nrow = length(step), ncol = N)
        
        
        # Diagonal matrix of log emission distributions for all observations (for all individuals).
        allprobs[1,] <- apply(cbind(mu,sigma),1,FUN=function(x)dgamma2(step[1],mean=x[1],sd=x[2]))*apply(cbind(mu.maxdepth,sigma.maxdepth),1,FUN=function(x)dgamma2(max_dep[1],mean=x[1],sd=x[2]))
        
        for(t in (2:length(step))){
          if(!is.na(step[t])){allprobs[t,] <- apply(cbind(mu,sigma),1,FUN=function(x)dgamma2(step[t],mean=x[1],sd=x[2]))}
          if(!is.na(angle[t])){allprobs[t,] <- allprobs[t,]*sapply(kappa,FUN=function(x)LaMa:::dvm(angle[t],mu=0,kappa=x))}
          if(!is.na(max_dep[t])){allprobs[t,] <- allprobs[t,]*apply(cbind(mu.maxdepth,sigma.maxdepth),1,FUN=function(x)dgamma2(max_dep[t],mean=x[1],sd=x[2]))}
        }
        
        delta_comp = matrix(rep(delta_comp, each = length(unique(data$ID))), ncol = length(delta_comp))
        
        -forward_g(delta = delta_comp, Gamma = Gamma_comp, allprobs=allprobs,trackID=data$ID) # simple forward algorithm
      }
      ii <- ii + 1
      
      if(b < 3){iteri <- 30}else{iteri = 300}
      obj <- MakeADFun(nll_move, par, silent = TRUE) # Create objective function
      opt <- nlminb(obj$par, obj$fn, obj$gr) # Optimization
      result = obj$report()
      nllk[ii+1] = opt$objective
      list(model = result, nllk =  opt$objective, code = opt$convergence)        
    }, error = function(e) {
      NULL
    })
    if(is.null(fit_result)){break}
    
    result <- fit_result
    par <- list(
      lbetaThr = result$model$lbetaThr,
      delta1 = c(0, 0),
      delta2 = c(0, 0),
      beta1 = result$model$beta1,
      beta2temp = result$model$beta2temp,
      logK = result$model$logK
    )
    dat = list(
      logmu = dat$logmu,                           # initial mean of step length
      logsigma = dat$logsigma,                     # initial standard deviation of step length
      logkappa = dat$logkappa,                     # initial concentration of turning angle
      logmu.maxdepth = dat$logmu.maxdepth,         # initial mean of maximum depth
      logsigma.maxdepth = dat$logsigma.maxdepth,   # initial standard deviation of maximum depth
      step = dat$step,                             # data of step length
      angle = dat$angle,                           # data of turning angle
      max_dep = dat$max_dep,                       # data of maximum depth
      N = dat$N,                                   # Number of states
      X1 = dat$X1,                                 # Vector of covariates in Gamma1 
      X2 = dat$X2,                                 # Vector of covariates in Gamma2
      U =  dat$U,                                  # Vector of covariates in the step function
      indexQ = dat$indexQ
    )
    
    # Increase shaprness parameter
    b <- b + ii/2
  }
  return(fit_result)
}

FitLassoTHMM_covinTPM <-  function(par0,dat0){ # Function to fit the Lasso THMM
  ii <- 1               
  nllk <- c(1e6, 0)     # Initialize neg-loglikelihood
  par <- par0           # Initial parameter values to be estimated during the iterative procedure.
  dat <- dat0           # Initial parameter values that are not estimated during the iterative procedure.
  
  # Identify indices corresponding to 'lbetaThr'
  p <- length(which(grepl("lbetaThr", names(unlist(par)))))
  
  while (abs(nllk[ii + 1] - nllk[ii]) > 1e-6) {
    if (ii %% 10 == 0) {
      cat("Iter", ii, ": nllk =", nllk[ii], "\n")
    }
    
    if (ii > 150){ return(list(model = result, nllk =  opt$objective, code = 2)) }
      fit_result <- tryCatch({
      nll_move = function(par) {
        getAll(par, dat)            # makes everything contained available without $
        b = 5e2                     # Sharpness parameter
        betaThr = exp(lbetaThr)     # Positive threshold paramter
        
        
        eta <- U %*% betaThr - 1
        
        # Define mixture probability nu: force the first quantile of the data in baseline 
        nu <- c(1 / (1 + exp(-b * eta)))
        nu[indexQ] <- 0
        n <- length(nu)

        delta1 = c(1, exp(delta1))
        delta1 = delta1 / sum(delta1)
        delta2 = c(1, exp(delta2))
        delta2 = delta2 / sum(delta2)
        delta_comp = matrix(delta1 * (1-nu[1]) + delta2 * nu[1], ncol = N)
        
        
        # Force a minimum different between the two components (only difference with SourceFunction.R). 
        K = 0.2+exp(logK) # Positive offset K to ensure separation; logK is a parameter to estimate
        
        el = sum(exp(beta1[1,(1:(N-1))]))-sum(exp(beta2temp[(1:(N-2))]))+K
        el = max0_smooth(el)+1e-5
        beta21 = log(el)
        beta2 = matrix(rbind(c(beta21,beta2temp),beta1[2,]),nrow=nrow(beta1))
        
        # Change dimensions of nu to make it appropriate for later use: NxNxT matrix
        nu <- array(rep(nu, each = N * N), dim = c(N, N, length(nu)))
        
        # Calculate t.p.m: Gamma1
        if (ncol(X1) == 1) { 
          Gamma1 <- tpm(beta1, byrow = TRUE)          # multinomial logit link
          Gamma1 <- array(Gamma1, dim = c(N, N, n))
          
        } else { # If covariate effects in t.p.m
          Gamma1 <- tpm_g(X1, t(beta1), byrow = TRUE) # multinomial logit link with covariate effects
        }
        
        # Calculate t.p.m: Gamma2
        if (ncol(X2) == 1) { 
          Gamma2 <- tpm(beta2, byrow = TRUE)          # multinomial logit link
          Gamma2 <- array(Gamma2, dim = c(N, N, n))
          
        } else {# If covariate effects in t.p.m
          Gamma2 <- tpm_g(X2, t(beta2), byrow = TRUE) # multinomial logit link with covariate effects
        }   
        # to change when multiple individuals
        Gamma_comp <- (1-nu) * Gamma1 + nu * Gamma2
        
        # exponentiating because all parameters strictly positive
        mu = exp(logmu)
        sigma = exp(logsigma)
        kappa = exp(logkappa)
        mu.maxdepth = exp(logmu.maxdepth)
        sigma.maxdepth = exp(logsigma.maxdepth)
        
        # reporting statements for later use
        REPORT(mu); ADREPORT(mu)
        REPORT(sigma); ADREPORT(sigma)
        REPORT(mu.maxdepth); ADREPORT(mu.maxdepth)
        REPORT(sigma.maxdepth); ADREPORT(sigma.maxdepth)
        REPORT(kappa); ADREPORT(kappa)
        REPORT(betaThr); ADREPORT(betaThr)
        REPORT(lbetaThr); ADREPORT(lbetaThr)
        REPORT(beta1); ADREPORT(beta1)
        REPORT(beta2); ADREPORT(beta2)
        REPORT(delta1); ADREPORT(delta1)
        REPORT(delta2); ADREPORT(delta2)
        REPORT(beta2temp); ADREPORT(beta2temp)
        REPORT(lambda_lasso); ADREPORT(lambda_lasso)
        REPORT(logK); ADREPORT(logK)
        
        # calculating all state-dependent densities
        allprobs = matrix(1, nrow = length(step), ncol = N)
        
        
        # Diagonal matrix of log emission distributions for all observations (for all individuals).
        allprobs[1,] <- apply(cbind(mu,sigma),1,FUN=function(x)dgamma2(step[1],mean=x[1],sd=x[2]))*apply(cbind(mu.maxdepth,sigma.maxdepth),1,FUN=function(x)dgamma2(max_dep[1],mean=x[1],sd=x[2]))
        
        for(t in (2:length(step))){
          if(!is.na(step[t])){allprobs[t,] <- apply(cbind(mu,sigma),1,FUN=function(x)dgamma2(step[t],mean=x[1],sd=x[2]))}
          if(!is.na(angle[t])){allprobs[t,] <- allprobs[t,]*sapply(kappa,FUN=function(x)LaMa:::dvm(angle[t],mu=0,kappa=x))}
          if(!is.na(max_dep[t])){allprobs[t,] <- allprobs[t,]*apply(cbind(mu.maxdepth,sigma.maxdepth),1,FUN=function(x)dgamma2(max_dep[t],mean=x[1],sd=x[2]))}
        }
        
        delta_comp = matrix(rep(delta_comp, each = length(unique(data$ID))), ncol = length(delta_comp))
        
        -forward_g(delta_comp, Gamma_comp, allprobs,data$ID) + lambda_lasso*sum(betaThr) # simple forward algorithm w Lasso penalty
        
      }
      ii <- ii + 1
      obj <- MakeADFun(nll_move, par, silent = TRUE) # Create objective function
      opt <- nlminb(obj$par, obj$fn, obj$gr) # Optimization
      result = obj$report()
      nllk[ii+1] = opt$objective
      
           list(model = result, nllk =  opt$objective, code = opt$convergence,He = obj$he()) 
    }, error = function(e) {
      NULL # Return NULL if an error occurs
    })
    
    if (is.null(fit_result)) {
      break
    }
    
    result <- fit_result
    lambda <- result$model$lambda_lasso
    betaThr<- exp(result$model$lbetaThr)
    
    par <- list(
      logmu = log(result$model$mu),     
      logsigma = log(c(result$model$sigma)),
      logkappa = log(c(result$model$kappa)), 
      logmu.maxdepth = log(c(result$model$mu.maxdepth)),
      logsigma.maxdepth = log(c(result$model$sigma.maxdepth)),
      lbetaThr = result$model$lbetaThr,
      delta1 = c(0, 0),
      delta2 = c(0, 0),
      beta1 = result$model$beta1,
      beta2temp = result$model$beta2temp,
      logK = result$model$logK
    )
    dat = list(
      step = dat$step,  
      angle = dat$angle,
      max_dep = dat$max_dep,
      N = dat$N,
      X1 = dat$X1,
      X2 = dat$X2,
      U =  dat$U,
      lambda_lasso = p/(sum(betaThr)),
      indexQ = dat$indexQ
    )
  }
  return(fit_result)
}



#### Fit Model ####

Nrep = 50       # Number of random initial values explored
bb = list()     # Initial list for Nrep models
nllk = array()  # Initial array for Nrep negative log-likelihoods
set.seed(1,kind="Mersenne-Twister",normal.kind = "Inversion")
# Define the function to be parallelized
fit_model <- function(i) {
  # We start by fitting a standard HMM to the data with momentuHMM to get estimates of state-dependent parameters and regression coefficients for TPMs.
  skip_to_next <- FALSE
  j <- 3
  beta01 <- beta02 <- array(NA, dim = c(j, j))
  
  # Generate beta matrices
  for (l in 1:j) {
    alpha <- rep(0.2 / (j - 1), j)
    alpha[l] <- 0.6
    alpha <- alpha * 10
    beta01[l, ] <- rdirichlet(1, alpha)
    beta02[l, ] <- rdirichlet(1, alpha)
  }
  
  # Process beta matrices to get initial parameters for the TPM on the multinomial logit scale
  beta01 <- logit(beta01)
  beta01 <- t(beta01)
  beta01 <- matrix(beta01[col(beta01) != row(beta01)], ncol = j * (j - 1), nrow = 1)
  beta01 <- rbind(beta01, do.call(rbind, replicate((ncol(X1) - 1), rep(0, j * (j - 1)), simplify = FALSE)))
  
  beta02 <- logit(beta02)
  beta02 <- t(beta02)
  beta02 <- matrix(beta02[col(beta02) != row(beta02)], ncol = (j * (j - 1)), nrow = 1)
  beta02 <- rbind(beta02, do.call(rbind, replicate((ncol(X2) - 1), rep(0, j * (j - 1)), simplify = FALSE)))
  
  # Define initial state-dependent parameters 
  mu0.step <-   c(rnorm(1,1,0.5), rnorm(1,2.5,0.5), rnorm(1,1,0.5))  # mean step-length, one per state
  shape0.step <- c(rgamma(j, shape = 1, scale = 2.5))                # shape (gamma distribution) step-length, one per state
  sigma0.step <- mu0.step / sqrt(shape0.step)                        # convert shape into standard deviation for the gamma distribution in momentuHMM
  
  mu0.max_dep <- abs(c(rnorm(1, mean = 30, sd = 0.5),rnorm(1, mean = 50, sd = 5), rnorm(1, mean = 350, sd = 50)))       # mean max depth (gamma distribution), one per state
  sigma0.max_dep <- abs(c(rnorm(1, mean = 30, sd = 0.5),rnorm(1, mean = 50, sd = 5), rnorm(1, mean = 150, sd = 30)))    # sd max depth (gamma distribution), one per state
   
  kappa0 <- c(rnorm(1,1,0.5), rnorm(1,5,0.5), rnorm(1,1,0.5))        # concentration parameter for von-mises, one per state
  
  # Put all initial parameters in Par0
  Par0 <- list(
    step = c(mu0.step, sigma0.step),
    angle = ifelse(abs(kappa0)<3.14,kappa0,rnorm(1,0.5,0.1)),
    max_dep = c(mu0.max_dep, sigma0.max_dep)
  )
  
  # Fit momentuHMM model with tryCatch to get initial paramters that fit the data well, with effect of covariate "dtoshore" in the TPM
  exData <- momentuHMM:::prepData(data %>% dplyr::select(step, angle, max_dep, ID, dtoshore),covNames = c("dtoshore"),coordNames = NULL)
  dist <- list(step = "gamma", angle = "vm", max_dep = "gamma")
  
  tryCatch(b <- momentuHMM:::fitHMM(exData, nbState = j,formula= ~dtoshore, dist = dist, Par0 = Par0, beta0 = beta01),error = function(e){skip_to_next<<- TRUE})
  
  if(!skip_to_next) {
    # As in the simulation study
    H0 = b
    
    # Initial state probabilities and threshold parmaeters (beta0)
    delta0 <- list(rep(1 / j, j), rep(1 / j, j))
    lbetaThr0= c(rnorm(1,-10,1),rnorm(1,0,2))  
    betaThr0 = exp(lbetaThr0)
    K = exp(-15) # very close to 0 but not 0 bc logarithm
    
    # Use the values obtained with momentuHMM as initial parameters for the THMM
    Par0 = list(mu = c(H0$mle$step[1,]),
                sigma = c(H0$mle$step[2,]),
                mu.maxdepth = c(as.numeric(H0$mle$max_dep[1,])),
                sigma.maxdepth = c(as.numeric(H0$mle$max_dep[2,])),
                kappa =  c(as.numeric(H0$mle$angle[2,])),
                lbetaThr= lbetaThr0,
                beta1 = H0$mle$beta,
                beta2temp = beta01[1,-1])
    Par0 = lapply(Par0,
                  FUN = function(x) ifelse(x > 1e3, abs(rnorm(1, 10, 10)), x)) # Control any outlier
    
    
    dat = list(
      logmu = log(Par0$mu),     
      logsigma = log(c(Par0$sigma)),
      logmu.maxdepth = log(c(as.numeric(Par0$mu.maxdepth))),
      logsigma.maxdepth = log(c(as.numeric(Par0$sigma.maxdepth))),
      logkappa = log(c(Par0$kappa)), 
      step = data$step,  
      angle = data$angle,
      max_dep = data$max_dep,
      N = j,
      X1 = X1,
      X2 = X2,
      U =  U,
      indexQ=indexQ
    )
    
    par = list(
      beta1 = Par0$beta1,        
      beta2temp = Par0$beta2temp,        
      lbetaThr = Par0$lbetaThr,
      delta1 = c(0,0),
      delta2 = c(0,0),
      logK = log(K)
    )   
    
  message("Getting initial value for betaThr..." )
  result <- getbetaThr_covinTPM(par,dat)
  if(!is.null(result)){
    # Minimise lasso likelihood
    dat0 = list(
      step = data$step, 
      angle = data$angle, 
      max_dep = data$max_dep,
      N = j,
      X1 = X1,
      X2 = X2,
      U =  U,
      lambda_lasso = abs(rt(1,0.8))
    )
    par0 = list(
      logmu = log(result$model$mu),     
      logsigma = log(c(result$model$sigma)),
      logkappa = log(c(result$model$kappa)), 
      logmu.maxdepth = log(c(result$model$mu.maxdepth)),
      logsigma.maxdepth = log(c(result$model$sigma.maxdepth)),
      beta1 = result$model$beta1,         
      beta2temp =  result$model$beta2temp,
      delta1 = c(0,0),
      delta2 = c(0,0),
      lbetaThr = result$model$lbetaThr,
      logK = result$model$logK
    ) 
    
    tryCatch({
      message("Fitting model ..." )
      result_lasso <- FitLassoTHMM_covinTPM(par0,dat0)
      
    }, error = function(e) {
      skip_to_next <<- TRUE
    })
    
    
    if(!skip_to_next){return(list(model = result_lasso$model, nllk =  result_lasso$nllk,He = result_lasso$He, code = result_lasso$code,betaThr0 = betaThr0)        
    )}else{
      return(list(model = NULL,He = NA , nllk = NA, code = NA , modH0 = NA,betaThr0=NA))
    }
  }else {
    return(list(model = NULL,He = NA, nllk = NA, code = NA , modH0 = NA,betaThr0=NA))
  }
  }else {
    return(list(model = NULL,He = NA, nllk = NA, code = NA , modH0 = NA,betaThr0=NA))
  }
}





# Run the parallelized code
results = mclapply(1:Nrep, fit_model, mc.cores = 4)


results <- lapply(results, function(res) {
  if (!"model" %in% names(res)) {
    res$model <- NA
    res$code <- NA
    res$He <- NA
    res$betaThr0 <- NA
  }
   res
})



bb = lapply(results, function(res) res$model)
He = lapply(results, function(res) res$He)
bb0 = lapply(results, function(res) res$betaThr0)
nllk = code  = array()

for(i in (1:Nrep)){
  if(is.null(bb[[i]])){
    nllk[i]= NA
    bb[[i]] = NA
    code[i] = NA
  }else{
    code[i] = results[[i]]$code
    if(is.na(code[i]) || code[i]!=0){ # only keep models that converged
      nllk[i] = NA
    }else{
      nllk[i] = results[[i]]$nllk
    }
  }
}

print(nllk)
print(code)
idx = which.min(nllk)


m = list(model = bb[[idx]], nllk = nllk[idx],He=He[[idx]],code = code[idx],betaThr0=bb0[idx])



# Get threshold from model fitted

U = cbind(data$exposure*(1-data$noland),data$exposure*data$noland)
U = ifelse(is.na(U),0,U) # NAs means no vessel so no exposure

print("exposure w land")
((1/exp(m$model$lbetaThr[1]))*(max(U[,1])-min(U[,1]))+min(U[,1]))^{-1} 

print("exposure w/o land")
((1/exp(m$model$lbetaThr[2]))*(max(U[,2])-min(U[,2]))+min(U[,2]))^{-1}  




