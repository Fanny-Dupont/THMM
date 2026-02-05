###############################################################################
#               Functions used for the simulation study                       #
#                     and THMM fitting for the paper:                         #
#       Estimating the distance at which narwhal respond to disturbance:      #
#               a penalized threshold hidden Markov model                     #
#                                                                             #
#                                                                             #
#     Authors: Fanny Dupont, Marianne Marcoux, Nigel E. Hussey,               #
#              Jackie Dawson and Marie Auger-Méthé                            #
###############################################################################



#' Simulates gamma HMM
#' 
#' @param TT Integer, length of the time-series to simulate
#' @param N Integer, number of states of the HMM
#' @param Gamma Numerical matrice of transition probability and dimension TTxNxN. For stationary HMMs, use replicate to obtain TTxNxN Gamma with the same value at each time.
#' @param delta Vector of initial state probabilities (sum to 1), of dimension N
#' @param mu.gamma Vector of (positive) mean of gamma state-dependent distribution, of dimension N
#' @param sd.gamma Vector of (positive) sd of gamma state-dependent distribution, of dimension N

#' @return  A data.frame with 3 columns and TT rows with the following components:
#'          Index: Index in the time-series 
#'          State: True underlying state
#'          Step:  Step-length observation
sim_HMM_gamma <- function(TT, N, Gamma, delta, mu.gamma, sd.gamma) {
  # initialize output
  output <- data.frame(
    Index = 1:TT,
    State = NA,
    step = NA            # "step" is the name of the observation that follows a gamma-HMM 
  )
  
  # possible states
  states <- 1:N
  
  # initial state
  output$State[1] <- sample(x = states, size = 1, prob = delta)
  
  # initial gamma step
  output$step[1] <- rgamma2(
    n = 1, 
    mean = mu.gamma[output$State[1]], 
    sd = sd.gamma[output$State[1]]
  )
  
  # simulate HMM sequence
  for (i in 2:TT) {
    transition <- Gamma[i,output$State[i-1],]
    output$State[i] <- sample(x = states, size = 1, prob = transition)
    output$step[i] <- rgamma2(
      n = 1, 
      mean = mu.gamma[output$State[i]], 
      sd = sd.gamma[output$State[i]]
    )
  }
  
  output$State <- as.factor(output$State)
  return(output)
}


#' Simulates movement from HMM with gamma step-length and vmises turning angles state-dependent distributions
#' 
#' @param TT Integer, length of the time-series to simulate
#' @param N Integer, number of states of the HMM
#' @param Gamma Numerical matrice of transition probability and dimension TTxNxN. For stationary HMMs, use replicate to obtain TTxNxN Gamma with the same value at each time.
#' @param delta Vector of initial state probabilities (sum to 1), of dimension N
#' @param mu.gamma Vector of (positive) mean of gamma state-dependent distribution, of dimension N
#' @param sd.gamma Vector of (positive) sd of gamma state-dependent distribution, of dimension N
#' @param kappa Vector of (positive) concentration of von Mises state-dependent distribution, of dimension N

#' @return  A data.frame with 6 columns and TT rows with the following components:
#'          Index: Index in the time-series 
#'          State: True underlying state
#'          Step:  Step-length observation
#'          Angle: Turning angle of observation
#'          x: x-axis location
#'          y: y-axis location

sim_HMM_move <- function(TT, N, Gamma, delta, mu.gamma, sd.gamma, kappa) {
  # initialize output
  output = data.frame("Index" = c(1:TT),
                      "ID" = rep(1,TT),
                      "State" = NA,
                      "step" = NA,
                      "angle" = NA,
                      "x" = NA,
                      "y" = NA)
  states = c(1:N)
  # get inital state from initial state distribution delta
  output$State[1] = sample(x=states, size=1, prob=delta)
  output$angle[1] = 0
  # get initial observation using initial state found and the state-dependent gamma distributions (tangle is NA)
  output$step[1] = rgamma2(n=1, mean=mu.gamma[output$State[1]], sd=sd.gamma[output$State[1]])
  output$x[1] = output$y[1] = 0
  
  for (i in 2:TT) {
    # using the previous state and the transition probability, get the following state and observation
    transition = Gamma[i,output$State[i-1],]
    output$State[i] = sample(x=states, size=1, prob=transition)
    output$step[i] = rgamma2(n=1, mean=mu.gamma[output$State[i]], sd=sd.gamma[output$State[i]])
    output$angle[i] = output$angle[i-1] + LaMa:::rvm(n=1, 0, kappa=kappa[output$State[i]])
    output$angle[i] = atan2(sin(output$angle[i]), cos(output$angle[i]))
    output$x[i] = output$x[i-1] +  output$step[i]*cos(output$angle[i]) 
    output$y[i] = output$y[i-1] +  output$step[i]*sin(output$angle[i]) 
    
  }
  
  output$State <- as.factor(output$State)
  return(output)
}


#' Simulate a Bernoulli sequence that is more likely to stay the same as previous than not (hence "sticky").
#'
#' This function generates a binary time series (`0`/`1`) where each element
#' tends to repeat the previous value with a given probability (p_same), with 
#' temporal dependence between consecutive values.
#'
#' @param n Integer, length of the sequence to simulate.
#' @param p_same Numeric (between 0 and 1), probability that the next value
#'   is the same as the previous one (default = 0.8).
#' @param x0 Integer (0 or 1), initial value of the sequence (default = 0).

simulate_sticky_bernoulli <- function(n, p_same = 0.8, x0 = 0) {
  x <- numeric(n)
  x[1] <- x0
  for (i in 2:n) {
    if (runif(1) < p_same) {
      x[i] <- x[i - 1]      # stay the same
    } else {
      x[i] <- 1 - x[i - 1]  # switch
    }
  }
  return(x)
}




#' Smooth approximation of the maximum function: max(x + eps, 0)
#' 
#' @param x Numeric vector or matrix
#' @param alpha numeric scalar. It controls the smoothness of the approximation. Larger values of `alpha` make the function closer to `max(0, x)`, while smaller values make it smoother.
#' @param eps Numeric scalar, helps with numerical stability.
#' @return  A numeric vector of the same length as `x`.

max0_smooth <- function(x, alpha = 10, eps = 0){
  1 / alpha * log(1 + exp(alpha * (x + eps)))
}






#' Progress sharpness procedure to estimate an appropriate threshold parameter
#'  
#' @param par0 Numeric vector.
#'   Initial parameter values to be estimated during the iterative procedure.
#'   Typically includes parameters associated with the hidden process.
#'
#' @param dat0 List or vector.
#'   Contains parameters and data that are not estimated during the procedure,
#'   such as state-dependent parameters, number of states, data streams of
#'   interest, and vectors of covariates used in the step and transition
#'   probability models.
#'
#' @param betaThr0 Numeric scalar.
#'   Initial threshold parameter (i.e., beta_0) for the sharpness initialization.
#'
#' @return A list with the following components:
#'   model: The final fitted model.
#'   nllk: The final negative log-likelihood value.
#'   code: Optimization convergence code.
#'   betaThr0: The initial threshold parameter value.


getbetaThr <- function(par0,dat0,betaThr0){
  ii <- 1
  nllk <- c(1e6, 0) # Initialize neg-loglikelihood
  b <- 1            # Initial sharpness parameter
  par <- par0       # Initial parameter values to be estimated during the iterative procedure.
  dat <- dat0       # Initial parameter values that are not estimated during the iterative procedure.
  
  while (abs(nllk[ii + 1] - nllk[ii]) > 1e-6) { # Stopping rule
    if (b == 1e3){ break}                       # Stopping rule
   
        fit_result =   tryCatch({
        nll_move = function(par) {
        getAll(par, dat)                        # makes everything contained available without $
        
        betaThr = exp(lbetaThr)                 # threshold parameter (i.e., beta_0)
       
        # Define mixture probability nu: force the first 10% quantile of the data in baseline 
        nu <- sapply(U,FUN=function(x)ifelse(x < quantile(U,0.1),0,c(1 / (1 + exp(-b *  (x %*% betaThr - 1))))))
        
        n <- length(nu)
        delta1 = c(1, exp(delta1))
        delta1 = delta1 / sum(delta1)
        delta2 = c(1, exp(delta2))
        delta2 = delta2 / sum(delta2)
        delta_comp = matrix(delta1 * (1-nu[1]) + delta2 * nu[1], ncol = N) # Vector of inital state probability
        
        
        # Force a minimum different between the two components. 
        K = 0.2+exp(logK) # Positive offset K to ensure separation between normalizing constants of 0.2 + exp(logK) for state 1; logK is a parameter to estimate

        # Compute the difference between the sum of the first row of beta1 and beta2temp components,
        # then add the offset K. 
        el = sum(exp(beta1[(1:(N-1))]))-sum(exp(beta2temp[(1:(N-2))]))+K
        
        # Use max(0, el) to ensure el is non-negative and numerically stable,
        el = max0_smooth(el)
        beta21 = log(el+1e-5) # avoid log(0) with 1e-5
        
        # Combine
        beta2 = c(beta21,beta2temp)
        
        # Change dimensions of nu to make it appropriate for later use: NxNxT matrix
        nu <- array(rep(nu, each = N * N), dim = c(N, N, length(nu)))
        
        if (ncol(X1) == 1) {
          Gamma1 <- tpm(beta1, byrow = TRUE)        # multinomial logit link
          Gamma1 <- array(Gamma1, dim = c(N, N, n)) # Gamma1 is N x N x T (here n) matrix
          
        } else {
          Gamma1 <- tpm_g(X1[-1, ], beta1, byrow = TRUE) # multinomial logit link, effect of covariate X1 in Gamma1
        }
        
        # Calculate Gamma2
        if (ncol(X2) == 1) {
          Gamma2 <- tpm(beta2, byrow = TRUE)        # multinomial logit link
          Gamma2 <- array(Gamma2, dim = c(N, N, n)) # Gamma2 is N by N by T (here n) matrix
          
        } else {
          Gamma2 <- tpm_g(X2[-1, ], beta2, byrow = TRUE) # multinomial logit link, effect of covariate X2 in Gamma2
        }
        
        Gamma_comp <- (1-nu) * Gamma1 + nu * Gamma2
        
        # exponentiating because all parameters strictly positive
        mu = exp(logmu)
        sigma = exp(logsigma)

        # reporting statements for later use
        REPORT(mu); ADREPORT(mu)
        REPORT(sigma); ADREPORT(sigma)
        REPORT(lbetaThr); ADREPORT(lbetaThr)
        REPORT(beta1); ADREPORT(beta1)
        REPORT(beta2); ADREPORT(beta2)
        REPORT(beta2temp); ADREPORT(beta2temp)
        REPORT(logK); ADREPORT(logK)
        
        # calculating all state-dependent densities
        allprobs = matrix(1, nrow = length(step), ncol = N)
        ind = which(!is.na(step)) # only for non-NA obs.
        
        for (j in 1:N) {
          allprobs[1, j] = dgamma2(step[1], mu[j], sigma[j])
          allprobs[ind, j] = dgamma2(step[ind], mu[j], sigma[j])
        }
        
        -forward_g(delta_comp, Gamma_comp, allprobs) # simple forward algorithm for non-stationary HMM, to adapt when there are multiple individuals, see LaMa documentation.
        
        }
        
      # increment iteration step
      ii <- ii + 1
      
      obj <- MakeADFun(nll_move, par, silent = TRUE) # Create objective function
      opt <- nlminb(obj$par, obj$fn, obj$gr)         # Optimization
      result = obj$report()                          # Get the results
      nllk[ii+1] = opt$objective                     # Store negative log-likelihood at the minimum
      list(model = result, nllk =  opt$objective, code = opt$convergence,betaThr0 = betaThr0)        
    }, error = function(e) {
      NULL
    })
    if(is.null(fit_result)){break}     # If the model didn't work, stop the procedure 
    
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
      logmu = dat0$logmu,              # initial mean of step length
      logsigma = dat0$logsigma,        # initial standard deviation of step length
      step = dat0$step,                # step lengths
      N = dat0$N,                      # Number of states
      X1 = dat0$X1,                    # Vector of covariates in Gamma1 
      X2 = dat0$X2,                    # Vector of covariates in Gamma2
      U =  as.matrix(dat0$U, ncol = 1) # Vector of covariates in the step function
    )
    
    # Increase shaprness parameter
    b <- b + ii/2
  }
  
  return(fit_result)
}



#' Fit the THMM with lasso regularization
#'  
#' @param par0 Numeric vector.
#'   Initial parameter values to be estimated during the iterative procedure, should include all parameters associated with the THMM + log(K) and lambda_lasso (Lasso penalty parameter)
#'
#' @param dat0 List or vector.
#'   Contains parameters and data that are not estimated during the procedure,
#'   such as number of states, data streams of
#'   interest, and vectors of covariates used in the step and transition
#'   probability models.
#'
#' @param betaThr0 Numeric scalar.
#'   Initial threshold parameter for the sharpness initialization.
#'
#' @return A list with the following components:
#'   model: The final fitted model, with a list of the estimated parameters and TPM.
#'   npllk: The final negative penalized log-likelihood value.
#'   code:  The optimization convergence code, default optimizer is nlminb()
#'   betaThr0: The initial threshold parameter value. 

FitLassoTHMM <-  function(par0,dat0){
  ii <- 1
  npllk <- c(1e6, 0) # Initialize neg-loglikelihood
  par <- par0        # Initial parameter values to be estimated during the iterative procedure.
  dat <- dat0        # Initial parameter values that are not estimated during the iterative procedure.
  
  # Identify indices corresponding to 'lbetaThr'
  p <- length(which(grepl("lbetaThr", names(unlist(par)))))
  
  while (abs(npllk[ii + 1] - npllk[ii]) > 1e-6) {
    fit_result <- tryCatch({
      nll_move = function(par) {
        getAll(par, dat)            # makes everything contained available without $
        b = 5e2                     # Sharpness parameter
        betaThr = exp(lbetaThr)     # Positive threshold paramter
       
        # Define mixture probability nu: force the first quantile of the data in baseline 
        nu <- sapply(U,FUN=function(x)ifelse(x < quantile(U,0.1),0,c(1 / (1 + exp(-b *  (x %*% betaThr - 1))))))
        n <- length(nu)

        delta1 = c(1, exp(delta1))
        delta1 = delta1 / sum(delta1)
        delta2 = c(1, exp(delta2))
        delta2 = delta2 / sum(delta2)
        delta_comp = matrix(delta1 * (1-nu[1]) + delta2 * nu[1], ncol = N) # Vector of inital state probability
        
        
        # Force a minimum different between the two components. 
        K = 0.2+exp(logK) # Positive offset K to ensure separation between normalizing constants of 0.2 + exp(logK) for state 1; logK is a parameter to estimate
        
        # Compute the difference between the sum of the first row of beta1 and beta2temp components,
        # then add the offset K.
        el = sum(exp(beta1[(1:(N-1))]))-sum(exp(beta2temp[(1:(N-2))]))+K
        
        # Use max(0, el) to ensure el is non-negative and numerically stable,
        # then add a small constant to avoid log(0)
        el = max0_smooth(el)+1e-5
        beta21 = log(el)
        
        # Combine
        beta2 = c(beta21,beta2temp)
        
        # Change dimensions of nu to make it appropriate for later use: NxNxTT matrix
        nu <- array(rep(nu, each = N * N), dim = c(N, N, length(nu)))
        
        if (ncol(X1) == 1) {
          Gamma1 <- tpm(beta1, byrow = TRUE)        # multinomial logit link with intercept only
          Gamma1 <- array(Gamma1, dim = c(N, N, n)) # Gamma1 is N x N x T matrix
          
        } else {
          Gamma1 <- tpm_g(X1[-1, ], beta1, byrow = TRUE) # multinomial logit link, effect of covariate X1 in Gamma1
        }
        
        # Calculate Gamma2
        if (ncol(X2) == 1) {
          Gamma2 <- tpm(beta2, byrow = TRUE)         # multinomial logit link w intercept only
          Gamma2 <- array(Gamma2, dim = c(N, N, n))  # Gamma2 is N by N by T matrix
          
        } else {
          Gamma2 <- tpm_g(X2[-1, ], beta2, byrow = TRUE) # multinomial logit link, effect of covariate X2 in Gamma2
        }
        
        Gamma_comp <- (1-nu) * Gamma1 + nu * Gamma2     
        
        # exponentiating because all parameters of gamma state-dependent distributions are strictly positive
        mu = exp(logmu)
        sigma = exp(logsigma)

        # reporting statements for later use
        REPORT(mu); ADREPORT(mu)
        REPORT(sigma); ADREPORT(sigma)
        REPORT(lbetaThr); ADREPORT(lbetaThr)
        REPORT(beta1); ADREPORT(beta1)
        REPORT(beta2); ADREPORT(beta2)
        REPORT(beta2temp); ADREPORT(beta2temp)
        REPORT(lambda_lasso); ADREPORT(lambda_lasso) 
        REPORT(logK); ADREPORT(logK)
        
        # calculating all state-dependent densities
        allprobs = matrix(1, nrow = length(step), ncol = N)
        ind = which(!is.na(step)) # only for non-NA obs.
        
        for (j in 1:N) {
          allprobs[1, j] = dgamma2(step[1], mu[j], sigma[j])
          allprobs[ind, j] = dgamma2(step[ind], mu[j], sigma[j])
        }
        
        
        -forward_g(delta_comp, Gamma_comp, allprobs) + lambda_lasso*sum(betaThr) # nllk + lasso penalty
        
      }
      ii <- ii + 1
      obj <- MakeADFun(nll_move, par, silent = TRUE) # Create objective function
      opt <- nlminb(obj$par, obj$fn, obj$gr)         # Optimization
      result = obj$report()                          # Get the results
      # print(opt$message)
      npllk[ii+1] = opt$objective                    # Store negative penalized log-likelihood at the minimum
      
      list(model = result, npllk =  opt$objective, code = opt$convergence)       
    }, error = function(e) {
      # print(opt$message)
      NULL                                           # Return NULL if an error occurs
    })
    
    if (is.null(fit_result)) {
      break
    }
    if (ii > 100){
      break                                          # max number of iterations
    } 
    
    result <- fit_result
    lambda <- result$model$lambda_lasso
    betaThr<- exp(result$model$lbetaThr)
    
    # Update parameters
    par <- list(
      logmu = log(result$model$mu),     
      logsigma = log(c(result$model$sigma)),
      lbetaThr = result$model$lbetaThr,
      delta1 = c(0, 0),
      delta2 = c(0, 0),
      beta1 = result$model$beta1,
      beta2temp = result$model$beta2temp,
      logK = result$model$logK
    )
    # Update lambda_lasso, using the update rule in equation (9) of the main manuscript
    dat = list(
      step = dat0$step,   
      N = dat0$N,
      X1 = dat0$X1,
      X2 = dat0$X2,
      U =  as.matrix(dat0$U, ncol = 1),
      lambda_lasso = p/(sum(betaThr))
    )
  }
  return(fit_result)
}

