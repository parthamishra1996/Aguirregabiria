u <- function(i, s, q, params){
  #print("Entered u")
  params$p_r*s - params$p_w*q*(params$Q - i + s) - params$alpha*(i - s) - params$eta_q*q
}

calc_log_CCP <- function(i, s, params, V){
  #print("Calculating CCP")
  force(V)
  ans <- u(i, s, 0, params) - u(i, s, 1, params) + params$beta*(V(i-s) - V(params$Q) ) #
  
  -log(1 + exp(ans))
}


## Bellman operator
bellman_contraction <- function(V.fn, params) {
  force(V.fn)
  v <- 
    function(i){
      #print("i")
      #print(i)
      #print(length(i))
      integrand_lower <- function(s) dlnorm(s, meanlog = params$mu, sdlog = params$sigma) * (u(i, s, 1, params) - calc_log_CCP(i, s, params, V.fn ))
      #print("done with lower")
      integration_upper <- 1 - plnorm(i, meanlog = params$mu, sdlog = params$sigma) * (u(i, i, 1, params) - calc_log_CCP(i, i, params, V.fn ))
      #print("done with upper")
      c(
        integrate( #lower part
          integrand_lower, 
          lower = 0,
          upper = i
        )$value,
        integration_upper
      ) %>% 
        sum %>% {
          params$beta*V.fn(i) + . 
        }
      #print("done here")
    }
  #print(v)
  return(Vectorize(v))
}
  
## Value function iteration, V~
value_fn_iteration <- function(params) {
  epsilon <- 10^-4
  delta <- 1 + epsilon
  eval.points <- seq(0, 21, length.out = 100)
  V <- function(x) x
  while(delta > epsilon){
    V. <- V
    V <- bellman_contraction(V, params)
    V <- 
      ipol(
        val = V,
        intervals = c(0, 10),
        dims = 20
      )
    delta <- max(abs(V(eval.points) - V.(eval.points)))
    print(delta)
  }
  V
}

simulation <- function(V, params, nsamples = 10^4){
  force(V)
  
  i <- c(params$Q)
  s <- c()
  logccp <- c()
  demand <- rlnorm(nsamples, log(params$mu), log(params$sigma))
  
  for(x in seq(nsamples)){
    s = c(s, min( i[x], demand[x] ) )
    lccp = calc_log_CCP(i[x], s[x], params, V)
    logccp <- c(lccp)
    if(lccp < -2){ 
      i_new = params$Q
    } else{
      i_new = i[x] - s[x]
    }
    if(x < nsamples) i = c(i, i_new)
  }
  sim_data_new <- data.frame("i" = i, "s" = s, "logccp" = logccp)
  
  sim_data_new %>% saveRDS('../data/sim_data.rds')
}

estimate_mu_sigma<- function(sim_data){
  mu = 3; sigma = 2
  theta <- c(mu, sigma)
  theta <- mle(theta, sim_data)
  return(theta)
}

log_like <- function(theta, sim_data){
  ll <- sim_data %>% group_by(q) %>% mutate(
    log_s = log(s),
    phi = dlnorm(
      log_s, 
      meanlog = log(theta[1]), 
      sdlog = log(theta[2])
      ),
    phi_ = 1 - phi) %>% 
    select(phi, phi_) %>% 
    colSums %>% ungroup
  ll[phi][1] + ll[phi_][2]
}

mle <- function(theta, sim_data){
  obj <- optim(par = unlist(theta),
               fn = log_like,
               sim_data = sim_data,
               control = list(fnscale = -1), 
               method = 'L-BFGS-B',
               lower = c(0,0)
  )
  print(obj$par)
}

param_estimation <- function(sim_data){
  #sim_data <- data.frame("i" = runif(10^4, 0, 20), "s" = runif(10^4, 0, 20))
  #sim_data["s"] = sim_data$i
  
  sim_data['q'] = as.integer(sim_data['i'] == sim_data['s'])
  sim_data["i-s"] = sim_data['i'] - sim_data['s']
  
  theta <- estimate_mu_sigma(sim_data)
  
  npr <- npreg(unlist(sim_data$q) ~ (sim_data$`i-s`))
  sim_data['ccp'] = predict(npr)
  
  theta_ <- c(0,0)
  
  mle2(theta_, sim_data)
  
}

func <- function(params){
  f_ <- 
    function(i){
      integrand_lower <- function(s) dlnorm(s, meanlog = params$mu, sdlog = params$sigma) * (u(i, s, 1, params) - log(calc_CCP(i, s, params, V.fn )))
      #print("done with lower")
      integration_upper <- 1 - plnorm(i, meanlog = params$mu, sdlog = params$sigma) * (u(i, i, 1, params) - log(calc_CCP(i, i, params, V.fn )))
      #print("done with upper")
      c(
        integrate( #lower part
          integrand_lower, 
          lower = 0,
          upper = i
        )$value,
        integration_upper
      ) %>% 
        sum
      #print("done here")
    }
  #print(v)
  return(Vectorize(f_))
}

## Value function iteration, V~
f_approximation <- function(params) {
  epsilon <- 10^-4
  delta <- 1 + epsilon
  eval.points <- seq(0, 10, length.out = 20)
  f <- function(x) x
  while(delta > epsilon){
    f. <- f
    f <- func(f, params)
    f <- 
      ipol(
        val = f,
        intervals = c(0, 10),
        dims = 20
      )
    delta <- max(abs(f(eval.points) - f.(eval.points)))
    print(delta)
  }
  f
}
mle2 <- function(theta, sim_data){
  obj <- optim(par = unlist(theta),
               fn = f_approximation,
               sim_data = sim_data,
               control = list(fnscale = -1), 
               method = 'L-BFGS-B',
               lower = c(0,0)
  )
  print(obj$par)
}