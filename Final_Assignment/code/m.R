source('header.R')

params = list(
  Q = 20,
  mu = 1,
  sigma = 1,
  p_r = 5,
  p_w = 1,
  alpha = 0.3,
  eta_q = 4,
  beta = 0.95
) %>% 
  value_fn_iteration %>% 
  simulation(params = params) %>%
  param_estimation