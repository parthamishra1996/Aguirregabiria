my.library <- 'C:/Users/mymai/Documents/R/win-library/3.6'
.libPaths(my.library)

library('tidyverse')
c('reshape2', 'stringr', 'magrittr', 'doParallel', 'chebpol', 'np') %>%
  walk(~library(., character.only=TRUE))

dir('modules') %>% 
  walk(~source(paste('./modules/', ., sep="")))

var_save <- '../variables/'

registerDoParallel(cores=28)