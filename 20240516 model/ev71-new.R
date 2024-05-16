
library(openxlsx)
library(readr)
library(dplyr)
library(odin)
library(odin.dust)
library(mcstate)
library(coda)

# ------

gen_sir <- odin.dust::odin_dust('sir.R')


# 2005-2019

data <- read.xlsx('ev71.xlsx') %>% 
  filter(year < 2020) %>%
  arrange(week.no) %>% select(week.no, value)

# plot(data, type = 'l')

# table(data$value)

data_p <- mcstate::particle_filter_data(data, time = "week.no", rate = 7, initial_time = 0)

# ------

pop <- read.xlsx('population.xlsx', sheet = 'Sheet2') %>% 
  filter(Year >= 2005) %>%
  mutate(birth_rate = Live.births/(Total*1000*52*7), death_rate = Deaths/(Total*1000*52*7),
         birth = Live.births/(52*7), death = Deaths/(52*7)) %>%         
  
  mutate(day = c(1 + seq(0, (17*52*7), (52*7)), (19*52*7))) %>%    # 2005-2023
  
  dplyr::select(Year, day, Total, birth_rate, death_rate, birth, death)

birth_rate_step <- as.data.frame(spline(pop$day, pop$birth_rate, n = (52*7*19)))$y
death_rate_step <- as.data.frame(spline(pop$day, pop$death_rate, n = (52*7*19)))$y
summary(birth_rate_step)
summary(death_rate_step)

# matplot(seq(1,52*7*19), birth_rate_step, type = "l",col = 'red', lty = 1, ylim = c(1.5e-05, 3.7e-05))
# matlines(seq(1,52*7*19), death_rate_step, col = 'blue', lty = 1)


# pop.total <- as.data.frame(spline(pop$day, pop$Total*1000, n = (52*7*19)))
# pop.birth <- as.data.frame(spline(pop$day, pop$birth, n = (52*7*19)))
# pop.death <- as.data.frame(spline(pop$day, pop$death, n = (52*7*19)))



index <- function(info) {
  list(run = c(I_case = info$index$I_case), 
       state = c(t = info$index$time,
                 N = info$index$N,
                 S = info$index$S,
                 I = info$index$I,
                 I_case = info$index$I_case))   
}


compare <- function(state, observed, pars = NULL) {
  exp_noise <- 1e6  # A small amount of noise is added to prevent zero expectations
  
  incidence_modelled <- state['I_case', , drop = TRUE]
  incidence_observed <- observed$value
  
  lambda <- incidence_modelled +
    rexp(n = length(incidence_modelled), rate = exp_noise)
  dpois(x = incidence_observed, lambda = lambda, log = TRUE)
}


filter <- mcstate::particle_filter$new(data = data_p, model = gen_sir,
                                       n_particles = 1000L, seed = 1L,
                                       compare = compare, index = index)

prior <- list(
  mcstate::pmcmc_parameter("beta_0", 1.4, min = 0, prior = function(p)
    dunif(p, min = 0, max = 10, log = TRUE)),
  
  mcstate::pmcmc_parameter("beta_s", 0.13, min = 0, max = 1, prior = function(p)
    dbeta(p, 2, 10, log = TRUE)),
  
  mcstate::pmcmc_parameter("i_0", 3e-05, min = 0, max = 1, prior = function(p)
    dbeta(p, 2, 100, log = TRUE)),
  
  mcstate::pmcmc_parameter("p", 2e-04, min = 0, max = 1, prior = function(p)
    dbeta(p, 2, 100, log = TRUE))
)


transform <- function(theta) {
  as.list(theta)
}

make_transform <- function(beta_0, beta_s, i_0, p) {
  function(theta) {
    list(beta_0 = theta[["beta_0"]],
         beta_s = theta[["beta_s"]],
         i_0 = theta[["i_0"]],
         p = theta[["p"]],
         birth_rate_step = birth_rate_step, 
         death_rate_step = death_rate_step)
  }
}

transform <- make_transform(beta_0, beta_s, i_0, p)


control <- mcstate::pmcmc_control(
  n_steps = 5000,
  n_threads_total = 8,
  save_state = TRUE,
  save_trajectories = TRUE,
  progress = TRUE,
  rerun_every = 50,
  rerun_random = TRUE
)

# (1.4, 0.13, 3e-05, 2e-04)
proposal_matrix <- diag(c(0.2, 0.01, 1e-09, 1e-08), 4)

mcmc_pars <- mcstate::pmcmc_parameters$new(prior, proposal_matrix, transform)
# mcmc_pars$model(mcmc_pars$initial())

pmcmc <- mcstate::pmcmc(mcmc_pars, filter, control = control)

# mcmc <- coda::as.mcmc(cbind(pmcmc$probabilities, pmcmc$pars))
# summary(mcmc)
# coda::effectiveSize(mcmc)
# 1 - coda::rejectionRate(mcmc)
# plot(mcmc)

# ----------------

proposal_matrix <- cov(pmcmc$pars)
proposal_matrix
mcmc_pars <- mcstate::pmcmc_parameters$new(prior, proposal_matrix, transform)

control <- mcstate::pmcmc_control(
  n_steps = 10000,
  save_state = TRUE,
  save_trajectories = TRUE,
  progress = TRUE,
  n_chains = 4,
  n_workers = 4,
  n_threads_total = 8,
  rerun_every = 50,
  rerun_random = TRUE
)


pmcmc1 <- mcstate::pmcmc(mcmc_pars, filter, control = control)




