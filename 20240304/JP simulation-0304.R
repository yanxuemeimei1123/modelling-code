
library(openxlsx)
library(readxl)
library(sf)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(odin)
library(odin.dust)
library(mcstate)


# 1.SIR model with odin.dust----------------------------------------------------

gen_sir <- odin.dust::odin_dust({
  
  # Parameters unit daily
  beta_0 <- user()     # average transmission rate (x people effectively contacted per person per day)
  beta_s <- user()     # amplitude of the transmission rate
  phi <- 0.75          # phase of the transmission rate (the time of annual seasonal peak, AUG or SEP)
  period <- 365        # period of seasonality
  gamma <- 1/7         # recovery rate
  
  birth_rate <- birth_rate_step[step]
  death_rate <- death_rate_step[step]
  
  birth_rate_step[] <- user()
  dim(birth_rate_step) <- user()
  death_rate_step[] <- user()
  dim(death_rate_step) <- user()
  
  
  # Definition of the time-step and output as "time"
  dt <- 1  
  initial(time) <- 0
  update(time) <- (step + 1) * dt
  
  # Initial values
  s_0 <- user()        # proportions of S
  i_0 <- user()        # proportions of I
  N_0 <- 127768000     # total population in 2005 year
  initial(N) <- N_0    
  initial(S) <- N_0 * s_0
  initial(I) <- N_0 * i_0
  initial(R) <- N_0 * (1 - s_0 - i_0)  
  
  # update for the next time step
  update(N) <- S + I + R
  update(S) <- S - n_deaths_S - n_infections_S + n_births
  update(I) <- I + n_infections_S - n_recoveries_I - n_deaths_I
  update(R) <- R + n_recoveries_I - n_deaths_R
  
  # Daily infections incidence
  update(I_inc) <- if (step %% 7 == 0) n_infections_S else I_inc + n_infections_S
  initial(I_inc) <- 0
  
  # reported case
  p <- user()                           # reporting probability per infection
  update(I_case) <- rpois(I_inc * p)    # observed case
  initial(I_case) <- 0
  
  # beta
  beta <- beta_0 * (1 + beta_s * cos(2 * 3.141593 * (time  + phi * period) / period))
  w <- user()   # imported case
  
  # Two types of events for S, so competing hazards.
  
  p_SI <- if (death_rate + beta * (I+w) / N > 1) 1 else (death_rate + beta * (I+w) / N)  # p < 1
  n_events_S <- if (p_SI > 0) rbinom(S, p_SI * dt) else 0                                # p > 0
  
  p_SID <- if (death_rate/p_SI < 0) 0 else
    if (death_rate/p_SI > 1) 1 else death_rate/p_SI
  
  n_deaths_S <- if (n_events_S > 1) rbinom(n_events_S, p_SID) else 0
  n_infections_S <- if (n_events_S - n_deaths_S > 0) (n_events_S - n_deaths_S) else 0
  
  # Two types of events for I, so competing hazards.
  p_IR <- death_rate + gamma                             # 0.0001626~0.0002408
  n_events_I <- if (I > 1) rbinom(I, p_IR * dt) else 0   # 0.1429
  n_deaths_I <- if (n_events_I > 1) rbinom(n_events_I, death_rate / p_IR) else 0
  n_recoveries_I <- if (n_events_I - n_deaths_I > 0) (n_events_I - n_deaths_I) else 0
  
  # one type of events for R
  n_deaths_R <- if (R > 1) rbinom(R, death_rate * dt) else 0 # (2.323~3.440)e-05
  
  # births and deaths
  n_births <- rbinom(N, birth_rate * dt)
  
})

# setting parameters

pop <- read.xlsx('clean data/Japan population.xlsx', sheet = 'Sheet2') %>% 
  filter(Year >= 2005) %>% 
  mutate(birth_rate = Live.births/(Total*1000*365), death_rate = Deaths/(Total*1000*365),
         time = c(1 + seq(0,17*365,365), 19*365),
         week = c(1 + seq(0,17*52,52), 19*52),
         birth = Live.births/365, death = Deaths/365) %>%         # 2005-2023
  dplyr::select(Year,time, week,Total, birth_rate, death_rate, birth, death)

birth_rate_step <- as.data.frame(spline(pop$time, pop$birth_rate, n = 365*19))$y
death_rate_step <- as.data.frame(spline(pop$time, pop$death_rate, n = 365*19))$y


# 2.Inferring parameters--------------------------------------------------------

# observed data

data0 <- read.xlsx('clean data/hfmd relate ev by week 4 main serotype 2005-2023.xlsx') %>% 
  # filter(virus == 'EV-A71') %>%
  filter(year <= 2019) %>% 
  select(-week) %>%
  rename(week = time)

data_eva71 <- mcstate::particle_filter_data(subset(data0,virus == 'EV-A71'), time = "week", rate = 7, initial_time = 0)
data_cva6 <- mcstate::particle_filter_data(subset(data0,virus == 'CVA6'), time = "week", rate = 7, initial_time = 312)
data_cva10 <- mcstate::particle_filter_data(subset(data0,virus == 'CVA10'), time = "week", rate = 7, initial_time = 0)
data_cva16 <- mcstate::particle_filter_data(subset(data0,virus == 'CVA16'), time = "week", rate = 7, initial_time = 0)


# defining the comparison function

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


# inferring parameters

filter_eva71 <- mcstate::particle_filter$new(data = data_eva71, model = gen_sir,
                                             n_particles = 10000, compare = compare, index = index)
filter_cva6 <- mcstate::particle_filter$new(data = data_cva6, model = gen_sir,
                                            n_particles = 10000, compare = compare, index = index)
filter_cva10 <- mcstate::particle_filter$new(data = data_cva10, model = gen_sir,
                                             n_particles = 10000, compare = compare, index = index)
filter_cva16 <- mcstate::particle_filter$new(data = data_cva16, model = gen_sir,
                                             n_particles = 10000, compare = compare, index = index)

# 3.Using MCMC to infer parameters--------------------------------------------

# define the prior function, default is a flat improper prior = uniform distribution

prior <- list(
  mcstate::pmcmc_parameter("beta_0", 2, min = 0, max = 10),
  mcstate::pmcmc_parameter("beta_s", 0.15, min = 0, max = 1),
  mcstate::pmcmc_parameter("s_0", 0.07, min = 0, max = 1),
  mcstate::pmcmc_parameter("i_0", 4e-05, min = 0, max = 1e-03),
  mcstate::pmcmc_parameter("p", 2e-04, min = 0, max = 1e-02),
  mcstate::pmcmc_parameter("w", 50, min = 0, max = 1000))

# The proposals for each MCMC step. drawn from a multi-variate normal distribution, 
# with the variance-covariance matrix. βn+1∼N(βn, 0.01)
proposal_matrix <- diag(c(1e-02,1e-03,1e-04,1e-07,1e-06,10), 6)

transform <- function(theta) {
  as.list(theta)
}

make_transform <- function(birth_rate_step,death_rate_step) {
  function(theta) {
    list(beta_0 = theta[["beta_0"]],
         beta_s = theta[["beta_s"]],
         s_0 = theta[["s_0"]],
         i_0 = theta[["i_0"]],
         p = theta[["p"]],
         w = theta[["w"]],
         birth_rate_step = birth_rate_step, 
         death_rate_step = death_rate_step)
  }
}

transform <- make_transform(birth_rate_step,death_rate_step)

# Run the pMCMC for the first time

mcmc_pars <- mcstate::pmcmc_parameters$new(prior, proposal_matrix, transform)

control <- mcstate::pmcmc_control(
  n_steps = 1000,
  save_state = TRUE,
  save_trajectories = TRUE,
  progress = TRUE)

pmcmc_eva71 <- mcstate::pmcmc(mcmc_pars, filter_eva71, control = control)

pmcmc_cva10 <- mcstate::pmcmc(mcmc_pars, filter_cva10, control = control)

pmcmc_cva16 <- mcstate::pmcmc(mcmc_pars, filter_cva16, control = control)

pmcmc_cva6 <- mcstate::pmcmc(mcmc_pars, filter_cva6, control = control)

# Run the pMCMC for the second time - tuning the pmcmc
# by modifying the proposal: use the covariance of the state as the proposal matrix

control_2 <- mcstate::pmcmc_control(
  n_steps = 1000000,
  save_state = TRUE,
  save_trajectories = TRUE,
  progress = TRUE)

# eva-71
proposal_matrix_eva71 <- cov(pmcmc_eva71$pars)
mcmc_pars_eva71 <- mcstate::pmcmc_parameters$new(prior, proposal_matrix_eva71, transform)
pmcmc_eva71_tuning <- mcstate::pmcmc(mcmc_pars_eva71, filter_eva71, control = control_2)
mcmc_eva71_tuning <- coda::as.mcmc(cbind(pmcmc_eva71_tuning$probabilities, pmcmc_eva71_tuning$pars))
write_rds(mcmc_eva71_tuning, "pMCMC_results_RDS/mcmc_eva71_tuning.rds")

# cva-10
proposal_matrix_cva10 <- cov(pmcmc_cva10$pars)
mcmc_pars_cva10 <- mcstate::pmcmc_parameters$new(prior, proposal_matrix_cva10, transform)
pmcmc_cva10_tuning <- mcstate::pmcmc(mcmc_pars_cva10, filter_cva10, control = control_2)
mcmc_cva10_tuning <- coda::as.mcmc(cbind(pmcmc_cva10_tuning$probabilities, pmcmc_cva10_tuning$pars))
write_rds(mcmc_cva10_tuning, "pMCMC_results_RDS/mcmc_cva10_tuning.rds")

# cva-16
proposal_matrix_cva16 <- cov(pmcmc_cva16$pars)
mcmc_pars_cva16 <- mcstate::pmcmc_parameters$new(prior, proposal_matrix_cva16, transform)
pmcmc_cva16_tuning <- mcstate::pmcmc(mcmc_pars_cva16, filter_cva16, control = control_2)
mcmc_cva16_tuning <- coda::as.mcmc(cbind(pmcmc_cva16_tuning$probabilities, pmcmc_cva16_tuning$pars))
write_rds(mcmc_cva16_tuning, "pMCMC_results_RDS/mcmc_cva16_tuning.rds")

# cva-6
proposal_matrix_cva6 <- cov(pmcmc_cva6$pars)
mcmc_pars_cva6 <- mcstate::pmcmc_parameters$new(prior, proposal_matrix_cva6, transform)
pmcmc_cva6_tuning <- mcstate::pmcmc(mcmc_pars_cva6, filter_cva6, control = control_2)
mcmc_cva6_tuning <- coda::as.mcmc(cbind(pmcmc_cva6_tuning$probabilities, pmcmc_cva6_tuning$pars))
write_rds(mcmc_cva6_tuning, "pMCMC_results_RDS/mcmc_cva6_tuning.rds")

