
library(odin)
library(odin.dust)
library(mcstate)


####---------seasonal stochastic SIR with birth and death-------------------####

sir.seasonal <- odin.dust::odin_dust({
  
  # Parameters
  beta0 <- 0.3           # baseline transmission rate
  gamma <- 1/7           # recovery rate
  N <- 12693.1430        # total population
  seasonality <- 0.2     # amplitude of seasonality
  period <- 365          # period of seasonality (assuming daily time steps)
  phi <- 0.7             # phase of transmission rate
  birth_rate <- 0.0085   # birth rate
  death_rate <- 0.0088   # death rate
  
  # # time-vary
  # death_rate <- if (time <= 365) 0.0085
  #               else if (time <= 365*2) 0.0075
  #               else if (time <= 365*3) 0.0080
  #               else 0.0090
  # 
  # birth_rate <- if (time <= 365) 0.0090
  #               else if (time <= 365*2) 0.0080
  #               else if (time <= 365*3) 0.0075
  #               else 0.0085
  
  # Definition of the time-step and output as "time"
  dt <- user(1)
  initial(time) <- 0
  update(time) <- (step + 1) * dt
  
  # seasonal beta
  beta <- beta0 * (1 + seasonality * sin(2 * 3.141593 * time / period + phi))
  
  # Initial values
  s_0 <- 0.7
  i_0 <- 0.01
  initial(S) <- N * s_0
  initial(I) <- N * i_0
  initial(R) <- (1 - s_0 - i_0) * N
  
  # Stochastic terms
  # Two types of events for S, so competing hazards.
  n_events_S <- rbinom(S, (beta * I / N + death_rate) * dt)
  n_deaths_S <- rbinom(n_events_S, death_rate/(beta * I / N + death_rate))
  n_infections_S <- n_events_S - n_deaths_S
  
  # Two types of events for I, so competing hazards.
  n_events_I <- rbinom(I, (gamma + death_rate) * dt)
  n_deaths_I <- rbinom(n_events_I, death_rate / (death_rate + gamma))
  n_recoveries_I <- n_events_I - n_deaths_I
  
  # one type of events for R
  n_deaths_R <- rbinom(R, death_rate * dt)
  
  # births
  n_births <- birth_rate * N
  
  # update for the next time step
  update(S) <- S - n_deaths_S - n_infections_S + n_births
  update(I) <- I + n_infections_S - n_recoveries_I - n_deaths_I
  update(R) <- R + n_recoveries_I - n_deaths_R
  
})

sir.seasonal_model <- sir.seasonal$new(pars = list(birth_rate <- 0.0085, 
                                                   death_rate <- 0.0088),
                                       time = 1,n_threads = 4L,
                                       seed = 1L,
                                       n_particles = 10L)


sir.seasonal_model$state()    # The initial state
sir.seasonal_model$info()     # index of (t, S, I, R)
# Run the particles (repeats) forward 10 time steps of length dt
# sir.seasonal_model$run(10)        #runs up to some time, returns final values(*)
# sir.seasonal_model$simulate(10)   #runs over a series of times, returning values at each time

n_times <- 365 *4
n_particles <- 10L   #repeat
x <- array(NA, dim = c(sir.seasonal_model$info()$len, n_particles, n_times))
head(x)

for (t in seq_len(n_times)) {
  x[ , , t] <- sir.seasonal_model$run(t)}

time <- x[1, 1, ]
x <- x[-1, , ]

par(mar = c(4.1, 5.1, 0.5, 0.5), las = 1)
cols <- c(S = "#8c8cd9", I = "#cc0044", R = "#999966")
matplot(time, t(x[1, , ]), type = "l",
        xlab = "Time", ylab = "Number of individuals",
        col = cols[["S"]], lty = 1, ylim = range(x))
matlines(time, t(x[2, , ]), col = cols[["I"]], lty = 1)
matlines(time, t(x[3, , ]), col = cols[["R"]], lty = 1)
legend("left", lwd = 1, col = cols, legend = names(cols), bty = "n")

matplot(time, t(x[2, , ]), type = "l",
        xlab = "Time", ylab = "Number of infected individuals",
        col = cols[["I"]], lty = 1, ylim = range(x[2, ,]))





####-------------- total simulation and fitting----------------- ####

# use simple SIR

incidence <- read.csv("incidence.csv")

sir <- odin.dust::odin_dust({
  
  N <- S + I + R
  p_SI <- 1 - exp(-(beta) * I / N)
  p_IR <- 1 - exp(-(gamma))
  n_IR <- rbinom(I, p_IR * dt)
  n_SI <- rbinom(S, p_SI * dt)
  
  update(time) <- (step + 1) * dt
  update(S) <- S - n_SI
  update(I) <- I + n_SI - n_IR
  update(R) <- R + n_IR
  update(cases_cumul) <- cases_cumul + n_SI    # accumulated infections include n_IR
  update(cases_inc) <- if (step %% freq == 0) n_SI else cases_inc + n_SI     # infections include n_IR
  
  initial(time) <- 0
  initial(S) <- 1000
  initial(R) <- 0
  initial(I) <- I0
  initial(cases_cumul) <- 0
  initial(cases_inc) <- 0
  
  beta <- user(0.2)
  gamma <- user(0.1)
  I0 <- user(10)
  
  freq <- user(4)
  dt <- 1.0 / freq
  
})


#  convert OBSERVED DATA into input for the particle filter

data <- mcstate::particle_filter_data(incidence, time = "day", rate = 4, initial_time = 0)
head(data)


pars <- list(beta = 0.25, gamma = 0.1)
mod <- sir$new(pars, 0, 20)
y <- mod$simulate(c(0, data$time_end))
i <- mod$info()$index[["time"]]
j <- mod$info()$index[["cases_inc"]]
matplot(y[i, 1, ], t(y[j, , ]), type = "l", col = "#00000055", lty = 1, las = 1,
        xlab = "Day", ylab = "Cases")
points(cases ~ day, incidence, col = "red", pch = 19)


# calculate the log-likelihood of the model run given the data.

index <- function(info) {
  list(run = c(incidence = info$index$cases_inc), # incidence = cases simulated
       state = c(t = info$index$time,
                 I = info$index$I,
                 cases = info$index$cases_inc))   # simulated cases
}
index(mod$info())


# the number of expected new cases (if all cases are observed) can be assumed to
# be Poisson distributed with mean yt (observed daily case counts)

compare <- function(state, observed, pars = NULL) {
  modelled <- state["incidence", , drop = TRUE]      # simulated cases (case_inc)
  lambda <- modelled + rexp(length(modelled), 1e6)   # ?? rexp-random exp number
  dpois(observed$cases, lambda, log = TRUE)          # d-density distribution
}


# index(mod$info())
# a <- mod$state()[6,, drop = TRUE]
# length(a)
# b <- a + rexp(length(a), 1e6)
# c <- dpois(data$cases, b, log = TRUE)
# plot(c)


# Only those trajectories consistent with the data are continued forward at each step,
# and a final log-likelihood of the model parameters given the data is produced.

filter <- mcstate::particle_filter$new(data, model = sir, n_particles = 100,
                                       compare = compare, index = index)

# log-likelihood
filter$run(pars = list(beta = 0.25, gamma = 0.1))


# run the filter while saving history (off by default)
filter <- mcstate::particle_filter$new(data, model = sir, n_particles = 100,
                                       compare = compare, index = index)
filter$run(save_history = TRUE)


times <- data$time_end
h <- filter$history()   # t I cases, 100 days, 100 repeats
head(h)

matplot(h["t", 1, ], t(h["cases", , ]), type = "l", col = "#00000011", 
        xlab = "Day", ylab = "Cases", las = 1)
points(cases ~ day, incidence, pch = 19, col = "red")

# unobserved states     # all case
matplot(h["t", 1, ], t(h["I", , ]), type = "l", col = "#00000011", 
        xlab = "Day", ylab = "Number of infecteds (I)", las = 1)



# estimate the posterior density for β and γ, are achieved by defining sampling
# distributions for the parameters, and running a set of MCMC chains.

beta <- pmcmc_parameter("beta", 0.2, min = 0)
gamma <- pmcmc_parameter("gamma", 0.1, min = 0,
                         prior = function(p)
                           dgamma(p, shape = 1, scale = 0.2, log = TRUE)
)

proposal_matrix <- diag(2) * 0.01^2
mcmc_pars <- pmcmc_parameters$new(list(beta = beta, gamma = gamma),
                                  proposal_matrix)




# Priors
priors <- list(
  mcstate::pmcmc_parameter("beta", 0.2, min = 0),
  mcstate::pmcmc_parameter("gamma", 0.1, min = 0, prior = function(p)
    dgamma(p, shape = 1, scale = 0.2, log = TRUE)))

# Proposal
vcv <- diag(0.01, 2)
vcv

# Transformation: Convert “MCMC parameters” into “model parameters”
transform <- function(theta) {
  as.list(theta)
}

# Final parameter object
mcmc_pars <- mcstate::pmcmc_parameters$new(priors, vcv, transform)




# Running PMCMC
control <- mcstate::pmcmc_control(
  n_steps = 500,
  progress = TRUE)
samples <- mcstate::pmcmc(mcmc_pars, filter, control = control)
samples


# PMCMC samples
plot(samples$probabilities[, "log_posterior"], type = "s",
     xlab = "Sample", ylab = "Log posterior")


hist(samples$pars[,'beta'])
hist(samples$pars[,'gamma'])
hist(samples$pars[,'beta']/samples$pars[,'gamma'])


