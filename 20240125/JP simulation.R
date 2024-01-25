
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

setwd('C:/Users/xy1123/OneDrive - Imperial College London/My PhD/phd')



data0 <- read.xlsx('data/JPdata/clean data/hfmd relate ev by week 4 main serotype 2005-2023.xlsx') %>% 
  filter(virus == 'EV-A71') %>% 
  # filter(year<=2019) %>%
  select(time, value) %>% 
  mutate(week = time) %>% select(-time)


## 1.Stochastic SIR model------------------------------------------------------

sir.seasonal <- odin.dust::odin_dust({

  # Parameters unit weekly
  beta0 <- 520/52     # average transmission rate (x people effectively contacted per person per day)
  betas <- 0.13        # amplitude of the transmission rate
  phi <- 0.75          # phase of the transmission rate (the time of annual seasonal peak, AUG or SEP)
  period <- 52        # period of seasonality
  gamma <- 1/1         # recovery rate


  birth_rate_step[] <- user()
  dim(birth_rate_step) <- user()
  birth_rate <- if (step >= length(birth_rate_step))
    birth_rate_step[length(birth_rate_step)] else birth_rate_step[step]
  
  death_rate_step[] <- user()
  dim(death_rate_step) <- user()
  # death_rate <- death_rate_step[step]
  death_rate <- if (step >= length(death_rate_step))
    death_rate_step[length(death_rate_step)] else death_rate_step[step]
  

  # Definition of the time-step and output as "time"  by day
  dt <- 1/10
  initial(time) <- 0
  update(time) <- (step + 1) * dt
  
  
  # seasonal beta
  
  beta <- beta0 * (1 + betas * cos(2 * 3.141593 * (time  + phi * period) / period))

  
  # Initial values
  N <- 127768000   # total population in 2005 year
  p <- 1.96E-04    # reporting probability per infection
  
  s_0 <- 0.0993
  i_0 <- 3.10e-05   # 0.000031
  initial(S) <- N * s_0
  initial(I) <- N * i_0
  initial(R) <- N * (1 - s_0 - i_0)
  initial(inc) <- 0
  initial(C) <- 0

  # Stochastic terms
  # Two types of events for S, so competing hazards.
  n_events_S <- rbinom(S, (beta * I / N + death_rate) * dt)
  n_deaths_S <- rbinom(n_events_S, death_rate/(beta * I / N + death_rate))
  n_infections_S <- n_events_S - n_deaths_S
  
  # Two types of events for I, so competing hazards.
  
  n_events_I <- rbinom(I, (gamma + death_rate) * dt)
  # n_events_I <- I # gamma = 1
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
  update(inc) <- n_infections_S     # infection
  update(C) <- rpois(n_infections_S * p)     # observed case

})


# prepare time vary death and birth rate by spline

pop0 <- read.xlsx('data/JPdata/clean data/Japan population.xlsx', sheet = 'Sheet1')

pop <- pop0 %>% 
  filter(Year >= 2005) %>%
  mutate(birth_rate = Live.births/(Total*52), death_rate = Deaths/(Total*52),
         time = c((1 + seq(0,16*52,52)), 52*18)) %>% 
  dplyr::select(time, birth_rate, death_rate)

birth_rate_step <- as.data.frame(spline(pop$time, pop$birth_rate, n=52*18))$y
death_rate_step <- as.data.frame(spline(pop$time, pop$death_rate, n=52*18))$y


model <- sir.seasonal$new(pars = list(birth_rate_step = birth_rate_step,
                                      death_rate_step = death_rate_step),
                          time = 1,n_threads = 4L,
                          seed = 1L,
                          n_particles = 10L)

dt <- 1/10
n_times <- 19*52/dt  # 2005-2023
n_particles <- 10L   # repeat
x <- array(NA, dim = c(model$info()$len, n_particles, n_times))


for (t in seq_len(n_times)) {
  x[ , , t] <- model$run(t)}

time <- x[1, 1, ]
head(time)
x <- x[-1, , ]
# head(x)

# ---------

cols <- c(C = "#8c8cd9", I = "#cc0044", S = "#999966")

matplot(time, t(x[5, , ]), type = "l", ylim = range(data0$value),
        xlab = "Time", ylab = "Number of case",
        col = cols[["C"]], lty = 1)
matlines(data0$week, data0$value, type = "l", col = "#cc0044") # pch = 19, cex = 0.4, 


# case.m <- as.data.frame(t(x[5, , ])) %>% 
#   mutate(week = seq(1,52*10*10,1)/10) %>%
#   melt(id.vars = "week",variable.name = "repeats")  
# 
# ggplot() +
#   geom_line(data = case.m, aes(x = week, y = value,fill = repeats),size = 1, color = 'gray')+
#   geom_line(data = data0, aes(x = week, y = value), size = 0.7, color = '#cc0044')+
#   scale_x_continuous(name="",expand = c(0,10), breaks = seq(26,(52*19),52),
#                      labels = seq(2005,2023,1)) +
#   theme(panel.background = element_blank(),
#         axis.line = element_line(colour = "black", size = 0.5),
#         axis.ticks = element_line(colour = "black", size = 0.5),
#         axis.ticks.length = unit(0.4,"lines"),
#         plot.margin = margin(0.5,0.2,0.2,0.5, "cm"),
#         axis.text = element_text(colour = "black",size = 11,hjust = 0.5,vjust = 0.5),
#         axis.title = element_text(size = 11,vjust = 0, hjust = 0.5,face = "bold"))



## 2.Inferring parameters---------------------------------------------------

# rate=1/dt
data <- mcstate::particle_filter_data(data0, time = "week", rate = 1/dt, initial_time = 0)
head(data)
str(data)


model$info()

index <- function(info) {
  list(run = c(Case = info$index$C), # incidence = cases simulated
       state = c(t = info$index$time,
                 I = info$index$I,
                 Inc = info$index$inc,
                 Case = info$index$C))   # simulated cases
}
index(model$info())

# likelihood, the probability of the simulation run given the data, by particle filter

# incidence_modelled <- model$state(6, , drop = TRUE)

compare <- function(state, observed, pars = NULL) {
  exp_noise <- 1e6  # A small amount of noise is added to prevent zero expectations
  
  incidence_modelled <- state['Case', , drop = TRUE]
  incidence_observed <- observed$value
  lambda <- incidence_modelled +
    rexp(n = length(incidence_modelled), rate = exp_noise)
  dpois(x = incidence_observed, lambda = lambda, log = TRUE)
}


filter <- mcstate::particle_filter$new(data = data,
                                       model = sir.seasonal,
                                       n_particles = 100,
                                       compare = compare,
                                       index = index,
                                       seed = 1L)

filter$run(pars = list(birth_rate_step = birth_rate_step,
                       death_rate_step = death_rate_step),
           save_history = TRUE) #save the history allows to plot the particle trajectories


# plot these along with the data，only trajectories consistent with the data are kept

h <- filter$history()
head(h)

par(mar = c(4, 4, 0.5, 0.5)) # 下左上右
matplot(data0$week, data0$value, type = 'l', col = "#cc0044",
        xlab = "Week", ylab = "Cases")
matlines(h[1,1,], t(h[4, , ]), type = "l",lwd = 0.2, col = "#00000011")


## 3.Using MCMC to infer parameters----------------------------------------------



