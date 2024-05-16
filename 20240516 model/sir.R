
# gen_sir <- odin.dust::odin_dust({

# Parameters unit daily
beta_0 <- user()     # average transmission rate (x people effectively contacted per person per day)
beta_s <- user()     # amplitude of the transmission rate
phi <- user(0.7)           # phase of the transmission rate (the time of annual seasonal peak, AUG or SEP)
period <- user(52*7)       # period of seasonality
gamma <- user(1/7)         # recovery rate

birth_rate <- birth_rate_step[step]
death_rate <- death_rate_step[step]

birth_rate_step[] <- user()
dim(birth_rate_step) <- user()
death_rate_step[] <- user()
dim(death_rate_step) <- user()


# Definition of the time-step and output as "time"  unit = 1 day
dt <- 1  
initial(time) <- 0
update(time) <- (step + 1) * dt

# Initial values
R0 <- beta/(death_rate+gamma)
s_star <- 1/R0
s_0 <- s_star
# s_0 <- user()              # proportions of S
i_0 <- user()              # proportions of I
N_0 <- user(127768000)     # total population in 2005 year

initial(N) <- N_0    
initial(S) <- N_0 * s_0
initial(I) <- N_0 * i_0
initial(R) <- N_0 * (1 - s_0 - i_0)  


# beta before and after covid
# 2020-3 to 2022-5 ---- (15*52+9)*7 to (18*52+21)*7 ---- 5523-6699
# cr <- user()
# beta <- if (time <= 6699 && time >= 5523) 
#   (cr * beta_0 * (1 + beta_s * cos(2 * 3.14159265358979323846 * (time  - phi * period) / period))) else
#     (beta_0 * (1 + beta_s * cos(2 * 3.14159265358979323846 * (time  - phi * period) / period)))

beta <- beta_0 * (1 + beta_s * cos(2 * 3.14159265358979323846 * (time  - phi * period) / period))
w <- user(0)                      # imported case
foi <- beta * (I + w) / N


# Two types of events for S, so competing hazards.
p_SI_event <- 1 - exp(-(foi + death_rate))  # transition probabilities of (death and I) from S
p_SI_death <- 1 - exp(-death_rate / (foi + death_rate))

n_events_S <- rbinom(S, p_SI_event * dt)
n_deaths_S <- rbinom(n_events_S, p_SI_death)
n_infections_S <- n_events_S - n_deaths_S

# Two types of events for I, so competing hazards.
p_IR_event <- 1 - exp(-(gamma + death_rate))  
p_IR_death <- 1 - exp(-death_rate / (death_rate + gamma))

n_events_I <- rbinom(I, p_IR_event * dt)
n_deaths_I <- rbinom(n_events_I, p_IR_death)
n_recoveries_I <- n_events_I - n_deaths_I

# one type of events for R
n_deaths_R <- rbinom(R, (1 - exp(-death_rate)) * dt) 

n_births <- rbinom(N, (1 - exp(-birth_rate)) * dt)


# update for the next time step
update(N) <- S + I + R
update(S) <- S + n_births - n_infections_S - n_deaths_S
update(I) <- I + n_infections_S - n_recoveries_I - n_deaths_I
update(R) <- R + n_recoveries_I - n_deaths_R

# Daily infections incidence
update(I_inc) <- if (step %% 7 == 0) n_infections_S else I_inc + n_infections_S
initial(I_inc) <- 0

# reported case
p <- user()                           # reporting probability per infection
update(I_case) <- rpois(I_inc * p)    # observed case
initial(I_case) <- 0

initial(births) <- rbinom(N, (1 - exp(-birth_rate)) * dt)
update(births) <- rbinom(N, (1 - exp(-birth_rate)) * dt)

# })


# par_week <- list(beta_0 = 520/52, beta_s = 0.13, phi = 0.75,
#                  birth_rate_step = birth_rate_step, death_rate_step = death_rate_step,
#                  s_0 = 0.0993, i_0 = 3.10e-05, N_0 = 126931430, w = 8, p = 1.96e-04)
# 
# model <- gen_sir$new(pars = par_week,
#                      time = 1,n_threads = 4L,
#                      seed = 1L,
#                      n_particles = 1000L)

# model$simulate(0)

# x <- model$simulate(seq(1, 52*7*15))  # 2005-2019
# time <- x[1,1,]
# 
# model$info()
# 
# cols <- c(S = "#8c8cd9", I = "#cc0044", R = "#999966", All = 'black')
# 
# matplot(time, t(x[2, , ]), type = "l", ylim = range(x),
#         xlab = "Time", ylab = "Number of individuals",
#         col = cols[["All"]], lty = 1)
# matlines(time, t(x[3, , ]), col = cols[["S"]], lty = 1)
# matlines(time, t(x[4, , ]), col = cols[["I"]], lty = 1)
# matlines(time, t(x[5, , ]), col = cols[["R"]], lty = 1)
# legend("left", lwd = 1, col = cols, legend = names(cols), bty = "n")
# 
# matplot(time, t(x[10, , ]), type = "l", col = cols[["I"]], lty = 1)  # , ylim = c(0,100000)
# summary(c(x[9, , ]))

