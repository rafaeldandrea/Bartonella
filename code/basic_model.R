library(deSolve)
library(tidyverse)
library(cowplot)
theme_set(theme_bw())


# Functions -------------------------------------------------------------------

# Apply smoothed step function to array: values less than 0 are set to 0;
# values between 0 and 1 are set to 10*n^3-15*n^4+6*n^5; and values larger
# than 1 are set to 1
# Input:
# - n: vector or arbitrary array of values
# Output:
# - array of values with smoothed step function applied to them
cutoff <- function(n) {
  return(ifelse(n<1, (1*(n>0))*(n*n*n*(10+n*(-15+6*n))), 1))
}

# Convert mobility values to peak infectivities
mobility_intensity <- function(nu) {
  return(1 / nu)
}

# Right hand side of host-pathogen community model
host_pathogen_model <- function(time, state, pars) {
  Sh <- length(pars$sigmah) # Number of host species
  Sp <- length(pars$sigmap) # Number of pathogen species
  cN <- 1e-8 # cutoff threshold for densities
  cT <- 1e-6 # cutoff threshold for trait evolution
  # Matrix N_ks where k indexes the host species and s the stage:
  Nks <- matrix(state[1:(Sh * (2 + Sp))], nrow = Sh, ncol = Sp + 2)
  S <- Nks[, 1] # Susceptibles for each species
  I <- matrix(Nks[, -c(1, Sp + 2)], Sh, Sp) # Infecteds per host per pathogen
  R <- Nks[, Sp + 2] # Recovereds for each species
  N <- rowSums(Nks) # Total population size per host species
  x <- state[Sh * (2 + Sp) + 1:Sh] # Host trait values
  z <- state[Sh * (3 + Sp) + 1:Sp] # Pathogen trait values
  dxz <- outer(x, z, FUN = `-`) ## Difference matrix of trait means (host - pathogen)
  eta <- mobility_intensity(pars$nu) # Intensity as a function of host mobility
  xi <- eta * exp(-dxz^2 / (2 * pars$omega^2)) / (pars$omega * sqrt(2 * pi)) ## infectivity rates, a matrix
  beta <- aperm(pars$nu %o% (pars$nu * xi), c(2, 1, 3)) # Transmission rates. beta_{kil} := nu_k * nu_i * xi_{kl}
  partial_beta <- sweep(beta, MARGIN = c(1, 3), STATS = dxz / pars$omega ^ 2, FUN = '*')
  SI = S %o% I ## SI_{kil} := S_k * I_{il}
  Iinv = ifelse(I > 0, 1 / I, 0)
  SII = sweep(SI, c(1, 3), Iinv, '*') ## SII_{kil} := S_k / I_{kl} * I_{il}
  # Effective fecundities (modified by infection status):
  phi_matrix <- pars$phi * cbind(1, 1 - exp(-dxz^2 / (2 * pars$theta^2)), 1)
  # ODEs:
  dSdt <- (rowSums(phi_matrix * Nks) - pars$alpha * N^2 - rowSums(beta * SI) -
             pars$mu * S) #* cutoff(S / cN) # Susceptibles
  dIdt <- (apply(beta * SI, c(1, 3), sum) - pars$gamma * I - pars$mu * I) #* cutoff(I / cN) # Infecteds
  dRdt <- (rowSums(pars$gamma * I) - pars$mu * R) # Recovereds
  dxdt <- pars$sigmah^2 * rowSums(pars$phi * dxz * exp(-dxz^2 / (2 * pars$theta^2)) /
                                    pars$theta^2) #* cutoff(N / cT) # Host traits
  dzdt <- pars$sigmap^2 * rowSums(aperm(SII * partial_beta, c(3, 1, 2))) #* cutoff(colSums(I) / cT) # Pathogen traits
  return(list(c(dSdt, dIdt, dRdt, dxdt, dzdt)))
}

# Solve model equations and organize them into a tidy tibble
solve_ode <- function(pars, integration_method = "lsoda") {
  # Rename data columns to match true variable names
  rename_columns <- function(sol, Sh, Sp) {
    names(sol)[1] <- "time" # The first column is "time"
    index <- 1 # Column counter index, to keep track of which columns we are renaming
    names(sol)[index + 1:Sh] <- paste0("S_", 1:Sh) # Susceptibles
    index <- index + Sh # Increment column counter index
    for (k in 1:Sp) { # For each pathogen species k:
      names(sol)[index + 1:Sh] <- paste0("I", k, "_", 1:Sh) # Infecteds
      index <- index + Sh # Increment column counter index
    }
    names(sol)[index + 1:Sh] <- paste0("R_", 1:Sh) # Recovereds
    index <- index + Sh # Increment column counter index
    names(sol)[index + 1:Sh] <- paste0("x_", 1:Sh) # Host trait values
    index <- index + Sh # Increment column counter index
    names(sol)[index + 1:Sp] <- paste0("z_", 1:Sp) # Pathogen trait values
    return(sol)
  }
  Sh <- length(pars$sigmah) # Number of host species
  Sp <- length(pars$sigmap) # Number of pathogen species
  # We now solve the equations and put them in a tidy form:
  sol <- ode(y = pars$initial_conditions, times = pars$time_sequence,
             func = pars$ode_system, parms = pars, method = integration_method) %>%
    as.data.frame() %>% # Convert to a data frame (needed for next step)
    as_tibble() %>% # Convert to a tibble
    rename_columns(Sh, Sp) # Rename the columns
  densities <- sol %>% # A tibble for the population densities
    select(!starts_with(c("x", "z"))) %>%
    pivot_longer(cols = !"time", names_to = "var", values_to = "density") %>%
    separate(col = "var", into = c("stage", "species"), sep = "_") %>%
    mutate(stage = fct_relevel(stage, c("S", paste0("I", 1:Sp), "R")))
  traits <- sol %>% # A tibble for host & pathogen trait values
    select(time, starts_with(c("x", "z"))) %>%
    pivot_longer(cols = !"time", names_to = "var", values_to = "trait") %>%
    separate(col = "var", into = c("role", "species"), sep = "_") %>%
    mutate(role = if_else(role == "x", "host", "pathogen"))
  # Return two tibbles (for densities and traits) in a list:
  return(list(densities = densities, traits = traits))
}


# Solve equations -------------------------------------------------------------

solution <- tibble( # Define table of parameters
  Sh = 2, # Number of host species
  Sp = 2, # Number of pathogen strains
  sigmah = list(rep(0, Sh)), # Vector of genetic variances of the host species
  sigmap = list(rep(0.1, Sp)), # Vector of genetic variances of the pathogen species
  nu = list(rep(1, Sh)), # Vector of mobilities
  phi = list(rep(5, Sh)), # Vector of host intrinsic growth rates
  mu = list(rep(1, Sh)), # Vector of host mortalities
  gamma = list(matrix(2, Sh, Sp)), # Matrix of recovery rates
  alpha = 0.1, # Competition coefficient
  omega = 0.15, # Width of infectivity curve
  theta = 0.5, # Width of the trait-matching curve for fecundity impact
  time_sequence = list(seq(0, 50, by = 0.1)),
  initial_conditions = list(c(rep(1, Sh), matrix(0.01, Sh, Sp), rep(0, Sh),
                              seq(0, 1, l = Sh), seq(0.1, 0.5, l = Sp))),
  ode_system = list(host_pathogen_model) # Function to evaluate in solving model
) %>%
  transmute(pars = mapply( # Gather all parameters in single column
    list, sigmah = sigmah, sigmap = sigmap, nu = nu, phi = phi, mu = mu,
    gamma = gamma, omega = omega, alpha = alpha, theta = theta,
    time_sequence = time_sequence, initial_conditions = initial_conditions,
    ode_system = ode_system, SIMPLIFY = FALSE)) %>%
  # Integrate model & put result into a new column, as a nested tibble:
  mutate(sol = map(pars, solve_ode, integration_method = "lsoda")) %>%
  select(sol) %>% # Keep only the column with the solution
  unnest(cols = sol) %>% # Unnest the nested data frame
  mutate(id = c("density", "trait")) %>% # Identify the two rows
  pivot_wider(names_from = "id", values_from = "sol")

plot_grid(
  solution$density[[1]] %>%
    mutate(species = paste("host species", species)) %>%
    ggplot() +
    aes(x = time, y = density, colour = stage) +
    geom_line() +
    facet_grid(. ~ species),
  solution$trait[[1]] %>%
    mutate(species = paste("host species", species)) %>%
    ggplot() +
    aes(x = time, y = trait, colour = role) +
    geom_line() +
    facet_grid(. ~ species),
  ncol = 1, align = "hv"
) %>%
  show
