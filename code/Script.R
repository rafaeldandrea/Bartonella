## Model on infectious diseases on a lattice with diversification of hosts and pathogens

## Determine whether working on SeaWulf (SBU clusters) or personal computer
seawulf = as.logical(Sys.info()['user'] == 'rdrocha')

## Read array argument from command line
index = as.numeric(commandArgs(TRUE))[1]
if(is.na(index)) index = 1

## load libraries
library(tidyverse)
library(magrittr)
library(furrr)
library(parallel)

## operational variables
do.dynamics = FALSE
do.save = seawulf
do.plots = FALSE
do.parallel = FALSE
do.wrangling = TRUE

## set multicore parallel computing if running on clusters
if(seawulf & do.parallel) 
  plan(multisession, workers = detectCores())


## simulation parameters
parameters = 
  expand_grid(
    J0 = 1e2,                              # size of rat community
    city_size = 1e3,                       # size of city in meters
    maxtime = 1e4,                         # simulation iterations
    infection_radius = c(100, 150, 200),   # radius of infection
    turnover_factor = .1,                  # proportion of deaths/births per iteration
    rat_mutation_rate = .01,               # probability of rat speciation per birth
    rat_mutation_amplitude = 5,            # mutation step ~ N(0, rma)
    infection_mutation_rate = .01,         # probability of disease speciation per infection
    infection_mutation_amplitude = .05,    # mutation step ~ N(0, ima)
    mortality = .1,                        # intrinsic rat mortality
    density_weight = c(0, 1),              # contribution of local density to mortality
    coef1 = 1,                             # contribution of the 1st power of virulence to mortality
    coef2 = c(0, 1),                       # contribution of the 2nd power of virulence to mortality
    virulence0 = 1,                        # initial virulence
    mobility0 = 40,                        # initial mobility
    seed = 1:5                             # seed generator
  ) %>%
  rowid_to_column(var = 'scenario')

parms = 
  parameters %>%
  split(.$infection_radius) %>%
  extract2(index)

## Functions
dist_torus = function(dtf, parms){
  
  rectangle_sides = with(parms, c(city_size, city_size))
  
  dx = rectangle_sides[1] / 2 - abs(as.matrix(dist(dtf$x)) - rectangle_sides[1] / 2)
  dy = rectangle_sides[2] / 2 - abs(as.matrix(dist(dtf$y)) - rectangle_sides[2] / 2)
  dists = sqrt(dx ^ 2 + dy ^ 2)
  return(dists)
}

move_rats = function(dtf, parms, distances){
  
  city_size = parms$city_size
  radius = parms$infection_radius
  
  distance_moved = abs(rnorm(nrow(dtf), mean = 0, sd = dtf$mobility))
  angle_moved = runif(nrow(dtf), min = 0, max = 2 * pi)
  
  foo =
    dtf %>%
    mutate(
      x = mod(x + distance_moved * cos(angle_moved), city_size),
      y = mod(y + distance_moved * sin(angle_moved), city_size)
    ) %>%
    mutate(
      local_density = rowSums(distances < radius)
    ) %>%
    return
  
}

infect_rats = function(dtf, parms, distances){
  
  if(all(dtf$infection == 0) | all(dtf$infection != 0)) return(dtf)
  
  infection_mutation_rate = parms$infection_mutation_rate
  infection_mutation_amplitude = parms$infection_mutation_amplitude
  infection_radius = parms$infection_radius
  
  dtf %<>% 
    rowid_to_column(var = 'id')
  
  susceptibles.id = 
    dtf %>% 
    filter(infection == 0) %>% 
    pull(id)
  
  infected.id = setdiff(dtf$id, susceptibles.id)
  
  foo = 
    expand_grid(
      id = susceptibles.id, 
      source = infected.id
    ) %>%
    mutate(distance = as.numeric(distances[infected.id, susceptibles.id])) %>%
    filter(distance > 0 & distance < infection_radius) %>%
    mutate(
      infection.source = dtf$infection[source],
      virulence.source = dtf$virulence[source],
      trial = rbernoulli(n(), p = virulence.source),
      infection = trial * infection.source,
      virulence = trial * virulence.source
    )
    
  if(sum(foo$trial) == 0){
    dtf %>%
      select(-id) %>%
      return
  }
  
  foo %<>%
    filter(trial == TRUE) %>%
    group_by(id) %>%
    slice_sample(n = 1) %>%
    ungroup %>%
    mutate(
      speciate = rbernoulli(n(), p = infection_mutation_rate)
    )
  
  if(any(foo$speciate == TRUE)){
    bar = 
      foo %>%
      filter(speciate == TRUE) %>%
      mutate(
        infection = max(dtf$infection) + seq(n()),
        virulence = pmax(0, rnorm(n(), mean = virulence, sd = infection_mutation_amplitude))
      ) 
    
    matches = match(bar$id, foo$id)
    foo$infection[matches] = bar$infection
    foo$virulence[matches] = bar$virulence
  }
  
  matches = match(foo$id, dtf$id)
  dtf$infection[matches] = foo$infection
  dtf$virulence[matches] = foo$virulence
  
  dtf %>% 
    select(-id) %>%
    return
  
}

kill_rats = function(dtf, parms){
  
  turnover_factor = parms$turnover_factor
  mortality = parms$mortality
  density_weight = parms$density_weight
  J0 = parms$J0
  coef1 = parms$coef1
  coef2 = parms$coef2
  
  deaths = 
    dtf %>%
    slice_sample(
      n = turnover_factor * J0,
      weight_by = 
        mortality + 
        density_weight * local_density / max(local_density) + 
        coef1 * virulence + coef2 * virulence ^ 2
    )
  
  dtf %>%
    anti_join(deaths, by = 'individual') %>%
    return
}

reproduce_rats = function(dtf, parms){
  
  turnover_factor = parms$turnover_factor
  rat_mutation_rate = parms$rat_mutation_rate
  rat_mutation_amplitude = parms$rat_mutation_amplitude
  J0 = parms$J0
  
  parents = 
    dtf %>%
    slice_sample(
      n = turnover_factor * J0, 
      weight_by = NULL
    )
  
  pups = 
    parents %>%
    mutate(
      x = jitter(x),
      y = jitter(y),
      # x = runif(turnover_factor * J0, min = 0, max = city_size),
      # y = runif(turnover_factor * J0, min = 0, max = city_size),
      infection = 0,
      virulence = 0,
      individual = max(dtf$individual) + seq(n()),
      speciate = rbernoulli(n(), p = rat_mutation_rate)
    )
  
  if(any(pups$speciate == TRUE)){
    bar = 
      pups %>%
      filter(speciate == TRUE) %>%
      mutate(
        species = max(dtf$species) + seq(n()),
        mobility = pmax(0, rnorm(n(), mean = mobility, sd = rat_mutation_amplitude))
      ) 
    
    matches = match(bar$individual, pups$individual)
    pups$species[matches] = bar$species
    pups$mobility[matches] = bar$mobility
    
  }
  
 
  dtf %>%
    bind_rows(
      pups %>%
        select(-speciate)
    ) %>%
    return
  
}

local_density = function(dtf, parms){
  distance_matrix = as.matrix(dist_torus(dtf, parms))
  radius = parms$infection_radius
  
  dtf %>%
    mutate(
      local_density = rowSums(distance_matrix < radius)
    ) %>%
    return

}

# simulation = function(parms, do.save = seawulf){
# 
#   list2env(as.list(parms), env = environment())

simulation = function(
  
  scenario,
  J0,                              # size of rat community
  city_size,                       # size of city in meters
  maxtime,                         # simulation iterations
  infection_radius,    # radius of infection
  turnover_factor,                  # proportion of deaths/births per iteration
  rat_mutation_rate,               # probability of rat speciation per birth
  rat_mutation_amplitude,            # mutation step ~ N(0, rma)
  infection_mutation_rate,         # probability of disease speciation per infection
  infection_mutation_amplitude,    # mutation step ~ N(0, ima)
  mortality,                        # intrinsic rat mortality
  density_weight,              # contribution of local density to mortality
  coef1,                             # contribution of the 1st power of virulence to mortality
  coef2,                       # contribution of the 2nd power of virulence to mortality
  virulence0,                        # initial virulence
  mobility0,                        # initial mobility
  seed                            # seed generator
  
){
  
  parms = 
    tibble(
      scenario,
      J0,                              # size of rat community
      city_size,                       # size of city in meters
      maxtime,                         # simulation iterations
      infection_radius,    # radius of infection
      turnover_factor,                  # proportion of deaths/births per iteration
      rat_mutation_rate,               # probability of rat speciation per birth
      rat_mutation_amplitude,            # mutation step ~ N(0, rma)
      infection_mutation_rate,         # probability of disease speciation per infection
      infection_mutation_amplitude,    # mutation step ~ N(0, ima)
      mortality,                        # intrinsic rat mortality
      density_weight,              # contribution of local density to mortality
      coef1,                             # contribution of the 1st power of virulence to mortality
      coef2,                       # contribution of the 2nd power of virulence to mortality
      virulence0,                        # initial virulence
      mobility0,                        # initial mobility
      seed      
    )
  
  set.seed(seed)
  
  # initial conditions
  simtime = 0
  
  rats = 
    tibble(
      x = runif(J0, min = 0, max = city_size),
      y = runif(J0, min = 0, max = city_size),
      time = 0,
      species = 1,
      mobility = mobility0,
      infection = 0,
      virulence = 0
    ) %>%
    rowid_to_column(var = 'individual') %>%
    local_density(parms)
  
  initial_infections = sample(J0, size = turnover_factor * J0)
  
  rats$infection[initial_infections] = 1
  rats$virulence[initial_infections] = virulence0
  
  # dynamics loop
  while(simtime < maxtime){
    
    simtime = simtime + 1
    
    if(!seawulf & simtime %% 10 == 0) 
      print(simtime)
    
    current = 
      rats %>% 
      slice_max(time)
    
    distance_matrix = as.matrix(dist_torus(current, parms))
    
    update =
      current %>%
      mutate(time = simtime) %>%
      kill_rats(parms) %>%
      reproduce_rats(parms) %>%
      infect_rats(parms, distance_matrix) %>%
      move_rats(parms, distance_matrix)
    
    rats %<>% bind_rows(update)
    
    if(sum(update$infection > 0) == 0){
      writeLines('All infections cleared')
      break
    } 
    
  }
  
  # save data
  if(do.save){
    data = list(parms = parms, rats = rats)
    saveRDS(data, file = paste0('~/RatDisease/data/rats_scenario_',scenario,'.rds'))
  }
  
  return(rats)
  
}


## dynamics
if(do.dynamics){
  
  if(!do.parallel){
    rats = 
      parms[1, ] %>% 
      simulation
  }
  
  if(do.parallel){
    rats =
      parms %>%
      future_pmap(
        .f = simulation,
        .options = furrr_options(seed = NULL)
      )
  }
  
}

if(do.wrangling){
  setwd('~/RatDisease/data/')
  
  files = list.files(pattern = 'rats_scenario_*')
  
  dat =
    files %>%
    map_dfr(function(char){
      
      dtf = readRDS(char)
      parms = dtf$parms
      rats = dtf$rats
      
      foo =
        rats %>% 
        filter(time > 0 & time %% 100 == 0) %>% 
        group_by(time) %>% 
        summarize_at(c('mobility', 'virulence'), mean) %>%
        left_join(
          rats %>%
            filter(time > 0 & time %% 100 == 0) %>% 
            group_by(time) %>% 
            summarize(infection_load = sum(virulence > 0) / 100, .groups = 'drop'),
          by = 'time'
        ) %>%
        bind_cols(parms)
      
      bar = 
        rats %>% 
        filter(time > 0 & time %% 100 == 0) %>% 
        group_by(time) %>% 
        count(virulence) %>% 
        ungroup %>% 
        group_by(time) %>% 
        summarize(
          v.richness = length(n), 
          v.shannon = -1 / log(100) * sum(n / 100 * log(n / 100))) %>% 
        ungroup
      
      hop = 
        rats %>% 
        filter(time > 0 & time %% 100 == 0) %>% 
        group_by(time) %>% 
        count(mobility) %>% 
        ungroup %>% 
        group_by(time) %>% 
        summarize(
          m.richness = length(n), 
          m.shannon = -1 / log(100) * sum(n / 100 * log(n / 100))) %>% 
        ungroup
      
      barhop = left_join(bar, hop, by = 'time')
      
      return(left_join(foo, barhop, by = 'time'))
      
    }) 
  
  res = 
    dat %>%
    group_by(infection_radius, density_weight, coef2, time) %>%
    summarize_at(
      c('mobility', 'virulence', 'v.richness', 'v.shannon', 'm.richness', 'm.shannon', 'infection_load'), 
      list(mean = mean, sd = sd)
    ) %>%
    ungroup
  
  res # %<>% filter(infection_radius == 100)
  
  plot_virulence = 
    res %>%
    mutate(infection_radius = factor(infection_radius)) %>%
    ggplot(aes(time, virulence_mean, group = infection_radius, color = infection_radius)) +
    geom_line() +
    facet_grid(density_weight ~ coef2, labeller = label_both)
  
  plot_mobility = 
    res %>%
    mutate(infection_radius = factor(infection_radius)) %>%
    ggplot(aes(time, mobility_mean, group = infection_radius, color = infection_radius)) +
    geom_line() +
    facet_grid(density_weight ~ coef2, labeller = label_both)
  
  plot_richness = 
    res %>%
    mutate(infection_radius = factor(infection_radius)) %>%
    pivot_longer(c(v.richness_mean, m.richness_mean), names_to = 'measure') %>%
    ggplot(aes(time, value, group = measure, color = measure)) +
    geom_line() +
    facet_grid(density_weight ~ coef2, labeller = label_both)
  
  plot_shannon = 
    res %>%
    mutate(infection_radius = factor(infection_radius)) %>%
    pivot_longer(c(v.shannon_mean, m.shannon_mean), names_to = 'measure') %>%
    ggplot(aes(time, value, group = measure, color = measure)) +
    geom_line() +
    facet_grid(density_weight ~ coef2, labeller = label_both)
  
  plot_infection_load = 
    res %>% 
    ggplot(aes(time, infection_load_mean)) + 
    geom_line() + 
    facet_grid(density_weight ~ coef2, labeller = label_both)
  
  
  gridExtra::grid.arrange(
    plot_virulence, 
    plot_mobility, 
    plot_richness,
    plot_shannon,
    plot_infected_load
  )
  
}

## plots
if(do.plots){
  
  ## set aesthetics
  theme_set(theme_bw())
  theme_update(
    panel.grid = element_blank(), 
    aspect.ratio = 1
  )
  
  longlived = 
    rats %>% 
    count(individual) %>% 
    arrange(desc(n)) %>% 
    head(9) %>% 
    pull(individual)
  
  plot_movement = 
    rats %>% 
    filter(individual %in% longlived) %>% 
    ggplot(aes(x, y)) + 
    geom_path() + 
    geom_point() + 
    facet_wrap(~individual, scales = 'free')
  
  
  plot_density = 
    rats %>% 
    filter(time %% (maxtime / 10) == 0) %>% 
    mutate(infection = factor(infection)) %>%
    ggplot(aes(x, y, color = infection)) + 
    geom_point() + 
    facet_wrap(~time)
  
  rats_sieved = 
    rats %>%
    filter(time %% (maxtime/100) == 0)
  
  
  plot_infections = 
    rats_sieved %>%
    mutate(infection = factor(infection)) %>%
    group_by(time, infection) %>%
    count %>%
    ungroup %>%
    ggplot(aes(time, n, color = infection)) +
    geom_line() +
    theme(legend.position = 'none') +
    ggtitle('Infection strains')
  
  plot_mean_virulence = 
    rats_sieved %>% 
    filter(infection > 0) %>%
    group_by(time) %>% 
    summarize(mean_virulence = mean(virulence)) %>% 
    ggplot(aes(time, mean_virulence)) + 
    geom_line() +
    ggtitle('Mean virulence')
  
  plot_species =
    rats_sieved %>%
    mutate(species = factor(species)) %>%
    group_by(time, species) %>%
    count %>%
    ungroup %>%
    ggplot(aes(time, n, color = species)) +
    geom_line() +
    theme(legend.position = 'none') +
    ggtitle('Rat species')
  
  plot_mean_mobility = 
    rats_sieved %>% 
    group_by(time) %>% 
    summarize(mean_mobility = mean(mobility)) %>% 
    ggplot(aes(time, mean_mobility)) + 
    geom_line() +
    ggtitle('Mean mobility')
  
  plot_movement %>%
    show
  
  plot_density %>%
    show
  
  gridExtra::grid.arrange(
    plot_infections, 
    plot_mean_virulence,
    plot_species, 
    plot_mean_mobility
  )
  
}

