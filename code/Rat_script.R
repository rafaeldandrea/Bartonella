## Model on infectious diseases on a lattice with diversification of hosts and pathogens

## Determine whether working on SeaWulf (SBU hpc) or personal computer
seawulf = as.logical(Sys.info()['user'] == 'rdrocha')

## Read array argument from command line
index = as.numeric(commandArgs(TRUE))[1]
if(is.na(index)) index = 1


library(tidyverse)
library(magrittr)

setwd('~/Krishnada_project/data/')

theme_set(theme_bw())
theme_update(
  panel.grid = element_blank(), 
  aspect.ratio = 1
)

do.dynamics = TRUE
do.save = FALSE
do.plots = TRUE

parameters = 
  expand_grid(
    J0 = 1e2,
    city_size = 1e3,
    maxtime = 1e4,
    infection_radius = 100,
    turnover_factor = .1,
    mortality = .1,
    rat_mutation_rate = .01,
    rat_mutation_amplitude = 5,
    infection_mutation_rate = .01,
    infection_mutation_amplitude = .05,
    virulence0 = .5,
    mobility0 = 10,
    seeder = 1:100
  )

parms = parameters[index, ]

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
  
  deaths = 
    dtf %>%
    slice_sample(
      n = turnover_factor * J0,
      weight_by = 1 * local_density / max(local_density) + mortality + 1 * virulence
    )
  
  dtf %>%
    anti_join(deaths, by = 'individual') %>%
    return
}


reproduce_rats = function(dtf, parms){
  
  turnover_factor = parms$turnover_factor
  rat_mutation_rate = parms$rat_mutation_rate
  rat_mutation_amplitude = parms$rat_mutation_amplitude
  
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


set.seed(parms$seeder)

## dynamics
if(do.dynamics){
  
  J0 = parms$J0
  maxtime = parms$maxtime
  seeder = parms$seeder
  
  rats = 
    tibble(
      x = runif(J0, min = 0, max = parms$city_size),
      y = runif(J0, min = 0, max = parms$city_size),
      time = 0,
      species = 1,
      mobility = parms$mobility0,
      infection = 0,
      virulence = 0
    ) %>%
    rowid_to_column(var = 'individual') %>%
    local_density(parms)
  
  initial_infections = sample(J0, size = parms$turnover_factor * J0)
  
  rats$infection[initial_infections] = 1
  rats$virulence[initial_infections] = parms$virulence0
  
  
  simtime = 0
  
  while(simtime < maxtime){
    
    simtime = simtime + 1
    if(simtime %% 10 == 0) print(simtime)
    
    current = 
      rats %>% 
      slice_max(time)
    
    distance_matrix = as.matrix(dist_torus(current, parms))
    seed = seeder + simtime
    
    update =
      current %>%
      mutate(time = simtime) %>%
      kill_rats(parms) %>%
      reproduce_rats(parms) %>%
      infect_rats(parms, distance_matrix) %>%
      move_rats(parms, distance_matrix)
    
    
    if(sum(update$infection > 0) == 0){
      writeLines('All infections cleared')
      break
    } 
    
    rats %<>% bind_rows(update)
  }
  
}

## save data
if(do.save){
  data = list(parameters = as_tibble(parameters), rats = rats)
  save(data, file = 'rats_data.RData')
}

## plots
if(do.plots){
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

