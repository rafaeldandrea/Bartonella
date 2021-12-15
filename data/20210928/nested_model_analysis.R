## Libraries

library(tidyverse)
library(pracma)
library(furrr)
library(parallel)
library(cowplot)
library(magrittr)
library(deSolve)

do.plots = 1

seawulf = as.logical(Sys.info()['user'] == 'rdrocha')
cores = if(seawulf) detectCores() else 4
plan(multisession, workers = cores)

filter = dplyr::filter

## Source functions
# source('https://github.com/dysordys/Krishnada_project/raw/main/code/nested_model_functions.R')
if(!seawulf) source('c:/users/rdand/Google Drive/GitHub/Krishnada_project/code/nested_model_functions.R')
if(seawulf) source('~/RatDisease/code/nested_model_functions.R')


## Plotting aesthetics

theme_set(theme_bw())
theme_update(
  panel.grid = element_blank(),
  aspect.ratio = 1
)

## Actions

do.equil_vs_dynamics = 0
do.evolution = 0
do.manual_evolution = 1

## Free parameters

parameters =
  expand_grid(
    number_of_species = 5,
    number_of_strains = 5,
    bmin = 1,
    bmax = 2 * bmin,
    beta0 = 5,
    omega0 = .5 * (bmax - bmin),
    # k = seq(1.2, 2, l = 10) / (bmax - bmin),
    k = 1.2 / (bmax - bmin),
    l = 0,
    lomega = c(.25 * 0:7, seq(2, 2.09, l = 10)),
    delta0 = 1,
    gamma0 = 1,
    sigma_delta = (bmax - bmin) * seq(.1, 2, l = 5),
    sigma_gamma = (bmax - bmin) *seq(.1, 10, l = 5),
    epsilon = .2,
    Kmin = 10
  ) %>%
  mutate(maxtime = ifelse(lomega < -.5, 20000, 1000)) %>%
  rowid_to_column(var = 'scenario')


## Tests and analyses

if(do.equil_vs_dynamics){
  ecoequil = 
    parameters %>%
    pmap_dfr(wrapper) %>%
    inner_join(parameters) %>%
    select(scenario, everything())
  
  ecoequil_summary =
    ecoequil %>%
    group_by(scenario) %>%
    summarize(
      minS = min(S),
      minIik = min(Iik),
      feasible = minS >=0 & minIik >= 0,
      .groups = 'drop'
    )
  
  scen = 
    ecoequil_summary %>%
    filter(feasible == TRUE) %>%
    pull(scenario) %>%
    pluck(1)
  
  parms_scenario = 
    parameters %>% 
    filter(scenario == scen)
  
  nspecies = 
    parms_scenario %>%
    pull(number_of_species)
  
  nstrains = 
    parms_scenario %>%
    pull(number_of_strains)
  
  Si_exact = 
    ecoequil %>%
    filter(scenario == scen) %>%
    select(species, S) %>%
    unique
  
  Iik_exact = 
    ecoequil %>%
    filter(scenario == scen) %>%
    select(species, strain, Iik)
  
  Ik_exact = 
    Iik_exact %>%
    group_by(strain) %>%
    summarize(
      Ik = sum(Iik), 
      .groups = 'drop'
    )
  
  set.seed(scen)
  
  ecosim = 
    ode(
      y = runif(nspecies * (1 + nstrains), min = 0, max = 2) * c(Si_exact$S, Iik_exact$Iik),
      times = seq(0, 1000, by = .1), 
      method = 'lsoda',
      func = eco_dynamics, 
      parms = parms_scenario
    ) %>%
    as.data.frame %>%
    as_tibble
  
  plot_Si_sim =
    ecosim[, 1:(nspecies + 1)] %>%
    filter(time > 0) %>%
    pivot_longer(-time) %>%
    drop_na(value) %>%
    mutate(species = factor(name, levels = seq(nspecies))) %>%
    ggplot(aes(time, value, group = species, color = species)) +
    geom_hline(yintercept = Si_exact$S, color = 'grey') +
    geom_line() +
    scale_x_log10() +
    ylab('Susceptibles') +
    ggtitle('Susceptibles')
  
  
  plot_Ik_sim =
    ecosim[, c(1, seq(2 + nspecies, ncol(ecosim)))] %>%
    as.data.frame %>%
    as_tibble %>%
    filter(time > 0) %>%
    pivot_longer(-time) %>%
    mutate(name = as.numeric(name) - nspecies) %>%
    mutate(
      species = arrayInd(name, .dim = c(nspecies, nstrains))[, 1],
      strain = arrayInd(name, .dim = c(nspecies, nstrains))[, 2]
    ) %>%
    group_by(time, strain) %>%
    summarize(
      Ik = sum(value), 
      .groups = 'drop'
    ) %>%
    mutate(strain = factor(strain, levels = seq(nstrains))) %>%
    ggplot(aes(time, Ik, group = strain, color = strain)) +
    geom_hline(yintercept = Ik_exact$Ik, color = 'grey') +
    geom_line() +
    scale_x_log10() +
    ylab('Total infections') +
    ggtitle('Infections')
  
  plot = 
    plot_grid(
      plot_Si_sim,
      plot_Ik_sim,
      nrow = 1,
      rel_heights = c(1, .9)
    )
  
  plot %>%
    show
}


if(do.evolution){
  evosim = 
    ode(
      y = dependent_parms(parms_scenario)$omega, 
      times = seq(0, 1, by = .0001), 
      method = 'lsoda',
      func = evo_dynamics, 
      parms = parms_scenario
    ) %>%
    as.data.frame %>%
    as_tibble
  
  evosim_manual = 
    evo_dynamics_manual(
      t = seq(0, 2, by = .0001),
      y = dependent_parms(parms_scenario)$omega,
      parms = parms_scenario,
      omegamin = .04
    )
  
  evosim_manual %>% 
    filter(time > 0) %>% 
    mutate(strain = factor(strain, levels = 1:nstrains)) %>% 
    ggplot(aes(time, omega, group = strain, color = strain)) + 
    geom_line() + 
    scale_x_log10()
  
}

if(do.manual_evolution){
  
  run =
    function(scenario, parms){
    
      parms_scenario = parms[scenario, ]
      
      result = 
        parms_scenario %>%
        bind_cols(
          evo_dynamics_fixedxbar_manual(
            t = seq(0, parms_scenario$maxtime, by = .01),
            y = 
              c(
                eco_equilibrium(parms_scenario) %>% 
                  select(species, S) %>% 
                  unique %>% 
                  pull(S),
                eco_equilibrium(parms_scenario) %>% 
                  select(species, strain, Iik) %>% 
                  pull(Iik)
              ),
            parms = parms_scenario,
            evolve_interval = 10,
            verbose = !seawulf
          ) 
        ) %>%
        mutate(
          species = factor(species, levels = seq(parms_scenario$number_of_species)),
          strain = factor(strain, levels = seq(parms_scenario$number_of_strains))
        )
      
      return(result)
      
    }
  
  indices = seq(nrow(parameters))
  
  result = 
    indices %>%
    future_map_dfr(
      .f = run,
      parms = parameters,
      .options = furrr_options(seed = TRUE)
    )
  
  if(!seawulf) filename = 'c:/users/rdand/Google Drive/GitHub/Krishnada_project/data/20210913/result.rds'
  if(seawulf) filename = '~/RatDisease/data/parameter_exploration.rds'
  
  if(seawulf) saveRDS(result, file = filename)
  
  if(!seawulf & do.plots){
    
    
    ## Time series plots
    
    dtf = 
      result %>%
      filter(time >= 100, lomega <= 2.04)
    
    dtf_host = 
      dtf %>%
      group_by(scenario, lomega, time, species) %>%
      summarize(
        N = S + sum(Iik),
        .groups = 'drop'
      ) 
    
    dtf_pathogen = 
      dtf %>%
      group_by(scenario, lomega, time, strain) %>%
      summarize(
        omega = unique(omega),
        Ik = sum(Iik),
        .groups = 'drop'
      )
    
    plot_abundance = 
      dtf_host %>%
      ggplot(aes(1 + time, N, group = species, color = species)) +
      geom_line() +
      # geom_point() +
      # scale_x_log10() +
      facet_wrap(~ lomega, scales = 'free')
    
    plot_infection = 
      dtf_pathogen %>%
      ggplot(aes(1 + time, Ik, group = strain, color = strain)) +
      geom_line() +
      # geom_point() +
      # scale_x_log10() +
      facet_wrap(~ lomega, scales = 'free')
    
    plot_omega = 
      dtf_pathogen %>%
      ggplot(aes(1 + time, omega, group = strain, color = strain)) +
      geom_line() +
      # geom_point() +
      # scale_x_log10() +
      facet_wrap(~ lomega, scales = 'free')
    
    gridExtra::grid.arrange(plot_abundance, plot_infection, plot_omega, nrow = 1)
    
    
    ## Equilibrium plots
    
    plot_omega_by_lomega = 
      result %>%
      group_by(scenario) %>%
      slice_max(time) %>%
      select(scenario, strain, omega, lomega, sigma_gamma, sigma_delta) %>%
      unique() %>%
      ungroup() %>%
      ggplot(aes(lomega, omega, group = strain, color = strain)) +
      geom_line() +
      geom_point() +
      labs(
        x = 'Penalty for niche width',
        y = 'Evolved niche width'
      )
    
    plot_infections_by_lomega = 
      result %>%
      group_by(scenario) %>%
      slice_max(time) %>%
      ungroup() %>%
      group_by(scenario, lomega, sigma_gamma, sigma_delta, strain) %>%
      summarize(Ik = sum(Iik), omega = unique(omega), .groups = 'drop') %>%
      ggplot(aes(lomega, Ik, group = strain, color = strain)) +
      geom_line() + 
      geom_point()
    
    plot_infections_by_omega = 
      result %>%
      group_by(scenario) %>%
      slice_max(time) %>%
      ungroup() %>%
      group_by(scenario, lomega, sigma_gamma, sigma_delta, strain) %>%
      summarize(Ik = sum(Iik), omega = unique(omega), .groups = 'drop') %>%
      ggplot(aes(omega, Ik, color = strain)) +
      geom_point() +
      facet_wrap(~lomega, scales = 'free')
    
    plot_abundance_by_species = 
      result %>%
      group_by(scenario) %>%
      slice_max(time) %>%
      ungroup() %>%
      group_by(scenario, lomega, sigma_gamma, sigma_delta, species) %>%
      summarize(
        Ii = sum(Iik), 
        Si = unique(S), 
        N = sum(Si + Ii), 
        .groups = 'drop'
      ) %>%
      ggplot(aes(species, N)) +
      geom_point() +
      facet_wrap(~lomega, scales = 'free')
    
    gridExtra::grid.arrange(
      plot_omega_by_lomega,
      plot_infections_by_lomega,
      plot_infections_by_omega,
      plot_abundance_by_species
    )
    
  }
  
  
    
}
