## Rodents vs pathogens model -- nested infectivity


## Libraries

library(tidyverse)
library(pracma)
library(furrr)
library(parallel)
library(cowplot)
library(magrittr)
library(deSolve)

do.explore = 1
do.analysis = 1
do.plot = 0


## Plotting aesthetics

theme_set(theme_bw())
theme_update(
  panel.grid = element_blank(),
  aspect.ratio = 1
)

## Free parameters

parameters =
  expand_grid(
    number_of_species = 1,
    number_of_strains = 1,
    bmin = .01,
    bmax = 1,
    beta0 = .5,
    omega0 = .5,
    k = seq(1.2, 2, l = 20),
    l = 0,
    delta0 = 1,
    gamma0 = 1,
    sigma_delta = seq(.01, .2, l = 20),
    sigma_gamma = seq(.01, 1, l = 20),
    epsilon = .05
  ) %>%
  rowid_to_column(var = 'scenario')

cutoff <- function(n) {
  return(ifelse(n<1, (1*(n>0))*(n*n*n*(10+n*(-15+6*n))), 1))
}

dynamics = 
  function(t, y, parms){
    
    list2env(parms, envir = environment())
    
    # y = pmax(0, y)
    
    S = y[1:number_of_species]
    
    Iik = 
      matrix(
        y[(number_of_species + 1): length(y)],
        number_of_species, 
        number_of_strains
      )
    
    Ik = colSums(Iik)
    SiIk = outer(S, Ik)
    
    dSdt = 
      (fecundity - mortality) * S + rowSums((fecundity + gamma) * Iik) -
      rowSums(beta * SiIk)
    
    dIdt = beta * SiIk - (mortality * (1 + delta) + gamma) * Iik
    
    return(list(c(as.numeric(dSdt), as.numeric(dIdt))))
  }

        
omega_dynamics = 
  function(t, y, parms){
    
    list2env(parms, envir = environment())
    
    fecundity = seq(bmin, bmax, l = number_of_species)
    mortality = fecundity * (1 - epsilon)
    
    omega = pmax(0, y)
    xbar = bmin - 1 / k * log(omega / omega0)
    xminusxbar = outer(fecundity, xbar, '-')
    omegainv = matrix(1 / omega, number_of_species, number_of_strains, byrow = TRUE)
    
    beta = beta0 * exp(l * fecundity) * exp(-xminusxbar ^ 2 * omegainv ^ 2)
    delta = delta0 * exp(-(xminusxbar / sigma_delta) ^ 2)
    gamma = gamma0 * exp(-(xminusxbar / sigma_gamma) ^ 2)
    
    zeta = 1 / (mortality * (1 + delta) + gamma)
    
    Q = t(zeta * beta)
    
    M = (1 - (fecundity + gamma) * zeta) * beta
    
    Si = pmax(0, solve(Q, rep(1, number_of_species)))
    
    Ik = pmax(0, solve(M, fecundity - mortality))
    
    Iik = zeta * beta * outer(Si, Ik)
    
    xbar_prime = 
      matrix(
        - 1 / k / omega, 
        number_of_species, 
        number_of_strains, 
        byrow = TRUE
      )
    
    beta_prime = 
      beta * 2 * xminusxbar * omegainv ^ 2 * (xbar_prime + xminusxbar * omegainv)
    
    zeta_prime = -zeta ^ 2 * 2 * xminusxbar * xbar_prime *
      (mortality * delta / sigma_delta ^ 2 + gamma / sigma_gamma ^ 2)
    
    domega_dt = 
      colSums(beta_prime * Si) +
      colSums(1 / zeta ^ 2 * zeta_prime * Iik %*% diag(1 / (1e-16 + Ik)))
    
    return(list(as.numeric(domega_dt)))
  }

omega_equil = 
  function(omega, parms){
    
    list2env(parms, envir = environment())
    
    fecundity = seq(bmin, bmax, l = number_of_species)
    
    mortality = fecundity * (1 - epsilon)
    
    xbar = bmin - 1 / k * log(omega / omega0)
    xminusxbar = outer(fecundity, xbar, '-')
    omegainv = matrix(1 / omega, number_of_species, number_of_strains, byrow = TRUE)
    
    beta = beta0 * exp(l * fecundity) * exp(-xminusxbar ^ 2 * omegainv ^ 2)
    
    delta = delta0 * exp(-(xminusxbar / sigma_delta) ^ 2)
    
    gamma = gamma0 * exp(-(xminusxbar / sigma_gamma) ^ 2)
    
    zeta = 1 / (mortality * (1 + delta) + gamma)
    
    Q = t(zeta * beta)
    
    M = (1 - (fecundity + gamma) * zeta) * beta
    
    # Si = try(solve(Q, rep(1, number_of_species)), silent = TRUE)
    Si = pmax(0, solve(Q, rep(1, number_of_species)))
    
    # if(class(Si) == 'try-error') return(NULL)
    
    # Ik = try(solve(M, fecundity - mortality), silent = TRUE)
    Ik = pmax(0, solve(M, fecundity - mortality))
    
    # if(class(Ik) == 'try-error') return(NULL)
    
    Iik = zeta * beta * outer(Si, Ik)
    
    xbar_prime = 
      matrix(
        - 1 / k / omega, 
        number_of_species, 
        number_of_strains, 
        byrow = TRUE
      )
    
    beta_prime = 
      beta * 2 * xminusxbar * omegainv ^ 2 * (xbar_prime + xminusxbar * omegainv)
    
    zeta_prime = -zeta ^ 2 * 2 * xminusxbar * xbar_prime *
      (mortality * delta / sigma_delta ^ 2 + gamma / sigma_gamma ^ 2)
    
    domega_dt = 
      colSums(beta_prime * Si) +
      colSums(1 / zeta ^ 2 * zeta_prime * Iik %*% diag(1 / (Ik + 1e-16)))
    
    return(domega_dt)
  }

# parms = 
#   parameters %>% 
#   filter(scenario == scen) %>%
#   as.list
#   
# sim = 
#   multiroot(
#     f = omega_equil, 
#     start = with(parms, exp(-k * seq(bmin, bmax, l = number_of_strains))),
#     positive = TRUE,
#     parms = parms
#   )

model = 
  function(
    number_of_species,
    number_of_strains,
    bmin,
    bmax,
    beta0,
    omega0,
    k,
    l,
    delta0,
    gamma0,
    sigma_delta,
    sigma_gamma,
    epsilon,
    ...
  ){
    
    ## Dependent parameters
    
    fecundity = seq(bmin, bmax, l = number_of_species)
    
    mortality = fecundity * (1 - epsilon)
    
    xbar = fecundity
    omega = omega0 * exp(-k * (xbar - bmin))

    xminusxbar = outer(fecundity, xbar, '-')
    
    omegainv = matrix(1 / omega, number_of_species, number_of_strains, byrow = TRUE)
    
    R = beta0 * exp(l * fecundity)
    E = exp(-xminusxbar ^ 2 * omegainv ^ 2)
    beta = E * R
    
    delta = delta0 * exp(-(xminusxbar / sigma_delta) ^ 2)
    
    gamma = gamma0 * exp(-(xminusxbar / sigma_gamma) ^ 2)
    
    ## Model Equilibrium
    
    zeta = 1 / (mortality * (1 + delta) + gamma)
    
    Q = t(zeta * beta)
    
    M = (1 - (fecundity + gamma) * zeta) * beta
    
    # Si = try(solve(Q, rep(1, number_of_species)), silent = TRUE)
    Si = solve(Q, rep(1, number_of_species))
    
    if(class(Si) == 'try-error') return(NULL)
    
    # Ik = try(solve(M, fecundity - mortality), silent = TRUE)
    Ik = solve(M, fecundity - mortality)
    
    if(class(Ik) == 'try-error') return(NULL)
    
    Iik = zeta * beta * outer(Si, Ik)
    
    xbar_prime = 
      matrix(
        - 1 / k / omega, 
        number_of_species, 
        number_of_strains, 
        byrow = TRUE
      )
    
    beta_prime = 
      beta * 2 * xminusxbar * omegainv ^ 2 * (xbar_prime + xminusxbar * omegainv)
    
    zeta_prime = -zeta ^ 2 * 2 * xminusxbar * xbar_prime *
      (mortality * delta / sigma_delta ^ 2 + gamma / sigma_gamma ^ 2)
    
    domega_dt = 
      colSums(beta_prime * Si) +
      colSums(
        1 / zeta ^ 2 * zeta_prime * 
          Iik %*% diag(1 / (1e-16 + Ik), nrow = length(Ik), ncol = length(Ik))
      )
    
    
    parms_tbl = 
      tibble(
        number_of_species,
        number_of_strains,
        bmin,
        bmax,
        beta0,
        omega0,
        k,
        l,
        delta0,
        gamma0,
        sigma_delta,
        sigma_gamma,
        epsilon
      )
    
    host_tbl = 
      tibble(
        species = seq(number_of_species),
        fecundity = fecundity,
        mortality = mortality,
        S = Si
      )
    
    pathogen_tbl = 
      tibble(
        strain = seq(number_of_strains),
        xbar = xbar,
        omega = omega,
        Ik = Ik,
        selection_on_niche_width = domega_dt
      )
      
    parms_tbl %>%
      bind_cols(host_tbl) %>%
      expand_grid(pathogen_tbl) %>%
      arrange(strain, species) %>%
      mutate(
        beta = as.numeric(beta),
        delta = as.numeric(delta),
        gamma = as.numeric(gamma),
        Iik = as.numeric(Iik)
      ) %>%
      return()
    
}

## Parameter exploration
if(do.explore){
  plan(multisession, workers = detectCores() - 1)
  
  results_raw = 
    parameters %>%
    future_pmap_dfr(.f = model)
  
  results = 
    results_raw %>%
    inner_join(parameters) %>%
    select(scenario, everything()) %>%
    group_by(scenario) %>%
    mutate(
      maxSel = max(selection_on_niche_width),
      minSel = min(selection_on_niche_width),
      minS = min(S),
      minI = min(Iik),
      minSI = minS * minI,
      feasible = minS >= 0 & minI >= 0
    ) %>%
    ungroup
}
  
plot_kl =
  results %>%
  group_by(k, l) %>%
  summarize(
    mean = mean(feasible), 
    .groups = 'drop'
  ) %>%
  ggplot(aes(k, l, fill = mean)) +
  geom_tile() +
  scale_fill_gradientn(colors = terrain.colors(100))

plot_gamma_delta =
  results %>%
  group_by(sigma_delta, sigma_gamma) %>%
  summarize(
    mean = mean(feasible), 
    .groups = 'drop'
  ) %>%
  ggplot(aes(sigma_delta, sigma_gamma, fill = mean)) +
  geom_tile() +
  scale_fill_gradientn(colors = terrain.colors(100))

gridExtra::grid.arrange(plot_kl, plot_gamma_delta)

if(do.analysis){
    
    scens = 
      results %>%
      filter(feasible == TRUE) %>%
      pull(scenario) %>%
      unique()
      
    scen = scens[2]
    
    parms = 
      with(
        results %>%
          filter(scenario == scen),
        list(
          number_of_species = unique(number_of_species),
          number_of_strains = unique(number_of_strains),
          fecundity = sort(unique(fecundity)),
          mortality = sort(unique(mortality)),
          beta = matrix(beta, number_of_species, number_of_strains),
          delta = matrix(delta, number_of_species, number_of_strains),
          gamma = matrix(gamma, number_of_species, number_of_strains)
        )
      )
    
    nspecies = unique(results$number_of_species)
    nstrains = unique(results$number_of_strains)
    y = c(rep(.9, nspecies), rep(.1, nspecies * nstrains))
    
    thetimes = exp(seq(0, log(2500), by = .001))
    thetimes = seq(0, 2500, by = .1)
    r = 
      ode(
        y = y, 
        times = thetimes, 
        method = 'vode',
        func = dynamics, 
        parms = parms,
        rtol = 1e-12,
        atol = 1e-12
      )
    
    S = r[nrow(r), 1 + seq(nspecies)]
    Iik = matrix(r[nrow(r), (2 + nspecies):ncol(r)], nspecies, nstrains)
    
    Sexact = 
      results %>%
      filter(scenario == scen) %>%
      pull(S) %>%
      unique()
    
    plot(parms$fecundity, Sexact, t = 'h', ylim = c(0, max(Sexact, S)), las = 1)
    points(parms$fecundity, S, pch = 20, col = 'red')
    
    r[, c(1, seq(2 + nspecies, ncol(r)))] %>%
      as.data.frame %>%
      as_tibble %>%
      pivot_longer(-time) %>%
      mutate(name = as.numeric(name) - nspecies) %>%
      mutate(
        species = arrayInd(name, .dim = c(nspecies, nstrains))[, 1],
        strain = arrayInd(name, .dim = c(nspecies, nstrains))[, 2]
      ) %>%
      select(-name) %>%
      group_by(time, species) %>%
      summarize(Ii = sum(value), .groups = 'drop') %>%
      mutate(species = as.character(species)) %>%
      inner_join(
        r[, seq(1 + nspecies)] %>%
          as.data.frame %>%
          as_tibble %>%
          # filter(time < 5e3) %>%
          pivot_longer(-time) %>%
          rename(species = name, Si = value)
      ) %>%
      mutate(abundance = Si + Ii) %>%
      ggplot(aes(time, abundance, color = species)) +
      geom_line()
    
    
    
    
    plot_time = 
      r[, seq(1 + nspecies)] %>%
      as.data.frame %>%
      as_tibble %>%
      # filter(time < 5e3) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, color = name, group = name)) +
      geom_line()
    
    plot_time +
      scale_x_log10() %>% 
      show
    
    y = 
      with(
        parameters %>% 
          filter(scenario == scen), 
        omega0 * exp(-k * (seq(bmin, bmax, l = number_of_species) - bmin))
      )
    w = 
      ode(
        y = y, 
        times = thetimes, 
        func = omega_dynamics, 
        parms = 
          parameters %>% 
          filter(scenario == scen) %>%
          as.list
      )
    
    results %>%
      filter(feasible == TRUE) %>%
      group_by(k, sigma_delta, sigma_gamma) %>%
      slice_max(omega, n = 1) %>%
      ungroup %>%
      ggplot(aes(k, selection_on_niche_width)) +
      geom_point() +
      facet_grid(sigma_delta ~ sigma_gamma)
    
    results %>%
      filter(feasible == TRUE) %>%
      group_by(k, sigma_delta, sigma_gamma) %>%
      slice_max(omega, n = 1) %>%
      ungroup %>%
      group_by(sigma_delta, sigma_gamma) %>%
      slice_max(selection_on_niche_width) %>%
      ungroup %>%
      select(sigma_delta, sigma_gamma, selection_on_niche_width) %>%
      unique() %>%
      ggplot(aes(sigma_delta, sigma_gamma, fill = selection_on_niche_width)) +
      geom_tile() +
      scale_fill_gradient(low = 'blue', high = 'red') +
      ggtitle('Selection on strain with largest niche width') +
      labs(
        fill = 'selection', 
        x = expression(paste('width of virulence curve (', sigma[delta], ')')), 
        y = expression(paste('width of recovery curve (', sigma[gamma], ')'))
      )
    
    results |> 
      filter(feasible == TRUE) |> 
      select(k, sigma_delta, sigma_gamma, feasible) |> 
      unique() |> 
      count(sigma_gamma, sigma_delta) |> 
      ggplot(aes(sigma_delta, sigma_gamma, fill = n)) + 
      geom_tile() +
      scale_fill_gradient(low = 'blue', high = 'red') +
      ggtitle('Number of feasible scenarios (by steepness of niche width gradient)') +
      labs(
        fill = 'selection', 
        x = expression(paste('width of virulence curve (', sigma[delta], ')')), 
        y = expression(paste('width of recovery curve (', sigma[gamma], ')'))
      )
    
    results |> 
      filter(feasible == TRUE, k == max(k)) |> 
      select(sigma_delta, sigma_gamma, strain, omega, k, Ik) |> 
      unique() |> 
      mutate(
        s_delta = round(sigma_delta, 2),
        s_gamma = round(sigma_gamma, 2)
      ) %>%
      ggplot(aes(strain, Ik)) + 
      geom_line() + 
      facet_grid(s_delta ~ s_gamma, labeller = label_both)
    
    results |> 
      filter(feasible == TRUE) |> 
      select(
        sigma_delta, 
        sigma_gamma, 
        k, 
        l,
        omega0,
        strain, 
        omega, 
        selection_on_niche_width
      ) |> 
      unique() |> 
      mutate(
        s_delta = round(sigma_delta, 2),
        s_gamma = round(sigma_gamma, 2)
      ) %>%
      mutate(k = factor(k)) %>%
      ggplot(aes(strain, selection_on_niche_width, group = k, color = k)) + 
      geom_hline(yintercept = 0, color = 'grey') +
      geom_line() + 
      facet_grid(s_delta ~ l + s_gamma, labeller = label_both)
    
    
    plot_kl =
      results %>%
      group_by(k, l) %>%
      summarize(
        mean = mean(feasible), 
        .groups = 'drop'
      ) %>%
      ggplot(aes(k, l, fill = mean)) +
      geom_tile() +
      scale_fill_gradientn(colors = terrain.colors(100))
    
    plot_gamma_delta =
      results %>%
      group_by(sigma_delta, sigma_gamma) %>%
      summarize(
        mean = mean(feasible), 
        .groups = 'drop'
      ) %>%
      ggplot(aes(sigma_delta, sigma_gamma, fill = mean)) +
      geom_tile() +
      scale_fill_gradientn(colors = terrain.colors(100))
    
    results_feasible =
      results %>%
      filter(feasible == TRUE) %>%
      select(scenario, minS, minI, minSI, maxSel, minSel) %>%
      unique()
    
    selected_scenarios = 
      with(
        results_feasible,
        tibble(
          maxS = scenario[which.max(minS)],
          maxI = scenario[which.max(minI)],
          maxSI = scenario[which.max(minSI)],
          maxmaxSel = scenario[which.max(maxSel)],
          maxminSel = scenario[which.max(minSel)],
          minmaxSel = scenario[which.min(maxSel)],
          minminSel = scenario[which.min(minSel)]
        )
      )

}

## Plots
if(do.plot){
  
  # index = 426 ## maximizes S: diseases are all specialized -- most infect a single host, 
                ## diseases have low and equitable presence
  
  # index = 110 ## maximizes Ik: no niche breadth gradient
                ## diseases have very high presence (Si make up 10% of community)
  
  # index = 810 ## maximizes Ik under constraint of negative niche breadth gradient
                ## rodent species abundances do not seem to follow an obvious pattern
  
  # index = 9410
               
  index = 310
  mod = 
    parameters[index, ] %>%
    pmap_dfr(.f = model) %>%
    mutate(
      species = factor(species),
      strain = factor(strain),
      selection = 
        factor(
          selection_on_niche_width > 0, 
          levels = c(TRUE, FALSE), 
          labels = c('positive', 'negative')
        ),
      sel = selection_on_niche_width / 
        sum(abs(selection_on_niche_width)) * 
        number_of_species * max(omega)
    )
  
  stopifnot(
    mod %>% pull(S) %>% min() >= 0,
    mod %>% pull(Ik) %>% min() >= 0
  )
  
  summary_host = 
    mod %>% 
    group_by(species) %>%
    summarize(
      fecundity = unique(fecundity),
      mortality = unique(mortality),
      `b - d` = fecundity - mortality,
      infection = sum(beta), 
      virulence = sum(delta),
      recovery  = sum(gamma),
      .groups = 'drop'
    ) %>%
    mutate(species = as.numeric(species))
  
  summary_strain = 
    mod %>% 
    group_by(strain) %>%
    summarize(
      omega = unique(omega),
      selection = unique(selection_on_niche_width),
      infection = sum(beta), 
      virulence = sum(delta),
      recovery  = sum(gamma),
      .groups = 'drop'
    ) %>%
    mutate(strain = as.numeric(strain))
  
  rescale = 
    function(x){
      if(sd(x) == 0) return(rep(0, length(x)))
      return(scale(x))
    }
  
  plot_summary_host = 
    summary_host %>% 
    select(-c(fecundity, mortality)) %>%
    pivot_longer(-species) %>%
    group_by(name) %>%
    mutate(value = rescale(value)) %>%
    ungroup %>%
    ggplot(aes(species, value, group = name, color = name)) + 
    geom_line() + 
    scale_x_continuous(breaks = 1:10) +
    labs(x = 'host species', y = 'scaled value') +
    ggtitle('Host statistics')
  
  plot_summary_strain = 
    summary_strain %>% 
    pivot_longer(-strain) %>%
    group_by(name) %>%
    mutate(value = rescale(value)) %>%
    ungroup %>%
    ggplot(aes(strain, value, group = name, color = name)) + 
    geom_line() + 
    scale_x_continuous(breaks = 1:10) +
    labs(x = 'strain', y = 'scaled value') +
    ggtitle('Strain statistics')
  
  plot_species = 
    mod %>%
    group_by(species) %>%
    summarize(
      susceptible = S,
      infected = sum(Iik), 
      .groups = 'drop'
    ) %>%
    unique() %>%
    ggplot(aes(x = species,  y = rep(0, length(species)))) +
    geom_segment(
      aes(xend = species, yend = susceptible),
      size = 1.5
    ) +
    geom_segment(
      aes(xend = species, y = susceptible, yend = susceptible + infected),
      color = 'red',
      size = 1.5
    ) +
    labs(
      x = 'host species', 
      y = 'abundance'
    ) +
    ggtitle('Host abundance: S + I')
    
  
  plot_Ik = 
    mod %>% 
    select(strain, Ik) %>%
    unique() %>%
    ggplot(aes(strain, Ik)) + 
    geom_segment(aes(xend = strain, yend = rep(0, length(strain)))) + 
    geom_point() +
    labs(y = 'total infections') +
    ggtitle(expression(paste('Infection abundance: ',Sigma[i], I[ik])))
  
  plot_niche_width = 
    mod %>%
    ggplot(aes(strain, omega)) +
    geom_segment(
      aes(
        xend = strain, 
        yend = omega + sel,
        group = selection,
        color = selection
      ),
      arrow = arrow(type = 'closed', length = unit(.1, 'inches'))
    ) + 
    geom_point() +
    labs(
      x = 'strain',
      y = 'niche width (omega)'   
    ) +
    theme(legend.position = 'none') +
    ggtitle(expression(paste('Selection vs ', omega)))
  
  plot_infectivity = 
    mod %>% 
    ggplot(aes(species, beta, group = strain, color = strain)) +
    geom_line() +
    labs(
      x = 'host species',
      y = 'infectivity rate (beta)'
    ) +
    theme(legend.position = 'none') +
    ggtitle('Infectivity rates')
  
  plot_virulence = 
    mod %>%
    ggplot(aes(species, delta)) +
    geom_line(aes(group = strain, color = strain)) +
    geom_point(aes(group = strain, color = strain)) +
    geom_line(
      aes(species, virulence), 
      data = summary_host,
      color = 'grey',
      size = 2
    ) +
    geom_line(
      aes(strain, virulence), 
      data = summary_strain,
      color = rgb(153 / 255, 0, 0),
      size = 2
    ) +
    labs(
      x = 'host species',
      y = 'mortality rate (delta)'
    ) +
    ggtitle('Virulence (disease-induced mortality) -- Total by host (grey) and by strain (red)')
  
  plot_recovery = 
    mod %>%
    ggplot(aes(species, gamma)) +
    geom_line(aes(group = strain, color = strain)) +
    geom_point(aes(group = strain, color = strain)) +
    geom_line(
      aes(species, recovery), 
      data = summary_host,
      color = 'grey',
      size = 2
    ) +
    geom_line(
      aes(strain, recovery), 
      data = summary_strain,
      color = rgb(153 / 255, 0, 0),
      size = 2
    ) +
    labs(
      x = 'host species',
      y = 'recovery rate (gamma)'
    ) +
    ggtitle('Recovery rates -- Total by host (grey) and by strain (red)')
  
  gridExtra::grid.arrange(
    plot_species,
    plot_Ik,
    plot_summary_host,
    plot_summary_strain,
    plot_infectivity,
    plot_niche_width
  )
  
}


## ================== Old Code ================
# plot_bars =
#   results |> 
#   mutate(
#     feasible = as.factor((Smin > 0) * (Iikmin > 0))
#   ) |> 
#   group_by(feasible) |> 
#   summarize_all(mean) |> 
#   select(-c(2:12)) |> 
#   pivot_longer(-feasible, names_to = 'index') |> 
#   ggplot(aes(index, value, group = feasible, fill = feasible)) + 
#   geom_col(position = 'dodge', color = 'black') +
#   theme(aspect.ratio = .5)
# 
# plot_bars2 = 
#   results |>
#   mutate(
#     feasible = as.factor((Smin > 0) * (Iikmin > 0))
#   ) |>
#   filter(feasible == 1) |>
#   select(-c(1:11)) |>
#   pivot_longer(-feasible, names_to = 'index') |>
#   group_by(index, value) |>
#   count(feasible) |>
#   ungroup() |>
#   mutate(value = factor(value)) |>
#   ggplot(aes(index, n, group = value, fill = value)) +
#   geom_col(position = 'dodge', color = 'black') +
#   theme(aspect.ratio = .5)
# 
# gridExtra::grid.arrange(plot_bars, plot_bars2)
# 
# 
# plot_scatter = 
#   results |>
#   filter(Smin > 0 & Iikmin > 0) |>
#   select(-c(1, 3, 5:11)) |>
#   pivot_longer(-c(Smin, Ikmin), names_to = 'parameter') |>
#   group_by(parameter, value) |>
#   summarise_at(c('Smin', 'Ikmin'), mean) |>
#   ungroup() |>
#   pivot_longer(c(Smin, Ikmin), names_to = 'min', values_to = 'outcome') |>
#   ggplot(aes(value, outcome, group = min, color = min)) +
#   geom_line() +
#   geom_point() +
#   facet_wrap(~parameter, scales = 'free') +
#   ggtitle('Smin > 0 & Iikmin > 0')
# 
# plot_scatter2 = 
#   results |>
#   filter(!(Smin > 0 & Iikmin > 0)) |>
#   select(-c(1, 3, 5:11)) |>
#   pivot_longer(-c(Smin, Ikmin), names_to = 'parameter') |>
#   group_by(parameter, value) |>
#   summarise_at(c('Smin', 'Ikmin'), mean) |>
#   ungroup() |>
#   pivot_longer(c(Smin, Ikmin), names_to = 'min', values_to = 'outcome') |>
#   ggplot(aes(value, outcome, group = min, color = min)) +
#   geom_line() +
#   geom_point() +
#   facet_wrap(~parameter, scales = 'free') +
#   ggtitle('Smin <= 0 | Iikmin <= 0')


