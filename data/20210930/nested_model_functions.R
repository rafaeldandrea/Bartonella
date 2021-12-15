## Rodents vs pathogens model -- nested infectivity


## Libraries

library(tidyverse)
library(magrittr)
library(rootSolve)

## Functions

parms_to_list = 
  function(...) list(...)
  
dependent_parms = 
  function(parms, omega = NULL, xbar = NULL){
    
    list2env(parms, envir = environment())
    
    fecundity = seq(bmin, bmax, l = number_of_species)
    
    mortality = fecundity * (1 - epsilon)
    
    if(is.null(xbar)) xbar = fecundity[1:min(number_of_species, number_of_strains)]
    
    if(is.null(omega)) omega = omega0 * exp(-k * xbar)
    
    xminusxbar = outer(fecundity, xbar, '-')
    
    omegamat = matrix(omega, number_of_species, number_of_strains, byrow = TRUE)
    
    omegainv = 1 / (1e-16 + omegamat)
    
    beta = 
      beta0 * 
      exp(l * fecundity - lomega * omegamat) * 
      exp(-xminusxbar ^ 2 * omegainv ^ 2) / 
      mean(exp(-xminusxbar ^ 2 * omegainv ^ 2))
    
    delta = 
      delta0 * 
      exp(-(xminusxbar / sigma_delta) ^ 2) / 
      mean(exp(-(xminusxbar / sigma_delta) ^ 2))
    
    gamma = 
      gamma0 * 
      exp(-(xminusxbar / sigma_gamma) ^ 2) / 
      mean(exp(-(xminusxbar / sigma_gamma) ^ 2))
    
    list(
      fecundity = fecundity,
      mortality = mortality,
      a = epsilon * bmin / Kmin,
      xbar = xbar,
      omega = omega,
      xminusxbar = xminusxbar,
      omegainv = omegainv,
      omegamat = omegamat,
      beta = beta,
      delta = delta,
      gamma = gamma
    ) %>%
      return()
    
  }

eco_dynamics = 
  function(t, y, parms){
    
    list2env(parms, envir = environment())
    
    ## Dependent parameters
    parms2 = dependent_parms(parms) 
    
    list2env(parms2, envir = environment())
    
    y = pmax(0, y)
    
    S = y[1:number_of_species]
    
    Iik = 
      matrix(
        y[(number_of_species + 1): length(y)],
        number_of_species, 
        number_of_strains
      )
    
    Ik = colSums(Iik)
    SiIk = outer(S, Ik)
    
    N = S + rowSums(Iik)
    
    dSdt = 
      (fecundity - mortality - a * N) * S + rowSums((fecundity - a * N + gamma) * Iik) -
      rowSums(beta * SiIk)
    
    dIdt = beta * SiIk - (mortality * (1 + delta) + gamma) * Iik
    
    return(list(c(as.numeric(dSdt), as.numeric(dIdt))))
  }

evo_dynamics = 
  function(t, y, parms){
    
    list2env(parms, envir = environment())
    
    ## Dependent parameters
    parms2 = dependent_parms(parms, omega = y) 
    
    list2env(parms2, envir = environment())
    
    zeta = 1 / (mortality * (1 + delta) + gamma)
    
    Q = t(zeta * beta)
    
    M = (1 - (fecundity + gamma) * zeta) * beta
    
    Si = try(pmax(0, solve(Q, rep(1, number_of_species))), silent = TRUE)
    
    Ik = try(pmax(0, solve(M, fecundity - mortality)), silent = TRUE)
    
    if(class(Si) == 'try-error') Si = rep(0, number_of_species)
    if(class(Ik) == 'try-error') Ik = rep(0, number_of_strains)
    
    Iik = zeta * beta * outer(Si, Ik)
    
    xbar_prime = - 1 / k * omegainv
    
    beta_prime = beta * 2 * xminusxbar * omegainv ^ 2 * (xbar_prime + xminusxbar * omegainv)
    
    zeta_prime = -zeta ^ 2 * 2 * xminusxbar * xbar_prime *
      (mortality * delta / sigma_delta ^ 2 + gamma / sigma_gamma ^ 2)
    
    dwdt =
      colSums(beta_prime * Si) +
      colSums(1 / zeta ^ 2 * zeta_prime * t(t(Iik) / (1e-16 + Ik)))
    
    return(list(dwdt))
  }

evo_dynamics_fixedxbar_no_self_regulation = 
  function(t, y, parms, xbar = NULL){
    
    list2env(parms, envir = environment())
    
    ## Dependent parameters
    parms2 = dependent_parms(parms, omega = y, xbar = xbar) 
    
    list2env(parms2, envir = environment())
    
    zeta = 1 / (mortality * (1 + delta) + gamma)
    
    Q = t(zeta * beta)
    
    Si = pmax(0, solve(Q, rep(1, number_of_species)))
    
    beta_prime = beta * (2 * xminusxbar ^ 2 * omegainv ^ 3 - lomega)
    
    dwdt = colSums(beta_prime * Si)
    
    return(list(dwdt))
  }

evo_equilibrium = 
  function(y, parms){
    multiroot(
      f = unlist(evo_dynamics(t = 0, y = y, parms = parms)),
      start = y,
      positive = TRUE,
      parms = parms
    ) 
  }

evo_dynamics_fixedxbar_manual = 
  function(t, y, parms, omegamin){
    w = y
    delta_t = diff(t)
    nstrains = parms$number_of_strains
    
    wrec = 
      tibble(
        strain = seq(nstrains),
        omega = w
      )
    
    for(dt in delta_t){
      w = pmax(omegamin, w + dt * unlist(evo_dynamics_fixedxbar(t = 0, y = w, parms)))
      wrec %<>%
        bind_rows(
          tibble(
            strain = seq(nstrains),
            omega = w
          )
        )
    }
  
    wrec %>%
      mutate(time = rep(t, each = nstrains)) %>%
      select(time, everything()) %>%
      return
  }

eco_equilibrium = 
  function(parms, ...){
    
    list2env(parms, envir = environment())
    
    ## Dependent parameters
    parms2 = dependent_parms(parms, ...) 
    
    list2env(parms2, envir = environment())
    
    ## Model Equilibrium
    zeta = 1 / (mortality * (1 + delta) + gamma)
    
    Q = t(zeta * beta)
    
    M = (1 - (fecundity + gamma) * zeta) * beta
    
    Si = solve(Q, rep(1, number_of_species))
    
    Ik = solve(M, (fecundity - mortality - a * Si))
    
    Iik = zeta * beta * outer(Si, Ik)
    
    parms_tbl = as_tibble(parms)
    
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
        Ik = Ik
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

wrapper = 
  function(
    number_of_species,
    number_of_strains,
    bmin,
    bmax,
    beta0,
    omega0,
    k,
    l,
    lomega,
    delta0,
    gamma0,
    sigma_delta,
    sigma_gamma,
    epsilon, 
    Kmin,
    ...
  ){
    parms = 
      parms_to_list(
        number_of_species = number_of_species,
        number_of_strains = number_of_strains,
        bmin = bmin,
        bmax = bmax,
        beta0 = beta0,
        omega0 = omega0,
        k = k,
        l = l,
        lomega = lomega,
        delta0 = delta0,
        gamma0 = gamma0,
        sigma_delta = sigma_delta,
        sigma_gamma = sigma_gamma,
        epsilon = epsilon,
        Kmin = Kmin
      )
    
    parms %>%
      eco_equilibrium %>%
      return()
    
  }


evo_dynamics_fixedxbar_manual = 
  function(t, y, parms, evolve_interval, verbose, omegamin = 0){
    list2env(parms, envir = environment())
    
    ## Dependent parameters
    parms2 = dependent_parms(parms) 
    
    list2env(parms2, envir = environment())
    
    S = y[1:number_of_species]
    
    Iik = 
      matrix(
        y[(number_of_species + 1): length(y)],
        number_of_species, 
        number_of_strains
      )
    
    Ik = colSums(Iik)
    SiIk = outer(S, Ik)
    
    dt = diff(t)
    
    time = 0 
    step = 0
    
    SI_dtf = 
      tibble(
        time = time, 
        species = seq(number_of_species), 
        S = S
      ) %>% 
      inner_join(
        expand_grid(
          time = time, 
          species = seq(number_of_species), 
          strain = seq(number_of_strains)
        ) %>% 
          arrange(strain, species) %>%
          mutate(Iik = as.numeric(Iik)),
        by = c('time', 'species')
      )
    
    omega_dtf = 
      tibble(
        time = time,
        strain = seq(number_of_strains),
        omega = omega
      )
    
    foo = 
      expand_grid(
        strain = seq(number_of_strains), 
        species = seq(number_of_species)
      )
    
    while(time < max(t) & !is.na(min(S)) & !is.na(min(Iik))){
      
      step = step + 1
      
      ## Check for evolution of niche breadth
      if(time %% evolve_interval == 0){
        
        if(verbose == TRUE) writeLines(paste('time = ', time))
        
        beta_prime = beta * (2 * xminusxbar ^ 2 * omegainv ^ 3 - lomega)
        
        dwdt = colSums(beta_prime * S)
        
        omega = pmax(omegamin, omega + dt[step] * dwdt)
        
        omegamat = matrix(omega, number_of_species, number_of_strains, byrow = TRUE)
        omegainv = 1 / (1e-16 + omegamat)
        
        beta = beta0 * exp(l * fecundity - lomega * omegamat) * exp(-xminusxbar ^ 2 * omegainv ^ 2)
        
        ## Update time
        time = time + dt[step]
        
        omega_dtf %<>%
          bind_rows(
            tibble(
              time = time, 
              strain = seq(number_of_strains),
              omega = omega
            )
          )
        
        SI_dtf %<>%
          bind_rows(
            tibble(
              time = time, 
              species = seq(number_of_species), 
              S = S
            ) %>% 
              inner_join(
                foo %>%
                  mutate(
                    time = time,
                    Iik = as.numeric(Iik)
                  ),
                by = c('time', 'species')
              )
          )
      } else{
        
        ## Update time
        time = time + dt[step]
      }
      
      
      ## Update abundances
      
      N = S + rowSums(Iik)
      
      dSdt = 
        (fecundity - mortality - a * N) * S + rowSums((fecundity - a * N + gamma) * Iik) -
        rowSums(beta * SiIk)
      
      dIdt = beta * SiIk - (mortality * (1 + delta) + gamma) * Iik
      
      S = pmax(0, S + dt[step] * dSdt)
      Iik = matrix(pmax(0, Iik + dt[step] * dIdt), number_of_species, number_of_strains)
      
      Ik = colSums(Iik)
      SiIk = outer(S, Ik)
      
    }
    
    SI_dtf %>%
      left_join(
        omega_dtf,
        by = c('time', 'strain')
      ) %>%
      return
    
  }