library(SimLT)
library(msm)
library(lubridate)
library(igraph)
library(rstan)
library(parallel)

source("models.R")
source("utils.R")
source("data_stan.R")
source("fit_stan.R")
source("distributions.R")
source("SIMULATION/population_hazard.R")
source("SIMULATION/excess_hazard.R")

LT <- as.data.frame(read.table("DATA/ENGLAND_LT_2010_2015.txt", header = T))

##################################################

model_data <- "PGWABST"
prop <- 50
dist <- gsub(pattern = "_", replacement = "", x = substring(text = model_data, first = c(1, 4), last = c(3, 7))[1])

sample_size <- 2000
year_init <- 2010
year_last <- 2014

precision_tilde <- 10
precision <- 10

last_day <- ymd(paste(year_last, "-12-31", sep = ""))
sample_size <- ifelse(test = (sample_size %% 2 == 0), yes = sample_size, no = sample_size + 1)
X_fem <- simDesMatrix(seed = 999, n = (sample_size / 2), admin.cens = last_day, scale.age = F, site = "lung", sex = "female")
X_mal <- simDesMatrix(seed = 999, n = (sample_size / 2), admin.cens = last_day, scale.age = F, site = "lung", sex = "male")
fixed_times <- c(seq(from = ymd(paste(year_init, "-01-01", sep = "")), to = ymd(paste(year_last, "-01-01", sep = "")), by = "year"), last_day)

hrates_fem <- compute_hazard_rates(X = X_mal, fixed_times = fixed_times, sex = 1)
hrates_mal <- compute_hazard_rates(X = X_fem, fixed_times = fixed_times, sex = 2) 

N_sim <- 50
data <- list()

progressbar <- txtProgressBar(min = 1, max = N_sim, initial = 1)

for (k in 1:N_sim) {
  set.seed(k + 999)
  
  sim_fem <- sim_pophaz(seed = k + 999, lst = hrates_fem)
  sim_mal <- sim_pophaz(seed = k + 999, lst = hrates_mal)
  
  sim_pop <- list()
  sim_pop$sim.pop <- c(sim_fem$sim.pop, sim_mal$sim.pop)
  sim_pop$status <- c(sim_fem$status, sim_mal$status)
  
  hist(sim_pop$sim.pop)
  mean(sim_pop$status)
  
  ##################################################
  
  if (model_data %in% c("PGWABST", "LN_ABST", "LL_ABST")) {
    s_re <- simulate_re(struc = "ICAR", precision_tilde = precision_tilde, precision = precision)
  } else if (model_data %in% c("PGWABCD", "LN_ABCD", "LL_ABCD")) {
    s_re <- simulate_re(struc = "IID",  precision_tilde = precision_tilde, precision = precision)
  } else if (model_data %in% c("PGWABXX", "LN_ABXX", "LL_ABXX")) {
    s_re <- simulate_re(struc = "NONE")
  }
  
  re_tilde <- s_re$re_tilde
  re <- s_re$re
    
  # Design Matrix
  desMat <- rbind(X_fem, X_mal)
  desMat <- cbind(desMat, sex = c(rep(x = 0, times = nrow(X_fem)), rep(x = 1, times = nrow(X_mal))))
  desMat$dep <- scale(desMat$dep)
  desMat$age <- scale(desMat$age)
  dep <- c(X_fem$dep, X_mal$dep)
  
  X_tilde <- matrix(data = desMat$age, ncol = 1, byrow = F) 
  X <- matrix(data = c(desMat$age, desMat$sex, desMat$dep), ncol = 3, byrow = F) 
  
  idxs <- list()
  for (i in 1:length(unique(desMat$gor))) {
    idxs[[i]] <- which(desMat$gor == i)
  }
  
  alphas <- list()
  betas <- list()
  
  for (i in 1:length(unique(dep))) {
    if (prop == 50) {
      alphas[[i]] <- c(2)
      betas[[i]] <-  c(2, - 1, - 2)
    } else if (prop == 75) {
      alphas[[i]] <- c(2)
      betas[[i]] <-  c(2, - 1, - 2)
    } else {
      stop("Choose a valid ratio for 'observed / total individuals.'")
    }
  }
  
  if (dist == "LN") {
    if (prop == 50) {
      pars <- list(mu = 1, sigma = 1)
    } else if (prop == 75) {
      pars <- list(mu = 0, sigma = 0.1)
    }
  } else if (dist == "PGW") {
    if (prop == 50) {
      pars <- list(eta = 1.5, nu = 0.75, theta = 1.5) 
    } else if (prop == 75) {
      pars <- list(eta = 1, nu = 2, theta = 1) 
    }
  }
  
  C_1 <- rep(x = 5, times = nrow(desMat))
  C_2 <- rexp(n = nrow(desMat), rate = 0.01)
  sim_exc <- list()
  sim_exc$sim.exc <- c()
  sim_exc$status <- c()
  tim_exc <- simulate_RS_MEGH(dist = dist, pars = pars, X_tilde = X_tilde, X = X, alphas = alphas, betas = betas, dep = dep, re_tilde = re_tilde, re = re, idxs = idxs)
  for (i in 1:nrow(desMat)) { sim_exc$sim.exc <- c(sim_exc$sim.exc, min(tim_exc[i], C_1[i], C_2[i])) }
  for (i in 1:nrow(desMat)) { sim_exc$status <- c(sim_exc$status, ifelse(test = ((tim_exc[i] < C_1[i]) & (tim_exc[i] < C_2[i])), yes = 1, no = 0)) }
    
  hist(sim_exc$sim.exc)
  mean(sim_exc$status)
  
  ##################################################
  
  # Final survival times and censoring status
  times <- vector()
  status <- vector()
  for (i in 1:nrow(desMat)) {
    if (sim_pop$status[i] == 1) {
      times[i] <- min(sim_pop$sim.pop[i], sim_exc$sim.exc[i])
      status[i] <- 1
    } else if ((sim_pop$status[i] == 0) & (sim_exc$sim.exc[i] < sim_pop$sim.pop[i])) {
      times[i] <- sim_exc$sim.exc[i]
      if (sim_exc$status[i] == 1) {
        status[i] <- 1  
      } else {
        status[i] <- 0
      }
    } else if ((sim_pop$status[i] == 0) & (sim_exc$sim.exc[i] >= sim_pop$sim.pop[i])) {
      times[i] <- sim_pop$sim.pop[i]
      status[i] <- 0
    }
  }
  times <- ifelse(test = (times < (1 / 365)), yes = (1 / 365), no = times) # Making survival at least one day
  
  # hist(times, breaks = 10)
  # mean(status)
  
  ##################################################
  
  day_out <- rbind(X_fem, X_mal)[, "date.diag"] + (times * 365)
  age_out <- rbind(X_fem, X_mal)[, "age"] + times
  MAT_ID <- as.data.frame(cbind(1:nrow(desMat), year(day_out), c(rep(x = 0, times = nrow(X_fem)), rep(x = 1, times = nrow(X_mal))) + 1, rbind(X_fem, X_mal)[, "dep"], floor(age_out), rbind(X_fem, X_mal)[, "gor"]))
  colnames(MAT_ID) <- c("index", "X_year", "sex", "dep", "age", "gor")
  pop_haz <- merge(x = MAT_ID, y = LT, by = colnames(MAT_ID)[-1], all.x = TRUE, sort = FALSE)
  pop_haz <- pop_haz[order(pop_haz$index), ]
  pop_haz <- pop_haz$rate
  
  data[[k]] <- data.frame(time = times, obs = status, region = desMat$gor, age = desMat$age, sex = desMat$sex, dep = desMat$dep, pop.haz = pop_haz)
  
  setTxtProgressBar(progressbar, k) 
}

close(progressbar)

save.image(file = paste("SIMULATION/DATA/", model_data, "_n_", sample_size, "_prop_", prop, ".RData", sep = ""))
