compute_ages <- function (dates, dob, ...) {
  time_length(difftime(ymd(dates), ymd(dob)), "years")
}

compute_hazard_rates <- function (X, fixed_times, sex = 1, ...) {
  times <- list()
  dates <- list() 
  hrate <- list()
  
  year_init <- as.numeric(year(fixed_times[1]))
  sample_size <- nrow(X)
  
  progressbar <- txtProgressBar(min = 1, max = sample_size, initial = 1)
  for(i in 1:sample_size){
    dob <- as.POSIXlt(ymd(X[i, "date.diag"]) - (X[i, "age"] * 365))
    if (day(dob) == 29 & month(dob) == 2) {
      day(dob) <- 1
      month(dob) <- 3
    }
    # Birthday on year_init
    birthday <- dob
    year(birthday) <- year_init
    birthday <- ymd(birthday)
    
    dod <- ymd(X[i, "date.diag"]) 
    var_times <- c(seq(from = ymd(birthday), to = last_day, by = "year"), dod)
    time_pts <- sort(unique(c(fixed_times, var_times)))
    age_tmpt <- compute_ages(dates = time_pts, dob = dob)
    
    # Removing dates before the date of diagnosis
    ind <- which(time_pts == dod)
    if (ind == 1) {
      times[[i]] <- c(0, as.numeric(diff(time_pts)))
      dates[[i]] <- time_pts
      age_tmpt <- floor(age_tmpt)
    } else if(ind > 1 ) {
      times[[i]] <- c(0, as.numeric(diff(time_pts))[-c(1:(ind - 1))])
      dates[[i]] <- time_pts[-c(1:(ind - 1))]
      age_tmpt <- floor(age_tmpt[-c(1:(ind - 1))])
    }
    
    MAT_ID <- as.data.frame(cbind(1:length(age_tmpt), year(dates[[i]]), rep(sex, length(age_tmpt)), rep(X[i, "dep"], length(age_tmpt)), age_tmpt, rep(X[i, "gor"], length(age_tmpt))))
    colnames(MAT_ID) <- c("index", "X_year", "sex", "dep", "age", "gor")
    hrate[[i]] <- merge(x = MAT_ID, y = LT, by = colnames(MAT_ID)[-1], all.x = TRUE, sort = FALSE)
    hrate[[i]] <- hrate[[i]][order(hrate[[i]]$index), ]
    hrate[[i]] <- hrate[[i]]$rate
    
    setTxtProgressBar(progressbar, i) 
  }
  close(progressbar)
  
  list(times = times, hrates = hrate, dates = dates)
}

sim_pophaz <- function (seed, lst, days_year = 365, ...) {
  sim.pop <- vector()
  ind.pop <- vector()
  n <- length(lst$times)
  for(i in 1:n){
    times.temp <- cumsum(lst$times[[i]]) / days_year
    rates.temp <- lst$hrates[[i]]
    temp.sim <- rsim.pwexp(seed = (seed + i), n = 1, rate = rates.temp, t = times.temp)
    sim.pop[i] <- temp.sim$sim.pop
    ind.pop[i] <- temp.sim$status.pop
  }
  list(sim.pop = sim.pop, status = ind.pop) 
}
