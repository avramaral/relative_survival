compile_model <- function (model) {
  
  if (!validate_model(model = model)) {
    stop("Select a valid model.")
  }
  
  original_model <- model
  name  <- substring(text = model, first = c(1, 4), last = c(3, 7))
  dist  <- name[1]
  model <- name[2]
  
  if (model == "ABDD") { model <- "ABCC" }
  if (model == "AADD" | model == "BBCC" | model == "BBDD") { model <- "AACC" }
  if (model == "ABTT") { model <- "ABSS" }
  if (model == "AATT" | model == "BBSS" | model == "BBTT") { model <- "AASS" }
  if (model == "ABZZ") { model <- "ABYY" }
  if (model == "AAZZ" | model == "BBYY" | model == "BBZZ") { model <- "AAYY" }
  if (model == "BBXX") { model <- "AAXX" }
  
  mod <- stan_model(file = paste("MODELS/", gsub(pattern = "_", replacement = "", x = dist), "/", dist, model, ".stan", sep = ""))
  list(model = original_model, mod = mod)
}

fit_stan <- function (mod, data, chains = 4, iter = 4e3, warmup = 2e3, adapt_delta = 0.8, max_treedepth = 10, ...) {

  chains <- chains
  iter <- iter
  warmup <- warmup
  
  start_time <- Sys.time()
  
  fit <- sampling(object = mod$mod,
                  data = data,
                  chains = chains,
                  iter = iter,
                  warmup = warmup,
                  pars = c("lp_tilde", "lp", "excessHaz", "cumExcessHaz"),
                  include = F,
                  control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
                  cores = getOption(x = "mc.cores", default = detectCores()))
  
  end_time <- Sys.time()
  time_taken <- end_time - start_time
  
  list(model = mod$model, chains = chains, iter = iter, warmup = warmup, adapt_delta = adapt_delta, max_treedepth = max_treedepth, fit = fit, time_taken = time_taken)
}

