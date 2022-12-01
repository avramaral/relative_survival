source("header.R")

##################################################

model <- "LN_ABST"
dist <- gsub(pattern = "_", replacement = "", x = substring(text = model, first = c(1, 4), last = c(3, 7))[1])

d <- readRDS(file = paste("FITTED_MODELS/", dist, "/d_", model, ".rds", sep = ""))
r <- readRDS(file = paste("FITTED_MODELS/", dist, "/",   model, ".rds", sep = ""))
fit <- r$fit

X_tilde <- d$X_tilde
X <- d$X
if (!substring(text = model, first = c(1, 6), last = c(5, 7))[2] == "XX") {
  region <- d$region
} else {
  region <- NULL
}

time <- seq(from = 0.0, to = 5, by = 0.2) # Change step-size, if needed

##################################################

####################
# Estimated curves #
####################

res <- result_processing(fit = fit, 
                         model = model, 
                         time = time, 
                         X_tilde = X_tilde, 
                         X = X,
                         region = region)

excHaz    <- res$excHaz
excCumHaz <- res$excCumHaz
netSur    <- res$netSur

##################################################

#################
# Grouping step #
#################

# If "region" and "XX" model, we have to load the data again
# data <- readRDS(file = "./path/file.rds")
# region <- data$region

group <- region # Assumes it goes from 1:N
grouped_res <- group_marginal_quantities(group = group, excHaz = excHaz, excCumHaz = excCumHaz, netSur = netSur)

