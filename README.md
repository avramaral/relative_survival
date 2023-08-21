## Extended Excess Hazard Models for Spatially Dependent Survival Data

> **Amaral, AVR**, Rubio, FJ, Quaresma, M, Rodríguez-Cortés, FJ, and Moraga, P (2023). *Extended Excess Hazard Models for Spatially Dependent Survival Data*. Submitted. arXiv preprint. [arXiv:2302.09392](https://arxiv.org/abs/2302.09392).


### How to fit the `RS-SGH` model

As an example, we can fit the RS-SGH model for leukemia-diagnosed patients (Henderson et al., 2002) as follows

``` r
source("header.R") # load libraries and needed functions

data <- readRDS(file = "DATA/leuk.rds") # load "leukemia" data

# Optional
data$age <- scale(data$age)
data$wbc <- scale(data$wbc)
data$dep <- scale(data$dep)

map <- readRDS(file = "DATA/nwengland_map.rds") # load England map
adj_info <- adj_list(map = map) # create an object with information about the neighborhood structure

model <- "LN_ABST" 
dist <- gsub(pattern = "_", replacement = "", x = substring(text = model, first = c(1, 4), last = c(3, 7))[1]) # extract the distribution code from "model"

d <- data_stan(data = data, model = model, cov.tilde = c("age"), cov = c("age", "wbc", "sex", "dep"), nonlinear = c(), adj_info = adj_info) # create the data object
m <- compile_model(model = model) # compile the Stan model
r <- fit_stan(mod = m, data = d) # fit the model
```

#### References

> Henderson, R., Shimakura, S. and Gorst, D. (2002). *Modeling spatial variation in leukemia survival data*. Journal of the American Statistical Association 97, 965–972.

