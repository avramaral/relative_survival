library(spBayesSurv)
library(sp)
library(R2BayesX)

data(LeukSurv)

n <- nrow(LeukSurv)
time <- LeukSurv$time / 365 

quintiles <- quantile(LeukSurv$tpi, c(0.2, 0.4, 0.6, 0.8))

for (i in 1:n) {
  if (LeukSurv$tpi[i] <= quintiles[1]) {
  LeukSurv$tpi_transf[i] <- 1 } else if (quintiles[1] < LeukSurv$tpi[i] & LeukSurv$tpi[i] <= quintiles[2]) {
  LeukSurv$tpi_transf[i] <- 2 } else if (quintiles[2] < LeukSurv$tpi[i] & LeukSurv$tpi[i] <= quintiles[3]) {
  LeukSurv$tpi_transf[i] <- 3 } else if (quintiles[3] < LeukSurv$tpi[i] & LeukSurv$tpi[i] <= quintiles[4]) {
  LeukSurv$tpi_transf[i] <- 4 } else {
  LeukSurv$tpi_transf[i] <- 5 }
}

leuk <- data.frame(index = 1:n,
                   time = time,
                   obs  = LeukSurv$cens,
                   region = LeukSurv$district,
                   age = LeukSurv$age,
                   sex = LeukSurv$sex,
                   wbc = LeukSurv$wbc,
                   dep = LeukSurv$tpi,
                   dep_transf = LeukSurv$tpi_transf)

# Life tables
lt <- as.data.frame(read.table("DATA/ENGLAND_LT_2010_2015.txt", header = T))
lt <- lt[lt$X_year == 2010, ] 
lt <- lt[lt$gor == 2, ]
lt <- lt[, c('sex', 'dep', 'age', 'rate')]
colnames(lt) <- c('sex', 'dep_transf', 'age', 'rate')
lt$sex <- lt$sex - 1

# Merging data
data <- merge(x = leuk, y = lt, by = c("sex", "dep_transf", "age"), all.x = T)
data <- data[order(data$index), ]
data <- cbind(data[, c(5:7, 3, 1, 8:10)])
rownames(data) <- 1:nrow(data)
colnames(data) <- c(colnames(data)[1:7], "pop.haz")

saveRDS(object = data, file = "DATA/leuk.rds")

# Geo Boundaries
map <- read.bnd(file = system.file("otherdata/nwengland.bnd", package = "spBayesSurv"))
map <- bnd2sp(bndObject = map)
proj4string(map) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

saveRDS(object = map, file = "DATA/nwengland_map.rds")
