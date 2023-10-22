rm(list = ls()); gc(reset = TRUE)
your_path <- "~"
setwd(your_path)

generate_spd_data <- function(value) {
  simplify2array(lapply(value, function(theta) {
    cov_mat <- 8 * diag(2)
    cov_mat[2, 1] <- cov_mat[1, 2] <- theta
    rWishart(n = 1, df = 35, Sigma = cov_mat)[, , 1]
  }))
}

joint_independence_simulator <- function(seed, case = "1-1", num = 3000) {
  set.seed(seed)
  if (case == "1-1") {
    x <- matrix(rnorm(num * 3), ncol = 3)
    z <- x[, 1]
    y <- x[, 2]
    x <- x[, 3]
  }
  if (case == "1-2") {
    cov_mat <- (0.25)^abs(outer(1:3, 1:3, "-"))
    x <- rmvnorm(num, mean = rep(0, 3), sigma = cov_mat)
    z <- x[, 1]
    y <- x[, 2]
    x <- x[, 3]
  }
  if (case == "2-1") {
    x <- matrix(rnorm(num * 3 * 2), nrow = num)
    x <- rbind(x[, 1:2], x[, 3:4], x[, 5:6])
    x <- t(apply(x, 1, function(theta) {
        rmovMF(1, theta, alpha = 1)
      }))
    z <- x[1:num,]
    y <- x[(num+1):(2*num),]
    x <- x[(2*num+1):(3*num),]
  }
  if (case == "2-2") {
    cov_mat <- (0.6)^abs(outer(1:6, 1:6, "-"))
    x <- rmvnorm(num, mean = rep(0, 6), sigma = cov_mat)
    x <- rbind(x[, 1:2], x[, 3:4], x[, 5:6])
    x <- t(apply(x, 1, function(theta) {
        rmovMF(1, theta, alpha = 1)
    }))
    z <- x[1:num,]
    y <- x[(num+1):(2*num),]
    x <- x[(2*num+1):(3*num),]
  }
  if (case == "3-1") {
    x <- matrix(rnorm(num * 3), ncol = 3)
    z <- generate_spd_data(x[, 1])
    y <- generate_spd_data(x[, 2])
    x <- generate_spd_data(x[, 3])
  }
  if (case == "3-2") {
    cov_mat <- (0.6)^abs(outer(1:3, 1:3, "-"))
    x <- rmvnorm(num, mean = rep(0, 3), sigma = cov_mat)
    z <- generate_spd_data(x[, 1])
    y <- generate_spd_data(x[, 2])
    x <- generate_spd_data(x[, 3])
  }
  
  dst_type <- strsplit(case, "-")[[1]][1]
  if (dst_type == "1") {
    x <- dist(x)
    y <- dist(y)
    z <- dist(z)
  } else if (dst_type == "2") {
    x <- as.dist(Ball::nhdist(x))
    y <- as.dist(Ball::nhdist(y))
    z <- as.dist(Ball::nhdist(z))
  } else if (dst_type == "3") {
    x <- as.dist(CovDist(x, method = "C"))
    y <- as.dist(CovDist(y, method = "C"))
    z <- as.dist(CovDist(z, method = "C"))
  }

  dx <- list(as.vector(x), as.vector(y), as.vector(z))
  dx2 <- list(x, y, z)
  dx3 <- list(as.matrix(x), as.matrix(y), as.matrix(z))
  
  permute_p <- bcov.test(dx3, distance=TRUE, num.permutations = 199)[["p.value"]]
  limit_p <- KBCovLimit(dx, n=num)[["p.value"]]

  c(permute_p, limit_p)
}


simulation_wrapper <- function(seed, case = "1-1", 
                               num_list = c(120, 240, 480, 960, 1920))
{
  pvalue_list <- list()
  for(i in 1:length(num_list)) {
    pvalue_list[[i]] = joint_independence_simulator(seed = seed, case = case, 
                                                    num = num_list[i])
  }
  pvalue_list <- do.call("rbind", pvalue_list)
  pvalue_list < 0.05
}

calculate_power <- function(pvalue_mat_list){
  cols <- ncol(pvalue_mat_list[[1]])
  pvalue_mat_list <- lapply(1:cols, function(i){
    sapply(pvalue_mat_list,function(x) {x[,i]})
  })
  power_tab <- sapply(pvalue_mat_list, rowMeans)
  power_tab
}

library(snowfall)
sfInit(parallel = TRUE, cpus = 50)
sfLibrary(Ball)
sfLibrary(movMF)
sfLibrary(CovTools)
sfLibrary(mvtnorm)
sfExportAll()

nrep <- 100

print(paste0("------------ Example 1-1 ------------"))
example_11 <- sfLapply(1:nrep, simulation_wrapper, case = "1-1")
res11 <- t(calculate_power(example_11)); print(res11)
save(res11, file="spectral_ji_1_1.rda")

print(paste0("------------ Example 1-2 ------------"))
example_12 <- sfLapply(1:nrep, simulation_wrapper, case = "1-2")
res12 <- t(calculate_power(example_12)); print(res12)
save(res12, file="spectral_ji_1_2.rda")

print(paste0("------------ Example 2-1 ------------"))
example_21 <- sfLapply(1:nrep, simulation_wrapper, case = "2-1")
res21 <- t(calculate_power(example_21)); print(res21)
save(res21, file="spectral_ji_2_1.rda")

print(paste0("------------ Example 2-2 ------------"))
example_22 <- sfLapply(1:nrep, simulation_wrapper, case = "2-2")
res22 <- t(calculate_power(example_22)); print(res22)
save(res22, file="spectral_ji_2_2.rda")

print(paste0("------------ Example 3-1 ------------"))
example_31 <- sfLapply(1:nrep, simulation_wrapper, case = "3-1")
res31 <- t(calculate_power(example_31)); print(res31)
save(res31, file="spectral_ji_3_1.rda")

print(paste0("------------ Example 3-2 ------------"))
example_32 <- sfLapply(1:nrep, simulation_wrapper, case = "3-2")
res32 <- t(calculate_power(example_32)); print(res32)
save(res32, file="spectral_ji_3_2.rda")

sfStop()
