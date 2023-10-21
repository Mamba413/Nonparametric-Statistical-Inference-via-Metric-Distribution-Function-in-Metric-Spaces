rm(list = ls()); gc(reset = TRUE)
your_path <- "/public/home/zhujin/bddistribution_great/sim_testing"
setwd(your_path)

homogeneity_simulator <- function(seed, case = "1-1", num = 3000, p = 2) {
  each_n <- num / 3
  if (case == "1-1") {
    x <- matrix(rnorm(num * p), ncol = p)
  }
  if (case == "1-2") {
    x1 <- matrix(rnorm(each_n * p), ncol = p)
    x2 <- matrix(rnorm(each_n * p, mean = -0.2), ncol = p)
    x3 <- matrix(rnorm(each_n * p, mean = 0.2), ncol = p)
    x <- rbind(x1, x2, x3)
  }
  if (case == "2-1") {
    x <- rmovMF(num, c(2, 0) / sqrt(4), alpha = 1)
  }
  if (case == "2-2") {
    x1 <- rmovMF(each_n, c(2, 0) / sqrt(4))
    x2 <- rmovMF(each_n, c(2, 1) / sqrt(5))
    x3 <- rmovMF(each_n, c(2, 2) / sqrt(8))
    x <- rbind(x1, x2, x3)
  }
  if (case == "3-1") {
    cov_mat <- 0.0^abs(outer(1:p, 1:p, "-"))
    x <- rWishart(n = num, df = 8, Sigma = cov_mat)
  }
  if (case == "3-2") {
    cov_mat1 <- (0.0)^abs(outer(1:p, 1:p, "-"))
    cov_mat2 <- (-0.1)^abs(outer(1:p, 1:p, "-"))
    cov_mat3 <- (0.1)^abs(outer(1:p, 1:p, "-"))
    x1 <- rWishart(n = each_n, df = 8, Sigma = cov_mat1)
    x2 <- rWishart(n = each_n, df = 8, Sigma = cov_mat2)
    x3 <- rWishart(n = each_n, df = 8, Sigma = cov_mat3)
    x <- abind::abind(x1, x2, x3)
  }
  
  dst_type <- strsplit(case, "-")[[1]][1]
  if (dst_type == "1") {
    x <- dist(x)
  } else if (dst_type == "2") {
    x <- as.dist(Ball::nhdist(x))
  } else if (dst_type == "3") {
    x <- as.dist(CovDist(x, method = "C"))
  }
  
  y <- matrix(0, nrow = num, ncol = 3)
  index <- rep(1:3, each = each_n)
  y[, 1] <- 1.0 * (index == 1)
  y[, 2] <- 1.0 * (index == 2)
  y[, 3] <- 1.0 * (index == 3)
  y <- dist(y)
  
  permute_p <- bcov.test(x, y, distance = TRUE, 
                         weight = "rbf")[["p.value"]]
  limit_p <- bcov.test(x, y, distance = TRUE, 
                       method = "limit", weight = "rbf")[["p.value"]]
  c(permute_p, limit_p)
}


simulation_wrapper <- function(seed, case = "1-1", 
                               num_list = c(210, 300, 420, 510, 600), 
                               p = 2)
{
  pvalue_list <- list()
  for(i in 1:length(num_list)) {
    pvalue_list[[i]] = homogeneity_simulator(seed = seed, case = case, 
                                             p = p, num = num_list[i])
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
sfInit(parallel = TRUE, cpus = 15)
sfLibrary(Ball)
sfLibrary(movMF)
sfLibrary(CovTools)
sfExportAll()

nrep <- 500

print(paste0("------------ Example 1-1 ------------"))
example_11 <- sfLapply(1:nrep, simulation_wrapper, case = "1-1")
res11 <- t(calculate_power(example_11)); print(res11)
save(res11, file="spectral_homo_1_1.rda")

print(paste0("------------ Example 1-2 ------------"))
example_12 <- sfLapply(1:nrep, simulation_wrapper, case = "1-2")
res12 <- t(calculate_power(example_12)); print(res12)
save(res12, file="spectral_homo_1_2.rda")

print(paste0("------------ Example 2-1 ------------"))
example_21 <- sfLapply(1:nrep, simulation_wrapper, case = "2-1")
res21 <- t(calculate_power(example_21)); print(res21)
save(res21, file="spectral_homo_2_1.rda")

print(paste0("------------ Example 2-2 ------------"))
example_22 <- sfLapply(1:nrep, simulation_wrapper, case = "2-2")
res22 <- t(calculate_power(example_22)); print(res22)
save(res22, file="spectral_homo_2_2.rda")

print(paste0("------------ Example 3-1 ------------"))
example_31 <- sfLapply(1:nrep, simulation_wrapper, case = "3-1")
res31 <- t(calculate_power(example_31)); print(res31)
save(res31, file="spectral_homo_3_1.rda")

print(paste0("------------ Example 3-2 ------------"))
example_32 <- sfLapply(1:nrep, simulation_wrapper, case = "3-2")
res32 <- t(calculate_power(example_32)); print(res32)
save(res32, file="spectral_homo_3_2.rda")

sfStop()


