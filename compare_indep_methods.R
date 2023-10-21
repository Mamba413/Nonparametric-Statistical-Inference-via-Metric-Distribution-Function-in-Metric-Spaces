rm(list = ls()); gc(reset = TRUE)
your_path <- "~/mdf/code/"
setwd(your_path)

measure.Choleksy.3d <- function(array3d){
  M        = dim(array3d)[3]
  outdist  = array(0, c(M,M))
  vec_chol = list()
  for (i in 1:M){
    vec_chol[[i]] = base::chol(array3d[,,i])
  }
  for (i in 1:(M-1)){
    cholA = vec_chol[[i]]
    for (j in (i+1):M){
      cholB = vec_chol[[j]]
      output = norm(cholA-cholB,"F")
      outdist[i,j] = output
      outdist[j,i] = output
    }
  }
  return(outdist)
}

band_matrix <- function(x, p) {
  mat <- diag(1, p, p)
  for(i in 1:p) {
    for (j in 1:p) {
      mat[i, j] <- (x)^(abs(i - j))
    }
  }
  mat
}

check_pd <- function(A){
  eigens <- eigen(A, symmetric = TRUE, only.values = TRUE)[["values"]]
  eigens <- min(eigens)
  if (eigens <= 0){
    return(FALSE)
  } else {
    return(TRUE)
  }
}

double_normalized_distance <- function(dx) {
  if (class(dx)[1] != "matrix") {
    dx <- as.matrix(dx)
  }
  size <- nrow(dx)
  C <- diag(size) - matrix(data = 1 / size, nrow = size, ncol = size)
  -C %*% dx %*% C
}

distance_multivariate <- function(data_list, R = 199) {
  data_list <- lapply(data_list, double_normalized_distance)
  tm <- multivariance::total.multivariance(x = data_list)
  # print(tm)
  tm.resampled <- numeric(R)
  for (r in 1:R) {
    set.seed(r)
    tm.resampled[r] <- multivariance::total.multivariance(multivariance::sample.cdms(data_list))
  }
  dcov_pvalue <- mean(c(tm.resampled >= tm, 1))
  dcov_pvalue
}

generate_spd_data <- function(value) {
  simplify2array(lapply(value, function(theta) {
    cov_mat <- matrix(c(1, 0.1, 0.1, 0.1, 1, 0.1, 0.1, 0.1, 1), nrow = 3)
    # cov_mat <- theta^abs(outer(1:6, 1:6, "-"))
    # diag(cov_mat) <- 8
    rWishart(n = 1, df = 3 + 8 * abs(theta), Sigma = cov_mat)[, , 1]
  }))
}

generate_spd_data2 <- function(value) {
  simplify2array(lapply(value, function(theta) {
    cov_mat <- matrix(c(1, 0.1, 0.1, 0.1, 1, 0.1, 0.1, 0.1, 1), nrow = 3)
    cov_mat <- cov_mat / (3 + 12 * abs(theta))
    rWishart(n = 1, df = (3 + 12 * abs(theta)), Sigma = cov_mat)[, , 1]
  }))
}

joint_independence_simulator <- function(seed, 
                                         case = "1-1", 
                                         num = 3000, 
                                         delta = NULL) {
  set.seed(seed = seed)
  if (case == "1-1") {
    x <- matrix(rnorm(num * 3), ncol = 3)
    z <- x[, 1]
    y <- x[, 2]
    x <- x[, 3]
  }
  if (case == "1-2") {
    delta <- ifelse(is.null(delta), 0.4, delta)
    cov_mat <- (delta)^abs(outer(1:3, 1:3, "-") != 0)
    x <- rmvnorm(num, mean = rep(0, 3), sigma = cov_mat)
    z <- x[, 1]
    y <- x[, 2]
    x <- x[, 3]
  }
  if (case == "1-3") {
    delta <- ifelse(is.null(delta), 0.6, delta)
    cov_mat1 <- (delta)^abs(outer(1:3, 1:3, "-") != 0)
    cov_mat2 <- (-delta)^abs(outer(1:3, 1:3, "-") != 0)
    x1 <- rmvnorm(num, mean = rep(0, 3), sigma = cov_mat1)
    x2 <- rmvnorm(num, mean = rep(0, 3), sigma = cov_mat2)
    mix_index <- rbinom(num, size = 1, prob = 0.5) 
    x <- x1 * mix_index + x2 * (1 - mix_index)
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
    delta <- ifelse(is.null(delta), 0.8, delta)
    p <- 3
    cov_mat <- (delta)^abs(outer(1:p, 1:p, "-") != 0)
    x <- rmvnorm(num, mean = rep(0, p), sigma = cov_mat)
    x <- c(x[, 1], x[, 2], x[, 3])
    x <- t(sapply(x, function(theta) {
      rmovMF(1, 2.8 * c(cos(theta), sin(theta), 0, 0), alpha = 1)
    }))
    z <- x[1:num,]
    y <- x[(num+1):(2*num),]
    x <- x[(2*num+1):(3*num),]
  }
  if (case == "2-3") {
    delta <- ifelse(is.null(delta), 0.9, delta)
    p <- 3
    cov_mat <- (delta)^abs(outer(1:p, 1:p, "-") != 0)
    x <- rmvnorm(num, mean = rep(0, p), sigma = cov_mat)
    x <- c(x[, 1], x[, 2], x[, 3])
    x <- t(sapply(x, function(theta) {
        rmovMF(1, 6.0 * c(abs(theta), 0, 0, 0), alpha = 1)
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
    delta <- ifelse(is.null(delta), 0.75, delta)
    cov_mat <- (delta)^abs(outer(1:3, 1:3, "-") != 0)
    x <- rmvnorm(num, mean = rep(0, 3), sigma = cov_mat)
    z <- generate_spd_data(x[, 1])
    y <- generate_spd_data(x[, 2])
    x <- generate_spd_data(x[, 3])
  }
  if (case == "3-3") {
    delta <- ifelse(is.null(delta), 0.8, delta)
    cov_mat1 <- (delta)^abs(outer(1:3, 1:3, "-") != 0)
    x <- rmvnorm(num, mean = rep(0, 3), sigma = cov_mat1)
    z <- generate_spd_data2(x[, 1])
    y <- generate_spd_data2(x[, 2])
    x <- generate_spd_data2(x[, 3])
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
    if (Sys.info()[1] == "Darwin") {
      x <- as.dist(measure.Choleksy.3d(x))
      y <- as.dist(measure.Choleksy.3d(y))
      z <- as.dist(measure.Choleksy.3d(z))
    } else {
      x <- as.dist(CovDist(x, method = "C"))
      y <- as.dist(CovDist(y, method = "C"))
      z <- as.dist(CovDist(z, method = "C"))
    }
  }
  
  dx1 <- list(as.vector(x), as.vector(y), as.vector(z))
  dx2 <- list(x, y, z)
  permute_p <- permute_p2 <- limit_p <- ed_value <- 1.0
  
  # limit_p <- KBCov(dx1, n=num)[["p.value"]]
  # permute_p <- bcov.test(dx2, distance=TRUE)[["p.value"]]
  
  # dx11 <- list(as.vector(x), as.vector(y), as.vector(z))
  # permute_p <- KBCov_permuted(dx11, n=num)[["p.value"]]
  
  dx12 <- list(as.matrix(x), as.matrix(y), as.matrix(z))
  permute_p2 <- Ball::bcov.test(dx12, num.permutations = 299, 
                                distance = TRUE)[["p.value"]]
  
  # TM test
  dx_tm <- list(as.matrix(x), as.matrix(y), as.matrix(z))
  ed_value <- distance_multivariate(dx_tm, R = 399)
  
  c(ed_value, permute_p, permute_p2, limit_p)
}

simulation_wrapper <- function(seed, case = "1-1", 
                               num_list = c(40, 60, 80, 100, 120), 
                               delta_list = c()
                               )
{
  pvalue_list <- list()
  if (length(delta_list) == 0) {
    for(i in 1:length(num_list)) {
      pvalue_list[[i]] = joint_independence_simulator(seed = seed, case = case, 
                                                      num = num_list[i])
    }
  } else {
    for(i in 1:length(delta_list)) {
      pvalue_list[[i]] = joint_independence_simulator(seed = seed, case = case, 
                                                      num = 80, 
                                                      delta = delta_list[i])
    }
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
sfInit(parallel = TRUE, cpus = 4)
sfLibrary(Ball)
sfLibrary(movMF)
sfLibrary(mvtnorm)
sfExportAll()

nrep <- 500

print(paste0("------------ Example 1-2-1 ------------"))
example_12 <- sfLapply(1:nrep, simulation_wrapper, case = "1-2", 
                       delta_list = c(0, 0.2, 0.4, 0.6, 0.8, 0.99))
res12_1 <- t(calculate_power(example_12)); print(res12_1)
save(res12_1, file="compare_indep_1_2_1.rda")

print(paste0("------------ Example 1-3-1 ------------"))
example_13 <- sfLapply(1:nrep, simulation_wrapper, case = "1-3", 
                       delta_list = c(0, 0.2, 0.4, 0.6, 0.8, 0.99))
res13_1 <- t(calculate_power(example_13)); print(res13_1)
save(res13_1, file="compare_indep_1_3_1.rda")

print(paste0("------------ Example 2-2-1 ------------"))
example_22 <- sfLapply(1:nrep, simulation_wrapper, case = "2-2", 
                       delta_list = c(0, 0.2, 0.4, 0.6, 0.8, 0.99))
res22_1 <- t(calculate_power(example_22)); print(res22_1)
save(res22_1, file="compare_indep_2_2_1.rda")

print(paste0("------------ Example 2-3-1 ------------"))
example_23 <- sfLapply(1:nrep, simulation_wrapper, case = "2-3", 
                       delta_list = c(0, 0.2, 0.4, 0.6, 0.8, 0.99))
res23_1 <- t(calculate_power(example_23)); print(res23_1)
save(res23_1, file="compare_indep_2_3_1.rda")

print(paste0("------------ Example 3-2-1 ------------"))
example_32 <- sfLapply(1:nrep, simulation_wrapper, case = "3-2", 
                       delta_list = c(0, 0.2, 0.4, 0.6, 0.8, 0.99))
res32_1 <- t(calculate_power(example_32)); print(res32_1)
save(res32_1, file="compare_indep_3_2_1.rda")

print(paste0("------------ Example 3-3-1 ------------"))
example_33 <- sfLapply(1:nrep, simulation_wrapper, case = "3-3", 
                       delta_list = c(0, 0.2, 0.4, 0.6, 0.8, 0.99))
res33_1 <- t(calculate_power(example_33)); print(res33_1)
save(res33_1, file="compare_indep_3_3_1.rda")

sfStop()


