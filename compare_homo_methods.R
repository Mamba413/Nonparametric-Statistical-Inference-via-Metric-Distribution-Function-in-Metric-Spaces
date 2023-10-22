rm(list = ls()); gc(reset = TRUE)
your_path <- "~"
setwd(your_path)
library(energy)
library(Ball)
your_python_path <- "/public/home/zhujin/miniconda3/envs/mdf/bin/python"   # change to your python path

homogeneity_simulator <- function(seed, case = "1-1", 
                                  methods = c("frechet", "ball", "energy"), 
                                  num = 3000, delta = NULL, p = 2) {
  set.seed(seed)
  each_n <- num / 3
  if (case == "1-1") {
    x <- matrix(rnorm(num * p), ncol = p)
  }
  if (case == "1-2") {
    delta <- ifelse(is.null(delta), 0.4, delta)
    x1 <- matrix(rnorm(each_n * p), ncol = p)
    x2 <- matrix(rnorm(each_n * p, mean = -delta), ncol = p)
    x3 <- matrix(rnorm(each_n * p, mean = delta), ncol = p)
    x <- rbind(x1, x2, x3)
  }
  if (case == "1-3") {
    delta <- ifelse(is.null(delta), 2.0, delta)
    x1 <- matrix(rt(each_n * p, df = 3.0 - delta), ncol = p)
    x2 <- matrix(rt(each_n * p, df = 3.0), ncol = p)
    x3 <- matrix(rt(each_n * p, df = 3.0 + delta), ncol = p)
    x <- rbind(x1, x2, x3)
  }
  if (case == "2-1") {
    x <- rmovMF(num, 2.5 * c(2, 0) / sqrt(4), alpha = 1)
  }
  if (case == "2-2") {
    delta <- tanh(delta * pi) * pi / 4
    x1 <- rmovMF(each_n, 2.5 * c(cos(pi/4 + delta), sin(pi/4 + delta)))
    x2 <- rmovMF(each_n, 2.5 * c(cos(pi/4), sin(pi/4)))
    x3 <- rmovMF(each_n, 2.5 * c(cos(pi/4 - delta), sin(pi/4 - delta)))
    x <- rbind(x1, x2, x3)
  }
  if (case == "2-3") {
    delta <- tanh(delta * pi) * pi / 4
    x1 <- rmovMF(each_n, 2.5 * rbind(c(cos(pi/4 + delta), sin(pi/4 + delta)), c(cos(pi/4 + delta), -sin(pi/4 + delta))), c(1, 1))
    x2 <- rmovMF(each_n, 2.5 * rbind(c(cos(pi/4), sin(pi/4)), c(cos(pi/4), -sin(pi/4))), c(1, 1))
    x3 <- rmovMF(each_n, 2.5 * rbind(c(cos(pi/4 - delta), sin(pi/4 - delta)), c(cos(pi/4 - delta), -sin(pi/4 - delta))), c(1, 1))
    x <- rbind(x1, x2, x3)
  }
  if (case == "3-1") {
    p <- 3
    cov_mat <- 0.0^abs(outer(1:p, 1:p, "-"))
    x <- rWishart(n = num, df = 8, Sigma = cov_mat)
  }
  if (case == "3-2") {
    delta <- ifelse(is.null(delta), 0.15, delta)
    p <- 3
    cov_mat1 <- (0.0)^abs(outer(1:p, 1:p, "-") != 0)
    cov_mat2 <- (-delta)^abs(outer(1:p, 1:p, "-") != 0)
    cov_mat3 <- (delta)^abs(outer(1:p, 1:p, "-") != 0)
    x1 <- rWishart(n = each_n, df = 8, Sigma = cov_mat1)
    x2 <- rWishart(n = each_n, df = 8, Sigma = cov_mat2)
    x3 <- rWishart(n = each_n, df = 8, Sigma = cov_mat3)
    # interestingly! this implementation works when using python
    x <- array(dim = c(p, p, 3 * each_n))
    for (i in 1:each_n) {
      x[, , i] <- x1[, , i]
      x[, , each_n + i] <- x2[, , i]
      x[, , 2 * each_n + i] <- x3[, , i]
    }
    # the following one-line implementation does not work
    # x <- abind::abind(x1, x2, x3)
  }
  if (case == "3-3") {
    delta <- ifelse(is.null(delta), 4.0, delta)
    p <- 3
    cov_mat <- (0.1)^abs(outer(1:p, 1:p, "-") != 0)
    x1 <- rWishart(n = each_n, df = 8.0, Sigma = cov_mat)
    x2 <- rWishart(n = each_n, df = 8.0 - delta, 
                   Sigma = (8.0 / (8.0 - delta)) * cov_mat)
    x3 <- rWishart(n = each_n, df = 8.0 + delta, 
                   Sigma = (8.0 / (8.0 + delta)) * cov_mat)
    x <- array(dim = c(p, p, 3 * each_n))
    for (i in 1:each_n) {
      x[, , i] <- x1[, , i]
      x[, , each_n + i] <- x2[, , i]
      x[, , 2 * each_n + i] <- x3[, , i]
    }
  }
  all_sample <- x
  
  dst_type <- strsplit(case, "-")[[1]][1]
  if (dst_type == "1") {
    x <- dist(x)
    data_type <- "multivariate"
  } else if (dst_type == "2") {
    x <- as.dist(Ball::nhdist(x))
    data_type <- "sphere"
  } else if (dst_type == "3") {
    x <- as.dist(CovDist(x, method = "C"))
    data_type <- "spd"
  }
  
  y <- matrix(0, nrow = num, ncol = 3)
  index <- rep(1:3, each = each_n)
  y[, 1] <- 1.0 * (index == 1)
  y[, 2] <- 1.0 * (index == 2)
  y[, 3] <- 1.0 * (index == 3)
  y <- dist(y)
  
  f_value <- ed_value <- permute_p <- limit_p <- 1.0
  
  if ("frechet" %in% methods) {
    group_size <- c(each_n, each_n, each_n)
    distance_type <- switch(dst_type,
                            "1" = "Euclidean", 
                            "2" = "Geodestic", 
                            "3" = "Cholesky")
    testing_data <- all_sample
    file_name <- sprintf("hypothesis_test_%s_%s_%s.RData", 
                         case, seed, each_n)
    save(testing_data, group_size, data_type, distance_type, 
         file = file_name)
    test_command <- "%s sim_test_one.py %s %s %s"
    test_command <- sprintf(test_command, 
                            your_python_path, case, seed, each_n)
    system(test_command)
    ftest_result_file_name <- "frechet_test_result_%s_%s_%s.txt"
    ftest_result_file_name <- sprintf(ftest_result_file_name, 
                                      case, seed, each_n)
    f_value <- read.table(ftest_result_file_name)
    f_value <- f_value[2, , drop = TRUE]
    file.remove(file_name)
    file.remove(ftest_result_file_name)
  }
  if ("ball" %in% methods) {
    permute_p <- bcov.test(x, y, distance = TRUE, 
                           weight = "rbf")[["p.value"]]
    # limit_p <- bcov.test(x, y, distance = TRUE, 
    #                      method = "limit", weight = "rbf")[["p.value"]]
  }
  if ("energy" %in% methods) {
    set.seed(1); 
    ed_value <- energy::eqdist.etest(x = x, sizes = rep(each_n, 3),
                                     distance = TRUE, R = 399)[["p.value"]]
  }
  
  c(f_value, ed_value, permute_p, limit_p)
}


simulation_wrapper <- function(seed, case = "1-1", 
                               num_list = c(30, 60, 90, 120, 150), 
                               delta_list = c(),
                               methods = c("frechet", "ball", "energy"), 
                               p = 3)
{
  pvalue_list <- list()
  if (length(delta_list) == 0) {
    for(i in 1:length(num_list)) {
      pvalue_list[[i]] = homogeneity_simulator(seed = seed, case = case, 
                                              methods = methods, 
                                              p = p, num = num_list[i])
    }
  } else {
    for(i in 1:length(delta_list)) {
      pvalue_list[[i]] = homogeneity_simulator(seed = seed, case = case, 
                                              methods = methods, 
                                              p = p, num = 90, 
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
sfInit(parallel = TRUE, cpus = 10)
sfLibrary(Ball)
sfLibrary(movMF)
sfLibrary(CovTools)
sfExportAll()

nrep <- 500

print(paste0("------------ Example 1-2 ------------"))
example_12 <- sfLapply(1:nrep, simulation_wrapper, case = "1-2", 
                       delta_list = c(0.0, 0.2, 0.4, 0.6, 0.8))
res12_1 <- t(calculate_power(example_12)); print(res12_1)
save(res12_1, file="compare_homo_1_2_1.rda")

print(paste0("------------ Example 1-3 ------------"))
example_13 <- sfLapply(1:nrep, simulation_wrapper, case = "1-3", 
                       delta_list = c(0.0, 0.2, 0.24, 0.36, 0.48) * 5)
res13_1 <- t(calculate_power(example_13)); print(res13_1)
save(res13_1, file="compare_homo_1_3_1.rda")

print(paste0("------------ Example 2-2 ------------"))
example_22 <- sfLapply(1:nrep, simulation_wrapper, case = "2-2", 
                       delta_list = seq(0, 0.6, length.out = 5))
res22_1 <- t(calculate_power(example_22)); print(res22_1)
save(res22_1, file="compare_homo_2_2_1.rda")

print(paste0("------------ Example 2-3 ------------"))
example_23 <- sfLapply(1:nrep, simulation_wrapper, case = "2-3", 
                       delta_list = seq(0, 0.7, length.out = 5))
res23_1 <- t(calculate_power(example_23)); print(res23_1)
save(res23_1, file="compare_homo_2_3_1.rda")

print(paste0("------------ Example 3-2 ------------"))
example_32 <- sfLapply(1:nrep, simulation_wrapper, case = "3-2", 
                       delta_list = c(0.0, 0.5, 1.0, 1.5, 2.0))
res32_1 <- t(calculate_power(example_32)); print(res32_1)
save(res32_1, file="compare_homo_3_2_1.rda")

print(paste0("------------ Example 3-3 ------------"))
example_33 <- sfLapply(1:nrep, simulation_wrapper, case = "3-3", 
                       delta_list = c(0.0, 1.0, 2.0, 3.0, 4.0))
res33_1 <- t(calculate_power(example_33)); print(res33_1)
save(res33_1, file="compare_homo_3_3_1.rda")

sfStop()
