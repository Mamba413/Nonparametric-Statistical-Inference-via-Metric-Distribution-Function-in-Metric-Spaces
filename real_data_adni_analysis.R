rm(list = ls()); gc(reset = TRUE)
your_path <- "~"
setwd(your_path)
library(fda.usc)
library(fda)
library(Ball)
library(multivariance)
library(energy)
library(lokern)
library(fdapace)

#####################################################
################# clinical variable #################
#####################################################
load("clinical.dat")
scoredata <- as.data.frame(scoredata)
colnames(scoredata) <- c("time", "status", "id", "gender", "handness", 
                         "widowed", "divorced", "neverMarried", 
                         "educationLength", "retirement", "age", 
                         "APOE2", "APOE3", "APOE4", 
                         "ADASScore")
dim(scoredata)
anyNA(scoredata)

dx_covariate <- list()
dx_covariate[["gender"]] <- dist(scoredata[["gender"]])
dx_covariate[["handness"]] <- dist(scoredata[["handness"]])
dx_covariate[["maritalStatus"]] <- dist(scoredata[, c("widowed", "divorced", "neverMarried")])
dx_covariate[["educationLength"]] <- dist(scoredata[["educationLength"]])
dx_covariate[["retirement"]] <- dist(scoredata[["retirement"]])
dx_covariate[["age"]] <- dist(scoredata[["age"]])
dx_covariate[["genetic1"]] <- dist(scoredata[["APOE2"]])
dx_covariate[["genetic2"]] <- dist(scoredata[["APOE3"]])
dx_covariate[["genetic3"]] <- dist(scoredata[["APOE4"]])
dx_covariate[["ADASScore"]] <- dist(scoredata[["ADASScore"]])

##############################################
################# Hippcampus #################
##############################################
load("hippocampussurface.dat")
hippcampus_l <- t(cleandata[1:15000, ])
hippcampus_r <- t(cleandata[15001:30000, ])

#############################################################################
######################### Plot hippocampus example ##########################
#############################################################################
library(reshape2)
library(ggplot2)
pdat <- cbind.data.frame(rbind(hippcampus_l[1, , drop = TRUE],
                               hippcampus_r[1, , drop = TRUE]),
                         "side" = c("Left", "Right"))
x_lab <- ""
pdat <- melt(pdat, id.vars = "side", variable.name = "x", value.name = "distance")
pdat[["x"]] <- as.numeric(as.character(pdat[["x"]]))
p <- ggplot(pdat, aes(x = x, y = distance)) + 
  facet_wrap(side ~ ., ncol = 1, strip.position = "right") + 
  geom_line(size = 0.1) + 
  scale_y_continuous() + 
  scale_x_continuous(expand = c(0.01, 0.01)) + 
  ylab("Radial distance") + 
  xlab(x_lab) + 
  theme_bw() + 
  theme(
    axis.title.x = theme_x_title,
    strip.text = element_text(face = "bold", size = 10),
    panel.grid = element_blank(), 
  )
p
ggsave(p, filename = "hippocampus_curves_example.jpg", width = 9.6, height = 3.66)

###################################################################
######################### Functional PCA ##########################
###################################################################
smooth_method <- "fpca"
fpca_smooth <- function(hippcampus, FVEthreshold = 0.99, max_K = 34) {
  num <- nrow(hippcampus)
  M <- ncol(hippcampus)
  s <- 1:ncol(hippcampus)
  L <- MakeFPCAInputs(IDs = rep(1:num, each = M), 
                      tVec = rep(s, num), t(hippcampus))
  FPCAdense <- FPCA(L$Ly, L$Lt, optns = list(plot = TRUE, methodMuCovEst = 'smooth', 
                                             maxK = max_K, FVEthreshold = FVEthreshold))
  select_K <- FPCAdense[["selectK"]]
  print(select_K)
  smooth_hippcampus <- fitted(FPCAdense, K = select_K)
  smooth_hippcampus
}
smooth_l <- fpca_smooth(hippcampus_l)
smooth_r <- fpca_smooth(hippcampus_r)
# save(smooth_l, smooth_r, file = "fpca_smooth_99.rda")
# load("fpca_smooth_99.rda")   

############################
##### compute distance #####
############################
hippcampus_l <- fdata(smooth_l)
hippcampus_r <- fdata(smooth_r)

## L2 distance
distance_method <- "L2"
dx_hippcampus_l <- metric.lp(hippcampus_l)
dx_hippcampus_r <- metric.lp(hippcampus_r)
dx_hippcampus_L2 <- sqrt(dx_hippcampus_l^2 + dx_hippcampus_r^2)

file_name <- sprintf("covariate_dist_%s_%s.rda", smooth_method, distance_method)
# save(dx_hippcampus_l, dx_hippcampus_r, dx_hippcampus, dx_covariate, scoredata, file = file_name)
# load(file_name)

#################################################################
######################### Testing ###############################
#################################################################
test_num <- length(names(dx_covariate))
tm_value <- numeric(test_num)
ma_value <- numeric(test_num)
ma_value2 <- numeric(test_num)
ma_value_spectrum <- numeric(test_num)
ma_value_spectrum2 <- numeric(test_num)
dx_hippcampus_list <- list("L2" = dx_hippcampus_L2)
dx_type <- names(dx_hippcampus_list)
for (type in dx_type) {
  cat("------------------------------------------\n")
  cat("----------------", type, "----------------\n")
  cat("------------------------------------------\n")
  dx_hippcampus <- dx_hippcampus_list[[type]]
  
  i <- 1
  for (cov_name in names(dx_covariate)) {
    cat("------------------------------------------\n")
    cat("--------------", cov_name, "--------------\n")
    cat("------------------------------------------\n")
    dx_cov_name <- dx_covariate[[cov_name]]
    set.seed(1234);
    tm_value[i] <- dcov.test(as.dist(dx_hippcampus), as.dist(dx_cov_name), R = 999)[["p.value"]]
    ma_value[i] <- bcov.test(dx_hippcampus, dx_cov_name, distance = TRUE,
                             num.permutations = 999, seed = 1234)[["p.value"]]
    dx_list <- list(as.vector(as.dist(dx_hippcampus)), 
                    as.vector(as.dist(dx_cov_name)))
    ma_value_spectrum[i] <- KBCovLimit(dx_list, n = nrow(dx_hippcampus))[["p.value"]]
    
    ma_value2[i] <- KBCov_permuted(dx_list, n = nrow(dx_hippcampus), 
                                   num.permutations = 999, 
                                   seed = 1234)[["p.value"]]
    ma_value_spectrum2[i] <- KBCov(dx_list, n = nrow(dx_hippcampus))[["p.value"]]
    i <- i + 1
  }
  
  print("Raw p-value: ")
  print(round(tm_value, 3))
  print(round(ma_value, 3))
  print(round(ma_value_spectrum, 3))
  print(round(ma_value2, 3))
  print(round(ma_value_spectrum2, 3))
  print("Adjusted p-value: ")
  print(round(p.adjust(tm_value), 3))
  print(round(p.adjust(ma_value), 3))
  print(round(p.adjust(ma_value_spectrum), 3))
  print(round(p.adjust(ma_value2), 3))
  print(round(p.adjust(ma_value_spectrum2), 3))
}
