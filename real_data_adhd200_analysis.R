library(Ball)
library(energy)
rm(list = ls()); gc(reset = TRUE)
path <- "~/mdf/code/"
setwd(path)
tmp_file_path <- "ADHD200-HO"

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

############################################
############# download dataset #############
############################################
source("adhd200_download.R")

phenotype <- read.csv("adhd200_preprocessed_phenotypics.tsv", sep = "\t")
## there are 162 observations because some subjects are filtered according to frequency
analysis_index <- match(as.numeric(file_id), phenotype[["ScanDir.ID"]])
demograph_all <- phenotype[analysis_index, ]
demograph_all[demograph_all == "-999"] <- NA
demograph_all[demograph_all == "N/A"] <- NA
demograph_all[["Handedness"]] <- as.numeric(as.character(demograph_all[["Handedness"]]))
disease_group <- ifelse(demograph_all[["DX"]] == 0, 1, 2)
demograph_all[["disease_group"]] <- disease_group
save(roi_ts, demograph_all, file = "ADHD200-HO/ts_covariate_data.rda")

###################### Pre-process ######################
cvglasso_batch <- function(ts) {
  res <- tryCatch(
    {
      CVglasso::CVglasso(X = ts, nlam = 10, K = 5, trace = "none", diagonal = TRUE)
    },
    error = function(cond) {
      list()
    }
  )
  res[c("Sigma", "Omega")]
}
library(snowfall)
sfInit(parallel = TRUE, cpus = 6)
sfLibrary(CVglasso)
sfExport("cvglasso_batch")
result <- sfLapply(roi_ts, cvglasso_batch)
sfStop()
# save(result, file = "ADHD200-HO/ts_glasso.rda")
# load("ADHD200-HO/ts_glasso.rda")

rm(list = setdiff(ls(), c("demograph_all", "result"))); gc(reset = TRUE)
graph_list <- lapply(result, function(x) {
  ((x[[2]]) + t(x[[2]])) / 2
})
null_graph <- which(sapply(graph_list, is.null))
if (length(null_graph) != 0) {
  graph_list <- graph_list[-null_graph]
  analysis_phenotype <- analysis_phenotype[-null_graph, ]
}
graph_num <- length(graph_list)
graph_array <- array(dim = c(dim(graph_list[[1]]), graph_num))
for (i in 1:graph_num) {
  graph_array[, , i] <- graph_list[[i]]
}
# save(graph_array, demograph_all, file = "ADHD200-HO/graph_covariate_data.rda")
# load("ADHD200-HO/graph_covariate_data.rda")

############################################
# Compute metric
############################################
dx <- measure.Choleksy.3d(graph_array)

############################################
double_normalized_distance <- function(dx) {
  if (!("matrix" %in% class(dx))) {
    dx <- as.matrix(dx)
  }
  size <- nrow(dx)
  C <- diag(size) - matrix(data = 1, nrow = size, ncol = size) / size
  -C %*% dx %*% C
}

distance_multivariate <- function(data_list, R = 199) {
  data_list <- lapply(data_list, double_normalized_distance)
  tm <- multivariance::total.multivariance(x = data_list)
  tm.resampled <- numeric(R)
  for (r in 1:R) {
    set.seed(r)
    tm.resampled[r] <- multivariance::total.multivariance(
      multivariance::sample.cdms(data_list))
  }
  dcov_pvalue <- (sum(tm.resampled >= tm) + 1) / (R + 1)
  dcov_pvalue
}

############################### Dependence Testing ######################################
summary(demograph_all)
demograph_all[["DX"]] <- NULL
table(demograph_all[["Secondary.Dx"]])

summary(demograph_all)
seed <- 1234

na_index <- which(is.na(demograph_all[["Gender"]]) | is.na(demograph_all[["Handedness"]]))
bcov.test(demograph_all[["Handedness"]][-na_index], demograph_all[["Gender"]][-na_index], num.permutations = 999)
set.seed(seed); dcov.test(demograph_all[["Handedness"]][-na_index], demograph_all[["Gender"]][-na_index], R = 999)
dst_list2 <- list(as.vector(dist(demograph_all[["Handedness"]][-na_index])), 
                  as.vector(dist(demograph_all[["Gender"]][-na_index])))
KBCovLimit(dst_list2, n = length(demograph_all[["Gender"]][-na_index]))[["p.value"]]
KBCov(dst_list2, n = length(demograph_all[["Gender"]][-na_index]))[["p.value"]]
KBCov_permuted(dst_list2, n = length(demograph_all[["Gender"]][-na_index]), 
               num.permutations = 999, seed = seed)[["p.value"]]

dxgh <- as.matrix(dist(demograph_all[-na_index, c("Gender", "Handedness")]))
dxfc <- dx[-na_index, -na_index]
bcov.test(dxgh, dxfc, num.permutations = 999, distance = TRUE)
set.seed(seed); dcov.test(as.dist(dxgh), as.dist(dxfc), R = 999)
dst_list3 <- list(dxgh[lower.tri(dxgh)], dxfc[lower.tri(dxfc)])
KBCovLimit(dst_list3, n = length(demograph_all[["Gender"]][-na_index]))[["p.value"]]
KBCov(dst_list3, n = length(demograph_all[["Gender"]][-na_index]))[["p.value"]]
KBCov_permuted(dst_list3, n = length(demograph_all[["Gender"]][-na_index]), 
               num.permutations = 999, seed = seed)[["p.value"]]

dxh <- as.matrix(dist(demograph_all[-na_index, c("Handedness")]))
dxfcg <- sqrt((dx[-na_index, -na_index])^2 + 
                (as.matrix(dist(demograph_all[["Gender"]][-na_index])))^2)
bcov.test(dxh, dxfcg, num.permutations = 999, distance = TRUE)
set.seed(seed); dcov.test(as.dist(dxh), as.dist(dxfcg), R = 999)
dst_list4 <- list(dxh[lower.tri(dxh)], dxfcg[lower.tri(dxfcg)])
KBCovLimit(dst_list4, n = length(demograph_all[["Gender"]][-na_index]))[["p.value"]]
KBCov(dst_list4, n = length(demograph_all[["Gender"]][-na_index]))[["p.value"]]
KBCov_permuted(dst_list4, n = length(demograph_all[["Gender"]][-na_index]), 
               num.permutations = 999, seed = seed)[["p.value"]]

dxg <- as.matrix(dist(demograph_all[-na_index, c("Gender")]))
dxfch <- sqrt((dx[-na_index, -na_index])^2 + 
                (as.matrix(dist(demograph_all[["Handedness"]][-na_index])))^2)
bcov.test(dxg, dxfch, num.permutations = 999, distance = TRUE)
set.seed(seed); dcov.test(as.dist(dxg), as.dist(dxfch), R = 999)
dst_list5 <- list(dxg[lower.tri(dxg)], dxfch[lower.tri(dxfch)])
KBCovLimit(dst_list5, n = length(demograph_all[["Gender"]][-na_index]))[["p.value"]]
KBCov(dst_list5, n = length(demograph_all[["Gender"]][-na_index]))[["p.value"]]
KBCov_permuted(dst_list5, n = length(demograph_all[["Gender"]][-na_index]), 
               num.permutations = 999, seed = seed)[["p.value"]]


res1 <- cbind(c(0.057, 0.02, 0.727), 
             c(0.139, 0.008, 0.876), 
             c(0.112117496469696, 0.00538956535891366, 0.9887666), 
             c(0.124407557204615, 0.00742301292891823, 0.9730179), 
             c(0.103, 0.006, 0.884))
t(apply(res1, 2, round, digits = 3))
apply(t(apply(res1, 2, p.adjust)), 2, round, digits = 3)

res2 <- cbind(c(0.018, 0.321, 0.037),
              c(0.005, 0.701, 0.017),
              c(0.005951589, 0.6417801, 0.01730777), 
              c(0.009864958, 0.6557691, 0.02238525),
              c(0.006, 0.721, 0.02))
t(apply(res2, 2, round, digits = 3))
apply(t(apply(res2, 2, p.adjust)), 2, round, digits = 3)
