####################
#Project script    # 
#cortex data       #
#PC Lenovo L13     #
#R version 4.1.1   #
#Windows x86_64-w64#
#Script by MEmam   #
####################



#Simulation of MOFA_data

#Disclaimer "THIS Part FOR SIMULATION METHOD IS THE ONLY PART THAT COME FROM MOFA+ OPEN SOURCE CODE FROM THEIR GITHUB"
make_example_data <- function(n_views=3, n_features=100, n_samples = 50, n_groups = 1,
                              n_factors = 5, likelihood = "gaussian",
                              lscales = 1, sample_cov = NULL, as.data.frame = FALSE) {
  
  # Sanity checks
  if (!all(likelihood %in% c("gaussian", "bernoulli", "poisson")))
    stop("Likelihood not implemented: Use either gaussian, bernoulli or poisson")
  
  if(length(lscales) == 1)
    lscales = rep(lscales, n_factors)
  if(!length(lscales) == n_factors)
    stop("Lengthscalces lscales need to be of length n_factors")
  if(all(lscales == 0)){
    sample_cov <- NULL
  }
  
  if (length(likelihood)==1) likelihood <- rep(likelihood, n_views) 
  if (!length(likelihood) == n_views) 
    stop("Likelihood needs to be a single string or matching the number of views!")
  
  if(!is.null(sample_cov)){
    if(sample_cov[1] == "equidistant") {
      sample_cov <- seq_len(n_samples)
    }
    if(is.null(dim(sample_cov))) sample_cov <- matrix(sample_cov, nrow = 1)
    if(ncol(sample_cov) != n_samples){
      stop("Number of columns in sample_cov must match number of samples n_samples.")
    }
    
    # Simulate covariance for factors
    Sigma = lapply(lscales, function(ls) {
      if(ls == 0) diag(1, n_samples)
      else (1) * exp(-as.matrix(stats::dist(t(sample_cov)))^2/(2*ls^2))
      # else (1-0.001) * exp(-as.matrix(stats::dist(t(sample_cov)))^2/(2*ls^2)) + diag(0.001, n_samples)
    })
    
    # simulate factors
    alpha_z <- NULL
    S_z <- lapply(seq_len(n_groups), function(vw) matrix(1, nrow=n_samples, ncol=n_factors))
    Z <-  vapply(seq_len(n_factors), function(fc) mvtnorm::rmvnorm(1, rep(0, n_samples), Sigma[[fc]]), numeric(n_samples))
    colnames(Z) <- paste0("simulated_factor_", 1:ncol(Z))
    Z <- lapply(seq_len(n_groups), function(gr) Z)
    sample_cov <- Reduce(cbind, lapply(seq_len(n_groups), function(gr) sample_cov))
  } else {
    # set sparsity for factors
    theta_z <- 0.5
    
    # set ARD prior for factors, each factor being active in at least one group
    alpha_z <- vapply(seq_len(n_factors), function(fc) {
      active_gw <- sample(seq_len(n_groups), 1)
      alpha_fc <- sample(c(1, 1000), n_groups, replace = TRUE)
      if(all(alpha_fc==1000)) alpha_fc[active_gw] <- 1
      alpha_fc
    }, numeric(n_groups))
    alpha_z <- matrix(alpha_z, nrow=n_factors, ncol=n_groups, byrow=TRUE)
    
    # simulate facors 
    S_z <- lapply(seq_len(n_groups), function(vw) matrix(rbinom(n_samples * n_factors, 1, theta_z),
                                                         nrow=n_samples, ncol=n_factors))
    Z <- lapply(seq_len(n_groups), function(vw) vapply(seq_len(n_factors), function(fc) rnorm(n_samples, 0, sqrt(1/alpha_z[fc,vw])), numeric(n_samples)))
  }
  
  # set sparsity for weights
  theta_w <- 0.5
  
  # set ARD prior, each factor being active in at least one view
  alpha_w <- vapply(seq_len(n_factors), function(fc) {
    active_vw <- sample(seq_len(n_views), 1)
    alpha_fc <- sample(c(1, 1000), n_views, replace = TRUE)
    if(all(alpha_fc==1000)) alpha_fc[active_vw] <- 1
    alpha_fc
  }, numeric(n_views))
  alpha_w <- matrix(alpha_w, nrow=n_factors, ncol=n_views, byrow=TRUE)
  
  # simulate weights 
  S_w <- lapply(seq_len(n_views), function(vw) matrix(rbinom(n_features*n_factors, 1, theta_w),
                                                      nrow=n_features, ncol=n_factors))
  W <- lapply(seq_len(n_views), function(vw) vapply(seq_len(n_factors), function(fc) rnorm(n_features, 0, sqrt(1/alpha_w[fc,vw])), numeric(n_features)))
  
  # set noise level (for gaussian likelihood)
  tau <- 10
  
  # pre-compute linear term and rbind groups
  mu <- lapply(seq_len(n_views), function(vw) lapply(seq_len(n_groups), function(gw)  (S_z[[gw]]*Z[[gw]]) %*% t(S_w[[vw]]*W[[vw]])))
  mu <- lapply(mu, function(l) Reduce(rbind, l))
  groups <- rep(paste("group",seq_len(n_groups), sep = "_"), each = n_samples)
  
  # simulate data according to the likelihood
  data <- lapply(seq_len(n_views), function(vw){
    lk <- likelihood[vw]
    if (lk == "gaussian"){
      dd <- t(mu[[vw]] + rnorm(length(mu[[vw]]),0,sqrt(1/tau)))
    }
    else if (lk == "poisson"){
      term <- log(1+exp(mu[[vw]]))
      dd <- t(apply(term, 2, function(tt) rpois(length(tt),tt)))
    }
    else if (lk == "bernoulli") {
      term <- 1/(1+exp(-mu[[vw]]))
      dd <- t(apply(term, 2, function(tt) rbinom(length(tt),1,tt)))
    }
    colnames(dd) <- paste0("sample_", seq_len(ncol(dd)))
    rownames(dd) <- paste0("feature_", seq_len(nrow(dd)),"_view", vw)
    dd
  })
  
  if(!is.null(sample_cov)) {
    colnames(sample_cov) <- colnames(data[[1]])
    rownames(sample_cov) <- paste0("covariate_", seq_len(nrow(sample_cov)))
  }
  
  names(data) <- paste0("view_", seq_len(n_views))
  
  if(as.data.frame){
    gr_df <- data.frame(group = groups, sample = colnames(data[[1]]))
    dat <- lapply(names(data), function(vw){
      tmp <- data[[vw]]
      df <- melt(tmp, varnames = c("feature", "sample"))
      df$view <- vw
      df
    })
    data <- bind_rows(dat)
    data <- dplyr::left_join(data, gr_df, by = "sample")
    
    sample_cov <- melt(sample_cov, varnames = c("covariate", "sample"))
  }
  return(list(data = data, groups = groups, alpha_w=alpha_w, alpha_z =alpha_z,
              lscales = lscales, sample_cov = sample_cov, Z = Z, W= W))
}


#Create_Mofa_data_object
MOFAexample <- make_example_data(n_views=2, n_features=1000, n_samples = 100, n_groups = 1,
                                 n_factors = 1, likelihood = "gaussian",
                                 lscales = 1, sample_cov = NULL, as.data.frame = FALSE)

omic1= as.data.frame(MOFAexample$data$view_1)
omic2= as.data.frame(MOFAexample$data$view_2)
Simulated_data_model = rbind(omic1, omic2)

#plot_ground_truth

#Sample_latent_variable(score)

sample_score = as.data.frame(MOFAexample$Z)
sample_score $index_s = rownames(sample_score)
colnames(sample_score )= c("sample_value", "sample_index")
plot(x=sample_score$sample_index, y= sample_score $sample_value, xlab="Sample_index", ylab="Latent_score")


#Features_loading(weight)

feature_weight= as.data.frame(MOFAexample$W)
omic1_weight = feature_weight[1]
colnames(omic1_weight) = "weight"
omic2_weight = feature_weight[2]
colnames(omic2_weight) = "weight"
feature_weight = rbind(omic1_weight, omic2_weight)
feature_weight$label = 1
feature_weight$label[1: 1000] = "omic1"
feature_weight$label[1001: 2000] = "omic2"
ggplot(feature_weight, aes(x= c(1:2000), y=weight, color=label)) + scale_colour_manual(values = c("lightcoral", "cyan3")) + 
  labs(y= "Feature weights", x = "Feature index") + geom_point(alpha=0.5, size=1) + theme_bw() + ggtitle("Ground Truth simulated dataset") + 
  theme(plot.title = element_text(hjust = 0.5, face = 'bold'), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = 'none', legend.background = element_rect(fill="white", size=0.6, linetype="solid", colour ="lemonchiffon4"))


#bulid_model_with_FABIA

library("fabia")

#Function_request_data_frame_&_number_of_factor_create_fabia_object
Fabia_function <- function(merged_fabia_data, num){
  set.seed(123)
  
  resFabia = fabia(as.matrix(merged_fabia_data), p=num, alpha=0, 1, cyc=1000, spl=0.5, spz=0.5, random=1.0, center=2, norm=2, scale=0.0, lap=1.0, nL=1)
  extractPlot(resFabia, which = 5)
  #rb1 contains the information about all the biclusters
  rb1 = extractBic(resFabia,thresZ=0.5,thresL=NULL)
  return(resFabia)
}

#create_feature_loading_data_frame_function
Fabia_loading <- function(Fabia_object, BC_num){
  
  loadings10_FABIA = Fabia_object@L[,BC_num]
  BC1_df_loading = as.data.frame(loadings10_FABIA)
  return(BC1_df_loading)
}
#create_sample_score_function
Fabia_score <- function(Fabia_object, BC_num){
  
  scores10_FABIA1 = Fabia_object@Z[BC_num,] 
  return(scores10_FABIA1)
}

#Create_object
Fabia_object = Fabia_function(merged_fabia_data = Simulated_data_model, num = 1)
#create_sample_score_object
Fabia_score_BC1 = Fabia_score(Fabia_object = Fabia_object, BC_num = 1)
#create_feature_loading_data_frame_object
BC1_df_loading = Fabia_loading(Fabia_object = Fabia_object, BC_num = 1)

#MOFA
library(tidyverse)
library(ggplot2)
library(MOFA2)
library(data.table)
library(ggplot2)

MOFA_function <- function(merged_MOFA_data, num){
  MOFAobject <- create_mofa(merged_MOFA_data$data)
  data_opts <- get_default_data_options(MOFAobject)
  data_opts$scale_views <- TRUE
  data_opts$scale_groups <- TRUE
  model_opts <- get_default_model_options(MOFAobject)
  train_opts <- get_default_training_options(MOFAobject)
  train_opts$maxiter <-25
  train_opts$convergence_mode <- "slow"
  train_opts$seed <- 42
  model_opts$num_factors= num
  
  MOFAobject <- prepare_mofa(
    object = MOFAobject,
    data_options = data_opts,
    model_options = model_opts,
    training_options = train_opts
  )
  MOFAobject.trained <- run_mofa(MOFAobject, outfile)
  model<-MOFAobject.trained
  plot_data_overview(model)
  return(model)
  
}
model = MOFA_function(MOFAexample, num = 1)
plot_variance_explained(model, x="group", y="factor")

#create_MOFA_feature_loading_data_frame_function
MOFA_score <- function(model, Factor_num){
  
  factors <- get_factors(model, factors =Factor_num, as.data.frame = T)
  
  factor1 = factors[1:100, 1:3]
  
  factor1_df_score = as.data.frame(factor1$value)
  rownames(factor1_df_score)= factor1$sample
  return(factor1_df_score)
}

MOFA_loading <- function(model, Factor_num){
  
  weights <- get_weights(model, factors = Factor_num, as.data.frame = T)
  weights$weight= weights$value
  df1_weights_1 = weights[c("feature","weight")] 
  colnames(df1_weights_1)= c("features","weights")
  weights1= df1_weights_1
  factor1_df_loading= as.data.frame(weights1$weights)
  
  substrRight <- function(x, n){
    substr(x, nchar(x)-n+1, nchar(x))
  }
  for (i in 1:dim(df2_weights_1)[1]) {
    if(substrRight(as.character(df2_weights_1$features[i]), 5)=="view1"){
      df2_weights_1$feature[i] = "omic1"}
    if(substrRight(as.character(df2_weights_1$features[i]), 5)!="view1"){
      df2_weights_1$feature[i] = "omic2"}
  }
  plot_df1_weights_1 <- ggplot(df2_weights_1, aes(x=feature_index, y=weights, color=feature)) + scale_colour_manual(values = c("lightcoral", "cyan3")) + 
    labs(y= "Feature weights", x = "Feature index") + geom_point(alpha=0.5, size=1) + theme_bw() + ggtitle("Feature weight") + 
    theme(plot.title = element_text(hjust = 0.5, face = 'bold'), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = 'none', legend.background = element_rect(fill="white", size=0.6, linetype="solid", colour ="lemonchiffon4")) 
  print(plot_df1_weights_1)
  return(factor1_df_loading)
}


factor1_df_score= MOFA_score(model = model, Factor_num = 1)
factor1_df_loading= MOFA_loading(model = model, Factor_num = 1)



library("GFA")

#Function_request_data_frame_&_number_of_factor_create_GFA_object
GFA_function <- function(merged_GFA_data, num){
  set.seed(123)
  model_option <- getDefaultOpts()                         
  GFA_object <- gfa(t(merged_GFA_data), K= num, opts=model_option) 
  normalized_data <- normalizeData(merged_GFA_data, type="center")
  visulization <- visualizeComponents(GFA_object, merged_GFA_data, normalized_data )
  return(GFA_object)
}

#Create_sample_score_GFA_function
GFA_score <- function(GFA_object, BC_num){
  
  GFA_scores = as.data.frame(GFA_object$X)
  g_factor1_df_score= data.frame(GFA_scores$V+BC_num)
  rownames(g_factor1_df_score)= rownames(GFA_scores) 
  return(g_factor1_df_score)
}

#Create_feature_loading_data_frame_GFA_function
GFA_loading <- function(GFA_object, BC_num){
  
  GFA_weight= as.data.frame(GFA_object$W)
  g_factor1_df_loading= data.frame(GFA_weight$V1+BC_num)
  rownames(g_factor1_df_loading)= rownames(GFA_weight)
  return(g_factor1_df_loading)
}
GFA_object= GFA_function(merged_GFA_data = list(t(omic1), t(omic2)), num = 1)
g_factor1_df_score= GFA_score(GFA_object= GFA_object, BC_num= 1)
g_factor1_df_loading = GFA_loading(GFA_object=GFA_object , BC_num= 1)


#Multiple_factorAnalysis
library(factoextra)
library(FactoMineR)

#Function_request_data_frame_&_number_of_factor_create_MFA_object
MFA_function <- function(merged_MFA_data, group,  num){
  set.seed(123)
  resMFA_scale_after = MFA(as.data.frame(merged_MFA_data),group=group,ncp = num)
  return(resMFA_scale_after)
}


#Create_sample_score_MFA_function
MFA_score <- function(MFA_object, num){
  scores1_scaled_after = data.frame(MFA_object$ind$coord[,num])
  colnames(scores1_scaled_after) = "MFA_score"
  return(scores1_scaled_after)
}

#Create_feature_loading_data_frame_MFA_function
MFA_loading <- function(MFA_object, num){
  fviz_contrib(MFA_object,choice="partial.axes", axes = num, top =10,palette= "jco")
  loadings1_scaled_after = data.frame(MFA_object$quanti.var$coord[,num])
  colnames(loadings1_scaled_after) = "MFA_loading"
  return(loadings1_scaled_after)
}
MFA_object = MFA_function(merged_MFA_data= cbind(t(omic1),t(omic2)), group = c(ncol(t((omic1))),ncol(t(omic2))), num = 2)
scores1_scaled_after= MFA_score(MFA_object= MFA_object, num = 1)
loadings1_scaled_after= MFA_loading(MFA_object= MFA_object, num = 1)
updated= as.data.frame(cbind(feature_weight$weight, factor1_df_loading,BC1_df_loading, loadings1_scaled_after, g_factor1_df_loading))

colnames(updated)= c("True_w", "MOFA_w","FABIA_w",  "MFA_w", "GFA_w")

updated_s= as.data.frame(cbind(sample_score$sample_value, scores10_FABIA1, factor1_df_score, scores1_scaled_after, g_factor1_df_score))
colnames(updated_s)= c("True_s", "FABIA_s", "MOFA_s", "MFA_s", "GFA_s")


library(ggpubr)
ggscatter(updated_s, x = "MOFA_s", y = "MFA_s", conf.int = TRUE, 
          ,add = "reg.line", cor.coef = TRUE, cor.method = "pearson",
          xlab = "MOFA_Factor1", ylab = "True")


#Take_quantile_according_to_each_tool_threshold_values


#Ground_truth_threshold_values
updated_Tru_s_80 = updated_s[which(updated_s$True_s > quantile(updated_s$True_s, probs = 0.80)),]
updated_Tru_w_80= updated[which(updated$True_w > quantile(updated$True_w, probs = 0.80)),]
updated_Tru_s_20 = updated_s[which(updated_s$True_s < quantile(updated_s$True_s, probs = 0.20)),]
updated_Tru_w_20= updated[which(updated$True_w < quantile(updated$True_w, probs = 0.20)),]
updated_Tru_top_s = rbind(updated_Tru_s_20, updated_Tru_s_80)
updated_Tru_top_w = rbind(updated_Tru_w_20, updated_Tru_w_80)

#MOFA
updated_MOFA_s_80 = updated_s[which(updated_s$MOFA_s > quantile(updated_s$MOFA_s, probs = 0.80)),]
updated_MOFA_w_80= updated[which(updated$MOFA_w > quantile(updated$MOFA_w, probs = 0.80)),]
updated_MOFA_s_20 = updated_s[which(updated_s$MOFA_s < quantile(updated_s$MOFA_s, probs = 0.20)),]
updated_MOFA_w_20= updated[which(updated$MOFA_w < quantile(updated$MOFA_w, probs = 0.20)),]
updated_MOFA_top_s = rbind(updated_MOFA_s_20, updated_MOFA_s_80)
updated_MOFA_top_w = rbind(updated_MOFA_w_20, updated_MOFA_w_80)

#MFA
updated_MFA_s_80 = updated_s[which(updated_s$MFA_s > quantile(updated_s$MFA_s, probs = 0.80)),]
updated_MFA_w_80= updated[which(updated$MFA_w > quantile(updated$MFA_w, probs = 0.80)),]
updated_MFA_s_20 = updated_s[which(updated_s$MFA_s < quantile(updated_s$MFA_s, probs = 0.20)),]
updated_MFA_w_20= updated[which(updated$MFA_w < quantile(updated$MFA_w, probs = 0.20)),]
updated_MFA_top_s = rbind(updated_MFA_s_20, updated_MFA_s_80)
updated_MFA_top_w = rbind(updated_MFA_w_20, updated_MFA_w_80)

#GFA
updated_GFA_s_80 = updated_s[which(updated_s$GFA_s > quantile(updated_s$GFA_s, probs = 0.80)),]
updated_GFA_w_80= updated[which(updated$GFA_w > quantile(updated$GFA_w, probs = 0.80)),]
updated_GFA_s_20 = updated_s[which(updated_s$GFA_s < quantile(updated_s$GFA_s, probs = 0.20)),]
updated_GFA_w_20= updated[which(updated$GFA_w < quantile(updated$GFA_w, probs = 0.20)),]
updated_GFA_top_s = rbind(updated_GFA_s_20, updated_GFA_s_80)
updated_GFA_top_w = rbind(updated_GFA_w_20, updated_GFA_w_80)

#FABIA
updated_FABIA_s = updated_s[which(updated_s$FABIA_s > quantile(updated_s$FABIA_s, probs = 0.80)),]
updated_FABIA_w= updated[which(updated$FABIA_w > quantile(updated$FABIA_w, probs = 0.80)),]

updated_FABIA_s_80 = updated_s[which(updated_s$FABIA_s > quantile(updated_s$FABIA_s, probs = 0.80)),]
updated_FABIA_w_80= updated[which(updated$FABIA_w > quantile(updated$FABIA_w, probs = 0.80)),]
updated_FABIA_s_20 = updated_s[which(updated_s$FABIA_s < quantile(updated_s$FABIA_s, probs = 0.20)),]
updated_FABIA_w_20= updated[which(updated$FABIA_w < quantile(updated$FABIA_w, probs = 0.20)),]
updated_FABIA_top_s = rbind(updated_FABIA_s_20, updated_FABIA_s_80)
updated_FABIA_top_w = rbind(updated_FABIA_w_20, updated_FABIA_w_80)

ggscatter(updated_MOFA_top_s, x = "MOFA_w", y = "True_w", conf.int = TRUE, 
          ,add = "reg.line", cor.coef = TRUE, cor.method = "pearson",
          xlab = "MOFA_Factor1", ylab = "True")
library(corrplot)
corrplot(cor(updated_MOFA_top_s),
         method = "number",       
         order = "hclust",         # Ordering method of the matrix
         hclust.method = "ward.D", # If order = "hclust", is the cluster method to be used
         addrect = 2,              # If order = "hclust", number of cluster rectangles
         rect.col = 3,             # Color of the rectangles
         rect.lwd = 3) 
