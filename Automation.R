####################
#Project script    # 
#multi_omics       #
#PC Lenovo L13     #
#R version 4.1.1   #
#Windows x86_64-w64#
#Script by MEmam   #
####################

#Simulated_function_to_simulate_multi_omics
simulate_data<- function(n_s, n_d){
  library("stringr") 
  #Sampled_Signal_start_and_end
  s_sig_s= ceiling(n_s/sample(3:5, 1))
  s_sig_e= ceiling(n_s/ sample(1:2, 1))
  alpha<-rnorm(n_s,0,0.05)
  alpha[s_sig_s:s_sig_e]<-rnorm(length(s_sig_s:s_sig_e),3,0.05)
  #FeaturesD_Signal_start_and_end
  d_sig_s= ceiling(n_d/sample(3:5, 1))
  d_sig_e= ceiling(n_d/ sample(1:2, 1))
  beta<-rnorm(n_d,0,0.05)
  beta[d_sig_s:d_sig_e]<-rnorm(length(d_sig_s:d_sig_e),5,0.05)
  data.1<-alpha%*%t(beta) #signal
  dim(data.1)
  eps<-rnorm(n_s*n_d,3,4) #noise
  data.1a<-data.1+matrix(eps,n_s,n_d) #signal + noise
  data.2<-data.1a[,c(1:(n_d/2))]
  data.3<-data.1a[,c(((n_d/2)+1):n_d)]
  omic1_so = as.data.frame(t(data.2))
  rownames(omic1_so) = paste0('Feature_omic1_', rownames(omic1_so))
  colnames(omic1_so) = str_replace_all(colnames(omic1_so), 'V', 'sample_')
  omic2_so = as.data.frame(t(data.3))
  rownames(omic2_so) = paste0('Feature_omic2_', rownames(omic2_so))
  colnames(omic2_so) = str_replace_all(colnames(omic2_so), 'V', 'sample_')
  results <- list(alpha=alpha,beta= beta, omic1_so = omic1_so, omic2_so = omic2_so)
  return(results)
}
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
#MOFA
library(tidyverse)
library(ggplot2)
library(MOFA2)
library(data.table)
library(ggplot2)
MOFA_function <- function(merged_MOFA_data, num){
  MOFAobject <- create_mofa(merged_MOFA_data)
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
  return(factor1_df_loading)
}
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
Integromics<- function(n_s, n_d){
  simulate_data = simulate_data(n_s,n_d)
  alpha <- simulate_data$alpha
  beta <- simulate_data$beta
  omic1_so <- as.data.frame(simulate_data$omic1_so)
  omic2_so <- as.data.frame(simulate_data$omic2_so)
  #creat_Mofa_data_structure
  data <- make_example_data(
    n_views = 2,
    n_samples = 2,
    n_features = 2,
    n_factors = 2
  )[[1]]
  names(data)<-c("omic1","omic2")
  lapply(data,dim)
  a<-data.matrix(omic1_so, rownames.force = NA)
  b<-data.matrix(omic2_so, rownames.force = NA)
  data[["omic1"]]<-a
  data[["omic2"]]<-b
  N = ncol(data[[1]])
  Simulated_data_model = as.data.frame(rbind(omic1_so, omic2_so))
  #Create_object
  Fabia_object = Fabia_function(merged_fabia_data = Simulated_data_model, num = 1)
  #create_sample_score_object
  Fabia_score_BC1 = Fabia_score(Fabia_object = Fabia_object, BC_num = 1)
  #create_feature_loading_data_frame_object
  BC1_df_loading = Fabia_loading(Fabia_object = Fabia_object, BC_num = 1)
  #data <- simulate_data(n_s, n_d)$data
  model = MOFA_function(data, num = 1)
  plot_variance_explained(model, x="group", y="factor")
  factor1_df_score= MOFA_score(model = model, Factor_num = 1)
  factor1_df_loading= MOFA_loading(model = model, Factor_num = 1)
  GFA_object= GFA_function(merged_GFA_data = list(t(omic1_so), t(omic2_so)), num = 1)
  g_factor1_df_score= GFA_score(GFA_object= GFA_object, BC_num= 1)
  g_factor1_df_loading = GFA_loading(GFA_object=GFA_object , BC_num= 1)
  MFA_object = MFA_function(merged_MFA_data= cbind(t(omic1_so),t(omic2_so)), group = c(ncol(t((omic1_so))),ncol(t(omic2_so))), num = 2)
  scores1_scaled_after= MFA_score(MFA_object= MFA_object, num = 1)
  loadings1_scaled_after= MFA_loading(MFA_object= MFA_object, num = 1)
  fviz_contrib(MFA_object,choice="partial.axes", axes = 1, top =10,palette= "jco")
  correlation_w= as.data.frame(cbind(beta,BC1_df_loading,factor1_df_loading, loadings1_scaled_after, g_factor1_df_loading))
  colnames(correlation_w)= c("True_w", "FABIA_w", "MOFA_w",  "MFA_w", "GFA_w")
  correlation_w= cor(correlation_w$"True_w", correlation_w)  
  correlation_s= as.data.frame(cbind(alpha, Fabia_score_BC1,factor1_df_score, scores1_scaled_after, g_factor1_df_score))
  colnames(correlation_s)= c("True_s", "FABIA_s", "MOFA_s", "MFA_s", "GFA_s")
  correlation_s= cor(correlation_s$"True_s", correlation_s)  
  correlation <- list(correlation_s=correlation_s,correlation_w = correlation_w)
  return(correlation)
}
#First_correlation_var
first_cor = Integromics(100, 2000)
#Internalize_iterator
i <- 1
#new_w_is_a_var_to_accumulate_wight_from_each_iteration
new_w <- as.data.frame(first_cor$correlation_w) 
#new_s_is_a_var_to_accumulate_score_from_each_iteration
new_s <- as.data.frame(first_cor$correlation_s) 
while (i < 4) {
  cor_i = Integromics(100, 2000)
  df_s=as.data.frame(cor_i$correlation_s)
  df_w=as.data.frame(cor_i$correlation_w)
  print(i)
  new_w = rbind(new_w, df_w)
  new_s = rbind(new_s, df_s)
  i = i+1
}


