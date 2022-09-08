####################
#Project script    # 
#multi_omics       #
#PC Lenovo L13     #
#R version 4.1.1   #
#Windows x86_64-w64#
#Script by MEmam   #
####################

#Integromics_function_return_directory_contains_all_tools_features_loading,samples_scores, correlation, and visualization_for_all_methods)
Integromics<- function(n_features, n_samples, n_groups,
                       n_factors){
  #' @param nbFeatures vector with number of features per dataset
  #' @param nbSamples integer, number of samples
  #' @param n_groups integer, number of groups
  #' @param nbFactors integer, number of factors
  
  ###Create_run_directory
  dir.create("Automation_MOFA_simulator")
  #Change_dir_to_run_directory
  setwd("Automation_MOFA_simulator")
  #Path_variable_to_store_main_run_folder
  path= getwd()
  #Simulate_data
  MOFAexample <- make_example_data(n_views=2, n_features=n_features, n_samples = n_samples, n_groups = 1,
                                   n_factors = n_factors, likelihood = "gaussian",
                                   lscales = 1, sample_cov = NULL, as.data.frame = FALSE)
  #Note_this_part_is_manually_changed_according_to_num_of_omics
  #Split_omic1
  omic1_so= as.data.frame(MOFAexample$data$view_1)
  #Split_omic2
  omic2_so= as.data.frame(MOFAexample$data$view_2)
  #Convert_simulated_data_object_into_data_frame
  Simulated_data_model = as.data.frame(rbind(omic1_so, omic2_so))
  #Export_merged_data_frame_into_run_directory
  write.csv(Simulated_data_model,"simulated_data_frame.csv" )
  #Ground_truth_sample_latent_variable(score)_stored_on_data_z_vector
  sample_score = as.data.frame(MOFAexample$Z)
  #Ground_truth_Features_loading(weight)_W_vector_is_a_weight_vector_per_each_omics
  feature_weight= as.data.frame(MOFAexample$W)
  #Create_feature_weight_variable_merge_all_omics_weight
  omic1_weight = feature_weight[1]
  colnames(omic1_weight) = "weight"
  omic2_weight = feature_weight[2]
  colnames(omic2_weight) = "weight"
  feature_weight = rbind(omic1_weight, omic2_weight)
  #Create_MOFA_model
  model = MOFA_function(MOFAexample$data, num = 1)
  #Factor1_sample_score_data frame
  factor1_df_score= MOFA_score(model = model, Factor_num = 1)
  #Factor1_feature_loading_data frame
  factor1_df_loading= MOFA_loading(model = model, Factor_num = 1)
  #Prepare_data_for_FABIA_model
  Simulated_data_model = as.data.frame(rbind(omic1_so, omic2_so))
  #Create_FABIA_object
  Fabia_object = Fabia_function(merged_fabia_data = Simulated_data_model, num = 1)
  #Create_FABIA_sample_score_object
  Fabia_score_BC1 = Fabia_score(Fabia_object = Fabia_object, BC_num = 1)
  #Create_FABIA_feature_loading_data_frame_object
  BC1_df_loading = Fabia_loading(Fabia_object = Fabia_object, BC_num = 1)
  #Create_GFA_object
  GFA_object= GFA_function(merged_GFA_data = list(t(omic1_so), t(omic2_so)), num = 1)
  #Factor1_GFA_sample_score_data frame
  g_factor1_df_score= GFA_score(GFA_object= GFA_object, BC_num= 1)
  #Factor1_GFA_feature_loading_data frame
  g_factor1_df_loading = GFA_loading(GFA_object=GFA_object , BC_num= 1)
  #Create_MFA_object
  MFA_object = MFA_function(merged_MFA_data= cbind(t(omic1_so),t(omic2_so)), group = c(ncol(t((omic1_so))),ncol(t(omic2_so))), num = 2)
  #Factor1_MFA_sample_score_data frame
  scores1_scaled_after= MFA_score(MFA_object= MFA_object, num = 1)
  #Factor1_MFA_feature_loading_data frame
  loadings1_scaled_after= MFA_loading(MFA_object= MFA_object, num = 1)
  #Combine_feature_loading_from_all_tools_in_one_data_frame
  correlation_w= as.data.frame(cbind(feature_weight,factor1_df_loading,BC1_df_loading, loadings1_scaled_after, g_factor1_df_loading))
  #Change_col_name_according_to_each_tool
  colnames(correlation_w)= c("True_w", "MOFA_w","FABIA_w",  "MFA_w", "GFA_w")
  #Export_feature_loading_data_frame (correlation_w)
  write.csv(correlation_w,"Feature_loadings_for_all_tools.csv" )
  #Perform_feature_loading_pearson_correlation_against_ground_true_loadings
  correlation_w_ground_true= cor(correlation_w$"True_w", correlation_w)  
  #Combine_sample_scores_from_all_tools_in_one_data_frame
  correlation_s= as.data.frame(cbind(sample_score, Fabia_score_BC1, factor1_df_score, scores1_scaled_after, g_factor1_df_score))
  #Change_col_name_according_to_each_tool
  colnames(correlation_s)= c("True_s", "FABIA_s", "MOFA_s", "MFA_s", "GFA_s")
  #Export_sample_score_data_frame (correlation_s)
  write.csv(correlation_s,"Sample_scores_for_all_tools.csv" )
  #erform_sample_score_pearson_correlation_against_ground_true_sample_score
  correlation_s_ground_true= cor(correlation_s$"True_s", correlation_s)  
  ###############
  ##ploting
  ##
  ###############
  #Create_MOFA_directory_to_save_MOFA_plots
  dir.create("MOFA+")
  #set_directory_to_MOFA+
  setwd("MOFA+")
  #Add_new_column_called_label_to_MOFA_loading_data_frame_to_color_each_omic_with_different_color
  factor1_df_loading$label = 1
  #Omic1_index_take_value_omic1_in_label_column
  factor1_df_loading$label[1: n_features] = "omic1"
  #Omic2_index_take_value_omic2_in_label_column
  factor1_df_loading$label[(n_features+1): (n_features * n_views)] = "omic2"
  #Plot_MOFA_feature_loading_different_color_per_omics x= c(1:(n_features * n_views)
  MOFA_feature_loading_p= ggplot(factor1_df_loading, aes(x= c(1:(n_features * n_views)), y=`weights1$weights`, color=label)) + scale_colour_manual(values = c("lightcoral", "cyan3")) + 
    labs(y= "Feature weights", x = "Feature index") + geom_point(alpha=0.5, size=1) + theme_bw() +  
    theme(plot.title = element_text(hjust = 0.5, face = 'bold'), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = 'none', legend.background = element_rect(fill="white", size=0.6, linetype="solid", colour ="lemonchiffon4"))
  #Export_MOFA_feature_loading_plot
  ggsave(filename="MOFA_feature_loading_p.jpg", plot=MOFA_feature_loading_p)
  #Plot_MOFA_denstiy_loading
  MOFA_denstiy_loading_p=ggplot(correlation_w, aes(x = `MOFA_w`)) + 
    geom_histogram(aes(y = ..density..),
                   colour = 4, fill = "white", bins = 30) +
    geom_density()
  #Export_MOFA_denstiy_loading_plot
  ggsave(filename="MOFA_denstiy_loading_plot.jpg", plot=MOFA_denstiy_loading_p)
  
  #Plot_MOFA_denstiy_samples_scores
  MOFA_denstiy_scores_p=ggplot(correlation_s, aes(x = `MOFA_s`)) + 
    geom_histogram(aes(y = ..density..),
                   colour = 4, fill = "white", bins = 30) +
    geom_density()
  #Export_MOFA_denstiy_sample_plot
  ggsave(filename="MOFA_denstiy_scores_p.jpg", plot=MOFA_denstiy_scores_p)
  #visualize_MOFA_model_variance
  Mofa_model_vairance= plot_variance_explained(model, x="group", y="factor")
  #export_mofa_model_variance_image
  ggsave(filename="Mofa_model_vairance.jpg", plot=Mofa_model_vairance)
  ##########################
  #Return_to_main_run_directory
  setwd(path)
  
  #Ground_truth_weight
  
  #Create_ground_truth_directory_to_save_ground_truth_plots
  dir.create("true")
  setwd("true")
  #Add_new_column_called_label_to_ground_truth_loading_data_frame_to_color_each_omic_with_different_color
  feature_weight$label = 1
  #Omic1_index_take_value_omic1_in_label_column
  feature_weight$label[1: n_features] = "omic1"
  #Omic2_index_take_value_omic2_in_label_column
  feature_weight$label[(n_features+1): (n_features * n_views)] = "omic2"
  #Plot_ground_truth_feature_loading_different_color_per_omics
  ground_truth_feature_loading= ggplot(feature_weight, aes(x= c(1:(n_features * n_views)), y=weight, color=label)) + scale_colour_manual(values = c("lightcoral", "cyan3")) + 
    labs(y= "Feature weights", x = "Feature index") + geom_point(alpha=0.5, size=1) + theme_bw() +  
    theme(plot.title = element_text(hjust = 0.5, face = 'bold'), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = 'none', legend.background = element_rect(fill="white", size=0.6, linetype="solid", colour ="lemonchiffon4"))
  #Export_ground_truth_feature_loading_plot
  ggsave(filename="ground_truth_feature_loading.jpg", plot=ground_truth_feature_loading)
  #plot_ground_truth_sample_score
  sample_score $index_s = rownames(sample_score)
  colnames(sample_score )= c("sample_value", "sample_index")
  png(width=600, height=350)
  plot(x=sample_score$sample_index, y= sample_score $sample_value, xlab="Sample_index", ylab="Latent_score")
  dev.off()
  ###########################
  #return_to_main_run_directory
  setwd(path)
  #FABIA
  #Create_FABIA_directory_to_save_FABIA_plots
  
  dir.create("FABIA")
  setwd("FABIA")
  #Add_new_column_called_label_to_FABIA_loading_data_frame_to_color_each_omic_with_different_color
  BC1_df_loading$label = 1
  #Omic1_index_take_value_omic1_in_label_column
  BC1_df_loading$label[1: n_features] = "omic1"
  #Omic2_index_take_value_omic2_in_label_column
  BC1_df_loading$label[(n_features+1): (n_features * n_views)] = "omic2"
  #Plot_FABIA_feature_loading_different_color_per_omics
  FABIA_feature_loading_plot=ggplot(BC1_df_loading, aes(x= c(1:(n_features * n_views)), y=loadings10_FABIA, color=label)) + scale_colour_manual(values = c("lightcoral", "cyan3")) + 
    labs(y= "Feature weights", x = "Feature index") + geom_point(alpha=0.5, size=1) + theme_bw() +  
    theme(plot.title = element_text(hjust = 0.5, face = 'bold'), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = 'none', legend.background = element_rect(fill="white", size=0.6, linetype="solid", colour ="lemonchiffon4"))
  
  FABIA_denstiy_loading_plot=ggplot(correlation_w, aes(x = `FABIA_w`)) + 
    geom_histogram(aes(y = ..density..),
                   colour = 4, fill = "white", bins = 30) +
    geom_density()
  #Export_FABIA_feature_loading_plot
  ggsave(filename="FABIA_feature_loading_p.jpg", plot=FABIA_feature_loading_plot)
  #Export_FABIA_denstiy_loading_plot
  ggsave(filename="FABIA_denstiy_loading_plot.jpg", plot=FABIA_denstiy_loading_plot)
  
  #Plot_FABIA_denstiy_samples_scores
  FABIA_denstiy_scores_p=ggplot(correlation_s, aes(x = `FABIA_s`)) + 
    geom_histogram(aes(y = ..density..),
                   colour = 4, fill = "white", bins = 30) +
    geom_density()
  #Export_FABIA_denstiy_sample_plot
  ggsave(filename="FABIA_denstiy_scores_p.jpg", plot=FABIA_denstiy_scores_p)
  
  
  #return_to_main_run_directory
  
  setwd(path)
  #MFA
  #Create_MFA_directory_to_save_MFA_plots
  dir.create("MFA")
  setwd("MFA")
  #Add_new_column_called_label_to_MOFA_loading_data_frame_to_color_each_omic_with_different_color
  loadings1_scaled_after$label = 1
  #Omic1_index_take_value_omic1_in_label_column
  loadings1_scaled_after$label[1: n_features] = "omic1"
  #Omic2_index_take_value_omic2_in_label_column
  loadings1_scaled_after$label[(n_features+1): (n_features * n_views)] = "omic2"
  #Plot_MFA_feature_loading_different_color_per_omics
  MFA_loading_plot= ggplot(loadings1_scaled_after, aes(x= c(1:(n_features * n_views)), y=MFA_loading, color=label)) + scale_colour_manual(values = c("lightcoral", "cyan3")) + 
    labs(y= "Feature weights", x = "Feature index") + geom_point(alpha=0.5, size=1) + theme_bw() +  
    theme(plot.title = element_text(hjust = 0.5, face = 'bold'), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = 'none', legend.background = element_rect(fill="white", size=0.6, linetype="solid", colour ="lemonchiffon4"))
  
  MFA_denstiy_loading_plot=ggplot(correlation_w, aes(x = `MFA_w`)) + 
    geom_histogram(aes(y = ..density..),
                   colour = 4, fill = "white", bins = 30) +
    geom_density()
  #Export_MFA_feature_loading_plot
  ggsave(filename="MFA_loading_plot.jpg", plot=MFA_loading_plot)
  #Export_MFA_denstiy_loading_plot
  ggsave(filename="MFA_denstiy_loading_plot.jpg", plot=MFA_denstiy_loading_plot)
  
  #Plot_MFA_denstiy_samples_scores
  MFA_denstiy_scores_p=ggplot(correlation_s, aes(x = `MFA_s`)) + 
    geom_histogram(aes(y = ..density..),
                   colour = 4, fill = "white", bins = 30) +
    geom_density()
  #Export_FABIA_denstiy_sample_plot
  ggsave(filename="MFA_denstiy_scores_p.jpg", plot=MFA_denstiy_scores_p)
  
  #visualize_variance_across_MFA_score_per_samples
  variance_across_samples= fviz_contrib(MFA_object,choice="partial.axes", axes = 1, top =10,palette= "jco")
  ggsave(filename="MFAvariance_across_samples.jpg", plot=variance_across_samples)
  #return_to_main_run_directory
  
  setwd(path)
  #GFA
  #Create_GFA_directory_to_save_GFA_plots
  dir.create("GFA")
  setwd("GFA")
  #Add_new_column_called_label_to_GFA_loading_data_frame_to_color_each_omic_with_different_color
  g_factor1_df_loading$label = 1
  #Omic1_index_take_value_omic1_in_label_column
  g_factor1_df_loading$label[1: n_features] = "omic1"
  #Omic2_index_take_value_omic2_in_label_column
  g_factor1_df_loading$label[(n_features+1): (n_features * n_views)] = "omic2"
  #Plot_GFA_feature_loading_different_color_per_omics
  GFA_loading_plot=ggplot(g_factor1_df_loading, aes(x= c(1:(n_features * n_views)), y=GFA_weight.V1...BC_num, color=label)) + scale_colour_manual(values = c("lightcoral", "cyan3")) + 
    labs(y= "Feature weights", x = "Feature index") + geom_point(alpha=0.5, size=1) + theme_bw() +  
    theme(plot.title = element_text(hjust = 0.5, face = 'bold'), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = 'none', legend.background = element_rect(fill="white", size=0.6, linetype="solid", colour ="lemonchiffon4"))
  
  GFA_denstiy_loading_plot=ggplot(correlation_w, aes(x = `GFA_w`)) + 
    geom_histogram(aes(y = ..density..),
                   colour = 4, fill = "white", bins = 30) +
    geom_density()
  #Export_GFA_feature_loading_plot
  ggsave(filename="GFA_loading_plot.jpg", plot=GFA_loading_plot)
  #Export_GFA_denstiy_loading_plot
  ggsave(filename="GFA_denstiy_loading_plot.jpg", plot=GFA_denstiy_loading_plot)
  
  #Plot_GFA_denstiy_samples_scores
  GFA_denstiy_scores_p=ggplot(correlation_s, aes(x = `GFA_s`)) + 
    geom_histogram(aes(y = ..density..),
                   colour = 4, fill = "white", bins = 30) +
    geom_density()
  #Export_GFA_denstiy_sample_plot
  ggsave(filename="GFA_denstiy_scores_p.jpg", plot=GFA_denstiy_scores_p)
  #return_to_main_run_directory
  #Correlation_plot_for_all_tools_latent(sample)score
  setwd(path)
  library(corrplot)
  png("sample_score_correlation.png",width=600, height=350)
  corrplot(cor(correlation_s),
           method = "number",       
           order = "hclust",         # Ordering method of the matrix
           hclust.method = "ward.D", # If order = "hclust", is the cluster method to be used
           addrect = 2,              # If order = "hclust", number of cluster rectangles
           rect.col = 3,             # Color of the rectangles
           rect.lwd = 3) 
  
  dev.off()
  #Correlation_plot_for_all_tools_feature_weight
  png("Feature_loading_correlation.png",width=600, height=350)
  
  corrplot(cor(correlation_w),
           method = "number",       
           order = "hclust",         # Ordering method of the matrix
           hclust.method = "ward.D", # If order = "hclust", is the cluster method to be used
           addrect = 2,              # If order = "hclust", number of cluster rectangles
           rect.col = 3,             # Color of the rectangles
           rect.lwd = 3) 
  
  dev.off()
  correlation <- list(correlation_w_ground_true=correlation_w_ground_true,correlation_s_ground_true = correlation_s_ground_true)
  return(correlation)
}

#Integration_model_with_FABIA
library("fabia")
#Function_request_two_arguments_and_return_FABIA_object
#merged_fabia_data = data_frame_merged_two_omics
#num= number_of_factor


Fabia_function <- function(merged_fabia_data, num){
  #set_seed
  set.seed(123)
  #Create_FBAIA_biclustering_object
  resFabia = fabia(as.matrix(merged_fabia_data), p=num, alpha=0, 1, cyc=1000, spl=0.5, spz=0.5, random=1.0, center=2, norm=2, scale=0.0, lap=1.0, nL=1)
  return(resFabia)
}
#Create_feature_loading_data_frame_function
#Two_arguments_reuired_1-Fabia_object,2- BC_num = number_of_bicluster"factor"
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

#Integration_model_with_MOFA
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
  head(model_opts)
  train_opts <- get_default_training_options(MOFAobject)
  train_opts$maxiter <-25
  train_opts$convergence_mode <- "slow"
  train_opts$seed <- 42
  head(train_opts)
  MOFAobject <- prepare_mofa(
    object = MOFAobject,
    data_options = data_opts,
    model_options = model_opts,
    training_options = train_opts
  )
  
  outfile = paste0(getwd(),"/model_run.hdf5")
  MOFAobject.trained <- run_mofa(MOFAobject, outfile)
  model= MOFAobject.trained
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


#Integration_model_with_GFA
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

#Integration_model_with_Multiple_factor_Analysis
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



#Correlation_on_feature_loading_and_sample_scores
Loadings_scores_correlation = Integromics(n_views=2, n_features=1000, n_samples = 100,                         n_factors = 1)


