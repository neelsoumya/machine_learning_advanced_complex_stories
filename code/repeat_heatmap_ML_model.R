##################################################
# a repeat analysis of a heatmap from ML model
#
#  run
#   python3 generic_cv_test.py df_med_suicide_features_lengthened_final_ml_f20_FULL_column_deleted.csv  drug_suidata_autoencoder_f20_MLP_FULL_column_deleted_05_relu_sigmoid_10  0.5  relu  sigmoid
#
##################################################

####################
# load library
####################
library(sqldf)
library(pheatmap)

setwd('~/cam_project/code/cpft_mortality')

#################################################
# Class-contrastive for deep learning models
#################################################


# ```{r, echo=FALSE, fig.cap="Visualization of the amount of change predicted in the probability of death by setting a particular metafeature to 1. Rows represent patients and columns represent metafeatures. Predictions are made using logistic regression on the test set."  }
# generate a heatmap for this row is patient and column is change in predicted value
# pheatmap::pheatmap(mat = df_visualization_class_contrast[,18:31],
#                    scale = "none",
#                    show_colnames = TRUE,
#                    show_rownames = FALSE,
#                    cluster_rows = TRUE,
#                    cluster_cols = TRUE)
# load data frame generated from R
# output modiofied with following line in Librecalc or excel as header
# dementia_alzheimer,delirium,mild_cognitive_disorder,abuse_alcohol_drugs,respiratory,cardiovascular,diabetes,self_harm,lack_family_support,personal_risk_factors,SGA,anti_depressant,dementia_drug,antimanic_drug,thyroid,FGA,diuretic,blood_pressure
df_visualization_class_contrast_noextra_ML_deep_learning_model <- 
                        read.csv('df_class_contrastive_ml_from_python_relu_sigmoid_05_MODnoantimanic.csv', 
                                sep = ',', header = TRUE, 
                                stringsAsFactors=FALSE, na.strings="..") # ,strip.white = TRUE)
                                # quote="" # for dealing with ' etc in files
#df_visualization_class_contrast_noextra <- cbind(df_visualization_class_contrast[,18:31], df_visualization_class_contrast[39:41])
#df_visualization_class_contrast_noextra <- cbind(df_visualization_class_contrast_noextra, df_visualization_class_contrast[43])
#df_visualization_class_contrast_noextra <- cbind(df_visualization_class_contrast_noextra, df_visualization_class_contrast[45:49])
pheatmap::pheatmap(mat = df_visualization_class_contrast_noextra_ML_deep_learning_model,
                   scale = "none",
                   show_colnames = TRUE,
                   show_rownames = FALSE,
                   cluster_rows = TRUE,
                   cluster_cols = TRUE)
# df_visualization_class_contrast_noextra_complement <- cbind(df_visualization_class_contrast_complement[,18:31], df_visualization_class_contrast_complement[39:41])
# df_visualization_class_contrast_noextra_complement <- cbind(df_visualization_class_contrast_noextra_complement, df_visualization_class_contrast_complement[45:49])
#pheatmap::pheatmap(mat = df_visualization_class_contrast_noextra_complement,
#                   scale = "none",
#                   show_colnames = TRUE,
#                   show_rownames = FALSE,
#                   cluster_rows = TRUE,
#                   cluster_cols = TRUE)

df_visualization_class_contrast_noextra_ML_deep_learning_model_2comb <-
                    read.csv('df_class_contrastive_ml_2comb_from_python_relu_sigmoid_05.csv',
                                sep = '\t', header = FALSE,
                                stringsAsFactors=FALSE, na.strings="..") # ,strip.white = TRUE)
                                # quote="" # for dealing with ' etc in files
pheatmap::pheatmap(mat = df_visualization_class_contrast_noextra_ML_deep_learning_model_2comb,
                   scale = "none",
                   show_colnames = FALSE,
                   show_rownames = FALSE,
                   cluster_rows = TRUE,
                   cluster_cols = TRUE)
# advanced paper with 3 change simulatenously
df_visualization_class_contrast_noextra_ML_deep_learning_model_3comb <-
                    read.csv('df_class_contrastive_ml_3comb_from_python_advanced_MOD.csv',
                                sep = '\t', header = FALSE,
                                stringsAsFactors=FALSE, na.strings="..") # ,strip.white = TRUE)
                                # quote="" # for dealing with ' etc in files
pheatmap::pheatmap(mat = df_visualization_class_contrast_noextra_ML_deep_learning_model_3comb,
                   scale = "none",
                   show_colnames = FALSE,
                   show_rownames = FALSE,
                   cluster_rows = TRUE,
                   cluster_cols = TRUE)
                   
# tuples are in list_tuples_3_advanced.csv
# column names are in df_med_suicide_features_lengthened_final_ml_f20_FULL_column_deleted.csv
# 1 age_now_scaled,
# 2 dementia_alzheimer,
# 3 dementia,
# 4 mild_cognitive_disorder,
# 5 abuse_alcohol_drugs,
# 6 respiratory,
# 7 cardiovascular,
# 8 diabetes,
# 9 self_harm,
# 10 lack_family_support,
# 11 personal_risk_factors,
# 12 SGA,
# 13 SSRI,
# 14 dementia_drug,
# 15 bipolar_drug, aso called antimanic_drug
# 16 thyroid,
# 17 FGA,
# 18 diuretic,
# 19 blood_pressure,
# 20 aspirin,
# 21 Death_Flag
