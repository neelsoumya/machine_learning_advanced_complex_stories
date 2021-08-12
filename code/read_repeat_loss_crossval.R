####################################################
# read in csv files with cross validation error
#     these are output from
#  
#   python3 generic_cv_test.py df_med_suicide_features_lengthened_final_ml_f20_FULL_column_deleted.csv  drug_suidata_autoencoder_f20_MLP_FULL_column_deleted_05_relu_sigmoid_1  0.5  relu  sigmoid
#
####################################################


###############################################
# load library
###############################################
library(sqldf)

setwd('~/cam_project/code/cpft_mortality/')
###############################################
# read all cross val outputs from ML model
###############################################
df_all_rep1 <- read.csv('drug_suidata_autoencoder_f20_MLP_FULL_column_deleted_05_relu_sigmoid_1xval_loss_test.csv', 
                                        sep = '\t', header = FALSE, 
                                        stringsAsFactors=FALSE, na.strings="..") # ,strip.white = TRUE)
                                        
                                        
df_all_rep2 <- read.csv('drug_suidata_autoencoder_f20_MLP_FULL_column_deleted_05_relu_sigmoid_2xval_loss_test.csv', 
                                        sep = '\t', header = FALSE, 
                                        stringsAsFactors=FALSE, na.strings="..") # ,strip.white = TRUE)

df_all_rep3 <- read.csv('drug_suidata_autoencoder_f20_MLP_FULL_column_deleted_05_relu_sigmoid_3xval_loss_test.csv', 
                                        sep = '\t', header = FALSE, 
                                        stringsAsFactors=FALSE, na.strings="..") # ,strip.white = TRUE)


df_all_rep4 <- read.csv('drug_suidata_autoencoder_f20_MLP_FULL_column_deleted_05_relu_sigmoid_4xval_loss_test.csv', 
                                        sep = '\t', header = FALSE, 
                                        stringsAsFactors=FALSE, na.strings="..") # ,strip.white = TRUE)

df_all_rep5 <- read.csv('drug_suidata_autoencoder_f20_MLP_FULL_column_deleted_05_relu_sigmoid_5xval_loss_test.csv', 
                                        sep = '\t', header = FALSE, 
                                        stringsAsFactors=FALSE, na.strings="..") # ,strip.white = TRUE)

df_all_rep6 <- read.csv('drug_suidata_autoencoder_f20_MLP_FULL_column_deleted_05_relu_sigmoid_6xval_loss_test.csv', 
                                        sep = '\t', header = FALSE, 
                                        stringsAsFactors=FALSE, na.strings="..") # ,strip.white = TRUE)

df_all_rep7 <- read.csv('drug_suidata_autoencoder_f20_MLP_FULL_column_deleted_05_relu_sigmoid_7xval_loss_test.csv', 
                                        sep = '\t', header = FALSE, 
                                        stringsAsFactors=FALSE, na.strings="..") # ,strip.white = TRUE)

df_all_rep8 <- read.csv('drug_suidata_autoencoder_f20_MLP_FULL_column_deleted_05_relu_sigmoid_8xval_loss_test.csv', 
                                        sep = '\t', header = FALSE, 
                                        stringsAsFactors=FALSE, na.strings="..") # ,strip.white = TRUE)

                        

###############################################
# combine all these runs
###############################################
df_combined <- rbind(df_all_rep1, df_all_rep2, df_all_rep3, df_all_rep4, df_all_rep5, df_all_rep6, df_all_rep7, df_all_rep8)



###############################################
# get metrics
###############################################
cat('cross validation error ... \n')
cat('Mean .... \n')
cat(  mean(df_combined$V1 ))
cat('\n')
cat('Std ....')
cat(  sd(df_combined$V1 ))

cat('\n')
