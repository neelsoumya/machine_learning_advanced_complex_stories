#############################################################
# do GLM LASSO on data CPFT
# run after running everything till line 
#     1792 on analysis_mortality_cpft_metafeatures.rmd
#############################################################

# https://glmnet.stanford.edu/articles/glmnet.html

library(glmnet)

# create data frame of X
df_x = cbind(df_f20_TRAIN$age_now_scaled, df_f20_TRAIN$dementia_alzheimer, df_f20_TRAIN$delirium, df_f20_TRAIN$mild_cognitive_disorder,
             df_f20_TRAIN$abuse_alcohol_drugs, df_f20_TRAIN$schizophrenia, df_f20_TRAIN$bipolar_disorder, df_f20_TRAIN$major_depression,
             df_f20_TRAIN$emotional_personality_disorder, df_f20_TRAIN$respiratory, df_f20_TRAIN$cardiovascular, df_f20_TRAIN$diabetes,
             df_f20_TRAIN$self_harm, df_f20_TRAIN$lack_family_support, df_f20_TRAIN$personal_risk_factors, df_f20_TRAIN$SGA,
             df_f20_TRAIN$anti_depressant, df_f20_TRAIN$suicide_any, df_f20_TRAIN$dementia_drug, df_f20_TRAIN$antimanic_drug,
             df_f20_TRAIN$thyroid, df_f20_TRAIN$FGA, df_f20_TRAIN$diuretic, df_f20_TRAIN$anti_hypertensive, df_f20_TRAIN$aspirin)


# call glmnet
mylogit_glm_all_f20_glmnet <- glmnet::glmnet(
                                             y = df_f20_TRAIN$Death_Flag, 
  x = df_x #age_now_scaled  + dementia_alzheimer + delirium + mild_cognitive_disorder + abuse_alcohol_drugs + schizophrenia + bipolar_disorder + major_depression + emotional_personality_disorder + respiratory + cardiovascular + diabetes + self_harm + lack_family_support + personal_risk_factors + SGA + anti_depressant + suicide_any + dementia_drug + antimanic_drug + thyroid + FGA + diuretic + anti_hypertensive + aspirin , data = df_f20_TRAIN, family = 'binomial'
                                             , family = 'binomial')
  
  #data = df_f20_TRAIN)

# cross validation
cvfit <- glmnet::cv.glmnet(
                                             y = df_f20_TRAIN$Death_Flag, 
  x = df_x #age_now_scaled  + dementia_alzheimer + delirium + mild_cognitive_disorder + abuse_alcohol_drugs + schizophrenia + bipolar_disorder + major_depression + emotional_personality_disorder + respiratory + cardiovascular + diabetes + self_harm + lack_family_support + personal_risk_factors + SGA + anti_depressant + suicide_any + dementia_drug + antimanic_drug + thyroid + FGA + diuretic + anti_hypertensive + aspirin , data = df_f20_TRAIN, family = 'binomial'
                                             , family = 'binomial')

# metrics on test set
df_x_test = cbind(df_f20_TEST$age_now_scaled, df_f20_TEST$dementia_alzheimer, df_f20_TEST$delirium, df_f20_TEST$mild_cognitive_disorder,
             df_f20_TEST$abuse_alcohol_drugs, df_f20_TEST$schizophrenia, df_f20_TEST$bipolar_disorder, df_f20_TEST$major_depression,
             df_f20_TEST$emotional_personality_disorder, df_f20_TEST$respiratory, df_f20_TEST$cardiovascular, df_f20_TEST$diabetes,
             df_f20_TEST$self_harm, df_f20_TEST$lack_family_support, df_f20_TEST$personal_risk_factors, df_f20_TEST$SGA,
             df_f20_TEST$anti_depressant, df_f20_TEST$suicide_any, df_f20_TEST$dementia_drug, df_f20_TEST$antimanic_drug,
             df_f20_TEST$thyroid, df_f20_TEST$FGA, df_f20_TEST$diuretic, df_f20_TEST$anti_hypertensive, df_f20_TEST$aspirin)

metrics_test <- glmnet::assess.glmnet(mylogit_glm_all_f20_glmnet, 
                      newy = df_f20_TEST$Death_Flag,
                      newx = df_x_test)

metrics_test$auc

# TODO: repeat few times and get 95% CI of AUC
list_auc_glmnet <- c(metrics_test$auc)
lb = quantile(list_auc_glmnet, 0.025)
ub = quantile(list_auc_glmnet, 0.975)
cat(lb)
cat(ub)
# TODO: PCA in R code and feed to ML
# TODO: PCA in R code and feed to LR model


###################
# Perform PCA
###################
pr.out = prcomp(df_x, scale = TRUE)
names(pr.out)


###################
# Biplot
###################
biplot(x = pr.out, scale = 0)

pr.out$rotation = -pr.out$rotation
pr.out$x  = -pr.out$x
biplot(x = pr.out, scale = 0)

ggbiplot(pr.out, labels =  rownames(USArrests))
fviz_pca_biplot(pr.out)


######################################
# Plot variance explained
######################################
pr.var = pr.out$sdev^2
proportion_variance_explained = pr.var/sum(pr.var)
plot(proportion_variance_explained, xlab="Principal Component", ylab="Proportion of Variance Explained", ylim=c(0,1),type='b')
plot(cumsum(proportion_variance_explained), xlab="Principal Component", ylab="Cumulative Proportion of Variance Explained", ylim=c(0,1),type='b')



