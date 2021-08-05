# do GLM LASSO on data CPFT

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


# TODO: repeat few times and get 95% CI of AUC
# TODO: PCA in R code and feed to ML
# TODO: PCA in R code and feed to LR model

# FUTURE WORK Penalizd Cox models in glmnet
# https://cran.r-project.org/web/packages/glmnet/glmnet.pdf
# family = Cox 
