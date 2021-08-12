#############################################################
# do glmnet, PCA + LR, PCA + RF, RF alone on data CPFT
# run after running everything till line 
#     1792 on analysis_mortality_cpft_metafeatures.rmd
#############################################################

# https://glmnet.stanford.edu/articles/glmnet.html

# repeat
# 10 repeats
i_num_repeats = 100

# lists to store 95% CI limits
list_super_glmnet  <- NULL
list_super_pca_lr  <- NULL
list_super_pca_rf  <- NULL
list_super_rf_solo <- NULL

for (i_temp_counter in c(1:i_num_repeats))
{  

    ############################################
    # split into test and train without seed
    ############################################
    TRAIN = sample(c(TRUE,FALSE),
                   nrow(df_f20),
                   replace = TRUE)
    TEST = (!TRAIN)
    df_f20_TRAIN = df_f20[TRAIN,]
    df_f20_TEST  = df_f20[TEST,]



    ####################
    # glmnet L1 LASSO
    ####################  
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

    # repeat few times and get 95% CI of AUC
    list_auc_glmnet <- c(metrics_test$auc)
    lb = quantile(list_auc_glmnet, 0.025)
    ub = quantile(list_auc_glmnet, 0.975)
    cat(lb)
    cat(ub)
   
    # append to list of AUC values for this model
    list_super_glmnet <- cbind( list_super_glmnet, c(metrics_test$auc) )


    ###################
    # Perform PCA
    ###################
    pr.out = prcomp(df_x)#, scale = TRUE)
    names(pr.out)


    ###################
    # Biplot
    ###################
    # biplot(x = pr.out, scale = 0)

    pr.out$rotation = -pr.out$rotation
    pr.out$x  = -pr.out$x
    # biplot(x = pr.out, scale = 0)

    # ggbiplot(pr.out, labels =  rownames(USArrests))
    # fviz_pca_biplot(pr.out)


    ######################################
    # Plot variance explained
    ######################################
    pr.var = pr.out$sdev^2
    proportion_variance_explained = pr.var/sum(pr.var)
    # plot(proportion_variance_explained, xlab="Principal Component", ylab="Proportion of Variance Explained", ylim=c(0,1),type='b')
    # plot(cumsum(proportion_variance_explained), xlab="Principal Component", ylab="Cumulative Proportion of Variance Explained", ylim=c(0,1),type='b')



    ##############################################
    # logistic regression on reduced features
    ##############################################

    # str_query = " mylogit_glm_all_f20 <- glm(formula = Death_Flag ~ age_now_scaled "
    pr.out$x # educed dimenisons
    # create a data frame
    df_pca_features_train <-
    data.frame(
     cbind( pr.out$x[,1], pr.out$x[,2], pr.out$x[,3], pr.out$x[,4], pr.out$x[,5], pr.out$x[,6], pr.out$x[,7], pr.out$x[,8], pr.out$x[,9], pr.out$x[,10] ),
     stringsAsFactors = FALSE)

    # add death flag
    df_pca_features_train$Death_Flag <- df_f20_TRAIN$Death_Flag

    # str_query = " mylogit_glm_all_f20_pca_lr <- glm(formula =  df_f20_TRAIN$Death_Flag ~  df_pca_features_train$X1 + df_pca_features_train$X2 + df_pca_features_train$X3, data = df_f20_TRAIN, family = 'binomial' ) "

    str_query = " mylogit_glm_all_f20_pca_lr <- glm(formula =  Death_Flag ~  X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, data = df_pca_features_train, family = 'binomial' ) "

    #for (str_temp_colname in  unique(str_metafeatures_present_temp) )
    #{
    #    
    #      str_query = paste0(str_query, " + ", str_temp_colname)

    #}
    #str_query = paste0(str_query, " , data = df_f20_TRAIN, family = 'binomial' ) ")
    # str_query 
    str_query
    # evaluate 

    eval(parse(text=str_query)) 

    #summary(mylogit_glm_all_f20_pca_lr)
    # log odds plots
    ###############################
    # Calculate ROC AUPR AUC etc
    ###############################
    # crteate new data TEST


    ############################
    # Perform PCA on TEST SET
    ############################
    pr_out_test = prcomp(df_x_test)#, scale = TRUE)
    names(pr_out_test)

    # create a data frame
    df_pca_features_test <-
    data.frame(
     cbind( pr_out_test$x[,1], pr_out_test$x[,2], pr_out_test$x[,3], pr_out_test$x[,4], pr_out_test$x[,5], pr_out_test$x[,6], pr_out_test$x[,7], pr_out_test$x[,8], pr_out_test$x[,9], pr_out_test$x[,10] ),
     stringsAsFactors = FALSE)
    # add Deth Flag
    df_pca_features_test$Death_Flag <- df_f20_TEST$Death_Flag

    prob = predict(mylogit_glm_all_f20_pca_lr, newdata = df_pca_features_test, type=c("response"))
    df_f20_TEST$prob = prob
    #g <- pROC::roc( paste0(str_exprn_roc_calculation), data=df_f20_TEST) #  Death_Flag ~ expression
    #plot(g)
    #cat("AUC is:", g$auc)
    #auc_atnf = g$auc
    # Precision recall curve
    #mmdata(df_t_all_gene_matched_withprobeid_osm_mod$expression, df_t_all_gene_matched_withprobeid_osm_mod$non_responder)
    #mmdata_atnf = mmdata(df_t_all_gene_matched_withprobeid_osm_mod$expression,
    #                     df_t_all_gene_matched_withprobeid_osm_mod$non_responder)
    #smcurves <- evalmod(mmdata_atnf, raw_curves = TRUE)
    # plot(smcurves, raw_curves = FALSE)
    fg <- prob[df_f20_TEST$Death_Flag == 1]
    bg <- prob[df_f20_TEST$Death_Flag == 0]
    # ROC Curve    
    roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = TRUE)
    # png("figures/roc_curve_PRROC_f20.png")
    # plot(roc)
    # dev.off()
    # PR Curve

    cat('AUC from PCA + LR: ... \n')
    cat(roc$auc)

    pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = TRUE)
    cat(pr$auc.davis.goadrich)

    # append to list of AUC values for this model
    list_super_pca_lr <- cbind( list_super_pca_lr, roc$auc )
    
    
    # png("figures/aupr_curve_PRROC_f20.png")
    # plot(pr)
    # dev.off()
    # mylogit_active <- glm( paste0(str_exprn_roc_calculation), #Death_Flag ~ expression, 
    #                      data = df_f20_TRAIN, 
    #                      family = "binomial")
    # prob_active = predict(mylogit_active, newdata=df_f20_TEST, type=c("response"))
    # df_f20_TEST$prob_active = prob_active
    #g_active <- pROC::roc( paste0(str_exprn_roc_calculation), #Death_Flag ~ prob_active,
    #                      data=df_f20_TEST)
    #plot(g_active)
    #cat("AUC is:", g_active$auc)
    #auc_active = g_active$auc
    # fg_active <- prob_active[df_f20_TEST$Death_Flag == 1]
    # bg_active <- prob_active[df_f20_TEST$Death_Flag == 0]
    # ROC Curve    
    # roc_active <- roc.curve(scores.class0 = fg_active, scores.class1 = bg_active, curve = TRUE)
    # png("figures/roc_curve_f20.png")
    # plot(roc_active)
    # dev.off()
    # PR Curve
    # pr_active <- pr.curve(scores.class0 = fg_active, scores.class1 = bg_active, curve = TRUE)



    
    
    ##################################
    # PCA + RF
    ##################################
    #    https://www.blopig.com/blog/2017/04/a-very-basic-introduction-to-random-forests-using-r/
    library(randomForest)

    # perform training
    rf_model_pca = randomForest::randomForest(formula = Death_Flag ~  X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, 
                                              data = df_pca_features_train, 
                                              ntree = 100,
                                              importance = TRUE)

    summary(rf_model_pca)

    # varImpPlot(rf_model_pca)

    # now predeict on test set
    prob_pca_rf = predict(rf_model_pca, newdata = df_pca_features_test, type=c("response"))

    # prediction_for_table <- predict(rf_classifier,validation1[,-5])
    # table(observed=validation1[,5],predicted=prediction_for_table)

    df_f20_TEST$prob_pca_rf = prob_pca_rf
    fg_pca_rf = prob_pca_rf[df_f20_TEST$Death_Flag == 1]
    bg_pca_rf = prob_pca_rf[df_f20_TEST$Death_Flag == 0]

    # roc curve
    roc_pca_rf = roc.curve(scores.class0 = fg_pca_rf, scores.class1 = bg_pca_rf, curve = TRUE)

    cat('AUC of PCA + RF ... \n')
    cat(roc_pca_rf$auc)


    # append to list of AUC values for this model
    list_super_pca_rf <- cbind( list_super_pca_rf, roc_pca_rf$auc )


    
    
    
    ##################################
    # RF alone on original features
    ##################################

    # perform training
    rf_model_solo = randomForest::randomForest(formula = Death_Flag ~ age_now_scaled  + dementia_alzheimer + delirium + mild_cognitive_disorder + abuse_alcohol_drugs + schizophrenia + bipolar_disorder + major_depression + emotional_personality_disorder + respiratory + cardiovascular + diabetes + self_harm + lack_family_support + personal_risk_factors + SGA + anti_depressant + suicide_any + dementia_drug + antimanic_drug + thyroid + FGA + diuretic + anti_hypertensive + aspirin , 
                                              data = df_f20_TRAIN, 
                                              ntree = 100,
                                              importance = TRUE)

    summary(rf_model_solo)

    # varImpPlot(rf_model_solo)

    # now predeict on test set
    prob_rf_solo = predict(rf_model_solo, newdata = df_f20_TEST, type=c("response"))

    # prediction_for_table <- predict(rf_classifier,validation1[,-5])
    # table(observed=validation1[,5],predicted=prediction_for_table)

    df_f20_TEST$prob_rf_solo = prob_rf_solo
    fg_rf_solo = prob_rf_solo[df_f20_TEST$Death_Flag == 1]
    bg_rf_solo = prob_rf_solo[df_f20_TEST$Death_Flag == 0]

    # roc curve
    roc_rf_solo = roc.curve(scores.class0 = fg_rf_solo, scores.class1 = bg_rf_solo, curve = TRUE)

    cat('AUC of RF alone on original features ... \n')
    cat(roc_rf_solo$auc)

    # append to list of AUC values for this model
    list_super_rf_solo <- cbind( list_super_rf_solo, roc_rf_solo$auc )

    # https://www.r-bloggers.com/2012/12/binary-classification-a-comparison-of-titanic-proportions-between-logistic-regression-random-forests-and-conditional-trees/

}


##############################################
# print all metrics
##############################################

cat(' ....................... \n')
cat('Metrics ....\n')
cat(' ....................... \n')

cat('GLMNET L1 LASSO 95% CI of AUC \n')
cat(  quantile(list_super_glmnet, 0.025) )
cat('\n')
cat(  quantile(list_super_glmnet, 0.975) )
cat('\n')
cat('GLMNET L1 LASSO Mean \n')
cat(  mean(list_super_glmnet) )
cat('\n')

cat('PCA + LR 95% CI of AUC \n')
cat(  quantile(list_super_pca_lr, 0.025) )
cat('\n')
cat(  quantile(list_super_pca_lr, 0.975) )
cat('\n')
cat('PCA + LR Mean \n')
cat(  mean(list_super_pca_lr) )
cat('\n')

cat('PCA + RF 95% CI of AUC \n')
cat(  quantile(list_super_pca_rf, 0.025) )
cat('\n')
cat(  quantile(list_super_pca_rf, 0.975) )
cat('\n')
cat('PCA + RF Mean \n')
cat(  mean(list_super_pca_rf) )
cat('\n')

cat('RF alone on unmodified features 95% CI of AUC \n')
cat(  quantile(list_super_rf_solo, 0.025) )
cat('\n')
cat(  quantile(list_super_rf_solo, 0.975) )
cat('\n')
cat('RF alone Mean \n')
cat(  mean(list_super_rf_solo) )
cat('\n')


cat(' ....................... \n')



