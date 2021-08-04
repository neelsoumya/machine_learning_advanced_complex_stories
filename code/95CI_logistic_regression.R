# run till line 1792 of analysis_mortality_cpft_metafeatures.rmd
# script to calculate 95% CI of logistic regression model


list_auc <- NULL # list to store AUC

for (x in c(1:100))
{

    TRAIN = sample(c(TRUE,FALSE),
                   nrow(df_f20),
                   replace = TRUE)
    TEST = (!TRAIN)
    df_f20_TRAIN = df_f20[TRAIN,]
    df_f20_TEST  = df_f20[TEST,]


    idx_to_keep_temp <- which( str_metafeatures %in% colnames(df_f20) )
    # which of these metafeatures are presnet stire them in a list
    str_metafeatures_present_temp <- str_metafeatures[idx_to_keep_temp]

    str_query = " mylogit_glm_all_f20 <- glm(formula = Death_Flag ~ age_now_scaled "
    for (str_temp_colname in  unique(str_metafeatures_present_temp) )
    {

      str_query = paste0(str_query, " + ", str_temp_colname)

    }

    str_query = paste0(str_query, " , data = df_f20_TRAIN, family = 'binomial' ) ")

    # str_query 
    # str_query

    # evaluate 
    eval(parse(text=str_query)) 

    #summary(mylogit_glm_all_f20)

    # log odds plots


    ###############################
    # Calculate ROC AUPR AUC etc
    ###############################

    prob = predict(mylogit_glm_all_f20, newdata=df_f20_TEST, type=c("response"))
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
    cat('... \n')
    cat(roc$auc) 
  
    # add to list
    list_auc <- cbind(list_auc, roc$auc)

}

# TODO: get 95% CI boostrap using summer school formula

lb = quantile(list_auc, 0.025)
ub = quantile(list_auc, 0.975)

cat('95% CI .... \n')
cat(lb)
cat('\n')
cat(ub)
