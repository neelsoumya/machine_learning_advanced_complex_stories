


This is a project that uses routinely collected clinical data to predict mortality in patients with serious mental illnesses.
We have developed a human-interpretable machine learning technique to make predictions for mortality. This repository has code and instructions for running code on clinical data which is private.

## Layout

    * code folder
        
      * has code 
    
    * manuscript folder
    
        * has manuscript and supplementary material
    
    * installation
    
        
        ```R
        R --no-save < INSTALL_MANY_MODULES.R
        ```

        ```R
        python -m pip install -r requirements.txt
        ```
        
    * data extraction from clinical database
    
        * data_saver.rmd
        
        * usage:
        
            * run or knit in R Studio
        
    * classical statistical models
    
        * analysis_mortality_cpft_metafeatures.rmd
        
        * usage:
        
            * run or knit in R Studio
            
        * advanced_analysis_mortality_cpft_metafeatures.rmd
        
                * advanced code for follow-up paper

        * 95CI_logistic_regression.R

                * script to calculate 95% CI of logistic regression model

        * glmnet.R

                * script to perform LASSO GLM on CPFT model

        * pca_logistic.R

                * do glmnet, PCA + LR, PCA + RF, RF alone on data CPFT
        
    * machine learning models
    
        * mixed_datatype_autoencoder_keras_MLP_multicategorical4_column_deleted_contrast.py
        
        * usage:
        
            * python mixed_datatype_autoencoder_keras_MLP_multicategorical4_column_deleted_contrast.py df_med_suicide_features_lengthened_final_ml_f20_FULL.csv  drug_suidata_autoencoder_f20_MLP_FULL  0.5
            
            * python mixed_datatype_autoencoder_keras_MLP_multicategorical4_column_deleted_contrast_bootstrap.py df_med_suicide_features_lengthened_final_ml_f20_FULL.csv  drug_suidata_autoencoder_f20_MLP_FULL  0.5
    
            * python advanced_mixed_datatype_autoencoder_keras_MLP_multicategorical4_column_deleted_contrast.py df_med_suicide_features_lengthened_final_ml_f20_FULL_column_deleted.csv  drug_suidata_autoencoder_f20_MLP_FULL_column_deleted_advanced  0.5
                
            * advanced_analysis_mortality_cpft_metafeatures.rmd
        
                * advanced code for follow-up paper
                
            * python ImmuneModel_nomovement_AIS_COMPLETELY_GENERIC_healthcare_bioAI.py tcell_survival_ImmuneModel_nomovement_AIS_COMPLETELY_GENERIC_healthcare_bioAI  tcell_matches_ImmuneModel_nomovement_AIS_COMPLETELY_GENERIC_healthcare_bioAI 1000 False False 100 False   50 1000 2 nocheckauto turnstile notcrseq
                
                * advanced code for follow-up paper data augmentation
                
            * parse_for_negative_selection.R
            
                * advanced code for follow-up paper data augmentation    
                
            * generic_cv_test.py

                * cross validation error for generic depth and activation functions   

            * read_repeat_loss_crossval.R
            
                * cross validation error over repeated calls to same ML model

            * repeat_heatmap_ML_model.R

                * run ML model and then generate class contrastive ML heatmap for repeat replicability 
                
    * synthetic_data
    
        * code for synthetic data
        
                * autoencoder and feedforward neural network on synthetic data
        
                * calculating and testing standardised mortality ratios
        
    * public_data
    
        * has public data on mortality rates
        
    * other functions and helper scripts (also in functions folder)
    
        * convert_aupr_to_ggplot.R
            * makes AUPR plots in ggplot
        * datetimefunc.R
            * miscellaneous datetime functions from rlib 
        * feature_to_metafeature_mapping.tsv
            * maps features to metafeatures
        * short_descriptions_icd10.tsv
            * stores short text descriptions for some ICD-10 codes
        * chebi_ontology.csv
            * stores drug categories and groupings
        * cam_project.bib
            * bibliography file
        
        

## Contact
    * Soumya Banerjee
    * sb2333@cam.ac.uk
    
