

University of Cambridge Mental Health Prediction Project.

This is a project that uses routinely collected clinical data to predict mortality in patients with serious mental illnesses.
We have developed a human-interpretable machine learning technique to make predictions for mortality. 

## Files and Layout

    * code folder
        
      * has code 
    
    * manuscript folder
    
        * has manuscript and supplementary material
    
    * installation
    
        * R --no-save < INSTALL_MANY_MODULES.R
        * python -m pip install -r requirements.txt
        
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

        
##  Installations

Install Anaconda (Python distribution)

and install all packages by running the following command at the Terminal

```R
pip install -r requirements.txt
```

For the R packages, install R and run the following command 

```R
R --no-save < INSTALL_MANY_MODULES.R
```

## Citation and manuscript

     Forthcoming        

## Contact
    * Soumya Banerjee
    * sb2333@cam.ac.uk
  
