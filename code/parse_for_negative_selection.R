# Program tio read output of Rmd and produceinput for negative selection proghram
#   output from analysis_mortality_cpft_medicine2.Rmd
# also parse for autoencoder


library(sqldf)


df_t_records_patients_using_topMeds_pca <- read.csv('t_records_patients_using_topMeds_pca.tsv', 
                                                    sep = '\t', header = TRUE, 
                                                    stringsAsFactors=FALSE, na.strings="..") # ,strip.white = TRUE)


#  all patients and binned age, categroized medicine
df_t_records_patients_using_topMeds_pca_ALL_binnedage_catmed <- 
  sqldf::sqldf(" SELECT *, 
               CASE 
               WHEN age_now < 10 THEN '0000'
               WHEN age_now >= 10 and age_now < 20 THEN '0001'
               WHEN age_now >= 20 and age_now < 30 THEN '0010'
               WHEN age_now >= 35 and age_now < 40 THEN '0011'
               WHEN age_now >= 40 and age_now < 45 THEN '0100'
               WHEN age_now >= 45 and age_now < 50 THEN '0101'
               WHEN age_now >= 50 and age_now < 55 THEN '0110'
               WHEN age_now >= 55 and age_now < 60 THEN '0111'
               WHEN age_now >= 60 and age_now < 65 THEN '1000'
               WHEN age_now >= 65 and age_now < 70 THEN '1001'
               WHEN age_now >= 70 and age_now < 75 THEN '1010'
               WHEN age_now >= 75 and age_now < 80 THEN '1011'
               WHEN age_now >= 80 and age_now < 85 THEN '1100'
               WHEN age_now >= 85 and age_now < 90 THEN '1101'
               WHEN age_now >= 90 and age_now < 95 THEN '1110'
               WHEN age_now >= 95                  THEN '1111'
               ELSE '0000' 
               END
               AS age_range
               FROM df_t_records_patients_using_topMeds_pca
               ")



idx_remove_2_temp <- which( colnames(df_t_records_patients_using_topMeds_pca_ALL_binnedage_catmed) == 'age_now' )
# idx_remove_3 <- which( colnames(df_t_records_patients_using_topMeds_pca_alive) == 'age_range' )

## RMPVE TGHSE COLUMNS
df_t_records_patients_using_topMeds_pca_ALL_binnedage_catmed_truncated <- 
  df_t_records_patients_using_topMeds_pca_ALL_binnedage_catmed[,-c(idx_remove_2_temp)]

# death_flag column needs to be at end; so add new column, remove old column, then rename column
df_t_records_patients_using_topMeds_pca_ALL_binnedage_catmed_truncated$Death_Flag2 <- df_t_records_patients_using_topMeds_pca_ALL_binnedage_catmed_truncated$Death_Flag
df_t_records_patients_using_topMeds_pca_ALL_binnedage_catmed_truncated$Death_Flag <- NULL

idx_temp <- which( colnames(df_t_records_patients_using_topMeds_pca_ALL_binnedage_catmed_truncated) == 'Death_Flag2' )
colnames(df_t_records_patients_using_topMeds_pca_ALL_binnedage_catmed_truncated)[idx_temp] <- "Death_Flag"

# save thi sno separators no clumn name for use in stochastic simulator
write.table(df_t_records_patients_using_topMeds_pca_ALL_binnedage_catmed_truncated,
            file='df_t_records_patients_using_topMeds_pca_ALL_binnedage_catmed_truncated.tsv',
            row.names = FALSE, quote=FALSE, append = FALSE, sep = "\t",
            col.names = FALSE)


#  first select only those that are alive 
#   NOTE: this could be control in the future
df_t_records_patients_using_topMeds_pca_alive <- sqldf::sqldf("select * 
                                                             from df_t_records_patients_using_topMeds_pca
                                                             where Death_Flag = 0  
                                                             ")

#df_t_records_patients_using_topMeds_pca_alive <- 
#  sqldf::sqldf(" SELECT *, 
#               CASE 
#               WHEN age_now < 10 THEN '0,0,0,0'
#               WHEN age_now >= 10 and age_now < 20 THEN '0,0,0,1'
#               WHEN age_now >= 20 and age_now < 30 THEN '0,0,1,0'
#               ELSE '0,0,0,0' 
#               END
#               AS age_range
#               FROM df_t_records_patients_using_topMeds_pca_alive
#               ")




df_t_records_patients_using_topMeds_pca_alive <- 
        sqldf::sqldf(" SELECT *, 
                          CASE 
                            WHEN age_now < 10 THEN '0000'
                            WHEN age_now >= 10 and age_now < 20 THEN '0001'
                            WHEN age_now >= 20 and age_now < 30 THEN '0010'
                            WHEN age_now >= 35 and age_now < 40 THEN '0011'
                            WHEN age_now >= 40 and age_now < 45 THEN '0100'
                            WHEN age_now >= 45 and age_now < 50 THEN '0101'
                            WHEN age_now >= 50 and age_now < 55 THEN '0110'
                            WHEN age_now >= 55 and age_now < 60 THEN '0111'
                            WHEN age_now >= 60 and age_now < 65 THEN '1000'
                            WHEN age_now >= 65 and age_now < 70 THEN '1001'
                            WHEN age_now >= 70 and age_now < 75 THEN '1010'
                            WHEN age_now >= 75 and age_now < 80 THEN '1011'
                            WHEN age_now >= 80 and age_now < 85 THEN '1100'
                            WHEN age_now >= 85 and age_now < 90 THEN '1101'
                            WHEN age_now >= 90 and age_now < 95 THEN '1110'
                            WHEN age_now >= 95                  THEN '1111'
                          ELSE '0000' 
                          END
                    AS age_range
                    FROM df_t_records_patients_using_topMeds_pca_alive
        ")


# save table fo nergatiove selection prhgram
write.table(df_t_records_patients_using_topMeds_pca_alive,
            file='df_t_records_patients_using_topMeds_pca_alive.csv',
            row.names = FALSE, quote=FALSE, append = FALSE, sep = ",")

# save in another format no heaers amd take out 
# colnames(df_t_records_patients_using_topMeds_pca_alive)[74] <- "age_range_final"

# npw remove unncessary columns lke age_range, Death_Flag, age_now
idx_remove_1 <- which( colnames(df_t_records_patients_using_topMeds_pca_alive) == 'Death_Flag' )
idx_remove_2 <- which( colnames(df_t_records_patients_using_topMeds_pca_alive) == 'age_now' )
# idx_remove_3 <- which( colnames(df_t_records_patients_using_topMeds_pca_alive) == 'age_range' )

## RMPVE TGHSE COLUMNS
df_t_records_patients_using_topMeds_pca_alive_truncated <- 
df_t_records_patients_using_topMeds_pca_alive[,-c(idx_remove_1,idx_remove_2)]

# save thi sno separators no clumn name for use in stochastic simulator
write.table(df_t_records_patients_using_topMeds_pca_alive_truncated,
            file='df_t_records_patients_using_topMeds_pca_alive_truncated.CSV',
            row.names = FALSE, quote=FALSE, append = FALSE, sep = "",
            col.names = FALSE)



##############################
# parse for autoencoder
##############################

df_t_records_patients_using_topMeds_pca <- 
  sqldf::sqldf(" SELECT *, 
               CASE 
               WHEN age_now < 10 THEN '0,0,0,0'
               WHEN age_now >= 10 and age_now < 20 THEN '0,0,0,1'
               WHEN age_now >= 20 and age_now < 30 THEN '0,0,1,0'
               WHEN age_now >= 35 and age_now < 40 THEN '0,0,1,1'
               WHEN age_now >= 40 and age_now < 45 THEN '0,1,0,0'
               WHEN age_now >= 45 and age_now < 50 THEN '0,1,0,1'
               WHEN age_now >= 50 and age_now < 55 THEN '0,1,1,0'
               WHEN age_now >= 55 and age_now < 60 THEN '0,1,1,1'
               WHEN age_now >= 60 and age_now < 65 THEN '1,0,0,0'
               WHEN age_now >= 65 and age_now < 70 THEN '1,0,0,1'
               WHEN age_now >= 70 and age_now < 75 THEN '1,0,1,0'
               WHEN age_now >= 75 and age_now < 80 THEN '1,0,1,1'
               WHEN age_now >= 80 and age_now < 85 THEN '1,1,0,0'
               WHEN age_now >= 85 and age_now < 90 THEN '1,1,0,1'
               WHEN age_now >= 90 and age_now < 95 THEN '1,1,1,0'
               WHEN age_now >= 95                  THEN '1,1,1,1'
               ELSE '0,0,0,0' 
               END
               AS age_range
               FROM df_t_records_patients_using_topMeds_pca
               ")


# save table fo nergatiove selection prhgram
#write.table(df_t_records_patients_using_topMeds_pca,
#            file='df_t_records_patients_using_topMeds_pca_alive.csv',
#            row.names = FALSE, quote=FALSE, append = FALSE, sep = ",")

# save in another format no heaers amd take out 
# colnames(df_t_records_patients_using_topMeds_pca_alive)[74] <- "age_range_final"

# npw remove unncessary columns lke age_range, Death_Flag, age_now
idx_remove_1 <- which( colnames(df_t_records_patients_using_topMeds_pca) == 'Death_Flag' )
idx_remove_2 <- which( colnames(df_t_records_patients_using_topMeds_pca) == 'age_now' )
# idx_remove_3 <- which( colnames(df_t_records_patients_using_topMeds_pca_alive) == 'age_range' )

## REMOVE THESE COLUMNS
df_t_records_patients_using_topMeds_pca_truncated <- 
df_t_records_patients_using_topMeds_pca[,-c(idx_remove_1,idx_remove_2)]

# save thi sno separators no clumn name for use in stochastic simulator
write.table(df_t_records_patients_using_topMeds_pca_truncated,
            file='df_t_records_patients_using_topMeds_pca_truncated.csv',
            row.names = FALSE, quote=FALSE, append = FALSE, sep = ",",
            col.names = FALSE)
