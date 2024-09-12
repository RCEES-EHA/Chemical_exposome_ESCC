##  Permutational multivariate analysis of variance，PERMANOVA，nonparametric multivariate analysis of variance）、or ADONIS analysis

library(vegan)

org_exp<-read.csv("",row.names = 1)
sample_info <- read.csv("",row.names = 1)


adonis_result_organics<- adonis2(org_exp ~location+gender+age+BMI+TNM_stage+tumor_location+family_cancer_history+education+ smoking_history+drinking_habit+
                         green_tea_habit+regular_staple_food+frenquency_on_pickled_food+frenquency_on_fresh_fruit_vegetable+preference_of_hot_food_and_water+preference_of_hard_food+daily_meal_times+ medication_history+LYM_pcent+MON_pcent+NEU_pcent+GLB+ALT_.AST+ALP+γ_GGT+UREA+CRE, data = sample_info, permutations = 999, method="bray",by="margin")
adonis_result_organics
##---------------------------------------------------------------------------------------------
inorg_exp<-read.csv("",row.names = 1)
adonis_result_inorganics<- adonis2(inorg_exp ~location+gender+age+BMI+TNM_stage+tumor_location+family_cancer_history+education+ smoking_history+drinking_habit+
                                   green_tea_habit+regular_staple_food+frenquency_on_pickled_food+frenquency_on_fresh_fruit_vegetable+preference_of_hot_food_and_water+preference_of_hard_food+daily_meal_times+ medication_history+LYM_pcent+MON_pcent+NEU_pcent+GLB+ALT_.AST+ALP+γ_GGT+UREA+CRE , data = sample_info, permutations = 999, method="bray",by="margin")
adonis_result_inorganics
##--------------------------------------------------------------------

endo_exp<-read.csv("",row.names = 1)
adonis_result_endo<- adonis2(endo_exp ~location+gender+age+BMI+TNM_stage+tumor_location+family_cancer_history+education+ smoking_history+drinking_habit+
                                     green_tea_habit+regular_staple_food+frenquency_on_pickled_food+frenquency_on_fresh_fruit_vegetable+preference_of_hot_food_and_water+preference_of_hard_food+daily_meal_times+ medication_history+LYM_pcent+MON_pcent+NEU_pcent+GLB+ALT_.AST+ALP+γ_GGT+UREA+CRE , data = sample_info, permutations = 999, method="bray",by="margin")
adonis_result_endo
##-------------------------------

exo_exp<-read.csv("",row.names = 1)
adonis_result_exo<- adonis2(exo_exp ~location+gender+age+BMI+TNM_stage+tumor_location+family_cancer_history+education+ smoking_history+drinking_habit+
                               green_tea_habit+regular_staple_food+frenquency_on_pickled_food+frenquency_on_fresh_fruit_vegetable+preference_of_hot_food_and_water+preference_of_hard_food+daily_meal_times+ medication_history+LYM_pcent+MON_pcent+NEU_pcent+GLB+ALT_.AST+ALP+γ_GGT+UREA+CRE , data = sample_info, permutations = 999, method="bray",by="margin")
adonis_result_exo
##-----------

EMs_exp<-read.csv("",row.names = 1)
adonis_result_EMs<- adonis2(EMs_exp ~location+gender+age+BMI+TNM_stage+tumor_location+family_cancer_history+education+ smoking_history+drinking_habit+
                              green_tea_habit+regular_staple_food+frenquency_on_pickled_food+frenquency_on_fresh_fruit_vegetable+preference_of_hot_food_and_water+preference_of_hard_food+daily_meal_times+ medication_history+LYM_pcent+MON_pcent+NEU_pcent+GLB+ALT_.AST+ALP+γ_GGT+UREA+CRE , data = sample_info, permutations = 999, method="bray",by="margin")
adonis_result_EMs

##-------------------

NEMs_exp<-read.csv("",row.names = 1)
adonis_result_NEMs<- adonis2(NEMs_exp ~location+gender+age+BMI+TNM_stage+tumor_location+family_cancer_history+education+ smoking_history+drinking_habit+
                              green_tea_habit+regular_staple_food+frenquency_on_pickled_food+frenquency_on_fresh_fruit_vegetable+preference_of_hot_food_and_water+preference_of_hard_food+daily_meal_times+ medication_history+LYM_pcent+MON_pcent+NEU_pcent+GLB+ALT_.AST+ALP+γ_GGT+UREA+CRE , data = sample_info, permutations = 999, method="bray",by="margin")
adonis_result_NEMs
##-------------------------------------

write.csv(adonis_result_organics, file = 'adonis_result_organics.csv', row.names = T, quote = FALSE, na = '')
write.csv(adonis_result_inorganics, file = 'adonis_result_inorganics.csv', row.names = T, quote = FALSE, na = '')
write.csv(adonis_result_endo, file = 'adonis_result_endo.csv', row.names = T, quote = FALSE, na = '')
write.csv(adonis_result_exo, file = 'adonis_result_exo.csv', row.names = T, quote = FALSE, na = '')
write.csv(adonis_result_EMs, file = 'adonis_result_EMs.csv', row.names = T, quote = FALSE, na = '')
write.csv(adonis_result_NEMs, file = 'adonis_result_NEMs.csv', row.names = T, quote = FALSE, na = '')
##-----------------------------------------------------------------
## dissimilarity
endo_bray_dis <-vegdist(endo_exp, method="bray")
endo_bray_dis=as.matrix(endo_bray_dis)
write.csv(endo_bray_dis, file = 'endo_bray_dis.csv', row.names = T, quote = FALSE, na = '')


exo_bray_dis <-vegdist(exo_exp, method="bray")
exo_bray_dis=as.matrix(exo_bray_dis)
write.csv(exo_bray_dis, file = 'exo_bray_dis.csv', row.names = T, quote = FALSE, na = '')

EMs_bray_dis <-vegdist(EMs_exp, method="bray")
EMs_bray_dis=as.matrix(EMs_bray_dis)
write.csv(EMs_bray_dis, file = 'EMs_bray_dis.csv', row.names = T, quote = FALSE, na = '')

NEMs_bray_dis <-vegdist(NEMs_exp, method="bray")
NEMs_bray_dis=as.matrix(NEMs_bray_dis)
write.csv(NEMs_bray_dis, file = 'NEMs_bray_dis.csv', row.names = T, quote = FALSE, na = '')
