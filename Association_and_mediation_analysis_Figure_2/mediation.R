

data <- read.csv("", check.names = T)

# covert factor for categorical variable

data[,c(2,3)]<-lapply(data[,c(2,3)], as.factor)
data[,c(6:19)]<-lapply(data[,c(6:19)], as.factor)


habit_name <- colnames(data[10:18])
chemical_name <- colnames(data[29:343])
outcome_name <- colnames(data[20:28])

##---------------------------------------------------------------------------------------------------------------------------------------------
# install.packages("bruceR")
library(bruceR)

## forward
result1 <- list()

set.seed(1)  # Setting seed outside the loop for varied bootstrap samples

for (i in habit_name) {
  result1[[i]] <- list()
  for (a in outcome_name) {
    result1[[i]][[a]] <- list()
    k <- 1
    for (j in chemical_name) {
      result1[[i]][[a]][[k]] <- PROCESS(data, y = a, x = i,
                                            meds = j,
                                            covs = c("group", "gender", "age", "BMI", "TNM_stage", "medication_history"),
                                            ci = "boot", nsim = 1000)
      k <- k + 1
    }
  }
}


## backward

result2 <- list()

set.seed(1)  # Setting seed outside the loop for varied bootstrap samples

for (i in habit_name) {
  result1[[i]] <- list()
  for (a in outcome_name) {
    result1[[i]][[a]] <- list()
    k <- 1
    for (j in chemical_name) {
      result1[[i]][[a]][[k]] <- PROCESS(data, y = j, x = i,
                                        meds = a,
                                        covs = c("group", "gender", "age", "BMI", "TNM_stage", "medication_history"),
                                        ci = "boot", nsim = 1000)
      k <- k + 1
    }
  }
}


##---------------------------------------------------------------------------------------------------------------------------------------------
# further using Mediation.package for mediation.result 


library(mediation)
library(lpSolve)

data <- read.csv("", check.names = T)
set.seed(123)




fit.totaleffect=lm(Androsterone.sulfate ~ smoking_history+group+gender+age+BMI+TNM_stage+medication_history,data)



fit.mediator=lm(ALP~smoking_history+group+gender+age+BMI+TNM_stage+medication_history,data)



fit.Y=lm(Androsterone.sulfate~smoking_history+ALP+group+gender+age+BMI+TNM_stage+medication_history,data)


results = mediate(fit.mediator, fit.Y, treat='smoking_history', mediator='ALP', robustSE = TRUE,sims = 1000)

summary(results)

medsen1<-medsens(results, rho.by = 0.01, sims = 100, effect.type = "indirect")

summary(medsen1)



fit.totaleffect=lm(N..Formylkynurenine~ frenquency_on_fresh_fruit_vegetable+group+gender+age+BMI+TNM_stage+medication_history,data)



fit.mediator=lm(NEU_pcent~frenquency_on_fresh_fruit_vegetable+group+gender+age+BMI+TNM_stage+medication_history,data)



fit.Y=lm(N..Formylkynurenine~frenquency_on_fresh_fruit_vegetable+NEU_pcent+group+gender+age+BMI+TNM_stage+medication_history,data)


results = mediate(fit.mediator, fit.Y, treat='frenquency_on_fresh_fruit_vegetable', mediator='NEU_pcent', robustSE = TRUE,sims = 1000)

summary(results)


## get rho for every significant mediation  

# for forward

#1 

fit.totaleffect=lm(ALT_.AST ~ daily_meal_times+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.mediator=lm(X1.2.3.Trihydroxybenzene~daily_meal_times+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.Y=lm(ALT_.AST~daily_meal_times+X1.2.3.Trihydroxybenzene+group+gender+age+BMI+TNM_stage+medication_history,data)

results = mediate(fit.mediator, fit.Y, treat='daily_meal_times', mediator='X1.2.3.Trihydroxybenzene', robustSE = TRUE,sims = 1000)

summary(results)

medsen1<-medsens(results, rho.by = 0.01, sims = 100, effect.type = "indirect")

summary(medsen1)

#2 

fit.totaleffect=lm(LYM_pcent ~ smoking_history+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.mediator=lm(L.Pipecolic.acid~smoking_history+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.Y=lm(LYM_pcent~ smoking_history+L.Pipecolic.acid+group+gender+age+BMI+TNM_stage+medication_history,data)

results = mediate(fit.mediator, fit.Y, treat='smoking_history', mediator='L.Pipecolic.acid', robustSE = TRUE,sims = 1000)

summary(results)

medsen2<-medsens(results, rho.by = 0.01, sims = 100, effect.type = "indirect")

summary(medsen2)

#3

fit.totaleffect=lm(GLB ~ smoking_history+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.mediator=lm(L.Pipecolic.acid~smoking_history+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.Y=lm(GLB~smoking_history+L.Pipecolic.acid+group+gender+age+BMI+TNM_stage+medication_history,data)

results = mediate(fit.mediator, fit.Y, treat='smoking_history', mediator='L.Pipecolic.acid', robustSE = TRUE,sims = 1000)

summary(results)

medsen3<-medsens(results, rho.by = 0.01, sims = 100, effect.type = "indirect")

summary(medsen3)

#4

fit.totaleffect=lm(LYM_pcent ~ drinking_habit+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.mediator=lm(L.Pipecolic.acid~drinking_habit+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.Y=lm(LYM_pcent~drinking_habit+L.Pipecolic.acid+group+gender+age+BMI+TNM_stage+medication_history,data)

results = mediate(fit.mediator, fit.Y, treat='drinking_habit', mediator='L.Pipecolic.acid', robustSE = TRUE,sims = 1000)

summary(results)

medsen4<-medsens(results, rho.by = 0.01, sims = 100, effect.type = "indirect")

summary(medsen4)

#5

fit.totaleffect=lm(NEU_pcent ~ drinking_habit+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.mediator=lm(L.Pipecolic.acid~drinking_habit+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.Y=lm(NEU_pcent~drinking_habit+L.Pipecolic.acid+group+gender+age+BMI+TNM_stage+medication_history,data)

results = mediate(fit.mediator, fit.Y, treat='drinking_habit', mediator='L.Pipecolic.acid', robustSE = TRUE,sims = 1000)

summary(results)

medsen5<-medsens(results, rho.by = 0.01, sims = 100, effect.type = "indirect")

summary(medsen5)

#6

fit.totaleffect=lm(GLB ~ drinking_habit+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.mediator=lm(L.Pipecolic.acid~drinking_habit+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.Y=lm(GLB~drinking_habit+L.Pipecolic.acid+group+gender+age+BMI+TNM_stage+medication_history,data)

results = mediate(fit.mediator, fit.Y, treat='drinking_habit', mediator='L.Pipecolic.acid', robustSE = TRUE,sims = 1000)

summary(results)

medsen6<-medsens(results, rho.by = 0.01, sims = 100, effect.type = "indirect")

summary(medsen6)

#7

fit.totaleffect=lm(ALT_.AST ~ drinking_habit+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.mediator=lm(Gallic.acid~drinking_habit+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.Y=lm(ALT_.AST~drinking_habit+Gallic.acid+group+gender+age+BMI+TNM_stage+medication_history,data)

results = mediate(fit.mediator, fit.Y, treat='drinking_habit', mediator='Gallic.acid', robustSE = TRUE,sims = 1000)

summary(results)

medsen7<-medsens(results, rho.by = 0.01, sims = 100, effect.type = "indirect")

summary(medsen7)

#8

fit.totaleffect=lm(CRE ~ green_tea_habit+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.mediator=lm(X2.Methylbutyroylcarnitine~green_tea_habit+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.Y=lm(CRE~green_tea_habit+X2.Methylbutyroylcarnitine+group+gender+age+BMI+TNM_stage+medication_history,data)

results = mediate(fit.mediator, fit.Y, treat='green_tea_habit', mediator='X2.Methylbutyroylcarnitine', robustSE = TRUE,sims = 1000)

summary(results)

medsen8<-medsens(results, rho.by = 0.01, sims = 100, effect.type = "indirect")

summary(medsen8)

#9

fit.totaleffect=lm(LYM_pcent ~ green_tea_habit+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.mediator=lm(LPE.20.3.~green_tea_habit+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.Y=lm(LYM_pcent~green_tea_habit+LPE.20.3.+group+gender+age+BMI+TNM_stage+medication_history,data)

results = mediate(fit.mediator, fit.Y, treat='green_tea_habit', mediator='LPE.20.3.', robustSE = TRUE,sims = 1000)

summary(results)

medsen9<-medsens(results, rho.by = 0.01, sims = 100, effect.type = "indirect")

summary(medsen9)

#10

fit.totaleffect=lm(NEU_pcent ~ green_tea_habit+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.mediator=lm(LPE.20.3.~green_tea_habit+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.Y=lm(NEU_pcent~green_tea_habit+LPE.20.3.+group+gender+age+BMI+TNM_stage+medication_history,data)

results = mediate(fit.mediator, fit.Y, treat='green_tea_habit', mediator='LPE.20.3.', robustSE = TRUE,sims = 1000)

summary(results)

medsen10<-medsens(results, rho.by = 0.01, sims = 100, effect.type = "indirect")

summary(medsen10)

#11

fit.totaleffect=lm(UREA ~ frenquency_on_pickled_food+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.mediator=lm(D....Glucose~frenquency_on_pickled_food+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.Y=lm(UREA~frenquency_on_pickled_food+D....Glucose+group+gender+age+BMI+TNM_stage+medication_history,data)

results = mediate(fit.mediator, fit.Y, treat='frenquency_on_pickled_food', mediator='D....Glucose', robustSE = TRUE,sims = 1000)

summary(results)

medsen11<-medsens(results, rho.by = 0.01, sims = 100, effect.type = "indirect")

summary(medsen11)

#12

fit.totaleffect=lm(ALT_.AST ~ daily_meal_times+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.mediator=lm(D....Glucose~daily_meal_times+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.Y=lm(ALT_.AST~daily_meal_times+D....Glucose+group+gender+age+BMI+TNM_stage+medication_history,data)

results = mediate(fit.mediator, fit.Y, treat='daily_meal_times', mediator='D....Glucose', robustSE = TRUE,sims = 1000)

summary(results)

medsen12<-medsens(results, rho.by = 0.01, sims = 100, effect.type = "indirect")

summary(medsen12)


#13

fit.totaleffect=lm(GLB ~ frenquency_on_fresh_fruit_vegetable+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.mediator=lm(ACar.11.1.~frenquency_on_fresh_fruit_vegetable+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.Y=lm(GLB~frenquency_on_fresh_fruit_vegetable+ACar.11.1.+group+gender+age+BMI+TNM_stage+medication_history,data)

results = mediate(fit.mediator, fit.Y, treat='frenquency_on_fresh_fruit_vegetable', mediator='ACar.11.1.', robustSE = TRUE,sims = 1000)

summary(results)

medsen13<-medsens(results, rho.by = 0.01, sims = 100, effect.type = "indirect")

summary(medsen13)

#14

fit.totaleffect=lm(NEU_pcent ~ drinking_habit+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.mediator=lm(DL.2.Aminooctanoic.acid~drinking_habit+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.Y=lm(NEU_pcent~drinking_habit+DL.2.Aminooctanoic.acid+group+gender+age+BMI+TNM_stage+medication_history,data)

results = mediate(fit.mediator, fit.Y, treat='drinking_habit', mediator='DL.2.Aminooctanoic.acid', robustSE = TRUE,sims = 1000)

summary(results)

medsen14<-medsens(results, rho.by = 0.01, sims = 100, effect.type = "indirect")

summary(medsen14)

#15

fit.totaleffect=lm(γ_GGT ~ frenquency_on_pickled_food+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.mediator=lm(Glutaric.acid~frenquency_on_pickled_food+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.Y=lm(γ_GGT~frenquency_on_pickled_food+Glutaric.acid+group+gender+age+BMI+TNM_stage+medication_history,data)

results = mediate(fit.mediator, fit.Y, treat='frenquency_on_pickled_food', mediator='Glutaric.acid', robustSE = TRUE,sims = 1000)

summary(results)

medsen15<-medsens(results, rho.by = 0.01, sims = 100, effect.type = "indirect")

summary(medsen15)

#16

fit.totaleffect=lm(UREA ~ frenquency_on_pickled_food+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.mediator=lm(Glutaric.acid~frenquency_on_pickled_food+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.Y=lm(UREA~frenquency_on_pickled_food+Glutaric.acid+group+gender+age+BMI+TNM_stage+medication_history,data)

results = mediate(fit.mediator, fit.Y, treat='frenquency_on_pickled_food', mediator='Glutaric.acid', robustSE = TRUE,sims = 1000)

summary(results)

medsen16<-medsens(results, rho.by = 0.01, sims = 100, effect.type = "indirect")

summary(medsen16)


#17

fit.totaleffect=lm(GLB ~ daily_meal_times+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.mediator=lm(X2.Hydroxyphenylacetic.acid~daily_meal_times+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.Y=lm(GLB~daily_meal_times+X2.Hydroxyphenylacetic.acid+group+gender+age+BMI+TNM_stage+medication_history,data)

results = mediate(fit.mediator, fit.Y, treat='daily_meal_times', mediator='X2.Hydroxyphenylacetic.acid', robustSE = TRUE,sims = 1000)

summary(results)

medsen17<-medsens(results, rho.by = 0.01, sims = 100, effect.type = "indirect")

summary(medsen17)

#18

fit.totaleffect=lm(LYM_pcent ~ smoking_history+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.mediator=lm(X3.4.Dihydro.2.2.5.7.8.pentamethyl.2H.1.benzopyran.6.ol~smoking_history+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.Y=lm(LYM_pcent~smoking_history+X3.4.Dihydro.2.2.5.7.8.pentamethyl.2H.1.benzopyran.6.ol+group+gender+age+BMI+TNM_stage+medication_history,data)

results = mediate(fit.mediator, fit.Y, treat='smoking_history', mediator='X3.4.Dihydro.2.2.5.7.8.pentamethyl.2H.1.benzopyran.6.ol', robustSE = TRUE,sims = 1000)

summary(results)

medsen18<-medsens(results, rho.by = 0.01, sims = 100, effect.type = "indirect")

summary(medsen18)

#19

fit.totaleffect=lm(LYM_pcent ~ green_tea_habit+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.mediator=lm(X3.4.Dihydro.2.2.5.7.8.pentamethyl.2H.1.benzopyran.6.ol~green_tea_habit+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.Y=lm(LYM_pcent~green_tea_habit+X3.4.Dihydro.2.2.5.7.8.pentamethyl.2H.1.benzopyran.6.ol+group+gender+age+BMI+TNM_stage+medication_history,data)

results = mediate(fit.mediator, fit.Y, treat='green_tea_habit', mediator='X3.4.Dihydro.2.2.5.7.8.pentamethyl.2H.1.benzopyran.6.ol', robustSE = TRUE,sims = 1000)

summary(results)

medsen19<-medsens(results, rho.by = 0.01, sims = 100, effect.type = "indirect")

summary(medsen19)


#20

fit.totaleffect=lm(NEU_pcent ~ green_tea_habit+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.mediator=lm(X3.4.Dihydro.2.2.5.7.8.pentamethyl.2H.1.benzopyran.6.ol~green_tea_habit+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.Y=lm(NEU_pcent~green_tea_habit+X3.4.Dihydro.2.2.5.7.8.pentamethyl.2H.1.benzopyran.6.ol+group+gender+age+BMI+TNM_stage+medication_history,data)

results = mediate(fit.mediator, fit.Y, treat='green_tea_habit', mediator='X3.4.Dihydro.2.2.5.7.8.pentamethyl.2H.1.benzopyran.6.ol', robustSE = TRUE,sims = 1000)

summary(results)

medsen20<-medsens(results, rho.by = 0.01, sims = 100, effect.type = "indirect")

summary(medsen20)

#21

fit.totaleffect=lm(MON_pcent ~ daily_meal_times+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.mediator=lm(PC.20.4.22.6.~daily_meal_times+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.Y=lm(MON_pcent~daily_meal_times+PC.20.4.22.6.+group+gender+age+BMI+TNM_stage+medication_history,data)

results = mediate(fit.mediator, fit.Y, treat='daily_meal_times', mediator='PC.20.4.22.6.', robustSE = TRUE,sims = 1000)

summary(results)

medsen21<-medsens(results, rho.by = 0.01, sims = 100, effect.type = "indirect")

summary(medsen21)

#22

fit.totaleffect=lm(NEU_pcent ~ drinking_habit+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.mediator=lm(Vulgarone.A~drinking_habit+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.Y=lm(NEU_pcent~drinking_habit+Vulgarone.A+group+gender+age+BMI+TNM_stage+medication_history,data)

results = mediate(fit.mediator, fit.Y, treat='drinking_habit', mediator='Vulgarone.A', robustSE = TRUE,sims = 1000)

summary(results)

medsen22<-medsens(results, rho.by = 0.01, sims = 100, effect.type = "indirect")

summary(medsen22)

#23

fit.totaleffect=lm(GLB ~ drinking_habit+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.mediator=lm(Vulgarone.A~drinking_habit+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.Y=lm(GLB~drinking_habit+Vulgarone.A+group+gender+age+BMI+TNM_stage+medication_history,data)

results = mediate(fit.mediator, fit.Y, treat='drinking_habit', mediator='Vulgarone.A', robustSE = TRUE,sims = 1000)

summary(results)

medsen23<-medsens(results, rho.by = 0.01, sims = 100, effect.type = "indirect")

summary(medsen23)

#24

fit.totaleffect=lm(ALT_.AST ~ daily_meal_times+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.mediator=lm(Sphingosine.1.phosphate~daily_meal_times+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.Y=lm(ALT_.AST~daily_meal_times+Sphingosine.1.phosphate+group+gender+age+BMI+TNM_stage+medication_history,data)

results = mediate(fit.mediator, fit.Y, treat='daily_meal_times', mediator='Sphingosine.1.phosphate', robustSE = TRUE,sims = 1000)

summary(results)

medsen24<-medsens(results, rho.by = 0.01, sims = 100, effect.type = "indirect")

summary(medsen24)

#25

fit.totaleffect=lm(MON_pcent ~ drinking_habit+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.mediator=lm(X1.2.3.Trihydroxybenzene~drinking_habit+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.Y=lm(MON_pcent~drinking_habit+X1.2.3.Trihydroxybenzene+group+gender+age+BMI+TNM_stage+medication_history,data)

results = mediate(fit.mediator, fit.Y, treat='drinking_habit', mediator='X1.2.3.Trihydroxybenzene', robustSE = TRUE,sims = 1000)

summary(results)

medsen25<-medsens(results, rho.by = 0.01, sims = 100, effect.type = "indirect")

summary(medsen25)

#26

fit.totaleffect=lm(GLB ~ drinking_habit+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.mediator=lm(X1.2.3.Trihydroxybenzene~drinking_habit+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.Y=lm(GLB~drinking_habit+X1.2.3.Trihydroxybenzene+group+gender+age+BMI+TNM_stage+medication_history,data)

results = mediate(fit.mediator, fit.Y, treat='drinking_habit', mediator='X1.2.3.Trihydroxybenzene', robustSE = TRUE,sims = 1000)

summary(results)

medsen26<-medsens(results, rho.by = 0.01, sims = 100, effect.type = "indirect")

summary(medsen26)



#27

fit.totaleffect=lm(ALT_.AST ~ drinking_habit+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.mediator=lm(X1.2.3.Trihydroxybenzene~drinking_habit+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.Y=lm(ALT_.AST~drinking_habit+X1.2.3.Trihydroxybenzene+group+gender+age+BMI+TNM_stage+medication_history,data)

results = mediate(fit.mediator, fit.Y, treat='drinking_habit', mediator='X1.2.3.Trihydroxybenzene', robustSE = TRUE,sims = 1000)

summary(results)

plot(results)

medsen27<-medsens(results, rho.by = 0.01, sims = 100, effect.type = "indirect")

summary(medsen27)

#28

fit.totaleffect=lm(UREA ~ green_tea_habit+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.mediator=lm(X2.Methylbutyroylcarnitine~green_tea_habit+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.Y=lm(UREA~green_tea_habit+X2.Methylbutyroylcarnitine+group+gender+age+BMI+TNM_stage+medication_history,data)

results = mediate(fit.mediator, fit.Y, treat='green_tea_habit', mediator='X2.Methylbutyroylcarnitine', robustSE = TRUE,sims = 1000)

summary(results)

medsen28<-medsens(results, rho.by = 0.01, sims = 100, effect.type = "indirect")

summary(medsen28)

#29

fit.totaleffect=lm(UREA ~ daily_meal_times+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.mediator=lm(D....Glucose~daily_meal_times+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.Y=lm(UREA~daily_meal_times+D....Glucose+group+gender+age+BMI+TNM_stage+medication_history,data)

results = mediate(fit.mediator, fit.Y, treat='daily_meal_times', mediator='D....Glucose', robustSE = TRUE,sims = 1000)

summary(results)

medsen29<-medsens(results, rho.by = 0.01, sims = 100, effect.type = "indirect")

summary(medsen29)

#30

fit.totaleffect=lm(LYM_pcent~ frenquency_on_fresh_fruit_vegetable+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.mediator=lm(N..Formylkynurenine~frenquency_on_fresh_fruit_vegetable+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.Y=lm(LYM_pcent~frenquency_on_fresh_fruit_vegetable+N..Formylkynurenine+group+gender+age+BMI+TNM_stage+medication_history,data)

results = mediate(fit.mediator, fit.Y, treat='frenquency_on_fresh_fruit_vegetable', mediator='N..Formylkynurenine', robustSE = TRUE,sims = 1000)

summary(results)

medsen30<-medsens(results, rho.by = 0.01, sims = 100, effect.type = "indirect")

summary(medsen30)

#31

fit.totaleffect=lm(NEU_pcent~ frenquency_on_fresh_fruit_vegetable+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.mediator=lm(N..Formylkynurenine~frenquency_on_fresh_fruit_vegetable+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.Y=lm(NEU_pcent~frenquency_on_fresh_fruit_vegetable+N..Formylkynurenine+group+gender+age+BMI+TNM_stage+medication_history,data)

results = mediate(fit.mediator, fit.Y, treat='frenquency_on_fresh_fruit_vegetable', mediator='N..Formylkynurenine', robustSE = TRUE,sims = 1000)

summary(results)

medsen31<-medsens(results, rho.by = 0.01, sims = 100, effect.type = "indirect")

summary(medsen31)



# backward

#1

fit.totaleffect=lm(Androsterone.sulfate ~ smoking_history+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.mediator=lm(ALP~smoking_history+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.Y=lm(Androsterone.sulfate~smoking_history+ALP+group+gender+age+BMI+TNM_stage+medication_history,data)

results = mediate(fit.mediator, fit.Y, treat='smoking_history', mediator='ALP', robustSE = TRUE,sims = 1000)

summary(results)

back_medsen1<-medsens(results, rho.by = 0.01, sims = 100, effect.type = "indirect")

summary(back_medsen1)

#2

fit.totaleffect=lm(L.Methionine~ smoking_history+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.mediator=lm(ALP~smoking_history+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.Y=lm(L.Methionine~smoking_history+ALP+group+gender+age+BMI+TNM_stage+medication_history,data)

results = mediate(fit.mediator, fit.Y, treat='smoking_history', mediator='ALP', robustSE = TRUE,sims = 1000)

summary(results)

back_medsen2<-medsens(results, rho.by = 0.01, sims = 100, effect.type = "indirect")

summary(back_medsen2)

#3

fit.totaleffect=lm(N1.Methyl.2.pyridone.5.carboxamide ~ preference_of_hard_food+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.mediator=lm(MON_pcent~preference_of_hard_food+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.Y=lm(N1.Methyl.2.pyridone.5.carboxamide~preference_of_hard_food+MON_pcent+group+gender+age+BMI+TNM_stage+medication_history,data)

results = mediate(fit.mediator, fit.Y, treat='preference_of_hard_food', mediator='MON_pcent', robustSE = TRUE,sims = 1000)

summary(results)

back_medsen3<-medsens(results, rho.by = 0.01, sims = 100, effect.type = "indirect")

summary(back_medsen3)

#4

fit.totaleffect=lm(LPE.18.1. ~ preference_of_hard_food+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.mediator=lm(ALT_.AST~preference_of_hard_food+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.Y=lm(LPE.18.1.~preference_of_hard_food+ALT_.AST+group+gender+age+BMI+TNM_stage+medication_history,data)

results = mediate(fit.mediator, fit.Y, treat='preference_of_hard_food', mediator='ALT_.AST', robustSE = TRUE,sims = 1000)

summary(results)

back_medsen4<-medsens(results, rho.by = 0.01, sims = 100, effect.type = "indirect")

summary(back_medsen4)


#5

fit.totaleffect=lm(LPC.16.0. ~ preference_of_hard_food+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.mediator=lm(MON_pcent~preference_of_hard_food+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.Y=lm(LPC.16.0.~preference_of_hard_food+MON_pcent+group+gender+age+BMI+TNM_stage+medication_history,data)

results = mediate(fit.mediator, fit.Y, treat='preference_of_hard_food', mediator='MON_pcent', robustSE = TRUE,sims = 1000)

summary(results)

back_medsen5<-medsens(results, rho.by = 0.01, sims = 100, effect.type = "indirect")

summary(back_medsen5)

#6

fit.totaleffect=lm(LPG.16.0. ~ preference_of_hard_food+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.mediator=lm(MON_pcent~preference_of_hard_food+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.Y=lm(LPG.16.0.~preference_of_hard_food+MON_pcent+group+gender+age+BMI+TNM_stage+medication_history,data)

results = mediate(fit.mediator, fit.Y, treat='preference_of_hard_food', mediator='MON_pcent', robustSE = TRUE,sims = 1000)

summary(results)

back_medsen6<-medsens(results, rho.by = 0.01, sims = 100, effect.type = "indirect")

summary(back_medsen6)

#7

fit.totaleffect=lm(LysoPE.18.1.9Z..0.0. ~ preference_of_hard_food+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.mediator=lm(ALT_.AST~preference_of_hard_food+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.Y=lm(LysoPE.18.1.9Z..0.0.~preference_of_hard_food+ALT_.AST+group+gender+age+BMI+TNM_stage+medication_history,data)

results = mediate(fit.mediator, fit.Y, treat='preference_of_hard_food', mediator='ALT_.AST', robustSE = TRUE,sims = 1000)

summary(results)

back_medsen7<-medsens(results, rho.by = 0.01, sims = 100, effect.type = "indirect")

summary(back_medsen7)


#8

fit.totaleffect=lm(Lysyl.Alanine ~ frenquency_on_fresh_fruit_vegetable+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.mediator=lm(ALT_.AST~frenquency_on_fresh_fruit_vegetable+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.Y=lm(Lysyl.Alanine~frenquency_on_fresh_fruit_vegetable+ALT_.AST+group+gender+age+BMI+TNM_stage+medication_history,data)

results = mediate(fit.mediator, fit.Y, treat='frenquency_on_fresh_fruit_vegetable', mediator='ALT_.AST', robustSE = TRUE,sims = 1000)

summary(results)

back_medsen8<-medsens(results, rho.by = 0.01, sims = 100, effect.type = "indirect")

summary(back_medsen8)


#9

fit.totaleffect=lm(LysoPC.18.1.9Z.. ~ preference_of_hard_food+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.mediator=lm(ALT_.AST~preference_of_hard_food+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.Y=lm(LysoPC.18.1.9Z..~preference_of_hard_food+ALT_.AST+group+gender+age+BMI+TNM_stage+medication_history,data)

results = mediate(fit.mediator, fit.Y, treat='preference_of_hard_food', mediator='ALT_.AST', robustSE = TRUE,sims = 1000)

summary(results)

back_medsen9<-medsens(results, rho.by = 0.01, sims = 100, effect.type = "indirect")

summary(back_medsen9)


## forward:about drinking habit and immune function


#4

fit.totaleffect1=lm(LYM_pcent ~ drinking_habit+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.mediator1=lm(L.Pipecolic.acid~drinking_habit+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.Y1=lm(LYM_pcent~drinking_habit+L.Pipecolic.acid+group+gender+age+BMI+TNM_stage+medication_history,data)

results1 = mediate(fit.mediator1, fit.Y1, treat='drinking_habit', mediator='L.Pipecolic.acid', robustSE = TRUE,sims = 1000)

summary(results1)


#5

fit.totaleffect2=lm(NEU_pcent ~ drinking_habit+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.mediator2=lm(L.Pipecolic.acid~drinking_habit+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.Y2=lm(NEU_pcent~drinking_habit+L.Pipecolic.acid+group+gender+age+BMI+TNM_stage+medication_history,data)

results2 = mediate(fit.mediator2, fit.Y2, treat='drinking_habit', mediator='L.Pipecolic.acid', robustSE = TRUE,sims = 1000)

summary(results2)


#6

fit.totaleffect3=lm(GLB ~ drinking_habit+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.mediator3=lm(L.Pipecolic.acid~drinking_habit+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.Y3=lm(GLB~drinking_habit+L.Pipecolic.acid+group+gender+age+BMI+TNM_stage+medication_history,data)

results3 = mediate(fit.mediator3, fit.Y3, treat='drinking_habit', mediator='L.Pipecolic.acid', robustSE = TRUE,sims = 1000)

summary(results3)



#14

fit.totaleffect4=lm(NEU_pcent ~ drinking_habit+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.mediator4=lm(DL.2.Aminooctanoic.acid~drinking_habit+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.Y4=lm(NEU_pcent~drinking_habit+DL.2.Aminooctanoic.acid+group+gender+age+BMI+TNM_stage+medication_history,data)

results4 = mediate(fit.mediator4, fit.Y4, treat='drinking_habit', mediator='DL.2.Aminooctanoic.acid', robustSE = TRUE,sims = 1000)

summary(results4)



#22

fit.totaleffect5=lm(NEU_pcent ~ drinking_habit+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.mediator5=lm(Vulgarone.A~drinking_habit+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.Y5=lm(NEU_pcent~drinking_habit+Vulgarone.A+group+gender+age+BMI+TNM_stage+medication_history,data)

results5 = mediate(fit.mediator5, fit.Y5, treat='drinking_habit', mediator='Vulgarone.A', robustSE = TRUE,sims = 1000)

plot(results5)

summary(results5)


#23

fit.totaleffect6=lm(GLB ~ drinking_habit+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.mediator6=lm(Vulgarone.A~drinking_habit+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.Y6=lm(GLB~drinking_habit+Vulgarone.A+group+gender+age+BMI+TNM_stage+medication_history,data)

results6 = mediate(fit.mediator6, fit.Y6, treat='drinking_habit', mediator='Vulgarone.A', robustSE = TRUE,sims = 1000)

plot(results6)

summary(results6)



#25

fit.totaleffect7=lm(MON_pcent ~ drinking_habit+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.mediator7=lm(X1.2.3.Trihydroxybenzene~drinking_habit+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.Y7=lm(MON_pcent~drinking_habit+X1.2.3.Trihydroxybenzene+group+gender+age+BMI+TNM_stage+medication_history,data)

results7 = mediate(fit.mediator7, fit.Y7, treat='drinking_habit', mediator='X1.2.3.Trihydroxybenzene', robustSE = TRUE,sims = 1000)

plot(results7)

summary(results7)

#26

fit.totaleffect8=lm(GLB ~ drinking_habit+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.mediator8=lm(X1.2.3.Trihydroxybenzene~drinking_habit+group+gender+age+BMI+TNM_stage+medication_history,data)

fit.Y8=lm(GLB~drinking_habit+X1.2.3.Trihydroxybenzene+group+gender+age+BMI+TNM_stage+medication_history,data)

results8 = mediate(fit.mediator8, fit.Y8, treat='drinking_habit', mediator='X1.2.3.Trihydroxybenzene', robustSE = TRUE,sims = 1000)

plot(results8)

summary(results8)


