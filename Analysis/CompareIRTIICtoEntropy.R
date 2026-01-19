### Compare IRT to my "entropy" regression approach:


# 0.0 Setup ####

include(psych)

if(!exists("spaAcc")) source("Analysis/SPA.DataImport.R")

spaAll = t(spaAcc %>% select(-ItemType, -Stimulus)) %>% 
  data.frame %>% 
  mutate(across(everything(), as.numeric))
colnames(spaAll) = spaAcc$Stimulus

itemCols = colnames(spaAll)

# 1.0 fit the IRT ####

spaAll.irt = irt.fa(spaAll)

# 1.1 use IRT difficulties as predictors for items ###

spaAll$ability = as.vector(factor.scores(spaAll, spaAll.irt$fa)$scores)

# 2.0 Entropy from abilities ####
spa.ability.glm = 
  apply(spaAll[,itemCols], 2
        , FUN=function(x) 
          glm(acc~ability, data=data.frame(acc=x, ability=spaAll$ability)
              , family=binomial))

## 2.1 "predict" the performance on each item for abilities of -2, -1.5,..., +2 ####
spaitems.pred = 
  lapply(spa.ability.glm
         , function(x) predict(x, data.frame(ability=-3:3)
                               , type="response")) %>% bind_cols()

## 2.2 calculate "info" from entropy ####
spaitems.info = 
  lapply(spaitems.pred
         , function(x) -x*log(x)-(1-x)*log(1-x)) %>% 
  bind_cols() %>% t() %>% data.frame

colnames(spaitems.info)=c(gsub("p-", "n", paste0("p",-3:3)))

# 3.0 Compare item Infos ####

irtIIC = plot(spaAll.irt)$sumInfo[[1]]

all(rownames(spaitems.info)==rownames(irtIIC))

cor(irtIIC, spaitems.info)

### Ok, so entropy does not correlate quite as well as I would have liked.

# 4.0 Use the Rasch model directly ####
rasch.pred = spaAll %>% slice(1:7)
rasch.pred$ability = seq(-3,3)

for(i in itemCols)
  rasch.pred[,i] =
  1/(
    1+
      exp(-1.7*spaAll.irt$irt$discrimination[i,]*
           (rasch.pred$ability - spaAll.irt$irt$difficulty[[1]][i])
          )
   )

rasch.pred = t(rasch.pred %>% select(-ability))
colnames(rasch.pred) = -3:3
rasch.info = rasch.pred
for(i in 1:ncol(rasch.pred))
  rasch.info[,i] = (1.7*spaAll.irt$irt$discrimination)^2*rasch.pred[,i]*(1-rasch.pred[,i])

cor(irtIIC, rasch.info)
