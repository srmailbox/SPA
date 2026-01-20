###############
# SPA Project - SPA test item analysis - Age Screeners ####
# 
# Final results here are just based on simple tetrachoric correlations by Age Group
#
# Created: 2024-10-07
# Changelog:
# 2025-05-14: Use 20th percentile instead of 16th.

# 0. Set up ####
include(CTT)
include(psych)
include(zoo)

if(!exists("spaAcc")) source("Analysis/SPA.DataImport.R")

# 1.0 reshape data and split by item type ####

## 1.1 all items, regardless of type ####

spaAll = t(spaAcc %>% select(-ItemType, -Stimulus)) %>% 
  data.frame %>% 
  mutate(across(everything(), as.numeric))
colnames(spaAll) = spaAcc$Stimulus

## 1.2 items by type ####
spaItemTypes = by(spaAcc, spaAcc$ItemType
                  , function(x) {
                    temp = t(x %>% select(-ItemType, -Stimulus)) %>% 
                      data.frame %>% mutate(across(everything(), as.numeric))
                    colnames(temp) = x$Stimulus
                    return(temp)
                  }
)


# # 2.0 Classical Item Analysis ####
# 
# ## 2.1 Full Set of Items ####
# 
# spaAll.ctt = itemAnalysis(spaAll)
# 
# 
# ## 2.2 Item Type ####
# 
# spaCTT = lapply(spaItemTypes, itemAnalysis)
# 
# 
# ### That's all fine, but it's not that informative since we are much more
# # concerned with items that discriminate at the lower end of the distribution
# 
# # 3.0 IRT ####
# 
# ## 3.1 Full Set of Items ####
# 
spaAll.irt = irt.fa(spaAll)
# 
# spaAll.f3 = irt.fa(spaAll, 3, plot=F)
# # I don't love how these are falling out
# 
# ## 3.2 Item Type ####
# 
# spaIRT = lapply(spaItemTypes, irt.fa, plot=F)
# 
# par(mfrow=c(2,2))
# spaIRT$Irreg %>% plot(main="Irr")
# spaIRT$Reg %>% plot(main="Reg")
# spaIRT$NW %>% plot(main="NW")
# spaAll.irt %>% plot(main="All")
# par(mfrow=c(1,1))
# 
### 3.2.1 IRT by Item Type vs All together ####
# 
# # How well do the Information matrices from the overall analysis match up
# # with the analysis by each item type. (e.g., can we just use the overall
# # IRT to identify useful items, or do we need to consider each itemType 
# # separately)
# spaIIC = 
#   merge(
#     lapply(spaIRT
#            , function(x) plot(x)$sumInfo[[1]] %>% data.frame) %>% 
#       bind_rows(.id="Type") %>% 
#       mutate(Item = rownames(.))
#     , plot(spaAll.irt)$sumInfo[[1]] %>% data.frame %>% mutate(Item=rownames(.))
#     , by="Item"
#     , suffixes=c(".byType", ".All")
#   )
# 
# cor(spaIIC %>% select(ends_with(".byType"))
#     , spaIIC %>% select(ends_with(".All"))) %>% diag
# 
## Huh! Turns out the overall analysis does really well, so we can maybe take
# an easier approach here.

# Let's add the kids' ability to their data.
itemCols = colnames(spaAll)
spaAll$ability =
  factor.scores(spaAll, spaAll.irt$fa)$scores
spaResults = merge(spaAll %>% mutate(newID = rownames(.)
                                     , across(all_of(itemCols), ~ifelse(is.na(.), 0, .)))
                   , pptDetails %>% select(newID, starts_with("sp_")) %>%
                     # rename(newID = Child_ID, sp_age=sp.age) %>%
                     mutate(across(starts_with("sp_"), as.numeric))
                   , by = "newID") %>%
  # rowwise() %>%
  # The pptDetails scores just flat out wrong, so this recalculates those scores
  # based on item level scoring.  
  mutate(sp_alt = rowSums(select(., all_of(itemCols)))
         , sp_irr_alt = rowSums(select(., all_of(itemCols[1:29])))
         , sp_reg_alt = rowSums(select(., all_of(itemCols[1:29+29])))
         , sp_nw_alt = rowSums(select(., all_of(itemCols[1:29+58])))
         )
#   
# 
# plot(spaResults$sp_age, spaResults$ability)
# # plot(spaResults$nst_total, spaResults$ability)
# # plot(spaResults$sp_age, spaResults$nst_total)
# ### The problem here is that by lumping everyone together, we assume they should
# # have similar abilities. But younger kids have lower abilities, as expected.
# 
# # 4.0 Adjusting for Age ####
# # Ok, so the only issue here is that ability is correlated with age in a way
# # that will undermine our results.
# 
# # I think we want to redo this with the following logic:
# # 1. filter out kids under 78 months (6.5yrs) because the really young kids
# #    just can't be evaluated with the same items as the 12 year olds
# # 2. Apply a full IRT model to determine a general ability score.
# # 3. Adjust the ability score for age in months with regression, to get 
# #    an age-corrected ability (standardized residuals)
# # 4. "predict" the performance on each item for abilities of -2, -1.5,..., +2
# #     based on the IRT model
# # 5. Calculate information based on the IRT model
# 
# # Then we can use this to determine a set of items that have high entropy for
# # each level of age-corrected ability.
# 
## 4.1 Age Filter ####
# Filtering kids under 6.5 yrs old, and removing a single outlier 18 yr old.
spa7plus = spaResults %>% filter(sp_age>=78, sp_age < 200)

# # Leaves us with 770 kids
# 
# ## 4.2 IRT on the restricted set to get abilities ####
spa7.irt = irt.fa(spa7plus[,itemCols])
# spa7plus$ability = factor.scores(spa7plus[,itemCols], spa7.irt$fa)$scores
# 
# ## 4.3 Age adjusted abilities ####
# spa7plus$ability.age= rstandard(lm(ability~sp_age, spa7plus))
# 
# # ggplot(spa7plus, aes(x=sp_age, y=ability))+
# #   theme_minimal()+
# #   geom_point(aes(y=ability.age), col="red", alpha=.5)+
# #   geom_smooth(aes(y=ability.age), col="red")
# # 
# # ggplot(spa7plus, aes(x=sp_age, y=ability))+
# #   theme_minimal()+
# #   geom_point(col="grey")+
# #   geom_smooth(col="grey")
# 
# ## 4.4 logistic regression for each item against ability ####
# # spa7.ability.glm = 
# #   apply(spa7plus[,itemCols], 2
# #         , FUN=function(x) 
# #           glm(acc~ability, data=data.frame(acc=x, ability=spa7plus$ability.age)
# #               , family=binomial))
# 
# ## 4.4 "predict" the performance on each item for abilities of -2, -1.5,..., +2 ###
# spa7items.pred = spa7plus %>% slice(1:length(seq(-2,2,.25))) %>% 
#   select(all_of(itemCols))
# abilities = seq(-2,2,.25)
# 
# for(i in itemCols) {
#   spa7items.pred[,i] =
#   1/(
#     1+
#       exp(-1.7*spa7.irt$irt$discrimination[i,]*
#             (abilities - spa7.irt$irt$difficulty[[1]][i])
#       )
#   )
# }
# 
# 
# ## 4.5 calculate "info" from entropy ####
# spa7items.info = 
#   lapply(spa7items.pred
#          , function(x) -x*log(x)-(1-x)*log(1-x)) %>% 
#   bind_cols() %>% t() %>% data.frame %>% 
#   mutate(Type=gl(3,29,labels=c("Irr", "Reg", "NW")))
# 
# colnames(spa7items.info)=c(gsub("p-", "n", paste0("p",seq(-2,2,.25))), "Type")
# # colnames(spa7items.info)[9]="n0"
# 
# ## 4.6 Item Informativeness ####
# 
# info.n.5 = spa7items.info %>% select(n0.5, Type) %>% rename(Info = n0.5) %>% 
#   arrange(desc(Info)) %>%
#   mutate(items = factor(rownames(.), levels=rownames(.)[order(-Info)])) %>% 
#   ggplot(aes(x=items, y=Info, fill=Type))+
#   geom_col()+
#   theme_minimal()+
#   theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust=1))+
#   ggtitle("How informative is each item around ability = -.5 SD?") +
#   ylab("Information")+xlab(NULL)
# 
# info.n1 = spa7items.info %>% select(n1, Type) %>% rename(Info = n1) %>% 
#   arrange(desc(Info)) %>%
#   mutate(items = factor(rownames(.), levels=rownames(.)[order(-Info)])) %>% 
#   ggplot(aes(x=items, y=Info, fill=Type))+
#   geom_col()+
#   theme_minimal()+
#   theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust=1))+
#   ggtitle("How informative is each item around ability = -1 SD?") +
#   ylab("Information")+xlab(NULL)
# 
# 
# gridExtra::grid.arrange(
#   info.n.5, info.n1
# )
# ### This just doesn't work


# 5.0 Item Performance by Age ####

spa7.item.glm = 
  apply(spa7plus[,itemCols], 2
        , function(x) {
          mod = glm(Acc~scale(sp_age), data=data.frame(Acc=x, sp_age=spa7plus$sp_age)
              , family=binomial)
          list(glm=mod, performance = data.frame(P = predict(mod, data.frame(sp_age=78:144), type="response")
                     , sp_age=78:144)
               , p50age = data.frame(p50=-coef(mod)[1]/coef(mod)[2]))
        })

p50s = lapply(spa7.item.glm, function(x) x$p50age) %>% bind_rows(.id="item")
rownames(p50s)=NULL

ageICC = data.frame(
  lapply(spa7.item.glm, function(x) x$performance[,1]) %>% bind_cols()
  , age=78:144
  )

ageICC %>% 
  pivot_longer(-age, names_to="item", values_to="pCorrect") %>% 
  ggplot(aes(x=age, y=pCorrect, group=item))+
  geom_line()+
  theme_minimal()+
  theme(panel.grid.major.x = element_blank()
        , panel.grid.minor.x = element_blank()
        )

itemCoefs = lapply(spa7.item.glm, function(x) coef(x$glm)) %>% 
  bind_rows(.id="item") %>% 
  rename(intrcpt="(Intercept)", age="scale(sp_age)") %>%  
  arrange(age) %>% 
  # filter(intrcpt<0) %>% 
  data.frame

low15slopes = itemCoefs %>% arrange(age) %>% slice(1:15) %>% select(item)

ageICC.shallowSlope=ageICC %>% 
  pivot_longer(-age, names_to="item", values_to="pCorrect") %>% 
  filter(item %in% low15slopes$item)

ggplot(ageICC.shallowSlope, aes(x=age, y=pCorrect, group=item, label=item))+
  geom_line()+
  geom_label(data=ageICC.shallowSlope %>% group_by(item) %>% sample_n(1)
             , alpha = .8)+
  theme_minimal()+
  theme(panel.grid.major.x = element_blank()
        , panel.grid.minor.x = element_blank()
  )


mostAverage = itemCoefs %>% arrange(intrcpt^2) %>% slice(1:15) %>% select(item)

ageICC.average=ageICC %>% 
  pivot_longer(-age, names_to="item", values_to="pCorrect") %>% 
  filter(item %in% mostAverage$item)

ggplot(ageICC.average, aes(x=age, y=pCorrect, group=item, label=item))+
  geom_line()+
  geom_label(data=ageICC.average %>% group_by(item) %>% sample_n(1)
             , alpha = .8)+
  theme_minimal()+
  theme(panel.grid.major.x = element_blank()
        , panel.grid.minor.x = element_blank()
  )


# 6.0 Item Performance by Ability (unadj) ####

spa7.ability.glm = 
  apply(spa7plus[,itemCols], 2
        , function(x) {
          mod = glm(Acc~ability
                    , data=data.frame(
                      Acc=x, ability=as.vector(spa7plus$ability)
                    )
                    , family=binomial)
          list(glm=mod
               , performance = 
                 data.frame(
                   P = predict(mod
                               , data.frame(
                                 ability=seq(min(spa7plus$ability)
                                             , max(spa7plus$ability)
                                             , length.out=100))
                               , type="response")
                   , ability=seq(
                     min(spa7plus$ability), max(spa7plus$ability)
                     , length.out=100)
                   )
               , p50ability = data.frame(p50=-coef(mod)[1]/coef(mod)[2]))
        })

p50s.ability = lapply(spa7.ability.glm, function(x) x$p50ability) %>% bind_rows(.id="item")
rownames(p50s.ability)=NULL

ageICC.ability = data.frame(
  lapply(spa7.ability.glm, function(x) x$performance[,1]) %>% bind_cols()
  , abilities = seq(
    min(spa7plus$ability), max(spa7plus$ability)
    , length.out=100)
)

ageICC.ability %>% 
  pivot_longer(-abilities, names_to="item", values_to="pCorrect") %>% 
  ggplot(aes(x=abilities, y=pCorrect, group=item))+
  geom_line()+
  theme_minimal()+
  theme(panel.grid.major.x = element_blank()
        , panel.grid.minor.x = element_blank()
  )

itemCoefs = lapply(spa7.ability.glm, function(x) coef(x$glm)) %>% 
  bind_rows(.id="item") %>% 
  rename(intrcpt="(Intercept)") %>%  
  arrange(ability) %>% 
  # filter(intrcpt<0) %>% 
  data.frame

low15slopes = itemCoefs %>% arrange(ability) %>% slice(1:15) %>% select(item)

ageICC.shallowSlope=ageICC.ability %>% 
  pivot_longer(-abilities, names_to="item", values_to="pCorrect") %>% 
  filter(item %in% low15slopes$item)

ggplot(ageICC.shallowSlope, aes(x=abilities, y=pCorrect, group=item, label=item))+
  geom_line()+
  geom_label(data=ageICC.shallowSlope %>% group_by(item) %>% sample_n(1)
             , alpha = .8)+
  theme_minimal()+
  theme(panel.grid.major.x = element_blank()
        , panel.grid.minor.x = element_blank()
  )


mostAverage = itemCoefs %>% arrange(intrcpt^2) %>% slice(1:15)

ageICC.average=ageICC.ability %>% 
  pivot_longer(-abilities, names_to="item", values_to="pCorrect") %>% 
  filter(item %in% mostAverage$item)

ggplot(ageICC.average, aes(x=abilities, y=pCorrect, group=item, label=item))+
  geom_line()+
  geom_label(data=ageICC.average %>% group_by(item) %>% sample_n(1)
             , alpha = .8)+
  theme_minimal()+
  theme(panel.grid.major.x = element_blank()
        , panel.grid.minor.x = element_blank()
  )

geoMean = itemCoefs %>% arrange(sqrt(abs(intrcpt)*abs(ability))) %>% slice(1:15)

ageICC.geo=ageICC.ability %>% 
  pivot_longer(-abilities, names_to="item", values_to="pCorrect") %>% 
  filter(item %in% geoMean$item)

ggplot(ageICC.geo, aes(x=abilities, y=pCorrect, group=item, label=item))+
  geom_line()+
  geom_label(data=ageICC.geo %>% group_by(item) %>% sample_n(1)
             , alpha = .8)+
  theme_minimal()+
  theme(
    # panel.grid.major.x = element_blank()
        , panel.grid.minor.x = element_blank()
  )

# 7.0 What if we do just use the IRT Results ####

itemCharacteristics = data.frame(plot(spa7.irt)$sumInfo %>% data.frame
           , Diff= spa7.irt$irt$difficulty[[1]]
           , Disc = as.vector(spa7.irt$irt$discrimination))

midDiffItems =itemCharacteristics %>% arrange(Diff) %>% slice(34:53) %>% rownames()

spaMidDiff = spa7plus[, c(midDiffItems, "sp_age")]
spaMidDiff$Score = rowSums(spaMidDiff %>% select(-sp_age))
spaMidDiff$logisScore = qlogis((spaMidDiff$Score+1)/22)

## 7.1 Quant Reg ####

include(quantreg)
spaMidDiff.qr = rq(logisScore~sp_age, data=spaMidDiff, tau=c(.1,.2, .3, .4))

spaMidDiff.qrpred = data.frame(age=78:144
                               ,  plogis(predict(spaMidDiff.qr
                                                 , data.frame(sp_age = 78:144))
                                         )
)


spaMidDiff.qrpred %>% pivot_longer(-age, values_to="P", names_to="tau") %>% 
  mutate(tau = gsub("tau..", "", tau)) %>% 
  ggplot(aes(x=age, y=P*20, color = tau))+
  geom_line()+
  theme_minimal()+
  labs(x="Age in Months", y="Items correct of 20")

# Looks like that does an ok job. even poorer young spellers can spell a bunch
# and poorer older readers are not at ceiling.

## 7.2 QR for IDing "20th percentile" ####

# ggplot(spa7plus, aes(x=sp_alt, y=ability, col=sp_age))+geom_point()
spa7plus$sp_logit = qlogis((spa7plus$sp_alt+1)/89) 
spa7.qr = rq(sp_logit~sp_age, spa7plus, tau=.20)

spa7plus$qr.fit = 87*plogis(predict(spa7.qr))

qrMod = data.frame(sp_age=78:144) %>% 
  mutate(sp_alt=87*plogis(predict(spa7.qr, data.frame(sp_age = 78:144))))

ggplot(spa7plus, aes(x=sp_age, y=sp_alt))+
  geom_point(col="darkgrey")+
  # geom_smooth()+
  geom_line(data=qrMod, aes(colour=NULL), linewidth=1, col="blue")+
  theme_minimal()+
  labs(x="Age in months", y="Score (of 87)"
       , title="20th Percentile Performance", subtitle="as a function of age")

# So, let's flag the poorest 20% of SPA spellers (age adjusted)

spa7plus = spa7plus %>% 
  mutate(sub20p = ifelse(sp_alt < qr.fit, 1, 0))
# table(spa7plus$sub20p)/nrow(spa7plus)

## 7.3 Calculate MidItems score ####

spa7plus = spa7plus %>% 
  mutate(screener = rowSums(select(., all_of(midDiffItems)))
         , ageCat = cut(sp_age, breaks=c(0,89, 101, 113, 125, 150)
                        , labels=c("<= 7;5", paste0("<= ",8:10,";",5), "10;6+")))

## 7.4 Assess performance of screeners ####

### 7.4.1 Overall, ignoring ages ####
nSub20 = sum(spa7plus$sub20p)
spa7roc.overall = spa7plus %>% 
  arrange(screener) %>% 
  mutate(detection = cumsum(sub20p)/nSub20) %>% 
  group_by(screener) %>% 
  summarise(detection = max(detection))
spa7roc.overall

ggplot(spa7roc.overall, aes(x=screener, y=detection))+
  geom_line(col="red")+
  ylab("Detect sub 20th percentile") + xlab("Screener threshold")+
  theme_minimal()

# Would need to apply a rule of 18 or less to catch 100% of sub20ps

### 7.4.2 Thresholds by Age ####
spa7detection.age = spa7plus %>% 
  group_by(ageCat, screener) %>% 
  summarise(detection=sum(sub20p)) %>% #slice(1:20) %>% data.frame()
  # ungroup(screener) %>% 
  mutate(detection=cumsum(detection)/sum(detection)) %>% 
  ungroup() %>% 
  arrange(ageCat, screener)

ggplot(spa7detection.age, aes(x=screener, y=detection, col=ageCat))+
  geom_line()+
  labs(title = "Screener's rate of detecting Full SPA 20th percentile kids"
       , subtitle = "by ages and threshold score (using IRT)"
       , colour = "Age Group"
       , x = "Screener threshold"
       , y = "Detect sub 20th percentile") +
  theme_minimal()+
  scale_x_continuous(breaks = 0:20)+
  scale_y_continuous(minor_breaks = 0:20/20)+
  theme(panel.grid.minor.x=element_blank())

### Ok, so with this we could apply different thresholds:
# sub 7;5: 8 or less indicates follow up
# 7;6 - 8;5: 9 or less
# 8;6 - 9;5: 10 or less
# 9;5 - 10;5: 12 or less
# 10;6 plus: 18 or less
thrshHlds = c(8,9,10,12,18)
# False positive rate
spa7fpr.age = spa7plus %>% 
  group_by(ageCat, screener) %>%
  filter(!sub20p) %>% 
  summarise(n = n()) %>% #slice(1:20) %>% data.frame()
  # ungroup(screener) %>%
  mutate(fpr=cumsum(n)) %>% 
  mutate(fpr=fpr/max(fpr)) %>% 
  ungroup() %>% 
  arrange(ageCat, screener)

ggplot(spa7fpr.age, aes(x=screener, y=fpr, col=ageCat))+
  geom_vline(data=data.frame(ageCat = levels(spa7plus$ageCat)
                             , threshold = thrshHlds)
             , mapping=aes(xintercept=threshold, col=ageCat)
             , size=.25, linetype="dashed")+
  geom_line()+
  labs(title = "False Positive Rate (Overdiagnosis)"
       , subtitle = "by ages and threshold score"
       , colour = "Age Group"
       , x = "Screener threshold"
       , y = "FPR") +
  theme_minimal()+
  scale_x_continuous(breaks = 0:length(midDiffItems))+
  scale_y_continuous(minor_breaks = 0:length(midDiffItems)/length(midDiffItems))+
  theme(panel.grid.minor.x=element_blank())

## the false positive rates are pretty high, especially for the youngest and
# oldest kids. Not too bad in the mid-range.
# age <= 7;5: 35% false positives
# - 8;5: 19.4%
# - 9;5: 10.8%
# - 10;5: 6.6%
# 10;6 +: 57.6% false positives
# so that approach produces good sensitivity, but very poor specificity.

# 8.0 Using tetrachoric cors to identify item lists ####

# In this case, the idea is to identify which items are most predictive of
# whether a kid is likely to fall in the "sub20p" category.

spa7tetrachoric = cor(spa7plus %>% select(all_of(itemCols)), spa7plus$sub20p) %>% data.frame %>% 
  rename(cor = ".") %>% arrange(cor) %>% slice(1:20)

## 8.1 by age ####
spa7tetrachoric.age = by(spa7plus, INDICES=spa7plus$ageCat
   , function(x) cor(x %>% select(all_of(itemCols)), x$sub20p) %>% data.frame %>% 
     rename(cor = ".")# %>% arrange(cor)
   , simplify = F) %>% 
  lapply(FUN=as.data.frame) %>% bind_cols() 

colnames(spa7tetrachoric.age) = paste0("cor", 7:11)

### 8.1.1 age-based screeners ####

spa7screeners=apply(
  spa7tetrachoric.age, 2
  , function(x) {
    spa7plus %>% 
      mutate(screener = rowSums(select(., all_of(names(x[order(x)])[1:20])))) %>% 
      select(newID, screener, ageCat, sub20p)
  }) %>% bind_cols() %>%
  rename(newID = newID...1, sub20p=sub20p...4, ageCat = ageCat...3
         , screener7 = screener...2, screener8 = screener...6
         , screener9 = screener...10, screener10 = screener...14
         , screener11 = screener...18) %>% 
  select(newID, ageCat, sub20p, starts_with("screener")
         , -starts_with("screener...")) %>% 
  mutate(screener.age = 
           case_when(ageCat == "<= 7;5"~screener7, ageCat=="<= 8;5"~screener8
                     , ageCat=="<= 9;5"~screener9, ageCat=="<= 10;5"~screener10
                     , ageCat=="10;6+"~screener11)
         , screener.IRTflag = as.numeric(screener.age <= thrshHlds[as.numeric(ageCat)])
  )

by(spa7screeners, INDICES=spa7screeners$ageCat
   , FUN=function(x) cor(x$sub20p, x$screener.age))

# Those aren't bad

## 8.2 assess age-based screeners ####

spa7detection.agescreeners = spa7screeners %>% 
  group_by(ageCat, screener.age) %>% 
  summarise(detection=sum(sub20p)) %>% #slice(1:20) %>% data.frame()
  # ungroup(screener) %>% 
  mutate(detection=cumsum(detection)/sum(detection)) %>% 
  ungroup() %>% 
  arrange(ageCat, screener.age)

### 8.2.1 Detection ####
ggplot(spa7detection.agescreeners, aes(x=screener.age, y=detection, col=ageCat))+
  geom_line()+
  labs(title = "Screener's rate of detecting Full SPA 20th percentile kids"
       , subtitle = "by ages and threshold score (using tetrachoric correlation)"
       , colour = "Age Group"
       , x = "Screener threshold"
       , y = "Detect sub 20th percentile") +
  theme_minimal()+
  scale_x_continuous(breaks = 0:20)+
  scale_y_continuous(minor_breaks = 0:20/20)+
  theme(panel.grid.minor.x=element_blank())

thrshHlds.ttrchrc = spa7detection.agescreeners %>%
  group_by(ageCat) %>% 
  filter(detection==1) %>% 
  summarize(thrshHld.tc = min(screener.age))
### Ok, so with different items for each age, we could use the following
# sub 7;5: 11 or less indicates follow up
# 7;6 - 8;5: 13
# 8;6 - 9;5: 14 or less
# 9;6 - 10;5: 17 or less
# 10;6 plus: 15 or less
# 9;6-10;5 has a higher threshold than the older 10;6+ 

### 8.2.2 False positive rate ####

spa7screeners = 
  merge(spa7screeners %>% data.frame, thrshHlds.ttrchrc %>% data.frame) %>% 
  mutate(
    screener.tcflag = 
      as.numeric(screener.age <= thrshHld.tc)
    )

spa7fpr.agescreener = spa7screeners %>% 
  group_by(ageCat, screener.age) %>%
  filter(!sub20p) %>% 
  summarise(n = n()) %>% #slice(1:20) %>% data.frame()
  # ungroup(screener) %>%
  mutate(fpr=cumsum(n)) %>% 
  mutate(fpr=fpr/max(fpr)) %>% 
  ungroup() %>% 
  arrange(ageCat, screener.age)

ggplot(spa7fpr.agescreener, aes(x=screener.age, y=fpr, col=ageCat))+
  geom_vline(data=thrshHlds.ttrchrc
             , mapping=aes(xintercept=thrshHld.tc, col=ageCat)
             , size=.25, linetype="dashed")+
  geom_line()+
  labs(title = "False Positive Rate (Overdiagnosis)"
       , subtitle = "by ages and threshold score - Age-specific Items"
       , colour = "Age Group"
       , x = "Screener threshold"
       , y = "FPR") +
  theme_minimal()+
  scale_x_continuous(breaks = 0:20)+
  scale_y_continuous(minor_breaks = 0:20/20)+
  theme(panel.grid.minor.x=element_blank())

# This keeps false alarms under 25% for everyone. 
# 7;5: 12.5%
# 8;5: 21.6%
# 9;5: 11.7%
# 10;5: 20.6%
# 10;6+: 14.6%

# Works much less well than the p16 version did.

### 8.2.3 "Saved" assessments ####

# How many kids would not meet the threshold for follow-up

table(spa7screeners$sub20p, spa7screeners$screener.tcflag)

# So across the dataset of 769 kids
# You would not follow up 512/614 (83%) kids who are not at risk (correct reject)
# You would correctly follow up on 155/155 (100%) kids who are at risk (hit)
# Of the kids you followed up, 102/257 (or 40%) would not be at risk (false alarm)
# Of the kids you did not follow up, 0/542 (or 0% would be misses)

# Overall you would reduce your follow up assessments by 512/769 or 66%

#### 8.2.3.1 Effect of Base Rate ####
# Note, this is in a general population where the proportion of at risk spellers
# is around 20%. If your proportion is much higher (say in a clinic where
# 2/3 of referred kids are actually at risk), the savings would be much 
# diminished.

baseRate = 0:100/100
corrRejRate = 512/614

pSaved = (1-baseRate)*corrRejRate

plot(baseRate, pSaved, type="l", xlim=c(0,1.1), ylim=c(0,1)
     , xlab = "% at risk in population"
     , ylab = "% reduction in # of full assessments"
     , main = "Efficiency of Screener")
points(baseRate[seq(1,101, length.out=11)],pSaved[seq(1,101, length.out=11)])
abline(v=.2, col="grey")
text(baseRate[seq(1,101, length.out=11)]+.05,pSaved[seq(1,101, length.out=11)]+.05, round(pSaved[seq(1,101, length.out=11)],2))
text(.2,0, "gen. pop'n", col="grey", pos=4)

## 8.3 Final list output ####
# Here are the items:
spa7agescreener.items =
  spa7tetrachoric.age %>% 
  apply(2
        , FUN= function(x) names(x[order(x)])[1:20])

write.csv(spa7agescreener.items, "spaScreener.p20.20Item.csv", row.names=F)

spa7items.long = spa7agescreener.items %>% data.frame() %>% 
  pivot_longer(cols=1:5, names_to = "age", values_to = "item") %>% 
  mutate(cor=NA)
for (i in 1:nrow(spa7items.long))
  spa7items.long$cor[i] = spa7tetrachoric.age[spa7items.long$item[i]
                                           , spa7items.long$age[i]]


(spa7items.long %>% 
  pivot_wider(id_cols = "item", names_from="age", values_from = "cor") %>% 
  data.frame() %>% 
    mutate(n = rowSums(!is.na(.))-1) %>% 
    arrange(desc(n), desc(cor7), desc(cor8), desc(cor9), desc(cor10), desc(cor11))
  ) %>% 
  write.csv("spaScreener.p20.20Item.ItemCors.csv", row.names=F)

# 9.0 Performance on other scales ####

# The above analysis only considers how well the screeners do at predicting
# overall SPA spelling. Here we consider how the screeners do when trying to
# identify kids with more specific spelling issues (irregs, regs, nonwords
# and even Morph) or on other spelling tests (WIAT)

### The focus will be on the Irregs, regs, and nonwords with the intention of
# perhaps expanding the number of those items in order to try and improve the
# performance on subtasks.

#### NOTE: I think in order to do this, the smarter approach might be to look
# at how many items would need to be included, to get perfect performance o
# each subtype, and then look at how the aggregate of those works

### See script SPA.ItemTypeAnalysis.R

## 9.1 SPA subscales ####
# The SPA is made up of three item-types. A child could be fine overall, but
# poor at one type (say nonword spelling). The screener items were not chosen
# with any consideration to item-type, though they do include a variety:

table(cut(which(order(spa7tetrachoric.age$cor7)<=20), breaks=c(0,29, 58, 87), labels = c("I", "R", "N")))
# 5 Irr, 4 Reg, 11 NW
table(cut(which(order(spa7tetrachoric.age$cor8)<=20), breaks=c(0,29, 58, 87), labels = c("I", "R", "N")))
# 9 Irr, 9 Reg, 2 NW
table(cut(which(order(spa7tetrachoric.age$cor9)<=20), breaks=c(0,29, 58, 87), labels = c("I", "R", "N")))
# 5 Irr, 10 Reg, 5 NW
table(cut(which(order(spa7tetrachoric.age$cor10)<=20), breaks=c(0,29, 58, 87), labels = c("I", "R", "N")))
# 3 Irr, 8 Reg, 9 NW
table(cut(which(order(spa7tetrachoric.age$cor11)<=20), breaks=c(0,29, 58, 87), labels = c("I", "R", "N")))
# 8 Irr, 8 Reg, 4 NW

# So for each subsale, we can identify the poorest readers using the same Quant 
# Regression technique, then apply the same ROC curve approach to see how well 
# the screeners do at catching kids who might be ok overall, but poor at one type

# This is just sensitivity (we don't care as much about specificity here, 
# because the screener will only be flagging them for follow up - their specific 
# deficits can be assessed more directly later)

### 9.1.1 Irregular Words ####

#### 9.1.1.1 QR for IDing "20th percentile" ####

# ggplot(spa7plus, aes(x=sp_alt, y=ability, col=sp_age))+geom_point()
spa7plus$sp_irr_logit = qlogis((spa7plus$sp_irr_alt+1)/31)
spa7.qr_irr = rq(sp_irr_logit~sp_age, spa7plus, tau=.20)

spa7plus$qr_irr.fit = 29*plogis(predict(spa7.qr_irr))

qrMod.irr = data.frame(sp_age=78:144) %>% 
  mutate(sp_irr_alt=29*plogis(predict(spa7.qr_irr, data.frame(sp_age = 78:144))))

ggplot(spa7plus, aes(x=sp_age, y=sp_irr_alt))+
  geom_point(col="darkgrey")+
  # geom_smooth()+
  geom_line(data=qrMod.irr, aes(colour=NULL), size=1, col="blue")+
  theme_minimal()+
  labs(x="Age in months", y="Score (of 29)"
       , title="20th Percentile Performance (Irregs)", subtitle="as a function of age")

# So, let's flag the poorest 20% of SPA spellers (age adjusted)

spa7plus = spa7plus %>% 
  mutate(sub20p.irr = ifelse(sp_irr_alt < qr_irr.fit, 1, 0))
# table(spa7plus$sub20p.irr)/nrow(spa7plus)

#### 9.1.1.2 Sensitivity ####

spa7screeners = spa7screeners %>% 
  merge(spa7plus %>% select(newID, sub20p.irr))

spa7detection.agescreeners.irr = spa7screeners %>% 
  group_by(ageCat, screener.age) %>% 
  summarise(detection.irr=sum(sub20p.irr)) %>% #slice(1:20) %>% data.frame()
  # ungroup(screener) %>% 
  mutate(detection.irr=cumsum(detection.irr)/sum(detection.irr)) %>% 
  ungroup() %>% 
  arrange(ageCat, screener.age)

ggplot(spa7detection.agescreeners.irr
       , aes(x=screener.age, y=detection.irr, col=ageCat))+
  geom_line()+
  labs(title = "Screener's rate of detecting Irregular SPA 20th percentile kids"
       , subtitle = "by ages and threshold score"
       , colour = "Age Group"
       , x = "Screener threshold"
       , y = "Detect sub 20th percentile") +
  theme_minimal()+
  scale_x_continuous(breaks = 0:20)+
  scale_y_continuous(minor_breaks = 0:20/20)+
  theme(panel.grid.minor.x=element_blank())

spa7detection.agescreeners.irr %>% 
  filter(detection.irr == 1 | (ageCat == "<= 7;5" & screener.age == 11) | 
  (ageCat == "<= 8;5" & screener.age == 9) |
  (ageCat == "<= 9;5" & screener.age == 12) |
  (ageCat == "<= 10;5" & screener.age == 14) | 
  (ageCat == "10;6+" & screener.age == 16)) %>% 
  data.frame %>% 
  group_by(ageCat) %>% slice(1:2)
# ok, so from this we have that the "full-scale thresholds" that did well,
# miss a lot of poor irreg readers.
# 7;5: threshold 11 catches only 65%, would need a threshold of 16 to get all poor 
#       Irreg spellers (5 irreg items)
# 8;5: threshold 9 catches 88%, would need a threshold of 12 to get 100% (9)
# 9;5: threshold 12 catches 80+%, would need 17 to get 100% (5)
# 10;5: threshold 14 catches 75+%, would need 18 to get 100% (3)
# 11;5: threshold 16 catches 95+%, would need 17 to get 100% (8)

# There is something of a link between the number of irreg items on the screener
# and how well it does, though not a perfect correlation.


### 9.1.2 Nonwords ####

#### 9.1.2.1 QR for IDing "20th percentile" ####

spa7plus$sp_nw_logit = qlogis((spa7plus$sp_nw_alt+1)/31)
(spa7.qr_nw = rq(sp_nw_logit~sp_age, spa7plus, tau=.20)) %>% summary

spa7plus$qr_nw.fit = 29*plogis(predict(spa7.qr_nw))

qrMod.nw = data.frame(sp_age=78:144) %>% 
  mutate(sp_nw_alt=29*plogis(predict(spa7.qr_nw, data.frame(sp_age = 78:144))))

ggplot(spa7plus, aes(x=sp_age, y=sp_nw_alt))+
  geom_point(col="darkgrey")+
  # geom_smooth()+
  geom_line(data=qrMod.nw, aes(colour=NULL), size=1, col="blue")+
  theme_minimal()+
  labs(x="Age in months", y="Score (of 29)"
       , title="20th Percentile Performance (Regs)", subtitle="as a function of age")

# So, let's flag the poorest 20% of SPA spellers (age adjusted)

spa7plus = spa7plus %>% 
  mutate(sub20p.nw = ifelse(sp_nw_alt < qr_nw.fit, 1, 0))
# table(spa7plus$sub20p.nw)/nrow(spa7plus)

#### 9.1.2.2 Sensitivity ####

spa7screeners = spa7screeners %>% 
  merge(spa7plus %>% select(newID, sub20p.nw))

spa7detection.agescreeners.nw = spa7screeners %>% 
  group_by(ageCat, screener.age) %>% 
  summarise(detection.nw=sum(sub20p.nw)) %>% #slice(1:20) %>% data.frame()
  # ungroup(screener) %>% 
  mutate(detection.nw=cumsum(detection.nw)/sum(detection.nw)) %>% 
  ungroup() %>% 
  arrange(ageCat, screener.age)

ggplot(spa7detection.agescreeners.nw
       , aes(x=screener.age, y=detection.nw, col=ageCat))+
  geom_line()+
  labs(title = "Screener's rate of detecting NW SPA 20th percentile kids"
       , subtitle = "by ages and threshold score"
       , colour = "Age Group"
       , x = "Screener threshold"
       , y = "Detect sub 20th percentile") +
  theme_minimal()+
  scale_x_continuous(breaks = 0:20)+
  scale_y_continuous(minor_breaks = 0:20/20)+
  theme(panel.grid.minor.x=element_blank())

spa7detection.agescreeners.nw %>% 
  filter(detection.nw == 1 | (ageCat == "<= 7;5" & screener.age == 11) | 
           (ageCat == "<= 8;5" & screener.age == 9) |
           (ageCat == "<= 9;5" & screener.age == 12) |
           (ageCat == "<= 10;5" & screener.age == 14) | 
           (ageCat == "10;6+" & screener.age == 16)) %>% 
  data.frame %>% 
  group_by(ageCat) %>% slice(1:2)
# NWs do not feature much in the screeners so it doesn't do very well here.
# the scale does pretty well here
# 7;5: threshold 11 catches 95%, would need a threshold of 14 to get all poor 
#       Irnw spellers (4 Regular items)
# 8;5: threshold 9 catches 75+%, would need a threshold of 19 to get 100% (9)
# 9;5: threshold 12 catches 75+%, would need 17 to get 100% (10)
# 10;5: threshold 14 catches 85+%, would need 17 to get 100% (8)
# 11;5: threshold 16 catches 75+%, need 20 to get 100% (8)

# There is something of a link between the number of Reg items on the screener
# and how well it does.

### 9.1.3 Regular Words ####

#### 9.1.3.1 QR for IDing "20th percentile" ####

spa7plus$sp_reg_logit = qlogis((spa7plus$sp_reg_alt+1)/31)
(spa7.qr_reg = rq(sp_reg_logit~sp_age, spa7plus, tau=.20)) %>% summary

spa7plus$qr_reg.fit = 29*plogis(predict(spa7.qr_reg))

qrMod.reg = data.frame(sp_age=78:144) %>% 
  mutate(sp_reg_alt=29*plogis(predict(spa7.qr_reg, data.frame(sp_age = 78:144))))

ggplot(spa7plus, aes(x=sp_age, y=sp_reg_alt))+
  geom_point(col="darkgrey")+
  # geom_smooth()+
  geom_line(data=qrMod.reg, aes(colour=NULL), size=1, col="blue")+
  theme_minimal()+
  labs(x="Age in months", y="Score (of 29)"
       , title="20th Percentile Performance (Regular)", subtitle="as a function of age")

# So, let's flag the poorest 20% of SPA spellers (age adjusted)

spa7plus = spa7plus %>% 
  mutate(sub20p.reg = ifelse(sp_reg_alt < qr_reg.fit, 1, 0))

#### 9.1.3.2 Sensitivity ####

spa7screeners = spa7screeners %>% 
  merge(spa7plus %>% select(newID, sub20p.reg))

spa7detection.agescreeners.reg = spa7screeners %>% 
  group_by(ageCat, screener.age) %>% 
  summarise(detection.reg=sum(sub20p.reg)) %>% #slice(1:20) %>% data.frame()
  # ungroup(screener) %>% 
  mutate(detection.reg=cumsum(detection.reg)/sum(detection.reg)) %>% 
  ungroup() %>% 
  arrange(ageCat, screener.age)

ggplot(spa7detection.agescreeners.reg
       , aes(x=screener.age, y=detection.reg, col=ageCat))+
  geom_line()+
  labs(title = "Screener's rate of detecting Regular SPA 20th percentile kids"
       , subtitle = "by ages and threshold score"
       , colour = "Age Group"
       , x = "Screener threshold"
       , y = "Detect sub 20th percentile") +
  theme_minimal()+
  scale_x_continuous(breaks = 0:20)+
  scale_y_continuous(minor_breaks = 0:20/20)+
  theme(panel.grid.minor.x=element_blank())

spa7detection.agescreeners.reg %>% 
  filter(detection.reg == 1 | (ageCat == "<= 7;5" & screener.age == 11) | 
           (ageCat == "<= 8;5" & screener.age == 9) |
           (ageCat == "<= 9;5" & screener.age == 12) |
           (ageCat == "<= 10;5" & screener.age == 14) | 
           (ageCat == "10;6+" & screener.age == 20)) %>% 
  data.frame %>% 
  group_by(ageCat) %>% slice(1:2) %>% 
  mutate(res = c("detection", "threshold")
         , val = ifelse(res=="detection", detection.reg, screener.age)) %>%
  select(ageCat, res, val) %>% 
  pivot_wider(id_cols=ageCat, names_from=res, values_from=val) %>% 
  mutate(value = sprintf("%.1f%% (%d)", detection*100, threshold))
# ok, so from this we have that the "full-scale thresholds" that did well,
# miss a lot of poor nw readers.
# 7;5: threshold 11 catches 90%, would need a threshold of 15 to get all poor 
#       Reg spellers (11 nw items)
# 8;5: threshold 9 catches 95%, would need a threshold of 15 to get 100% (2)
# 9;5: threshold 12 catches 100%, only need 11 to get 100% (9)
# 10;5: threshold 14 catches 95%, would need 16 to get 100% (5)
# 11;5: threshold 16 catches 100%, need 16 to get 100% (4)

# There is something of a link between the number of nw items on the screener
# and how well it does.

## 9.2 WIAT scale ####

### 9.2.1 import WIAT and Morph data ####
# pptDetails includes the scores for Morph and WIAT

spaMorphWiat = merge(spa7screeners
                     , pptDetails %>% 
                       select(newID, sp_age
                              , `wiat.sp.ss mean 100`
                              , Morph.tot
                              , `Morph_spell test_raw PS`
                              , `Morph_spell test_raw RW`) %>% 
                       rename(WIAT.ss = `wiat.sp.ss mean 100`
                              , Morph.ps = `Morph_spell test_raw PS`
                              , Morph.rw = `Morph_spell test_raw RW`) %>% 
                       mutate(sp_age = as.numeric(sp_age))
                     , by = "newID", all.x=T)

### 9.2.3 WIAT Sensitivity ####
# The WIAT is already age standardized, so we don't need to do anything to
# identify the 20th percentile. Instead we can just use 1 SD below the mean,
# which in this case is any score of 85 or less.

spaMorphWiat = spaMorphWiat %>% 
  mutate(WIAT.sub20p = ifelse(WIAT.ss <=85, 1, 0))

WIATdetection.agescreeners = spaMorphWiat %>% 
  group_by(ageCat, screener.age) %>% 
  summarise(detection=sum(WIAT.sub20p, na.rm=T)) %>% #slice(1:20) %>% data.frame()
  # ungroup(screener) %>% 
  mutate(detection=cumsum(detection)/sum(detection)) %>% 
  ungroup() %>% 
  arrange(ageCat, screener.age)

ggplot(WIATdetection.agescreeners
       , aes(x=screener.age, y=detection, col=ageCat))+
  geom_line()+
  labs(title = "Screener's rate of detecting WIAT -1SD kids"
       , subtitle = "by ages and threshold score"
       , colour = "Age Group"
       , x = "Screener threshold"
       , y = "Detect -1SD") +
  theme_minimal()+
  scale_x_continuous(breaks = 0:20)+
  scale_y_continuous(minor_breaks = 0:20/20)+
  theme(panel.grid.minor.x=element_blank())

WIATdetection.agescreeners %>% 
  filter(detection == 1 | (ageCat == "<= 7;5" & screener.age == 11) | 
           (ageCat == "<= 8;5" & screener.age == 9) |
           (ageCat == "<= 9;5" & screener.age == 12) |
           (ageCat == "<= 10;5" & screener.age == 14) | 
           (ageCat == "10;6+" & screener.age == 16)) %>% 
  data.frame %>% 
  group_by(ageCat) %>% slice(1:2)
# ouch. Very poor, other than for kids in the 10;6 - 12;0 range.

## 9.3 Morph ####
# Morph does not have a standardized scoring yet, so we can only use this 
# sample. As above, I'll rely on Quantile Regression to identify age based
# poor readers.

# Note that the sample is actually quite small for this - N = 400

### 9.3.1 Morph Total ####

#### 9.3.1.1 QR for IDing "20th percentile" ####

## TODO: I don't know for sure how many items there are on the morph
# but kids do as well as 27 on the words, and 26 on the pseudowords
# so it must be at least 53.

spaMorphWiat$morph_logit = qlogis((spaMorphWiat$Morph.tot+1)/55)
(morph.qr_tot = rq(morph_logit~sp_age, spaMorphWiat, tau=.20)) %>% summary

spaMorphWiat$morph.fit = NA
spaMorphWiat[rownames(morph.qr_tot$x),"morph.fit"] = 53*plogis(predict(morph.qr_tot))

# Age-based threshold table
qrMod.morph.tot = data.frame(sp_age=78:144) %>% 
  mutate(Morph.tot=53*plogis(predict(morph.qr_tot, data.frame(sp_age = 78:144))))

ggplot(spaMorphWiat, aes(x=sp_age, y=Morph.tot))+
  geom_point(col="darkgrey")+
  # geom_smooth()+
  geom_line(data=qrMod.morph.tot, aes(colour=NULL), size=1, col="blue")+
  theme_minimal()+
  labs(x="Age in months", y="Score (of 29)"
       , title="20th Percentile Performance (Morph Total)", subtitle="as a function of age")

# So, let's flag the poorest 20% of SPA spellers (age adjusted)

spaMorphWiat = spaMorphWiat %>% 
  mutate(sub20p.morph.tot = ifelse(Morph.tot < morph.fit, 1, 0))

ggplot(spaMorphWiat, aes(x=sp_age, y=Morph.tot, col=sub20p.morph.tot))+
  geom_point()+
  # geom_smooth()+
  geom_line(data=qrMod.morph.tot, aes(colour=NULL), size=1, col="blue")+
  theme_minimal()+
  labs(x="Age in months", y="Score (of 53)"
       , title="20th Percentile Performance (Morph Total)", subtitle="as a function of age")

#### 9.3.1.2 Sensitivity ####
# Problem - some scores don't show up in the Morph data.
morphScreeners = spa7screeners %>% 
  filter(newID %in% spaMorphWiat$newID[!is.na(spaMorphWiat$Morph.tot)]) %>% 
  select(newID, ageCat, starts_with("screener")) %>% 
  merge(spaMorphWiat %>% select(newID, sub20p.morph.tot))

Morphdetection.agescreeners = expand.grid(
  ageCat = c('<= 7;5', '<= 8;5', '<= 9;5', '<= 10;5', '10;6+')
  , screener.age = 0:20
) %>% data.frame()

Morphdetection.agescreeners.morph.tot = morphScreeners %>% 
  group_by(ageCat, screener.age) %>% 
  summarise(detection.morph=sum(sub20p.morph.tot)) %>% #slice(1:20) %>% data.frame()
  # ungroup(screener) %>% 
  mutate(detection.morph=cumsum(detection.morph)/sum(detection.morph, na.rm=T)) %>% 
  merge(Morphdetection.agescreeners, all.y=T) %>% 
  group_by(ageCat) %>% 
  mutate(
    detection.morph = na.locf(detection.morph, na.rm=FALSE)
    , detection.morph = ifelse(is.na(detection.morph), 0, detection.morph)
  ) %>% 
  ungroup() %>% 
  arrange(ageCat, screener.age)

ggplot(Morphdetection.agescreeners.morph.tot
       , aes(x=screener.age, y=detection.morph, col=ageCat))+
  geom_line()+
  labs(title = "Screener's rate of detecting Morph Total 20th percentile kids"
       , subtitle = "by ages and threshold score"
       , colour = "Age Group"
       , x = "Screener threshold"
       , y = "Detect sub 20th percentile") +
  theme_minimal()+
  scale_x_continuous(breaks = 0:20)+
  scale_y_continuous(minor_breaks = 0:20/20)+
  theme(panel.grid.minor.x=element_blank())

Morphdetection.agescreeners.morph.tot %>%
  arrange(detection.morph, screener.age) %>% 
  filter(detection.morph == 1 | (ageCat == "<= 7;5" & screener.age == 11) | 
           (ageCat == "<= 8;5" & screener.age == 9) |
           (ageCat == "<= 9;5" & screener.age == 12) |
           (ageCat == "<= 10;5" & screener.age == 14) | 
           (ageCat == "10;6+" & screener.age == 16)) %>% 
  data.frame %>% 
  group_by(ageCat) %>% slice(1:2) %>% 
  mutate(res = c("detection", "threshold")
         , val = ifelse(res=="detection", detection.morph, screener.age)) %>%
  select(ageCat, res, val) %>% 
  pivot_wider(id_cols=ageCat, names_from=res, values_from=val) %>% 
  mutate(value = sprintf("%.1f%% (%d)", detection*100, threshold))
# ok, so from this we have that the "full-scale thresholds" that did well,
# miss a lot of poor nw readers.
# 7;5: threshold 11 catches 57%, would need a threshold of 17 to get all poor 
#       morph spellers (11 nw items)
# 8;5: threshold 9 catches 92%, would need a threshold of 19 to get 100% (2)
# 9;5: threshold 12 catches 57%, would need 20 to get 100% (9)
# 10;5: threshold 14 catches 75%, would need 20 to get 100% (5)
# 11;5: threshold 16 catches 46%, need 20 to get 100% (4)

### 9.3.2 Morph Words ####

## TODO: I don't know for sure how many items there are on the morph
# but kids do as well as 27 on the words.

#### 9.3.2.1 QR for IDing "20th percentile" ####

spaMorphWiat$morph.rw_logit = qlogis((spaMorphWiat$Morph.rw+1)/29)
(morph.qr_rw = rq(morph.rw_logit~sp_age, spaMorphWiat, tau=.20)) %>% summary

spaMorphWiat$morph.rw.fit = NA
spaMorphWiat[rownames(morph.qr_rw$x),"morph.rw.fit"] = 27*plogis(predict(morph.qr_rw))

# Age-based threshold table
qrMod.morph.rw = data.frame(sp_age=78:144) %>% 
  mutate(Morph.rw=27*plogis(predict(morph.qr_rw, data.frame(sp_age = 78:144))))

ggplot(spaMorphWiat, aes(x=sp_age, y=Morph.rw))+
  geom_point(col="darkgrey")+
  # geom_smooth()+
  geom_line(data=qrMod.morph.rw, aes(colour=NULL), size=1, col="blue")+
  theme_minimal()+
  labs(x="Age in months", y="Score (of 27)"
       , title="20th Percentile Performance (Morph Words)", subtitle="as a function of age")

# So, let's flag the poorest 20% of SPA spellers (age adjusted)

spaMorphWiat = spaMorphWiat %>% 
  mutate(sub20p.morph.rw = ifelse(Morph.rw < morph.rw.fit, 1, 0))

ggplot(spaMorphWiat, aes(x=sp_age, y=Morph.rw, col=sub20p.morph.rw))+
  geom_point()+
  # geom_smooth()+
  geom_line(data=qrMod.morph.rw, aes(colour=NULL), size=1, col="blue")+
  theme_minimal()+
  labs(x="Age in months", y="Score (of 29)"
       , title="20th Percentile Performance (Morph Words)", subtitle="as a function of age")

#### 9.3.2.2 Sensitivity ####

morphScreeners = morphScreeners %>% 
  filter(newID %in% spaMorphWiat$newID[!is.na(spaMorphWiat$Morph.rw)]) %>% 
  #select(newID, ageCat, starts_with("screener")) %>% 
  merge(spaMorphWiat %>% select(newID, sub20p.morph.rw))

# Morphdetection.agescreeners.morph.rw = morphScreeners %>% 
#   group_by(ageCat, screener.age) %>% 
#   summarise(detection.morph.rw=sum(sub20p.morph.rw)) %>% #slice(1:20) %>% data.frame()
#   # ungroup(screener) %>% 
#   mutate(detection.morph=cumsum(detection.morph.rw)/sum(detection.morph.rw, na.rm=T)) %>% 
#   ungroup() %>% 
#   arrange(ageCat, screener.age)

Morphdetection.agescreeners.morph.rw = morphScreeners %>% 
  group_by(ageCat, screener.age) %>% 
  summarise(detection.morph.rw=sum(sub20p.morph.rw)) %>% #slice(1:20) %>% data.frame()
  # ungroup(screener) %>% 
  mutate(detection.morph.rw=cumsum(detection.morph.rw)/sum(detection.morph.rw, na.rm=T)) %>% 
  merge(Morphdetection.agescreeners, all.y=T) %>% 
  group_by(ageCat) %>% 
  mutate(
    detection.morph.rw = na.locf(detection.morph.rw, na.rm=FALSE)
    , detection.morph.rw = ifelse(is.na(detection.morph.rw), 0, detection.morph.rw)
  ) %>% 
  ungroup() %>% 
  arrange(ageCat, screener.age)

ggplot(Morphdetection.agescreeners.morph.rw
       , aes(x=screener.age, y=detection.morph.rw, col=ageCat))+
  geom_line()+
  labs(title = "Screener's rate of detecting Morph Words 20th percentile kids"
       , subtitle = "by ages and threshold score"
       , colour = "Age Group"
       , x = "Screener threshold"
       , y = "Detect sub 20th percentile") +
  theme_minimal()+
  scale_x_continuous(breaks = 0:20)+
  scale_y_continuous(minor_breaks = 0:20/20)+
  theme(panel.grid.minor.x=element_blank())

Morphdetection.agescreeners.morph.rw %>% 
  filter(detection.morph.rw == 1 | (ageCat == "<= 7;5" & screener.age == 11) | 
           (ageCat == "<= 8;5" & screener.age == 9) |
           (ageCat == "<= 9;5" & screener.age == 12) |
           (ageCat == "<= 10;5" & screener.age == 14) | 
           (ageCat == "10;6+" & screener.age == 16)) %>% 
  data.frame %>% 
  group_by(ageCat) %>% slice(1:2) %>% 
  mutate(res = c("detection", "threshold")
         , val = ifelse(res=="detection", detection.morph.rw, screener.age)) %>%
  select(ageCat, res, val) %>% 
  pivot_wider(id_cols=ageCat, names_from=res, values_from=val) %>% 
  mutate(value = sprintf("%.1f%% (%d)", detection*100, threshold))
# ok, so from this we have that the "full-scale thresholds" that did well,
# miss a lot of poor nw readers.
# 7;5: threshold 11 catches 58%, would need a threshold of 17 to get all poor 
#       morph spellers (11 nw items)
# 8;5: threshold 9 catches 86%, would need a threshold of 19 to get 100% (2)
# 9;5: threshold 12 catches 50%, would need 20 to get 100% (9)
# 10;5: threshold 14 catches 67%, would need 20 to get 100% (5)
# 11;5: threshold 16 catches 70%, need 20 to get 100% (4)

### 9.3.3 Morph Pseudowords ####

## TODO: I don't know for sure how many items there are on the morph
# but kids do as well as 26 on the pseudowords.

#### 9.3.3.1 QR for IDing "20th percentile" ####

spaMorphWiat$morph.ps_logit = qlogis((spaMorphWiat$Morph.ps+1)/28)
(morph.qr_ps = rq(morph.ps_logit~sp_age, spaMorphWiat, tau=.20)) %>% summary

spaMorphWiat$morph.ps.fit = NA
spaMorphWiat[rownames(morph.qr_ps$x),"morph.ps.fit"] = 26*plogis(predict(morph.qr_ps))

# Age-based threshold table
qrMod.morph.ps = data.frame(sp_age=78:144) %>% 
  mutate(Morph.ps=26*plogis(predict(morph.qr_ps, data.frame(sp_age = 78:144))))

ggplot(spaMorphWiat, aes(x=sp_age, y=Morph.ps))+
  geom_point(col="darkgrey")+
  # geom_smooth()+
  geom_line(data=qrMod.morph.ps, aes(colour=NULL), size=1, col="blue")+
  theme_minimal()+
  labs(x="Age in months", y="Score (of 27)"
       , title="20th Percentile Performance (Morph Words)", subtitle="as a function of age")

# So, let's flag the poorest 20% of SPA spellers (age adjusted)

spaMorphWiat = spaMorphWiat %>% 
  mutate(sub20p.morph.ps = ifelse(Morph.ps < morph.ps.fit, 1, 0))

ggplot(spaMorphWiat, aes(x=sp_age, y=Morph.ps, col=sub20p.morph.ps))+
  geom_point()+
  # geom_smooth()+
  geom_line(data=qrMod.morph.ps, aes(colour=NULL), size=1, col="blue")+
  theme_minimal()+
  labs(x="Age in months", y="Score (of 29)"
       , title="20th Percentile Performance (Morph Pseudowords)", subtitle="as a function of age")

#### 9.3.3.2 Sensitivity ####

morphScreeners = morphScreeners %>% 
  filter(newID %in% spaMorphWiat$newID[!is.na(spaMorphWiat$Morph.ps)]) %>% 
  #select(newID, ageCat, starts_with("screener")) %>% 
  merge(spaMorphWiat %>% select(newID, sub20p.morph.ps))

# Morphdetection.agescreeners.morph.ps = morphScreeners %>% 
#   group_by(ageCat, screener.age) %>% 
#   summarise(detection.morph.ps=sum(sub20p.morph.ps)) %>% #slice(1:20) %>% data.frame()
#   # ungroup(screener) %>% 
#   mutate(detection.morph=cumsum(detection.morph.ps)/sum(detection.morph.ps, na.rm=T)) %>% 
#   ungroup() %>% 
#   arrange(ageCat, screener.age)

Morphdetection.agescreeners.morph.ps = morphScreeners %>% 
  group_by(ageCat, screener.age) %>% 
  summarise(detection.morph.ps=sum(sub20p.morph.ps)) %>% #slice(1:20) %>% data.frame()
  # ungroup(screener) %>% 
  mutate(detection.morph.ps=cumsum(detection.morph.ps)/sum(detection.morph.ps, na.rm=T)) %>% 
  merge(Morphdetection.agescreeners, all.y=T) %>% 
  group_by(ageCat) %>% 
  mutate(
    detection.morph.ps = na.locf(detection.morph.ps, na.rm=FALSE)
    , detection.morph.ps = ifelse(is.na(detection.morph.ps), 0, detection.morph.ps)
  ) %>% 
  ungroup() %>% 
  arrange(ageCat, screener.age)

ggplot(Morphdetection.agescreeners.morph.ps
       , aes(x=screener.age, y=detection.morph.ps, col=ageCat))+
  geom_line()+
  labs(title = "Screener's rate of detecting Morph Pseudowords 20th percentile kids"
       , subtitle = "by ages and threshold score"
       , colour = "Age Group"
       , x = "Screener threshold"
       , y = "Detect sub 20th percentile") +
  theme_minimal()+
  scale_x_continuous(breaks = 0:20)+
  scale_y_continuous(minor_breaks = 0:20/20)+
  theme(panel.grid.minor.x=element_blank())

Morphdetection.agescreeners.morph.ps %>% 
  filter(detection.morph.ps == 1 | (ageCat == "<= 7;5" & screener.age == 11) | 
           (ageCat == "<= 8;5" & screener.age == 9) |
           (ageCat == "<= 9;5" & screener.age == 12) |
           (ageCat == "<= 10;5" & screener.age == 14) | 
           (ageCat == "10;6+" & screener.age == 16)) %>% 
  data.frame %>% 
  group_by(ageCat) %>% slice(1:2) %>% 
  mutate(res = c("detection", "threshold")
         , val = ifelse(res=="detection", detection.morph.ps, screener.age)) %>%
  select(ageCat, res, val) %>% 
  pivot_wider(id_cols=ageCat, names_from=res, values_from=val) %>% 
  mutate(value = sprintf("%.1f%% (%d)", detection*100, threshold))

# ok, so from this we have that the "full-scale thresholds" that did well,
# miss a lot of poor nw readers.
# 7;5: threshold 11 catches 60%, would need a threshold of 17 to get all poor 
#       morph spellers (11 nw items)
# 8;5: threshold 9 catches 80%, would need a threshold of 19 to get 100% (2)
# 9;5: threshold 12 catches 31%, would need 20 to get 100% (9)
# 10;5: threshold 14 catches 58%, would need 20 to get 100% (5)
# 11;5: threshold 16 catches 29%, need 20 to get 100% (4)





## 9.4 Morph (assumes 29 items per scale) ####
# Morph does not have a standardized scoring yet, so we can only use this 
# sample. As above, I'll rely on Quantile Regression to identify age based
# poor readers.

# Note that the sample is actually quite small for this - N = 400

### 9.4.1 Morph Total ####

#### 9.4.1.1 QR for IDing "20th percentile" ####

## TODO: I don't know for sure how many items there are on the morph
# but kids do as well as 27 on the words, and 26 on the pseudowords
# so it must be at least 53.

# Per Saskia's email, I'll do this with 29 and 29 for 58 until we get
# a confirmed value.

spaMorphWiat$morph_logit = qlogis((spaMorphWiat$Morph.tot+1)/60)
(morph.qr_tot = rq(morph_logit~sp_age, spaMorphWiat, tau=.20)) %>% summary

spaMorphWiat$morph.fit = NA
spaMorphWiat[rownames(morph.qr_tot$x),"morph.fit"] = 58*plogis(predict(morph.qr_tot))

# Age-based threshold table
qrMod.morph.tot = data.frame(sp_age=78:144) %>% 
  mutate(Morph.tot=58*plogis(predict(morph.qr_tot, data.frame(sp_age = 78:144))))

ggplot(spaMorphWiat, aes(x=sp_age, y=Morph.tot))+
  geom_point(col="darkgrey")+
  # geom_smooth()+
  geom_line(data=qrMod.morph.tot, aes(colour=NULL), size=1, col="blue")+
  theme_minimal()+
  labs(x="Age in months", y="Score (of 58)"
       , title="20th Percentile Performance (Morph Total)", subtitle="as a function of age")

# So, let's flag the poorest 20% of SPA spellers (age adjusted)

spaMorphWiat = spaMorphWiat %>% 
  mutate(sub20p.morph.tot = ifelse(Morph.tot < morph.fit, 1, 0))

ggplot(spaMorphWiat, aes(x=sp_age, y=Morph.tot, col=sub20p.morph.tot))+
  geom_point()+
  # geom_smooth()+
  geom_line(data=qrMod.morph.tot, aes(colour=NULL), size=1, col="blue")+
  theme_minimal()+
  labs(x="Age in months", y="Score (of 58)"
       , title="20th Percentile Performance (Morph Total)", subtitle="as a function of age")

#### 9.4.1.2 Sensitivity ####
# Problem - some scores don't show up in the Morph data.
morphScreeners = spa7screeners %>% 
  filter(newID %in% spaMorphWiat$newID[!is.na(spaMorphWiat$Morph.tot)]) %>% 
  select(newID, ageCat, starts_with("screener")) %>% 
  merge(spaMorphWiat %>% select(newID, sub20p.morph.tot))

Morphdetection.agescreeners = expand.grid(
  ageCat = c('<= 7;5', '<= 8;5', '<= 9;5', '<= 10;5', '10;6+')
  , screener.age = 0:20
) %>% data.frame()

Morphdetection.agescreeners.morph.tot = morphScreeners %>% 
  group_by(ageCat, screener.age) %>% 
  summarise(detection.morph=sum(sub20p.morph.tot)) %>% #slice(1:20) %>% data.frame()
  # ungroup(screener) %>% 
  mutate(detection.morph=cumsum(detection.morph)/sum(detection.morph, na.rm=T)) %>% 
  merge(Morphdetection.agescreeners, all.y=T) %>% 
  group_by(ageCat) %>% 
  mutate(
    detection.morph = na.locf(detection.morph, na.rm=FALSE)
    , detection.morph = ifelse(is.na(detection.morph), 0, detection.morph)
  ) %>% 
  ungroup() %>% 
  arrange(ageCat, screener.age)

ggplot(Morphdetection.agescreeners.morph.tot
       , aes(x=screener.age, y=detection.morph, col=ageCat))+
  geom_line()+
  labs(title = "Screener's rate of detecting Morph Total 20th percentile kids"
       , subtitle = "by ages and threshold score"
       , colour = "Age Group"
       , x = "Screener threshold"
       , y = "Detect sub 20th percentile") +
  theme_minimal()+
  scale_x_continuous(breaks = 0:20)+
  scale_y_continuous(minor_breaks = 0:20/20)+
  theme(panel.grid.minor.x=element_blank())

Morphdetection.agescreeners.morph.tot %>%
  arrange(detection.morph, screener.age) %>% 
  filter(detection.morph == 1 | (ageCat == "<= 7;5" & screener.age == 11) | 
           (ageCat == "<= 8;5" & screener.age == 9) |
           (ageCat == "<= 9;5" & screener.age == 12) |
           (ageCat == "<= 10;5" & screener.age == 14) | 
           (ageCat == "10;6+" & screener.age == 16)) %>% 
  data.frame %>% 
  group_by(ageCat) %>% slice(1:2) %>% 
  mutate(res = c("detection", "threshold")
         , val = ifelse(res=="detection", detection.morph, screener.age)) %>%
  select(ageCat, res, val) %>% 
  pivot_wider(id_cols=ageCat, names_from=res, values_from=val) %>% 
  mutate(value = sprintf("%.1f%% (%d)", detection*100, threshold))


### 9.4.2 Morph Words ####

## TODO: I don't know for sure how many items there are on the morph
# but kids do as well as 27 on the words. Assume 29 per Saskia.

#### 9.4.2.1 QR for IDing "20th percentile" ####

spaMorphWiat$morph.rw_logit = qlogis((spaMorphWiat$Morph.rw+1)/31)
(morph.qr_rw = rq(morph.rw_logit~sp_age, spaMorphWiat, tau=.20)) %>% summary

spaMorphWiat$morph.rw.fit = NA
spaMorphWiat[rownames(morph.qr_rw$x),"morph.rw.fit"] = 29*plogis(predict(morph.qr_rw))

# Age-based threshold table
qrMod.morph.rw = data.frame(sp_age=78:144) %>% 
  mutate(Morph.rw=29*plogis(predict(morph.qr_rw, data.frame(sp_age = 78:144))))

ggplot(spaMorphWiat, aes(x=sp_age, y=Morph.rw))+
  geom_point(col="darkgrey")+
  # geom_smooth()+
  geom_line(data=qrMod.morph.rw, aes(colour=NULL), size=1, col="blue")+
  theme_minimal()+
  labs(x="Age in months", y="Score (of 29)"
       , title="20th Percentile Performance (Morph Words)", subtitle="as a function of age")

# So, let's flag the poorest 20% of SPA spellers (age adjusted)

spaMorphWiat = spaMorphWiat %>% 
  mutate(sub20p.morph.rw = ifelse(Morph.rw < morph.rw.fit, 1, 0))

ggplot(spaMorphWiat, aes(x=sp_age, y=Morph.rw, col=sub20p.morph.rw))+
  geom_point()+
  # geom_smooth()+
  geom_line(data=qrMod.morph.rw, aes(colour=NULL), size=1, col="blue")+
  theme_minimal()+
  labs(x="Age in months", y="Score (of 29)"
       , title="20th Percentile Performance (Morph Words)", subtitle="as a function of age")

#### 9.4.2.2 Sensitivity ####

morphScreeners = morphScreeners %>% 
  filter(newID %in% spaMorphWiat$newID[!is.na(spaMorphWiat$Morph.rw)]) %>% 
  #select(newID, ageCat, starts_with("screener")) %>% 
  merge(spaMorphWiat %>% select(newID, sub20p.morph.rw))

# Morphdetection.agescreeners.morph.rw = morphScreeners %>% 
#   group_by(ageCat, screener.age) %>% 
#   summarise(detection.morph.rw=sum(sub20p.morph.rw)) %>% #slice(1:20) %>% data.frame()
#   # ungroup(screener) %>% 
#   mutate(detection.morph=cumsum(detection.morph.rw)/sum(detection.morph.rw, na.rm=T)) %>% 
#   ungroup() %>% 
#   arrange(ageCat, screener.age)

Morphdetection.agescreeners.morph.rw = morphScreeners %>% 
  group_by(ageCat, screener.age) %>% 
  summarise(detection.morph.rw=sum(sub20p.morph.rw)) %>% #slice(1:20) %>% data.frame()
  # ungroup(screener) %>% 
  mutate(detection.morph.rw=cumsum(detection.morph.rw)/sum(detection.morph.rw, na.rm=T)) %>% 
  merge(Morphdetection.agescreeners, all.y=T) %>% 
  group_by(ageCat) %>% 
  mutate(
    detection.morph.rw = na.locf(detection.morph.rw, na.rm=FALSE)
    , detection.morph.rw = ifelse(is.na(detection.morph.rw), 0, detection.morph.rw)
  ) %>% 
  ungroup() %>% 
  arrange(ageCat, screener.age)

ggplot(Morphdetection.agescreeners.morph.rw
       , aes(x=screener.age, y=detection.morph.rw, col=ageCat))+
  geom_line()+
  labs(title = "Screener's rate of detecting Morph Words 20th percentile kids"
       , subtitle = "by ages and threshold score"
       , colour = "Age Group"
       , x = "Screener threshold"
       , y = "Detect sub 20th percentile") +
  theme_minimal()+
  scale_x_continuous(breaks = 0:20)+
  scale_y_continuous(minor_breaks = 0:20/20)+
  theme(panel.grid.minor.x=element_blank())

Morphdetection.agescreeners.morph.rw %>% 
  filter(detection.morph.rw == 1 | (ageCat == "<= 7;5" & screener.age == 11) | 
           (ageCat == "<= 8;5" & screener.age == 9) |
           (ageCat == "<= 9;5" & screener.age == 12) |
           (ageCat == "<= 10;5" & screener.age == 14) | 
           (ageCat == "10;6+" & screener.age == 16)) %>% 
  data.frame %>% 
  group_by(ageCat) %>% slice(1:2) %>% 
  mutate(res = c("detection", "threshold")
         , val = ifelse(res=="detection", detection.morph.rw, screener.age)) %>%
  select(ageCat, res, val) %>% 
  pivot_wider(id_cols=ageCat, names_from=res, values_from=val) %>% 
  mutate(value = sprintf("%.1f%% (%d)", detection*100, threshold))

### 9.4.3 Morph Pseudowords ####

## TODO: I don't know for sure how many items there are on the morph
# but kids do as well as 26 on the pseudowords. Assume 30 here.

#### 9.4.3.1 QR for IDing "20th percentile" ####

spaMorphWiat$morph.ps_logit = qlogis((spaMorphWiat$Morph.ps+1)/31)
(morph.qr_ps = rq(morph.ps_logit~sp_age, spaMorphWiat, tau=.20)) %>% summary

spaMorphWiat$morph.ps.fit = NA
spaMorphWiat[rownames(morph.qr_ps$x),"morph.ps.fit"] = 29*plogis(predict(morph.qr_ps))

# Age-based threshold table
qrMod.morph.ps = data.frame(sp_age=78:144) %>% 
  mutate(Morph.ps=29*plogis(predict(morph.qr_ps, data.frame(sp_age = 78:144))))

ggplot(spaMorphWiat, aes(x=sp_age, y=Morph.ps))+
  geom_point(col="darkgrey")+
  # geom_smooth()+
  geom_line(data=qrMod.morph.ps, aes(colour=NULL), size=1, col="blue")+
  theme_minimal()+
  labs(x="Age in months", y="Score (of 29)"
       , title="20th Percentile Performance (Morph Words)", subtitle="as a function of age")

# So, let's flag the poorest 20% of SPA spellers (age adjusted)

spaMorphWiat = spaMorphWiat %>% 
  mutate(sub20p.morph.ps = ifelse(Morph.ps < morph.ps.fit, 1, 0))

ggplot(spaMorphWiat, aes(x=sp_age, y=Morph.ps, col=sub20p.morph.ps))+
  geom_point()+
  # geom_smooth()+
  geom_line(data=qrMod.morph.ps, aes(colour=NULL), size=1, col="blue")+
  theme_minimal()+
  labs(x="Age in months", y="Score (of 29)"
       , title="20th Percentile Performance (Morph Pseudowords)", subtitle="as a function of age")

#### 9.4.3.2 Sensitivity ####

morphScreeners = morphScreeners %>% 
  filter(newID %in% spaMorphWiat$newID[!is.na(spaMorphWiat$Morph.ps)]) %>% 
  #select(newID, ageCat, starts_with("screener")) %>% 
  merge(spaMorphWiat %>% select(newID, sub20p.morph.ps))

# Morphdetection.agescreeners.morph.ps = morphScreeners %>% 
#   group_by(ageCat, screener.age) %>% 
#   summarise(detection.morph.ps=sum(sub20p.morph.ps)) %>% #slice(1:20) %>% data.frame()
#   # ungroup(screener) %>% 
#   mutate(detection.morph=cumsum(detection.morph.ps)/sum(detection.morph.ps, na.rm=T)) %>% 
#   ungroup() %>% 
#   arrange(ageCat, screener.age)

Morphdetection.agescreeners.morph.ps = morphScreeners %>% 
  group_by(ageCat, screener.age) %>% 
  summarise(detection.morph.ps=sum(sub20p.morph.ps)) %>% #slice(1:20) %>% data.frame()
  # ungroup(screener) %>% 
  mutate(detection.morph.ps=cumsum(detection.morph.ps)/sum(detection.morph.ps, na.rm=T)) %>% 
  merge(Morphdetection.agescreeners, all.y=T) %>% 
  group_by(ageCat) %>% 
  mutate(
    detection.morph.ps = na.locf(detection.morph.ps, na.rm=FALSE)
    , detection.morph.ps = ifelse(is.na(detection.morph.ps), 0, detection.morph.ps)
  ) %>% 
  ungroup() %>% 
  arrange(ageCat, screener.age)

ggplot(Morphdetection.agescreeners.morph.ps
       , aes(x=screener.age, y=detection.morph.ps, col=ageCat))+
  geom_line()+
  labs(title = "Screener's rate of detecting Morph Pseudowords 20th percentile kids"
       , subtitle = "by ages and threshold score"
       , colour = "Age Group"
       , x = "Screener threshold"
       , y = "Detect sub 20th percentile") +
  theme_minimal()+
  scale_x_continuous(breaks = 0:20)+
  scale_y_continuous(minor_breaks = 0:20/20)+
  theme(panel.grid.minor.x=element_blank())

Morphdetection.agescreeners.morph.ps %>% 
  filter(detection.morph.ps == 1 | (ageCat == "<= 7;5" & screener.age == 11) | 
           (ageCat == "<= 8;5" & screener.age == 9) |
           (ageCat == "<= 9;5" & screener.age == 12) |
           (ageCat == "<= 10;5" & screener.age == 14) | 
           (ageCat == "10;6+" & screener.age == 16)) %>% 
  data.frame %>% 
  group_by(ageCat) %>% slice(1:2) %>% 
  mutate(res = c("detection", "threshold")
         , val = ifelse(res=="detection", detection.morph.ps, screener.age)) %>%
  select(ageCat, res, val) %>% 
  pivot_wider(id_cols=ageCat, names_from=res, values_from=val) %>% 
  mutate(value = sprintf("%.1f%% (%d)", detection*100, threshold))

# ok, so from this we have that the "full-scale thresholds" that did well,
# miss a lot of poor nw readers.
# 7;5: threshold 11 catches 60%, would need a threshold of 17 to get all poor 
#       morph spellers (11 nw items)
# 8;5: threshold 9 catches 80%, would need a threshold of 19 to get 100% (2)
# 9;5: threshold 12 catches 31%, would need 20 to get 100% (9)
# 10;5: threshold 14 catches 58%, would need 20 to get 100% (5)
# 11;5: threshold 16 catches 29%, need 20 to get 100% (4)

### 9.4.4 Conclusion ####
# The number of items does seem to matter, at least a bit.
# values assuming 27+26=53 items differ somewhat from results assuming
# 30+30=60 items.

# 10 Correlations among subscales ####

# If we're going to use the screener to assess performance on various subscales
# it is worth having a rough sense of how well each test/subscale "predicts"
# the other subscales. For example, if the WIAT and the SPA total are not
# highly correlated, then we would not expect a screener derived from the SPA to
# do well on the WIAT or vice versa.

## 10.1 raw scores ####
# For this to work, we need some kind of age-adjusted scores for the
# (sub)scales that don't have standardized norms. Right now, this can only be
# based on the sample.

# I'll have to give some thought to the best method here, but I think maybe
# the standardized residuals from a model regressing age out of the scores.

## 10.2 tetrachoric "20th percentile" ####

# it's also worth considering the degree to which being a poor performer one one
# (sub)scale associates with being a poor performer on other (sub)scales. This
# one is easier, since the QR analysis creates an age-adjusted categorization.

ageAdj20thP = spaMorphWiat %>% 
  select(newID, sp_age, contains("sub20p")) %>% 
  rename(
    Age = sp_age
    , SPA.total = sub20p, SPA.reg = sub20p.reg, SPA.irr = sub20p.irr
    , SPA.nw = sub20p.nw
    , Morph.total = sub20p.morph.tot, Morph.rw = sub20p.morph.rw
    , Morph.ps = sub20p.morph.ps
    , WIAT = WIAT.sub20p
  ) %>% 
  select(newID, Age, SPA.total, SPA.reg, SPA.irr, SPA.nw
         , Morph.total, Morph.rw, Morph.ps, WIAT)
cor(ageAdj20thP %>% select(-newID), use="pairwise.complete")
