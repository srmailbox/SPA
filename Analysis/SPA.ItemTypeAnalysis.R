###############
# SPA Project - SPA test item analysis - Item types
#
# In the first pass at creating a screener, I focused on a screener that would
# identify "overall" at risk kids. This led to poor performance on Item subtypes
# e.g., kids who were poor at Nonwords, but not other word types. 
# (SPA.ItemAnalysis.R)
#
# In this version, I develop screeners for each item type (Irregs, Regs, NW)
# and then aggregate these into overall screeners.
# Final approach - identify each item subscale for each age, but choose a number 
# of items and threshold that keeps false positives low.
# 
# In this version, I order the items by how well they predict poor performance 
# on that item type for that age group. Then for each age group and item type, 
# I consider all lengths and threshold scores to calculate a) the detection rate
# for at risk (20th percentile) on that item type and age, and b) if the 
# detection rate is perfect, I calculate the false positive rate for that list 
# length and threshold.
# 
# Then, using my eyeballs, I choose a length and threshold that keeps the false 
# positive rate below 10% or as low as possible (sometimes, if one extra item 
# improves the score, I do that).
#
# Created: 2025-05-15
# Changelog:

# 0. Set up ####
# include(CTT)
# include(psych)
include(zoo)
include(quantreg)

if(!exists("spaAcc")) source("Analysis/SPA.DataImport.R")

# 1.0 reshape data and split by item type ####

## 1.1 all items, regardless of type ####

spaAll = t(spaAcc %>% select(-ItemType, -Stimulus)) %>% 
  data.frame %>% 
  mutate(across(everything(), as.numeric)
         , score = rowSums(.)
         , score_logis = qlogis((score+1)/89)
         , newID = rownames(.)
  )%>% 
  merge(pptDetails %>% 
          select(newID, sp_age) %>% 
          mutate(sp_age = as.numeric(sp_age)
                 , ageCat = cut(sp_age, breaks=c(0,89, 101, 113, 125, 150)
                                , labels=c("<= 7;5", paste0("<= ",8:10,";",5), "10;6+"))
          )
  )
colnames(spaAll)[2:88] = spaAcc$Stimulus

# Limit to kids between 6;5 and 12;0
spaAnalysis = spaAll %>% filter(sp_age < 200, sp_age >=78)

### 1.1.1 Identify poor spellers (overall) ####
spaAnalysis.qr.2 = rq(score_logis ~ sp_age, spaAnalysis, tau=.2)
spaAnalysis$qr.fit = 87*plogis(predict(spaAnalysis.qr.2))
spaAnalysis$sub20p_tot = ifelse(spaAnalysis$qr.fit>spaAnalysis$score, 1, 0)


## 1.2 items by type ####
spaItemTypes = 
  list(Irregs = spaAnalysis[,c(1, 2:30, 89:94)]
       , Regs = spaAnalysis[,c(1, 2:30+29, 89:94)]
       , NWs = spaAnalysis[,c(1, 2:30+58, 89:94)]
  ) %>% 
  lapply(function(x) {
    x %>% 
      mutate(subscale_score = rowSums(across(all_of(2:30)))
             , subscale_logis = qlogis((subscale_score+1)/31)
             )
  })

# 769 kids

## 1.3 QR analysis to identify 20th percentiles ####
spaItemTypes.qr.2 = lapply(spaItemTypes
                           , function(x) {
                             mod = rq(subscale_logis~sp_age, x, tau=.20)
                             x$qr.fit = 29*plogis(predict(mod))
                             x$sub20p_subscale = ifelse(x$subscale_score < x$qr.fit, 1, 0)
                             list(qr = mod, data=x)
                           })

spaAnalysis = 
  merge(
    spaAnalysis
    , spaItemTypes.qr.2$Irregs$data %>% select(newID, subscale_score, sub20p_subscale) %>% 
      rename(score_irr = subscale_score, sub20p_irr = sub20p_subscale)
    , by="newID"
  ) %>% 
  merge(
    spaItemTypes.qr.2$Regs$data %>% select(newID, subscale_score, sub20p_subscale) %>% 
      rename(score_reg = subscale_score, sub20p_reg = sub20p_subscale)
    , by="newID"
  ) %>% 
  merge(
    spaItemTypes.qr.2$NWs$data %>% select(newID, subscale_score, sub20p_subscale) %>% 
      rename(score_nws = subscale_score, sub20p_nws = sub20p_subscale)
    , by="newID"
  ) %>% 
  mutate(sub20p = rowMaxes(across(starts_with("sub20p"))))


## 1.4 Comorbidity ####

spaAnalysis %>% group_by(sub20p, sub20p_tot, sub20p_irr, sub20p_nws, sub20p_reg) %>% 
  summarize(n())

### OK, so any child that is "poor" overall, is also "poor" according
# to at least one of the three subscales (and most are poor on all 3)

# This means that we only need to identify the poor performers on each subscale
# This will "catch" all of the "overall" poor spellers as well.

# 2. Using tetrachoric cors to identify item lists ####

spa.tetrachoric = function(x) {
  itemCols = 2:30
  byAge = by(x$data
             , INDICES=x$data$ageCat
             , function(y) cor(y %>% select(all_of(itemCols))
                               , y$sub20p_subscale) %>% 
               data.frame %>% 
               rename(cor = ".") %>% arrange(cor)
             , simplify = F)
  list(data=x$data, qr=x$qr, itemCors = byAge)
  
}

spaItemTypes.qr.2 = lapply(spaItemTypes.qr.2, spa.tetrachoric)

## 2.1 Cumulative scores and thresholds ####
spa.detection.age = function(x) {
  detection = lapply(x$itemCors, function(y) y)
  fpr = detection
  for(age in names(x$itemCors)) {
    dat = merge(x$data, spaAnalysis %>% select(newID, sub20p)) %>% filter(ageCat == age)
    nSub20p = sum(dat$sub20p_subscale)
    nOk = nrow(dat)-nSub20p
    detection[[age]] = data.frame(detection[[age]], matrix(NA, 29, 30))
    fpr[[age]] = data.frame(fpr[[age]], matrix(NA, 29, 30))
    colnames(detection[[age]])[-1]=paste0("l", 0:29)
    colnames(fpr[[age]])[-1]=paste0("l", 0:29)
    for(sclLength in 1:29) {
      itms = rownames(detection[[age]])[1:sclLength]
      itmsDat = dat %>% select(all_of(itms), sub20p, sub20p_subscale) %>% 
        mutate(score = rowSums(across(all_of(itms))))
      for(thrsh in 1:sclLength-1){
        nDetect = itmsDat %>% 
          filter(sub20p_subscale==1,score<=thrsh) %>% nrow
        detection[[age]][sclLength,thrsh+2] = itmsDat %>% 
          filter(sub20p_subscale==1,score<=thrsh) %>% nrow
        if(nDetect == nSub20p) fpr[[age]][sclLength,thrsh+2] = round((itmsDat %>% 
          filter(sub20p==0,score<=thrsh) %>% nrow)/nOk,3)*100
      }
    }
      minItems = sum(rowMaxes(detection[[age]], na.rm=T)<=max(detection[[age]], na.rm=T)-1)+1
      detection[[age]] = list(
        itemDetection = detection[[age]]
        , itemFPR = fpr[[age]]
        , length = minItems
        , items = rownames(detection[[age]])[1:minItems]
      )
    
  }
  
  # This structure is terrible. I really want the screeners to go in the age
  # summaries.
  list(data=x$data, qr=x$qr, detection = detection)
}

spaItemTypesScreeners = lapply(spaItemTypes.qr.2, spa.detection.age)

### Ok, so looking at the "FPR" for each item type:

spaItemTypesScreeners$Irregs$detection$`<= 7;5`$itemFPR
spaItemTypesScreeners$Regs$detection$`<= 7;5`$itemFPR
spaItemTypesScreeners$NWs$detection$`<= 7;5`$itemFPR
# Needs 4/11 irregs, we'd get 8.5% FPR
# 7/11 regs 6.4%
# 3/7 NWs to get to 6.1% (3/6 gives 9.8% FPR)
# That's a 29 item scale.

screeners = vector("list", 5)
names(screeners) = names(spaItemTypesScreeners$Irregs$detection)
screeners[["<= 7;5"]] = 
  list(
    Age = "<= 7;5"
    , Irreg = 
      list(
        items = rownames(spaItemTypesScreeners$Irregs$detection$`<= 7;5`$itemFPR)[1:11]
        , thresh = 4
        )
    , Regs = 
      list(
        items = rownames(spaItemTypesScreeners$Regs$detection$`<= 7;5`$itemFPR)[1:11]
        , thresh = 7
      )
    , NWs = 
      list(
        items = rownames(spaItemTypesScreeners$NWs$detection$`<= 7;5`$itemFPR)[1:7]
        , thresh = 3
      )
  )

spaItemTypesScreeners$Irregs$detection$`<= 8;5`$itemFPR
spaItemTypesScreeners$Regs$detection$`<= 8;5`$itemFPR
spaItemTypesScreeners$NWs$detection$`<= 8;5`$itemFPR
# Needs 2/6 irregs, we'd get 6.9% FPR
# 4/6 regs 7.5%
# 3/8 NWs to get to 9.8
# That's a 20 item scale.

screeners[["<= 8;5"]] = 
  list(
    Age = "<= 8;5"
    , Irreg = 
      list(
        items = rownames(spaItemTypesScreeners$Irregs$detection$`<= 8;5`$itemFPR)[1:6]
        , thresh = 2
      )
    , Regs = 
      list(
        items = rownames(spaItemTypesScreeners$Regs$detection$`<= 8;5`$itemFPR)[1:6]
        , thresh = 4
      )
    , NWs = 
      list(
        items = rownames(spaItemTypesScreeners$NWs$detection$`<= 8;5`$itemFPR)[1:8]
        , thresh = 3
      )
  )

spaItemTypesScreeners$Irregs$detection$`<= 9;5`$itemFPR[1:15,1:15]
spaItemTypesScreeners$Regs$detection$`<= 9;5`$itemFPR[1:15,1:15]
spaItemTypesScreeners$NWs$detection$`<= 9;5`$itemFPR[1:15,1:15]
# Needs 2/5 irregs, we'd get 4.2% FPR
# 4/6 regs 2.5%
# 5/9 NWs to get to 5%
# That's a 20 item scale.

screeners[["<= 9;5"]] = 
  list(
    Age = "<= 9;5"
    , Irreg = 
      list(
        items = rownames(spaItemTypesScreeners$Irregs$detection$`<= 9;5`$itemFPR)[1:5]
        , thresh = 2
      )
    , Regs = 
      list(
        items = rownames(spaItemTypesScreeners$Regs$detection$`<= 9;5`$itemFPR)[1:6]
        , thresh = 4
      )
    , NWs = 
      list(
        items = rownames(spaItemTypesScreeners$NWs$detection$`<= 9;5`$itemFPR)[1:9]
        , thresh = 5
      )
  )

spaItemTypesScreeners$Irregs$detection$`<= 10;5`$itemFPR[1:15,1:15]
spaItemTypesScreeners$Regs$detection$`<= 10;5`$itemFPR[1:15,1:15]
spaItemTypesScreeners$NWs$detection$`<= 10;5`$itemFPR[,10:30]
# Needs 2/5 irregs, we'd get 2.1% FPR 
# 10/12 regs 7.4%
# 11/18 NWs to get to 9.8% (11/19 gives 3%)
# That's a 35 item scale.

screeners[["<= 10;5"]] = 
  list(
    Age = "<= 10;5"
    , Irreg = 
      list(
        items = rownames(spaItemTypesScreeners$Irregs$detection$`<= 10;5`$itemFPR)[1:5]
        , thresh = 2
      )
    , Regs = 
      list(
        items = rownames(spaItemTypesScreeners$Regs$detection$`<= 10;5`$itemFPR)[1:12]
        , thresh = 10
      )
    , NWs = 
      list(
        items = rownames(spaItemTypesScreeners$NWs$detection$`<= 10;5`$itemFPR)[1:18]
        , thresh = 11
      )
  )

spaItemTypesScreeners$Irregs$detection$`10;6+`$itemFPR[,1:20]
spaItemTypesScreeners$Regs$detection$`10;6+`$itemFPR[,1:25]
spaItemTypesScreeners$NWs$detection$`10;6+`$itemFPR[,1:25]
# Needs 8/11 irregs, gets 10.1% (as good as it gets)
# 15/20 regs to get to 8.3% (10/14 gives 11.7%)
# 10/15 NWs to get to 9.5% (11/19 gives 3%)
# That's a 46 item scale.

screeners[["10;6+"]] = 
  list(
    Age = "10;6+"
    , Irreg = 
      list(
        items = rownames(spaItemTypesScreeners$Irregs$detection$`10;6+`$itemFPR)[1:11]
        , thresh = 8
      )
    , Regs = 
      list(
        items = rownames(spaItemTypesScreeners$Regs$detection$`10;6+`$itemFPR)[1:20]
        , thresh = 15
      )
    , NWs = 
      list(
        items = rownames(spaItemTypesScreeners$NWs$detection$`10;6+`$itemFPR)[1:15]
        , thresh = 10
      )
  )

spaScreened = 
  lapply(
    screeners
    , function(x) {
      dat = spaAnalysis %>% filter(ageCat==x$Age)
      dat %>% 
        mutate(scr_irr_score = rowSums(across(x$Irreg$items))
               , scr_irr_sub20p = as.numeric(scr_irr_score <= x$Irreg$thresh)
               , scr_reg_score = rowSums(across(x$Regs$items))
               , scr_reg_sub20p = as.numeric(scr_reg_score <= x$Regs$thresh)
               , scr_nws_score = rowSums(across(x$NWs$items))
               , scr_nws_sub20p = as.numeric(scr_nws_score <= x$NWs$thresh)
               , scr_sub20p = rowMaxes(across(contains("_sub20p")))
               )
      
    }
  ) %>% bind_rows

### FPR
table(spaScreened$sub20p, spaScreened$scr_sub20p)

# Overall this would lead to 103/522 people being flagged for follow up ~ 20%

### 2.2 Output ####
lapply(screeners, function(x) lapply(x[-1], function(y) y$items) %>% unlist %>% data.frame(Item=.))

# 3.0 ####
spaScreened %>% group_by(sub20p_tot, ageCat) %>% summarise(across(contains("sub20p"), mean))
# Just out of curiosity, I wanted to see which subscales are most associated with "total" at risk
# and it turns out that overall poor performance predicts a bad Regular item score more
# than poor Irregs or NWs.

#################################################

spaScreenerItems = lapply(spaItemTypesScreeners, 
       function(x) lapply(x$detection, function(y) y$items))

itmsByAge = vector("list", length(levels(spaItemTypes$Irregs$ageCat)))
names(itmsByAge) = levels(spaItemTypes$Irregs$ageCat)
for (age in names(itmsByAge)) {
  itmsByAge[[age]] = c(
    spaScreenerItems$Irregs[[age]], spaScreenerItems$Regs[[age]]
    , spaScreenerItems$NWs[[age]]
  )
}

# 3.0 assess item type derived screeners against "total" ####

qrMod = data.frame(sp_age=78:144) %>% 
  mutate(score=87*plogis(predict(spaAnalysis.qr.2, data.frame(sp_age = 78:144))))

ggplot(spaAnalysis, aes(x=sp_age, y=score, colour=sub20p))+
  geom_point(alpha=.5)+
  # geom_smooth()+
  geom_line(data=qrMod, aes(x=sp_age, colour=NULL), size=1, col="darkred")+
  theme_minimal()+
  labs(x="Age in months", y="Score (of 87)"
       , title="20th Percentile Performance", subtitle="as a function of age")

for(age in names(itmsByAge)){
  spaAnalysis$screener.score[spaAnalysis$ageCat==age] = 
    rowSums(spaAnalysis[spaAnalysis$ageCat==age, itmsByAge[[age]]])
}

table(spaAnalysis$screener.score, spaAnalysis$ageCat, spaAnalysis$sub20p)[,,2]
# so for each subscale, we can apply a threshold as follows:
# 7;5: 6/9 (15/80 false alarms, 18.8%) worse than total (12.5)
# 8;5: 4/7 (19/134 false alarms, 14.2%) better (21.6)
# 9;5: 9/13 (14/120 false alarms, 11.7%) same (11.7)
# 10;5: 13/16 (32/136 false alarms, 23.5%) worse (20.6)
# 11;0: 15/18 (56/144 false alarms, 38.9%... oof.) much worse (14.6)

## 3.1 Additional items for "total" ####

# For each age group, take the items that are not in the screener yet
# and choose the ones that correlate strongest with the "20th percentile"

# get the item correlations
spaItemCorsByAge = lapply(
  spaItemTypes.qr.2
  , function(x) {
    itms = x$itemCors
    class(itms)="list"
    itms %>% bind_cols()
  }
  ) %>% bind_rows()

colnames(spaItemCorsByAge) = c(paste0("<= ", c("7;5", "8;5", "9;5", "10;5")), "10;6+")

# get the items that are already on the scale.
spaScreenerItemsByAge = spaScreenerItems %>% unlist(recursive = T, use.names = T)
spaScreenerItemsByAge = data.frame(Item = spaScreenerItemsByAge
           , label = names(spaScreenerItemsByAge)) %>% 
  mutate(
    Scale = gsub("\\..*", "", label)
    , Age = ifelse(
      gsub(".*\\.", "", label) == "<= 7;5"
      , gsub(".*\\.", "", label)
      , substr(gsub(".*\\.", "", label), 1, str_length(gsub(".*\\.", "", label))-1)
    )
  ) %>% select(-label)


## Get the items that were in the original screener (ignoring Item type)
addItems = vector("list", 5)
names(addItems) = colnames(spaItemCorsByAge)
for(age in names(addItems)) {
  usedItems = spaScreenerItemsByAge %>% filter(Age==age) %>% select(Item) %>% unlist
  addItems[[age]] = spaItemCorsByAge %>% select(age) %>% 
    rename(Cor = 1) %>% 
    mutate(Item = rownames(.), Scale = "Total", Age = age) %>% 
    filter(!Item %in% usedItems) %>% 
    arrange(Cor) %>% 
    slice(1:10)
}

# Final screeners
spaScreenerLong = rbind(
  spaScreenerItemsByAge
  , addItems %>% bind_rows %>% select(Item, Scale, Age)
)

write.csv(spaScreenerItemsByAge, "spaItemTypeScreener.p20.csv", row.names = F)
write.csv(spaScreenerLong, "spaLongScreener.p20.csv", row.names = F)

### Check quality of the "Total" screeners
for(age in names(addItems)){
  spaAnalysis$longScreener.score[spaAnalysis$ageCat==age] = 
    rowSums(
      spaAnalysis[spaAnalysis$ageCat==age
                  , spaScreenerLong$Item[spaScreenerLong$Age==age]]
      )
}

table(spaAnalysis$longScreener.score, spaAnalysis$ageCat, spaAnalysis$sub20p)[,,2]
# so for each subscale, we can apply a threshold as follows:
# 7;5: 11/19 (10/80 false alarms, 12.5%) (same "Total")
# 8;5: 13/17 (30/134 false alarms, 22.4%) (worse, keep Items)
# 9;5: 16/23 (6/120 false alarms, 5%) (much better than either)
# 10;5: 23/26 (42/136 false alarms, 30.8%) (keep Items)
# 11;0: 21/28 (12/144 false alarms, 8.3%) (much better than either)

#### original version ignoring subtypes:
# This keeps false alarms under 25% for everyone. 
# 7;5: 12.5%
# 8;5: 21.6%
# 9;5: 11.7%
# 10;5: 20.6%
# 10;6+: 14.6%

