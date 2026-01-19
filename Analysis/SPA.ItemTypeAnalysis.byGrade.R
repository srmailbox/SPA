###############
# SPA Project - SPA test item analysis - Item types
#  Grade -based screeners
#
# In the first pass at creating a screener, I focused on a screener that would
# identify "overall" at risk kids. This led to poor performance on Item subtypes
# e.g., kids who were poor at Nonwords, but not other word types. 
# (SPA.ItemAnalysis.R). Recently, I did the same but by grade 
# (SPA.ItemAnalysis.ByGrade.R).
#
# We had ultimately decided to do a screener that tested for each item type 
# (Irregs, Regs, NW) and then aggregated these into overall screeners
# (SPA.ItemTypeAnalysis.R).
# Here I repeat that process but Grade-based again.
# 
# REST OF THE DOCUMENTATION OF THE ORIGINAL
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
# Created: 2026-01-19
# Changelog:

# 0. Set up ####
# include(CTT)
# include(psych)
# include(zoo)
# include(quantreg)

if(!exists("spaAcc")) source("Analysis/SPA.DataImport.R")

# 1.0 reshape data and split by item type ####

## 1.1 all items, regardless of type ####

spaYrAll = t(spaAcc %>% select(-ItemType, -Stimulus)) %>% 
  data.frame %>% 
  mutate(across(everything(), as.numeric)
         , score = rowSums(.)
         # , score_logis = qlogis((score+1)/89)
         , newID = rownames(.)
  )%>% 
  merge(pptDetails %>% 
          select(newID, YEAR, sp_age) %>% 
          mutate(
            YEAR = ordered(YEAR, labels=paste0("Yr", 0:6))
            , sp_age = as.numeric(sp_age)
          )
  )
colnames(spaYrAll)[2:88] = spaAcc$Stimulus

# Limit to kids between 6;5 and 12;0
spaYrAnalysis = spaYrAll %>% filter(YEAR != "Yr0", sp_age < 200) %>% 
  mutate(YEAR = ordered(YEAR))
# ggplot(spaYrAnalysis, aes(x=sp_age, y=jitter(YEAR)))+geom_point()

### 1.1.1 Identify poor spellers (overall) ####
# In the age based analysis, this used Quantile Regression, but here we will
# just use simple percentile ranks.

spaYrAnalysis = spaYrAnalysis %>% 
  group_by(YEAR) %>% 
  mutate(sub20p_tot=rank(score, ties.method="min")/n()
         , sub20p_tot = ifelse(sub20p_tot>.2, 0, 1)
         ) %>% 
  ungroup()

## 1.2 items by type ####
spaYrAnalysis =
  spaYrAnalysis %>% 
  rowwise %>% 
  mutate(
    score_irr =  sum(across(2:30))
    , score_reg = sum(across(31:59))
    , score_nws = sum(across(60:88))
    ) %>% 
  ungroup()

## 1.3 Identify 20th percentiles by Item type and Grade ####

spaYrAnalysis =
  spaYrAnalysis %>% 
  group_by(YEAR) %>% 
  mutate(
    across(
      starts_with("score_")
      , ~ifelse(rank(., ties.method="min")/n()>.2, 0,1)
      , .names="sub20p_{.col}"
           )
    ) %>% 
  ungroup
colnames(spaYrAnalysis)=gsub("_score_", "_", colnames(spaYrAnalysis))

## 1.4 Comorbidity ####

spaYrAnalysis %>% group_by(sub20p_tot, sub20p_irr, sub20p_nws, sub20p_reg) %>% 
  summarize(n())


# 2. Using tetrachoric cors to identify item lists ####

# first create a list of datasets by type.

spaYrItemTypes = 
  list(Irregs = spaYrAnalysis[,c(1, 2:30, 90,93,96)]
       , Regs = spaYrAnalysis[,c(1, 2:30+29, 90,94,97)]
       , NWs = spaYrAnalysis[,c(1, 2:30+58, 90,95,98)]
  )

## 2.1 calculate correlations for each outcome ####
spaYrItemTypes = 
  lapply(spaYrItemTypes
         , function(x) list(
           data=x, 
           itemCors=by(x, x$YEAR
              , function(x) cor(x %>% select(2:30)
                                , x %>% select(starts_with("sub20p"))
              )
           )
         )
  )
  
## 2.3 Cumulative scores & thresholds ####
spaYr.detection = function(x) {
  detection = lapply(x$itemCors, function(y) y %>% data.frame %>% 
                       arrange(pick(starts_with("sub20"))))
  fpr = detection
  for(yr in names(x$itemCors)) {
    dat = x$data %>% filter(YEAR == yr)
    nSub20p = sum(dat[,33])
    nOk = nrow(dat)-nSub20p
    detection[[yr]] = data.frame(detection[[yr]], matrix(NA, 29, 30))
    colnames(detection[[yr]])[-1]=paste0("l", 0:29)
    fpr[[yr]] = detection[[yr]]
    for(sclLength in 1:29) {
      itms = rownames(detection[[yr]])[1:sclLength]
      itmsDat = dat %>% select(all_of(itms), starts_with("sub20p")) %>% 
        mutate(score = rowSums(across(all_of(itms))))
      for(thrsh in 1:sclLength-1){
        nDetect = itmsDat %>% 
          filter(if_any(starts_with("sub20p"), ~.==1),score<=thrsh) %>% nrow
        detection[[yr]][sclLength,thrsh+2] = itmsDat %>% 
          filter(if_any(starts_with("sub20p"), ~.==1),score<=thrsh) %>% nrow
        if(nDetect == nSub20p) fpr[[yr]][sclLength,thrsh+2] = 
          round(
            (itmsDat %>% filter(if_any(starts_with("sub20p"), ~.==0)
                                , score<=thrsh) %>% nrow)/nOk
            ,3)*100
      }
    }
    minItems = sum(rowMaxes(detection[[yr]], na.rm=T)<=max(detection[[yr]], na.rm=T)-1)+1
    detection[[yr]] = list(
      itemDetection = detection[[yr]]
      , itemFPR = fpr[[yr]]
      , length = minItems
      , items = rownames(detection[[yr]])[1:minItems]
    )
    
  }
  
  # This structure is terrible. 
  list(data=x$data, detection = detection)
}

spaYrItemTypesScreeners = lapply(spaYrItemTypes, spaYr.detection)

## OK, I need to find all the combos that produce perfect detection,
## but relatively low "false positives"

# Heatmaps seem like a quick way to visualize this.
# For a given itemFPR matrix, this is what I want to generate:

spaYrItemTypeFPRs = lapply(spaYrItemTypesScreeners
       , function (type)
         lapply(type$detection,
                function(yr) {
                  yr$itemFPR %>% 
                  select(-starts_with("sub20p")) %>% mutate(testLength=1:29) %>% 
                    pivot_longer(-testLength, values_to = "FPR", names_to="thresh") %>% 
                    mutate(thresh = ordered(thresh, levels=paste0("l", 0:29))) %>% 
                    ggplot(aes(thresh, testLength, fill=FPR))+
                    geom_tile(na.rm=T)+
                    theme_minimal()+
                    scale_y_continuous(breaks = 0:30)+
                    scale_fill_binned(breaks = 0:5*20, na.value="grey90")
                }
         )
)

par(mfrow=c(2,3))
png(filename="FPRs.Irregs.png", width=960)
spaYrItemTypeFPRs$Irregs$Yr1+ggtitle("Irregs Yr1")
spaYrItemTypeFPRs$Irregs$Yr2+ggtitle("Irregs Yr2")
spaYrItemTypeFPRs$Irregs$Yr3+ggtitle("Irregs Yr3")
spaYrItemTypeFPRs$Irregs$Yr4+ggtitle("Irregs Yr4")
spaYrItemTypeFPRs$Irregs$Yr5+ggtitle("Irregs Yr5")
spaYrItemTypeFPRs$Irregs$Yr6+ggtitle("Irregs Yr6")
dev.off()

png(filename="FPRs.Regs.png", width=960)
spaYrItemTypeFPRs$Regs$Yr1+ggtitle("Regs Yr1")
spaYrItemTypeFPRs$Regs$Yr2+ggtitle("Regs Yr2")
spaYrItemTypeFPRs$Regs$Yr3+ggtitle("Regs Yr3")
spaYrItemTypeFPRs$Regs$Yr4+ggtitle("Regs Yr4")
spaYrItemTypeFPRs$Regs$Yr5+ggtitle("Regs Yr5")
spaYrItemTypeFPRs$Regs$Yr6+ggtitle("Regs Yr6")
dev.off()

png(filename="FPRs.NWs.png", width=960)
spaYrItemTypeFPRs$NWs$Yr1+ggtitle("NWs Yr1")
spaYrItemTypeFPRs$NWs$Yr2+ggtitle("NWs Yr2")
spaYrItemTypeFPRs$NWs$Yr3+ggtitle("NWs Yr3")
spaYrItemTypeFPRs$NWs$Yr4+ggtitle("NWs Yr4")
spaYrItemTypeFPRs$NWs$Yr5+ggtitle("NWs Yr5")
spaYrItemTypeFPRs$NWs$Yr6+ggtitle("NWs Yr6")
dev.off()
