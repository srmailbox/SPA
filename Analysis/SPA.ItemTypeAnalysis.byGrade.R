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
include(gridExtra)

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
           itemCors=by(
             x, x$YEAR
             , function(x) 
               cor(x %>% select(2:30)
                   , x %>% select(starts_with("sub20p"))
               )
             , simplify=F
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
                    geom_text(aes(label=round(FPR,0)), size=2.5, color="white")+
                    theme_minimal()+
                    scale_y_continuous(
                      name="Items", breaks = 1:29, labels = rownames(yr$itemFPR)
                      )+
                    scale_fill_binned(
                      breaks = c(0,5,10,15,20,100), na.value="grey90"
                      , transform="pseudo_log"
                      )
                }
         )
)

png(filename="FPRs.Irregs.png", height=960, width=1920)
grid.arrange(
  spaYrItemTypeFPRs$Irregs$Yr1+ggtitle("Irregs Yr1")
  ,spaYrItemTypeFPRs$Irregs$Yr2+ggtitle("Irregs Yr2")
  ,spaYrItemTypeFPRs$Irregs$Yr3+ggtitle("Irregs Yr3")
  ,spaYrItemTypeFPRs$Irregs$Yr4+ggtitle("Irregs Yr4")
  ,spaYrItemTypeFPRs$Irregs$Yr5+ggtitle("Irregs Yr5")
  ,spaYrItemTypeFPRs$Irregs$Yr6+ggtitle("Irregs Yr6")
  , nrow=2)
dev.off()

png(filename="FPRs.Regs.png", height=960, width=1920)
grid.arrange(
  spaYrItemTypeFPRs$Regs$Yr1+ggtitle("Regs Yr1")
,spaYrItemTypeFPRs$Regs$Yr2+ggtitle("Regs Yr2")
,spaYrItemTypeFPRs$Regs$Yr3+ggtitle("Regs Yr3")
,spaYrItemTypeFPRs$Regs$Yr4+ggtitle("Regs Yr4")
,spaYrItemTypeFPRs$Regs$Yr5+ggtitle("Regs Yr5")
,spaYrItemTypeFPRs$Regs$Yr6+ggtitle("Regs Yr6")
,nrow=2)
dev.off()

png(filename="FPRs.NWs.png", height=960, width=1920)
grid.arrange(
  spaYrItemTypeFPRs$NWs$Yr1+ggtitle("NWs Yr1")
,spaYrItemTypeFPRs$NWs$Yr2+ggtitle("NWs Yr2")
,spaYrItemTypeFPRs$NWs$Yr3+ggtitle("NWs Yr3")
,spaYrItemTypeFPRs$NWs$Yr4+ggtitle("NWs Yr4")
,spaYrItemTypeFPRs$NWs$Yr5+ggtitle("NWs Yr5")
,spaYrItemTypeFPRs$NWs$Yr6+ggtitle("NWs Yr6")
,nrow=2)
dev.off()


## let's look at rankings by year ####
# I don't love that items pop in and out of "importance" with this method. Treating
# each grade as fully distinct allows for that in a way that using QR didn't as much
# in the age-based version. e.g., it's odd to me that "nature" is one of the best
# predictors for Yrs 4 and 6, but not even top 10 for year 5.

lapply(spaYrItemTypes
       , function (x) lapply(x$itemCors, function(y) y %>% data.frame) %>% bind_cols
       )

# OK, welp, this is just kind of the reality, I guess. It's either accept this
# or make highly subjective changes guided by the feeling that items should
# behave well.

### 2.3.1 OK, let's set the scales ####
# Irregs:
# Yr1: No Options - Irregs are at the floor.
# Yr2: 5 items, threshold 2 (e.g., must score higher than 2)
#   said, come, was, have, other (4% FPR)
# Yr3: deadly young dairy please other people search country shoulder (station)
#   9 (10) items, threshold 3 (must score 4 or better - FPR 5% (4%))
# Yr4: search nature dairy young remember country other flavour fortune station
#   10 items, threshold 6 (must score 7 - FPR 13%)
#   (could add people then please, increasing threshold by 1 and then 1 more)
#   15 items (people please deadly aspire shoulder) would give 7% FPR (9 thresh)
# Yr5: need a lot - 14 (16) items
#   search country remember flavour science station please shoulder fortune 
#   language young deadly nature neighbour (dairy wheat)
#   thresh 10 (11) FPR 15% (12%)
# Yr6: fortune nature search country shoulder neighbour (aspire young dairy)
#   6 (9) items with threshold 4 (6) gets FPR 8% (6%)

# Regs:
# Yr1: 5 items (best cup clip stem grand) with a threshold of 1 gets 2% FPR
# Yr2: grand strong best clip fond life - 6 items, thresh 3, fpr 3%
# Yr3: understand sporty blade subside shelf ground prize yesterday selfishly scarper
#   10 items, threshold 5, fpr 10%
# Yr4: strong sporty blade shelf understand clip cup best
#   8 items threshold 6 fpr 3%
# Yr5: check subside property selfishly prize yesterday scarper demonstrate 
#   shelf judgement impact trombone swept scram
#   14 items, thresh 9, fpr 10%
# Yr6: property impact judgement subside blade demonstrate lantern selfishly
#   8 items, thresh 5, fpr 8%

# NWs:
# Yr1: swem clep lep prond scrand - 5 items, thresh 1, fpr 0%
# Yr2: prond swem greem chust linst = 5 items, thresh 2, fpr 3%
# Yr3:  chust scade swem underbably jubside bife strofect trelfishly clep fize 
#       lorty scrand withound prostand
#   items 14, thresh 6, fpr 12%
# Yr4:  clep swem chust shilf greem lep scade prond strofect prostand jubside 
#       bife lorty spantern scrand
#   items 15, thresh 8, fpr 6% (thresh 9, 15%)
# Yr5:  scade shilf trelfishly strofect prostand greem impabit sandstrate 
#       sombone lorty fize withound prond lirst
#   14 items thresh 9 fpr 9%
# Yr 6: underbably strofect sandstrate jubside shilf lorty greem prostand scade 
#       grong chust bife (sombone blarper spantern)
#   12 items (15), thresh 9 (10), fpr 15% (12%)


### Set up the rules
scrnrStruct = 
data.frame(
  Type = gl(3, 6, labels=c("Irreg", "Regs", "NWs"))
  , Year = paste0("Yr", 1:6)
  , nItems = c(0,5,11,10,14,6
               ,5,6,10,8,14,8
               ,5,5,14,15,14,12)
  , thrshhld = c(NA,2,3,6,10,4
                 , 1,3,5,6,9,5
                 , 1,2,6,8,9,9)
)

scrnrItems = lapply(
  spaYrItemTypes
  , function(typ) 
    lapply(
      typ$itemCors
      , function(yr) rownames(yr)[order(yr)]
    ) %>% bind_cols
)

scrnrStructList = vector("list", nrow(scrnrStruct))
for (i in 1:nrow(scrnrStruct)) {
  if(scrnrStruct$nItems[i]) {
    Type = scrnrStruct$Type[i]
    Year = scrnrStruct$Year[i]
    nItems = scrnrStruct$nItems[i]
    Threshhold = scrnrStruct$thrshhld[i]
    scrnrStructList[[i]] = 
      list(Type = Type, Year = Year
           , nItems = nItems, Threshhold=Threshhold
           , items = unlist(scrnrItems[[Type]][1:nItems, Year]))
  }
}

scrnrScores = lapply(
  scrnrStructList
  , function(x) {
    if(is.null(x)) return(x)
    varName = paste0("scrnr20p.",x$Type)
    dat = spaYrAnalysis %>% 
      filter(YEAR==x$Year) %>% 
      select(newID, all_of(x$items)) %>% 
      mutate(subscale = ifelse(rowSums(pick(-newID))<=x$Threshhold, 1,0))
    colnames(dat)[ncol(dat)]=varName
    list(Year = x$Year, scrnr = dat %>% select(newID, ncol(dat)))
  }
)

scrnrScoresYr6 = 
  lapply(scrnrScores
         , function (x) if (!is.null(x) && x$Year=="Yr6") x$scrnr else NULL
         ) %>% bind_cols

scrnrScoresYr5 = 
  lapply(scrnrScores
         , function (x) if (!is.null(x) && x$Year=="Yr5") x$scrnr else NULL
  ) %>% bind_cols

scrnrScoresYr4 = 
  lapply(scrnrScores
         , function (x) if (!is.null(x) && x$Year=="Yr4") x$scrnr else NULL
  ) %>% bind_cols

scrnrScoresYr3 = 
  lapply(scrnrScores
         , function (x) if (!is.null(x) && x$Year=="Yr3") x$scrnr else NULL
  ) %>% bind_cols

scrnrScoresYr2 = 
  lapply(scrnrScores
         , function (x) if (!is.null(x) && x$Year=="Yr2") x$scrnr else NULL
  ) %>% bind_cols

scrnrScoresYr1 = 
  lapply(scrnrScores
         , function (x) if (!is.null(x) && x$Year=="Yr1") x$scrnr else NULL
  ) %>% bind_cols %>% mutate(scrnr20p.Irreg = NA)

scrnrScoresAll = 
  rbind(
    scrnrScoresYr1 %>% select(newID...1, starts_with("scrnr")) %>% 
      rename(newID = newID...1)
    , scrnrScoresYr2 %>% select(newID...1, starts_with("scrnr")) %>% 
      rename(newID = newID...1)
    , scrnrScoresYr3 %>% select(newID...1, starts_with("scrnr")) %>% 
      rename(newID = newID...1)
    , scrnrScoresYr4 %>% select(newID...1, starts_with("scrnr")) %>% 
      rename(newID = newID...1)
    , scrnrScoresYr5 %>% select(newID...1, starts_with("scrnr")) %>% 
      rename(newID = newID...1)
    , scrnrScoresYr6 %>% select(newID...1, starts_with("scrnr")) %>% 
      rename(newID = newID...1)
  )

detectionFPRCheck = 
  merge(
    spaYrAnalysis %>% select(newID, YEAR, starts_with("sub20p"))
    , scrnrScoresAll
  ) %>% 
  mutate(
    anySub20p = rowSums(pick(starts_with("sub20p")), na.rm=T)>0
    , anyScrnr20p = rowSums(pick(starts_with("scrnr")), na.rm=T)>0
  )


Contingencies = by(detectionFPRCheck, detectionFPRCheck$YEAR
   , function(x) 
     list(
       Irregs = table(x$sub20p_irr, x$scrnr20p.Irreg)
       , Regs = table(x$sub20p_reg, x$scrnr20p.Regs)
       , NWs = table(x$sub20p_nws, x$scrnr20p.NWs)
       , All = table(x$anySub20p, x$anyScrnr20p)
       )
   )
