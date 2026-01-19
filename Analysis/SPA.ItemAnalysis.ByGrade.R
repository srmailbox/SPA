###############
# SPA Project - SPA test item analysis - Age Screeners ####
# 
# Final results here are just based on simple tetrachoric correlations by Age Group
#
# Created: 2026-01-16
# Changelog:

# 0. Set up ####
include(psych)

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

itemCols = colnames(spaAll)
spaResults = merge(spaAll %>% mutate(newID = rownames(.)
                                     , across(all_of(itemCols), ~ifelse(is.na(.), 0, .)))
                   , pptDetails %>% select(newID, YEAR, starts_with("sp_")) %>%
                     # rename(newID = Child_ID, sp_age=sp.age) %>%
                     mutate(across(starts_with("sp_"), as.numeric))
                   , by = "newID") %>%
  # rowwise() %>%
  # The pptDetails scores are not accurate based on the item data, so this corrects for
  # that.
  mutate(sp_alt = rowSums(select(., all_of(itemCols)))
         , sp_irr_alt = rowSums(select(., all_of(itemCols[1:29])))
         , sp_reg_alt = rowSums(select(., all_of(itemCols[1:29+29])))
         , sp_nw_alt = rowSums(select(., all_of(itemCols[1:29+58])))
         )

# remove YEAR 0 (kindy) kids, and the one older kid with special needs
spaYear = spaResults %>% filter(YEAR!=0, sp_age<200) %>% 
  select(-sp_tot, -sp_irr, -sp_reg, -sp_nw) %>% 
  mutate(YEAR = ordered(YEAR))

# 2.0 ID 20th percentile on each subscale, by year ####
# Also identifies kids in the bottom 20p

spaYear = spaYear %>% 
  group_by(YEAR) %>% 
  mutate(across(c(starts_with("sp_"), -sp_age), ~quantile(., probs=.2, type=1)
                , .names = "{.col}.20p")) %>% 
  ungroup %>% 
  mutate(
    sp_alt_sub20 = sp_alt <= sp_alt.20p
    , sp_irr_sub20 = sp_irr_alt <= sp_irr_alt.20p
    , sp_reg_sub20 = sp_reg_alt <= sp_reg_alt.20p
    , sp_nw_sub20 = sp_nw_alt <= sp_nw_alt.20p
  ) %>% 
  mutate(across(ends_with("sub20"), ~ ifelse(., 1, 0))) %>% 
  data.frame

## 2.1 Thresholds ####
spaYear %>% select(YEAR, ends_with("20p")) %>% unique %>% 
  arrange(YEAR)
# Floor performance for Year 1 across the board, and Year 2 irregs

## 2.2 "prevalence" ####
spaYear %>% group_by(YEAR) %>% summarize(across(ends_with("sub20"), mean))

# Massive overestimation of "20th percentile" for Year 1 kids on NWs.
# Mainly because they are at floor - the threshold is 0, and 48% of kids did not
# spell a single nonword correctly.
# Irregs also produce somewhat elevated values for years 1 and 2
# Oddly regs don't pose an issue for Year 1, but that is probably volatile given
# the small number of kids in that Year

# 3.0 Using tetrachoric cors to identify item lists ####

# In this case, the idea is to identify which items are most predictive of
# whether a kid is likely to fall in the "sub20" category.

spaYearCor = cor(spaYear %>% select(all_of(itemCols)), spaYear$sp_alt_sub20) %>% data.frame %>% 
  rename(cor = ".") %>% arrange(cor) %>% slice(1:20)

## 3.1 by Year ####
spaYearCor.Year = by(spaYear, INDICES=spaYear$YEAR
   , function(x) cor(x %>% select(all_of(itemCols)), x$sp_alt_sub20) %>% data.frame %>% 
     rename(cor = ".")# %>% arrange(cor)
   , simplify = F) %>% 
  lapply(FUN=as.data.frame) %>% bind_cols() 

colnames(spaYearCor.Year) = paste0("cor", 1:6)

### 3.1.1 age-based screeners ####

spaYear.screeners=apply(
  spaYearCor.Year, 2
  , function(x) {
    spaYear %>% 
      mutate(screener = rowSums(select(., all_of(names(x[order(x)])[1:20])))) %>% 
      select(newID, screener, YEAR, sp_alt_sub20)
  }) %>% bind_cols() %>%
  rename(newID = newID...1, sub20=sp_alt_sub20...8, Year = YEAR...3
         , screenerY1 = screener...2, screenerY2 = screener...6
         , screenerY3 = screener...10, screenerY4 = screener...14
         , screenerY5 = screener...18, screenerY6 = screener...22) %>% 
  select(newID, Year, sub20, starts_with("screenerY")) %>%
  mutate(screener.yr =
           case_when(Year==1~screenerY1, Year==2~screenerY2
                     , Year==3~screenerY3, Year==4~screenerY4
                     , Year==5~screenerY5, Year==6~screenerY6)
         # , screener.IRTflag = as.numeric(screener.yr <= thrshHlds[as.numeric(Year)])
  )

by(spaYear.screeners, INDICES=spaYear.screeners$Year
   , FUN=function(x) cor(x$sub20, x$screener.yr))

# At the total level, the correlation between the derived screener scores and
# the probability of being a 20th percentile kid is reasonably high
# running from .74 (Year 1) to .87 (Year 4)

## 3.2 assess age-based screeners ####

spaYear.detection = spaYear.screeners %>% 
  group_by(Year, screener.yr) %>% 
  summarise(detection=sum(sub20)) %>% #slice(1:20) %>% data.frame()
  # ungroup(screener) %>% 
  mutate(detection=cumsum(detection)/sum(detection)) %>% 
  ungroup() %>% 
  arrange(Year, screener.yr)

### 3.2.1 Detection ####
ggplot(spaYear.detection, aes(x=screener.yr, y=detection, col=Year))+
  geom_line()+
  labs(title = "Screener's rate of detecting Full SPA 20th percentile kids"
       , subtitle = "by Year and threshold score (using tetrachoric correlation)"
       , colour = "Year"
       , x = "Screener threshold"
       , y = "Detect sub 20th percentile") +
  theme_minimal()+
  scale_x_continuous(breaks = 0:20)+
  scale_y_continuous(minor_breaks = 0:20/20)+
  theme(panel.grid.minor.x=element_blank())

thrshHlds.ttrchrc = spaYear.detection %>%
  group_by(Year) %>% 
  filter(detection==1) %>% 
  summarize(thrshHld.tc = min(screener.yr))
### Ok, so with different items for each age, we could use the following
# Year  thrshHld.tc
# <ord>       <dbl>
# 1 1               3
# 2 2              10
# 3 3              12
# 4 4              15
# 5 5              16
# 6 6              16

### 3.2.2 False positive rate ####

spaYear.screeners = 
  merge(spaYear.screeners %>% data.frame, thrshHlds.ttrchrc %>% data.frame) %>% 
  mutate(
    screener.tcflag = 
      as.numeric(screener.yr <= thrshHld.tc)
    )

spaYear.fpr = spaYear.screeners %>% 
  group_by(Year, screener.yr) %>%
  filter(!sub20) %>% 
  summarise(n = n()) %>% #slice(1:20) %>% data.frame()
  # ungroup(screener) %>%
  mutate(fpr=cumsum(n)) %>% 
  mutate(fpr=fpr/max(fpr)) %>% 
  ungroup() %>% 
  arrange(Year, screener.yr)

ggplot(spaYear.fpr, aes(x=screener.yr, y=fpr, col=Year))+
  geom_vline(data=thrshHlds.ttrchrc
             , mapping=aes(xintercept=thrshHld.tc, col=Year)
             , size=.25, linetype="dashed")+
  geom_line()+
  labs(title = "False Positive Rate (Overdiagnosis)"
       , subtitle = "by Year and threshold score - Year-specific Items"
       , colour = "Year"
       , x = "Screener threshold"
       , y = "FPR") +
  theme_minimal()+
  scale_x_continuous(breaks = 0:20)+
  scale_y_continuous(minor_breaks = 0:20/20)+
  theme(panel.grid.minor.x=element_blank())

merge(thrshHlds.ttrchrc, spaYear.fpr, by.y=c("Year", "screener.yr")
      , by.x=c("Year", "thrshHld.tc"), all.x=T, all.y=F)
# This keeps false alarms under 25% for everyone. 
# Year thrshHld.tc  n        fpr
# 1    1           3  1 0.04761905
# 2    2          10  4 0.07446809
# 3    3          12  7 0.20540541
# 4    4          15  6 0.09836066
# 5    5          16 18 0.32773109
# 6    6          16  9 0.23300971
# 
# Not bad for Year 1, 2, and 4. A bit high for 3 and 6
# Rough for 5.

### 3.2.3 "Saved" assessments ####

# How many kids would not meet the threshold for follow-up

table(spaYear.screeners$sub20, spaYear.screeners$screener.tcflag)

# So across the dataset of 769 kids
# You would not follow up 543/665 (82%) kids who are not at risk (correct reject)
# You would correctly follow up on 176/176 (100%) kids who are at risk (hit)**
# Of the kids you followed up, 122/298 (or 41%) would not be at risk (false alarm)
# Of the kids you did not follow up, 0/542 (or 0% would be misses)**
# ** these values are essentially forced to be perfect by the process.

# Overall you would reduce your follow up assessments by 543/841 or 62%

#### 3.2.3.1 Effect of Base Rate ####
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

## 3.3 Final list output ####
# Here are the items:
spaYearagescreener.items =
  spaYeartetrachoric.age %>% 
  apply(2
        , FUN= function(x) names(x[order(x)])[1:20])

write.csv(spaYearagescreener.items, "spaScreener.p20.20Item.csv", row.names=F)

spaYearitems.long = spaYearagescreener.items %>% data.frame() %>% 
  pivot_longer(cols=1:5, names_to = "age", values_to = "item") %>% 
  mutate(cor=NA)
for (i in 1:nrow(spaYearitems.long))
  spaYearitems.long$cor[i] = spaYeartetrachoric.age[spaYearitems.long$item[i]
                                           , spaYearitems.long$age[i]]


(spaYearitems.long %>% 
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

table(cut(which(order(spaYeartetrachoric.age$cor7)<=20), breaks=c(0,29, 58, 87), labels = c("I", "R", "N")))
# 5 Irr, 4 Reg, 11 NW
table(cut(which(order(spaYeartetrachoric.age$cor8)<=20), breaks=c(0,29, 58, 87), labels = c("I", "R", "N")))
# 9 Irr, 9 Reg, 2 NW
table(cut(which(order(spaYeartetrachoric.age$cor9)<=20), breaks=c(0,29, 58, 87), labels = c("I", "R", "N")))
# 5 Irr, 10 Reg, 5 NW
table(cut(which(order(spaYeartetrachoric.age$cor10)<=20), breaks=c(0,29, 58, 87), labels = c("I", "R", "N")))
# 3 Irr, 8 Reg, 9 NW
table(cut(which(order(spaYeartetrachoric.age$cor11)<=20), breaks=c(0,29, 58, 87), labels = c("I", "R", "N")))
# 8 Irr, 8 Reg, 4 NW

# So for each subsale, we can identify the poorest readers using the same Quant 
# Regression technique, then apply the same ROC curve approach to see how well 
# the screeners do at catching kids who might be ok overall, but poor at one type

# This is just sensitivity (we don't care as much about specificity here, 
# because the screener will only be flagging them for follow up - their specific 
# deficits can be assessed more directly later)

### 9.1.1 Irregular Words ####

#### 9.1.1.1 QR for IDing "20th percentile" ####

# ggplot(spaYear, aes(x=sp_alt, y=ability, col=sp_age))+geom_point()
spaYear$sp_irr_logit = qlogis((spaYear$sp_irr_alt+1)/31)
spaYear.qr_irr = rq(sp_irr_logit~sp_age, spaYear, tau=.20)

spaYear$qr_irr.fit = 29*plogis(predict(spaYear.qr_irr))

qrMod.irr = data.frame(sp_age=78:144) %>% 
  mutate(sp_irr_alt=29*plogis(predict(spaYear.qr_irr, data.frame(sp_age = 78:144))))

ggplot(spaYear, aes(x=sp_age, y=sp_irr_alt))+
  geom_point(col="darkgrey")+
  # geom_smooth()+
  geom_line(data=qrMod.irr, aes(colour=NULL), size=1, col="blue")+
  theme_minimal()+
  labs(x="Age in months", y="Score (of 29)"
       , title="20th Percentile Performance (Irregs)", subtitle="as a function of age")

# So, let's flag the poorest 20% of SPA spellers (age adjusted)

spaYear = spaYear %>% 
  mutate(sub20.irr = ifelse(sp_irr_alt < qr_irr.fit, 1, 0))
# table(spaYear$sub20.irr)/nrow(spaYear)

#### 9.1.1.2 Sensitivity ####

spaYearscreeners = spaYearscreeners %>% 
  merge(spaYear %>% select(newID, sub20.irr))

spaYeardetection.agescreeners.irr = spaYearscreeners %>% 
  group_by(Year, screener.yr) %>% 
  summarise(detection.irr=sum(sub20.irr)) %>% #slice(1:20) %>% data.frame()
  # ungroup(screener) %>% 
  mutate(detection.irr=cumsum(detection.irr)/sum(detection.irr)) %>% 
  ungroup() %>% 
  arrange(Year, screener.yr)

ggplot(spaYeardetection.agescreeners.irr
       , aes(x=screener.yr, y=detection.irr, col=Year))+
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

spaYeardetection.agescreeners.irr %>% 
  filter(detection.irr == 1 | (Year == "<= 7;5" & screener.yr == 11) | 
  (Year == "<= 8;5" & screener.yr == 9) |
  (Year == "<= 9;5" & screener.yr == 12) |
  (Year == "<= 10;5" & screener.yr == 14) | 
  (Year == "10;6+" & screener.yr == 16)) %>% 
  data.frame %>% 
  group_by(Year) %>% slice(1:2)
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

spaYear$sp_nw_logit = qlogis((spaYear$sp_nw_alt+1)/31)
(spaYear.qr_nw = rq(sp_nw_logit~sp_age, spaYear, tau=.20)) %>% summary

spaYear$qr_nw.fit = 29*plogis(predict(spaYear.qr_nw))

qrMod.nw = data.frame(sp_age=78:144) %>% 
  mutate(sp_nw_alt=29*plogis(predict(spaYear.qr_nw, data.frame(sp_age = 78:144))))

ggplot(spaYear, aes(x=sp_age, y=sp_nw_alt))+
  geom_point(col="darkgrey")+
  # geom_smooth()+
  geom_line(data=qrMod.nw, aes(colour=NULL), size=1, col="blue")+
  theme_minimal()+
  labs(x="Age in months", y="Score (of 29)"
       , title="20th Percentile Performance (Regs)", subtitle="as a function of age")

# So, let's flag the poorest 20% of SPA spellers (age adjusted)

spaYear = spaYear %>% 
  mutate(sub20.nw = ifelse(sp_nw_alt < qr_nw.fit, 1, 0))
# table(spaYear$sub20.nw)/nrow(spaYear)

#### 9.1.2.2 Sensitivity ####

spaYearscreeners = spaYearscreeners %>% 
  merge(spaYear %>% select(newID, sub20.nw))

spaYeardetection.agescreeners.nw = spaYearscreeners %>% 
  group_by(Year, screener.yr) %>% 
  summarise(detection.nw=sum(sub20.nw)) %>% #slice(1:20) %>% data.frame()
  # ungroup(screener) %>% 
  mutate(detection.nw=cumsum(detection.nw)/sum(detection.nw)) %>% 
  ungroup() %>% 
  arrange(Year, screener.yr)

ggplot(spaYeardetection.agescreeners.nw
       , aes(x=screener.yr, y=detection.nw, col=Year))+
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

spaYeardetection.agescreeners.nw %>% 
  filter(detection.nw == 1 | (Year == "<= 7;5" & screener.yr == 11) | 
           (Year == "<= 8;5" & screener.yr == 9) |
           (Year == "<= 9;5" & screener.yr == 12) |
           (Year == "<= 10;5" & screener.yr == 14) | 
           (Year == "10;6+" & screener.yr == 16)) %>% 
  data.frame %>% 
  group_by(Year) %>% slice(1:2)
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

spaYear$sp_reg_logit = qlogis((spaYear$sp_reg_alt+1)/31)
(spaYear.qr_reg = rq(sp_reg_logit~sp_age, spaYear, tau=.20)) %>% summary

spaYear$qr_reg.fit = 29*plogis(predict(spaYear.qr_reg))

qrMod.reg = data.frame(sp_age=78:144) %>% 
  mutate(sp_reg_alt=29*plogis(predict(spaYear.qr_reg, data.frame(sp_age = 78:144))))

ggplot(spaYear, aes(x=sp_age, y=sp_reg_alt))+
  geom_point(col="darkgrey")+
  # geom_smooth()+
  geom_line(data=qrMod.reg, aes(colour=NULL), size=1, col="blue")+
  theme_minimal()+
  labs(x="Age in months", y="Score (of 29)"
       , title="20th Percentile Performance (Regular)", subtitle="as a function of age")

# So, let's flag the poorest 20% of SPA spellers (age adjusted)

spaYear = spaYear %>% 
  mutate(sub20.reg = ifelse(sp_reg_alt < qr_reg.fit, 1, 0))

#### 9.1.3.2 Sensitivity ####

spaYearscreeners = spaYearscreeners %>% 
  merge(spaYear %>% select(newID, sub20.reg))

spaYeardetection.agescreeners.reg = spaYearscreeners %>% 
  group_by(Year, screener.yr) %>% 
  summarise(detection.reg=sum(sub20.reg)) %>% #slice(1:20) %>% data.frame()
  # ungroup(screener) %>% 
  mutate(detection.reg=cumsum(detection.reg)/sum(detection.reg)) %>% 
  ungroup() %>% 
  arrange(Year, screener.yr)

ggplot(spaYeardetection.agescreeners.reg
       , aes(x=screener.yr, y=detection.reg, col=Year))+
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

spaYeardetection.agescreeners.reg %>% 
  filter(detection.reg == 1 | (Year == "<= 7;5" & screener.yr == 11) | 
           (Year == "<= 8;5" & screener.yr == 9) |
           (Year == "<= 9;5" & screener.yr == 12) |
           (Year == "<= 10;5" & screener.yr == 14) | 
           (Year == "10;6+" & screener.yr == 20)) %>% 
  data.frame %>% 
  group_by(Year) %>% slice(1:2) %>% 
  mutate(res = c("detection", "threshold")
         , val = ifelse(res=="detection", detection.reg, screener.yr)) %>%
  select(Year, res, val) %>% 
  pivot_wider(id_cols=Year, names_from=res, values_from=val) %>% 
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

spaMorphWiat = merge(spaYearscreeners
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
  mutate(WIAT.sub20 = ifelse(WIAT.ss <=85, 1, 0))

WIATdetection.agescreeners = spaMorphWiat %>% 
  group_by(Year, screener.yr) %>% 
  summarise(detection=sum(WIAT.sub20, na.rm=T)) %>% #slice(1:20) %>% data.frame()
  # ungroup(screener) %>% 
  mutate(detection=cumsum(detection)/sum(detection)) %>% 
  ungroup() %>% 
  arrange(Year, screener.yr)

ggplot(WIATdetection.agescreeners
       , aes(x=screener.yr, y=detection, col=Year))+
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
  filter(detection == 1 | (Year == "<= 7;5" & screener.yr == 11) | 
           (Year == "<= 8;5" & screener.yr == 9) |
           (Year == "<= 9;5" & screener.yr == 12) |
           (Year == "<= 10;5" & screener.yr == 14) | 
           (Year == "10;6+" & screener.yr == 16)) %>% 
  data.frame %>% 
  group_by(Year) %>% slice(1:2)
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
  mutate(sub20.morph.tot = ifelse(Morph.tot < morph.fit, 1, 0))

ggplot(spaMorphWiat, aes(x=sp_age, y=Morph.tot, col=sub20.morph.tot))+
  geom_point()+
  # geom_smooth()+
  geom_line(data=qrMod.morph.tot, aes(colour=NULL), size=1, col="blue")+
  theme_minimal()+
  labs(x="Age in months", y="Score (of 53)"
       , title="20th Percentile Performance (Morph Total)", subtitle="as a function of age")

#### 9.3.1.2 Sensitivity ####
# Problem - some scores don't show up in the Morph data.
morphScreeners = spaYearscreeners %>% 
  filter(newID %in% spaMorphWiat$newID[!is.na(spaMorphWiat$Morph.tot)]) %>% 
  select(newID, Year, starts_with("screener")) %>% 
  merge(spaMorphWiat %>% select(newID, sub20.morph.tot))

Morphdetection.agescreeners = expand.grid(
  Year = c('<= 7;5', '<= 8;5', '<= 9;5', '<= 10;5', '10;6+')
  , screener.yr = 0:20
) %>% data.frame()

Morphdetection.agescreeners.morph.tot = morphScreeners %>% 
  group_by(Year, screener.yr) %>% 
  summarise(detection.morph=sum(sub20.morph.tot)) %>% #slice(1:20) %>% data.frame()
  # ungroup(screener) %>% 
  mutate(detection.morph=cumsum(detection.morph)/sum(detection.morph, na.rm=T)) %>% 
  merge(Morphdetection.agescreeners, all.y=T) %>% 
  group_by(Year) %>% 
  mutate(
    detection.morph = na.locf(detection.morph, na.rm=FALSE)
    , detection.morph = ifelse(is.na(detection.morph), 0, detection.morph)
  ) %>% 
  ungroup() %>% 
  arrange(Year, screener.yr)

ggplot(Morphdetection.agescreeners.morph.tot
       , aes(x=screener.yr, y=detection.morph, col=Year))+
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
  arrange(detection.morph, screener.yr) %>% 
  filter(detection.morph == 1 | (Year == "<= 7;5" & screener.yr == 11) | 
           (Year == "<= 8;5" & screener.yr == 9) |
           (Year == "<= 9;5" & screener.yr == 12) |
           (Year == "<= 10;5" & screener.yr == 14) | 
           (Year == "10;6+" & screener.yr == 16)) %>% 
  data.frame %>% 
  group_by(Year) %>% slice(1:2) %>% 
  mutate(res = c("detection", "threshold")
         , val = ifelse(res=="detection", detection.morph, screener.yr)) %>%
  select(Year, res, val) %>% 
  pivot_wider(id_cols=Year, names_from=res, values_from=val) %>% 
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
  mutate(sub20.morph.rw = ifelse(Morph.rw < morph.rw.fit, 1, 0))

ggplot(spaMorphWiat, aes(x=sp_age, y=Morph.rw, col=sub20.morph.rw))+
  geom_point()+
  # geom_smooth()+
  geom_line(data=qrMod.morph.rw, aes(colour=NULL), size=1, col="blue")+
  theme_minimal()+
  labs(x="Age in months", y="Score (of 29)"
       , title="20th Percentile Performance (Morph Words)", subtitle="as a function of age")

#### 9.3.2.2 Sensitivity ####

morphScreeners = morphScreeners %>% 
  filter(newID %in% spaMorphWiat$newID[!is.na(spaMorphWiat$Morph.rw)]) %>% 
  #select(newID, Year, starts_with("screener")) %>% 
  merge(spaMorphWiat %>% select(newID, sub20.morph.rw))

# Morphdetection.agescreeners.morph.rw = morphScreeners %>% 
#   group_by(Year, screener.yr) %>% 
#   summarise(detection.morph.rw=sum(sub20.morph.rw)) %>% #slice(1:20) %>% data.frame()
#   # ungroup(screener) %>% 
#   mutate(detection.morph=cumsum(detection.morph.rw)/sum(detection.morph.rw, na.rm=T)) %>% 
#   ungroup() %>% 
#   arrange(Year, screener.yr)

Morphdetection.agescreeners.morph.rw = morphScreeners %>% 
  group_by(Year, screener.yr) %>% 
  summarise(detection.morph.rw=sum(sub20.morph.rw)) %>% #slice(1:20) %>% data.frame()
  # ungroup(screener) %>% 
  mutate(detection.morph.rw=cumsum(detection.morph.rw)/sum(detection.morph.rw, na.rm=T)) %>% 
  merge(Morphdetection.agescreeners, all.y=T) %>% 
  group_by(Year) %>% 
  mutate(
    detection.morph.rw = na.locf(detection.morph.rw, na.rm=FALSE)
    , detection.morph.rw = ifelse(is.na(detection.morph.rw), 0, detection.morph.rw)
  ) %>% 
  ungroup() %>% 
  arrange(Year, screener.yr)

ggplot(Morphdetection.agescreeners.morph.rw
       , aes(x=screener.yr, y=detection.morph.rw, col=Year))+
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
  filter(detection.morph.rw == 1 | (Year == "<= 7;5" & screener.yr == 11) | 
           (Year == "<= 8;5" & screener.yr == 9) |
           (Year == "<= 9;5" & screener.yr == 12) |
           (Year == "<= 10;5" & screener.yr == 14) | 
           (Year == "10;6+" & screener.yr == 16)) %>% 
  data.frame %>% 
  group_by(Year) %>% slice(1:2) %>% 
  mutate(res = c("detection", "threshold")
         , val = ifelse(res=="detection", detection.morph.rw, screener.yr)) %>%
  select(Year, res, val) %>% 
  pivot_wider(id_cols=Year, names_from=res, values_from=val) %>% 
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
  mutate(sub20.morph.ps = ifelse(Morph.ps < morph.ps.fit, 1, 0))

ggplot(spaMorphWiat, aes(x=sp_age, y=Morph.ps, col=sub20.morph.ps))+
  geom_point()+
  # geom_smooth()+
  geom_line(data=qrMod.morph.ps, aes(colour=NULL), size=1, col="blue")+
  theme_minimal()+
  labs(x="Age in months", y="Score (of 29)"
       , title="20th Percentile Performance (Morph Pseudowords)", subtitle="as a function of age")

#### 9.3.3.2 Sensitivity ####

morphScreeners = morphScreeners %>% 
  filter(newID %in% spaMorphWiat$newID[!is.na(spaMorphWiat$Morph.ps)]) %>% 
  #select(newID, Year, starts_with("screener")) %>% 
  merge(spaMorphWiat %>% select(newID, sub20.morph.ps))

# Morphdetection.agescreeners.morph.ps = morphScreeners %>% 
#   group_by(Year, screener.yr) %>% 
#   summarise(detection.morph.ps=sum(sub20.morph.ps)) %>% #slice(1:20) %>% data.frame()
#   # ungroup(screener) %>% 
#   mutate(detection.morph=cumsum(detection.morph.ps)/sum(detection.morph.ps, na.rm=T)) %>% 
#   ungroup() %>% 
#   arrange(Year, screener.yr)

Morphdetection.agescreeners.morph.ps = morphScreeners %>% 
  group_by(Year, screener.yr) %>% 
  summarise(detection.morph.ps=sum(sub20.morph.ps)) %>% #slice(1:20) %>% data.frame()
  # ungroup(screener) %>% 
  mutate(detection.morph.ps=cumsum(detection.morph.ps)/sum(detection.morph.ps, na.rm=T)) %>% 
  merge(Morphdetection.agescreeners, all.y=T) %>% 
  group_by(Year) %>% 
  mutate(
    detection.morph.ps = na.locf(detection.morph.ps, na.rm=FALSE)
    , detection.morph.ps = ifelse(is.na(detection.morph.ps), 0, detection.morph.ps)
  ) %>% 
  ungroup() %>% 
  arrange(Year, screener.yr)

ggplot(Morphdetection.agescreeners.morph.ps
       , aes(x=screener.yr, y=detection.morph.ps, col=Year))+
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
  filter(detection.morph.ps == 1 | (Year == "<= 7;5" & screener.yr == 11) | 
           (Year == "<= 8;5" & screener.yr == 9) |
           (Year == "<= 9;5" & screener.yr == 12) |
           (Year == "<= 10;5" & screener.yr == 14) | 
           (Year == "10;6+" & screener.yr == 16)) %>% 
  data.frame %>% 
  group_by(Year) %>% slice(1:2) %>% 
  mutate(res = c("detection", "threshold")
         , val = ifelse(res=="detection", detection.morph.ps, screener.yr)) %>%
  select(Year, res, val) %>% 
  pivot_wider(id_cols=Year, names_from=res, values_from=val) %>% 
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
  mutate(sub20.morph.tot = ifelse(Morph.tot < morph.fit, 1, 0))

ggplot(spaMorphWiat, aes(x=sp_age, y=Morph.tot, col=sub20.morph.tot))+
  geom_point()+
  # geom_smooth()+
  geom_line(data=qrMod.morph.tot, aes(colour=NULL), size=1, col="blue")+
  theme_minimal()+
  labs(x="Age in months", y="Score (of 58)"
       , title="20th Percentile Performance (Morph Total)", subtitle="as a function of age")

#### 9.4.1.2 Sensitivity ####
# Problem - some scores don't show up in the Morph data.
morphScreeners = spaYearscreeners %>% 
  filter(newID %in% spaMorphWiat$newID[!is.na(spaMorphWiat$Morph.tot)]) %>% 
  select(newID, Year, starts_with("screener")) %>% 
  merge(spaMorphWiat %>% select(newID, sub20.morph.tot))

Morphdetection.agescreeners = expand.grid(
  Year = c('<= 7;5', '<= 8;5', '<= 9;5', '<= 10;5', '10;6+')
  , screener.yr = 0:20
) %>% data.frame()

Morphdetection.agescreeners.morph.tot = morphScreeners %>% 
  group_by(Year, screener.yr) %>% 
  summarise(detection.morph=sum(sub20.morph.tot)) %>% #slice(1:20) %>% data.frame()
  # ungroup(screener) %>% 
  mutate(detection.morph=cumsum(detection.morph)/sum(detection.morph, na.rm=T)) %>% 
  merge(Morphdetection.agescreeners, all.y=T) %>% 
  group_by(Year) %>% 
  mutate(
    detection.morph = na.locf(detection.morph, na.rm=FALSE)
    , detection.morph = ifelse(is.na(detection.morph), 0, detection.morph)
  ) %>% 
  ungroup() %>% 
  arrange(Year, screener.yr)

ggplot(Morphdetection.agescreeners.morph.tot
       , aes(x=screener.yr, y=detection.morph, col=Year))+
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
  arrange(detection.morph, screener.yr) %>% 
  filter(detection.morph == 1 | (Year == "<= 7;5" & screener.yr == 11) | 
           (Year == "<= 8;5" & screener.yr == 9) |
           (Year == "<= 9;5" & screener.yr == 12) |
           (Year == "<= 10;5" & screener.yr == 14) | 
           (Year == "10;6+" & screener.yr == 16)) %>% 
  data.frame %>% 
  group_by(Year) %>% slice(1:2) %>% 
  mutate(res = c("detection", "threshold")
         , val = ifelse(res=="detection", detection.morph, screener.yr)) %>%
  select(Year, res, val) %>% 
  pivot_wider(id_cols=Year, names_from=res, values_from=val) %>% 
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
  mutate(sub20.morph.rw = ifelse(Morph.rw < morph.rw.fit, 1, 0))

ggplot(spaMorphWiat, aes(x=sp_age, y=Morph.rw, col=sub20.morph.rw))+
  geom_point()+
  # geom_smooth()+
  geom_line(data=qrMod.morph.rw, aes(colour=NULL), size=1, col="blue")+
  theme_minimal()+
  labs(x="Age in months", y="Score (of 29)"
       , title="20th Percentile Performance (Morph Words)", subtitle="as a function of age")

#### 9.4.2.2 Sensitivity ####

morphScreeners = morphScreeners %>% 
  filter(newID %in% spaMorphWiat$newID[!is.na(spaMorphWiat$Morph.rw)]) %>% 
  #select(newID, Year, starts_with("screener")) %>% 
  merge(spaMorphWiat %>% select(newID, sub20.morph.rw))

# Morphdetection.agescreeners.morph.rw = morphScreeners %>% 
#   group_by(Year, screener.yr) %>% 
#   summarise(detection.morph.rw=sum(sub20.morph.rw)) %>% #slice(1:20) %>% data.frame()
#   # ungroup(screener) %>% 
#   mutate(detection.morph=cumsum(detection.morph.rw)/sum(detection.morph.rw, na.rm=T)) %>% 
#   ungroup() %>% 
#   arrange(Year, screener.yr)

Morphdetection.agescreeners.morph.rw = morphScreeners %>% 
  group_by(Year, screener.yr) %>% 
  summarise(detection.morph.rw=sum(sub20.morph.rw)) %>% #slice(1:20) %>% data.frame()
  # ungroup(screener) %>% 
  mutate(detection.morph.rw=cumsum(detection.morph.rw)/sum(detection.morph.rw, na.rm=T)) %>% 
  merge(Morphdetection.agescreeners, all.y=T) %>% 
  group_by(Year) %>% 
  mutate(
    detection.morph.rw = na.locf(detection.morph.rw, na.rm=FALSE)
    , detection.morph.rw = ifelse(is.na(detection.morph.rw), 0, detection.morph.rw)
  ) %>% 
  ungroup() %>% 
  arrange(Year, screener.yr)

ggplot(Morphdetection.agescreeners.morph.rw
       , aes(x=screener.yr, y=detection.morph.rw, col=Year))+
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
  filter(detection.morph.rw == 1 | (Year == "<= 7;5" & screener.yr == 11) | 
           (Year == "<= 8;5" & screener.yr == 9) |
           (Year == "<= 9;5" & screener.yr == 12) |
           (Year == "<= 10;5" & screener.yr == 14) | 
           (Year == "10;6+" & screener.yr == 16)) %>% 
  data.frame %>% 
  group_by(Year) %>% slice(1:2) %>% 
  mutate(res = c("detection", "threshold")
         , val = ifelse(res=="detection", detection.morph.rw, screener.yr)) %>%
  select(Year, res, val) %>% 
  pivot_wider(id_cols=Year, names_from=res, values_from=val) %>% 
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
  mutate(sub20.morph.ps = ifelse(Morph.ps < morph.ps.fit, 1, 0))

ggplot(spaMorphWiat, aes(x=sp_age, y=Morph.ps, col=sub20.morph.ps))+
  geom_point()+
  # geom_smooth()+
  geom_line(data=qrMod.morph.ps, aes(colour=NULL), size=1, col="blue")+
  theme_minimal()+
  labs(x="Age in months", y="Score (of 29)"
       , title="20th Percentile Performance (Morph Pseudowords)", subtitle="as a function of age")

#### 9.4.3.2 Sensitivity ####

morphScreeners = morphScreeners %>% 
  filter(newID %in% spaMorphWiat$newID[!is.na(spaMorphWiat$Morph.ps)]) %>% 
  #select(newID, Year, starts_with("screener")) %>% 
  merge(spaMorphWiat %>% select(newID, sub20.morph.ps))

# Morphdetection.agescreeners.morph.ps = morphScreeners %>% 
#   group_by(Year, screener.yr) %>% 
#   summarise(detection.morph.ps=sum(sub20.morph.ps)) %>% #slice(1:20) %>% data.frame()
#   # ungroup(screener) %>% 
#   mutate(detection.morph=cumsum(detection.morph.ps)/sum(detection.morph.ps, na.rm=T)) %>% 
#   ungroup() %>% 
#   arrange(Year, screener.yr)

Morphdetection.agescreeners.morph.ps = morphScreeners %>% 
  group_by(Year, screener.yr) %>% 
  summarise(detection.morph.ps=sum(sub20.morph.ps)) %>% #slice(1:20) %>% data.frame()
  # ungroup(screener) %>% 
  mutate(detection.morph.ps=cumsum(detection.morph.ps)/sum(detection.morph.ps, na.rm=T)) %>% 
  merge(Morphdetection.agescreeners, all.y=T) %>% 
  group_by(Year) %>% 
  mutate(
    detection.morph.ps = na.locf(detection.morph.ps, na.rm=FALSE)
    , detection.morph.ps = ifelse(is.na(detection.morph.ps), 0, detection.morph.ps)
  ) %>% 
  ungroup() %>% 
  arrange(Year, screener.yr)

ggplot(Morphdetection.agescreeners.morph.ps
       , aes(x=screener.yr, y=detection.morph.ps, col=Year))+
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
  filter(detection.morph.ps == 1 | (Year == "<= 7;5" & screener.yr == 11) | 
           (Year == "<= 8;5" & screener.yr == 9) |
           (Year == "<= 9;5" & screener.yr == 12) |
           (Year == "<= 10;5" & screener.yr == 14) | 
           (Year == "10;6+" & screener.yr == 16)) %>% 
  data.frame %>% 
  group_by(Year) %>% slice(1:2) %>% 
  mutate(res = c("detection", "threshold")
         , val = ifelse(res=="detection", detection.morph.ps, screener.yr)) %>%
  select(Year, res, val) %>% 
  pivot_wider(id_cols=Year, names_from=res, values_from=val) %>% 
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
  select(newID, sp_age, contains("sub20")) %>% 
  rename(
    Age = sp_age
    , SPA.total = sub20, SPA.reg = sub20.reg, SPA.irr = sub20.irr
    , SPA.nw = sub20.nw
    , Morph.total = sub20.morph.tot, Morph.rw = sub20.morph.rw
    , Morph.ps = sub20.morph.ps
    , WIAT = WIAT.sub20
  ) %>% 
  select(newID, Age, SPA.total, SPA.reg, SPA.irr, SPA.nw
         , Morph.total, Morph.rw, Morph.ps, WIAT)
cor(ageAdj20thP %>% select(-newID), use="pairwise.complete")
