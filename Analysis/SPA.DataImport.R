###############
# SPA Project - Data Importing and cleaning.
#
# Created: 2024-09-20

# 0. Set up ####

include(readxl)

## 0.1 all PPT Data ###

pptDetails = read_xlsx("spelling qualitative.cleaned.xlsx"
                   , sheet="SubjDets", range="A1:Y1054"
                   , na = c("sr_na", "nr", "NR")
                   )

# 1.0 Read Data Files ####

## 1.1 SPA data ####

spaRaw = read_xlsx("spelling qualitative.cleaned.xlsx"
                   , sheet="SPA", range="A1:ANF88", na = "sr_na") %>% 
  select(-starts_with("sr_DROP"))

# grep("\\.\\.\\.", colnames(spaRaw), value = T) %>% sort


## 1.2 SPA data "integrity" ####

spa.uniqueValsTables = apply(spaRaw %>% select(-ItemType, -Simulus) %>% 
                         mutate_all(tolower)
                       , 1, function(x) sort(table(x), decreasing = T))
# There are up to 510 different responses
maxLen = max(unlist(lapply(spa.uniqueValsTables, length)))
minLen = min(unlist(lapply(spa.uniqueValsTables, length)))
spa.uniqueValsDFs = lapply(spa.uniqueValsTables
                        , function(x) 
                          rbind(data.frame(x)# %>% mutate(Var1=as.character(Var1))
                                , data.frame(x=rep(NA, maxLen-length(x))
                                             , Freq=rep(NA, maxLen-length(x))
                                )
                          )
                        )

spa.uniqueVals = spa.uniqueValsDFs %>% bind_cols()
colnames(spa.uniqueVals) = paste(
  gl(length(spaRaw$Simulus), 2, labels = spaRaw$Simulus)
  , c("resp", "cnt"), sep="."
)

# write.table(spa.uniqueVals, 'SPA.uniqueResponses.csv', sep=";", row.names = F)

## Ok, that gives us a complete list of all codings (responses) associated with
# each variable.
### Need to check responses to make sure the values are aligned with the right
# stimuli.


### There seem to be different coding schemes as well:

spa.coding = apply(spaRaw %>% select(-ItemType, -Simulus) %>% 
                     mutate_all(tolower)
                   , 2, function(x) sort(unique(x)))

## I have replaced "." with 1, since it appears that at least 1 coder used it 
# instead.
## I have replaced ` with 1, since it seems likely to just be a typo with the 
# ` right next to 1 on the top row of many keyboards.
## I am leaving "-" as is, because it is unclear to me if it should mean 
# illegible, or no response.
## nr/NR etc... are coded as (no response)
## ? is coded as (illegible) - which means when embedded in responses it could
# look like a single illegible letter: la(illegible)(illegible)oon
# when it appears alone, it may indicate a string of illegible letters, or a 
# single illegible letter - I'm not sure.

## OK, I think that is "final" for the response data from SPA

## 1.3 SPA Accuracy coding ####

spaAcc = spaRaw %>% 
  rename(Stimulus = Simulus) %>% 
  mutate(across(!c("ItemType", "Stimulus"), ~ .==1)
         , across(!c("ItemType", "Stimulus"), ~ifelse(is.na(.), F, .))
         )

# summary(t(spaAcc %>% select(-ItemType, -Stimulus)))
### 1.3.1 Write to file ####
spaRaw %>% 
  rename(Stimulus = Simulus) %>% 
  mutate(across(!c("ItemType", "Stimulus"), ~ .==1)) %>% 
  mutate(across(c(-ItemType, -Stimulus), ~as.numeric(.))) %>% 
  write.csv(file="SPA.Accuracy.csv")

pptDetails %>% 
  write.csv("ppt.details.csv")

## 1.4 SPA scoring ####

# Scoring scheme for SPA is fairly simple, each participant is meant to be
# presented with all of the items from all 3 lists, and the score is the total #
# correct from each list.

SPA.scores = by(spaAcc %>% select(-ItemType, -Stimulus), spaAcc$ItemType
   , FUN = colSums, na.rm=T)
SPA.nas = by(spaAcc %>% select(-ItemType, -Stimulus), spaAcc$ItemType
             , FUN = function(x) colSums(!is.na(x)))

SPA = data.frame(SPA.irr = SPA.scores$Irreg, SPA.irr.resps = SPA.nas$Irreg
                 , SPA.reg = SPA.scores$Reg, SPA.reg.resps = SPA.nas$Reg
                 , SPA.nw = SPA.scores$NW, SPA.nw.resps = SPA.nas$NW) %>% 
  mutate(SPA.irr = ifelse(SPA.irr.resps == 0, NA, SPA.irr)
         , SPA.reg = ifelse(SPA.reg.resps == 0, NA, SPA.reg)
         , SPA.nw = ifelse(SPA.nw.resps == 0, NA, SPA.nw)
         , SPA.total = SPA.irr+SPA.reg+SPA.nw
         , SPA.total.resps = SPA.irr.resps+SPA.reg.resps+SPA.nw.resps)


# Ok, so there are 21 participants missing at least one of the sub scales

# This still leaves us with 1023 participants with complete data.


## OK, so for looking at item features etc... I can use all of the SPA
# data assuming each column is a unique ppt.

## However, for comparing against Morph and WIAT, I will have to limit things
# to those participants that we can feasibly match up.

