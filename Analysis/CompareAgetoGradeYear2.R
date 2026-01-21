ageDx = spaAnalysis %>% select(newID, sp_age, starts_with("sub20p"))
yrDx = spaYrAnalysis %>% select(newID, YEAR, starts_with("sub20p"))

compareDx = merge(ageDx, yrDx, by="newID", suffixes=c(".age", ".yr"))

compareDx %>% #group_by(YEAR) %>% 
  select(sub20p_tot.age, sub20p_tot.yr, YEAR) %>% table


table(spaYrAnalysis$sp_age, spaYrAnalysis$YEAR)

compareDx %>% filter(YEAR=="Yr2", sub20p_tot.age==1, sub20p_tot.yr==0) %>% 
  summarize(range(sp_age))
