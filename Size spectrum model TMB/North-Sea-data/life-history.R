require(rfishbase)

load('ICES_catch.Rdata')
spp <- unique(ices_catch$SPECIES_NAME)[-2]
lw  <- rfishbase::length_weight(spp)
g   <- rfishbase::popgrowth(spp)

sp_lw <- lw %>% filter(grepl('North Sea', Locality)) %>%
  group_by(sciname) %>% 
  summarise(a = mean(a, na.rm =T),
            b = mean(b, na.rm=T))

sp_g <- g %>% filter(grepl('North Sea', Locality)) %>%
  group_by(sciname) %>% 
  summarise(Linf = mean(Loo, na.rm =T),
            K = mean(K, na.rm=T))

life_history <- inner_join(sp_lw, sp_g) %>% mutate(Winf = a*Linf^b)

write.csv(life_history,'life_history.csv')
