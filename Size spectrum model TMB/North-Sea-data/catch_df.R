# Catch data to use

require(fisheryO)
icc <- ices_catch_data()
head(icc)
dim(icc)

icc_ag <- icc %>% filter(ECOREGION == 'Greater North Sea Ecoregion') %>%
  group_by(YEAR, SPECIES_NAME, SPECIES_CODE) %>%
  summarise(VALUE = sum(VALUE, na.rm=T))
  
icc_sp <- icc_ag %>% group_by(YEAR) %>% 
  arrange(desc(VALUE)) %>% 
  mutate(CUM_PROP_VAL = cumsum(VALUE)/sum(VALUE,na.rm=T)) %>%
  filter(CUM_PROP_VAL<=0.9) %>% 
  arrange(YEAR) %>%
  pull(SPECIES_NAME) %>% unique(.)

ices_catch <- icc_ag %>% 
  group_by(YEAR) %>% 
  filter(!SPECIES_NAME %in% c("Crangon crangon",
                              "Mytilus edulis",
                              "Cerastoderma edule",
                              "Nephrops norvegicus"),
         SPECIES_NAME %in% c(icc_sp,
                             'Eutrigla gurnardus')) %>%
  arrange(YEAR, desc(VALUE))

ices_catch_ag <- ices_catch %>% group_by(YEAR) %>%
  summarise(sum(VALUE))

save(ices_catch, ices_catch_ag, file = 'ICES_catch.Rdata')

