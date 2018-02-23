#install.packages('DATRAS',repos='http://www.rforge.net/',type='source')
require(DATRAS)
#datras <-  downloadExchange(survey = 'NS-IBTS')
#datras_cpue <-  downloadExchange(survey = 'CPUE per length per Hour and Swept Area', years = '2001')

lh <- read.csv('life_history.csv')

specset <- c('Sprattus sprattus',
             "Ammodytes",
             "Ammodytes marinus",
             "Ammodytes tobianus",
             "Ammodytidae",
             'Trisopterus esmarkii',
             'Clupea harengus',
             'Merlangius merlangus',
             'Solea solea',
             'Solea vulgaris',
             'Pleuronectes platessa',
             'Melanogrammus aeglefinus',
             'Pollachius virens',
             'Gadus morhua',
             'Eutrigla gurnardus')


NS <- readExchange(sprintf('surveys/NS-IBTS_%d.zip', 1965:2016))
NSset <- subset(NS, Species %in% specset)

NSset[[3]]$Species <- as.character(NSset[[3]]$Species)
NSset[[3]]$Species[grepl('Solea',NSset[[3]]$Species)] <- "Solea solea"
NSset[[3]]$Species[grepl('Ammodyt',NSset[[3]]$Species)] <- "Ammodytes marinus"

NSset[[3]] <- inner_join(NSset[[3]],lh,by=c('Species'='sciname'))

NSset[[3]] <- NSset[[3]] %>% group_by(Species) %>%
  mutate(Weight=a*LngtCm^b)

w.breaks <- 10^seq(-4,5,0.5)

nspec <- unique(NSset[[3]]$Species)
spec <- vector('list', length(nspec)); i=0
names(spec) <- nspec
for (s in nspec){
  i=i+1
  
  x <- subset(NSset,Species == s)
  x[[3]]$sizeGroup <- cut(x[[3]]$Weight, breaks = w.breaks, 
                          right = FALSE)
  N <- xtabs(Count ~ haul.id + sizeGroup, data = x[[3]])
  N <- round(N)
  x[[2]]$N <- N[as.character(x[[2]]$haul.id), , drop = FALSE]
  attr(x, "w.breaks") <- w.breaks
  
  try <- x[[2]] %>% dplyr::select(Year, N) 
  nexttry <- data.frame(try$Year,as.data.frame(try$N))
  
  midpts <- (w.breaks[1:(length(w.breaks)-1)] + w.breaks[2:(length(w.breaks))])/2
  
  spec[[i]] <- nexttry %>% group_by(try.Year,sizeGroup) %>%
    summarise(N=mean(Freq)) %>% 
    mutate(midpt = midpts,
           scale = diff(w.breaks)) %>%
    filter(try.Year!='1965')
}

meltspec <- reshape2::melt(spec, id.vars = c('try.Year','sizeGroup','midpt',"scale"))
names(meltspec) <-c('Year', 'sizeGroup', 'mid',"scale",'n','N','Species')
spectra <- meltspec %>% dplyr::select(-n,-sizeGroup)

spectra %>% 
  filter(N>0) %>%
  ggplot() +
  geom_line(aes(x=mid,y=N/scale,col=Year))+
  facet_wrap(~Species,scales = 'free_y') + 
  scale_y_log10() + 
  scale_x_log10()

spectra %>% 
  group_by(Species, Year) %>% 
  mutate(year = as.numeric(as.character(Year))) %>%
  filter(year>1972) %>%
  ggplot() + 
  scale_y_log10()+
  geom_line(aes(x=year,y=biomass)) + 
  facet_wrap(~Species, scales='free')

# group spectra of genus/species confusion

save(spectra, file='NS_survey_spectra.Rdata')

write.csv(specsum,'NS_numberspec.csv')
