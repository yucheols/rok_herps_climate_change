##### process coordinates for the project

# prevent encoding errors
Sys.getlocale()
Sys.setlocale("LC_ALL","korean")
Sys.getlocale()

# load packages
library(readxl)
library(plyr)
library(dplyr)
library(readr)

#####  process NES data
nes.coords <- list.files(path = 'occs_compiled/NES_raw', pattern = 'xlsx', full.names = T) %>%
  lapply(read_xlsx) %>%
  plyr::rbind.fill() %>% 
  dplyr::select(3,5,6,8,9)

colnames(nes.coords) = c('scientific_name', 'lat', 'long', 'year', 'region')
head(nes.coords)

# check species names
unique(nes.coords$scientific_name)

# recode some names
nes.coords$scientific_name <- dplyr::recode(nes.coords$scientific_name, 
                                            'Mauremys sinensis (Gray, 1834)' = 'Mauremys sinensis',
                                            'Pseudemys peninsularis Carr, 1938' = 'Pseudemys peninsularis',
                                            'Pseudemys concinna (Le Conte, 1830)' = 'Pseudemys concinna')

unique(nes.coords$scientific_name)

# export per species
list_sp <- unique(nes.coords$scientific_name)
length(list_sp)

for (i in 1:length(list_sp)) {
  per_sp <- nes.coords %>% dplyr::filter(scientific_name == list_sp[[i]])
  file_name <- paste0('occs_compiled/NES_per_sp/', list_sp[[i]], '.csv')
  write.csv(per_sp, file_name)
}
