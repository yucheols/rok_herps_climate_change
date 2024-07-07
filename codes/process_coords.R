##### process coordinates for the project

# prevent encoding errors
Sys.getlocale()
Sys.setlocale("LC_CTYPE", ".1251")
Sys.getlocale()

# load packages
library(readxl)
library(plyr)
library(dplyr)
library(readr)
library(megaSDM)

#####  process NES data
nes.coords <- list.files(path = 'occs_compiled/NES_raw', pattern = 'xlsx', full.names = T) %>%
  lapply(read_xlsx) %>%
  plyr::rbind.fill() %>% 
  dplyr::select(3,5,6,8,9)

nes.coords$source <- 'NES'
colnames(nes.coords) = c('scientific_name', 'lat', 'long', 'year', 'region', 'source')
head(nes.coords)

# check species names
unique(nes.coords$scientific_name)

# recode some names
nes.coords$scientific_name <- dplyr::recode(nes.coords$scientific_name, 
                                            'Hyla japonica' = 'Dryophytes japonicus',
                                            'Rhabdophis lateralis' = 'Rahbdophis tigrinus',
                                            'Glandirana rugosa' = 'Glandirana emeljanovi',
                                            'Amphiesma vibakari' = 'Hebius vibakari',
                                            'Dinodon rufozonatum' = 'Lycodon rufozonatus',
                                            'Gloydius brevicaudus' = 'Gloydius brevicauda',
                                            'Gloydius saxatilis' = 'Gloydius intermedius',
                                            'Rana dybowskii' = 'Rana uenoi',
                                            'Hierophis spinalis' = 'Orientocoluber spinalis',
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


#####  Get GBIF data via the megaSDM package

# get species list == 240707 22:05 PM
list_sp_gbif <- c('Bombina orientalis',
                  'Bufo gargarizans',
                  'Bufo stejnegeri',
                  'Dryophytes japonicus',
                  'Dryophytes suweonensis',
                  'Dryophytes flaviventris',
                  'Glandirana emeljanovi',
                  'Hynobius leechii',
                  'Hynobius quelpaertensis',
                  'Hynobius yangi',
                  'Hynobius unisacculus',
                  'Hynobius notialis', 
                  'Hynobius geojeensis',
                  'Hynobius perplicatus',
                  'Karsenia koreana',
                  'Kaloula borealis',
                  'Lithobates catesbeianus',
                  'Onychodactylus koreanus',
                  'Onychodactylus sillanus',
                  'Pelophylax nigromaculatus',
                  'Pelophylax chosenicus',
                  'Rana coreana',
                  'Rana uenoi',
                  'Rana huanrenensis',
                  'Elaphe schrenckii',
                  'Elaphe anomala',
                  'Elaphe dione',
                  'Oocatochus rufodorsatus',
                  'Rhabdophis tigrinus',
                  'Hebius vibakari',
                  'Orientocoluber spinalis',
                  'Lycodon rufozonatus',
                  'Sibynophis chinensis',
                  'Gloydius brevicauda',
                  'Gloydius intermedius',
                  'Gloydius ussuriensis',
                  'Takydromus amurensis',
                  'Takydromus wolteri',
                  'Eremias argus',
                  'Scincella huanrenensis',
                  'Scincella vandenburghi',
                  'Gekko japonicus',
                  'Pelodiscus sinensis',
                  'Pelodiscus maackii',
                  'Mauremys sinensis')

# collect occurrence points
OccurrenceCollection(spplist = list_sp_gbif, 
                     output = 'occs_compiled/GBIF_raw',
                     trainingarea = c(125.0818, 130.9404, 33.11208, 38.61215))

# sort out coordinates
gbif.coords <- list.files(path = 'occs_compiled/GBIF_raw', pattern = '.csv', full.names = T) %>%
  lapply(read_csv) %>%
  plyr::rbind.fill() %>%
  dplyr::select(4,5,6)

gbif.coords$source = 'GBIF'
colnames(gbif.coords) = c('scientific_name', 'lat', 'long', 'source')
head(gbif.coords)

# export per species
for (i in 1:length(list_sp_gbif[-c(6,19)])) {
  per_sp <- gbif.coords %>% dplyr::filter(scientific_name == list_sp_gbif[-c(6,19)][[i]])
  file_name <- paste0('occs_compiled/GBIF_per_sp/', list_sp_gbif[-c(6,19)][[i]], '.csv')
  write.csv(per_sp, file_name)
}
