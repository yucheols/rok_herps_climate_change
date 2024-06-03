### test workflow == probably need to update the ENMwrap package at some point to remove raster and rgdal dependencies.....
library(ENMwrap)
library(megaSDM)
library(raster)
library(dismo)
library(plyr)
library(dplyr)
library(readr)

# clear working environment
rm(list = ls(all.names = T))
gc()

# set random seed for reproducibility of random sampling elements
set.seed(333)

# prevent encoding error
Sys.getlocale()
Sys.setlocale("LC_CTYPE", ".1251")
Sys.getlocale()

##### part 1 ::: get climate data ---------------------------------------------------

# mask polygon == Republic of Korea
poly <- terra::vect('poly/KOR_adm0.shp')

# load climate data
envs <- terra::rast(list.files(path = 'E:/env layers/worldclim', pattern = '.tif$', full.names = T))
names(envs) = c('bio1', 'bio10', 'bio11', 'bio12', 'bio13', 'bio14', 'bio15', 'bio16', 'bio17',
                'bio18', 'bio19', 'bio2', 'bio3', 'bio4', 'bio5', 'bio6', 'bio7', 'bio8', 'bio9')

# layer masking and processing
envs <- terra::crop(envs, terra::ext(poly))
envs <- terra::mask(envs, poly) %>% raster::stack()
raster::plot(envs[[1]])

# export processed
for (i in 1:nlayers(envs)) {
  r <- envs[[i]]
  name <- paste0('envs/processed/', names(envs)[i], '.bil')
  writeRaster(r, filename = name, overwrite = T)
}


##### part 2 ::: get occurrence data ---------------------------------------------------

# make a list of species for draft modeling
# note B.gargarizans is used for B.sachalinensis and, D.immaculatus is used for D.suweonensis....just to follow the nomenclature recognized by the package
spplist <- c('Bombina orientalis',
             'Bufo stejnegeri',
             'Bufo gargarizans',
             'Dryophytes japonicus',
             'Dryophytes immaculatus',
             'Glandirana emeljanovi',
             'Hynobius leechii',
             'Hynobius quelpaertensis',
             'Hynobius yangi',
             'Kaloula borealis',
             'Karsenia koreana',
             'Onychodactylus koreanus',
             'Pelophylax chosenicus',
             'Pelophylax nigromaculatus',
             'Rana coreana',
             'Rana huanrenensis',
             'Rana uenoi')

# get data
#OccurrenceCollection(spplist = spplist,
#                     output = 'occs_test_workflow/raw',
#                     trainingarea = envs[[1]])

# compile raw data
occs_all <- list.files(path = 'occs_test_workflow/raw', pattern = '.csv', full.names = T) %>%
  lapply(read_csv) %>%
  plyr::rbind.fill() %>%
  dplyr::select('species', 'decimalLongitude', 'decimalLatitude')

colnames(occs_all) = c('species', 'long', 'lat')

# sort compiled data into two column (long, lat) dataframe by species
occs_list <- list(occs_all %>% filter(species == spplist[[1]]) %>% select(2,3),
                  occs_all %>% filter(species == spplist[[2]]) %>% select(2,3),
                  occs_all %>% filter(species == spplist[[3]]) %>% select(2,3),
                  occs_all %>% filter(species == spplist[[4]]) %>% select(2,3),
                  occs_all %>% filter(species == spplist[[5]]) %>% select(2,3),
                  occs_all %>% filter(species == spplist[[6]]) %>% select(2,3),
                  occs_all %>% filter(species == spplist[[7]]) %>% select(2,3),
                  occs_all %>% filter(species == spplist[[8]]) %>% select(2,3),
                  occs_all %>% filter(species == spplist[[9]]) %>% select(2,3),
                  occs_all %>% filter(species == spplist[[10]]) %>% select(2,3),
                  occs_all %>% filter(species == spplist[[11]]) %>% select(2,3),
                  occs_all %>% filter(species == spplist[[12]]) %>% select(2,3),
                  occs_all %>% filter(species == spplist[[13]]) %>% select(2,3),
                  occs_all %>% filter(species == spplist[[14]]) %>% select(2,3),
                  occs_all %>% filter(species == spplist[[15]]) %>% select(2,3),
                  occs_all %>% filter(species == spplist[[16]]) %>% select(2,3),
                  occs_all %>% filter(species == spplist[[17]]) %>% select(2,3))

# thin
thin <- occs_thinner(occs_list = occs_list, envs = envs[[1]], long = 'long', lat = 'lat', spp_list = spplist)

# export thinned occurrence data per species
for (i in 1:length(thin)) {
  file <- thin[[i]]
  write.csv(thin[[i]], paste0('occs_test_workflow/thinned/', spplist[[i]], '.csv'))
}


##### part 3 ::: get background data ---------------------------------------------------
# get random background points
bg <- randomPoints(mask = envs[[1]], n = 10000) %>% as.data.frame()
colnames(bg) = c('long', 'lat')
head(bg)


##### part 4 ::: select environmental variables ---------------------------------------------------
# may need to revise this part later on to better reflect species-specific climatic requirements....but this should be enough for now to test the workflow

# extract pixel values
vals <- raster::extract(envs, bg) %>% as.data.frame()
head(vals)

# generate correlation matrix
cormat <- cor(vals, method = 'pearson')
print(cormat)

# run correlation test
testcor <- caret::findCorrelation(cormat, cutoff = abs(0.7))
print(testcor)

# reduce the env dataset
envs.subs <- raster::dropLayer(envs, testcor) 
print(envs.subs)


##### part 5 ::: test candidate models per species  ---------------------------------------------------
testsp <- test_multisp(taxon.list = spplist,
                       occs.list = occs_list,
                       envs = envs.subs,
                       bg = bg,
                       tune.args = list(fc = c('LQ','LQHP'), rm = c(1, 1.5)),
                       partitions = 'block',
                       partition.settings = list(orientation = 'lat_lon'),
                       type = 'type1')

