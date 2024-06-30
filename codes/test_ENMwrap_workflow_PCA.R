### test workflow == probably need to update the ENMwrap package at some point to remove raster and rgdal dependencies.....
### use PCA instead of selecting variables from the correlation test

# clear working environment
rm(list = ls(all.names = T))
gc()

# turn off scientific notation
options(scipen = 999)

# load packages
library(ENMwrap)
library(ENMTools)
library(megaSDM)
library(raster)
library(dismo)
library(plyr)
library(dplyr)
library(readr)

# set random seed for reproducibility of random sampling elements
set.seed(333)

# prevent encoding error
Sys.getlocale()
Sys.setlocale("LC_CTYPE", ".1251")
Sys.getlocale()


##### part 1 ::: get environmental data ---------------------------------------------------
# mask polygon == Republic of Korea
poly <- terra::vect('poly/KOR_adm0.shp')

# load climate data
clim <- terra::rast(list.files(path = 'E:/env layers/worldclim', pattern = '.tif$', full.names = T))
names(clim) = c('bio1', 'bio10', 'bio11', 'bio12', 'bio13', 'bio14', 'bio15', 'bio16', 'bio17',
                'bio18', 'bio19', 'bio2', 'bio3', 'bio4', 'bio5', 'bio6', 'bio7', 'bio8', 'bio9')

clim <- terra::crop(clim, terra::ext(poly))
clim <- terra::mask(clim, poly)
plot(clim[[1]])

# topography
topo <- terra::rast(list.files(path = 'E:/env layers/topo_1km_earthenv', pattern = '.tif$', full.names = T))
names(topo) = c('elevation', 'slope')

topo <- terra::crop(topo, terra::ext(poly))
topo <- terra::mask(topo, poly)
plot(topo[[1]])

# land cover
land <- terra::rast(list.files(path = 'E:/env layers/land cover', pattern = '.tif$', full.names = T))
land <- subset(land, c(2,5))
names(land) = c('agriculture', 'forest')

land <- terra::crop(land, terra::ext(poly))
land <- terra::mask(land, poly)
plot(land[[1]])

# layer masking and processing
envs <- c(clim, topo, land)
plot(envs[[1]])

# export processed
for (i in 1:terra::nlyr(envs)) {
  r <- envs[[i]]
  terra::writeRaster(r, filename = paste0('envs/processed_full/', names(envs)[i], '.tif'), overwrite = T)
}


##### part 2 ::: conduct PCA on environmental variables ---------------------------------------------------
# run raster PCA
envs.pca <- raster.pca(raster::stack(envs), n = 7)

# get eigenvalues == Dim 1 ~ Dim 5 have eigenvalue > 1
factoextra::get_eigenvalue(envs.pca$pca.object)

# get 5 rasters to be used
envs <- raster::stack(subset(envs.pca$rasters, c(1:5)))
print(envs)

# check loadings


##### part 3 ::: get occurrence data ---------------------------------------------------
# make a list of species for draft modeling
# note B.gargarizans is used for B.sachalinensis and, D.immaculatus is used for D.suweonensis....just to follow the nomenclature recognized by the package
# name matching is required for data filtering and downstream model testing. If the names dont match the model testing step will throw errors
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


##### part 4 ::: get background data ---------------------------------------------------
# get random background points
bg <- randomPoints(mask = envs[[1]], n = 10000) %>% as.data.frame()
colnames(bg) = c('long', 'lat')
head(bg)


##### part 5 ::: test candidate models per species  ---------------------------------------------------
# fit models per species == may need to increase the feature class (fc) and regularization (rm) combinations for actual application
# for the code below (testing two features and two regularizations), it takes aproximately 2~3 minutes per species
testsp <- test_multisp(taxon.list = spplist,
                       occs.list = occs_list,
                       envs = envs,
                       bg = bg,
                       tune.args = list(fc = c('L','LQ','H','LQH','LQHP','LQHPT'), rm = seq(1,4, by = 0.5)),
                       partitions = 'block',
                       partition.settings = list(orientation = 'lat_lon'),
                       type = 'type1')

# check results
print(testsp$metrics)
print(testsp$models)
print(testsp$preds)
print(testsp$contrib)
print(testsp$taxon.list)

# plot predictions
# should fix the package code at some point to remove the rgdal and raster dependencies....so that I can directly use the SpatVector object loaded earlier....
# but I dont have enough time now so this will have to do
rok <- rgdal::readOGR('poly/KOR_adm0.shp')
plot_preds(preds = testsp$preds, poly = rok, pred.names = spplist)

# make predictions under current climate
