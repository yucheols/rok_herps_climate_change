##### Pick one reptile and one amphibian species to test the ENM workflows...for reptile, use G. brevicauda and for amphibian use Bufo stejnegeri

# clear working environment
rm(list = ls(all.names = T))
gc()

# turn off scientific notation
options(scipen = 999)

# set random seed
set.seed(100)

# prevent encoding errors
Sys.getlocale()
Sys.setlocale("LC_CTYPE", ".1251")
Sys.getlocale()

# load packages
library(ENMwrap)
library(ENMTools)
library(dismo)
library(ecospat)
library(raster)
library(dplyr)
library(factoextra)


##### part 1 ::: environmental data ----------
# use climate and topography
envs <- raster::stack(list.files(path = 'envs/processed_full', pattern = '.tif$', full.names = T))
envs <- raster::dropLayer(envs, c('agriculture', 'forest'))
print(envs)

# conduct raster PCA
envs.pca <- ENMTools::raster.pca(env = envs, n = 10)
print(envs.pca$pca.object)

# get PCA eigenvalues == drop PC layers with eigen < 1.0 following the Kaiser's criterion
get_eigenvalue(envs.pca$pca.object)

# get loading contributions == contribution of each variable to PCs
get.res <- get_pca_var(envs.pca$pca.object)
get_contrib <- as.data.frame(get.res$contrib[, c(1:5)])
print(get_contrib)

rev(sort(get_contrib[[1]]))  # bio11 (9.71 %) ... bio6 (9.67 %) == coldest temp
rev(sort(get_contrib[[2]]))  # bio19 (11.02 %) ... bio17 (10.97 %) == precip of coldest/driest months
rev(sort(get_contrib[[3]]))  # bio16 (20.64 %) ... bio18 (17.66 %) == precip of warmest/wettest months
rev(sort(get_contrib[[4]]))  # bio14 (17.84 % ) ... bio2 (16.23 %) == precip of driest month & mean diurnal range
rev(sort(get_contrib[[5]]))  # bio3 (34.08 %) ... slope (16.62 %) == isothermality and slope

# set new env object
envs <- envs.pca$rasters[[1:5]]
print(envs)
plot(envs[[1]])


##### part 2 ::: occurrence data ----------
# NES data
amp_nes <- read.csv('occs_compiled/NES_per_sp/Bufo stejnegeri.csv') %>% dplyr::select(2,4,3)
head(amp_nes)

rep_nes <- read.csv('occs_compiled/NES_per_sp/Gloydius brevicauda.csv') %>% dplyr::select(2,4,3) 
head(rep_nes)

# GBIF data
amp_gbif <- read.csv('occs_compiled/GBIF_per_sp/Bufo stejnegeri.csv') %>% dplyr::select(2,4,3)
head(amp_gbif)

rep_gbif <- read.csv('occs_compiled/GBIF_per_sp/Gloydius brevicauda.csv') %>% dplyr::select(2,4,3)
head(rep_gbif)

# merge
amp <- rbind(amp_nes, amp_gbif)
rep <- rbind(rep_nes, rep_gbif)

# thin occurrence points
thin <- occs_thinner(occs_list = list(amp[, c(2,3)], rep[, c(2,3)]), envs = envs[[1]], long = 'long', lat = 'lat', 
                     spp_list = c('Bufo stejnegeri', 'Gloydius brevicauda'))

glimpse(thin)

#amp_thin <- thin[[1]]
#rep_thin <- thin[[2]]

#head(amp_thin)
#head(rep_thin)


##### part 3 ::: background data ----------
bg <- randomPoints(mask = envs[[1]], n = 10000) %>% as.data.frame()
colnames(bg) = c('long', 'lat')
head(bg)

##### part 4 ::: test candidate models per species  ----------
test.mod <- test_multisp(taxon.list = c('Bufo stejnegeri', 'Gloydius brevicauda'),
                         occs.list = thin,
                         envs = envs, 
                         bg = bg,
                         tune.args = list(fc = c('L', 'LQ', 'H', 'LQH', 'LQHP', 'LQHPT'),
                                          rm = seq(1, 4, by = 0.5)),
                         partitions = 'block',
                         partition.settings = list(orientation = 'lat_lon'),
                         type = 'type1')
