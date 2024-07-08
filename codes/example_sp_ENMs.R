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
rev(sort(get_contrib[[4]]))  # bio14 (17.84 %) ... bio2 (16.23 %) == precip of driest month & mean diurnal range
rev(sort(get_contrib[[5]]))  # bio3 (34.08 %) ... slope (16.62 %) == isothermality and slope

# set new env object
envs <- envs.pca$rasters[[1:5]]
print(envs)
plot(envs[[1]])

# export PC rasters
for (i in 1:nlayers(envs)) {
  r <- envs[[i]]
  file_name <- paste0('envs/PC_rasters/', names(envs)[i], '.bil')
  writeRaster(r, filename = file_name, overwrite = T)
}


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
# model tuning
test.mod <- test_multisp(taxon.list = c('Bufo stejnegeri', 'Gloydius brevicauda'),
                         occs.list = thin,
                         envs = envs, 
                         bg = bg,
                         tune.args = list(fc = c('L', 'LQ', 'H', 'LQH', 'LQHP', 'LQHPT'),
                                          rm = seq(1, 4, by = 0.5)),
                         partitions = 'block',
                         partition.settings = list(orientation = 'lat_lon'),
                         type = 'type1')

# check outputs
print(test.mod$metrics)
print(test.mod$models)
print(test.mod$preds)

# save output
saveRDS(test.mod, 'output_rds/example_sp_ENMs_20240708.rds')

# import saved model (shortcut)
#test.mod <- readRDS('output_rds/example_sp_ENMs_20240708.rds')


##### part 5 ::: variable importance  ----------
# B. stejnegeri
print(test.mod$contrib[[1]])

# G. brevicauda 
print(test.mod$contrib[[2]])


##### part 6 ::: response curves  ----------
## get plot data
# B.stejnegeri
amp.resp.data <- respDataPull(sp.name = 'B.stejnegeri', model = test.mod$models[[1]], names.var = names(envs))
head(amp.resp.data)

# G.brevicauda
rep.resp.data <- respDataPull(sp.name = 'G.brevicauda', model = test.mod$models[[2]], names.var = names(envs))
head(rep.resp.data)


##### part 7 ::: future layer prep == HadGEM3.GC31.LL_ssp585  ----------
# import future climate layers
hadgem <- raster::stack('future_ssp/wc2.1_30s_bioc_HadGEM3-GC31-LL_ssp585_2081-2100.tif')
names(hadgem) = c('bio1','bio2','bio3','bio4','bio5', 'bio6','bio7','bio8','bio9','bio10','bio11', 
                  'bio12','bio13', 'bio14','bio15','bio16','bio17','bio18','bio19')

print(hadgem)

# process future climate layers
hadgem <- raster::crop(hadgem, extent(envs))
hadgem <- raster::mask(hadgem, envs[[1]])
plot(hadgem[[1]])

# add topo layers
elev <- raster('envs/processed_full/elevation.tif')
slope <- raster('envs/processed_full/slope.tif')

fut.envs <- raster::stack(hadgem, elev, slope)

# export future layers
for (i in 1:nlayers(fut.envs)) {
  r <- fut.envs[[i]]
  file_name <- paste0('future_ssp/processed/', names(fut.envs)[i], '.bil')
  writeRaster(r, filename = file_name, overwrite = T)
}

# conduct raster PCA == get 5 raster layers regardless of PC eigenvalues...bcz nlayers should match between current and future layer sets
fut.envs.pca <- ENMTools::raster.pca(env = fut.envs, n = 5)
print(fut.envs.pca$pca.object)

# but lets at least look at variable contribution to each PC
get_res_fut <- get_pca_var(fut.envs.pca$pca.object)
get_contrib_fut <- as.data.frame(get_res_fut$contrib[, c(1:5)])
print(get_contrib_fut)

rev(sort(get_contrib_fut[[1]]))  # bio11 (10.55 %) ... bio6 (10.47 %) == coldest temp
rev(sort(get_contrib_fut[[2]]))  # bio5 (13.74 %) ... bio8 (12.97 %) == max temp of warmest month & mean temp of wettest months
rev(sort(get_contrib_fut[[3]]))  # bio18 (17.35 %) ... bio16 (16.70 %) == precip of warmest/wettest months
rev(sort(get_contrib_fut[[4]]))  # bio14 (33.56 %) ... bio17 (21.47 %) == precip of driest months
rev(sort(get_contrib_fut[[5]]))  # bio3 (38.83 %) ... slope (22.14 %) == isothermality and annual precip

# set new future env object
fut.envs <- fut.envs.pca$rasters
print(fut.envs)

# export PC rasters
for (i in 1:nlayers(fut.envs)) {
  r <- fut.envs[[i]]
  file_name <- paste0('future_ssp/PC_rasters/', names(fut.envs)[i], '.bil')
  writeRaster(r, filename = file_name, overwrite = T)
}


##### part 8 ::: future projections  ----------
# B. stejnegeri
amp.fut.pred <- dismo::predict(object = test.mod$models[[1]], x = fut.envs)
plot(amp.fut.pred)

# G. brevicauda
rep.fut.pred <- dismo::predict(object = test.mod$models[[2]], x = fut.envs)
plot(rep.fut.pred)
