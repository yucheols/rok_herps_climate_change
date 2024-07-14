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

# export loading contributions
write.csv(get_contrib, 'output_other/current_envs_raster_PC_var_contrib.csv')

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

# export
write.csv(bg, 'bg/bg.csv')


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
write.csv(test.mod$contrib[[1]], 'output_other/amp_current_ENMs_varimp.csv')

# G. brevicauda 
print(test.mod$contrib[[2]])
write.csv(test.mod$contrib[[2]], 'output_other/rep_current_ENMs_varimp.csv')


##### part 6 ::: response curves  ----------
### B.stejnegeri
# get plot data
amp.resp.data <- respDataPull(sp.name = 'B.stejnegeri', model = test.mod$models[[1]], names.var = names(envs))
head(amp.resp.data)

# plot
amp.resp.plot <- plot_response(resp.data = amp.resp.data)
amp.resp.plot +
  theme(axis.title = element_text(size = 14, face = 'bold'),
        axis.title.x = element_text(margin = margin(t = 20)),
        axis.title.y = element_text(margin = margin(r = 20)),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12),
        legend.title = element_text(size = 14, face = 'bold'),
        legend.text = element_text(size = 14, face = 'italic'),
        legend.position = 'top')

# save response plot
ggsave('plots/example_sp_ENMs/amp_resp.png', width = 30, height = 22, dpi = 800, units = 'cm')

### G.brevicauda
# get plot data
rep.resp.data <- respDataPull(sp.name = 'G.brevicauda', model = test.mod$models[[2]], names.var = names(envs))
head(rep.resp.data) 

# plot
rep.resp.plot <- plot_response(resp.data = rep.resp.data)
rep.resp.plot +
  theme(axis.title = element_text(size = 14, face = 'bold'),
        axis.title.x = element_text(margin = margin(t = 20)),
        axis.title.y = element_text(margin = margin(r = 20)),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12),
        legend.title = element_text(size = 14, face = 'bold'),
        legend.text = element_text(size = 14, face = 'italic'),
        legend.position = 'top')

# save response plot
ggsave('plots/example_sp_ENMs/rep_resp.png', width = 30, height = 22, dpi = 800, units = 'cm')

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

# export loading contributions
write.csv(get_contrib_fut, 'output_other/future_envs_raster_PC_var_contrib.csv')

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

##### part 9 ::: plot current & future side by side  ----------
# polygon
rok <- rgdal::readOGR('poly/KOR_adm1.shp')

### plot amphibian preds
amp.preds <- raster::stack(test.mod$preds[[1]], amp.fut.pred)
names(amp.preds) = c('Current', 'ssp585_2090')

amp.preds.plot <- plot_preds(preds = amp.preds, poly = rok, pred.names = names(amp.preds))
print(amp.preds.plot)

# save plot
ggsave('plots/example_sp_ENMs/amp_pred.png', width = 30, height = 14, dpi = 800, units = 'cm')

### plot reptile preds
rep.preds <- raster::stack(test.mod$preds[[2]], rep.fut.pred)
names(rep.preds) = c('Current', 'ssp585_2090')

rep.preds.plot <- plot_preds(preds = rep.preds, poly = rok, pred.names = names(rep.preds))
print(rep.preds.plot)

# save plot
ggsave('plots/example_sp_ENMs/rep_pred.png', width = 30, height = 14, dpi = 800, units = 'cm')


##### part 10 ::: threshold calculation, binary maps, area calculations  ----------
##### function to calculate MaxEnt thresholds == from https://babichmorrowc.github.io/post/2019-04-12-sdm-threshold/
sdm_threshold <- function(sdm, occs, type = "mtp", binary = FALSE){
  occPredVals <- raster::extract(sdm, occs)
  if(type == "mtp"){
    thresh <- min(na.omit(occPredVals))
  } else if(type == "p10"){
    if(length(occPredVals) < 10){
      p10 <- floor(length(occPredVals) * 0.9)
    } else {
      p10 <- ceiling(length(occPredVals) * 0.9)
    }
    thresh <- rev(sort(occPredVals))[p10]
  }
  sdm_thresh <- sdm
  sdm_thresh[sdm_thresh < thresh] <- NA
  if(binary){
    sdm_thresh[sdm_thresh >= thresh] <- 1
  }
  return(sdm_thresh)
}

##### B. stejnegeri
# calculate threshold
amp.th <- sdm_threshold(sdm = amp.preds[[1]], occs = thin[[1]], type = 'p10', binary = F)
print(minValue(amp.th))

# get binary map
amp.bins <- bin_maker(preds = amp.preds, th = rep(minValue(amp.th), nlayers(amp.preds)))
plot(amp.bins)

# plot
amp.bins.plot <- plot_preds(preds = amp.bins, poly = rok, colors = rev(terrain.colors(1000)), pred.names = names(amp.bins))
amp.bins.plot + theme(legend.position = 'none')

# save plot
ggsave('plots/example_sp_ENMs/amp_bin.png', width = 30, height = 14, dpi = 800, units = 'cm')

# calculate area
calc_ranges(bin.stack = amp.bins, bin.labs = c('current', 'future'), digits = 0)


##### G. brevicauda
# calculate threshold
rep.th <- sdm_threshold(sdm = rep.preds[[1]], occs = thin[[2]], type = 'p10', binary = F)
print(minValue(rep.th))

# get binary map
rep.bins <- bin_maker(preds = rep.preds, th = rep(minValue(rep.th), nlayers(rep.preds)))
plot(rep.bins)

# plot
rep.bins.plot <- plot_preds(preds = rep.bins, poly = rok, colors = rev(terrain.colors(1000)), pred.names = names(rep.bins))
rep.bins.plot + theme(legend.position = 'none')

# save plot
ggsave('plots/example_sp_ENMs/rep_bin.png', width = 30, height = 14, dpi = 800, units = 'cm')

# calculate area
calc_ranges(bin.stack = rep.bins, bin.labs = c('current', 'future'), digits = 0)


##### part 11 ::: export rasters  ----------

### B.stejnegeri
# continuous
for (i in 1:nlayers(amp.preds)) {
  r <- amp.preds[[i]]
  file_name <- paste0('output_rasters/example_sp/B.stejnegeri/cont_', names(amp.preds)[i], '.tif')
  writeRaster(r, filename = file_name, overwrite = T)
}

# binary
for (i in 1:nlayers(amp.bins)) {
  r <- amp.bins[[i]]
  file_name <- paste0('output_rasters/example_sp/B.stejnegeri/bin_', names(amp.bins)[i], '.tif')
  writeRaster(r, filename = file_name, overwrite = T)
}

### G.brevicauda
# continuous
for (i in 1:nlayers(rep.preds)) {
  r <- rep.preds[[i]]
  file_name <- paste0('output_rasters/example_sp/G.brevicauda/cont_', names(rep.preds)[i], '.tif')
  writeRaster(r, filename = file_name, overwrite = T)
}

# binary
for (i in 1:nlayers(rep.bins)) {
  r <- rep.bins[[i]]
  file_name <- paste0('output_rasters/example_sp/G.brevicauda/bin_', names(rep.bins)[i], '.tif')
  writeRaster(r, filename = file_name, overwrite = T)
}
