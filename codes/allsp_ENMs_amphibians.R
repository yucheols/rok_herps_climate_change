##### run ENMs for all Korean amphibian species

# clear working environment 
rm(list = ls(all.names = T))
gc()

# set random seed
set.seed(777)

# turn off scientific notation
options(scipen = 999)

# create directory to store outputs for this test
#dir.create('all_amphibians')
#dir.create('all_amphibians/occs_processed')
#dir.create('all_amphibians/output')
#dir.create('all_amphibians/output/contrib')
#dir.create('all_amphibians/output/models')
#dir.create('all_amphibians/output/preds_current')
#dir.create('all_amphibians/output/preds_future')
#dir.create('all_amphibians/output/metrics')
#dir.create('all_amphibians/output/other')
#dir.create('all_amphibians/envs')
#dir.create('all_amphibians/envs/PC_rasters')
#dir.create('all_amphibians/bg')

# load packages
library(ENMwrap)
library(ENMTools)
library(dismo)
library(factoextra)
library(raster)
library(plyr)
library(dplyr)
library(readr)


#####  part 1 ::: prep environmental data ----------
# use climate and topography
envs <- raster::stack(list.files(path = 'envs/processed_full', pattern = '.tif$', full.names = T))
envs <- raster::dropLayer(envs, c('agriculture', 'forest'))
print(envs)

# conduct raster PCA
envs.pca <- raster.pca(env = envs, n = 10)
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
write.csv(get_contrib, 'all_amphibians/output/other/current_envs_raster_PC_var_contrib.csv')

# set new env object
envs <- envs.pca$rasters[[1:5]]
print(envs)
plot(envs[[1]])

# export PC rasters
for (i in 1:nlayers(envs)) {
  r <- envs[[i]]
  file_name <- paste0('all_amphibians/envs/PC_rasters/', names(envs)[i], '.bil')
  writeRaster(r, filename = file_name, overwrite = T)
}


#####  part 2 ::: load occurrence points ----------
# GBIF data
gbif <- list.files(path = 'occs_compiled/GBIF_per_sp/amphibian', pattern = '.csv', full.names = T) %>%
  lapply(read_csv) %>%
  rbind.fill %>%
  select(2,3,4,5)

head(gbif)

# NES data
nes <- list.files(path = 'occs_compiled/NES_per_sp/amphibian', pattern = '.csv', full.names = T) %>%
  lapply(read_csv) %>%
  rbind.fill %>%
  select(2,3,4,7)

head(nes)

# merge data
occs <- rbind(gbif[, c(1,2,3)], nes[, c(1,2,3)])
occs$lat <- as.numeric(occs$lat)
occs$long <- as.numeric(occs$long)

# make into a list == exclude some species == e.g. new Hynobius species
unique(occs$scientific_name)

occs.list <- list(occs %>% filter(scientific_name == unique(scientific_name)[1]) %>% select(3,2),
                  occs %>% filter(scientific_name == unique(scientific_name)[2]) %>% select(3,2),
                  occs %>% filter(scientific_name == unique(scientific_name)[3]) %>% select(3,2),
                  occs %>% filter(scientific_name == unique(scientific_name)[4]) %>% select(3,2),
                  occs %>% filter(scientific_name == unique(scientific_name)[5]) %>% select(3,2),
                  occs %>% filter(scientific_name == unique(scientific_name)[7]) %>% select(3,2),
                  occs %>% filter(scientific_name == unique(scientific_name)[10]) %>% select(3,2),
                  occs %>% filter(scientific_name == unique(scientific_name)[11]) %>% select(3,2),
                  occs %>% filter(scientific_name == unique(scientific_name)[12]) %>% select(3,2),
                  occs %>% filter(scientific_name == unique(scientific_name)[13]) %>% select(3,2),
                  occs %>% filter(scientific_name == unique(scientific_name)[14]) %>% select(3,2),
                  occs %>% filter(scientific_name == unique(scientific_name)[15]) %>% select(3,2),
                  occs %>% filter(scientific_name == unique(scientific_name)[16]) %>% select(3,2),
                  occs %>% filter(scientific_name == unique(scientific_name)[17]) %>% select(3,2),
                  occs %>% filter(scientific_name == unique(scientific_name)[18]) %>% select(3,2),
                  occs %>% filter(scientific_name == unique(scientific_name)[19]) %>% select(3,2),
                  occs %>% filter(scientific_name == unique(scientific_name)[20]) %>% select(3,2),
                  occs %>% filter(scientific_name == unique(scientific_name)[21]) %>% select(3,2))

# thin occurrence points
occs.list.thin <- occs_thinner(occs_list = occs.list, envs = envs, long = 'long', lat = 'lat', spp_list = unique(occs$scientific_name))
glimpse(occs.list.thin)


# create a new taxon list since we dropped a few from the original list
taxon.list <- list(unique(occs$scientific_name)[1],
                   unique(occs$scientific_name)[2],
                   unique(occs$scientific_name)[3],
                   unique(occs$scientific_name)[4],
                   unique(occs$scientific_name)[5],
                   unique(occs$scientific_name)[7],
                   unique(occs$scientific_name)[10],
                   unique(occs$scientific_name)[11],
                   unique(occs$scientific_name)[12],
                   unique(occs$scientific_name)[13],
                   unique(occs$scientific_name)[14],
                   unique(occs$scientific_name)[15],
                   unique(occs$scientific_name)[16],
                   unique(occs$scientific_name)[17],
                   unique(occs$scientific_name)[18],
                   unique(occs$scientific_name)[19],
                   unique(occs$scientific_name)[20],
                   unique(occs$scientific_name)[21])


# export processed occurrence points
for (i in 1:length(occs.list.thin)) {
  o <- occs.list.thin[[i]]
  name <- paste0('all_amphibians/occs_processed/', taxon.list[[i]], '.csv')
  write.csv(o, name)
}

#####  part 3 ::: sample background points ----------
bg <- randomPoints(mask = envs[[1]], n = 10000) %>% as.data.frame()
colnames(bg) = c('long', 'lat')
head(bg)

write.csv(bg, 'all_amphibians/bg/bg.csv')


#####  part 4a ::: test candidate models per species ---------- 
test.mod <- test_multisp(taxon.list = taxon.list,
                         occs.list = occs.list.thin,
                         envs = envs,
                         bg = bg,
                         tune.args = list(fc = c('L', 'LQ', 'H', 'LQH', 'LQHP', 'LQHPT'),
                                          rm = seq(1, 4, by = 0.5)),
                         partitions = 'block',
                         partition.settings = list(orientation = 'lat_lon'),
                         type = 'type1')

# save the model object
#saveRDS(test.mod, 'all_amphibians/output/models/model_tuning.rds')

# import saved model object
#test.mod <- readRDS('all_amphibians/output/models/model_tuning.rds')

#####  part 4b ::: view results ---------- 
# print results
print(test.mod$metrics)
print(test.mod$models)
print(test.mod$preds)
print(test.mod$contrib)
print(test.mod$taxon.list)

# export metrics
write.csv(test.mod$metrics, 'all_amphibians/output/metrics/model_metrics.csv')

# export current prediction maps
for (i in 1:nlayers(test.mod$preds)) {
  r <- test.mod$preds[[i]]
  name <- paste0('all_amphibians/output/preds_current/', test.mod$taxon.list[[i]], '_preds.tif')
  writeRaster(r, filename = name, overwrite = T)
}

# export variable contribution 
for (i in 1:length(test.mod$contrib)) {
  c <- test.mod$contrib[[i]]
  name <- paste0('all_amphibians/output/contrib/', test.mod$taxon.list[[i]], '_contrib.csv')
  write.csv(c, name)
}


#####  part 5 ::: binary maps & ranges for current climate predictions ---------- 

# get thresholds
th <- get_thresh(preds = test.mod$preds, occs.list = occs.list.thin, type = 'p10')

# get binary maps
current.bin <- bin_maker(preds = test.mod$preds, th = th)
names(current.bin) = taxon.list

print(current.bin)
plot(current.bin)

# export current binary
#dir.create('all_amphibians/output/bins_current')

for (i in 1:nlayers(current.bin)) {
  r <- current.bin[[i]]
  name <- paste0('all_amphibians/output/bins_current/', names(current.bin)[[i]], '_bin.tif')
  writeRaster(r, filename = name, overwrite = T)
}


# calculate current ranges
current.range <- calc_ranges(bin.stack = current.bin, bin.labs = names(current.bin), digits = 0)
print(current.range)


#####  part 6 ::: climate change predictions // scenario == HadGEM3-GC31-LL_ssp585_2081-2100 ---------- 

# load future layers
envs.fut <- raster::stack(list.files(path = 'future_ssp/PC_rasters', pattern = '.bil$', full.names = T))
print(envs.fut)

# make future predictions
future.pred <- model_predictr(model = test.mod$models, preds.list = envs.fut, pred.names = taxon.list, method = 'multi2single')
print(future.pred)
plot(future.pred)

# export future predictions
for (i in 1:nlayers(future.pred)) {
  r <- future.pred[[i]]
  name <- paste0('all_amphibians/output/preds_future/', names(future.pred)[[i]], '_preds_future.tif')
  writeRaster(r, filename = name, overwrite = T)
}


#####  part 7 ::: binary maps & ranges for future climate predictions ---------- 

# get  binary maps
future.bin <- bin_maker(preds = future.pred, th = th)

print(future.bin)
plot(future.bin)

# export future binary
#dir.create('all_amphibians/output/bins_future')

for (i in 1:nlayers(future.bin)) {
  r <- future.bin[[i]]
  name <- paste0('all_amphibians/output/bins_future/', names(future.bin)[[i]], '_bin_future.tif')
  writeRaster(r, filename = name, overwrite = T)
}

# calculate future ranges
future.range <- calc_ranges(bin.stack = future.bin, bin.labs = names(future.bin), digits = 0)
print(future.range)


#####  part 8 ::: range calculation within National Parks ---------- 

# import Korea NP borders file
knp <- rgdal::readOGR('poly/KNP_boundaries/NLPRK_BNDRY.shp')

### current
# mask to NP borders
current.knp <- raster::mask(current.bin, knp)
plot(current.knp)

# calculate current ranges within NP
current.range.knp <- calc_ranges(bin.stack = current.knp, bin.labs = names(current.knp), digits = 0) 
print(current.range.knp)

### future
# mask to NP borders
future.knp <- raster::mask(future.bin, knp)
plot(future.knp)

# calculate future ranges within NP
future.range.knp <- calc_ranges(bin.stack = future.knp, bin.labs = names(future.knp), digits = 0)
print(future.range.knp)


#####  part 9 ::: GIS calculations ---------- 
# current annual mean temperature
cur.amt <- raster('envs/processed_full/bio1.tif')

# current annual precipitation
cur.map <- raster('envs/processed_full/bio12.tif')

# future annual mean temperature
fut.amt <- raster('future_ssp/processed/bio1.bil')

# future annual precipitation
fut.map <- raster('future_ssp/processed/bio12.bil')

# elevation
elev <- raster('envs/processed_full/elevation.tif')

# function to convert binary suitable pixels to coordinates
bin2coords <- function(bin) {
  output <- list()
  
  for (i in 1:nlayers(bin)) {
    r <- bin[[i]]
    locs.which <- which(r[] == 1)
    locs <- raster::xyFromCell(r, locs.which) %>% as.data.frame() %>% na.omit()
    colnames(locs) = c('long', 'lat')
    output[[i]] <- locs
  }
  return(output)
}

# convert current binary suitable pixels to coords
current.xy <- bin2coords(bin = current.bin)
print(current.xy)

# convert future binary suitable pixels to coords
future.xy <- bin2coords(bin = future.bin)
print(future.xy)

# function for batch extraction of environmental values
extractr <- function(env.var, coords, var.name) {
  output <- list()
  
  for (i in 1:length(coords)) {
    ext <- raster::extract(env.var, coords[[i]]) %>% as.data.frame() %>% na.omit()
    colnames(ext) = var.name
    output[[i]] <- ext
  }
  return(output)
}

##### extract environmental values

### 1. current annual mean temp
current.amt <- extractr(env.var = cur.amt, coords = current.xy, var.name = 'AMT_current') 
print(current.amt)
head(current.amt[[1]])

# check mean values
for (i in 1:length(current.amt)) {
  print(mean(current.amt[[i]]$AMT_current))
}

### 2. current annual precipitation
current.pr <- extractr(env.var = cur.map, coords = current.xy, var.name = 'MAP_current')
print(current.pr)
head(current.pr[[1]])

# check mean value
for (i in 1:length(current.pr)) {
  print(mean(current.pr[[i]]$MAP_current))
}


### 3. future annual mean temperature
future.amt <- extractr(env.var = fut.amt, coords = future.xy, var.name = 'AMT_future')
print(future.amt)
head(future.amt[[1]])

# check mean value
for (i in 1:length(future.amt)) {
  print(mean(future.amt[[i]]$AMT_future))
}


### 4. future annual precipitation
future.pr <- extractr(env.var = fut.map, coords = future.xy, var.name = 'MAP_future')
print(future.pr)
head(future.pr[[1]])

# check mean value
for (i in 1:length(future.pr)) {
  print(mean(future.pr[[i]]$MAP_future))
}


### 5. current elevation
current.elev <- extractr(env.var = elev, coords = current.xy, var.name = 'current_elev')
print(current.elev)
head(current.elev[[1]])

# check mean value
for (i in 1:length(current.elev)) {
  print(mean(current.elev[[i]]$current_elev))
}


### 6. future elevation
future.elev <- extractr(env.var = elev, coords = future.xy, var.name = 'future_elev')
print(future.elev)
head(future.elev[[1]])

# check mean value
for (i in 1:length(future.elev)) {
  print(mean(future.elev[[i]]$future_elev))
}

