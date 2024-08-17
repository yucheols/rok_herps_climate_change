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
#dir.create('all_amphibians/output_other')
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
write.csv(get_contrib, 'all_amphibians/output_other/current_envs_raster_PC_var_contrib.csv')

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

# make into a list
unique(occs$scientific_name)

occs.list <- list(occs %>% filter(scientific_name == unique(scientific_name)[1]) %>% select(3,2),
                  occs %>% filter(scientific_name == unique(scientific_name)[2]) %>% select(3,2),
                  occs %>% filter(scientific_name == unique(scientific_name)[3]) %>% select(3,2),
                  occs %>% filter(scientific_name == unique(scientific_name)[4]) %>% select(3,2),
                  occs %>% filter(scientific_name == unique(scientific_name)[5]) %>% select(3,2),
                  occs %>% filter(scientific_name == unique(scientific_name)[6]) %>% select(3,2),
                  occs %>% filter(scientific_name == unique(scientific_name)[7]) %>% select(3,2),
                  occs %>% filter(scientific_name == unique(scientific_name)[8]) %>% select(3,2),
                  occs %>% filter(scientific_name == unique(scientific_name)[9]) %>% select(3,2),
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


#####  part 3 ::: sample background points ----------
bg <- randomPoints(mask = envs[[1]], n = 10000) %>% as.data.frame()
colnames(bg) = c('long', 'lat')
head(bg)

write.csv(bg, 'all_amphibians/bg/bg.csv')


#####  part 4 ::: test candidate models per species ---------- 
test.mod <- test_multisp(taxon.list = unique(occs$scientific_name),
                         occs.list = occs.list.thin,
                         envs = envs,
                         bg = bg,
                         tune.args = list(fc = c('L', 'LQ', 'H', 'LQH', 'LQHP', 'LQHPT'),
                                          rm = seq(1, 4, by = 0.5)),
                         partitions = 'block',
                         partition.settings = list(orientation = 'lat_lon'),
                         type = 'type1')
