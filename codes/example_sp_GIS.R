##### Pick one reptile and one amphibian species to test the ENM workflows...for reptile, use G. brevicauda and for amphibian use Bufo stejnegeri
# extract GIS data

# clear working environment
rm(list = ls(all.names = T))
gc()

# turn off scientific notation
options(scipen = 999)

# load packages
library(raster)
library(dplyr)

### current annual mean temperature
cur.amp <- raster('envs/processed_full/bio1.tif')

### current annual precipitation
cur.mat <- raster('envs/processed_full/bio12.tif')

### future annual mean temperature
fut.amp <- raster('future_ssp/processed/bio1.bil')

### future annual precipitation
fut.mat <- raster('future_ssp/processed/bio12.bil')

### elevation
elev <- raster('envs/processed_full/elevation.tif')


##### part 1 ::: Bufo stejnegeri  ----------
# extract current "presence" coordinates from the binary map
amp.cur.bin <- raster('output_rasters/example_sp/B.stejnegeri/bin_Current.tif')
amp.cur.bin[amp.cur.bin < 1] <- NA
plot(amp.cur.bin)

amp.cur.locs <- raster::xyFromCell(amp.cur.bin, cell = 1:ncell(amp.cur.bin)) %>% as.data.frame()
colnames(amp.cur.locs) = c('long', 'lat')
head(amp.cur.locs)

# extract "presence" coordinates from the binary map
amp.fut.bin <- raster('output_rasters/example_sp/B.stejnegeri/bin_ssp585_2090.tif')
amp.fut.bin[amp.fut.bin < 1] <- NA
plot(amp.fut.bin)

amp.fut.locs <- raster::xyFromCell(amp.fut.bin, cell = 1:ncell(amp.fut.bin)) %>% as.data.frame()
colnames(amp.fut.locs) = c('long', 'lat')
head(amp.fut.locs)

# current annual mean temp
amp.cur.amp <- raster::extract(cur.amp, amp.cur.locs) %>% as.data.frame() %>% na.omit()
sd(amp.cur.amp[[1]])

# future annual mean temp
amp.fut.amp <- raster::extract(fut.amp, amp.fut.locs) %>% as.data.frame() %>% na.omit()
sd(amp.fut.amp[[1]])

# current annual precipitation
amp.cur.mat <- raster::extract(cur.mat, amp.cur.locs) %>% as.data.frame() %>% na.omit()
sd(amp.cur.mat[[1]])

# future annual precipitation
amp.fut.mat <- raster::extract(fut.mat, amp.fut.locs) %>% as.data.frame() %>% na.omit()
sd(amp.fut.mat[[1]])

# current elevation
amp.cur.elev <- raster::extract(elev, amp.cur.locs) %>% as.data.frame() %>% na.omit()
range(amp.cur.elev)
sd(amp.cur.elev[[1]])

# future elevation
amp.fut.elev <- raster::extract(elev, amp.fut.locs) %>% as.data.frame() %>% na.omit()
range(amp.fut.elev)
sd(amp.fut.elev[[1]])


##### part 2 ::: Gloydius brevicauda  ----------
# extract current "presence" coordinates from the binary map
rep.cur.bin <- raster('output_rasters/example_sp/G.brevicauda/bin_Current.tif')
rep.cur.bin[rep.cur.bin < 1] <- NA
plot(rep.cur.bin)

rep.cur.locs <- raster::xyFromCell(rep.cur.bin, cell = 1:ncell(rep.cur.bin)) %>% as.data.frame()
colnames(rep.cur.locs) = c('long', 'lat')
head(rep.cur.locs)

# extract "presence" coordinates from the binary map
rep.fut.bin <- raster('output_rasters/example_sp/G.brevicauda/bin_ssp585_2090.tif')
rep.fut.bin[rep.fut.bin < 1] <- NA
plot(rep.fut.bin)

rep.fut.locs <- raster::xyFromCell(rep.fut.bin, cell = 1:ncell(rep.fut.bin)) %>% as.data.frame()
colnames(rep.fut.locs) = c('long', 'lat')
head(rep.fut.locs)

# current annual mean temp
rep.cur.amp <- raster::extract(cur.amp, rep.cur.locs) %>% as.data.frame() %>% na.omit()
sd(rep.cur.amp[[1]])

# future annual mean temp
rep.fut.amp <- raster::extract(fut.amp, rep.fut.locs) %>% as.data.frame() %>% na.omit()
sd(rep.fut.amp[[1]])

# current annual precipitation
rep.cur.mat <- raster::extract(cur.mat, rep.cur.locs) %>% as.data.frame() %>% na.omit()
sd(rep.cur.mat[[1]])

# future annual precipitation
rep.fut.mat <- raster::extract(fut.mat, rep.fut.locs) %>% as.data.frame() %>% na.omit()
sd(rep.fut.mat[[1]])

# current elevation
rep.cur.elev <- raster::extract(elev, rep.cur.locs) %>% as.data.frame() %>% na.omit()
range(rep.cur.elev)
sd(rep.cur.elev[[1]])

# future elevation
rep.fut.elev <- raster::extract(elev, rep.fut.locs) %>% as.data.frame() %>% na.omit()
range(rep.fut.elev)
sd(rep.fut.elev[[1]])
