### test workflow
library(ENMwrap)
library(megaSDM)
library(raster)

# clear working environment
rm(list = ls(all.names = T))
gc()

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
spplist <- c('Bombina orientalis',
             'Bufo stejnegeri',
             'Bufo sachalinensis',
             'Dryophytes japonicus',
             'Dryophytes suweonensis',
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
OccurrenceCollection(spplist = spplist,
                     output = 'occs/raw',
                     trainingarea = envs[[1]])


##### part 3 ::: get background data ---------------------------------------------------
