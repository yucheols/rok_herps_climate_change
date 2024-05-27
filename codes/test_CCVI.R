### clear working environment before starting
rm(list = ls(all.names = T))
gc()

###  follow the tutorial
library(ccviR)
library(sf)
library(dplyr)
library(terra)
library(ggplot2)

# launch the vignette                                
vignette('app_vignette', package = 'ccviR')           # using demo data
vignette('app_details_vignette', package = 'ccviR')   # using app in practice
vignette('data_prep_vignette', package = 'ccviR')     # prepare custom climate data
vignette('package_vignette', package = 'ccviR')       # directly in R

##### run ccvi app with demo data == facing some errors here....  --------------------------------------------------
run_ccvi_app('demo')


##### run ccvi directly in R  --------------------------------------------------

# load example data
data_pth <- system.file('extdata', package = 'ccviR')


##### prepare climate data  --------------------------------------------------
# get climate data file names
list.files(file.path(data_pth, 'clim_files/raw'))

# MAT: mean annual temperature for the historical normal period (required)
# MAT_2050: mean annual temperature for the future under climate change. It can be any number eg 2050, 2100 (required)
# CMD: climate moisture deficit for the historical normal period (required)
# CMD_2050: climate moisture deficit for the future under climate change it can be any number eg 2050, 2100 (required)
# CCEI: Climate Change Exposure Index from NatureServe website
# MAP: mean annual precipitation for the historical normal period
# MWMT: mean warmest month temperature for the historical normal period
# MCMT: mean coldest month temperature for the historical normal period
# clim_poly: An optional shapefile with a polygon of the extent of the climate data. 
#            It will be created from the climate data if it is missing but it is faster to provide it.

# prepare the data and save it in out_folder
# RCP4.5 prep - saves breaks as brks for use in 8.5 prep
brks <- prep_clim_data(
  mat_norm = file.path(data_pth, 'clim_files/raw', 'NB_norm_MAT.tif'),
  mat_fut = file.path(data_pth, 'clim_files/raw', 'NB_RCP.4.5_MAT.tif'),
  cmd_norm = file.path(data_pth, 'clim_files/raw', 'NB_norm_CMD.tif'),
  cmd_fut = file.path(data_pth, 'clim_files/raw', 'NB_RCP.4.5_CMD.tif'),
  map = file.path(data_pth, 'clim_files/raw', 'NB_norm_MAP.tif'),
  mwmt = file.path(data_pth, 'clim_files/raw', 'NB_norm_MWMT.tif'),
  mcmt = file.path(data_pth, 'clim_files/raw', 'NB_norm_MCMT.tif'),
  clim_poly = file.path(data_pth, 'clim_files/processed', 'clim_poly.shp'),
  out_folder = file.path(data_pth, 'clim_files/processed'),
  overwrite = T,
  scenario_name = 'RCP 4.5'
  )

# RCP8.5 - using breaks from 4.5
# map, mwmt and mcmt only need to be processed once
prep_clim_data(
  mat_norm = file.path(data_pth, 'clim_files/raw', 'NB_norm_MAT.tif'),
  mat_fut = file.path(data_pth, 'clim_files/raw', 'NB_RCP.8.5_MAT.tif'),
  cmd_norm = file.path(data_pth, 'clim_files/raw', 'NB_norm_CMD.tif'),
  cmd_fut = file.path(data_pth, 'clim_files/raw', 'NB_RCP.8.5_CMD.tif'),
  out_folder = file.path(data_pth, 'clim_files/processed'),
  clim_poly = file.path(system.file('extdata', package = 'ccviR'),
                        'assess_poly.shp'),
  overwrite = T,
  scenario_name = 'RCP 8.5',
  brks_mat = brks$brks_mat,
  brks_cmd = brks$brks_cmd,
  brks_ccei = brks$brks_ccei
)


##### explore the prepared data  --------------------------------------------------
# select color ramp for plots
pal = c('#FFF9CA', '#FEE697', '#FEC24D', '#F88B22', '#D85A09', '#A33803')

# create rasters for the historical normal MAT and the future MAT
mat_norm <- rast(file.path(data_pth, 'clim_files/raw/NB_norm_MAT.tif'))
mat_fut <- rast(file.path(data_pth, 'clim_files/raw/NB_RCP.4.5_MAT.tif'))

# plot MAT change
plot(mat_norm - mat_fut,
     col = colorRampPalette(rev(pal))(50))

# breaks used to classify MAT data (determined when preparing RCP4.5 data above)
brks$brks_mat

# create raster for the classified MAT. Note that the plot only shows five classes as none of the data falls within the sixth class.
mat_classified <- rast(file.path(data_pth, 'clim_files/processed/MAT_reclassRCP_4.5.tif'))

# plot the classified MAT
plot(mat_classified, col = pal)

# create a csv readme file to record where the data came from and any relevant metadata
write.csv(
  data.frame(Scenario_name = c('RCP 4.5', 'RCP 8.5'),
             GCM_or_Ensemble_name = c('AdaptWest 15 CMIP5 AOGCM Ensemble'),
             Historical_normal_period = '1961-1990',
             Future_period = '2050s',
             Emissions_scenario = c('RCP 4.5', 'RCP 8.5'),
             Link_to_source = 'https://adaptwest.databasin.org/pages/adaptwest-climatena-cmip5/',
             brks_mat = brks$brks_mat %>% brks_to_txt(),
             brks_cmd = brks$brks_cmd %>% brks_to_txt(),
             brks_ccei = brks$brks_ccei %>% brks_to_txt()),
  file.path(data_pth, 'clim_files/processed/', 'climate_data_readme.csv'),
  row.names = F
)


##### Load species specific spatial data  --------------------------------------------------
# The following spatial data sets can be input for each species: 
  
# Species North American or global range polygon (required)
# Assessment area polygon (required)
# Non-breeding range polygon
# Projected range change raster and a matrix to reclassify it as 0 = unsuitable, 1 = lost, 2 = maintained and 3 = gained
# Physiological thermal niche (PTN) polygon. PTN polygon should include cool or cold environments that the species occupies 
# that may be lost or reduced in the assessment area as a result of climate change.

# input the data
rng_poly <- read_sf(file.path(data_pth, 'rng_poly.shp'), agr = 'constant')                        # species range polygon
assess_poly <- read_sf(file.path(data_pth, 'assess_poly.shp'), agr = 'constant')                  # assessment area polygon
rng_chg <- rast(c(file.path(data_pth, 'rng_chg_45.tif'), file.path(data_pth, 'rng_chg_85.tif')))  # range change data
PTN_poly <- read_sf(file.path(data_pth, 'PTN_poly.shp'), agr = 'constant')                        # physiological thermal niche (PTN) polygon

# the range change raster has values from -1, 0 and 1, this matrix is used to
# convert them to the 4 classes described above.
hs_rcl_mat <- matrix(c(-1:1, c(1, 2 ,3)), ncol = 2)


##### Load climate data  --------------------------------------------------
clim_dat <- get_clim_vars(file.path(data_pth, 'clim_files/processed'),
                          scenario_names = c('RCP 4.5', 'RCP 8.5'))

str(clim_dat, max.level = 1)


##### Run spatial data analysis  --------------------------------------------------
spat_res <- analyze_spatial(range_poly = rng_poly, scale_poly = assess_poly,
                            clim_vars_lst = clim_dat, ptn_poly = PTN_poly,
                            hs_rast = rng_chg, hs_rcl = hs_rcl_mat,
                            scenario_names = c('RCP 4.5', 'RCP 8.5'))

spat_res$spat_table[, 1:7]

# create a map to visualize the exposure to changes in MAT for the RCP 4.5 scenario
plot(clim_dat$mat$RCP_4.5, col = pal)
plot(st_geometry(spat_res$range_poly_assess), add = T)


##### Answer the vulnerability questions  --------------------------------------------------
# create a blank vulnerability question table
vuln <- make_vuln_df('sp_name')

# you can interactively edit the table
vuln <- edit(vuln)

# or use code (recommended for reproducibility)
vuln$Value1[3:19] <- c(0, 0, 1, 0, -1, -1, -1, -1, 0, 0, 1, 0, 0, 1, 0, 0, 0)
vuln$Value1[26:19] <- c(0, -1, -1, 0)

# include a second value to reflect uncertainty and trigger a monte carlo to
# determine confidence
vuln$Value2[3:5] <- c(2, 0, 0)

# editing using a csv
# save the table as a csv 
#write.csv(vuln, "path/to/write", row.names = FALSE)

# edit the csv and save it 
# load it back in to R
#vuln <- read.csv("path/to/write", stringsAsFactors = FALSE)


##### Calculate the index  --------------------------------------------------
# run calculation
index_res <- calc_vulnerability(spat_df = spat_res$spat_table, vuln_df = vuln, tax_grp = 'Bird', n_rnds = 20)

# check results
glimpse(index_res)
index_res$index

# See how the index values 
# (that are calculated separately for the modeled response to climate change and sensitivity and adaptive capacity modified by exposure sections)
# are combined into the overall index value.

# set plot theme
my_theme <- theme_classic() +
  theme(text = element_text(size = 12),
        strip.background = element_blank())

theme_set(my_theme)

plot_score_index(index_res)

# look at scores for each factor
q_score <- bind_rows(index_res$vuln_df %>% `names<-`(index_res$scenario_name),
                     .id  = 'scenario_name') %>%
  plot_q_score()

print(q_score)

