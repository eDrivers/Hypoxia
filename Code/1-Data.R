# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                     LIBRARIES
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(sf)
library(sp)
library(magrittr)
library(tidyverse)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                   DOWNLOAD DATA
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The data used to characterize hypoxia comes from DFO's annual surveys.
# We also use the benthic habitats description from Jean-Denis Dutil to get
# bathymetry data to interpolate dissolved oxygen.
# For more information read the repo's README.md document.

# Output location for downloaded data
output <- './Data/RawData'

# Will have to upload the data on zenodo and eventually get the data from SLGO.
# For now, I'm using the data downloaded manually from the website.

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                   IMPORT DATA
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ---------------- #
# Dissolved Oxygen #
# ---------------- #
# Import, spatial object &
hypoxia <- read.csv('./Data/RawData/Btm_O2_GSL_2013_2017.txt', sep = '\t') %>%
           st_as_sf(coords = c('LOND','LATD'), crs = 4326) %>%
           st_transform(crs = 32198)

# ---------- #
# Bathymetry #
# ---------- #
# Unzip file
unzip(zipfile = paste0(output, '/megahab.zip'), exdir = output)

# File name
fileName <- dir(output, pattern = 'shp$', full.names = T)


# Import and select bathymetry
bathy <- st_read(fileName) %>%
         select(Bathy_Mean) %>%
         mutate(Bathy_Mean = abs(round(Bathy_Mean))) %>%
         rename(DEPH = Bathy_Mean)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                               FORMAT OXYGEN DATA
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Divide by years and remove duplicated coordinates
ys <- unique(hypoxia$YEAR) %>% sort()
hyp <- vector('list', length(ys))

for(i in 1:length(ys)) {
  # Year ID
  id <- hypoxia$YEAR == ys[i]

  # Select only year
  hyp[[i]] <- hypoxia[id, ]

  # Add coordinates
  hyp[[i]] <- cbind(hyp[[i]], st_coordinates(hyp[[i]]))

  # Sort by coordinates and depth
  hyp[[i]] <- hyp[[i]][order(hyp[[i]]$X, hyp[[i]]$Y, hyp[[i]]$DEPH), ]

  # Remove duplicated coordinates
  hyp[[i]] <- hyp[[i]][!duplicated(hyp[[i]][, c('X','Y'), drop = T]), ]

  # Remove zero distance points
  temp <- as(hyp[[i]], 'Spatial')
  temp <- zerodist(temp)[,1]
  if(length(temp) > 1) hyp[[i]] <- hyp[[i]][-temp,]
}

# Data point in 2014 messing up the interpolation, likely because it's at the maximum observed depth in the Dutil data
hyp[[2]]$DEPH[hyp[[2]]$DEPH == 521] <- 520

# Even though the interpolation works, the results are somehow very weird.
#I will remove 2014 from the data for now, but this should be further investigated
hyp[[2]] <- NULL
ys <- ys[-2]

# Name hypoxia list
names(hyp) <- ys

# Change name
hypoxia <- hyp


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                  EXPORT DATA
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Export object as .RData
save(hypoxia, file = './data/rawData/hypoxia.RData')
save(bathy, file = './data/rawData/bathy.RData')
