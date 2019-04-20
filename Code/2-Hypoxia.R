# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                  LIBRARIES
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(sf)
library(sp)
library(magrittr)
library(tidyverse)
library(graphicsutils)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                    DATA
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load('./data/rawData/bathy.RData')
load('./data/rawData/hypoxia.RData')
ys <- names(hypoxia)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                CHECK OUTLIERS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Visualize outliers
graphicsutils::plot0(xlim = c(0,5), ylim = c(0, 150))
mtext('Outliers', 3, font = 2)
for(i in 1:length(ys)) boxplot(hypoxia[[i]][, 'sat', drop = T], add = T, at = i, frame = F, range = 3)
axis(1, at = 1:length(ys), labels = ys)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                 INTERPOLATION
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Get centroid of grid cells
coords <- st_centroid(bathy)

# Transform to sp object
coords <- as(coords, 'Spatial')
for(i in 1:length(ys)) hypoxia[[i]] <- as(hypoxia[[i]], 'Spatial')

# Kriging Formula
# varioForm <- as.formula("DOXY ~ DEPH")
varioForm <- as.formula("sat ~ DEPH")

# Interpolation
kriged <- vector('list', length(ys))
for(i in 1:length(ys)) {
  # Variogram
  varioFit <- automap::autofitVariogram(varioForm, hypoxia[[i]])
  plot(varioFit)

  # Krigging model
  kriged[[i]] <- gstat::krige(varioForm, hypoxia[[i]], coords, model = varioFit$var_model) %>%
                 st_as_sf(kriged)
}

# Single object
hypoxia <- bathy[, 'geometry']
for(i in 1:length(ys)) hypoxia <- cbind(hypoxia, kriged[[i]][, 1, drop = T])
colnames(hypoxia)[1:length(ys)] <- ys

# Annual mean of O2 saturation in deep waters
hypoxia$sat <- rowMeans(hypoxia[, as.character(ys), drop = T])


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                             HYPOXIC STRESS INDEX
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Oxygen saturation sigmoidal function (%)
sig <- function(x, L = 1, Q = 2e2, k = .15) {
  dat <- -L / (.99 + (Q * exp(-k * x))) + L
  dat[dat < 0] <- 0
  dat
}
x <- seq(0, 100, by = .1)
# png('./figures/hypFunction.png', width = 1000, height = 800, res = 300, pointsize = 6)
plot(x, sig(x), pch = 20, cex = .5, type = 'l', lwd = 1, xlab = 'Oxygen saturation (%)', ylab = 'Hypoxia index', frame = F)
# dev.off()

# Hypoxia
hypoxia$Hypoxia <- sig(hypoxia$sat)
hypoxia <- hypoxia[, 'Hypoxia']

# Remove cells with values = 0
hypoxia <- hypoxia[hypoxia$Hypoxia > 0, ]


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                  EXPORT DATA
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Export object as .RData
save(hypoxia, file = './Data/Driver/Hypoxia.RData')


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                 VISUALIZE DATA
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
png('./Figures/Hypoxia.png', width = 1280, height = 1000, res = 200, pointsize = 6)
plot(hypoxia[, 'Hypoxia'], border = 'transparent')
dev.off()
