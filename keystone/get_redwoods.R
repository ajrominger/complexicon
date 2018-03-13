library(rgbif)
library(rgdal)

setwd('~/Dropbox/Research/complexicon/keystone')
source('~/Dropbox/api_auth.R')

rwPoly <- readOGR('redwoods.kml', 'redwoods.kml')
rwCoord <- rwPoly@polygons[[1]]@Polygons[[1]]@coords
rwGBIFPoly <-
    sprintf('POLYGON((%s))', paste(paste(rwCoord[, 1], rwCoord[, 2], sep = ' '), collapse = ', '))

rwDownload <- occ_download(sprintf('geometry = %s', rwGBIFPoly))

ready <- FALSE
allTime <- 0

while (!ready) {
    Sys.sleep(30)
    m <- occ_download_meta(rwDownload)
    ready <- m$status == 'SUCCEEDED'
    allTime <- allTime + 30
    print(paste('time waited:', allTime / 30))
    
    if (allTime > 45 * 60)
        break
}

if (ready) {
    rwOcc <- occ_download_import(occ_download_get(m$key))
    write.csv(rwOcc, file = 'rwOcc.csv', row.names = FALSE)
} else {
    print('failed to download')
}
