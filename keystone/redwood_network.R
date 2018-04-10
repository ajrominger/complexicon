library(igraph)
library(pika)
library(png)

## variables to control appearances
rwCol <- hsv(0.05, 1, 0.7)
nodeCol <- hsv(0.6, 0.7, 0.8)
deadCol <- 'black'
edgeCol <- hsv(0.9, 0.5, 0.9)
bgCol <- rgb(143, 157, 62, alpha = 200, maxColorValue = 255)

## images to add to network
rw <- readPNG('~/Dropbox/Research/complexicon/keystone/web_images/redwood.png')
microbe <- readPNG('~/Dropbox/Research/complexicon/keystone/web_images/shutterstock_537846901.png')
bryophyte <- readPNG('~/Dropbox/Research/complexicon/keystone/web_images/bryophyte.png')
insect <- readPNG('~/Dropbox/Research/complexicon/keystone/web_images/2831012.png')
bird <- readPNG('~/Dropbox/Research/complexicon/keystone/web_images/shutterstock_141320173.png')


## set up the network
set.seed(120)
x <- sample_pa(1000, 0.1, directed = TRUE)

## add extra edges for the redwood
xadj <- as_adj(x)
redwood <- which.max(rowSums(as.matrix(xadj)))
rwAdd <- sample(0:1, nrow(xadj), rep = TRUE, prob = c(0.1, 0.9))
xadj[, redwood] <- xadj[, redwood] + rwAdd
xadj[xadj > 1] <- 1
diag(xadj) <- 0
xnew <- graph_from_adjacency_matrix(xadj, mode = 'directed')

## node colors
vcol <- rep(nodeCol, nrow(xadj))
vcol[redwood] <- rwCol

## node sizes
vsize <- rep(15, nrow(xadj))
vsize[redwood] <- 3.5 * max(vsize)

## make a fixed layout for the network
set.seed(1)
l <- layout_with_drl(xnew)
l[, 1] <- l[, 1] - l[redwood, 1]
l[, 2] <- l[, 2] - l[redwood, 2]

## add some space around 0
d <- sqrt(l[, 1]^2 + l[, 2]^2)
a <- acos(l[, 1] / d) * sign(l[, 2])
l <- cbind((d + 0.05*max(d)) * cos(a), (d + 0.05*max(d)) * sin(a))
l[redwood, ] <- c(0, 0)

## annimation code kindly borrowed from 
## https://davetang.org/muse/2015/02/12/animated-plots-using-r

## set-up temporary place to save individual frames
system('mkdir ~/Dropbox/Research/complexicon/keystone/temp')
setwd('~/Dropbox/Research/complexicon/keystone/temp')


## number of frames
frames <- 50

## the different zooms for each plot (a little extra at the beginning)
zooms <- c(rep(0.038, 4), seq(0.038, 1, length.out = frames - 4))


## loop through plots
for(i in 1:frames){
    ## save the plot as a .png file in the working directory
    png(paste0('frame_', paste0(rep(0, 4 - nchar(as.character(i))), collapse = ''), i, 
               '.png', collapse = ''))
    
    par(mar = rep(1.5, 4), xpd = NA)
    plot(l * 1.1, type = 'n', xlim = zooms[i] * range(l[, 1]), 
         ylim = zooms[i] * range(l[, 2]), 
         asp = 1, axes = FALSE)
    # polygon(1.1 * max(abs(l)) * cos(seq(0, 1.99 * pi, length.out = 100)), 
    #         1.1 * max(abs(l)) * sin(seq(0, 1.99 * pi, length.out = 100)), 
    #         col = bgCol, border = NA)
    plot(xnew, vertex.size = vsize, vertex.label = NA, vertex.color = vcol, 
         vertex.frame.color = NA,
         edge.width = 0.5, edge.color = edgeCol, edge.arrow.size = 0,
         layout = l, rescale = FALSE, add = TRUE)
    
    ## add images
    rasterImage(rw, -0.25, -0.25, 0.25, 0.25)
    rasterImage(microbe, -0.51, 0.005, -0.37, 0.145)
    rasterImage(bird, 0.2, -0.44, 0.335, -0.305)
    rasterImage(bryophyte, -0.605, -0.33, -0.465, -0.19)
    rasterImage(insect, 0.48, 0.105, 0.62, 0.245)
    
    dev.off()
}

## run ImageMagick
system('convert *.png -delay 3 -loop 0 ../redwood_zoom.gif')

## clean up
system('rm -f -r ~/Dropbox/Research/complexicon/keystone/temp')



## simulate extinction

## set-up temporary place to save individual frames
system('mkdir ~/Dropbox/Research/complexicon/keystone/temp')
setwd('~/Dropbox/Research/complexicon/keystone/temp')

## make adj mat to keep track of key interactions (without which, sp dies)
xadj <- as.matrix(as_adj(xnew, type = 'both'))

xadjKey <- sapply(1:nrow(xadj), function(i) {
    x <- xadj[i, ]
    n <- rtpois(1, 1)
    if(n > sum(x)) n <- sum(x)
    out <- numeric(length(x))
    out[sample(which(x == 1), n)] <- 1
    
    return(out)
})

## make adj mat to track extinctions
rownames(xadjKey) <- colnames(xadjKey) <- 1:nrow(xadjKey)
xadjExt <- xadjKey

## re-set values for trouble-shooting
vcol <- rep(nodeCol, nrow(xadj))
vcol[redwood] <- rwCol
xext <- xnew
xext <- set_edge_attr(xext, 'color', value = edgeCol)

## loop through extinctions
for(i in 1:5) {
    if(i == 1) {
        theseDead <- redwood
    } else {
        theseDead <- which(colSums(as.matrix(xadjExt)) == 0)
    }
    
    theseDeadInd <- which(!(1:nrow(xadjKey) %in% rownames(xadjExt)))
    vcol[theseDeadInd] <- deadCol
    # xext <- delete_edges(xext, E(xext)[to(theseDeadInd)|from(theseDeadInd)])
    xext <- set_edge_attr(xext, 'color', 
                          index = E(xext)[to(theseDeadInd)|from(theseDeadInd)], 
                          value = deadCol)
    
    ## multiple frames to slow it down
    for(j in 1:5) {
        png(paste0('frame_', paste0(rep(0, 4 - nchar(as.character(i))), collapse = ''), 
                   i, '-', j, '.png', collapse = ''))
        
        par(mar = rep(1.5, 4), xpd = NA)
        plot(l * 1.1, type = 'n', 
             xlim = zooms[length(zooms)] * range(l[, 1]), 
             ylim = zooms[length(zooms)] * range(l[, 2]), 
             asp = 1, axes = FALSE)
        # polygon(1.1 * max(abs(l)) * cos(seq(0, 1.99 * pi, length.out = 100)), 
        #         1.1 * max(abs(l)) * sin(seq(0, 1.99 * pi, length.out = 100)), 
        #         col = bgCol, border = NA)
        plot(xext, vertex.size = vsize, vertex.label = NA, vertex.color = vcol, 
             vertex.frame.color = NA,
             edge.width = 0.5, 
             # edge.color = edgeCol, 
             edge.arrow.size = 0,
             layout = l, rescale = FALSE, add = TRUE)
        
        legend(par('usr')[2], par('usr')[4], legend = 'extinctions', col = deadCol, 
               pch = 16, cex = 1.3,
               bg = 'transparent', box.col = 'transparent', 
               xjust = 0.8, yjust = 0.5)
        
        dev.off()
    }
    
    if(length(theseDead) == 0) break
    xadjExt <- xadjExt[-theseDead, -theseDead]
}

## run ImageMagick
system('convert *.png -delay 3 -loop 0 ../redwood_extinction.gif')

## clean up
system('rm -f -r ~/Dropbox/Research/complexicon/keystone/temp')
