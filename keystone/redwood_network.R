library(igraph)
library(pika)

## variables to control appearances
rwCol <- 'red'
nodeCol <- 'blue'
deadCol <- 'yellow'
extCol <- 'white'
bgCol <- rgb(143, 157, 62, alpha = 200, maxColorValue = 255)

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

## make a fixed layout for the network
set.seed(1)
l <- layout_with_drl(xnew)
l[, 1] <- l[, 1] - l[redwood, 1]
l[, 2] <- l[, 2] - l[redwood, 2]


## annimation code kindly borrowed from 
## https://davetang.org/muse/2015/02/12/animated-plots-using-r

## set-up temporary place to save individual frames
system('mkdir ~/Dropbox/Research/complexicon/keystone/temp')
setwd('~/Dropbox/Research/complexicon/keystone/temp')


## number of frames
frames <- 50

## the different zooms for each plot
zooms <- seq(0.05, 1, length.out = frames)


## loop through plots
for(i in 1:frames){
    ## save the plot as a .png file in the working directory
    png(paste0('frame_', paste0(rep(0, 4 - nchar(as.character(i))), collapse = ''), i, '.png', collapse = ''))
    
    par(mar = rep(1.5, 4), xpd = NA)
    plot(l * 1.1, type = 'n', xlim = zooms[i] * range(l[, 1]), ylim = zooms[i] * range(l[, 2]), asp = 1, axes = FALSE)
    polygon(1.1 * max(abs(l)) * cos(seq(0, 1.99 * pi, length.out = 100)), 
            1.1 * max(abs(l)) * sin(seq(0, 1.99 * pi, length.out = 100)), 
            col = bgCol, border = NA)
    plot(xnew, vertex.size = 10, vertex.label = NA, vertex.color = vcol, vertex.frame.color = NA,
         edge.width = 0.2, edge.color = 'black', edge.arrow.size = 0,
         layout = l, rescale = FALSE, add = TRUE)
    
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

## loop through extinctions
for(i in 1:50) {
    if(i == 1) {
        theseDead <- redwood
    } else {
        theseDead <- which(colSums(as.matrix(xadjExt)) == 0)
    }
    
    theseDeadInd <- which(!(1:nrow(xadjKey) %in% rownames(xadjExt)))
    vcol[theseDeadInd] <- deadCol
    xext <- delete_edges(xext, E(xext)[to(theseDeadInd)|from(theseDeadInd)])
    
    png(paste0('frame_', paste0(rep(0, 4 - nchar(as.character(i))), collapse = ''), i, '.png', collapse = ''))
    
    par(mar = rep(1.5, 4), xpd = NA)
    plot(l * 1.1, type = 'n', 
         xlim = zooms[length(zooms)] * range(l[, 1]), 
         ylim = zooms[length(zooms)] * range(l[, 2]), 
         asp = 1, axes = FALSE)
    polygon(1.1 * max(abs(l)) * cos(seq(0, 1.99 * pi, length.out = 100)), 
            1.1 * max(abs(l)) * sin(seq(0, 1.99 * pi, length.out = 100)), 
            col = bgCol, border = NA)
    plot(xext, vertex.size = 10, vertex.label = NA, vertex.color = vcol, vertex.frame.color = NA,
         edge.width = 0.2, edge.color = 'black', edge.arrow.size = 0,
         layout = l, rescale = FALSE, add = TRUE)
    
    legend(par('usr')[2], par('usr')[4], legend = 'extinctions', col = deadCol, pch = 16, bg = bgCol, 
           box.col = bgCol, xjust = 0.7, yjust = 0.5)
    
    dev.off()
    
    if(length(theseDead) == 0) break
    xadjExt <- xadjExt[-theseDead, -theseDead]
}

## run ImageMagick
system('convert *.png -delay 30 -loop 0 ../redwood_extinction.gif')

## clean up
system('rm -f -r ~/Dropbox/Research/complexicon/keystone/temp')
