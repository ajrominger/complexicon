library(igraph)
library(pika)

## set up the network
set.seed(2)
foo <- rstick(1000, r = 0.5)

if(sum(foo) %% 2) foo[1] <- foo[1] + 1
x <- sample_degseq(foo, method = 'vl')

## add extra edges for the redwood
xadj <- as_adj(x)
redwood <- which.max(rowSums(as.matrix(xadj)))
rwAdd <- sample(0:1, nrow(xadj), rep = TRUE, prob = c(0.1, 0.9))
xadj[redwood, ] <- xadj[, redwood] <- xadj[redwood, ] + rwAdd
xadj[xadj > 1] <- 1
diag(xadj) <- 0
xnew <- graph_from_adjacency_matrix(xadj, mode = 'undirected')

## node colors
vcol <- rep('blue', nrow(xadj))
vcol[redwood] <- 'red'

## make a fixed layout for the network
set.seed(1)
l <- layout_with_drl(xnew)
l[, 1] <- l[, 1] - l[redwood, 1]
l[, 2] <- l[, 2] - l[redwood, 2]
l[redwood, ]


## annimation code kindly borrowed from 
## https://davetang.org/muse/2015/02/12/animated-plots-using-r

setwd('./keystone/temp')

## number of frames
frames <- 50

## the different zooms for each plot
zooms <- seq(0.1, 1, length.out = frames)

## function for creating file name with leading zeros
## makes it easier to process them sequentially
rename <- function(x){
    if (x < 10) {
        return(name <- paste('000',i,'plot.png',sep=''))
    }
    if (x < 100 && i >= 10) {
        return(name <- paste('00',i,'plot.png', sep=''))
    }
    if (x >= 100) {
        return(name <- paste('0', i,'plot.png', sep=''))
    }
}

## loop through plots
for(i in 1:frames){
    name <- rename(i)
    
    ## saves the plot as a .png file in the working directory
    png(name)
    
    par(mar = rep(0.1, 4))
    plot(xnew, vertex.size = 10, vertex.label = NA, vertex.color = vcol, vertex.frame.color = NA,
         edge.width = 0.2, edge.color = 'black',
         layout = l, rescale = FALSE,
         xlim = zooms[i] * range(l[, 1]), ylim = zooms[i] * range(l[, 2]))
    box()
    
    dev.off()
}

## run ImageMagick
system('convert *.png -delay 3 -loop 0 redwood_zoom.gif')

