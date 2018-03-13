library(rgbif)
foo <- occ_data(geometry=c(-125.0, 38.4, -121.8, 40.9), limit = 20)
foo <- occ_data(geometry='POLYGON((30.1 10.1, 10 20, 20 40, 40 40, 30.1 10.1))', limit = 200)

m <- matrix(c(30.1, 10.1, 10, 20, 20, 40, 40, 40, 30.1, 10.1), ncol = 2, byrow = TRUE)
plot(m, col = 'red')
polygon(m[, 1], m[, 2])
points(foo$data[, c('decimalLongitude', 'decimalLatitude')])
