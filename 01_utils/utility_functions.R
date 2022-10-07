options(stringsAsFactors = F)

library(umap)
library(scales)
library(MASS)
library(RColorBrewer)
library(colorRamps)
library(vegan)
library(igraph)
library(entropy)
library(jpeg)
library(repr)
library(shape)
library(gtools)
library(entropy)

###########################
#####Utility functions#####
###########################
#Function to split a vector into n chunks with x amount of overlap
splitWithOverlap <- function(vec, seg.length, overlap) {
  starts = seq(1, nrow(vec), by = seg.length - overlap)
  ends   = starts + seg.length - 1
  ends[ends > nrow(vec)] = nrow(vec)
  
  lapply(1:length(starts), function(i)
    vec[starts[i]:ends[i],])
}

#Function to darken a given color, useful for plotting
darken_color = function(color, factor = 1.4) {
  col <- col2rgb(color)
  col <- col / factor
  col <- rgb(t(col), maxColorValue = 255)
  col
}

#Function to fit ellipse to points
fit.ellipse <- function (x, y = NULL) {

  EPS <- 1.0e-8
  dat <- xy.coords(x, y)
  
  D1 <- cbind(dat$x * dat$x, dat$x * dat$y, dat$y * dat$y)
  D2 <- cbind(dat$x, dat$y, 1)
  S1 <- t(D1) %*% D1
  S2 <- t(D1) %*% D2
  S3 <- t(D2) %*% D2
  T <- -solve(S3) %*% t(S2)
  M <- S1 + S2 %*% T
  M <- rbind(M[3,] / 2, -M[2,], M[1,] / 2)
  evec <- eigen(M)$vec
  cond <- 4 * evec[1,] * evec[3,] - evec[2,] ^ 2
  a1 <- evec[, which(cond > 0)]
  f <- c(a1, T %*% a1)
  names(f) <- letters[1:6]
  
  # calculate the center and lengths of the semi-axes
  A <-
    matrix(
      c(2 * f[1], f[2], f[2], 2 * f[3]),
      nrow = 2,
      ncol = 2,
      byrow = T
    )
  b <- matrix(c(-f[4], -f[5]),
              nrow = 2,
              ncol = 1,
              byrow = T)
  soln <- solve(A) %*% b
  b2 <- f[2] ^ 2 / 4
  
  center <- c(soln[1], soln[2])
  names(center) <- c("x", "y")
  
  num  <-
    2 * (f[1] * f[5] ^ 2 / 4 + f[3] * f[4] ^ 2 / 4 + f[6] * b2 - f[2] * f[4] *
           f[5] / 4 - f[1] * f[3] * f[6])
  den1 <- (b2 - f[1] * f[3])
  den2 <- sqrt((f[1] - f[3]) ^ 2 + 4 * b2)
  den3 <- f[1] + f[3]
  
  semi.axes <-
    sqrt(c(num / (den1 * (den2 - den3)),  num / (den1 * (-den2 - den3))))
  
  # calculate the angle of rotation
  term <- (f[1] - f[3]) / f[2]
  angle <- atan(1 / term) / 2
  
  list(
    coef = f,
    center = center,
    major = max(semi.axes),
    minor = min(semi.axes),
    angle = unname(angle)
  )
}

#Function calculate Euclidean distance
euc.dist <- function(x1, x2)
  sqrt(sum((x1 - x2) ^ 2))
