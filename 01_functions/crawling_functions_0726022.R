#rm(list=ls());
options(stringsAsFactors=F);
library(umap)
library(scales)
library(MASS)
library(RColorBrewer)
library(colorRamps)
library(vegan)
library(igraph)
library(entropy)
library(jpeg)

####################
#####Functions######
####################
#Function load cell image .jpegs
load_images = function(image_directory){
  setwd(image_directory)
  
  #Get directories
  files = list.files()
  
  dat = list()
  names = c()
  for(i in 1:length(files)){
    
    setwd(paste(files[i], '/', 'masks_normalized/', sep = ''))
    f = list.files(full.names = FALSE)
    f = gtools::mixedsort(sort(f))
    
    imgs = list()
    for(j in 1:length(f)){
      imgs[[f[j]]] = unlist(as.data.frame(readJPEG(f[j])))
    }
    
    dat[[as.character(files[i])]] = imgs
    names = c(names, rep(files[i], length(imgs)))
    
    setwd('../../')
  }
  
  l = list(dat, names)
  names(l) = c('images', 'names')
  return(l)
}

#Function to split a vector into n chunks with x amount of overlap
splitWithOverlap <- function(vec, seg.length, overlap) {
  starts = seq(1, nrow(vec), by=seg.length-overlap)
  ends   = starts + seg.length - 1
  ends[ends > nrow(vec)] = nrow(vec)
  
  lapply(1:length(starts), function(i) vec[starts[i]:ends[i],])
}

#Function to darken a given color, useful for plotting
darken_color = function(color, factor=1.4){
  col <- col2rgb(color)
  col <- col/factor
  col <- rgb(t(col), maxColorValue=255)
  col
}

#Function to calculate per bin entropies
calculate_bin_entropies = function(bin,
                                   strain_coords,
                                   strain_xy,
                                   window = 30,
                                   n_bins,
                                   plot_trajectories = FALSE,
                                   calculate_entropies = FALSE,
                                   n_trajectories = NULL){
  
  #Choose bin
  z = bin
  tmp = strain_coords
  
  #Split on bin and extract bouts in between
  y = which(!tmp == z)
  idx <- c(0, cumsum(abs(diff(y)) > 1))
  indices = split(y, idx)
  
  print("Extracting inter-bin bouts")
  #Get bin name
  series = lapply(indices, function(x) tmp[x])
  
  #Remove first bout (doesn't originate at bin of interest)
  series = series[-1]
  
  #Require bout of length n
  series = series[lapply(series, length)>=window]
  series = lapply(series, function(x) x[1:window])
  
  #Make empty vector for saving output
  entropies = c()
  
  if(length(series)>1){
    
    #Combine into dataframe
    series = do.call(cbind, series)
    print(ncol(series))
    
    #Add row corresponding to starting bin
    series = rbind(rep(z, ncol(series)), series)
    
    #Plot if desired
    if(plot_trajectories == TRUE){
      plot(series[,1],
           type = "l",
           lwd = 2,
           col = alpha("grey40", 0.5),
           ylim = c(min(unlist(series)), max(unlist(series))),
           ylab = "Bin",
           xlab = "Time",
           cex.axis = 1.5,
           cex.lab = 1.5,
           bty = 'n')
      for(i in 2:ncol(series)){
        lines(series[,i],
              lwd = 2,
              col = alpha("grey40", 0.5))}}
    
    if(!is.null(n_trajectories)){
      
      if(ncol(series)>n_trajectories){
        print("Calculating entropies")
        
        series = series[,sample(seq(1, ncol(series), 1), n_trajectories)]
        print(ncol(series))
        
        ##Calculate probability density functions
        for(h in 1:nrow(series)){
          
          #print(i)
          #Extract row(s) of interest
          row = unlist(series[h,])
          
          #Convert to coords
          row = cbind(as.numeric(unlist(lapply(strsplit(apply(as.data.frame(row), 1, function(x) unique(strain_xy[strain_coords == x])), "_"), function(v){v[1]}))),
                      as.numeric(unlist(lapply(strsplit(apply(as.data.frame(row), 1, function(x) unique(strain_xy[strain_coords == x])), "_"), function(v){v[2]}))))
          
          w = MASS::kde2d(row[,1],
                          row[,2],
                          h = c(1, 1),
                          n = n_bins,
                          lims = c(1, n_bins+1, 1, n_bins+1))$z
          
          pdf = as.numeric(unlist(as.data.frame(w)))
          entropies = c(entropies, entropy::entropy(pdf, unit = 'log2'))
        }
        #plot(entropies,
        #type = "l",
        #lwd = 1.5,
        #cex.lab = 1.5,
        #cex.axis = 1.5,
        #xlab = "Time",
        #bty = 'n',
        #ylim = c(0, 4))
      }
      #else{
      #entropies = c(entropies, NA)
      #}
    }
    
    if(calculate_entropies == TRUE){
      print("Calculating entropies")
      ##Calculate probability density functions
      for(h in 1:nrow(series)){
        #print(i)
        #Extract row(s) of interest
        row = unlist(series[h,])
        
        #Convert to coords
        row = cbind(as.numeric(unlist(lapply(strsplit(apply(as.data.frame(row), 1, function(x) unique(strain_xy[strain_coords == x])), "_"), function(v){v[1]}))),
                    as.numeric(unlist(lapply(strsplit(apply(as.data.frame(row), 1, function(x) unique(strain_xy[strain_coords == x])), "_"), function(v){v[2]}))))
        
        w = MASS::kde2d(row[,1],
                        row[,2],
                        h = c(1, 1),
                        n = n_bins,
                        lims = c(1, n_bins+1, 1, n_bins+1))$z
        
        pdf = as.numeric(unlist(as.data.frame(w)))
        entropies = c(entropies, entropy::entropy(pdf, unit = 'log2'))
      }
      #plot(entropies,
      #type = "l",
      #lwd = 1.5,
      #cex.lab = 1.5,
      #cex.axis = 1.5,
      #xlab = "Time",
      #bty = 'n',
      #ylim = c(0, 4))
    }
    
    l = list(series, entropies)
    names(l) = c("bouts", "entropies")
    return(l)
  }
}

#Function to fit ellipse to points
fit.ellipse <- function (x, y = NULL) {
  # from:
  # http://r.789695.n4.nabble.com/Fitting-a-half-ellipse-curve-tp2719037p2720560.html
  #
  # Least squares fitting of an ellipse to point data
  # using the algorithm described in:
  #   Radim Halir & Jan Flusser. 1998.
  #   Numerically stable direct least squares fitting of ellipses.
  #   Proceedings of the 6th International Conference in Central Europe
  #   on Computer Graphics and Visualization. WSCG '98, p. 125-132
  #
  # Adapted from the original Matlab code by Michael Bedward (2010)
  # michael.bedward@gmail.com
  #
  # Subsequently improved by John Minter (2012)
  #
  # Arguments:
  # x, y - x and y coordinates of the data points.
  #        If a single arg is provided it is assumed to be a
  #        two column matrix.
  #
  # Returns a list with the following elements:
  #
  # coef - coefficients of the ellipse as described by the general
  #        quadratic:  ax^2 + bxy + cy^2 + dx + ey + f = 0
  #
  # center - center x and y
  #
  # major - major semi-axis length
  #
  # minor - minor semi-axis length
  #
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
  cond <- 4 * evec[1,] * evec[3,] - evec[2,]^2
  a1 <- evec[, which(cond > 0)]
  f <- c(a1, T %*% a1)
  names(f) <- letters[1:6]
  
  # calculate the center and lengths of the semi-axes
  #
  # see http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2288654/
  # J. R. Minter
  # for the center, linear algebra to the rescue
  # center is the solution to the pair of equations
  # 2ax +  by + d = 0
  # bx  + 2cy + e = 0
  # or
  # | 2a   b |   |x|   |-d|
  # |  b  2c | * |y| = |-e|
  # or
  # A x = b
  # or
  # x = Ainv b
  # or
  # x = solve(A) %*% b
  A <- matrix(c(2*f[1], f[2], f[2], 2*f[3]), nrow=2, ncol=2, byrow=T )
  b <- matrix(c(-f[4], -f[5]), nrow=2, ncol=1, byrow=T)
  soln <- solve(A) %*% b
  b2 <- f[2]^2 / 4
  
  center <- c(soln[1], soln[2])
  names(center) <- c("x", "y")
  
  num  <- 2 * (f[1] * f[5]^2 / 4 + f[3] * f[4]^2 / 4 + f[6] * b2 - f[2]*f[4]*f[5]/4 - f[1]*f[3]*f[6])
  den1 <- (b2 - f[1]*f[3])
  den2 <- sqrt((f[1] - f[3])^2 + 4*b2)
  den3 <- f[1] + f[3]
  
  semi.axes <- sqrt(c( num / (den1 * (den2 - den3)),  num / (den1 * (-den2 - den3)) ))
  
  # calculate the angle of rotation
  term <- (f[1] - f[3]) / f[2]
  angle <- atan(1 / term) / 2
  
  list(coef=f, center = center, major = max(semi.axes), minor = min(semi.axes), angle = unname(angle))
}

#Function to plot markov on umap
plot_umap_markov = function(layout,
                            chord,
                            centers,
                            plot_umap_points = TRUE,
                            plot_self = FALSE,
                            cols = RColorBrewer::brewer.pal(nrow(centers), 'Paired'),
                            umap_cols = 'grey90',
                            plot_dynamic_points = FALSE,
                            arr.length = 0,
                            arr.width = 0){
  
  #Plotting in umap
  #Reformat chord diagram coords
  chord$fromx = centers[match(chord$rn, rownames(centers)), 1]
  chord$fromy = centers[match(chord$rn, rownames(centers)), 2]
  chord$tox = centers[match(chord$cn, rownames(centers)), 1]
  chord$toy = centers[match(chord$cn, rownames(centers)), 2]
  
  #Split into self and not
  chord_s = chord[chord$rn == chord$cn,]
  chord = chord[!chord$rn == chord$cn,]
  
  #Make colors not alpha
  chord_s$col = substr(chord_s$col, 1, 7)
  chord$col = substr(chord$col, 1, 7)
  
  #Remove zeroes
  chord = chord[!chord$value1==0,]
  
  if(plot_umap_points == FALSE){
    plot(layout$x,
         layout$y,
         pch = 20,
         #col = "grey90",
         col = NULL,
         bty = 'n',
         xaxt = 'n',
         yaxt = 'n',
         xlab = "",
         ylab = "")
  }else{
    plot(layout$x,
         layout$y,
         pch = 20,
         col = umap_cols,
         bty = 'n',
         cex = 0.5,
         xaxt = 'n',
         yaxt = 'n',
         xlab = "",
         ylab = "")
    
  }
  
  for(i in 1:nrow(chord)){
    diagram::curvedarrow(from = c(chord$fromx[i], chord$fromy[i]),
                         to = c(chord$tox[i], chord$toy[i]),
                         lcol = chord$col[i],
                         lwd = chord$value1[i]*50,
                         arr.pos = 0.9,
                         curve = 0.2,
                         arr.type = "triangle",
                         arr.length = arr.length,
                         arr.width = arr.width
    )}
  
  if(plot_self == TRUE){
    for(i in 1:nrow(chord)){
      diagram::selfarrow(c(chord$fromx[i], chord$fromy[i]),
                         lcol = chord$col[i],
                         lwd = chord$value1[i]*50,
                         curve = 1,
                         arr.length = 0,
                         arr.width = 0)
    }
  }
  
  if(plot_dynamic_points == TRUE){
    f = table(k$cluster)
    f = f/sum(f)
    
    for(i in 1:nrow(centers)){
      points(k$centers[i,1],
             k$centers[i,2],
             pch = 21,
             bg = cols[i],
             col = darken_color(cols[i]),
             cex = 4+(2*f[i]))
    }
  }else{
    points(centers,
           pch = 21,
           bg = darken_color(cols),
           col = darken_color(cols),
           cex = 3)
  }
}

#Function to plot parameters by behavior space position
plot_parameter = function(parameter, layout, n_bins, return = FALSE, ...){
  
  layout = bin_umap(layout, n_bins = n_bins)$layout
  
  a = parameter
  a = split(a[1:nrow(layout)], layout$xy_new)
  a = unlist(lapply(a, function(x) mean(x, na.rm = TRUE)))
  
  p = expand.grid(seq(1, n_bins+1, 1), seq(1, n_bins+1, 1))
  p = paste(p[,1], p[,2], sep = '_')
  a = a[match(p, names(a))]
  names(a) = p
  
  image(matrix(a, nrow = n_bins+1),
        xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', bty = 'n', ...)
  
  if(return == TRUE){
    return(a)
  }
}

#Function to extract windows of user defined features, size x and stepsize y
get_windows = function(features,
                       window_size = 10,
                       step_size = 1,
                       name = NULL){
  
  #Get windows
  peaks1 = splitWithOverlap(features, window_size, window_size-step_size)
  
  #Clean and combine into matrix
  peaks1 = peaks1[1:(length(peaks1)-window_size)]
  peaks1 = lapply(peaks1, function(x) unlist(as.data.frame(x)))
  peaks1 = do.call(cbind, peaks1)
  
  if(is.null(name) == FALSE){
    colnames(peaks1) = paste(name, '_', rownames(features)[1:ncol(peaks1)], sep = '')
  }else{
    colnames(peaks1) = rownames(features)[1:ncol(peaks1)]
  }
  
  return(peaks1)
}

#Function to extract velocity windows and (optionally) regularize, size x and stepsize y
#Velocities are expected to be in columns named: 'translational_velocity', 'angular_velocity', 'sideslip'
get_velocity_windows = function(features,
                                include_sideslip = FALSE,
                                return_xy_windows = FALSE,
                                window_size = 1,
                                step_size = 1,
                                symm = FALSE,
                                verbose = FALSE,
                                name = NULL){
  
  #Initialize window list
  frags = list()
  xys = list()
  
  #Get trajectories
  for(i in seq(1, (nrow(features)-window_size), step_size)){
    
    if(verbose == TRUE){
      if(i %% 10000 == TRUE){
        print(i)
      }
    }
    
    #Get velocity vectors
    vt = features$translational_velocity[seq(i, i+window_size, 1)]
    vr = features$angular_velocity[seq(i, i+window_size, 1)]
    if(include_sideslip == TRUE){
      vs = features$sideslip[seq(i, i+window_size, 1)]
    }
    x = features$x[seq(i, i+window_size, 1)]
    y = features$y[seq(i, i+window_size, 1)]
    
    time = features$time[seq(i, i+window_size, 1)]
    
    #Subtract first frame make t0 = 0
    vr = vr-vr[1]
    if(include_sideslip == TRUE){
      vs = vs-vs[1]
    }
    
    #Multiply by sign of second frame to normalize turn direction
    if(symm == TRUE){
      if(vr[2]<0){
        vr = vr*(-1)
      }
      
      if(include_sideslip == TRUE){
        if(vs[2]<0){
          vs = vs*(-1)
        }      }
      
      vr = abs(vr)
      
      if(include_sideslip == TRUE){
        vs = abs(vs)
      }
    }
    
    #Look at time difference
    t_diff = diff(time)
    
    if(include_sideslip == TRUE){
      frags[[paste(time[1], "_", time[length(time)], "_", name, sep = "")]] = c(vt, vr, vs)
      xys[[paste(time[1], "_", time[length(time)], "_", name, sep = "")]] = c(x, y) 
    }else{
      frags[[paste(time[1], "_", time[length(time)], "_", name, sep = "")]] = c(vt, vr) 
      xys[[paste(time[1], "_", time[length(time)], "_", name, sep = "")]] = c(x, y)
    }
  }
  
  #Remove elements with NAs
  frags = frags[lapply(frags, function(x) sum(is.na(x)))<1]
  xys = xys[lapply(xys, function(x) sum(is.na(x)))<1]
  
  #Combine into df
  if(return_xy_windows == TRUE){
    df = do.call(cbind, frags)
    xys_df = do.call(cbind, xys)
    
    l = list(df, xys_df)
    names(l) = c("velocity_windows", "xy_windows")
    return(l)
  }else{
    df = do.call(cbind, frags)
    return(df)
  }
}

#Function to bin a umap space into a n x n grid (umap coordinates are provided as the 'layout' object)
bin_umap = function(layout,
                    n_bins){
  
  n = n_bins
  
  #Split x into n bins
  x1 = seq(min(layout[,1]),
           max(layout[,1]),
           (max(layout[,1])-min(layout[,1]))/n)
  names(x1) = seq(1, n+1, 1)
  
  xnew = apply(layout, 1, function(x) names(x1)[which.min(abs(as.numeric(x[1]) - x1))])
  layout$xnew = xnew
  
  #Split y into n bins
  y1 = seq(min(layout[,2]),
           max(layout[,2]),
           (max(layout[,2])-min(layout[,2]))/n)
  names(y1) = seq(1, n+1, 1)
  
  ynew = apply(layout, 1, function(x) names(y1)[which.min(abs(as.numeric(x[2]) - y1))])
  layout$ynew = ynew
  
  #Paste xy to get unique bin combos (this will be input to sling shot as 'clusters')
  xy_new = paste(xnew, ynew, sep = "_")
  layout$xy_new = xy_new
  
  #Convert coordinates to numeric, sort first
  m = unique(xy_new)
  
  y = as.numeric(unlist(lapply(strsplit(m, "_"), function(v){v[2]})))
  m = m[order(y)]
  
  x = as.numeric(unlist(lapply(strsplit(m, "_"), function(v){v[1]})))
  m = m[order(x)]
  
  #Get names
  names(m) = seq(1, length(m), 1)
  names(xy_new) = names(m)[match(xy_new, m)]
  
  #Get vector of coords
  layout$coords = as.numeric(names(xy_new))
  
  #Return
  l = list(layout, xy_new)
  names(l) = c("layout", "new_coords")
  
  return(l)
}

#Function to iteratively run umap on windows of a desired size
iterative_umap = function(features,
                          velocity_windows = FALSE,
                          verbose = FALSE,
                          plot = FALSE,
                          step_size = 1,
                          window_size = 30,
                          n_bins = 32,
                          run_umap = TRUE,
                          return_windows = FALSE,
                          ...){
  
  if(verbose == TRUE){
    print("Getting windows")
    }
  #Get windows
  windows = list()
  
  #Set up plots
  if(plot == TRUE){
    n = length(features)*2
    x = ceiling(sqrt(n))
    y = floor(sqrt(n))
    par(mfrow = c(x,y), mar = c(1,1,1,1))
    rm(n, x, y)
  }
  
  for(i in 1:length(features)){
    if(velocity_windows == TRUE){
      #print(i)
      windows[[i]] = get_velocity_windows(features[[i]],
                                          window_size = window_size, 
                                          step_size = step_size,
                                          ...)
    }else{
      windows[[i]] = get_windows(features[[i]],
                                 window_size = window_size, 
                                 step_size = step_size,
                                 ...)
    }
  }
  
  if(run_umap == FALSE){
    #Return
    l = list(features, windows)
    names(l) = c("features", "windows")
    return(l)
  }else{
    if(verbose == TRUE){
      print("Running UMAP")
    }

    #Run UMAP
    umaps = list()
    for(i in 1:length(features)){
      if(verbose == TRUE){
        print(paste("umap", i, "out of", length(features)))
      }
      
      if(verbose == TRUE){
        umaps[[i]] = umap(t(windows[[i]]), verbose = TRUE)
      }else{
        umaps[[i]] = umap(t(windows[[i]]))
      }
      
      if(plot == TRUE){
        plot(umaps[[i]]$layout[,1:2],
             bty = 'n',
             xaxt = 'n',
             yaxt = 'n',
             ylab = "",
             pch = 20,
             xlab = "",
             col = alpha('grey50', 0.5))
        plot(umaps[[i]]$layout[,1:2],
             type = 'l',
             bty = 'n',
             xaxt = 'n',
             yaxt = 'n',
             ylab = "",
             xlab = "",
             col = alpha('grey50', 0.5))
      }
    }
    
    #Extract layouts
    umaps = lapply(umaps, function(z) data.frame(x = z$layout[,1], y = z$layout[,2]))
    
    #Bin
    umaps = lapply(umaps, function(x) bin_umap(x, n_bins = n_bins)$layout)
    
    #Return
    if(return_windows == TRUE){
      l = list(vel, windows, umaps)
      names(l) = c("features", "windows", "umaps")
    }else{
      l = list(umaps)
      names(l) = 'umaps'
    }
    return(l)
  }
}

#Function calculate Euclidean distance
euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))

#Function to calculate distance between umap layouts using Procrustes and Euclidean distances
run_procrustes = function(umaps,
                          run_protest = FALSE){
  
  #Get all combinations of umaps to compare
  x = combn(seq(1, length(umaps), 1), 2)
  
  #Run
  pr_res = c()
  pr_sig = list()
  dists = c()
  for(i in 1:ncol(x)){
    #print("running procrustes")
    pr = procrustes(umaps[[x[1,i]]][,1:2], 
                    umaps[[x[2,i]]][,1:2])
    pr_res = c(pr_res, summary(pr)$rmse)
    
    dists = c(dists, euc.dist(umaps[[x[1,i]]][,1:2], 
                              umaps[[x[2,i]]][,1:2]))
    
    if(run_protest == TRUE){
      #print("running protest")
      pr_sig[[i]] = protest(umaps[[x[1,i]]][,1:2], 
                            umaps[[x[2,i]]][,1:2])
    }
  }
  
  if(run_protest == TRUE){
    res = list(pr_res, pr_sig)
    names(res) = c("procrustes", "protest")
  }else{
    res = list(pr_res, dists)
    names(res) = c("procrustes", "euclidean_distances")
  }
  return(res)
}

#Function to calculate the amount and timing of recurrence in a behavior space
calculate_recurrence = function(umaps,
                                filter_outliers = FALSE,
                                n_bins = 16,
                                plot = FALSE, 
                                verbose = FALSE,
                                threshold = 0.05){
  results = list()
  for(h in 1:length(umaps)){
    
    if(verbose == TRUE){
      print(paste(h, "out of", length(umaps)))
    }
    
    u = umaps[[h]]
    
    if(filter_outliers == TRUE){
      u$x[u$x>30] = 30
      u$x[u$x<(-30)] = -30
      u$y[u$y>30] = 30
      u$y[u$y<(-30)] = -30
    }
    
    #Get distances
    l = bin_umap(u,
                 n_bins = n_bins)$layout
    res = list()
    pos = unique(l$xy_new)
    
    dists = list()
    for(i in 1:length(pos)){
      #if(i%%100 == TRUE){
      # print(paste(i, "out of", length(pos)))
      #}
      #x = as.numeric(as.numeric(unlist(lapply(strsplit(pos[i], "_"), function(v){v[1]}))),
      #               as.numeric(unlist(lapply(strsplit(pos[i], "_"), function(v){v[2]}))))
      x = c(as.numeric(unlist(lapply(strsplit(pos[i], "_"), function(v){v[1]}))),
            as.numeric(unlist(lapply(strsplit(pos[i], "_"), function(v){v[2]}))))
      z = apply(l, 1, function(y) euc.dist(x, c(as.numeric(y[3]), as.numeric(y[4]))))
      dists[[pos[i]]] = z
    }
    
    #Calculate distance distribution
    thresh = quantile(unlist(dists), probs = threshold)
    #10% 
    #5.09902 
    
    #Extract recurrences usins 10% threshold
    recs = lapply(dists, function(x){
      rs = which(x<thresh)
      ds = diff(rs)
      ds[ds>thresh]
    })
    
    if(plot == TRUE){
      histogram = hist(unlist(recs), 
                       breaks = seq(1, max(unlist(recs), na.rm = TRUE), 1), 
                       xlim = c(0,200))
    }else{
      histogram = hist(unlist(recs), 
                       breaks = seq(1, max(unlist(recs), na.rm = TRUE), 1),
                       plot = FALSE)
    }
    
    #Label points recurrent in 1 second bins
    prop_recurrent = list()
    for(i in 1:200){
      z = lapply(recs, function(x) which(x==i))
      z = sum(lapply(z, function(x) length(x))>0)/length(recs)
      prop_recurrent[[i]] = z
    }
    
    #barplot(unlist(prop_recurrent))
    total_recurrent = sum(lapply(recs, function(x) length(x))>0)/length(recs)
    
    l = list(dists, unlist(recs), histogram, unlist(prop_recurrent), total_recurrent)
    names(l) = c("distances", 
                 "recurrences", 
                 "histogram",
                 "proportion_recurrent_in_bins",
                 "total_proportion_recurrent")
    results[[h]] = l
  }
  
  if(plot == TRUE){
    image(do.call(cbind, lapply(results, function(x) x$histogram$counts[1:200]/max(x$histogram$counts[1:200]))),
          col = colorRamps::matlab.like(32),
          xaxt = 'n',
          yaxt = 'n',
          xlab = 'Time',
          ylab = 'Replicate',
          cex.lab = 1.5)
    axis(1,
         at = seq(0, 1, 1/(200-1)),
         labels = seq(1, 200, 1),
         cex.axis = 1.5)
  }
  
  return(results)
}

#Function to plot results of iterative tests
#Plot results
plot_results = function(res_list,
                        ylim = c(4,9),
                        ylab = NULL,
                        xlab = NULL,
                        plot_as_lines = FALSE,
                        col1 = 'grey50',
                        col2 = 'grey80',
                        ...){
  
  means = unlist(lapply(res_list, function(x) mean(x)))
  error = lapply(res_list, function(x) boxplot.stats(x)$stats)
  
  if(plot_as_lines == TRUE){
    
    means = unlist(lapply(res_list, function(x) median(x)))
    plot(means,
         xaxt = 'n',
         cex.axis = 1.5,
         cex.lab = 1.5,
         ylab = ylab,
         ylim = ylim,
         bty = 'n',
         las = 2,
         type = 'n',
         lwd = 3,
         col = col1,
         xlab = xlab,
         ...)
    polygon(c(seq(1, length(means),1),
              rev(seq(1, length(means),1))),
            c(lapply(error, function(x) x[5]),
              rev(lapply(error, function(x) x[1]))),
            col = col2, border = FALSE)
    lines(means,
          type = "l", 
          col = col1, 
          lwd = 3)
    axis(1, 
         at = seq(1, length(means),1),
         labels = names(means),
         cex.lab = 1.5,
         cex.axis = 1.5,
         las = 2)
  }else{
    plot(means,
         xaxt = 'n',
         cex.axis = 1.5,
         cex.lab = 1.5,
         ylab = ylab,
         ylim = ylim,
         pch = 21,
         cex = 2,
         bty = 'n',
         las = 2,
         type = 'n',
         bg = col2,
         col = col1,
         xlab = xlab,
         ...)
    axis(1, 
         at = seq(1, length(means),1),
         labels = names(means),
         cex.lab = 1.5,
         cex.axis = 1.5,
         las = 2)
    
    for(i in 1:length(res_list)){
      points(jitter(rep(i, length(res_list[[i]])), 0.25),
             res_list[[i]],
             pch = 20,
             col = alpha(col2, 0.5))}
    
    points(means,
           pch = 21,
           cex = 1.5,
           bg = col1,
           col = col1)
  }
  
  l = list(means, error)
  names(l) = c('means', 'error')
  return(l)
}

#Function to plot results of iterative tests as normalized variance
plot_variance = function(res_list,
                         ylim = c(0,1),
                         ylab = NULL,
                         xlab = NULL,
                         return = FALSE,
                         ...){
  
  variance = unlist(lapply(res_list, function(x) sd(x)/mean(x)))
  
  plot(variance,
       xaxt = 'n',
       cex.axis = 1.5,
       cex.lab = 1.5,
       ylab = ylab,
       ylim = ylim,
       pch = 21,
       cex = 2,
       bty = 'n',
       las = 2,
       type = 'n',
       bg = 'grey80',
       col = 'grey60',
       xlab = xlab,
       ...)
  axis(1, 
       at = seq(1, length(variance),1),
       labels = names(variance),
       cex.lab = 1.5,
       cex.axis = 1.5,
       las = 2)
  
  points(variance,
         pch = 21,
         cex = 1.5,
         bg = 'grey40',
         col = 'grey40')
  
  if(return == TRUE){
    return(variance)
  }
}

#Function to plot recurrence results
plot_recurrence = function(recurrences,
                           mar = c(0.5, 0.5, 0.5, 0.5)){
  
  #Analyze distribution of recurrences
  par(mfrow = c(length(recurrences), 1))
  par(mar = mar)
  
  for(i in 1:length(recurrences)){
    image(do.call(cbind, lapply(recurrences[[i]], function(y) y$proportion_recurrent_in_bins)),
          col = colorRampPalette(hcl.colors(12, "YlOrRd", rev = TRUE))(100),
          xaxt = 'n',
          yaxt = 'n')
    title(main = paste(names(recurrences)[i], 'frames'),
          cex.main = 1.5,
          font.main = 1)}
  
  axis(1,
       at = seq(0, 1, 0.125),
       labels = seq(0, max(as.numeric(names(recurrences))), max(as.numeric(names(recurrences)))/8),
       cex.axis = 1.5)
}

#Function to plot as vector field
plot_vector_field = function(layout,
                             bin_umap = FALSE,
                             n_bins = 32,
                             color_by_theta = FALSE,
                             arrow_color = 'grey50',
                             arrow_length = 0.05,
                             return = FALSE){
  
  if(bin_umap == TRUE){
    layout = bin_umap(layout, n_bins = n_bins)$layout
  }
  
  layout$dx = c(0, diff(layout$x))
  layout$dy = c(0, diff(layout$y))
  
  bins = split(layout, layout$xy)
  dx_mean = lapply(bins, function(x) mean(x$dx))
  dy_mean = lapply(bins, function(x) mean(x$dy))
  
  df = data.frame(x = as.numeric(unlist(lapply(strsplit(names(dx_mean), "_"), function(v){v[1]}))),
                  y = as.numeric(unlist(lapply(strsplit(names(dx_mean), "_"), function(v){v[2]}))),
                  dx = unlist(dx_mean),
                  dy = unlist(dy_mean))
  df$theta = rep(NA, nrow(df))
  for(i in 1:nrow(df)){
    x1 = df[i,1] 
    y1 = df[i,2] 
    x2 = df[i,1] + df[i,3]
    y2 = df[i,2] + df[i,4]
    
    df$theta[i] = atan2(y2-y1, x2-x1)*(180/pi)
  }
  
  df$dist = rep(NA, nrow(df))
  for(i in 1:nrow(df)){
    x1 = df[i,1] 
    y1 = df[i,2] 
    x2 = df[i,1] + df[i,3]
    y2 = df[i,2] + df[i,4]
    
    df$dist[i] = euc.dist(c(x1, y1), c(x2, y2))
  }
  
  par(mar = c(1,1,1,1))
  
  if(color_by_theta == TRUE){
    x = round(df$dy, 2)
    cols = colorRampPalette(c('cyan4', 'grey90', 'orangered3'))(length(seq(min(x), max(x), 0.01)))
    names(cols) = round(seq(min(x), max(x), 0.01), 2)
    cols = cols[match(as.numeric(x), 
                      as.numeric(names(cols)))]
    
    theta = round(df$theta)
    s = seq(-180, 180, 1)
    cols = c(colorRampPalette(c('midnightblue', 'cyan4'))(90),
             colorRampPalette(c('cyan4', 'lightgoldenrod1'))(90),
             colorRampPalette(c('lightgoldenrod1', 'sienna2'))(90),
             colorRampPalette(c('sienna2', 'orangered3'))(91))
    names(cols) = s
    cols = cols[match(theta, names(cols))]
    
    plot(df$x, 
         df$y, 
         xlim = c(min(df$x)-2, max(df$x)+2),
         ylim = c(min(df$y)-2, max(df$y)+2),
         type = "n",
         #col = alpha('grey50', 0.5),
         pch = 20,
         xlab = "",
         ylab = "",
         bty = 'n',
         xaxt = 'n',
         yaxt = 'n')
    shape::Arrows(df[,1], 
                  df[,2], 
                  df[,1] + df[,3]/2, 
                  df[,2] + df[,4]/2,
                  arr.length = arrow_length,
                  col = cols,
                  arr.type = "triangle")
  }else{
    
    plot(df$x, 
         df$y, 
         xlim = c(min(df$x)-2, max(df$x)+2),
         ylim = c(min(df$y)-2, max(df$y)+2),
         type = "n",
         #col = alpha('grey50', 0.5),
         pch = 20,
         xlab = "",
         ylab = "",
         bty = 'n',
         xaxt = 'n',
         yaxt = 'n')
    shape::Arrows(df[,1], 
                  df[,2], 
                  df[,1] + df[,3]/2, 
                  df[,2] + df[,4]/2,
                  arr.length = arrow_length,
                  col = arrow_color,
                  arr.type = "triangle")
  }
  
  if(return == TRUE){
    return(df)
  }
}

#Function to plot with features colored
plot_umap_features = function(layout,
                              windows,
                              bin_umap = FALSE,
                              n_bins = 32,
                              n_features = NULL,
                              feature_names = NULL,
                              colors = brewer.pal(11, 'Spectral'),
                              plot_points = FALSE,
                              return = FALSE,
                              ...){
  
  #Bin UMAP if desired
  if(bin_umap == TRUE){
    layout = bin_umap(layout,
                      n_bins = n_bins)$layout
  }
  
  #Get vector of rows to split windows on (as a function of feature number)
  tosplit = rep(1:n_features, 
                each=(nrow(windows)/n_features))
  
  #Split windows on features
  feat = split(as.data.frame(windows), tosplit)
  
  #Get colors
  cols = colorRampPalette(colors)(n_features)
  
  #Set up plotting aesthetics
  par(mfrow = c(1, n_features), bty = 'n', xaxt = 'n', yaxt = 'n', mar = c(2,2,2,2))
  
  #Loop through features and plot
  for(i in 1:length(feat)){
    
    #Calculate mean feature value per window
    m = colMeans(feat[[i]])
    
    #Match to layout
    m = m[1:nrow(layout)]
    
    #Add to layout
    if(is.null(feature_names) == FALSE){
      layout[,feature_names[i]] = m
    }else{
      layout = cbind(layout, m)
    }
    
    if(plot_points == TRUE){
      
      #Round mean feature value
      m = round(m, 2)
      
      #Take absolute value
      m = abs(m)
      
      #Get colors
      p = colorRampPalette(c('grey90', cols[i]))(length(seq(0, max(m), 0.01)))
      names(p) = seq(0, max(m), 0.01)
      p = p[match(m, names(p))]
      
      #Plot
      plot(layout[,1:2],
           pch = 20,
           col = p,
           ylab = '',
           xlab = '',
           ...)
      
      if(is.null(feature_names) == FALSE){
        title(main = feature_names[i],
              cex.main = 1.5,
              font.main = 1)}
      
    }else{
      
      #Split on bin
      m = split(m, layout$xy_new)
      
      #Get mean per bin
      m = lapply(m, function(x) mean(x, na.rm = TRUE))
      
      #Unlist
      m = unlist(m)
      
      #Round mean feature value
      m = round(m, 2)
      
      #Take absolute value
      m = abs(m)
      
      #Get colors
      p = colorRampPalette(c('grey90', cols[i]))(length(seq(0, max(m), 0.01)))
      names(p) = seq(0, max(m), 0.01)
      p = p[match(m, names(p))]
      
      #Plot
      plot(unlist(lapply(strsplit(names(m), "_"), function(v){v[1]})),
           unlist(lapply(strsplit(names(m), "_"), function(v){v[2]})),
           pch = 20,
           col = p,
           ylab = '',
           xlab = '',
           ...)
      
      if(is.null(feature_names) == FALSE){
        title(main = feature_names[i],
              cex.main = 1.5,
              font.main = 1)}
    }
  }
  if(return == TRUE){
    return(layout)
  }
}

#Function plot as a probability density function
plot_umap_pdf = function(layout,
                         h = 1,
                         n = 100, 
                         colors = matlab.like(100),
                         return = FALSE){
  
  #Get pdf
  pdf = kde2d(layout$x,
              layout$y,
              h = h,
              n = n)
  
  #Plot
  par(mar = c(1,1,1,1), bty = 'n', xaxt = 'n', yaxt = 'n')
  image(pdf$z,
        xlab = '',
        ylab = '',
        col = colors)
  
  #Return if desired
  if(return == TRUE){
    return(pdf)
  }
}

#Function to compare UMAP distributions across trials/individuals via bin-wise Fisher's test
run_umap_fishers = function(layout,
                            individuals_vector,
                            bin_umap = FALSE,
                            n_bins = 32,
                            odds_cutoff = 2,
                            cex = 0.5,
                            adjust_ps = FALSE,
                            verbose = FALSE,
                            return = FALSE){
  
  #Bin layout if desired
  if(bin_umap == TRUE){
    layout = bin_umap(layout,
                      n_bins = n_bins)$layout
  }
  
  #Split layout
  individuals = split(layout,
                      individuals_vector)
  
  #Set up plots
  n = length(individuals)*2
  x = ceiling(sqrt(n))
  y = floor(sqrt(n))
  par(mfrow = c(x,y), mar = c(2,2,2,2))
  rm(n, x, y)
  
  #Set up empty lists to save results
  all_odds = list()
  all_ps = list()
  
  #Loop through and run test on each individual
  for(h in 1:length(individuals)){
    
    if(verbose == TRUE){
      print(paste('individual', h, 'out of', length(individuals)))
    }
    r = individuals[[h]]$xy_new
    t = table(r)
    xy_table = table(layout$xy_new)
    
    t = t[match(names(xy_table), names(t))]
    names(t) = names(xy_table)
    t[is.na(t)] = 0
    
    odds = c()
    ps = c()
    
    for(i in 1:length(t)){
      w1 = t[i]
      x1 = xy_table[grep(paste("^", names(t[i]), "$", sep = ""), names(xy_table))]
      w2 = sum(t)-w1
      x2 = sum(xy_table)-x1
      
      out = fisher.test(as.matrix(
        rbind(
          c(w1, w2),
          c(x1, x2))))
      
      odds = c(odds, out$estimate)
      ps = c(ps, out$p.value)
    }
    
    #Adjust ps if desired
    if(adjust_ps == TRUE){
      ps = ps*length(ps)
    }
    
    #Add to list
    all_odds[[as.character(h)]] = odds
    all_ps[[as.character(h)]] = ps
    
    ##Plot odds ratios
    #Round and change names
    o = round(odds, 2)
    names(o) = names(t)
    
    #Get colors
    ints = seq(0, odds_cutoff, 0.01)
    o[o>odds_cutoff] = odds_cutoff
    cols = c(colorRampPalette(c("midnightblue", "grey90"))(length(seq(0, 1, 0.01))),
             colorRampPalette(c("grey90", "darkred"))(length(seq(1.01, odds_cutoff, 0.01))))
    names(cols) = ints
    cols = cols[match(o, names(cols))]
    
    #Plot
    plot(unlist(lapply(strsplit(names(o), "_"), function(v){v[1]})),
         unlist(lapply(strsplit(names(o), "_"), function(v){v[2]})),
         pch = 20,
         cex = cex,
         col = cols,
         cex.axis = 1.5,
         cex.lab = 1.5,
         xaxt = 'n',
         yaxt = 'n',
         ylab = "",
         xlab = "",
         bty = 'n')
    title(main = 'Odds ratios',
          cex.main = 1.5,
          font.main = 1)
    
    ##Plot p values
    #Set up colors
    cols = rep('grey90', length(ps))
    cols[ps<0.05] = 'red'
    
    #Plot
    plot(unlist(lapply(strsplit(names(o), "_"), function(v){v[1]})),
         unlist(lapply(strsplit(names(o), "_"), function(v){v[2]})),
         pch = 20,
         cex = cex,
         col = cols,
         cex.axis = 1.5,
         cex.lab = 1.5,
         xaxt = 'n',
         yaxt = 'n',
         ylab = "",
         xlab = "",
         bty = 'n')
    title(main = 'p-values',
          cex.main = 1.5,
          font.main = 1)
  }
  
  #Return if desired
  if(return == TRUE){
    l = list(all_odds, all_ps)
    names(l) = c('odds_ratios', 'p_values')
    return(l)
  }
}

#Function to plot markov on umap
plot_umap_markov = function(layout,
                            bin_umap = FALSE,
                            n_bins = 16,
                            chord,
                            plot_umap_points = TRUE,
                            plot_self = FALSE,
                            plot_dynamic_points = FALSE,
                            lwd_scale = 50,
                            curve = 0.2,
                            arr.width = 0,
                            point_col = NULL,
                            arr.length = 0){
  
  #Bin if desired
  if(bin_umap == TRUE){
    layout = bin_umap(layout, n_bins = n_bins)$layout
  }
  
  #Get markov
  probs = fitHigherOrder(layout$xy_new, 2)$Q[[1]]
  probs[probs<0.1] = 0
  
  #Get chord
  chord = circlize::chordDiagram(probs, 
                                 directional = TRUE,
                                 #direction.type = "arrows",
                                 grid.col =  alpha('grey60', 0.5),
                                 annotationTrack = c("name", "grid"));
  
  #Plotting in umap
  #Reformat chord diagram coords
  fromx = as.numeric(unlist(lapply(strsplit(chord$rn, "_"), function(v){v[1]})))
  fromy = as.numeric(unlist(lapply(strsplit(chord$rn, "_"), function(v){v[2]})))
  tox = as.numeric(unlist(lapply(strsplit(chord$cn, "_"), function(v){v[1]})))
  toy = as.numeric(unlist(lapply(strsplit(chord$cn, "_"), function(v){v[2]})))
  
  chord$fromx = fromx
  chord$fromy = fromy
  chord$tox = tox
  chord$toy = toy
  
  #Split into self and not
  chord_s = chord[chord$rn == chord$cn,]
  chord = chord[!chord$rn == chord$cn,]
  
  #Make colors not alpha
  chord_s$col = substr(chord_s$col, 1, 7)
  chord$col = substr(chord$col, 1, 7)
  
  #Remove zeroes
  chord = chord[!chord$value1==0,]
  
  if(plot_umap_points == FALSE){
    plot(layout$xnew,
         layout$ynew,
         pch = 20,
         #col = "grey90",
         col = NULL,
         bty = 'n',
         xaxt = 'n',
         yaxt = 'n',
         xlab = "",
         ylab = "")
  }else{
    plot(layout$xnew,
         layout$ynew,
         pch = 20,
         col = point_col,
         #col = NULL,
         cex = 3,
         bty = 'n',
         xaxt = 'n',
         yaxt = 'n',
         xlab = "",
         ylab = "")
    
  }
  
  for(i in 1:nrow(chord)){
    diagram::curvedarrow(from = c(chord$fromx[i], chord$fromy[i]),
                         to = c(chord$tox[i], chord$toy[i]),
                         lcol = chord$col[i],
                         lwd = chord$value1[i]*lwd_scale,
                         arr.pos = 0.9,
                         curve = curve,
                         arr.type = "triangle",
                         arr.length = arr.length,
                         arr.width = arr.width
    )}
  
  if(plot_self == TRUE){
    for(i in 1:nrow(chord)){
      diagram::selfarrow(c(chord$fromx[i], chord$fromy[i]),
                         lcol = chord$col[i],
                         lwd = chord$value1[i]*lwd_scale,
                         curve = 0.2,
                         arr.length = 0,
                         arr.width = 0)
    }
  }
  
  if(plot_dynamic_points == TRUE){
    f = table(k$cluster)
    f = f/sum(f)
    
    for(i in 1:nrow(k$centers)){
      points(k$centers[i,1],
             k$centers[i,2],
             pch = 21,
             bg = cols[i],
             col = darken_color(cols[i]),
             cex = 4+(2*f[i]))
    }
  }else{
    points(k$centers,
           pch = 21,
           bg = cols,
           col = darken_color(cols),
           cex = 4)
  }
}