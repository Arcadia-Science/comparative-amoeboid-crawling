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
library(swaRm)
library(dunn.test)
library(vioplot)
library(lsa)
library(shape)
library(gtools)
library(entropy)
library(pracma)
library(fractaldim)

##########################
#####TREBLE functions#####
##########################
#Function to load cell image .jpegs
load_images = function(image_directory) {
  setwd(image_directory)
  
  #Get directories
  files = list.files()
  
  dat = list()
  names = c()
  for (i in 1:length(files)) {
    setwd(paste(files[i], '/', 'masks_normalized/', sep = ''))
    f = list.files(full.names = FALSE)
    f = mixedsort(sort(f))
    
    imgs = list()
    for (j in 1:length(f)) {
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

#Function to extract windows of user defined features, size x and stepsize y
get_windows = function(features,
                       window_size = 10,
                       step_size = 1,
                       name = NULL) {
  #Get windows
  peaks1 = splitWithOverlap(features, window_size, window_size - step_size)
  
  #Clean and combine into matrix
  peaks1 = peaks1[1:(length(peaks1) - window_size)]
  peaks1 = lapply(peaks1, function(x)
    unlist(as.data.frame(x)))
  peaks1 = do.call(cbind, peaks1)
  
  if (is.null(name) == FALSE) {
    colnames(peaks1) = paste(name, '_', rownames(features)[1:ncol(peaks1)], sep = '')
  } else{
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
                                name = NULL) {
  #Initialize window list
  frags = list()
  xys = list()
  
  #Get trajectories
  for (i in seq(1, (nrow(features) - window_size), step_size)) {
    if (verbose == TRUE) {
      if (i %% 10000 == TRUE) {
        print(i)
      }
    }
    
    #Get velocity vectors
    vt = features$translational_velocity[seq(i, i + window_size, 1)]
    vr = features$angular_velocity[seq(i, i + window_size, 1)]
    if (include_sideslip == TRUE) {
      vs = features$sideslip[seq(i, i + window_size, 1)]
    }
    x = features$x[seq(i, i + window_size, 1)]
    y = features$y[seq(i, i + window_size, 1)]
    
    time = features$time[seq(i, i + window_size, 1)]
    
    #Subtract first frame make t0 = 0
    vr = vr - vr[1]
    if (include_sideslip == TRUE) {
      vs = vs - vs[1]
    }
    
    #Multiply by sign of second frame to normalize turn direction
    if (symm == TRUE) {
      if (vr[2] < 0) {
        vr = vr * (-1)
      }
      
      if (include_sideslip == TRUE) {
        if (vs[2] < 0) {
          vs = vs * (-1)
        }
      }
      
      vr = abs(vr)
      
      if (include_sideslip == TRUE) {
        vs = abs(vs)
      }
    }
    
    #Look at time difference
    t_diff = diff(time)
    
    if (include_sideslip == TRUE) {
      frags[[paste(time[1], "_", time[length(time)], "_", name, sep = "")]] = c(vt, vr, vs)
      xys[[paste(time[1], "_", time[length(time)], "_", name, sep = "")]] = c(x, y)
    } else{
      frags[[paste(time[1], "_", time[length(time)], "_", name, sep = "")]] = c(vt, vr)
      xys[[paste(time[1], "_", time[length(time)], "_", name, sep = "")]] = c(x, y)
    }
  }
  
  #Remove elements with NAs
  frags = frags[lapply(frags, function(x)
    sum(is.na(x))) < 1]
  xys = xys[lapply(xys, function(x)
    sum(is.na(x))) < 1]
  
  #Combine into df
  if (return_xy_windows == TRUE) {
    df = do.call(cbind, frags)
    xys_df = do.call(cbind, xys)
    
    l = list(df, xys_df)
    names(l) = c("velocity_windows", "xy_windows")
    return(l)
  } else{
    df = do.call(cbind, frags)
    return(df)
  }
}

#Function to bin a umap space into a n x n grid (umap coordinates are provided as the 'layout' object)
bin_umap = function(layout,
                    n_bins) {
  n = n_bins
  
  #Split x into n bins
  x1 = seq(min(layout[, 1]),
           max(layout[, 1]),
           (max(layout[, 1]) - min(layout[, 1])) / n)
  names(x1) = seq(1, n + 1, 1)
  
  xnew = apply(layout, 1, function(x)
    names(x1)[which.min(abs(as.numeric(x[1]) - x1))])
  layout$xnew = xnew
  
  #Split y into n bins
  y1 = seq(min(layout[, 2]),
           max(layout[, 2]),
           (max(layout[, 2]) - min(layout[, 2])) / n)
  names(y1) = seq(1, n + 1, 1)
  
  ynew = apply(layout, 1, function(x)
    names(y1)[which.min(abs(as.numeric(x[2]) - y1))])
  layout$ynew = ynew
  
  #Paste xy to get unique bin combos (this will be input to sling shot as 'clusters')
  xy_new = paste(xnew, ynew, sep = "_")
  layout$xy_new = xy_new
  
  #Convert coordinates to numeric, sort first
  m = unique(xy_new)
  
  y = as.numeric(unlist(lapply(strsplit(m, "_"), function(v) {
    v[2]
  })))
  m = m[order(y)]
  
  x = as.numeric(unlist(lapply(strsplit(m, "_"), function(v) {
    v[1]
  })))
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
                          ...) {
  if (verbose == TRUE) {
    print("Getting windows")
  }
  #Get windows
  windows = list()
  
  #Set up plots
  if (plot == TRUE) {
    n = length(features) * 2
    x = ceiling(sqrt(n))
    y = floor(sqrt(n))
    par(mfrow = c(x, y), mar = c(1, 1, 1, 1))
    rm(n, x, y)
  }
  
  for (i in 1:length(features)) {
    if (velocity_windows == TRUE) {
      windows[[i]] = get_velocity_windows(features[[i]],
                                          window_size = window_size,
                                          step_size = step_size,
                                          ...)
    } else{
      windows[[i]] = get_windows(features[[i]],
                                 window_size = window_size,
                                 step_size = step_size,
                                 ...)
    }
  }
  
  if (run_umap == FALSE) {
    #Return
    l = list(features, windows)
    names(l) = c("features", "windows")
    return(l)
  } else{
    if (verbose == TRUE) {
      print("Running UMAP")
    }
    
    #Run UMAP
    umaps = list()
    for (i in 1:length(features)) {
      if (verbose == TRUE) {
        print(paste("umap", i, "out of", length(features)))
      }
      
      if (verbose == TRUE) {
        umaps[[i]] = umap(t(windows[[i]]), verbose = TRUE)
      } else{
        umaps[[i]] = umap(t(windows[[i]]))
      }
      
      if (plot == TRUE) {
        plot(
          umaps[[i]]$layout[, 1:2],
          bty = 'n',
          xaxt = 'n',
          yaxt = 'n',
          ylab = "",
          pch = 20,
          xlab = "",
          col = alpha('grey50', 0.5)
        )
        plot(
          umaps[[i]]$layout[, 1:2],
          type = 'l',
          bty = 'n',
          xaxt = 'n',
          yaxt = 'n',
          ylab = "",
          xlab = "",
          col = alpha('grey50', 0.5)
        )
      }
    }
    
    #Extract layouts
    umaps = lapply(umaps, function(z)
      data.frame(x = z$layout[, 1], y = z$layout[, 2]))
    
    #Bin
    umaps = lapply(umaps, function(x)
      bin_umap(x, n_bins = n_bins)$layout)
    
    #Return
    if (return_windows == TRUE) {
      l = list(vel, windows, umaps)
      names(l) = c("features", "windows", "umaps")
    } else{
      l = list(umaps)
      names(l) = 'umaps'
    }
    return(l)
  }
}

#Function to calculate distance between umap layouts using Procrustes and Euclidean distances
run_procrustes = function(umaps,
                          run_protest = FALSE) {
  #Get all combinations of umaps to compare
  x = combn(seq(1, length(umaps), 1), 2)
  
  #Run
  pr_res = c()
  pr_sig = list()
  dists = c()
  for (i in 1:ncol(x)) {
    pr = procrustes(umaps[[x[1, i]]][, 1:2],
                    umaps[[x[2, i]]][, 1:2])
    pr_res = c(pr_res, summary(pr)$rmse)
    
    dists = c(dists, euc.dist(umaps[[x[1, i]]][, 1:2],
                              umaps[[x[2, i]]][, 1:2]))
    
    if (run_protest == TRUE) {
      pr_sig[[i]] = protest(umaps[[x[1, i]]][, 1:2],
                            umaps[[x[2, i]]][, 1:2])
    }
  }
  
  if (run_protest == TRUE) {
    res = list(pr_res, pr_sig)
    names(res) = c("procrustes", "protest")
  } else{
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
                                threshold = 0.05) {
  results = list()
  for (h in 1:length(umaps)) {
    if (verbose == TRUE) {
      print(paste(h, "out of", length(umaps)))
    }
    
    u = umaps[[h]]
    
    if (filter_outliers == TRUE) {
      u$x[u$x > 30] = 30
      u$x[u$x < (-30)] = -30
      u$y[u$y > 30] = 30
      u$y[u$y < (-30)] = -30
    }
    
    #Get distances
    l = bin_umap(u,
                 n_bins = n_bins)$layout
    res = list()
    pos = unique(l$xy_new)
    
    dists = list()
    for (i in 1:length(pos)) {
      x = c(as.numeric(unlist(lapply(strsplit(pos[i], "_"), function(v) {
        v[1]
      }))),
      as.numeric(unlist(lapply(strsplit(pos[i], "_"), function(v) {
        v[2]
      }))))
      z = apply(l, 1, function(y)
        euc.dist(x, c(as.numeric(y[3]), as.numeric(y[4]))))
      dists[[pos[i]]] = z
    }
    
    #Calculate distance distribution
    thresh = quantile(unlist(dists), probs = threshold)
    
    #Extract recurrences usins 10% threshold
    recs = lapply(dists, function(x) {
      rs = which(x < thresh)
      ds = diff(rs)
      ds[ds > thresh]
    })
    
    if (plot == TRUE) {
      histogram = hist(unlist(recs),
                       breaks = seq(1, max(unlist(recs), na.rm = TRUE), 1),
                       xlim = c(0, 200))
    } else{
      histogram = hist(unlist(recs),
                       breaks = seq(1, max(unlist(recs), na.rm = TRUE), 1),
                       plot = FALSE)
    }
    
    #Label points recurrent in 1 second bins
    prop_recurrent = list()
    for (i in 1:200) {
      z = lapply(recs, function(x)
        which(x == i))
      z = sum(lapply(z, function(x)
        length(x)) > 0) / length(recs)
      prop_recurrent[[i]] = z
    }
    
    #barplot(unlist(prop_recurrent))
    total_recurrent = sum(lapply(recs, function(x)
      length(x)) > 0) / length(recs)
    
    l = list(dists,
             unlist(recs),
             histogram,
             unlist(prop_recurrent),
             total_recurrent)
    names(l) = c(
      "distances",
      "recurrences",
      "histogram",
      "proportion_recurrent_in_bins",
      "total_proportion_recurrent"
    )
    results[[h]] = l
  }
  
  if (plot == TRUE) {
    image(
      do.call(
        cbind,
        lapply(results, function(x)
          x$histogram$counts[1:200] / max(x$histogram$counts[1:200]))
      ),
      col = matlab.like(32),
      xaxt = 'n',
      yaxt = 'n',
      xlab = 'Time',
      ylab = 'Replicate',
      cex.lab = 1.5
    )
    axis(
      1,
      at = seq(0, 1, 1 / (200 - 1)),
      labels = seq(1, 200, 1),
      cex.axis = 1.5
    )
  }
  
  return(results)
}

#Function to calculate per bin entropies
calculate_bin_entropies = function(bin,
                                   strain_coords,
                                   strain_xy,
                                   window = 30,
                                   n_bins,
                                   plot_trajectories = FALSE,
                                   calculate_entropies = FALSE,
                                   n_trajectories = NULL) {
  #Choose bin
  z = bin
  tmp = strain_coords
  
  #Split on bin and extract bouts in between
  y = which(!tmp == z)
  idx <- c(0, cumsum(abs(diff(y)) > 1))
  indices = split(y, idx)
  
  print("Extracting inter-bin bouts")
  #Get bin name
  series = lapply(indices, function(x)
    tmp[x])
  
  #Remove first bout (doesn't originate at bin of interest)
  series = series[-1]
  
  #Require bout of length n
  series = series[lapply(series, length) >= window]
  series = lapply(series, function(x)
    x[1:window])
  
  #Make empty vector for saving output
  entropies = c()
  
  if (length(series) > 1) {
    #Combine into dataframe
    series = do.call(cbind, series)
    print(ncol(series))
    
    #Add row corresponding to starting bin
    series = rbind(rep(z, ncol(series)), series)
    
    #Plot if desired
    if (plot_trajectories == TRUE) {
      plot(
        series[, 1],
        type = "l",
        lwd = 2,
        col = alpha("grey40", 0.5),
        ylim = c(min(unlist(series)), max(unlist(series))),
        ylab = "Bin",
        xlab = "Time",
        cex.axis = 1.5,
        cex.lab = 1.5,
        bty = 'n'
      )
      for (i in 2:ncol(series)) {
        lines(series[, i],
              lwd = 2,
              col = alpha("grey40", 0.5))
      }
    }
    
    if (!is.null(n_trajectories)) {
      if (ncol(series) > n_trajectories) {
        print("Calculating entropies")
        
        series = series[, sample(seq(1, ncol(series), 1), n_trajectories)]
        print(ncol(series))
        
        ##Calculate probability density functions
        for (h in 1:nrow(series)) {
          row = unlist(series[h, ])
          
          #Convert to coords
          row = cbind(as.numeric(unlist(lapply(strsplit(apply(as.data.frame(row), 1, function(x)
            unique(strain_xy[strain_coords == x])), "_"), function(v) {
              v[1]
            }))),
            as.numeric(unlist(lapply(strsplit(apply(as.data.frame(row), 1, function(x)
              unique(strain_xy[strain_coords == x])), "_"), function(v) {
                v[2]
              }))))
          
          w = kde2d(
            row[, 1],
            row[, 2],
            h = c(1, 1),
            n = n_bins,
            lims = c(1, n_bins + 1, 1, n_bins + 1)
          )$z
          
          pdf = as.numeric(unlist(as.data.frame(w)))
          entropies = c(entropies, entropy(pdf, unit = 'log2'))
        }
        
      }
    }
    
    if (calculate_entropies == TRUE) {
      print("Calculating entropies")
      ##Calculate probability density functions
      for (h in 1:nrow(series)) {
        row = unlist(series[h, ])
        
        #Convert to coords
        row = cbind(as.numeric(unlist(lapply(strsplit(apply(as.data.frame(row), 1, function(x)
          unique(strain_xy[strain_coords == x])), "_"), function(v) {
            v[1]
          }))),
          as.numeric(unlist(lapply(strsplit(apply(as.data.frame(row), 1, function(x)
            unique(strain_xy[strain_coords == x])), "_"), function(v) {
              v[2]
            }))))
        
        w = kde2d(
          row[, 1],
          row[, 2],
          h = c(1, 1),
          n = n_bins,
          lims = c(1, n_bins + 1, 1, n_bins + 1)
        )$z
        
        pdf = as.numeric(unlist(as.data.frame(w)))
        entropies = c(entropies, entropy(pdf, unit = 'log2'))
      }
    }
    
    l = list(series, entropies)
    names(l) = c("bouts", "entropies")
    return(l)
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
                            return = FALSE) {
  #Bin layout if desired
  if (bin_umap == TRUE) {
    layout = bin_umap(layout,
                      n_bins = n_bins)$layout
  }
  
  #Split layout
  individuals = split(layout,
                      individuals_vector)
  
  #Set up plots
  n = length(individuals) * 2
  x = ceiling(sqrt(n))
  y = floor(sqrt(n))
  par(mfrow = c(x, y), mar = c(2, 2, 2, 2))
  rm(n, x, y)
  
  #Set up empty lists to save results
  all_odds = list()
  all_ps = list()
  
  #Loop through and run test on each individual
  for (h in 1:length(individuals)) {
    if (verbose == TRUE) {
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
    
    for (i in 1:length(t)) {
      w1 = t[i]
      x1 = xy_table[grep(paste("^", names(t[i]), "$", sep = ""), names(xy_table))]
      w2 = sum(t) - w1
      x2 = sum(xy_table) - x1
      
      out = fisher.test(as.matrix(rbind(c(w1, w2),
                                        c(x1, x2))))
      
      odds = c(odds, out$estimate)
      ps = c(ps, out$p.value)
    }
    
    #Adjust ps if desired
    if (adjust_ps == TRUE) {
      ps = ps * length(ps)
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
    o[o > odds_cutoff] = odds_cutoff
    cols = c(colorRampPalette(c("midnightblue", "grey90"))(length(seq(0, 1, 0.01))),
             colorRampPalette(c("grey90", "darkred"))(length(seq(
               1.01, odds_cutoff, 0.01
             ))))
    names(cols) = ints
    cols = cols[match(o, names(cols))]
    
    #Plot
    plot(
      unlist(lapply(strsplit(names(
        o
      ), "_"), function(v) {
        v[1]
      })),
      unlist(lapply(strsplit(names(
        o
      ), "_"), function(v) {
        v[2]
      })),
      pch = 20,
      cex = cex,
      col = cols,
      cex.axis = 1.5,
      cex.lab = 1.5,
      xaxt = 'n',
      yaxt = 'n',
      ylab = "",
      xlab = "",
      bty = 'n'
    )
    title(main = 'Odds ratios',
          cex.main = 1.5,
          font.main = 1)
    
    ##Plot p values
    #Set up colors
    cols = rep('grey90', length(ps))
    cols[ps < 0.05] = 'red'
    
    #Plot
    plot(
      unlist(lapply(strsplit(names(
        o
      ), "_"), function(v) {
        v[1]
      })),
      unlist(lapply(strsplit(names(
        o
      ), "_"), function(v) {
        v[2]
      })),
      pch = 20,
      cex = cex,
      col = cols,
      cex.axis = 1.5,
      cex.lab = 1.5,
      xaxt = 'n',
      yaxt = 'n',
      ylab = "",
      xlab = "",
      bty = 'n'
    )
    title(main = 'p-values',
          cex.main = 1.5,
          font.main = 1)
  }
  
  #Return if desired
  if (return == TRUE) {
    l = list(all_odds, all_ps)
    names(l) = c('odds_ratios', 'p_values')
    return(l)
  }
}
