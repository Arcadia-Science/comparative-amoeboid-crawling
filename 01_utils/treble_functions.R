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
library(gplots)

##########################
##### TREBLE functions#####
##########################
##' load_images': Loads processed cell images from imoto et al. 2021 into R and
# compiles into a list of trials
load_images <- function(image_directory) {
  # Set work directory to where images are contained
  setwd(image_directory)

  # Get directories (corresponding to cell types)
  files <- list.files()

  # Generate empty list for save images into
  dat <- list()

  # Generate empty vector for saving trial names into
  names <- c()

  # Loop through cell type directories
  for (i in 1:length(files)) {
    # Set working directory to cell type directory
    setwd(paste(files[i], "/", "masks_normalized/", sep = ""))

    # List files in directories
    f <- list.files(full.names = FALSE)

    # Sort to be numerically increasing
    f <- mixedsort(sort(f))

    # Create empty list to save images into
    imgs <- list()

    # Loop through all trial directories and load images using 'readJPEG'
    for (j in 1:length(f)) {
      imgs[[f[j]]] <- unlist(as.data.frame(readJPEG(f[j])))
    }

    # Add images to data list
    dat[[as.character(files[i])]] <- imgs

    # Add names to names vector
    names <- c(names, rep(files[i], length(imgs)))

    # Return to root directory
    setwd("../../")
  }

  # Return list of cell images and names
  l <- list(dat, names)
  names(l) <- c("images", "names")
  return(l)
}

##' get_windows': Extracts windows from a matrix of data features using a
# desired window size ('window_size') and step size ('step_size)
get_windows <- function(features,
                        window_size = 10,
                        step_size = 1,
                        name = NULL) {
  # Use 'split_with_overlap' to generate windows of desired size
  win <- split_with_overlap(features, window_size, window_size - step_size)

  # Remove windows that contain NAs (i.e. are timepoints at the end of the
  # matrix that are within the window size distance from the end)
  win <- win[1:(length(win) - window_size)]

  # Unlist and combine into a matrix
  win <- lapply(win, function(x) {
    unlist(as.data.frame(x))
  })
  win <- do.call(cbind, win)

  if (is.null(name) == FALSE) {
    # If no names are given, name columns with time points
    colnames(win) <- paste(name, "_", rownames(features)[1:ncol(win)], sep = "")
  } else {
    # If names are given, use as column names
    colnames(win) <- rownames(features)[1:ncol(win)]
  }

  # Return windows
  return(win)
}

##' bin_umap': Bins a umap space into a nxn gird (umap coordinates are to be contained in the 'layout' object)
bin_umap <- function(layout,
                     n_bins) {
  # Generate a vector of all possible binned x coordinates
  x1 <- seq(
    min(layout[, 1]),
    max(layout[, 1]),
    (max(layout[, 1]) - min(layout[, 1])) / n_bins
  )

  # Name the binned x coordinates
  names(x1) <- seq(1, n_bins + 1, 1)

  # Calculate distance of all x points to the binned x coordinates, select the binned x coordinate that is closest
  xnew <- apply(layout, 1, function(x) {
    names(x1)[which.min(abs(as.numeric(x[1]) - x1))]
  })

  # Add binned coordinates to the umap layout as 'xnew'
  layout$xnew <- xnew

  # Generate a vector of all possible binned y coordinates
  y1 <- seq(
    min(layout[, 2]),
    max(layout[, 2]),
    (max(layout[, 2]) - min(layout[, 2])) / n_bins
  )

  # Name the binned y coordinates
  names(y1) <- seq(1, n_bins + 1, 1)

  # Calculate distance of all y points to the binned y coordinates, select the
  # binned y coordinate that is closest
  ynew <- apply(layout, 1, function(x) {
    names(y1)[which.min(abs(as.numeric(x[2]) - y1))]
  })

  # Add binned coordinates to the umap layout as 'ynew'
  layout$ynew <- ynew

  # Paste xy to get unique bin combos
  xy_new <- paste(xnew, ynew, sep = "_")

  # Add to the umap layout as 'xy_new'
  layout$xy_new <- xy_new

  # Exctract all unique binned values
  m <- unique(xy_new)

  # Sort x and y values
  y <- as.numeric(unlist(lapply(strsplit(m, "_"), function(v) {
    v[2]
  })))
  m <- m[order(y)]

  x <- as.numeric(unlist(lapply(strsplit(m, "_"), function(v) {
    v[1]
  })))
  m <- m[order(x)]

  # Add numeric names to the sorted, binned xy coordinates (these will be used
  # as numeric representations of the binned coordinates)
  names(m) <- seq(1, length(m), 1)

  # Match to the real data
  names(xy_new) <- names(m)[match(xy_new, m)]

  # Add the unique numeric coordinates to the umap layout as 'coords'
  layout$coords <- as.numeric(names(xy_new))

  # Return
  l <- list(layout, xy_new)
  names(l) <- c("layout", "new_coords")

  return(l)
}

##' iterative_umap': Generates a behavior space(s) from a given window size
# Acts as the central function of the iterative window search procedure part
# of the TREBLE framework
# The function expects a list of data feature matrices where list elements
# corresponding to individual trials
# Can be run over a range of window sizes by looping over window size values
iterative_umap <- function(features,
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
  # Create empty list to save windows into
  windows <- list()

  # Set up plots
  if (plot == TRUE) {
    n <- length(features) * 2
    x <- ceiling(sqrt(n))
    y <- floor(sqrt(n))
    par(mfrow = c(x, y), mar = c(1, 1, 1, 1))
    rm(n, x, y)
  }

  # Loop through trials and extract windows of desired size
  for (i in 1:length(features)) {
    windows[[i]] <- get_windows(features[[i]],
      window_size = window_size,
      step_size = step_size,
      ...
    )
  }

  if (run_umap == FALSE) {
    # If UMAP is not being run, return the features and windows
    l <- list(features, windows)
    names(l) <- c("features", "windows")
    return(l)
  } else {
    if (verbose == TRUE) {
      print("Running UMAP")
    }

    # Generate empty list to save UMAP embeddings into
    umaps <- list()

    # Loop through and generate UMAP embeddings from feature windows
    for (i in 1:length(features)) {
      if (verbose == TRUE) {
        print(paste("umap", i, "out of", length(features)))
      }

      # Run UMAP
      if (verbose == TRUE) {
        umaps[[i]] <- umap(t(windows[[i]]), verbose = TRUE)
      } else {
        umaps[[i]] <- umap(t(windows[[i]]))
      }

      # If desired, plot behavior space e,bedding
      if (plot == TRUE) {
        # As points
        plot(
          umaps[[i]]$layout[, 1:2],
          bty = "n",
          xaxt = "n",
          yaxt = "n",
          ylab = "",
          pch = 20,
          xlab = "",
          col = alpha("grey50", 0.5)
        )
        # As lines
        plot(
          umaps[[i]]$layout[, 1:2],
          type = "l",
          bty = "n",
          xaxt = "n",
          yaxt = "n",
          ylab = "",
          xlab = "",
          col = alpha("grey50", 0.5)
        )
      }
    }

    # Extract layouts
    umaps <- lapply(umaps, function(z) {
      data.frame(x = z$layout[, 1], y = z$layout[, 2])
    })

    # Bin
    umaps <- lapply(umaps, function(x) {
      bin_umap(x, n_bins = n_bins)$layout
    })

    # Return
    if (return_windows == TRUE) {
      l <- list(vel, windows, umaps)
      names(l) <- c("features", "windows", "umaps")
    } else {
      l <- list(umaps)
      names(l) <- "umaps"
    }
    return(l)
  }
}

##' run_procrustes': Calculates Euclidean and Procrustest distances between
# multiple UMAP embeddings generated by 'iterative_umap'
# This function can be used as a diagnostic to identify the optimal window size
# to use for a given data set
# Lower Procrustes values reflect similar structure in the UMAP embeddings
# across trials
# Lower Euclidean distance values, along with decreased variance in distance,
# reflect 'smoother' UMAP embeddings
run_procrustes <- function(umaps,
                           run_protest = FALSE) {
  # Get all combinations of UMAP embeddings to compare
  x <- combn(seq(1, length(umaps), 1), 2)

  # Generate empty lists and vectors to save results
  pr_res <- c()
  pr_sig <- list()
  dists <- c()

  # Loop through all combinations of UMAP embeddings and calculate procrustes distance
  for (i in 1:ncol(x)) {
    # Calculate procrustes distance of embeddings
    pr <- procrustes(
      umaps[[x[1, i]]][, 1:2],
      umaps[[x[2, i]]][, 1:2]
    )

    # Extract results
    pr_res <- c(pr_res, summary(pr)$rmse)

    # Calculate euclidean distance of all points in the two UMAP embeddings
    dists <- c(dists, euc_dist(
      umaps[[x[1, i]]][, 1:2],
      umaps[[x[2, i]]][, 1:2]
    ))

    # If desired, run the 'protest' function to calculate the significance of
    # the procrustes distance calculations
    if (run_protest == TRUE) {
      pr_sig[[i]] <- protest(
        umaps[[x[1, i]]][, 1:2],
        umaps[[x[2, i]]][, 1:2]
      )
    }
  }

  # Return
  if (run_protest == TRUE) {
    res <- list(pr_res, pr_sig)
    names(res) <- c("procrustes", "protest")
  } else {
    res <- list(pr_res, dists)
    names(res) <- c("procrustes", "euclidean_distances")
  }
  return(res)
}

##' calculate_recurrence': Measures the recurrence times (time it takes to
# return to a specific point) for all points in a behavior space
# Useful for identifying the optimal window size to use for a given data set
# Expects a list of UMAP embeddings as input ('umaps'); is often run with the
# results of multiple runs of 'iterative_umap'
calculate_recurrence <- function(umaps,
                                 filter_outliers = FALSE,
                                 n_bins = 16,
                                 plot = FALSE,
                                 verbose = FALSE,
                                 threshold = 0.05) {
  # Generate empty to list to contain results
  results <- list()

  # Loop through and calculate recurrence for all UMAP embeddings
  for (h in 1:length(umaps)) {
    if (verbose == TRUE) {
      print(paste(h, "out of", length(umaps)))
    }

    # Extract embedding of interest
    u <- umaps[[h]]

    # Filter outliers based on xy coordinates, if desired
    if (filter_outliers == TRUE) {
      u$x[u$x > 30] <- 30
      u$x[u$x < (-30)] <- -30
      u$y[u$y > 30] <- 30
      u$y[u$y < (-30)] <- -30
    }

    # Bin umap to desired resolution
    l <- bin_umap(u,
      n_bins = n_bins
    )$layout

    # Generate empty lists to save results
    res <- list()
    dists <- list()

    # Extract all unique points in binned UMAP embedding (will be used for
    # calculating recurrence times)
    pos <- unique(l$xy_new)

    # Loop through and calculate all returns to each binned position
    for (i in 1:length(pos)) {
      x <- c(
        as.numeric(unlist(lapply(strsplit(pos[i], "_"), function(v) {
          v[1]
        }))),
        as.numeric(unlist(lapply(strsplit(pos[i], "_"), function(v) {
          v[2]
        })))
      )
      z <- apply(l, 1, function(y) {
        euc_dist(x, c(as.numeric(y[3]), as.numeric(y[4])))
      })
      dists[[pos[i]]] <- z
    }

    # Calculate distance distribution
    thresh <- quantile(unlist(dists), probs = threshold)

    # Extract recurrences usins 10% threshold
    recs <- lapply(dists, function(x) {
      rs <- which(x < thresh)
      ds <- diff(rs)
      ds[ds > thresh]
    })

    # Calculate histogram of recurrences and plot if desired
    if (plot == TRUE) {
      histogram <- hist(unlist(recs),
        breaks = seq(1, max(unlist(recs), na.rm = TRUE), 1),
        xlim = c(0, 200)
      )
    } else {
      histogram <- hist(unlist(recs),
        breaks = seq(1, max(unlist(recs), na.rm = TRUE), 1),
        plot = FALSE
      )
    }

    # Label points recurrent in 1 second bins
    prop_recurrent <- list()
    for (i in 1:200) {
      z <- lapply(recs, function(x) {
        which(x == i)
      })
      z <- sum(lapply(z, function(x) {
        length(x)
      }) > 0) / length(recs)
      prop_recurrent[[i]] <- z
    }

    # Calculate total amount of points that are recurrent at timepoint x
    total_recurrent <- sum(lapply(recs, function(x) {
      length(x)
    }) > 0) / length(recs)

    # Save into results list
    l <- list(
      dists,
      unlist(recs),
      histogram,
      unlist(prop_recurrent),
      total_recurrent
    )
    names(l) <- c(
      "distances",
      "recurrences",
      "histogram",
      "proportion_recurrent_in_bins",
      "total_proportion_recurrent"
    )
    results[[h]] <- l
  }

  # Plot if desired
  if (plot == TRUE) {
    image(
      do.call(
        cbind,
        lapply(results, function(x) {
          x$histogram$counts[1:200] / max(x$histogram$counts[1:200])
        })
      ),
      col = matlab.like(32),
      xaxt = "n",
      yaxt = "n",
      xlab = "Time",
      ylab = "Replicate",
      cex.lab = 1.5
    )
    axis(
      1,
      at = seq(0, 1, 1 / (200 - 1)),
      labels = seq(1, 200, 1),
      cex.axis = 1.5
    )
  }

  # Return
  return(results)
}

##' calculate_bin_entropies': Bins a UMAP embedding and calculates the overall
# entropy
calculate_bin_entropies <- function(bin,
                                    strain_coords,
                                    strain_xy,
                                    window = 30,
                                    n_bins,
                                    plot_trajectories = FALSE,
                                    calculate_entropies = FALSE,
                                    n_trajectories = NULL) {
  # Choose bin
  z <- bin
  tmp <- strain_coords

  # Split on bin and extract bouts in between
  y <- which(!tmp == z)
  idx <- c(0, cumsum(abs(diff(y)) > 1))
  indices <- split(y, idx)

  print("Extracting inter-bin bouts")
  # Get bin name
  series <- lapply(indices, function(x) {
    tmp[x]
  })

  # Remove first bout (doesn't originate at bin of interest)
  series <- series[-1]

  # Require bout of length n
  series <- series[lapply(series, length) >= window]
  series <- lapply(series, function(x) {
    x[1:window]
  })

  # Make empty vector for saving output
  entropies <- c()

  if (length(series) > 1) {
    # Combine into dataframe
    series <- do.call(cbind, series)
    print(ncol(series))

    # Add row corresponding to starting bin
    series <- rbind(rep(z, ncol(series)), series)

    # Plot if desired
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
        bty = "n"
      )
      for (i in 2:ncol(series)) {
        lines(series[, i],
          lwd = 2,
          col = alpha("grey40", 0.5)
        )
      }
    }

    if (!is.null(n_trajectories)) {
      if (ncol(series) > n_trajectories) {
        print("Calculating entropies")

        series <- series[, sample(seq(1, ncol(series), 1), n_trajectories)]
        print(ncol(series))

        ## Calculate probability density functions
        for (h in 1:nrow(series)) {
          row <- unlist(series[h, ])

          # Convert to coords
          row <- cbind(
            as.numeric(unlist(lapply(strsplit(apply(as.data.frame(row), 1, function(x) {
              unique(strain_xy[strain_coords == x])
            }), "_"), function(v) {
              v[1]
            }))),
            as.numeric(unlist(lapply(strsplit(apply(as.data.frame(row), 1, function(x) {
              unique(strain_xy[strain_coords == x])
            }), "_"), function(v) {
              v[2]
            })))
          )

          w <- kde2d(
            row[, 1],
            row[, 2],
            h = c(1, 1),
            n = n_bins,
            lims = c(1, n_bins + 1, 1, n_bins + 1)
          )$z

          pdf <- as.numeric(unlist(as.data.frame(w)))
          entropies <- c(entropies, entropy(pdf, unit = "log2"))
        }
      }
    }

    if (calculate_entropies == TRUE) {
      print("Calculating entropies")
      ## Calculate probability density functions
      for (h in 1:nrow(series)) {
        row <- unlist(series[h, ])

        # Convert to coords
        row <- cbind(
          as.numeric(unlist(lapply(strsplit(apply(as.data.frame(row), 1, function(x) {
            unique(strain_xy[strain_coords == x])
          }), "_"), function(v) {
            v[1]
          }))),
          as.numeric(unlist(lapply(strsplit(apply(as.data.frame(row), 1, function(x) {
            unique(strain_xy[strain_coords == x])
          }), "_"), function(v) {
            v[2]
          })))
        )

        w <- kde2d(
          row[, 1],
          row[, 2],
          h = c(1, 1),
          n = n_bins,
          lims = c(1, n_bins + 1, 1, n_bins + 1)
        )$z

        pdf <- as.numeric(unlist(as.data.frame(w)))
        entropies <- c(entropies, entropy(pdf, unit = "log2"))
      }
    }

    l <- list(series, entropies)
    names(l) <- c("bouts", "entropies")
    return(l)
  }
}
