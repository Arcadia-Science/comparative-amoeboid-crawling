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

############################
##### Plotting functions#####
############################
## 'plot_parameter': Function to plot parameters by behavior space position
# This is useful for diagnostic plotting to understand how behavior space is
# separating data (e.g. by coloring a space with an organism's average speed)
# 'plot_parameter' works by splitting the parameter of choice into a list
# corresponding to regions of behavior space (which are calculated with the
# function 'bin_umap' and are contained in 'layout$xy_new')
# Mean parameter values are calculated as a function of behavior space
# position and plotted using the 'image' function

plot_parameter <-
  function(parameter, layout, n_bins, return = FALSE, ...) {
    # Bin layout onto a discrete grid of nxn size
    # (determined by the 'n_bins' argument)
    layout <- bin_umap(layout, n_bins = n_bins)$layout

    # Split the parameter of choice based on layout position
    a <- parameter
    a <- split(a[1:nrow(layout)], layout$xy_new)

    # Calculate mean parameter values as a function of behavior space position
    a <- unlist(lapply(a, function(x) {
      mean(x, na.rm = TRUE)
    }))

    # Make a matrix containing all possible behavior space positions
    # (indexed as in 'xy_new' i.e. x coord, '_', y coord)
    p <- expand.grid(seq(1, n_bins + 1, 1), seq(1, n_bins + 1, 1))

    # Combine using an underscore to match the indices of 'xy_new'
    p <- paste(p[, 1], p[, 2], sep = "_")

    # Match the mean parameter value names (which correspond to 'xy_new') with
    # the full list of positions
    a <- a[match(p, names(a))]
    names(a) <- p

    # Plot using image, converting the mean parameter values into a matrix of
    #   size n_bins x n_bins
    image(
      matrix(a, nrow = n_bins + 1),
      xaxt = "n",
      yaxt = "n",
      xlab = "",
      ylab = "",
      bty = "n",
      ...
    )

    # Return parameter values if desired
    if (return == TRUE) {
      return(a)
    }
  }

## 'plot_results': Function to plot distributions of points as raw data with
# boxplot statistics, used extensively in the iterative window search analyses
# Takes a list of data 'res_list' wherein each entry corresponds to a different
# group of data
plot_results <- function(res_list,
                         ylim = c(0, 9),
                         ylab = NULL,
                         xlab = NULL,
                         median = FALSE,
                         plot_as_lines = FALSE,
                         col1 = "grey50",
                         col2 = "grey80",
                         ...) {
  # Calculate the mean values of all data groups
  if (median == TRUE) {
    means <- unlist(lapply(res_list, function(x) {
      median(x)
    }))
  } else {
    means <- unlist(lapply(res_list, function(x) {
      mean(x)
    }))
  }

  # Calculate boxplot statistics for all data groups
  error <- lapply(res_list, function(x) {
    boxplot.stats(x)$stats
  })

  # If 'plot_as_lines' is desired, plot the resulting means as connected lines
  # with error represented as polygons
  if (plot_as_lines == TRUE) {
    # Generate empty plot
    plot(
      means,
      xaxt = "n",
      cex.axis = 1.5,
      cex.lab = 1.5,
      ylab = ylab,
      ylim = ylim,
      bty = "n",
      las = 2,
      type = "n",
      lwd = 3,
      col = col1,
      xlab = xlab,
      ...
    )
    # Add polygons corresponding to error
    polygon(c(
      seq(1, length(means), 1),
      rev(seq(
        1, length(means), 1
      ))
    ),
    c(
      lapply(error, function(x) {
        x[5]
      }),
      rev(lapply(error, function(x) {
        x[1]
      }))
    ),
    col = col2,
    border = FALSE
    )
    # Add mean lines
    lines(means,
      type = "l",
      col = col1,
      lwd = 3
    )
    # Add axis
    axis(
      1,
      at = seq(1, length(means), 1),
      labels = names(means),
      cex.lab = 1.5,
      cex.axis = 1.5,
      las = 2
    )
  } else {
    # If not plotting as lines, plot as point distributions
    # Initiate empty plot
    plot(
      means,
      xaxt = "n",
      cex.axis = 1.5,
      cex.lab = 1.5,
      ylab = ylab,
      ylim = ylim,
      pch = 21,
      cex = 2,
      bty = "n",
      las = 2,
      type = "n",
      bg = col2,
      col = col1,
      xlab = xlab,
      ...
    )
    # Add axis
    axis(
      1,
      at = seq(1, length(means), 1),
      labels = names(means),
      cex.lab = 1.5,
      cex.axis = 1.5,
      las = 2
    )
    # Add data
    for (i in 1:length(res_list)) {
      points(jitter(rep(i, length(res_list[[i]])), 0.25),
        res_list[[i]],
        pch = 20,
        col = alpha(col2, 0.5)
      )
    }
    # Add mean points
    points(
      means,
      pch = 21,
      cex = 1.5,
      bg = col1,
      col = col1
    )
  }

  # Return results
  l <- list(means, error)
  names(l) <- c("means", "error")
  return(l)
}

## 'plot_variance': Like 'plot_results' but instead calculates and plots mean
# normalized variance for a list of data
plot_variance <- function(res_list,
                          ylim = c(0, 1),
                          ylab = NULL,
                          xlab = NULL,
                          return = FALSE,
                          ...) {
  # Calculate mean normalized variance
  variance <- unlist(lapply(res_list, function(x) {
    sd(x) / mean(x)
  }))

  # Initiate plot
  plot(
    variance,
    xaxt = "n",
    cex.axis = 1.5,
    cex.lab = 1.5,
    ylab = ylab,
    ylim = ylim,
    pch = 21,
    cex = 2,
    bty = "n",
    las = 2,
    type = "n",
    bg = "grey80",
    col = "grey60",
    xlab = xlab,
    ...
  )
  # Add axis
  axis(
    1,
    at = seq(1, length(variance), 1),
    labels = names(variance),
    cex.lab = 1.5,
    cex.axis = 1.5,
    las = 2
  )
  # Add points to plot
  points(
    variance,
    pch = 21,
    cex = 1.5,
    bg = "grey40",
    col = "grey40"
  )
  # Return variance if desired
  if (return == TRUE) {
    return(variance)
  }
}

## 'plot_reccurence': Plots the distribution of recurrence times between points
# in behavior spaces calculated from a variety of window sizes.
# Takes the output of the function 'calculate_recurrence' as input and plots a
# series of heatmaps corresponding to recurrence distributions
plot_recurrence <- function(recurrences,
                            mar = c(0.5, 0.5, 0.5, 0.5)) {
  # Initiate plot parameters
  par(mfrow = c(length(recurrences), 1))
  par(mar = mar)

  # Loop through and plot recurrence distributions for each window size
  for (i in 1:length(recurrences)) {
    # Plot heatmap of recurrence distributions using the 'image' function
    image(
      # cbind together the recurrence distributions for all individuals in the
      # data set
      do.call(
        cbind,
        lapply(recurrences[[i]], function(y) {
          y$proportion_recurrent_in_bins
        })
      ),
      # Set color values
      col = colorRampPalette(hcl.colors(12, "YlOrRd", rev = TRUE))(100),
      xaxt = "n",
      yaxt = "n"
    )
    # Add title
    title(
      main = paste(names(recurrences)[i], "frames"),
      cex.main = 1.5,
      font.main = 1
    )
  }
  # Add axis
  axis(
    1,
    at = seq(0, 1, 0.125),
    labels = seq(0, max(as.numeric(
      names(recurrences)
    )), max(as.numeric(
      names(recurrences)
    )) / 8),
    cex.axis = 1.5
  )
}

## 'plot_vector_field': Calculates the mean positional vectors originating from
# each grid in a binned behavior space and plots
plot_vector_field <- function(layout,
                              bin_umap = FALSE,
                              n_bins = 32,
                              color_by_theta = FALSE,
                              arrow_color = "grey50",
                              arrow_length = 0.05,
                              return = FALSE) {
  # If desired, bin the behavior space to a given resolution
  # (denoted by 'n_bins')
  if (bin_umap == TRUE) {
    layout <- bin_umap(layout, n_bins = n_bins)$layout
  }

  # Calculate the change in x and y coordinates for all points in time contained
  # in the behavior space
  layout$dx <- c(0, diff(layout$x))
  layout$dy <- c(0, diff(layout$y))

  # Split the behavior space layout on unique xy coordinates
  bins <- split(layout, layout$xy)

  # Calculate mean x and y coordinate changes for each bin
  dx_mean <- lapply(bins, function(x) {
    mean(x$dx)
  })
  dy_mean <- lapply(bins, function(x) {
    mean(x$dy)
  })

  # Generate matrix of mean x and y coordinate changes
  df <- data.frame(
    x = as.numeric(unlist(lapply(strsplit(names(dx_mean), "_"), function(v) {
      v[1]
    }))),
    y = as.numeric(unlist(lapply(strsplit(names(dx_mean), "_"), function(v) {
      v[2]
    }))),
    dx = unlist(dx_mean),
    dy = unlist(dy_mean)
  )

  # Calculate theta (angle between point of origin and the average xy
  # coordinates for all moments exiting that position)
  df$theta <- rep(NA, nrow(df))
  for (i in 1:nrow(df)) {
    x1 <- df[i, 1]
    y1 <- df[i, 2]
    x2 <- df[i, 1] + df[i, 3]
    y2 <- df[i, 2] + df[i, 4]

    df$theta[i] <- atan2(y2 - y1, x2 - x1) * (180 / pi)
  }

  # Calculate the euclidean distance between the point of origin and the mean
  # coordinates of subsequent time points
  df$dist <- rep(NA, nrow(df))
  for (i in 1:nrow(df)) {
    x1 <- df[i, 1]
    y1 <- df[i, 2]
    x2 <- df[i, 1] + df[i, 3]
    y2 <- df[i, 2] + df[i, 4]

    df$dist[i] <- euc_dist(c(x1, y1), c(x2, y2))
  }

  # Set up plot parameters
  par(mar = c(1, 1, 1, 1))

  # If desired, plot vector field with arrows colored by theta (angle)
  if (color_by_theta == TRUE) {
    # Round theta values
    theta <- round(df$theta)

    # Generate all possible angle values to assign colors to
    s <- seq(-180, 180, 1)

    # Generate color ramp palette
    cols <- c(
      colorRampPalette(c("midnightblue", "cyan4"))(90),
      colorRampPalette(c("cyan4", "lightgoldenrod1"))(90),
      colorRampPalette(c("lightgoldenrod1", "sienna2"))(90),
      colorRampPalette(c("sienna2", "orangered3"))(91)
    )

    # Match colors to angles
    names(cols) <- s
    cols <- cols[match(theta, names(cols))]

    # Initiate plot
    plot(
      df$x,
      df$y,
      xlim = c(min(df$x) - 2, max(df$x) + 2),
      ylim = c(min(df$y) - 2, max(df$y) + 2),
      type = "n",
      pch = 20,
      xlab = "",
      ylab = "",
      bty = "n",
      xaxt = "n",
      yaxt = "n"
    )
    # Add arrows colored by theta
    Arrows(
      df[, 1],
      df[, 2],
      df[, 1] + df[, 3] / 2,
      df[, 2] + df[, 4] / 2,
      arr.length = arrow_length,
      col = cols,
      arr.type = "triangle"
    )
  } else {
    # Initiate plain vanilla plot
    plot(
      df$x,
      df$y,
      xlim = c(min(df$x) - 2, max(df$x) + 2),
      ylim = c(min(df$y) - 2, max(df$y) + 2),
      type = "n",
      pch = 20,
      xlab = "",
      ylab = "",
      bty = "n",
      xaxt = "n",
      yaxt = "n"
    )
    # Add arrows corresponding to vector distances and angles
    Arrows(
      df[, 1],
      df[, 2],
      df[, 1] + df[, 3] / 2,
      df[, 2] + df[, 4] / 2,
      arr.length = arrow_length,
      col = arrow_color,
      arr.type = "triangle"
    )
  }
  # If desired, return the angle and distance values
  if (return == TRUE) {
    return(df)
  }
}

## 'plot_umap_features': Plots a behavior space layout colored by a desired
# feature/parameter
# Expects that the behavior space layout provided contains named columns
# corresponding to the desired features (called by 'feature_names')
# Feature must therefore be values that correspond to all time points
# contained with the behavior space (e.g. velocity values)
plot_umap_features <- function(layout,
                               windows,
                               bin_umap = FALSE,
                               n_bins = 32,
                               n_features = NULL,
                               feature_names = NULL,
                               colors = brewer.pal(11, "Spectral"),
                               plot_points = FALSE,
                               return = FALSE,
                               ...) {
  # Bin UMAP if desired
  if (bin_umap == TRUE) {
    layout <- bin_umap(layout,
      n_bins = n_bins
    )$layout
  }

  # Get vector of rows to split windows on (as a function of feature number)
  tosplit <- rep(1:n_features,
    each = (nrow(windows) / n_features)
  )

  # Split windows on features
  feat <- split(as.data.frame(windows), tosplit)

  # Get colors
  cols <- colorRampPalette(colors)(n_features)

  # Set up plotting aesthetics
  par(
    mfrow = c(1, n_features),
    bty = "n",
    xaxt = "n",
    yaxt = "n",
    mar = c(2, 2, 2, 2)
  )

  # Loop through features and plot
  for (i in 1:length(feat)) {
    # Calculate mean feature value per window
    m <- colMeans(feat[[i]])

    # Match to layout
    m <- m[1:nrow(layout)]

    # Add to layout
    if (is.null(feature_names) == FALSE) {
      layout[, feature_names[i]] <- m
    } else {
      layout <- cbind(layout, m)
    }

    if (plot_points == TRUE) {
      # Round mean feature value
      m <- round(m, 2)

      # Take absolute value
      m <- abs(m)

      # Get colors
      p <-
        colorRampPalette(c("grey90", cols[i]))(length(seq(0, max(m), 0.01)))
      names(p) <- seq(0, max(m), 0.01)
      p <- p[match(m, names(p))]

      # Plot
      plot(
        layout[, 1:2],
        pch = 20,
        col = p,
        ylab = "",
        xlab = "",
        ...
      )

      if (is.null(feature_names) == FALSE) {
        title(
          main = feature_names[i],
          cex.main = 1.5,
          font.main = 1
        )
      }
    } else {
      # Split on bin
      m <- split(m, layout$xy_new)

      # Get mean per bin
      m <- lapply(m, function(x) {
        mean(x, na.rm = TRUE)
      })

      # Unlist
      m <- unlist(m)

      # Round mean feature value
      m <- round(m, 2)

      # Take absolute value
      m <- abs(m)

      # Get colors
      p <-
        colorRampPalette(c("grey90", cols[i]))(length(seq(0, max(m), 0.01)))
      names(p) <- seq(0, max(m), 0.01)
      p <- p[match(m, names(p))]

      # Plot
      plot(
        unlist(lapply(strsplit(names(
          m
        ), "_"), function(v) {
          v[1]
        })),
        unlist(lapply(strsplit(names(
          m
        ), "_"), function(v) {
          v[2]
        })),
        pch = 20,
        col = p,
        ylab = "",
        xlab = "",
        ...
      )

      if (is.null(feature_names) == FALSE) {
        title(
          main = feature_names[i],
          cex.main = 1.5,
          font.main = 1
        )
      }
    }
  }
  if (return == TRUE) {
    return(layout)
  }
}

##' plot_umap_pdf': Calculates a probability density function of points in a
# behavior space and plots
plot_umap_pdf <- function(layout,
                          h = 1,
                          n = 100,
                          colors = matlab.like(100),
                          return = FALSE) {
  # Calcuate probability density function using 'kde2d' from MASS, height ('h')
  # and resolution ('n') are determined by calls to the main function
  pdf <- kde2d(layout$x,
    layout$y,
    h = h,
    n = n
  )

  # Initiate plot
  par(
    mar = c(1, 1, 1, 1),
    bty = "n",
    xaxt = "n",
    yaxt = "n"
  )

  # Plot with image
  image(pdf$z,
    xlab = "",
    ylab = "",
    col = colors
  )

  # Return if desired
  if (return == TRUE) {
    return(pdf)
  }
}
