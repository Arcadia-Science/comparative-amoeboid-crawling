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

############################
#####Plotting functions#####
############################
# Function to plot parameters by behavior space position
plot_parameter <-
  function(parameter, layout, n_bins, return = FALSE, ...) {
    layout <- bin_umap(layout, n_bins = n_bins)$layout
    
    a <- parameter
    a <- split(a[1:nrow(layout)], layout$xy_new)
    a <- unlist(lapply(a, function(x)
      mean(x, na.rm = TRUE)))
    
    p <- expand.grid(seq(1, n_bins + 1, 1), seq(1, n_bins + 1, 1))
    p <- paste(p[, 1], p[, 2], sep = "_")
    a <- a[match(p, names(a))]
    names(a) <- p
    
    image(
      matrix(a, nrow = n_bins + 1),
      xaxt = "n",
      yaxt = "n",
      xlab = "",
      ylab = "",
      bty = "n",
      ...
    )
    
    if (return == TRUE) {
      return(a)
    }
  }

# Function to plot results of iterative tests
plot_results <- function(res_list,
                         ylim = c(4, 9),
                         ylab = NULL,
                         xlab = NULL,
                         plot_as_lines = FALSE,
                         col1 = "grey50",
                         col2 = "grey80",
                         ...) {
  means <- unlist(lapply(res_list, function(x)
    mean(x)))
  error <- lapply(res_list, function(x)
    boxplot.stats(x)$stats)
  
  if (plot_as_lines == TRUE) {
    means <- unlist(lapply(res_list, function(x)
      median(x)))
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
    polygon(c(seq(1, length(means), 1),
              rev(seq(
                1, length(means), 1
              ))),
            c(lapply(error, function(x)
              x[5]),
              rev(lapply(error, function(x)
                x[1]))),
            col = col2,
            border = FALSE)
    lines(means,
          type = "l",
          col = col1,
          lwd = 3)
    axis(
      1,
      at = seq(1, length(means), 1),
      labels = names(means),
      cex.lab = 1.5,
      cex.axis = 1.5,
      las = 2
    )
  } else {
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
    axis(
      1,
      at = seq(1, length(means), 1),
      labels = names(means),
      cex.lab = 1.5,
      cex.axis = 1.5,
      las = 2
    )
    
    for (i in 1:length(res_list)) {
      points(jitter(rep(i, length(res_list[[i]])), 0.25),
             res_list[[i]],
             pch = 20,
             col = alpha(col2, 0.5))
    }
    
    points(
      means,
      pch = 21,
      cex = 1.5,
      bg = col1,
      col = col1
    )
  }
  
  l <- list(means, error)
  names(l) <- c("means", "error")
  return(l)
}

# Function to plot results of iterative tests as normalized variance
plot_variance <- function(res_list,
                          ylim = c(0, 1),
                          ylab = NULL,
                          xlab = NULL,
                          return = FALSE,
                          ...) {
  variance <- unlist(lapply(res_list, function(x)
    sd(x) / mean(x)))
  
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
  axis(
    1,
    at = seq(1, length(variance), 1),
    labels = names(variance),
    cex.lab = 1.5,
    cex.axis = 1.5,
    las = 2
  )
  
  points(
    variance,
    pch = 21,
    cex = 1.5,
    bg = "grey40",
    col = "grey40"
  )
  
  if (return == TRUE) {
    return(variance)
  }
}

# Function to plot recurrence results
plot_recurrence <- function(recurrences,
                            mar = c(0.5, 0.5, 0.5, 0.5)) {
  # Analyze distribution of recurrences
  par(mfrow = c(length(recurrences), 1))
  par(mar = mar)
  
  for (i in 1:length(recurrences)) {
    image(
      do.call(
        cbind,
        lapply(recurrences[[i]], function(y)
          y$proportion_recurrent_in_bins)
      ),
      col = colorRampPalette(hcl.colors(12, "YlOrRd", rev = TRUE))(100),
      xaxt = "n",
      yaxt = "n"
    )
    title(
      main = paste(names(recurrences)[i], "frames"),
      cex.main = 1.5,
      font.main = 1
    )
  }
  
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

# Function to plot as vector field
plot_vector_field <- function(layout,
                              bin_umap = FALSE,
                              n_bins = 32,
                              color_by_theta = FALSE,
                              arrow_color = "grey50",
                              arrow_length = 0.05,
                              return = FALSE) {
  if (bin_umap == TRUE) {
    layout <- bin_umap(layout, n_bins = n_bins)$layout
  }
  
  layout$dx <- c(0, diff(layout$x))
  layout$dy <- c(0, diff(layout$y))
  
  bins <- split(layout, layout$xy)
  dx_mean <- lapply(bins, function(x)
    mean(x$dx))
  dy_mean <- lapply(bins, function(x)
    mean(x$dy))
  
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
  df$theta <- rep(NA, nrow(df))
  for (i in 1:nrow(df)) {
    x1 <- df[i, 1]
    y1 <- df[i, 2]
    x2 <- df[i, 1] + df[i, 3]
    y2 <- df[i, 2] + df[i, 4]
    
    df$theta[i] <- atan2(y2 - y1, x2 - x1) * (180 / pi)
  }
  
  df$dist <- rep(NA, nrow(df))
  for (i in 1:nrow(df)) {
    x1 <- df[i, 1]
    y1 <- df[i, 2]
    x2 <- df[i, 1] + df[i, 3]
    y2 <- df[i, 2] + df[i, 4]
    
    df$dist[i] <- euc.dist(c(x1, y1), c(x2, y2))
  }
  
  par(mar = c(1, 1, 1, 1))
  
  if (color_by_theta == TRUE) {
    x <- round(df$dy, 2)
    cols <-
      colorRampPalette(c("cyan4", "grey90", "orangered3"))(length(seq(min(x), max(x), 0.01)))
    names(cols) <- round(seq(min(x), max(x), 0.01), 2)
    cols <- cols[match(as.numeric(x),
                       as.numeric(names(cols)))]
    
    theta <- round(df$theta)
    s <- seq(-180, 180, 1)
    cols <- c(
      colorRampPalette(c("midnightblue", "cyan4"))(90),
      colorRampPalette(c("cyan4", "lightgoldenrod1"))(90),
      colorRampPalette(c("lightgoldenrod1", "sienna2"))(90),
      colorRampPalette(c("sienna2", "orangered3"))(91)
    )
    names(cols) <- s
    cols <- cols[match(theta, names(cols))]
    
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
  
  if (return == TRUE) {
    return(df)
  }
}

# Function to plot with features colored
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
                       n_bins = n_bins)$layout
  }
  
  # Get vector of rows to split windows on (as a function of feature number)
  tosplit <- rep(1:n_features,
                 each = (nrow(windows) / n_features))
  
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
        title(main = feature_names[i],
              cex.main = 1.5,
              font.main = 1)
      }
    } else {
      # Split on bin
      m <- split(m, layout$xy_new)
      
      # Get mean per bin
      m <- lapply(m, function(x)
        mean(x, na.rm = TRUE))
      
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
        title(main = feature_names[i],
              cex.main = 1.5,
              font.main = 1)
      }
    }
  }
  if (return == TRUE) {
    return(layout)
  }
}

# Function plot as a probability density function
plot_umap_pdf <- function(layout,
                          h = 1,
                          n = 100,
                          colors = matlab.like(100),
                          return = FALSE) {
  # Get pdf
  pdf <- kde2d(layout$x,
               layout$y,
               h = h,
               n = n)
  
  # Plot
  par(
    mar = c(1, 1, 1, 1),
    bty = "n",
    xaxt = "n",
    yaxt = "n"
  )
  image(pdf$z,
        xlab = "",
        ylab = "",
        col = colors)
  
  # Return if desired
  if (return == TRUE) {
    return(pdf)
  }
}
