plot_combined_psmc <- function(
    x1,  # BBCC02_psmc (will be blue)
    x2,  # T_psmc (will be red)
    type = "step", 
    xlim = NULL, 
    ylim = NULL,
    xlab = "Generations from present", 
    ylab = "Ne",
    show_present = TRUE, 
    mutation_rate = 1e-8, 
    generation_time = 1,
    scaled = FALSE, 
    bin_size = 100,
    show_bootstrap = TRUE,
    bootstrap_alpha = 0.2,
    line_size = 1.5,
    ...
) {
  # Process first PSMC (BBCC02 - blue)
  RS1 <- x1$RS[x1$RS[, "iter"] == x1$niters, ]
  theta0_1 <- x1$parameters[nrow(x1$parameters), "theta0"]
  
  if (scaled) {
    xx1 <- RS1[, "t_k"] / (theta0_1 * bin_size)
    yy1 <- theta0_1 * RS1[, "lambda_k"] / bin_size
  } else {
    denom <- 4 * mutation_rate * generation_time * bin_size
    N0_1 <- theta0_1 / denom
    xx1 <- 2 * N0_1 * RS1[, "t_k"]
    yy1 <- N0_1 * RS1[, "lambda_k"]
  }
  
  df_main1 <- data.frame(time = xx1, Ne = yy1, group = "BBCC02")
  
  # Process second PSMC (T - red)
  RS2 <- x2$RS[x2$RS[, "iter"] == x2$niters, ]
  theta0_2 <- x2$parameters[nrow(x2$parameters), "theta0"]
  
  if (scaled) {
    xx2 <- RS2[, "t_k"] / (theta0_2 * bin_size)
    yy2 <- theta0_2 * RS2[, "lambda_k"] / bin_size
  } else {
    denom <- 4 * mutation_rate * generation_time * bin_size
    N0_2 <- theta0_2 / denom
    xx2 <- 2 * N0_2 * RS2[, "t_k"]
    yy2 <- N0_2 * RS2[, "lambda_k"]
  }
  
  df_main2 <- data.frame(time = xx2, Ne = yy2, group = "T")
  
  # Combine main trajectories
  df_main <- rbind(df_main1, df_main2)
  
  # Calculate minimum effective population size
  min_Ne <- min(df_main$Ne)
  cat("Minimum Ne in data set:", min_Ne, "\n")
  
  # Create base plot
  p <- ggplot(df_main, aes(x = time, y = Ne, color = group)) +
    geom_step(aes(group = group), 
              linewidth = line_size,
              linetype = "solid") +
    geom_hline(yintercept = min_Ne, lty = 2) +
    scale_color_manual(values = c("BBCC02" = "blue", "T" = "red")) +
    labs(x = xlab, y = ylab, color = "Group") +
    theme_classic() +
    theme(legend.position = "top")
  
  # Add bootstrap trajectories if available
  if (show_bootstrap) {
    # Process bootstrap for first PSMC
    if (!is.null(x1$bootstrap)) {
      boot1 <- x1$bootstrap
      THETA0_1 <- rep(boot1$theta0, each = x1$n)
      
      if (scaled) {
        Tk_boot1 <- boot1$tk / (THETA0_1 * bin_size)
        Nk_boot1 <- THETA0_1 * boot1$lk / bin_size
      } else {
        denom <- 4 * mutation_rate * generation_time * bin_size
        N0_boot1 <- boot1$theta0 / denom
        Tk_boot1 <- 2 * N0_boot1 * boot1$tk
        Nk_boot1 <- N0_boot1 * boot1$lk
      }
      
      df_boot1 <- data.frame(
        time = as.vector(Tk_boot1),
        Ne = as.vector(Nk_boot1),
        rep = rep(1:ncol(Tk_boot1), each = nrow(Tk_boot1)),
        group = "BBCC02"
      )
      
      p <- p + 
        geom_step(data = df_boot1, 
                  aes(group = interaction(rep, group)),
                  color = "lightblue",
                  alpha = bootstrap_alpha,
                  linewidth = line_size * 0.7)
    }
    
    # Process bootstrap for second PSMC
    if (!is.null(x2$bootstrap)) {
      boot2 <- x2$bootstrap
      THETA0_2 <- rep(boot2$theta0, each = x2$n)
      
      if (scaled) {
        Tk_boot2 <- boot2$tk / (THETA0_2 * bin_size)
        Nk_boot2 <- THETA0_2 * boot2$lk / bin_size
      } else {
        denom <- 4 * mutation_rate * generation_time * bin_size
        N0_boot2 <- boot2$theta0 / denom
        Tk_boot2 <- 2 * N0_boot2 * boot2$tk
        Nk_boot2 <- N0_boot2 * boot2$lk
      }
      
      df_boot2 <- data.frame(
        time = as.vector(Tk_boot2),
        Ne = as.vector(Nk_boot2),
        rep = rep(1:ncol(Tk_boot2), each = nrow(Tk_boot2)),
        group = "T"
      )
      
      p <- p + 
        geom_step(data = df_boot2, 
                  aes(group = interaction(rep, group)),
                  color = "pink",
                  alpha = bootstrap_alpha,
                  linewidth = line_size * 0.7)
    }
  }
  
  # Add present time annotation
  if (show_present) {
    p <- p + 
      annotate("text", x = 0, y = max(df_main$Ne) * 0.05, 
               label = "", fontface = "italic", hjust = 0)
  }
  
  # Set axis limits
  if (!is.null(xlim)) p <- p + xlim(xlim)
  if (!is.null(ylim)) {
    p <- p + ylim(ylim)
  } else {
    y_max <- max(df_main$Ne)
    if (show_bootstrap) {
      if (!is.null(x1$bootstrap)) y_max <- max(y_max, df_boot1$Ne)
      if (!is.null(x2$bootstrap)) y_max <- max(y_max, df_boot2$Ne)
    }
    p <- p + expand_limits(y = c(0, y_max * 1.05))
  }
  
  return(p)
}