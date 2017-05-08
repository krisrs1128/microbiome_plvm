#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
## Functions to facilitate posterior predictive checks in the antibiotics study.
##
## author: kriss1@stanford.edu

#' Make histograms comparing simulated and true values
#'
#' @param mx [data.frame] A data.frame of the melted original counts data
#' @param m_sim [data.frame] A data.frame of the melted simulation array.
#'   Includes an extra column for the iteration number.
#' @param n_vis [int] The number of simulated samples to show
#' @return p [ggplot] The histogram plot object
compare_histograms <- function(mx, m_sim, n_vis = 4) {

  ## bind simulated and true data, in order to plot
  iter_max <- max(m_sim$iteration)
  hist_data <- rbind(
    cbind(
      iteration = NA,
      mx %>%
      rename(value = truth) %>%
      filter(method == "lda") %>%
      mutate(method = "truth")
    ),
    m_sim %>%
      filter(iteration %in% round(seq(1, iter_max, length.out = n_vis))) %>%
      rename(value = sim_value)
  )

  ## make histograms
  ggplot(hist_data) +
    geom_histogram(aes(x = asinh(value), fill = method), position = "dodge", bins = 100) +
    facet_grid(iteration ~ .) +
    scale_fill_manual(values = c("#86B8B1", "#000000", "#b186b8")) +
  theme(
    panel.border = element_rect(fill = "transparent", size = 0.5),
    legend.position = "bottom"
  )
}

#' Plot quantiles of true vs. simulated data
#'
#' @param mx [data.frame] A data.frame of the melted original counts data
#' @param m_sim [data.frame] A data.frame of the melted simulation array.
#'   Includes an extra column for the iteration number.
#' @param n_vis [int] The number of simulated samples to show
#' @return p [ggplot] A plot comparing the true quantiles in mx to those in the
#'   simulations m_sim.
compare_quantiles <- function(mx, q_sim) {
  q_sim$q_ix <- 0.01 * as.numeric(gsub("\\%", "", q_sim$q_ix))
  q_true <- data_frame(
    q_ix = seq(0, 1, 0.01),
    q_true = quantile(asinh(mx$truth), seq(0, 1, 0.01))
  )

  q_joined <- q_sim %>%
    full_join(q_true)

  ggplot(q_joined) +
    geom_abline(slope = 1, alpha = 0.6, size = 0.4, col = "#000000") +
    geom_point(
      aes(x = q_true, y = q, col = method),
      alpha = 0.05, size = 0.1,
      position = position_jitter(w = 0.2, h = 0.2)
    ) +
    coord_fixed() +
    scale_color_manual(values = c("#86B8B1", "#b186b8")) +
    guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2))) +
    labs(
      "x" = "Observed Quantiles",
      "y" = "Posterior Predictive Quantiles"
    )
}

compare_margins <- function(mx, m_sim) {
  group_totals <- mx %>%
    filter(method == "lda") %>%
    group_by(rsv) %>%
    summarise(group_total = sum(asinh(truth)))
  group_totals$rank <- rank(group_totals$group_total)

  sim_group_totals <- m_sim %>%
    group_by(iteration, method, rsv) %>%
    summarise(sim_total = sum(asinh(sim_value))) %>%
    left_join(group_totals)

  ggplot() +
    geom_boxplot(
      data = sim_group_totals,
      aes(y = as.factor(rank), x = sim_total, fill = method, col = method),
      alpha = 0.1, size = 0.1
    ) +
    geom_step(
      data = group_totals,
      aes(y = rank, x = group_total),
      col = "#9C964A"
    ) +
    labs(
      "x" = "x",
      "y" = "Prob(transformed microbe sum < x)"
    ) +
    theme(
      axis.text.y = element_blank()
    )
}

scores_summary <- function(data_list, supp_cols) {
  library("vegan")
  get_scores <- function(x) {
    princomp(scale(x, scale = FALSE))$scores[, seq_len(data_list$K)]
  }
  true_scores <- get_scores(data_list$x)
  scores <- get_scores(data_list$x_sim)
  aligned_scores <- procrustes(true_scores, scores)$Yrot
  dimnames(aligned_scores) <- NULL

  cbind(
    as_data_frame(aligned_scores),
    as_data_frame(supp_cols)
  )
}

loadings_summary <- function(data_list, supp_cols) {
  library("vegan")
  get_loadings <- function(x) {
    princomp(scale(x, scale = FALSE))$loadings[, seq_len(data_list$K)]
  }

  true_loadings <- get_loadings(data_list$x)
  loadings <- get_loadings(data_list$x_sim)
  aligned_loadings <- procrustes(true_loadings, loadings)$Yrot
  dimnames(aligned_loadings) <- NULL

  cbind(
    as_data_frame(aligned_loadings),
    as_data_frame(supp_cols)
  )
}

evals_summary <- function(data_list, supp_cols) {
  evals <- princomp(scale(data_list$x_sim, scale = FALSE))$sdev

  cbind(
    as_data_frame(evals),
    as_data_frame(supp_cols)
  )
}

sample_summary_fun <- function(x, x_sim, summary_fun, data_opts) {
  stat_list <- vector(
    length = nrow(x_sim),
    mode = "list"
  )

  stat_list[[1]] <- summary_fun(
    c(list("x" = x, "x_sim" = x), data_opts),
    list("iteration" = NA, "type" = "true")
  )
  stat_list[[1]]$row_ix <- seq_len(nrow(stat_list[[1]]))

  for (i in seq_along(stat_list)[-1]) {
    if (i %% 50 == 0) {
      cat(sprintf("processing %s\n", i))
    }

    stat_list[[i]] <- summary_fun(
      c(list("x" = x, "x_sim" = x_sim[i - 1,, ]), data_opts),
      list("iteration" = i - 1, "type" = "sim")
    )
    stat_list[[i]]$row_ix <- seq_len(nrow(stat_list[[i]]))
  }

  do.call(rbind, stat_list)
}

summary_contours <- function(summary_data, plot_opts, text_size = 3) {
  library("ggscaffold")
  ggcontours(
    summary_data %>% filter(type != "true"),
    plot_opts
  ) +
  geom_text(
    data = summary_data %>% filter(type == "true"),
    aes(x = V1, y = V2, label = row_ix),
    col = "#000000", size = text_size, alpha = 0.9
  ) +
  geom_text(
    data = summary_data %>%
      filter(type != "true") %>%
      group_by(row_ix, method) %>%
      summarise(V1 = mean(V1), V2 = mean(V2)),
    aes(x = V1, y = V2, label = row_ix, col = method),
    size = text_size
  ) +
  scale_color_manual(values = c("#3d6862", "#714678"), guide = FALSE) +
  guides(fill = guide_legend(override.aes = list(alpha = 1), ncol = 1)) +
  facet_grid(method ~ .)
}

posterior_checks_input <- function(x, x_sim, file_basename = NULL) {
  q_sim <- apply(asinh(x_sim), 1, quantile, seq(0, 1, 0.01))  %>%
    melt(
      varnames = c("q_ix", "iteration"),
      value.name = "q"
    ) %>%
    as_data_frame()

  mx <- x %>%
    melt(
      varnames = c("sample", "rsv"),
      value.name = "truth"
    ) %>%
    as_data_frame()

  mx_samples <- mx
  mx_samples$sample_id  <- sample_names(abt)[mx_samples$sample]

  ## show taxa chosen randomly among those present in >= 45% of samples
  keep_taxa <- c(160, 343, 891, 1036, 1045, 1086, 1116, 1131, 1417, 1463, 1659, 1895)
  m_sim <- x_sim[,, keep_taxa] %>%
    melt(
      varnames = c("iteration", "sample", "rsv"),
      value.name = "sim_value"
    ) %>%
    as_data_frame()
  m_sim$rsv <- keep_taxa[m_sim$rsv]

  mx_samples <- mx_samples %>%
    filter(rsv %in% keep_taxa) %>%
    left_join(
      cbind(
        sample_id = sample_names(abt),
        as_data_frame(sample_data(abt))
      )
    ) %>%
    left_join(m_sim)

  ## prefilter taxa used in PCA, so that the loadings plot is readable
  pca_taxa <- order(apply(asinh(x), 2, var), decreasing = TRUE)[1:1000]
  scores_data <- sample_summary_fun(
    asinh(t(x[, pca_taxa])),
    aperm(asinh(x_sim[,, pca_taxa]), c(1, 3, 2)),
    scores_summary,
    list("K" = 2)
  )

  loadings_data <- sample_summary_fun(
    asinh(t(x[, pca_taxa])),
    aperm(asinh(x_sim[,, pca_taxa]), c(1, 3, 2)),
    loadings_summary,
    list("K" = 2)
  )

  evals_data <- sample_summary_fun(
    asinh(t(x[, pca_taxa])),
    aperm(asinh(x_sim[,, pca_taxa]), c(1, 3, 2)),
    evals_summary,
    list()
  )

  input_data <- list(
    "q_sim" = q_sim,
    "mx" = mx,
    "mx_samples" = mx_samples,
    "scores_data" = scores_data,
    "loadings_data" = loadings_data,
    "evals_data" = evals_data
  )

  if (!is.null(file_basename)) {
    for (i in seq_along(input_data)) {
      write_feather(
        input_data[[i]],
        sprintf(
          "%s-%s.feather",
          file_basename,
          names(input_data[i])
        )
      )
    }
  }
  input_data
}

posterior_checks_plots <- function(input_data) {
  library("ggplot2")
  all_plots <- list()
  all_plots[["quantiles"]] <- compare_quantiles(input_data$mx, input_data$q_sim)
  input_data$m_sim <- NULL

  all_plots[["ts"]] <- ggplot() +
    geom_point(
      data = input_data$mx_samples,
      aes(x = time, y = asinh(sim_value), col = method),
      alpha = 0.01, size = 0.1
    ) +
    geom_line(
      data = input_data$mx_samples %>% filter(iteration == 1),
      aes(x = time, y = asinh(truth), group = rsv),
      size = 0.4, col = "#000000"
    ) +
    scale_color_manual(values = c("#86B8B1", "#b186b8")) +
    scale_y_continuous(breaks = scales::pretty_breaks(3)) +
    labs(x = "time", y = "asinh(abundance)") +
    guides(colour = guide_legend(nrow = 2, override.aes = list(alpha = 1, size = 2))) +
    facet_wrap(~rsv, scales = "free", ncol = 3)

  plot_opts <- list(
    "x" = "V1",
    "y" = "V2",
    "group" = "row_ix",
    "fill" = "method",
    "fill_type" = "category",
    "fill_colors" = c("#86B8B1","#b186b8"),
    "h" = 2.5,
    "theme_opts" = list("text_size" = 10, "subtitle_size" = 11, "key_width" = 1,
                        "border_size" = 1)
  )

  all_plots[["scores"]] <- summary_contours(input_data$scores_data, plot_opts, 2.3) +
    labs(x = "Axis 1", y = "Axis 2") +
    coord_fixed(0.8)

  plot_opts$h <- 0.01
  plot_opts$theme_opts$legend_position <- "none"
  all_plots[["loadings"]] <- summary_contours(input_data$loadings_data, plot_opts, 2.3) +
    scale_x_continuous(breaks = c(0.05, 0, 0.05)) +
    labs(x = "Axis 1", y = "Axis 2") +
    coord_fixed(0.4)

  all_plots[["evals"]] <- ggplot() +
   geom_point(
      data = input_data$evals_data %>%
        filter(type == "sim"),
      aes(x = as.factor(row_ix), y = log(value, 10), col = method),
      alpha = 0.05, size = 0.05, position = position_jitter(h = 0, w = 0.25)
    ) +
    geom_point(
      data = input_data$evals_data %>%
        filter(type == "true"),
      aes(x = as.factor(row_ix), y = log(value, 10)),
      col = "#000000", size = 0.9
    ) +
    ylim(log(0.25, 10), log(11, 10)) +
    scale_color_manual(values = c("#86B8B1", "#b186b8")) +
    scale_y_continuous(breaks = pretty_breaks(n = 2)) +
    guides(colour = guide_legend(override.aes = list(alpha = 1, size = 1), ncol = 1)) +
    theme(axis.text.x = element_blank()) +
    labs(
      "x" = expression(i),
      "y" = expression(log[10](lambda[i]))
    )

  all_plots
}
