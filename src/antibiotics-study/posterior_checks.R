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
      mx %>%
      rename(value = truth),
      iteration = NA,
      type = "truth"
    ),
    cbind(
      m_sim %>%
      filter(iteration %in% round(seq(1, iter_max, length.out = n_vis))) %>%
      rename(value = sim_value),
      type = "simulated"
    )
  )

  ## make histograms
  ggplot(hist_data) +
    geom_histogram(aes(x = asinh(value), fill = type), bins = 100) +
    facet_grid(iteration ~ .) +
    scale_fill_manual(values = c("#377eb8", "#4daf4a")) +
    theme(
      panel.border = element_rect(fill = "transparent", size = 0.5)
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
compare_quantiles <- function(mx, m_sim, q_probs = NULL) {
  if (is.null(q_probs)) {
    q_probs <- seq(0, 1, 0.01)
  }

  quantiles_comp <- m_sim %>%
    group_by(iteration) %>%
    do(
      data.frame(
        type = "sim",
        q_ix = q_probs,
        q = quantile(asinh(.$sim_value), q_probs)
      )
    )

  ggplot(quantiles_comp) +
    geom_step(
      aes(x = q, y = q_ix, group = iteration),
      alpha = 0.1, position = position_jitter(h = 0.005),
    ) +
    geom_step(
      data = data.frame(
        q_ix = q_probs,
        q = quantile(asinh(mx$truth), q_probs)
      ),
      aes(x = q, y = q_ix),
      col = "#79B5B7",
      size = 0.5
    ) +
    labs(
      "x" = "x",
      "y" = "Pr(asinh(count) < x)"
    )
}

compare_margins <- function(mx, m_sim, group_col) {
  group_totals <- mx %>%
    group_by_(group_col) %>%
    summarise(group_total = sum(asinh(truth)))
  group_totals$rank <- rank(group_totals$group_total)

  sim_group_totals <- m_sim %>%
    group_by_("iteration", group_col) %>%
    summarise(sim_total = sum(asinh(sim_value))) %>%
    left_join(group_totals)

  ggplot() +
    geom_boxplot(
      data = sim_group_totals,
      aes(y = as.factor(rank), x = sim_total),
      alpha = 0.1, size = 0.1
    ) +
    geom_step(
      data = sim_group_totals %>% filter(iteration == 1),
      aes(y = rank, x = group_total),
      col = "#79B5B7"
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

  data.frame(aligned_scores, supp_cols)
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

  data.frame(aligned_loadings, supp_cols)
}

evals_summary <- function(data_list, supp_cols) {
  evals <- princomp(scale(data_list$x_sim, scale = FALSE))$sdev
  data.frame(evals, supp_cols)
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

  rbindlist(stat_list)
}

summary_contours <- function(summary_data, plot_opts) {
  ggcontours(
    summary_data,
    plot_opts
  ) +
    geom_text(
      data = summary_data %>%
        filter(type != "true") %>%
        group_by(row_ix) %>%
        summarise(X1 = mean(X1), X2 = mean(X2)),
      aes(x = X1, y = X2, label = row_ix),
      size = 4
    ) +
    geom_text(
      data = summary_data %>% filter(type == "true"),
      aes(x = X1, y = X2, label = row_ix),
      col = "#79B5B7", size = 4
    )
}

counts_data_checker <- function(x, x_sim, output_dir = ".") {
  ## ---- posterior-data ----
  m_sim <- x_sim %>%
    melt(
      varnames = c("iteration", "sample", "rsv"),
      value.name = "sim_value"
    ) %>%
    as.data.table()

  mx <- x %>%
    melt(
      varnames = c("sample", "rsv"),
      value.name = "truth"
    ) %>%
    as.data.table()

  all_plots <- list()
  all_plots[["hists"]] <- compare_histograms(mx, m_sim)
  all_plots[["quantiles"]] <- compare_quantiles(mx, m_sim)
  all_plots[["margins"]] <- compare_margins(mx, m_sim, "rsv")

  mx_samples <- mx
  mx_samples$sample_id  <- sample_names(abt)[mx_samples$sample]
  mx_samples <- mx_samples %>%
    left_join(
      data.frame(
        sample_id = sample_names(abt),
        sample_data(abt)
      )
    ) %>%
    filter(rsv %in% sample(seq_len(ntaxa(abt)), 12)) %>%
    left_join(m_sim)

  all_plots[["ts"]] <- ggplot() +
    geom_point(
      data = mx_samples,
      aes(x = time, y = asinh(sim_value), group = interaction(iteration, rsv)),
      alpha = 0.01, size = 0.1
    ) +
    geom_line(
      data = mx_samples %>% filter(iteration == 1),
      aes(x = time, y = asinh(truth), group = rsv),
      size = 0.5, col = "#79B5B7"
    ) +
    facet_wrap(~rsv, scales = "free", ncol = 4)

  scores_data <- sample_summary_fun(
    asinh(t(x)),
    aperm(asinh(x_sim), c(1, 3, 2)),
    scores_summary,
    list("K" = 2)
  )

  loadings_data <- sample_summary_fun(
    asinh(t(x)),
    aperm(asinh(x_sim), c(1, 3, 2)),
    loadings_summary,
    list("K" = 2)
  )

  evals_data <- sample_summary_fun(
    asinh(t(x)),
    aperm(asinh(x_sim), c(1, 3, 2)),
    evals_summary,
    list()
  )

  plot_opts <- list(
    "x" = "X1",
    "y" = "X2",
    "group" = "row_ix",
    "h" = 1.5
  )

  all_plots[["scores"]] <- summary_contours(scores_data, plot_opts) +
    coord_fixed(0.5)

  plot_opts$h <- 0.01
  all_plots[["loadings"]] <- summary_contours(loadings_data, plot_opts) +
    coord_fixed(0.5)

  all_plots[["evals"]] <- ggplot() +
    geom_boxplot(
      data = evals_data %>%
        filter(type == "sim"),
      aes(x = as.factor(row_ix), y = evals),
      outlier.size = 0.1, size = 0.1
    ) +
    geom_point(
      data = evals_data %>%
        filter(type == "true"),
      aes(x = as.factor(row_ix), y = evals),
      col = "#79B5B7", size = 0.9
    ) +
    ylim(0, 11) +
    scale_y_log10() +
    theme(
      axis.text.x = element_blank()
    ) +
    labs(
      "x" = "Index",
      "y" = "Eigenvalue"
    )

  for (i in seq_along(all_plots)) {
    dir.create(output_dir, recursive = TRUE)
    ggsave(
      file = sprintf("%s/figure-%s.png", output_dir, names(all_plots)[i]),
      all_plots[[i]]
    )
  }
  all_plots
}
