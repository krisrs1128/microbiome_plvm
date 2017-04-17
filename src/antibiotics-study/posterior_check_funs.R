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
    scale_fill_manual(values = wes_palette("Moonrise3", 3)) +
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
    group_by(iteration, method) %>%
    do(
      data_frame(
        type = "sim",
        q_ix = q_probs,
        q = quantile(asinh(.$sim_value), q_probs)
      )
    )

  ggplot(quantiles_comp) +
    geom_step(
      aes(x = q, y = q_ix, col = method, group = iteration),
      alpha = 0.1, position = position_jitter(h = 0.005),
    ) +
    geom_step(
      data = data_frame(
        q_ix = q_probs,
        q = quantile(asinh(mx$truth), q_probs)
      ),
      aes(x = q, y = q_ix, col = "#85D4E3"),
      size = 0.5
    ) +
    scale_color_manual(values = wes_palette("Moonrise3", 3)) +
    labs(
      "x" = "x",
      "y" = "Pr(asinh(count) < x)"
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

summary_contours <- function(summary_data, plot_opts) {
  library("ggscaffold")
  summary_data$method[summary_data$type == "true"] <- "truth"
  ggcontours(
    summary_data %>% filter(type != "true"),
    plot_opts
  ) +
    geom_text(
      data = summary_data %>%
        filter(type != "true") %>%
        group_by(row_ix, method) %>%
        summarise(V1 = mean(V1), V2 = mean(V2)),
      aes(x = V1, y = V2, label = row_ix, col = method),
      size = 4
    ) +
    geom_text(
      data = summary_data %>% filter(type == "true"),
      aes(x = V1, y = V2, label = row_ix),
      col = "#85D4E3", size = 4
    ) +
    scale_color_manual(values = wes_palette("Moonrise3", 3), guide = FALSE) +
    facet_grid(. ~ method)
}

posterior_checks_input <- function(x, x_sim, file_basename = NULL) {
  m_sim <- x_sim %>%
    melt(
      varnames = c("iteration", "sample", "rsv"),
      value.name = "sim_value"
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
  mx_samples <- mx_samples %>%
    left_join(
      cbind(
        sample_id = sample_names(abt),
        as_data_frame(sample_data(abt))
      )
    ) %>%
    filter(rsv %in% sample(seq_len(ntaxa(abt)), 12)) %>%
    left_join(m_sim)

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

  input_data <- list(
    "m_sim" = m_sim,
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

posterior_checks_plots <- function(input_data, output_dir = ".", ...) {
  library("wesanderson")
  all_plots <- list()
  all_plots[["hists"]] <- compare_histograms(input_data$mx, input_data$m_sim)
  all_plots[["quantiles"]] <- compare_quantiles(input_data$mx, input_data$m_sim)
  all_plots[["margins"]] <- compare_margins(input_data$mx, input_data$m_sim)

  all_plots[["ts"]] <- ggplot() +
    geom_point(
      data = input_data$mx_samples,
      aes(x = time, y = asinh(sim_value), col = method),
      alpha = 0.01, size = 0.1
    ) +
    geom_line(
      data = input_data$mx_samples %>% filter(iteration == 1),
      aes(x = time, y = asinh(truth), group = rsv),
      size = 0.5, col = "#85D4E3"
    ) +
    scale_color_manual(values = wes_palette("Moonrise3", 3)) +
    labs(x = "time", y = "asinh(abundance)") +
    facet_wrap(~rsv, scales = "free", ncol = 6)

  plot_opts <- list(
    "x" = "V1",
    "y" = "V2",
    "group" = "row_ix",
    "fill" = "method",
    "fill_type" = "category",
    "fill_cols" = wes_palette("Moonrise3", 3),
    "h" = 1.5
  )

  all_plots[["scores"]] <- summary_contours(input_data$scores_data, plot_opts) +
    coord_fixed(0.5)

  plot_opts$h <- 0.01
  all_plots[["loadings"]] <- summary_contours(input_data$loadings_data, plot_opts) +
    coord_fixed(0.5)

  all_plots[["evals"]] <- ggplot() +
    geom_point(
      data = input_data$evals_data %>%
        filter(type == "true"),
      aes(x = as.factor(row_ix), y = value),
      col = "#85D4E3", size = 0.9
    ) +
    geom_boxplot(
      data = input_data$evals_data %>%
        filter(type == "sim"),
      aes(x = as.factor(row_ix), y = value, col = method, fill = method),
      outlier.size = 0.1, size = 0.1
    ) +
    ylim(0, 11) +
    scale_y_log10() +
    scale_color_manual(values = wes_palette("Moonrise3", 3), guide = FALSE) +
    scale_fill_manual(values = wes_palette("Moonrise3", 3)) +
    theme(
      axis.text.x = element_blank()
    ) +
    labs(
      "x" = "Index",
      "y" = "log(Eigenvalue)"
    )

  for (i in seq_along(all_plots)) {
    cat(sprintf("saving %s\n", names(all_plots)[i]))
    dir.create(output_dir, recursive = TRUE)
    ggsave(
      file = sprintf("%s/figure-%s.png", output_dir, names(all_plots)[i]),
      all_plots[[i]],
      ...
    )
  }
  all_plots
}
