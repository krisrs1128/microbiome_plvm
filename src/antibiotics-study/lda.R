#! /usr/bin/env Rscript

# File description -------------------------------------------------------------
# RStan code to run LDA on antibiotics dataset.
# Based on
# https://github.com/stan-dev/stan/releases/download/v2.12.0/stan-reference-2.12.0.pdf
# page 157

library("argparser")
parser <- arg_parser("Perform LDA on the antibiotics dataset")
parser <- add_argument(parser, "--subject", help = "Subject on which to perform analysis", default = "F")
argv <- parse_args(parser)

## ---- setup ----
library("rstan")
library("reshape2")
library("plyr")
library("dplyr")
library("ggplot2")
library("phyloseq")
library("RColorBrewer")
library("ggscaffold")
library("feather")
source("./posterior_check_funs.R")
dir.create("../../data/fits/", recursive = TRUE)
dir.create("../../data/figure-input/", recursive = TRUE)
dir.create("../../doc/figure/", recursive = TRUE)
set.seed(11242016)

## ---- get-data ----
abt <- get(load("../../data/antibiotics-study/abt.rda"))
abt <- abt %>%
  subset_samples(ind == argv$subject)

releveled_sample_data <- abt %>%
  sample_data %>%
  mutate(
    condition = revalue(
      condition,
      c("Pre Cp" = "Pre",
        "1st Cp" = "1st Course",
        "1st WPC" = "1st Course",
        "2nd Cp" = "2nd Course",
        "2nd WPC" = "2nd Course",
        "Post Cp" = "Post")
    )
  )
rownames(releveled_sample_data) <- abt %>%
  sample_data %>%
  rownames
sample_data(abt) <- releveled_sample_data

## ---- histograms ---
transformed_counts <- data_frame(
  count = c(get_taxa(abt), asinh(get_taxa(abt))),
  transformation = c(
    rep("original", ntaxa(abt) * nsamples(abt)),
    rep("asinh", ntaxa(abt) * nsamples(abt))
  )
)

p <- ggplot(transformed_counts) +
  geom_histogram(aes(x = count)) +
  facet_grid(. ~ transformation, scale = "free_x") +
  min_theme(list(text_size = 8, subtitle_size = 12))
ggsave("../../doc/figure/histograms-1.pdf", p)

## ---- heatmaps ----
x_order <- names(sort(taxa_sums(abt)))
y_order <- names(sort(sample_sums(abt)))
ordered_map <- function(x) {
  ggheatmap(
    x %>%
    melt(value.name = "fill", varnames = c("x", "y")),
    list("x_order" = x_order, "y_order" = y_order)
  ) +
    min_theme(list(text_size = 0)) +
    labs(x = "Sample", y = "Microbe")
}

p <- ordered_map(get_taxa(abt)) + ggtitle("Raw")
ggsave("../../doc/figure/heatmaps-1.pdf", p)

p <- ordered_map(asinh(get_taxa(abt))) + ggtitle("asinh")
ggsave("../../doc/figure/heatmaps-2.pdf", p)

## ---- lda ----
x <- t(get_taxa(abt))
dimnames(x) <- NULL
stan_data <- list(
  K = 4,
  V = ncol(x),
  D = nrow(x),
  n = x,
  alpha = rep(1, 4),
  gamma = rep(0.5, ncol(x))
)

m <- stan_model(file = "../stan/lda_counts.stan")
n_iter <- 1000
stan_fit <- vb(m, stan_data, iter = 2 * n_iter)
save(
  stan_fit,
  file = sprintf("../../data/fits/lda-%s-%s.rda", argv$subject, gsub("[:|| ||-]", "", Sys.time()))
)
samples <- rstan::extract(stan_fit)
rm(stan_fit)

## ---- extract_beta ----
# underlying RSV distributions
beta_logit <- samples$beta

for (i in seq_len(n_iter)) {
  for (k in seq_len(stan_data$K)) {
    beta_logit[i, k, ] <- log(beta_logit[i, k, ])
    beta_logit[i, k, ] <- beta_logit[i, k, ] - mean(beta_logit[i, k, ])
  }
}

beta_hat <- beta_logit %>%
  melt(
    varnames = c("iterations", "topic", "rsv_ix"),
    value.name = "beta_logit"
  ) %>%
  as_data_frame()

beta_hat$rsv <- rownames(tax_table(abt))[beta_hat$rsv_ix]
taxa <- as_data_frame(tax_table(abt)@.Data)
taxa$rsv <- rownames(tax_table(abt))
taxa$Taxon_5[which(taxa$Taxon_5 == "")] <- taxa$Taxon_4[which(taxa$Taxon_5 == "")]

beta_hat <- beta_hat %>%
  left_join(taxa) %>%
  mutate(
    topic = paste("Topic", topic),
    Taxon_5 = stringr::str_extract(Taxon_5, "[^_]+")
  )

sorted_taxa <- names(sort(table(beta_hat$Taxon_5), decreasing = TRUE))
beta_hat$Taxon_5 <- factor(beta_hat$Taxon_5, levels = sorted_taxa)
beta_hat$rsv <- factor(beta_hat$rsv, levels = taxa$rsv)

## ---- extract_theta ----
theta_logit <- samples$theta
for (i in seq_len(n_iter)) {
  for (d in seq_len(stan_data$D)) {
    theta_logit[i, d, ] <- log(theta_logit[i, d, ])
    theta_logit[i, d, ] <- theta_logit[i, d, ] - mean(theta_logit[i, d, ])
  }
}

theta_hat <- theta_logit %>%
  melt(
    varnames = c("iteration", "sample", "topic"),
    value.name = "theta_logit"
  )

theta_hat$sample <- sample_names(abt)[theta_hat$sample]
sample_info <- sample_data(abt)
sample_info$sample <- rownames(sample_info)
theta_hat$topic <- paste("Topic", theta_hat$topic)

theta_hat <- theta_hat %>%
  left_join(sample_info, by = "sample")

## ---- visualize_lda_theta_heatmap ----
plot_opts <- list(
  "x" = "time",
  "y" = "topic",
  "fill" = "mean_theta",
  "y_order" = paste("Topic", stan_data$K:1)
)

p <- ggheatmap(
  theta_hat %>%
  group_by(topic, time) %>%
  summarise(mean_theta = mean(theta_logit, na.rm = TRUE)) %>%
  as.data.frame(),
  plot_opts
) +
  labs(fill = "g(theta)")
ggsave(
  sprintf("../../doc/figure/visualize_lda_theta_heatmap-%s.pdf"),
  p, width = 7, height = 0.9
)

## ---- visualize_lda_theta_boxplot ----
p <- ggplot(theta_hat) +
  geom_boxplot(
    aes(x = as.factor(time), y = theta_logit),
    fill = "#C9C9C9",
    outlier.size = 0.05,
    size = 0.1,
    notchwidth = 0.1
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(3)) +
  min_theme(list(border_size = 0.7, text_size = 10, subtitle_size = 11)) +
  facet_grid(topic ~ condition, scales = "free_x", space = "free_x") +
  geom_hline(yintercept = 0, alpha = 0.4, size = 0.5, col = "#999999") +
  labs(x = "Time", y = expression(paste("g(", theta[k], ")"))) +
  theme(legend.position = "none") +
  scale_x_discrete(breaks = seq(1, 60, by = 10) - 1)
ggsave(
  sprintf("../../doc/figure/visualize_lda_theta_boxplot-%s.pdf", argv$subject),
  p, width = 6, height = 2.9
)

## ---- visualize_lda_beta ----
beta_summary <- beta_hat %>%
  group_by(rsv_ix, topic) %>%
  summarise(
    beta_median = median(beta_logit),
    Taxon_5 = Taxon_5[1],
    beta_upper = quantile(beta_logit, 0.975),
    beta_lower = quantile(beta_logit, 0.025)
  ) %>%
  filter(Taxon_5 %in% levels(beta_hat$Taxon_5)[1:5])
beta_summary$rsv_ix <- rep(seq_len(nrow(beta_summary) / 4), each = 4)

p <- ggplot(beta_summary) +
  geom_hline(yintercept = 0, alpha = 0.4, size = 0.5, col = "#999999") +
  geom_point(aes(x = rsv_ix, y = beta_median, col = Taxon_5), size = 0.1) +
  geom_errorbar(
    aes(x = rsv_ix, alpha = beta_upper, ymax = beta_upper, ymin = beta_lower, col = Taxon_5),
    size = 0.4
  ) +
  scale_color_brewer(palette = "Set2") +
  scale_alpha(range = c(0.01, 1), breaks = c(1, 2, 3), guide = FALSE) + ## larger values have darker CIs
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(breaks = scales::pretty_breaks(3), limits = c(-5, 12)) +
  facet_grid(topic ~ .) +
  labs(x = "Species", y = expression(paste("g(", beta[k], ")")), col = "Family") +
  theme(
    panel.border = element_rect(fill = "transparent", size = 0.75),
    axis.text.x = element_blank(),
    strip.text.x = element_blank(),
    legend.position = "bottom"
  )
ggsave(
  sprintf("../../doc/figure/visualize_lda_beta-%s.pdf", argv$subject),
  p, width = 6, height = 3.5
)

## ---- posterior-checks ----
checks_data <- posterior_checks_input(
  x,
  samples$x_sim,
  sprintf("../../data/figure-input/lda-%s", argv$subject)
)

## ---- js-input ----
colnames(beta_summary) <- c("ix", "topic", "median", "fill", "upper", "lower")
cat(
  sprintf("var beta = %s", jsonlite::toJSON(beta_summary, auto_unbox = TRUE)),
  file = sprintf("../../data/antibiotics-study/lda_beta-%s.js", argv$subject)
)
