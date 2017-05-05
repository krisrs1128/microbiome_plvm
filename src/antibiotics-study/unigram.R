#! /usr/bin/env Rscript

# File description -------------------------------------------------------------
# This is an application of the dynamic unigram model to the antibiotics data

library("argparser")
parser <- arg_parser("Perform dynamic unigrams on the antibiotics dataset")
parser <- add_argument(parser, "--subject", help = "Subject on which to perform analysis", default = "F")
argv <- parse_args(parser)

## ---- setup ----
library("rstan")
library("reshape2")
library("plyr")
library("dplyr")
library("ggplot2")
library("RColorBrewer")
library("phyloseq")
library("feather")
library("ggscaffold")
source("./posterior_check_funs.R")
dir.create("../../data/fits/", recursive = TRUE)
dir.create("../../doc/figure/", recursive = TRUE)
set.seed(11242016)

softmax <- function(x) {
  exp(x) / sum(exp(x))
}

## ---- get-data ----
abt <- get(load("../../data/antibiotics-study/abt.rda"))
abt <- abt %>%
  subset_samples(ind == argv$subject)

## ---- run-model ----
times <- sample_data(abt)$time
x <- t(get_taxa(abt))
dimnames(x) <- NULL
stan_data <- list(
  N = nrow(x),
 V = ncol(x),
  T = length(times),
  times = times,
  times_mapping = times,
  x = x,
  a0 = 0.5,
  b0 = 0.5
)

m <- stan_model("../stan/unigram.stan")
stan_fit <- vb(m, data = stan_data)
save(
  stan_fit,
  file = sprintf("../../data/fits/unigram-%s-%s.rda", argv$subject, gsub("[:|| ||-]", "", Sys.time()))
)
samples <- rstan::extract(stan_fit)
rm(stan_fit)

## ---- prepare-mu ----
taxa <- cbind(
  rsv = rownames(tax_table(abt)),
  as_data_frame(tax_table(abt)@.Data)
)
taxa$Taxon_5[which(taxa$Taxon_5 == "")] <- taxa$Taxon_4[which(taxa$Taxon_5 == "")]
taxa$Taxon_5 <- taxa$Taxon_5 %>%
  revalue(c("Bacteroidaceae_Bacteroides" = "Bacteroides"))

## center the mus
mu <- samples$mu
for (i in seq_len(stan_data$T)) {
  mu[, i,] <- mu[, i, ] - mean(mu[, i,])
}

mu_hat <- samples$mu %>%
  melt(
    varnames = c("iteration", "time", "rsv_ix"),
    value.name = "mu"
  )

mu_hat$rsv <- rownames(otu_table(abt))[mu_hat$rsv_ix]
mu_hat$time <- times[mu_hat$time]
mu_hat <- mu_hat %>%
  left_join(taxa) %>%
  left_join(sample_data(abt)[, c("time", "condition")]) %>%
  group_by(time) %>%
  mutate(
    prob = softmax(mu),
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

group_order <- sort(table(taxa$Taxon_5), decreasing = TRUE)
mu_hat$Taxon_5 <- factor(mu_hat$Taxon_5, levels = names(group_order))
mu_hat$rsv <- factor(
  taxa[mu_hat$rsv_ix, ]$rsv,
  levels = rownames(tax_table(abt))
)

## ---- unigramseries ----
plot_opts <- list(
  "x" = "time",
  "y" = "mean_mu",
  "col" = "Taxon_5",
  "facet_terms" = c("Taxon_5", "."),
  "alpha" = 0.4,
  "group" = "rsv"
)
p <- gglines(
  mu_hat %>%
  filter(Taxon_5 %in% levels(mu_hat$Taxon_5)[1:4]) %>%
  group_by(rsv, time) %>%
  summarise(mean_mu = mean(mu), Taxon_5 = Taxon_5[1]) %>%
  as.data.frame(),
  plot_opts
) +
  scale_y_continuous(breaks = scales::pretty_breaks(3)) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 1))) +
  theme(
    strip.text.y = element_blank(),
    legend.position = "bottom"
  )
ggsave(
  sprintf("../../doc/figure/unigramseries-%s.pdf", argv$subject),
  p
)

## ---- unigramboxplots ----
plot_opts <- list(
  "x" = "rsv",
  "y" = "mu",
  "fill" = "Taxon_5",
  "col" = "Taxon_5",
  "outlier.shape" = NA,
  "alpha" = 1,
  "size" = 0.4,
  "col_colors" = brewer.pal(8, "Set2"),
  "fill_colors" = brewer.pal(8, "Set2"),
  "theme_opts" = list(border_size = 0.7, text_size = 10, subtitle_size = 10)
)

mu_summary <- mu_hat %>%
  group_by(rsv, topic) %>%
  summarise(
    mu_median = median(mu),
    Taxon_5 = Taxon_5[1],
    mu_upper = quantile(mu, 0.975),
    mu_lower = quantile(mu, 0.025)
  ) %>%
  filter(
    Taxon_5 %in% levels(mu_hat$Taxon_5)[1:5],
    time %in% seq(10, 26, by = 5)
  )
mu_summary$rsv_ix <- rep(seq_len(nrow(mu_summary) / 4), each = 4)

theme_set(min_theme(list(text_size = 10, subtitle_size = 10)))
p <- ggplot(mu_summary) +
  geom_hline(yintercept = 0, alpha = 0.4, size = 0.5, col = "#999999") +
  geom_point(aes(x = rsv_ix, y = mu_median, col = Taxon_5), size = 0.1) +
  geom_errorbar(
    aes(x = rsv, alpha = mu_upper, ymax = mu_upper, ymin = mu_lower, col = Taxon_5),
    size = 0.4
  ) +
  scale_color_brewer(palette = "Set2") +
  scale_alpha(range = c(0.01, 1), breaks = c(1, 2, 3), guide = FALSE) + ## larger values have darker CIs
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(breaks = scales::pretty_breaks(3), limits = c(-9, 7)) +
  facet_grid(condition + time ~ Taxon_5, scales = "free_x", space = "free_x") +
  labs(x = "Species", y = expression(mu[t]), fill = "Family") +
  theme(
    panel.border = element_rect(fill = "transparent", size = 0.75),
    axis.text.x = element_blank(),
    strip.text.x = element_blank(),
    legend.position = "bottom"
  )

ggsave(
  sprintf("../../doc/figure/antibiotics_unigram_mu-%s.pdf", argv$subject),
  p, width = 6, height = 3.5
)

## ---- posterior-checks ----
checks_data <- posterior_checks_input(
  x,
  samples$x_sim,
  sprintf("../../data/figure-input/unigram-%s", argv$subject)
)
