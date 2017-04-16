#! /usr/bin/env Rscript

# File description -------------------------------------------------------------
# This is an application of the dynamic unigram model to the antibiotics data

# Code block ------------------------------------------------------------------

## ---- setup ----
library("rstan")
library("data.table")
library("reshape2")
library("plyr")
library("dplyr")
library("ggplot2")
library("RColorBrewer")
library("phyloseq")
library("treelapse")
library("feather")
library("ggscaffold")
source("./posterior_checks.R")
set.seed(11242016)

 theme_set(
  min_theme(list(text_size = 8, subtitle_size = 12))
)

softmax <- function(x) {
  exp(x) / sum(exp(x))
}

## ---- get-data ----
data(abt)
abt <- abt %>%
  filter_taxa(function(x) sum(x != 0) > .45 * nsamples(abt), prune = TRUE) %>%
  subset_samples(ind == "F")

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
  file = sprintf("../../data/fits/unigram-%s.rda", gsub("[:|| ||-]", "", Sys.time()))
)
samples <- rstan::extract(stan_fit)
rm(stan_fit)

## ---- prepare-mu ----
taxa <- data.table(
  rsv = rownames(tax_table(abt)),
  tax_table(abt)@.Data
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
  taxa[mu_hat$rsv_ix]$rsv,
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
ggsave("../../doc/figure/unigramseries-1.pdf", p)

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
p <- ggboxplot(
  mu_hat %>%
  filter(
    Taxon_5 %in% levels(mu_hat$Taxon_5)[1:4],
    time %in% seq(10, 26, by = 5)
  ) %>%
  as.data.frame(),
  plot_opts
) +
  geom_hline(yintercept = 0, alpha = 0.4, size = 0.5, col = "#999999") +
  scale_y_continuous(breaks = scales::pretty_breaks(3), limits = c(-9, 7)) +
  facet_grid(condition + time ~ Taxon_5, scales = "free_x", space = "free_x") +
  labs(y = expression(mu[t]), fill = "Family") +
  theme(
    strip.text.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title = element_text(size = 11),
    legend.position = "bottom",
  )
ggsave("../../doc/figure/antibiotics_unigram_mu.pdf", p, width = 6, height = 3.5)

## ---- posterior-checks ----
counts_data_checker(x, samples$x_sim, "../../doc/figure/unigram_post_checks")
