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
library("tidyr")
library("ggplot2")
library("RColorBrewer")
library("phyloseq")
library("feather")
library("ggscaffold")
source("./posterior_check_funs.R")
dir.create("../../data/fits/", recursive = TRUE)
dir.create("../../doc/figure/", recursive = TRUE)
theme_set(min_theme(list(text_size = 7, subtitle_size = 9)))
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
stan_fit <- vb(
  m,
  data = stan_data,
  iter = 6000,
  output_samples = 1000,
  eta = 0.1,
  adapt_engaged = FALSE
)
save(
  stan_fit,
  file = sprintf("../../data/fits/unigram-%s-%s.rda", argv$subject, gsub("[:|| ||-]", "", Sys.time()))
)
samples <- rstan::extract(stan_fit)
rm(stan_fit)

## ---- prepare-mu ----
taxa <- as_data_frame(tax_table(abt)@.Data)
taxa$rsv <- rownames(tax_table(abt))
taxa$Taxon_5[which(taxa$Taxon_5 == "")] <- taxa$Taxon_4[which(taxa$Taxon_5 == "")]
taxa$Taxon_5 <- taxa$Taxon_5 %>%
  revalue(
    c("Bacteroidaceae_Bacteroides" = "Bacteroidaceae",
      "Peptostreptococcaceae_1" = "Peptostreptococcaceae")
  )

## subset to small number of times, and center the mus
keep_times <- seq(10, 26, by = 5)
mu <- samples$mu[, keep_times, ]
for (i in seq_len(ncol(mu))) {
  mu[, i,] <- mu[, i, ] - mean(mu[, i,])
}

## ---- plot-data ----
## prepare summary statistics to plot
mu_summary <- apply(mu, c(2, 3), quantile, c(0.025, 0.5, 0.975)) %>%
  melt(
    varnames = c("quantile", "time", "rsv_ix"),
    value.name = "mu"
  ) %>%
  spread(quantile, mu) %>%
  mutate(
    time = keep_times[time]
  ) %>%
  left_join(sample_data(abt)[, c("time", "condition")]) %>%
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

colnames(mu_summary) <- c("time", "rsv_ix", "lower", "median", "upper", "condition")
mu_summary$rsv <- factor(
  taxa[mu_summary$rsv_ix, ]$rsv,
  levels = rownames(tax_table(abt))
)

mu_summary <- mu_summary %>%
  left_join(taxa[, c("rsv", "Taxon_5")])
sorted_taxa <- names(sort(table(mu_summary$Taxon_5), decreasing = TRUE))
mu_summary$Taxon_5 <- factor(
  mu_summary$Taxon_5,
  levels = c(sorted_taxa, "other")
)

mu_summary$Taxon_5[!(mu_summary$Taxon_5 %in% sorted_taxa[1:7])] <- "other"
p <- ggplot(mu_summary) +
  geom_hline(yintercept = 0, alpha = 0.4, size = 0.5, col = "#999999") +
  geom_point(aes(x = rsv_ix, y = median, col = Taxon_5), size = 0.1) +
  geom_errorbar(
    aes(x = rsv_ix, alpha = upper, ymax = upper, ymin = lower, col = Taxon_5),
    size = 0.4
  ) +
  scale_color_brewer(palette = "Set2") +
  scale_alpha(range = c(0.01, 1), breaks = c(1, 2, 3), guide = FALSE) + ## larger values have darker CIs
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(breaks = scales::pretty_breaks(3), limits = c(-3, 10)) +
  facet_grid(condition + time ~ ., scales = "free_x", space = "free_x") +
  labs(x = "Species", y = expression(mu[t]), col = "Family") +
  theme(
    panel.border = element_rect(fill = "transparent", size = 0.75),
    axis.text.x = element_blank(),
    strip.text.x = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 6.1)
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
