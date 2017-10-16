library("jsonlite")
df <- expand.grid(
  D = c(20, 100),
  V = c(325, 650),
  K = 2,
  N = c(1625, 3250, 6500),
  a0 = 1,
  b0 = 1
  sigma0 = 1
)

df$id <- seq_len(nrow(df))
unigram_dir <- file.path(
  Sys.getenv("MICROBIOME_PLVM_DIR"),
  "src",
  "sim",
  "unigram"
)

cat(
  toJSON(df, auto_unbox = TRUE),
  file = file.path(unigram_dir, "pipeline", "experiment.json")
)
