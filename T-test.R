################################################################################
################################################################################
################################################################################
################################################################################
#                        RESAMPLING - IS THE DATA OK?                         #
################################################################################
################################################################################
################################################################################
################################################################################


# To assess whether there was any bias in our data due to the different number 
# of sampled individuals per species, we developed a resampling protocol to 
# randomly sample four individuals from each species and calculated the mean 
# temperature for each one of the nine measured timepoints (morning, midday and
# evening each day), and isotope values (δ13C and δ15N). For each species, one 
# thousand resampling iterations were performed. A t-test was employed to compare
# the original mean values for each variable with the corresponding resampled
# means for each species. Since there were
# no differences, all the subsequent analyses only used the original values. 



# ----------------- Reproducibility setup -----------------
set.seed(42)                                # fixed seed
if (getRversion() >= "3.6.0") RNGkind(sample.kind = "Rounding")

# ----------------- Load data -----------------
holocron <- read.table(file.choose(), header = TRUE)

# ----------------- Load packages -----------------
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
library(tibble)

# ----------------- Coerce types -----------------
holocron <- holocron %>%
  mutate(
    Spp    = as.factor(Spp),
    I15N   = as.numeric(I15N),
    I13C   = as.numeric(I13C),
    Body_g = as.numeric(Body_g),
    B_mm   = as.numeric(B_mm),
    Abd_mm = as.numeric(Abd_mm),
    Cth_mm = as.numeric(Cth_mm),
    TA_H1  = as.numeric(TA_H1),
    TA_H2  = as.numeric(TA_H2),
    TA_H3  = as.numeric(TA_H3),
    TB_H1  = as.numeric(TB_H1),
    TB_H2  = as.numeric(TB_H2),
    TB_H3  = as.numeric(TB_H3),
    TC_H1  = as.numeric(TC_H1),
    TC_H2  = as.numeric(TC_H2),
    TC_H3  = as.numeric(TC_H3)
  )

# =========================================================
# RESAMPLING
# =========================================================

# Number of individuals per species per iteration
n_draw_per_spp <- 4L
n_iterations   <- 20000L    # set to 19000

vars <- c("TA_H1","TA_H2","TA_H3",
          "TB_H1","TB_H2","TB_H3",
          "TC_H1","TC_H2","TC_H3",
          "I13C","I15N")

# Split dataset by species
dat_by_spp <- split(holocron, holocron$Spp)

# Resolve draw sizes per species, capping at available n
draw_sizes <- imap_int(dat_by_spp, function(df, spp) {
  want <- if (length(n_draw_per_spp) == 1L) n_draw_per_spp else n_draw_per_spp[[spp]]
  min(nrow(df), want)
})

# ----------------- Create or load frozen plan -----------------
plan_path <- "resampling_plan.rds"

if (!file.exists(plan_path)) {
  index_draws <- imap(dat_by_spp, function(df, spp) {
    n_samp <- draw_sizes[[spp]]
    replicate(n_iterations, sample.int(nrow(df), n_samp), simplify = FALSE)
  })
  saveRDS(list(seed = .Random.seed,
               draw_sizes = draw_sizes,
               n_iterations = n_iterations,
               index_draws = index_draws),
          plan_path)
} else {
  plan <- readRDS(plan_path)
  index_draws <- plan$index_draws
  stopifnot(plan$n_iterations == n_iterations)
  stopifnot(identical(plan$draw_sizes, draw_sizes))
}

# ----------------- Run resampling using frozen indices -----------------
resampling_results <- imap_dfr(dat_by_spp, function(df, spp) {
  map_dfr(index_draws[[spp]], function(idx) {
    df[idx, , drop = FALSE] %>%
      summarise(across(all_of(vars), ~ mean(.x, na.rm = TRUE), .names = "{.col}_mean"))
  }) %>%
    mutate(Spp = spp, .before = 1)
})

# =========================================================
# Original species means
# =========================================================
original_means <- holocron %>%
  group_by(Spp) %>%
  summarise(across(all_of(vars), ~ mean(.x, na.rm = TRUE), .names = "{.col}_mean"),
            .groups = "drop")

# Average resampled means
resampling_means <- resampling_results %>%
  group_by(Spp) %>%
  summarise(across(ends_with("_mean"), ~ mean(.x, na.rm = TRUE),
                   .names = "{.col}_resampled"),
            .groups = "drop")

# Comparison (resampled – original)
comparison <- original_means %>%
  left_join(resampling_means, by = "Spp") %>%
  mutate(across(ends_with("_resampled"),
                ~ .x - get(sub("_resampled$", "", cur_column())),
                .names = "{.col}_diff"))

print(comparison)

# =========================================================
# Plots (change to any other if needed)
# =========================================================

# Example: TA_H1
ggplot(resampling_results, aes(x = TA_H1_mean, fill = Spp)) +
  geom_density(alpha = 0.5) +
  geom_vline(data = original_means,
             aes(xintercept = TA_H1_mean, color = Spp),
             linetype = "dashed", linewidth = 1) +
  labs(title = "Density of resampled TA_H1 means vs. original",
       x = "TA_H1 mean", y = "Density") +
  theme_minimal()

# Example: I13C
ggplot(resampling_results, aes(x = I13C_mean, fill = Spp)) +
  geom_density(alpha = 0.5) +
  geom_vline(data = original_means,
             aes(xintercept = I13C_mean, color = Spp),
             linetype = "dashed", linewidth = 1) +
  labs(title = "Density of resampled I13C means vs. original",
       x = "I13C mean", y = "Density") +
  theme_minimal()

# =========================================================
# T-tests
# =========================================================
do_ttest <- function(df, var, orig_mean) {
  x <- df[[paste0(var, "_mean")]]
  x <- x[!is.na(x)]
  if (length(x) < 2) return(tibble(statistic=NA, parameter=NA, p.value=NA, conf.low=NA, conf.high=NA))
  tt <- t.test(x, mu = orig_mean)
  tibble(statistic = unname(tt$statistic),
         parameter = unname(tt$parameter),
         p.value   = tt$p.value,
         conf.low  = unname(tt$conf.int[1]),
         conf.high = unname(tt$conf.int[2]))
}

ttests <- map_dfr(names(dat_by_spp), function(spp) {
  df_resamp <- resampling_results %>% filter(Spp == spp)
  df_orig   <- original_means %>% filter(Spp == spp)
  map_dfr(vars, function(v) {
    do_ttest(df_resamp, v, df_orig[[paste0(v, "_mean")]]) %>%
      mutate(Spp = spp, Variable = v, .before = 1)
  })
})

ttests_summary <- ttests %>%
  mutate(significant_0.05 = p.value < 0.05) %>%
  arrange(Spp, Variable)

print(ttests_summary)
