library(ggplot2)
library(tidyr)
library(dplyr)
library(purrr)

#load("Data/eval_outer_NO_1_1_75_1_2_1_1_100.RData"); results_list=eval_outer

# Data wrangling
score_long <- map2_dfr(score_items, names(score_items), ~ {
  df <- as.data.frame(.x)
  df$row <- seq_len(nrow(df))
  df_long <- pivot_longer(df, -row, names_to = "method", values_to = "value")
  df_long$score_type <- .y
  df_long
})

score_long <- score_long %>% select(-row)

method_order <- c("glm", "gee", "re_nosig", "re_np", "lme4", "gamm",
                  "cop_n", "cop", "cop_j", "cop_g", "cop_f")
rename_model <- function(x) {
  main_map <- c(
    glm = "GLM",
    gee = "GEE",
    re_nosig = "GAMLSS",
    re_np = "GAMLSS NP",
    lme4 = "LME4",
    gamm = "GAMM"
  )
  cop_map <- c(
    cop = "GJRM (C)",
    cop_n = "GJRM (N)",
    cop_j = "GJRM (J)",
    cop_g = "GJRM (G)",
    cop_f = "GJRM (F)"
  )
  if(x %in% names(main_map)) return(main_map[x])
  if(x %in% names(cop_map)) return(cop_map[x])
  return(x)
}
method_label_order <- sapply(method_order, rename_model)
score_long <- score_long %>%
  mutate(method_label = sapply(as.character(method), rename_model),
         method_label = factor(method_label, levels = method_label_order))

# Rename score_type for facets
score_type_map <- c(
  es = "Energy Score",
  vs1 = "Variogram Score (p=1)",
  vs2 = "Variogram Score (p=2)",
  vs2_wt = "Variogram Score (p=2, Weighted)",
  vs2_wt_coronly = "Variogram Score (p=2, Correlated Obs Only)"
)
score_long <- score_long %>%
  mutate(score_type_label = score_type_map[score_type])

# Calculate lowest median per score_type
meds <- score_long %>%
  group_by(score_type, method_label) %>%
  summarize(med = median(value), .groups = "drop") %>%
  group_by(score_type) %>%
  summarize(min_median = min(med), .groups = "drop") %>%
  mutate(score_type_label = score_type_map[score_type])

# Plot
ggplot(score_long, aes(x = method_label, y = value, fill = score_type_label)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  geom_hline(data = meds, aes(yintercept = min_median), linetype = "dashed", color = "red", inherit.aes = FALSE) +
  facet_wrap(~score_type_label, scales = "free", ncol = 1) +
  theme_bw() +
  labs(title = "Boxplot for each score item by method",
       y = "Score Value", x = "Method", fill="Score Type")