# load aggregated data firstly
# ......

# read in "gold standard"
parent_directory <- dirname(normalizePath("."))
gold_annotation <- readRDS(file.path(parent_directory, "dat/gold_annotation.rds")) %>% as.factor()

thresholds <- seq(0.15, 1, by = 0.05)
num_points <- length(thresholds)
mv_sen_spe_acc <- data.frame(
  threshold = thresholds, sensitivity = numeric(num_points),
  specificity = numeric(num_points), accuracy = numeric(num_points))
base_sen_spe_acc <- data.frame(
  threshold = thresholds, sensitivity = numeric(num_points),
  specificity = numeric(num_points), accuracy = numeric(num_points))
base_hierarchy_sen_spe_acc <- data.frame(
  threshold = thresholds, sensitivity = numeric(num_points),
  specificity = numeric(num_points), accuracy = numeric(num_points))
dpbmm_sen_spe_acc <- data.frame(
  threshold = thresholds, sensitivity = numeric(num_points),
  specificity = numeric(num_points), accuracy = numeric(num_points))
dpbmm_hierarchy_sen_spe_acc <- data.frame(
  threshold = thresholds, sensitivity = numeric(num_points),
  specificity = numeric(num_points), accuracy = numeric(num_points))

for (i in 1:length(thresholds)) {
  threshold <- thresholds[i]
  # all methods' prediction
  mv_predict <- ifelse(mv_probs >= threshold, 1, 0) %>% factor(levels = c(0, 1))
  base_predict <- ifelse(base_probs >= threshold, 1, 0) %>% factor(levels = c(0, 1))
  base_hierarchy_predict <- ifelse(base_hierarchy_probs >= threshold, 1, 0) %>% factor(levels = c(0, 1))
  dpbmm_predict <- ifelse(dpbmm_probs >=  threshold, 1, 0) %>% factor(levels = c(0, 1))
  dpbmm_hierarchy_predict <- ifelse(dpbmm_hierarchy_probs >= threshold, 1, 0) %>% factor(levels = c(0, 1))
  
  mv_sen_spe_acc[i, "specificity"] <- sensitivity(mv_predict, gold_annotation)
  base_sen_spe_acc[i, "specificity"] <- sensitivity(base_predict, gold_annotation)
  base_hierarchy_sen_spe_acc[i, "specificity"] <- sensitivity(base_hierarchy_predict, gold_annotation)
  dpbmm_sen_spe_acc[i, "specificity"] <- sensitivity(dpbmm_predict, gold_annotation)
  dpbmm_hierarchy_sen_spe_acc[i, "specificity"] <- sensitivity(dpbmm_hierarchy_predict, gold_annotation)
  
  mv_sen_spe_acc[i, "sensitivity"] <- specificity(mv_predict, gold_annotation)
  base_sen_spe_acc[i, "sensitivity"] <- specificity(base_predict, gold_annotation)
  base_hierarchy_sen_spe_acc[i, "sensitivity"] <- specificity(base_hierarchy_predict, gold_annotation)
  dpbmm_sen_spe_acc[i, "sensitivity"] <- specificity(dpbmm_predict, gold_annotation)
  dpbmm_hierarchy_sen_spe_acc[i, "sensitivity"] <- specificity(dpbmm_hierarchy_predict, gold_annotation)
  
  mv_sen_spe_acc[i, "accuracy"] <- (mv_predict == gold_annotation) %>% mean()
  base_sen_spe_acc[i, "accuracy"] <- (base_predict == gold_annotation) %>% mean()
  base_hierarchy_sen_spe_acc[i, "accuracy"] <- (base_hierarchy_predict == gold_annotation) %>% mean()
  dpbmm_sen_spe_acc[i, "accuracy"] <- (dpbmm_predict == gold_annotation) %>% mean()
  dpbmm_hierarchy_sen_spe_acc[i, "accuracy"] <- (dpbmm_hierarchy_predict == gold_annotation) %>% mean()
}

all_sen_spe_acc <- rbind(mv_sen_spe_acc,
                         base_sen_spe_acc,
                         base_hierarchy_sen_spe_acc,
                         dpbmm_sen_spe_acc,
                         dpbmm_hierarchy_sen_spe_acc)
all_sen_spe_acc["Algorithm"] <- c(
  rep("MV", num_points), rep("Base", num_points), rep("Base-Hierarchical", num_points),
  rep("DP-BMM", num_points), rep("DP-BMM-Hierarchical", num_points))


library(reshape2)

p <- all_sen_spe_acc %>% melt(id.vars = c("threshold", "Algorithm")) %>%
  ggplot(aes(x = threshold, y = value, color = Algorithm)) + 
  geom_line(size = 1.2) + 
  facet_wrap(~ variable, scales = "free_y", ncol = 3, nrow = 1) + 
  labs(x = "Threshold", y = "Value") + 
  theme(legend.title = element_blank())
ggsave("sen_spe_acc.pdf", plot = p, width = 36, height = 15, units = "cm")

