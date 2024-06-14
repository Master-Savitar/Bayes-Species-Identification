# load aggregated data firstly

all_roc_data <- rbind(mv_roc_data,
                      base_roc_data,
                      base_hierarchy_roc_data,
                      dpbmm_roc_data,
                      dpbmm_hierarchy_roc_data)

all_roc_data[all_roc_data$Algorithm == "DPBMM", "Algorithm"] <- "DP-BMM"
all_roc_data[all_roc_data$Algorithm == "DPBMM-Hierarchy", "Algorithm"] <- "DP-BMM-Hierarchy"

ggplot(all_roc_data, aes(x = FPR, y = TPR, color = Algorithm)) + 
  geom_line() + 
  labs(x = "False positive rate", y = "True positive rate") + 
  ylim(0.60, 1.0)
  
ggsave("roc.pdf", plot = p, width = 18, height = 12, units = "cm")
