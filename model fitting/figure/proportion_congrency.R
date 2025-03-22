
if (!require(R.matlab)) {
  install.packages("R.matlab")
  library(R.matlab)
}

if (!require(ggplot2)) {
  install.packages("ggplot2")
  library(ggplot2)
}

if (!require(ggsignif)) {
  install.packages("ggsignif")
  library(ggsignif)
}


file_path <- "F:/behavior_result/graph_data/RT_proportion_congruency.mat"
data_mat <- readMat(file_path)

alpha_vio <- data_mat$diffs.portition0 * 1000
alpha_stable <- data_mat$diffs.portition1 * 1000
data <- data.frame(
  Value = c(alpha_vio, alpha_stable),
  Condition = rep(c("20%Incongruent\nproportion", "80%Incongruent\nproportion"), each = length(alpha_vio))
)

p <- ggplot(data, aes(x = Condition, y = Value, fill = Condition)) +
  geom_violin(trim = TRUE, color = "black", fill = NA, width = 0.8) +  
  geom_jitter(width = 0.2, aes(color = Condition), size = 3, alpha = 0.8) +  
  geom_boxplot(width = 0.2, fill = NA, color = "black", outlier.shape = NA) +  
  scale_fill_manual(values = c("20%Incongruent\nproportion" = NA, "80%Incongruent\nproportion" = NA)) +  
  scale_color_manual(values = c("20%Incongruent\nproportion" = "#E69F00", "80%Incongruent\nproportion" = "#56B4E9")) +  
  labs(x = NULL,  
       y = "Stroop Effect (ms)") +  
  theme_minimal() +  
  theme(
    text = element_text(size = 12),
    legend.position = "none",  
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    axis.line.x = element_line(color = "black", size = 1),  
    axis.line.y = element_line(color = "black", size = 1),  
    axis.ticks = element_line(color = "black", size = 1.5),  
    axis.title.y = element_text(size = 20, family = "serif", face = "bold", margin = margin(r = 10)),  
    axis.text.y = element_text(size = 17, family = "serif", face = "bold"),  
    axis.text.x = element_text(size = 17, family = "serif", face = "bold"),  
    plot.title = element_blank()  
  ) +
  scale_y_continuous(breaks = seq(-0.02*1000, 0.16 * 1000, by = 0.02 * 1000), limits = c(-0.02 * 1000-0.02*1000/7, 0.16 * 1000+0.02*1000/7), expand = c(0, 0)) +  
  geom_signif(comparisons = list(c("20%Incongruent\nproportion", "80%Incongruent\nproportion")),
              map_signif_level = TRUE, textsize = 8, vjust = -0.2)  

print(p)

ggsave(filename = "C:/Users/wcy13/graph/Incongruent_proportion.png", 
       plot = p, 
       dpi = 1500,   
       width = 5.2,   
       height = 6,  
       units = "in",  
       bg = "white")  
