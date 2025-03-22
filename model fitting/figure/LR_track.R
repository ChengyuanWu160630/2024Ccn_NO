library(R.matlab)
library(ggplot2)


file_path <- "F:/behavior_result/graph_data/LR/LR_block.mat"
data <- readMat(file_path)
stab_block <- data$stab.block
vio_block <- data$vio.block
stab_df <- as.data.frame(stab_block)
vio_df <- as.data.frame(vio_block)
stab_mean <- apply(stab_df, 2, mean)
stab_sd <- apply(stab_df, 2, sd)
vio_mean <- apply(vio_df, 2, mean)
vio_sd <- apply(vio_df, 2, sd)
learning_rates <- 1:20
stab_data <- data.frame(LearningRate = learning_rates, Mean = stab_mean, SD = stab_sd, Group = "Stable block")
vio_data <- data.frame(LearningRate = learning_rates, Mean = vio_mean, SD = vio_sd, Group = "Volatile block")
plot_data <- rbind(stab_data, vio_data)

p <- ggplot(plot_data, aes(x = LearningRate, y = Mean, color = Group, shape = Group)) +
  geom_line(linetype = "solid", size = 1.2) +
  geom_errorbar(aes(ymin = Mean - SD/5, ymax = Mean + SD/5), width = 0.2, size = 0.8, color = "black") +
  geom_point(size = 2) +
  scale_color_manual(values = c("Stable block" = "#1F77B4", "Volatile block" = "#D62728")) +
  scale_shape_manual(values = c("Stable block" = 16, "Volatile block" = 16)) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black", size = 1),
    axis.ticks = element_line(color = "black", size = 1.3),
    axis.ticks.length = unit(0.17, "cm"),
    axis.title = element_text(size = 20, family = "serif", face = "bold"),
    axis.text = element_text(size = 17,family = "serif", face = "bold"),
    axis.title.y = element_text(size = 20, family = "serif", face = "bold", margin = margin(r = 10)),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(1, 1, 1, 1, "cm"),
    legend.position = c(0.85, 0.94),
    legend.direction = "vertical",
    legend.title = element_blank(),
    legend.text = element_text(size = 15, family = "serif", hjust = 0.5),  # 将图例文字居中对齐
    legend.key.size = unit(1.5, "cm")
  ) +
  guides(
    color = guide_legend(
      override.aes = list(size = 2),
      keyheight = unit(0.7, "cm"),
      default.unit = "cm"
    ),
    shape = guide_legend(
      override.aes = list(size = 2),
      keyheight = unit(0.7, "cm"),
      default.unit = "cm"
    )
  ) +
  labs(x = "Trials", y = "Learning Rate") +
  scale_y_continuous(breaks = seq(0.03, 0.038, by = 0.002), limits = c(0.03 - 0.002 / 7 , 0.038 + 0.002 / 7), expand = c(0, 0)) +
  scale_x_continuous(breaks = seq(2, 20, by = 2))


print(p)
ggsave(filename = "C:/Users/wcy13/graph/LR_track.png", 
       plot = p, 
       dpi = 1500,   
       width = 8,   
       height = 5,  
       units = "in",  
       bg = "white")  