library(ggplot2)
library(R.matlab)

wavelet_data_path <- "F:/ningbo/EEG/finished"
behavior_data_path <- "F:/ningbo/behavior"
behavior_result_path <- "F:/ningbo/behavior_result"

files <- list.files(wavelet_data_path, pattern = "*.set")
subjects_results <- vector("list", length(files))

precent_stable80 <- rep(0.8, 81)
trials_stable80 <- 0:80
precent_stable20 <- rep(0.2, 81)
trials_stable20 <- 0:80
precent_vio80 <- c(rep(0.8, 21), rep(0.2, 21), rep(0.8, 21), rep(0.2, 21))
trials_vio80 <- c(0:20, 20:40, 40:60, 60:80)
precent_vio20 <- c(rep(0.2, 21), rep(0.8, 21), rep(0.2, 21), rep(0.8, 21))
trials_vio20 <- c(0:20, 20:40, 40:60, 60:80)
# choose any subject
i <- 2
subject <- substr(files[i], 1, 3)
id <- as.numeric(subject)

load_file_by_id <- function(folder_path, id) {
  files <- list.files(folder_path, pattern = "\\.mat$", full.names = TRUE)
  for (file in files) {
    filename <- basename(file)
    if (grepl(as.character(id), filename)) {
      data <- readMat(file)
      return(data)
    }
  }
  warning("No file found for ID: ", id)
  return(NULL)
}

data <- load_file_by_id(behavior_data_path, id)

if (is.null(data)) {
  stop("Failed to load data for ID: ", id)
}

data$runtype <- unlist(data$runtype)
alpha_data <- readMat(file.path(behavior_result_path, subject, "alpha.mat"))
alpha <- alpha_data$alpha
f <- alpha_data$f
congruency <- alpha_data$congruency
condition <- numeric(0)
trialNO <- numeric(0)

for (k in 1:length(data$runtype)) {
  runtype <- data$runtype[k]
  if (runtype == "0.8" & sum(congruency[(80*(k-1)+1):(80*k)]) != 40) {
    condition <- c(condition, precent_stable80)
    trialNO <- c(trialNO, trials_stable80 + 80 * (k - 1))
  } else if (runtype == "0.2" & sum(congruency[(80*(k-1)+1):(80*k)]) != 40) {
    condition <- c(condition, precent_stable20)
    trialNO <- c(trialNO, trials_stable20 + 80 * (k - 1))
  } else if (runtype == "0.8" & sum(congruency[(80*(k-1)+1):(80*k)]) == 40) {
    condition <- c(condition, precent_vio20)
    trialNO <- c(trialNO, trials_vio20 + 80 * (k - 1))
  } else if (runtype == "0.2" & sum(congruency[(80*(k-1)+1):(80*k)]) == 40) {
    condition <- c(condition, precent_vio80)
    trialNO <- c(trialNO, trials_vio80 + 80 * (k - 1))
  }
}

realtrialNO <- 1:(length(alpha) / 2)
trialNO <- trialNO[1:(length(trialNO) / 2)]
congruency[congruency == -1] <- 0
PE <- abs(congruency - f)
data_f <- data.frame(realtrialNO = realtrialNO, f = f[1:320])
data_condition <- data.frame(trialNO = trialNO, condition = condition[1:length(trialNO)])
custom_linetype <- "42"  

p <- ggplot() +
  geom_line(data = data_f, aes(x = realtrialNO, y = f, color = "Predicted Conflict Level"), size = 2) +
  geom_line(data = data_condition, aes(x = trialNO, y = condition, color = "Proportion of Incongruent Trials"), 
            linetype = custom_linetype, size = 2) +  
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1.02), expand = c(0, 0)) +
  scale_x_continuous(
    breaks = seq(0, max(c(realtrialNO, trialNO)) + 10, by = 80),
    limits = c(0, max(c(realtrialNO, trialNO)) + max(c(realtrialNO, trialNO)) * 0.0095),
    expand = c(0, 0)
  ) +
  scale_color_manual(values = c("Predicted Conflict Level" = "black", "Proportion of Incongruent Trials" = "#5DADE2")) +
  theme_minimal() +
  theme(
    axis.line.x = element_line(color = "black", size = 1),
    axis.line.y = element_line(color = "black", size = 1),
    axis.ticks = element_line(color = "black", size = 1.3),
    axis.ticks.length = unit(0.17, "cm"),
    axis.title = element_text(size = 20, family = "serif", face = "bold"),
    axis.title.y = element_text(margin = margin(r = 15), size = 20, family = "serif", face = "bold"),
    axis.text = element_text(size = 17, family = "serif", face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(1, 1, 1, 1, "cm"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.title = element_blank(),
    legend.key.width = unit(3, "cm"),
    legend.text = element_text(size = 15, family = "serif", face = "bold"),
    legend.spacing.x = unit(1, 'cm')
  ) +
  guides(color = guide_legend(override.aes = list(
    linetype = c("solid", custom_linetype),  
    size = c(2, 2)
  ))) +
  labs(x = "Trials", y = "Predicted Conflict Level")

print(p)

ggsave(filename = "C:/Users/wcy13/graph/f_track.png", 
       plot = p, 
       width = 9.98,   
       height = 5.51,  
       dpi = 1500)    
