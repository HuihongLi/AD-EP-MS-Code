library(ggplot2)
library(dplyr)

data <- read.csv("plotJ20.csv", row.names = 1)

data_avg <- data %>%
  group_by(Age_months, Genotype) %>%
  summarise(mean_Cst7 = mean(Cst7), .groups = 'drop')


p <- ggplot(data, aes(x = Age_months, y = Cst7)) +
  geom_point(aes(color = Genotype), size = 3) +
  geom_line(data = data_avg, aes(x = Age_months, y = mean_Cst7, color = Genotype), linetype = "dashed", size = 1) +
  scale_color_manual(values = c("J20" = "#00A087", "WT" = "black")) +
  theme_minimal() +
  labs(x = "Age (months)", y = "Normalised counts", title = "Cst7") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = c(0.1, 0.9),  
    legend.direction = "vertical",
    legend.background = element_rect(fill = "white", color = NA),
    legend.key = element_rect(fill = "white", color = NA),
    legend.key.size = unit(1.5, "lines"),  
    legend.title = element_text(size = 10, face = "bold"),
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    panel.background = element_blank(),  
    axis.line.x = element_line(color = "black", size = 0.8),  
    axis.line.y = element_line(color = "black", size = 0.8),  
    axis.ticks.length = unit(-0.15, "cm"),  
    axis.ticks = element_line(color = "black", size = 0.8),  
    axis.text = element_text(color = "black")  # 轴文字颜色
  ) +
  guides(color = guide_legend(override.aes = list(shape = 15, linetype = 0, size = 5))) +  
  labs(color = "Genotype")

# 显示图像
print(p)
