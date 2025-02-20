library(ggplot2)

##Development of insular woodiness over time
# Simulate data
set.seed(42)
time1 <- seq(0, 10, length.out = 100)  # Time variable
y1 <- 1 / (1 + exp(- (time - 5))) + rnorm(100, 0, 0)  # Logistic function with slight noise

# Create the plot
ggplot(data.frame(time1, y1), aes(x = time1, y = y1)) +
  geom_line(color = "darkgreen", size = 1.2) +  # Smooth line
  theme_minimal() +
  labs(x = "Time",
       y = "Degree of insular woodiness",
       title = "Development of insular woodiness over time") +
  scale_x_continuous(breaks = c(0, 10), labels = c("Early", "Late")) +  # Custom X-axis labels
  scale_y_continuous(breaks = c(0, 1), labels = c("Low", "High")) +  # Custom Y-axis labels
  theme(plot.title = element_text(hjust = 0, size = 14, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.position = "none")

#Save with ggplot in 600 dpi and with a white background
ggsave("C:/Users/daank/OneDrive - University of Twente/Documents/Github/island_biology/plots/insular_woodiness.png", width = 6, height = 3.5, dpi = 600, bg = "white")


##Competitive strength with respect to native plant species
# Simulate data
set.seed(42)
time2 <- seq(0, 10, length.out = 100)  # Time variable
y2 <- 1 - (1 / (1 + exp(- (time - 5)))) + rnorm(100, 0, 0)  # Reversed logistic function

# Create the plot
ggplot(data.frame(time2, y2), aes(x = time2, y = y2)) +
  geom_line(color = "red", size = 1.2) +  # Smooth line
  theme_minimal() +
  labs(x = "Degree of insular woodiness",
       y = "Competitive strength",
       title = "Competitive strength of invasive VS native species") +
  scale_x_continuous(breaks = c(0, 10), labels = c("Low", "High")) +  # Custom X-axis labels
  scale_y_continuous(breaks = c(0, 1), labels = c("Low", "High")) +  # Custom Y-axis labels
  theme(plot.title = element_text(hjust = 0, size = 14, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.position = "none")

#Save with ggplot in 600 dpi and with a white background
ggsave("C:/Users/daank/OneDrive - University of Twente/Documents/Github/island_biology/plots/competitive_strength.png", width = 6, height = 3.5, dpi = 600, bg = "white")
