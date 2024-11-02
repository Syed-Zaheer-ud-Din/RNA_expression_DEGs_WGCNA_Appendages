
library(ggplot2)
library(tidyr)
library(readxl)
library(ggpubr)

setwd("C:/Users/szd_z/OneDrive/Desktop/desktop_22_07_2023/Thesis")
# Load data from file
data1 <- read_excel("coldll.xlsx")

# Gather the Tibia, Femur, Mandibular Horn Length, and Mandibular Horn Width into a single column
data_long <- gather(data1, Measurement, Values, -Group)

# Remove rows with missing or non-finite values
data_long <- na.omit(data_long)
data_long <- data_long[is.finite(data_long$Values), ]

# Create the boxplot with color and significance stars
p <- ggplot(data_long, aes(x = Group, y = Values, fill = Measurement)) +
  geom_boxplot() +
  ggtitle("Dll_RNAi Knockdown") +
  xlab("Groups") +
  ylab("Measurement Value") +
  facet_wrap(~ Measurement, ncol = 4) +  # Set ncol to 4 for a single row with 4 columns
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Perform t-test and add significance stars

p <- p + stat_compare_means(comparisons = list(c("Group1", "Group2"), c("Group1", "Group3"), c("Group2", "Group3")), label = "p.signif", method = "t.test")

# Print the plot
print(p)



##################3
