install.packages(c("tidyverse", "readxl", "ggplot2", "pheatmap", "caret", "randomForest"))
library(readxl)
library(tidyverse)

# Cargar archivo
data <- read_excel("ASD meta abundance (1).xlsx")

# Ver primeras filas
head(data)

# Ver estructura
glimpse(data)

