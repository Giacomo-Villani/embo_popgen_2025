library(tidyverse)

setwd("D:/git/embo_popgen_2025/Giacomo_Villani")

Salamandra <- read.table("Salamander.txt", header=T)

for (j in 3:5) {
for (i in 1:8){
if (Salamandra[i, j] == Salamandra[9, j]) {
Salamandra[i, j] <- 0
}
else {
Salamandra[i, j] <- 1
}
}
}

Salamandra <- Salamandra[1:8,] %>%
  mutate(across(Locus1:Locus3, as.numeric))

Pop1 <- Salamandra %>% filter(FAMID == "Pop1") 
Pop2 <- Salamandra %>% filter(FAMID == "Pop2")

freq_pop1 <- colSums(Pop1[, 3:5]) / nrow(Pop1)
freq_pop2 <- colSums(Pop2[, 3:5]) / nrow(Pop2)

## to finish with the other parts of the exercise