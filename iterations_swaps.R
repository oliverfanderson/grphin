library(tidyverse)
library(ggsci)
setwd('Desktop/GitHub_Repos/motifs/')
df <- read_tsv("iterations_swaps.txt")
df <- df %>% mutate(Randomization = Swaps/Edges)
sp_labels <- c("B. subtilis", "S. cerevisiae", "C. elegans", "D. melanogaster", "D. rerio")
ggplot(df, aes(Iterations, Randomization, color=Species))+
  geom_line()+
  theme_minimal()+
  scale_color_jama(labels=sp_labels)+
  scale_y_continuous(labels=scales::percent)+
  scale_x_continuous(label=scales::comma)+
  # scale_x_log10(label=scales::comma)+
  labs(title="Swap iterations and percent of edges changed in \nrandomized oxidative stress response networks")

ggplot(df, aes(Iterations, Swaps, color=Species))+
  geom_line()+
  theme_minimal()+
  scale_color_jama(labels=sp_labels)+
  scale_x_continuous(label=scales::comma)+
  scale_y_log10(label=scales::comma)+
  labs(title="Swaps iterations and changed edges in randomized \noxidative stress response networks")
  
