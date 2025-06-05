library(tidyverse)
library(ggsci)
setwd('Desktop/GitHub_Repos/motifs/')
df <- read_tsv("iterations_swaps.txt")
df <- df %>% mutate(Randomization = Swaps/Edges)
sp_labels <- c("B. subtilis", "S. cerevisiae", "C. elegans", "D. melanogaster", "D. rerio")

highlight_points <- tibble(
  Iterations = c(500, 500, 5000, 10000, 25000),
  Swaps = c(31, 47, 582, 1152, 2984),
  Randomization = c(31/68, 47/106, 582/1653, 1152/3092, 2984/8138),
  Species = c("txid224308", "txid7955", "txid7227", "txid6239", "txid559292")  # Replace with actual species names
)

ggplot()+
  geom_line(data=df, aes(Iterations, Randomization, color=Species))+
  geom_point(data=highlight_points, aes(x=Iterations, y=Randomization, color=Species))+
  theme_minimal()+
  scale_color_jama(labels=sp_labels)+
  scale_y_continuous(labels=scales::percent)+
  scale_x_continuous(label=scales::comma)+
  # scale_x_log10(label=scales::comma)+
  labs(title="Percent of edges changed in randomized oxidative \nstress response networks with threshold marked")

ggplot()+
  geom_line(data=df, aes(Iterations, Swaps, color=Species))+
  geom_point(data=highlight_points, aes(x=Iterations, y=Swaps, color=Species))+
  theme_minimal()+
  scale_color_jama(labels=sp_labels)+
  scale_x_continuous(label=scales::comma)+
  scale_y_log10(label=scales::comma)+
  labs(title="Number of edges changed in randomized oxidative \nstress response networks with threshold marked")
  
