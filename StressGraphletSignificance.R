library(tidyverse)
setwd("Desktop/GitHub_Repos/motifs/data/oxidative_stress/")

# Significance results
txid224308 <- read_tsv('txid224308/graphlet_significance.csv')
txid7955 <- read_tsv('txid7955/graphlet_significance.csv')
txid6239 <- read_tsv('txid6239/graphlet_significance.csv')
txid7227 <- read_tsv('txid7227/graphlet_significance.csv')
txid559292 <- read_tsv('txid559292/graphlet_significance.csv')

# results for a two-tailed test
txid224308_twotail <- txid224308 %>% filter(Significance >= 0.95 | Significance <= 0.05)
txid7955_twotail <- txid7955 %>% filter(Significance >= 0.95 | Significance <= 0.05)
txid6239_twotail <- txid6239 %>% filter(Significance >= 0.95 | Significance <= 0.05)
txid7227_twotail <- txid7227 %>% filter(Significance >= 0.95 | Significance <= 0.05)
txid559292_twotail <- txid559292 %>% filter(Significance >= 0.95 | Significance <= 0.05)

# results for a one-tailed test
txid224308_onetail <- txid224308 %>% filter(Significance <= 0.05)
txid7955_onetail <- txid7955 %>% filter(Significance <= 0.05)
txid6239_onetail <- txid6239 %>% filter(Significance <= 0.05)
txid7227_onetail <- txid7227 %>% filter(Significance <= 0.05)
txid559292_onetail <- txid559292 %>% filter(Significance <= 0.05)

# results for a one-tailed test at 0.01 significance
txid224308_onetail_01 <- txid224308 %>% filter(Significance <= 0.01)
txid7955_onetail_01 <- txid7955 %>% filter(Significance <= 0.01)
txid6239_onetail_01 <- txid6239 %>% filter(Significance <= 0.01)
txid7227_onetail_01 <- txid7227 %>% filter(Significance <= 0.01)
txid559292_onetail_01 <- txid559292 %>% filter(Significance <= 0.01)

# Get graphlets common to all 5 species
matrix1 <- rbind(txid224308, txid7955) %>% 
  rbind(txid6239) %>% rbind(txid7227) %>% rbind(txid559292)

all_5_sp <- matrix1 %>% count(Graphlet) %>% filter(n > 4)

# Generate a matrix to count occurences of significance graphlets across species
sig_matrix1 <- rbind(txid224308_onetail, txid7955_onetail) %>% 
  rbind(txid6239_onetail) %>% rbind(txid7227_onetail) %>% rbind(txid559292_onetail)

sig_matrix2 <- left_join(txid224308_onetail, txid7955_onetail, by = 'Graphlet') %>%
  left_join(txid6239_onetail, by = 'Graphlet') %>% 
  left_join(txid7227_onetail, by = 'Graphlet') %>%
  left_join(txid559292_onetail, by = 'Graphlet')

sig_3_plus_sp <- sig_matrix1 %>% count(Graphlet) %>% filter(n > 2)

# Get the number of unique graphlets found to be significant
length(unique(c(txid224308_onetail$Graphlet, txid7955_onetail$Graphlet, txid6239_onetail$Graphlet, txid7227_onetail$Graphlet, txid559292_onetail$Graphlet)))

# Print out the species that the common graphlets are found in
for (g in sig_3_plus_sp$Graphlet) {
  if (g %in% txid224308_onetail$Graphlet) {
    print(paste(g, "found in txid224308"))
  } else {print("not found")}
}

for (g in sig_3_plus_sp$Graphlet) {
  if (g %in% txid7955_onetail$Graphlet) {
    print(paste(g, "found in txid7955"))
  } else {print("not found")}
}

for (g in sig_3_plus_sp$Graphlet) {
  if (g %in% txid6239_onetail$Graphlet) {
    print(paste(g, "found in txid6239"))
  } else {print("not found")}
}

for (g in sig_3_plus_sp$Graphlet) {
  if (g %in% txid7227_onetail$Graphlet) {
    print(paste(g, "found in txid7227"))
  } else {print("not found")}
}

for (g in sig_3_plus_sp$Graphlet) {
  if (g %in% txid559292_onetail$Graphlet) {
    print(paste(g, "found in txid559292"))
  } else {print("not found")}
}

# Get the 3 most frequent significant graphlets for each species
txid224308_onetail %>% sort_by(desc(txid224308_onetail$Count)) %>% head(3)
txid6239_onetail %>% sort_by(desc(txid6239_onetail$Count)) %>% head(3)
txid7227_onetail %>% sort_by(desc(txid7227_onetail$Count)) %>% head(3)
txid7955_onetail %>% sort_by(desc(txid7955_onetail$Count)) %>% head(3)
txid559292_onetail %>% sort_by(desc(txid559292_onetail$Count)) %>% head(3)

# Get most count of unique graphlest in most frequent
top_3_most_common <- bind_rows(txid224308_onetail %>% sort_by(desc(txid224308_onetail$Count)) %>% head(3),
          txid6239_onetail %>% sort_by(desc(txid6239_onetail$Count)) %>% head(3),
          txid7227_onetail %>% sort_by(desc(txid7227_onetail$Count)) %>% head(3),
          txid7955_onetail %>% sort_by(desc(txid7955_onetail$Count)) %>% head(3),
          txid559292_onetail %>% sort_by(desc(txid559292_onetail$Count)) %>% head(3)
)

# Create key to identify graphlets
keyList <- c()

for (i in 1:98) {
  keyList <- c(keyList, paste0("G", i))
}

graphletList <- c(
  "((0, 1), (0, 1), (1, 1))",
  "((0, 1), (0, 4), (1, 5))",
  "((0, 1), (0, 5), (1, 4))",
  "((0, 4), (0, 4), (5, 5))",
  "((0, 5), (0, 5), (4, 4))",
  "((0, 4), (0, 5), (4, 5))",
  "((0, 1), (0, 7), (1, 7))",
  "((0, 5), (0, 7), (4, 7))",
  "((0, 4), (0, 7), (5, 7))",
  "((0, 7), (0, 7), (7, 7))",
  "((0, 1), (0, 3), (1, 2))",
  "((0, 1), (0, 2), (1, 3))",
  "((0, 1), (0, 6), (1, 6))",
  "((0, 2), (0, 4), (3, 5))",
  "((0, 2), (0, 5), (3, 4))",
  "((0, 3), (0, 5), (2, 4))",
  "((0, 3), (0, 4), (2, 5))",
  "((0, 4), (0, 6), (5, 6))",
  "((0, 5), (0, 6), (4, 6))",
  "((0, 3), (0, 7), (2, 7))",
  "((0, 2), (0, 7), (3, 7))",
  "((0, 6), (0, 7), (6, 7))",
  "((0, 2), (0, 2), (3, 3))",
  "((0, 3), (0, 3), (2, 2))",
  "((0, 2), (0, 3), (2, 3))",
  "((0, 3), (0, 6), (2, 6))",
  "((0, 2), (0, 6), (3, 6))",
  "((0, 6), (0, 6), (6, 6))",
  "((1, 1), (1, 1), (1, 1))",
  "((1, 1), (1, 4), (1, 5))",
  "((1, 4), (1, 4), (5, 5))",
  "((1, 5), (1, 5), (4, 4))",
  "((1, 4), (1, 5), (4, 5))",
  "((4, 5), (4, 5), (4, 5))",
  "((4, 4), (4, 5), (5, 5))",
  "((1, 1), (1, 7), (1, 7))",
  "((1, 5), (1, 7), (4, 7))",
  "((1, 4), (1, 7), (5, 7))",
  "((1, 7), (1, 7), (7, 7))",
  "((4, 4), (5, 7), (5, 7))",
  "((4, 7), (4, 7), (5, 5))",
  "((4, 5), (4, 7), (5, 7))",
  "((4, 7), (5, 7), (7, 7))",
  "((7, 7), (7, 7), (7, 7))",
  "((1, 1), (1, 2), (1, 3))",
  "((1, 2), (1, 4), (3, 5))",
  "((1, 2), (1, 5), (3, 4))",
  "((1, 3), (1, 4), (2, 5))",
  "((1, 3), (1, 5), (2, 4))",
  "((2, 5), (3, 4), (4, 5))",
  "((2, 5), (3, 5), (4, 4))",
  "((2, 4), (3, 4), (5, 5))",
  "((2, 4), (3, 5), (4, 5))",
  "((1, 1), (1, 6), (1, 6))",
  "((1, 4), (1, 6), (5, 6))",
  "((1, 5), (1, 6), (4, 6))",
  "((4, 5), (4, 6), (5, 6))",
  "((4, 4), (5, 6), (5, 6))",
  "((4, 6), (4, 6), (5, 5))",
  "((1, 6), (1, 7), (6, 7))",
  "((4, 6), (5, 7), (6, 7))",
  "((4, 7), (5, 6), (6, 7))",
  "((6, 7), (6, 7), (7, 7))",
  "((1, 3), (1, 7), (2, 7))",
  "((1, 2), (1, 7), (3, 7))",
  "((2, 5), (3, 7), (4, 7))",
  "((2, 7), (3, 4), (5, 7))",
  "((2, 7), (3, 5), (4, 7))",
  "((2, 4), (3, 7), (5, 7))",
  "((2, 7), (3, 7), (7, 7))",
  "((1, 2), (1, 2), (3, 3))",
  "((1, 2), (1, 3), (2, 3))",
  "((1, 3), (1, 3), (2, 2))",
  "((1, 2), (1, 6), (3, 6))",
  "((1, 3), (1, 6), (2, 6))",
  "((1, 6), (1, 6), (6, 6))",
  "((2, 4), (2, 5), (3, 3))",
  "((2, 2), (3, 4), (3, 5))",
  "((2, 3), (2, 4), (3, 5))",
  "((2, 3), (2, 5), (3, 4))",
  "((2, 5), (3, 6), (4, 6))",
  "((2, 6), (3, 5), (4, 6))",
  "((2, 4), (3, 6), (5, 6))",
  "((2, 6), (3, 4), (5, 6))",
  "((4, 6), (5, 6), (6, 6))",
  "((2, 7), (2, 7), (3, 3))",
  "((2, 3), (2, 7), (3, 7))",
  "((2, 2), (3, 7), (3, 7))",
  "((2, 7), (3, 6), (6, 7))",
  "((2, 6), (3, 7), (6, 7))",
  "((6, 6), (6, 7), (6, 7))",
  "((2, 3), (2, 3), (2, 3))",
  "((2, 2), (2, 3), (3, 3))",
  "((2, 3), (2, 6), (3, 6))",
  "((2, 2), (3, 6), (3, 6))",
  "((2, 6), (2, 6), (3, 3))",
  "((2, 6), (3, 6), (6, 6))",
  "((6, 6), (6, 6), (6, 6))"
)

graphlet_key <- tibble(
  Graphlet = graphletList, Key = keyList
)

# Join key back to datasets
txid224308_onetail <- inner_join(txid224308_onetail, graphlet_key)
txid6239_onetail <- inner_join(txid6239_onetail, graphlet_key)
txid559292_onetail <- inner_join(txid559292_onetail, graphlet_key)
txid7227_onetail <- inner_join(txid7227_onetail, graphlet_key)
txid7955_onetail <- inner_join(txid7955_onetail, graphlet_key)