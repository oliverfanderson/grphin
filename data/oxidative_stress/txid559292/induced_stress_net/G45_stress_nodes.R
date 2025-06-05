library(tidyverse)
orbits <- read_csv("Desktop/GitHub_Repos/motifs/data/oxidative_stress/txid559292/induced_stress_net/node_orbit.csv", col_names=FALSE)
ids <- read_csv("Desktop/GitHub_Repos/motifs/data/oxidative_stress/txid559292/induced_stress_net/protein_id_mapper.csv", col_names=FALSE)

colnames(ids)<-c("protein", "id")
?read_csv
orbits <- rowid_to_column(orbits, "id")

orbits <- orbits %>% filter(X114 > 0 | X115 > 0 | X116 > 0) %>% 
  inner_join(ids)

orbits$protein
