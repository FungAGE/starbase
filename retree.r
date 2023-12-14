library(ggtree)
library(treeio)
library(tidyverse)

jplace_tree <- read.jplace("tmp/epa_result.jplace")
tree<-read.tree()

tip_data <- data.frame(
  tip_label = jplace_tree$tip.label,
  highlight = ifelse(jplace_tree$tip.label == "QUERY", TRUE,FALSE)
)
ggtree(jplace_tree) + #%<+% tip_data +
  geom_tippoint(aes(color = ifelse(highlight, "red", "black")), size = 3)
