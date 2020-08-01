
# installing and loading packages
packages <- c("readstata13", "peperr", "QuantPsyc", "lme4", "bain", "tidyverse")

installed <- installed.packages()

sapply(packages, function(x) {
  if(!x %in% installed) install.packages(x, dep = TRUE) 
  library(x, character.only = TRUE)
  }
)

