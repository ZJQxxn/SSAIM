install.packages('fANCOVA')
install.packages('limSolve')
install.packages('matlab')
install.packages('VGAM')
install.packages('MTS')
install.packages('glmnet')
install.packages('fastclime')
install.packages("scar")
install.packages("mboost")

install.packages("BiocManager")
BiocManager::install(version = "3.12") # For 4.0
BiocManager::install("GEOquery")


# 3.4
source("http://bioconductor.org/biocLite.R")
BiocManager::install(version = "3.6") # For 3.4
biocLite("GEOquery")