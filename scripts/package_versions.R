suppressPackageStartupMessages (library(tximport, quietly = T))
suppressPackageStartupMessages (library(consensusDE, quietly = T))
suppressPackageStartupMessages (library(sva, quietly = T))
suppressPackageStartupMessages (library(DOSE, quietly = T))

pkgs <- c("tximport", "consensusDE", "sva", "DOSE")
versions <- sapply(pkgs, function(pkg) paste(pkg, as.character(packageVersion(pkg))))
cat(versions, sep = "\n")
