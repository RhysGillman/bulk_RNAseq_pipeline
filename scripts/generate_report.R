#!/usr/bin/env Rscript
suppressPackageStartupMessages (library(optparse, quietly = T))
suppressPackageStartupMessages (library(rmarkdown, quietly = T))

opt_list <- list(
  make_option(c("-d","--output_dir"),    type="character", default=NULL),
  make_option(c("-m","--metadata_file"), type="character", default=NULL)
)
opt <- parse_args(OptionParser(option_list = opt_list))
output_dir <- opt$output_dir
#pandoc dependency

pandoc_dir <- normalizePath(file.path(getwd(),"bin","pandoc-3.7.0.2","bin"))
Sys.setenv(RSTUDIO_PANDOC = pandoc_dir,
           PATH          = paste(pandoc_dir, Sys.getenv("PATH"), sep=":"))

report_dir <- file.path(output_dir,"report")

dir.create(report_dir, showWarnings = F)

if (dir.exists(file.path(output_dir, "qc"))) {
  # copy needed dirs to report dir
  invisible(file.copy(
    from      = file.path(output_dir, "qc"),
    to        = report_dir,
    recursive = TRUE
  ))
}
if (dir.exists(file.path(output_dir, "plots"))) {
  invisible(file.copy(
    from      = file.path(output_dir, "plots"),
    to        = report_dir,
    recursive = TRUE
  ))
}
message(paste0("Updating report in ",report_dir))

render(
  input       = "scripts/report_template.Rmd",
  output_file = file.path(report_dir, "Report.html"),
  knit_root_dir = output_dir,
  params      = list(
    output_dir    = output_dir,
    metadata_file = opt$metadata_file
  ),
  quiet=T
)
