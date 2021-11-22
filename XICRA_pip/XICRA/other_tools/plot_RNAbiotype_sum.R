#!/usr/bin/env Rscript
library("optparse")
library("XICRA.stats") ## Get this from conda: 

## get options
option_list = list(
  make_option(c("-f", "--file"), type="character", help="biotypes file", metavar="character"),
  make_option(c("-o", "--output"), type="character", help="name", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

## get message
if (is.null(opt$file)){
  print_help(opt_parser)
  stop("No arguments provided", call.=FALSE)
}

## get data
data_biotypes <- XICRA.stats::get_data(opt$file)

## plot into biotypes-plot.pdf in provided folder
XICRA.stats::plot_biotypes(data_biotypes, out_folder = opt$output)