#!/usr/bin/env Rscript
library("optparse")

## get options
option_list = list(
  make_option(c("-l", "--lib"),type="character",help="Library to check", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

library(opt$lib, character.only = TRUE)