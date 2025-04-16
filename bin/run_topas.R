#!/usr/bin/env Rscript


source('TOPAS.R')

library("optparse")
option_list = list(
  make_option(c("-n", "--network"), type="character", default=NULL,
              help ="path to network file", metavar="character"),
  make_option(c("-s","--seeds"), type="character", default=NULL,
              help="path to seed file", metavar="character"),
  make_option(c("-o","--output"), type="character", default="out.txt",
              help="output file name [default=%default]", metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
if(is.null(opt$network)){
  print_help(opt_parser)
  stop("provide the networkfile", call.=FALSE)
}
if(is.null(opt$seeds)){
  print_help(opt_parser)
  stop("provide the seed file", call.=FALSE)
}
if(is.null(opt$output)){
  print_help(opt_parser)
  stop("provide the seed file", call.=FALSE)
}

## Load network
network_example = 
  readr::read_tsv(file = opt$network, 
                  col_names = c('V1','V2'), 
                  col_types = 'cc')
## Seed gene set
seeds_example = 
  readr::read_table(file = opt$seeds, 
                    col_names = 'V',
                    col_types = 'c')

module_example = 
  TOPAS(
    network = network_example,
    seeds = seeds_example$V,
    expansion_steps = 2,
    cores = 4
    )

readr::write_tsv(module_example, file = file.path(path.expand("./"), opt$output))
