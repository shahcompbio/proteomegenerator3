#!/usr/bin/env Rscript
library(bambu)
library(optparse)
## load in commandline args
option_list = list(
  make_option("--se", type="character", default=NULL,
              help="bambu read classes", metavar="character"),
  make_option("--merge", type="logical", default=FALSE,
              help="if multisample merge was performed", metavar="logical"),
  make_option("--min_cpm", type="double", default=1,
              help="minimum cpm for detected transcript")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
# load summarized experiment object
load(opt$se)
# filter on cpms
min_cpm <- opt$min_cpm
if (opt$merge) {
  CPM <- assays(se)$CPM
  rows_satisfying_conditions <- apply(CPM, 1, function(row) any(row >= min_cpm))
  # Subset the SummarizedExperiment object based on the condition
  se.detected <- se[rows_satisfying_conditions, ]
} else {
  se.detected <- se[assays(se)$CPM >= min_cpm,]
}
detectgtf <- rowRanges(se.detected)
writeToGTF(detectgtf, "detected_transcripts.gtf")
