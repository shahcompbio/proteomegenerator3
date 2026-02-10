#!/usr/bin/env Rscript
library(bambu)
library(optparse)
## load in commandline args
option_list <- list(
  make_option("--se", type="character", default=NULL,
              help="bambu read classes", metavar="character"),
  make_option("--merge", type="logical", default=FALSE,
              help="if multisample merge was performed", metavar="logical"),
  make_option("--min_cts", type="double", default=1,
              help="min full-length reads for detected transcript")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
# load summarized experiment object
load(opt$se)
# filter on full-length read counts
min_cts <- opt$min_cts
if (opt$merge) {
  fullLengthCounts <- assays(se)$fullLengthCounts
  rows_satisfying_conditions <- apply(fullLengthCounts, 1, function(row) any(row >= min_cts))
  # Subset the SummarizedExperiment object based on the condition
  se.detected <- se[rows_satisfying_conditions, ]
} else {
  se.detected <- se[assays(se)$fullLengthCounts >= min_cts,]
}
detectgtf <- rowRanges(se.detected)
writeToGTF(detectgtf, "detected_transcripts.gtf")
