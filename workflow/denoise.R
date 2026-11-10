#!/usr/bin/env Rscript
# DADA2 forward-read workflow for ONE community (one FASTQ).
# Usage: Rscript denoise.R fastq outdir community silva truncLen trimLeft errorModel threads minReads
# Strategy (from the original denoise_retry.R):
#   attempt 1: truncLen, maxEE 2, truncQ 2, minLen 20, chosen error model
#   attempt 2: no truncation, maxEE 5, truncQ 4, minLen 50, noqualErrfun
suppressPackageStartupMessages(library(dada2))
a <- commandArgs(TRUE)
if (length(a) != 9) stop('Usage: denoise.R fastq outdir community silva truncLen trimLeft errorModel threads minReads')
fastq <- a[1]; out <- a[2]; name <- a[3]; silva <- a[4]
trunc_len <- as.integer(a[5]); trim_left <- as.integer(a[6]); model <- a[7]
threads <- as.integer(a[8]); min_reads <- as.integer(a[9])
mt <- if (threads > 1) threads else FALSE
dir.create(out, recursive = TRUE, showWarnings = FALSE)
unlink(file.path(out, c('error.txt', 'error_first.txt', 'status.txt', 'asvs.csv', 'counts.csv', 'qc.csv')))

run_one <- function(truncLen, maxEE, truncQ, minLen, errname, label) {
  set.seed(1729)
  errfun <- if (errname == 'noqual') noqualErrfun else loessErrfun
  filt <- file.path(out, 'forward.filtered.fastq.gz')
  qc <- filterAndTrim(fastq, filt, truncLen = truncLen, trimLeft = trim_left, maxN = 0,
                      maxEE = maxEE, truncQ = truncQ, minLen = minLen, rm.phix = TRUE,
                      multithread = FALSE, compress = TRUE)
  if (qc[1, 2] < min_reads) stop(sprintf('Fewer than %d reads after filtering', min_reads))
  err <- learnErrors(filt, multithread = mt, randomize = FALSE, nbases = 1e8,
                     verbose = FALSE, errorEstimationFunction = errfun)
  if (any(!is.finite(getErrors(err)))) stop('Error learning failed')
  d <- dada(filt, err = err, multithread = mt)
  st <- makeSequenceTable(d); rownames(st) <- name
  nc <- removeBimeraDenovo(st, method = 'consensus', multithread = mt)
  if (!ncol(nc) || sum(nc) == 0) stop('No nonchimeric ASVs')
  tx <- assignTaxonomy(nc, silva, multithread = mt, tryRC = TRUE, verbose = FALSE)
  saveRDS(nc, file.path(out, 'seqtab_nochim.rds')); saveRDS(tx, file.path(out, 'taxa_raw.rds'))
  ids <- sprintf('%s_ASV%04d', name, seq_len(ncol(nc)))
  rank <- function(r) if (r %in% colnames(tx)) tx[, r] else rep(NA_character_, nrow(tx))
  genus <- rank('Genus'); kingdom <- tx[, 1]
  eligible <- !is.na(genus) & !is.na(kingdom) & kingdom == 'Bacteria'
  write.csv(data.frame(ASV = ids, Sequence = colnames(nc), Kingdom = kingdom,
                       Phylum = rank('Phylum'), Class = rank('Class'), Order = rank('Order'),
                       Family = rank('Family'), Genus = genus,
                       Eligible = ifelse(eligible, 'TRUE', 'FALSE'),
                       Reason = ifelse(eligible, '', 'unclassified_or_nonbacterial')),
            file.path(out, 'asvs.csv'), row.names = FALSE, na = '')
  write.csv(data.frame(sample = name, ASV = ids, count = as.integer(nc[1, ])),
            file.path(out, 'counts.csv'), row.names = FALSE)
  write.csv(data.frame(sample = name, input = qc[1, 1], filtered = qc[1, 2], denoised = sum(st),
                       nonchim = sum(nc), asvs = ncol(nc), error_model = errname, attempt = label),
            file.path(out, 'qc.csv'), row.names = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(out, 'R_session.txt'))
  writeLines(label, file.path(out, 'status.txt'))
  TRUE
}

tryCatch({
  png(file.path(out, 'quality.png'), width = 1200, height = 700)
  print(plotQualityProfile(fastq)); dev.off()
}, error = function(e) try(dev.off(), silent = TRUE))

ok <- FALSE
tryCatch({ run_one(trunc_len, 2, 2, 20, model, 'complete'); ok <- TRUE },
         error = function(e) writeLines(conditionMessage(e), file.path(out, 'error_first.txt')))
if (!ok) tryCatch({ run_one(0, 5, 4, 50, 'noqual', 'complete_relaxed'); ok <- TRUE },
                  error = function(e) {
                    writeLines(conditionMessage(e), file.path(out, 'error.txt'))
                    writeLines('failed', file.path(out, 'status.txt'))
                  })
message(name, ' ok=', ok)
quit(status = if (ok) 0 else 3)
