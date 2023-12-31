#' Enrich for input motifs in a set of sequences.
#'
#' Given a set of target and background sequences, test if the input motifs
#' are significantly enriched in the targets sequences relative to the
#' background sequences. See the "Sequence manipulation and scanning" vignette.
#'
#' @param motifs See [convert_motifs()] for acceptable motif formats.
#' @param bkg.sequences \code{\link{XStringSet}} Optional. If missing,
#'    [shuffle_sequences()] is used to create background sequences from
#'    the input sequences.
#' @param max.p `numeric(1)` P-value threshold.
#' @param max.q `numeric(1)` Adjusted P-value threshold. This is only useful
#'    if multiple motifs are being enriched for.
#' @param max.e `numeric(1)`. The E-value is calculated by multiplying the
#'    P-value with the number of input motifs times two
#'    (McLeay and Bailey 2010).
#' @param qval.method `character(1)` See [stats::p.adjust()].
#' @param verbose `numeric(1)` 0 for no output, 4 for max verbosity.
#' @param shuffle.k `numeric(1)` The k-let size to use when shuffling input
#'    sequences. Only used if no background sequences are input. See
#'    [shuffle_sequences()].
#' @param shuffle.method `character(1)` One of `c('euler', 'markov', 'linear')`.
#'    See [shuffle_sequences()].
#' @param return.scan.results `logical(1)` Return output from
#'    [scan_sequences()]. For large jobs, leaving this as
#'    `FALSE` can save a small amount time by preventing construction of the complete
#'    results `data.frame` from [scan_sequences()].
#' @param nthreads `numeric(1)` Run [scan_sequences()] in parallel with `nthreads`
#'    threads. `nthreads = 0` uses all available threads.
#'    Note that no speed up will occur for jobs with only a single motif and
#'    sequence.
#' @param rng.seed `numeric(1)` Set random number generator seed. Since shuffling
#'    can occur simultaneously in multiple threads using C++, it cannot communicate
#'    with the regular `R` random number generator state and thus requires an
#'    independent seed. Each individual sequence in an \code{\link{XStringSet}} object
#'    will be given the following seed: `rng.seed * index`. See [shuffle_sequences()].
#' @param motif_pvalue.k `numeric(1)` Control [motif_pvalue()] approximation.
#'    See [motif_pvalue()].
#' @param use.gaps `logical(1)` Set this to `FALSE` to ignore motif gaps, if
#'    present.
#' @param allow.nonfinite `logical(1)` If `FALSE`, then apply a pseudocount if
#'    non-finite values are found in the PWM. Note that if the motif has a
#'    pseudocount greater than zero and the motif is not currently of type PWM,
#'    then this parameter has no effect as the pseudocount will be
#'    applied automatically when the motif is converted to a PWM internally. This
#'    value is set to `FALSE` by default in order to stay consistent with
#'    pre-version 1.8.0 behaviour. A message will be printed if a pseudocount
#'    is applied. To disable this, set `options(pseudocount.warning=FALSE)`.
#' @param warn.NA `logical(1)` Whether to warn about the presence of non-standard
#'    letters in the input sequence, such as those in masked sequences.
#' @param scan_sequences.qvals.method `character(1)` One of
#'    `c("fdr", "BH", "bonferroni")`. The method for calculating adjusted
#'    P-values for individual motif hits. These are described in depth in the
#'    Sequence Searches vignette.
#' @param mode `character(1)` One of `c("total.hits", "seq.hits")`. The former
#'    enriches for the total count of motif hits across all sequences,
#'    whereas the latter only counts motif hits once per sequence (useful for
#'    cases where there are many small sequences).
#' @param pseudocount `integer(1)` Add a pseudocount to the motif hit counts
#'    when performing the Fisher test.
#'
#' @return `DataFrame` Enrichment results in a `DataFrame`. Function args and
#'    (optionally) scan results are stored in the `metadata` slot.
#'
#' @details
#' To find enriched motifs, [scan_sequences()] is run on both
#' target and background sequences.
#' [stats::fisher.test()] is run to test for enrichment.
#'
#' See [scan_sequences()] for more info on scanning parameters.
#'
#' @examples
#' data(ArabidopsisPromoters)
#' data(ArabidopsisMotif)
#' if (R.Version()$arch != "i386") {
#' enrich_motifs(ArabidopsisMotif, ArabidopsisPromoters, threshold = 0.01)
#' }
#'
#' @references
#'
#' McLeay R, Bailey TL (2010). “Motif Enrichment Analysis: A unified
#' framework and method evaluation.” *BMC Bioinformatics*, **11**.
#'
#' @author Benjamin Jean-Marie Tremblay \email{benjamin.tremblay@@uwaterloo.ca}
#' @seealso [scan_sequences()], [shuffle_sequences()],
#'    [add_multifreq()], [motif_pvalue()]
#' @inheritParams scan_sequences
#' @export
enrich_motifs <- function(motifs, sequences, bkg.sequences,
  max.p = 10e-4, max.q = 10e-4, max.e = 10e-4, qval.method = "fdr",
  threshold = 0.0001, threshold.type = "pvalue", verbose = 0, RC = TRUE,
  use.freq = 1, shuffle.k = 2, shuffle.method = "euler",
  return.scan.results = FALSE, nthreads = 1, rng.seed = sample.int(1e4, 1),
  motif_pvalue.k = 8, use.gaps = TRUE, allow.nonfinite = FALSE,
  warn.NA = TRUE, no.overlaps = TRUE, no.overlaps.by.strand = FALSE,
  no.overlaps.strat = "score", respect.strand = FALSE,
  motif_pvalue.method = c("dynamic", "exhaustive"),
  scan_sequences.qvals.method = c("BH", "fdr", "bonferroni"),
  mode = c("total.hits", "seq.hits"), pseudocount = 1) {

  motif_pvalue.method <- match.arg(motif_pvalue.method)
  scan_sequences.qvals.method <- match.arg(scan_sequences.qvals.method)
  mode <- match.arg(mode)

  # param check --------------------------------------------
  args <- as.list(environment())
  all_checks <- character(0)
  if (!qval.method %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
                          "fdr", "none")) {
    qval.method_check <- paste0(" * Incorrect 'qval.method': expected `holm`, ",
                                "`hochberg`, `hommel`, `bonferroni`, `BH`, `BY`, ",
                                "`fdr` or `none`; got `",
                                qval.method, "`")
    qval.method_check <- wmsg2(qval.method_check, 4, 2)
    all_checks <- c(all_checks, qval.method_check)
  }
  if (!threshold.type %in% c("logodds", "pvalue", "logodds.abs", "qvalue")) {
    threshold.type_check <- paste0(" * Incorrect 'threshold.type': expected ",
                                   "`logodds`, `logodds.abs`, `pvalue` or `qvalue`; got `",
                                   threshold.type, "`")
    threshold.type_check <- wmsg2(threshold.type_check, 4, 2)
    all_checks <- c(all_checks, threshold.type_check)
  }
  if (!shuffle.method %in% c("euler", "markov", "linear", "random")) {
    shuffle.method_check <- paste0(" * Incorrect 'shuffle.method': expected ",
                                   "`markov`, `linear` or `random`; got `",
                                   shuffle.method, "`")
    shuffle.method_check <- wmsg2(shuffle.method_check, 4, 2)
    all_checks <- c(all_checks, shuffle.method_check)
  }
  char_check <- check_fun_params(list(qval.method = args$qval.method,
                                      threshold.type = args$threshold.type,
                                      shuffle.method = args$shuffle.method),
                                 numeric(), logical(), TYPE_CHAR)
  num_check <- check_fun_params(list(max.p = args$max.p, max.q = args$max.q,
                                     max.e = args$max.e,
                                     threshold = args$threshold,
                                     verbose = args$verbose, use.freq = args$use.freq,
                                     shuffle.k = args$shuffle.k,
                                     nthreads = args$nthreads,
                                     motif_pvalue.k = args$motif_pvalue.k),
                                c(1, 1, 1, 0, 1, 1, 1, 1, 1), logical(), TYPE_NUM)
  logi_check <- check_fun_params(list(RC = args$RC, use.gaps = args$use.gaps,
                                      return.scan.results = args$return.scan.results),
                                 numeric(), logical(), TYPE_LOGI)
  s4_check <- check_fun_params(list(sequences = args$sequences,
                                    bkg.sequences = args$bkg.sequences),
                               numeric(), c(FALSE, TRUE), TYPE_S4)
  all_checks <- c(all_checks, char_check, num_check, logi_check, s4_check)
  if (length(all_checks) > 0) stop(all_checks_collapse(print(all_checks)))
  pseudocount <- as.integer(pseudocount)[1]
  if (is.na(pseudocount))
    stop(" * Incorrect 'pseudocount': got `NA`", call. = FALSE)
  #---------------------------------------------------------

  if (verbose > 2) {
    message(" > Input parameters")
    message("   > motifs:              ", deparse(substitute(motifs)))
    message("   > sequences:           ", deparse(substitute(sequences)))
    message("   > bkg.sequences:       ", ifelse(!missing(bkg.sequences),
                                                 deparse(substitute(bkg.sequences)),
                                                 "none"))
    if (missing(bkg.sequences)) {
      message("   > shuffle.k:           ", shuffle.k)
      message("   > shuffle.method:      ", shuffle.method)
    }
    message("   > max.p:               ", max.p)
    message("   > max.q:               ", max.q)
    message("   > max.e:               ", max.e)
    message("   > qval.method:         ", qval.method)
    message("   > threshold:           ", threshold)
    message("   > threshold.type:      ", threshold.type)
    message("   > verbose:             ", verbose)
    message("   > RC:                  ", RC)
    message("   > respect.strand:      ", respect.strand)
    message("   > use.freq:            ", use.freq)
    message("   > use.gaps:            ", use.gaps)
    message("   > no.overlaps:         ", no.overlaps)
    message("   > return.scan.results: ", return.scan.results)
  }

  if (mode == "seq.hits" && no.overlaps) {
    if (verbose > 2) {
      message(" ! ignoring `no.overlaps=TRUE`, not needed when `mode=\"seq.hits\"`")
    }
    no.overlaps <- FALSE
  }

  motifs <- convert_motifs(motifs)
  motifs <- convert_type_internal(motifs, "PWM")

  if (missing(bkg.sequences)) {
    if (verbose > 0) message(" > Shuffling input sequences")
    bkg.sequences <- shuffle_sequences(sequences, shuffle.k, shuffle.method,
                                       nthreads = nthreads, rng.seed = rng.seed)
  }

  if (!is.list(motifs)) motifs <- list(motifs)
  motcount <- length(motifs)

  if (threshold.type == "pvalue") {
    if (verbose > 0)
      message(" > Converting P-values to logodds thresholds")
    threshold <- suppressMessages(motif_pvalue(motifs, pvalue = threshold,
        method = motif_pvalue.method,
        use.freq = use.freq, k = motif_pvalue.k, allow.nonfinite = allow.nonfinite))
    max.scores <- vapply(motifs, function(x)
      suppressMessages(motif_score(x, 1, use.freq, threshold.type = "fromzero",
          allow.nonfinite = allow.nonfinite)),
      numeric(1))
    if (verbose > 3) {
      mot.names <- vapply(motifs, function(x) x@name, character(1))
      for (i in seq_along(threshold)) {
        message("   > Motif ", mot.names[i], ": max.score = ", max.scores[i],
                ", threshold.score = ", threshold[i])
      }
    }
    threshold.type <- "logodds.abs"
  }

  res.all <- enrich_mots2(motifs, sequences, bkg.sequences, threshold,
    verbose, RC, use.freq, threshold.type, motcount, return.scan.results,
    nthreads, args[-(1:3)], use.gaps, allow.nonfinite, warn.NA,
    no.overlaps, no.overlaps.by.strand, no.overlaps.strat, respect.strand,
    motif_pvalue.method, scan_sequences.qvals.method, mode, pseudocount)

  if (nrow(res.all) == 0) {
    message(" ! No enriched motifs")
    return(res.all)
  }

  res.all$Qval <- p.adjust(res.all$Pval, method = qval.method)
  res.all$Eval <- res.all$Pval * length(motifs) * 2

  res.all <- res.all[!is.na(res.all$motif), ]

  res.all <- res.all[order(res.all$Qval), ]
  res.all <- res.all[res.all$Pval < max.p, ]
  res.all <- res.all[res.all$Qval < max.q, ]
  res.all <- res.all[res.all$Eval < max.e, ]

  res.all$pct.target.seq.hits <- 100 * (res.all$target.seq.hits / res.all$target.seq.count)
  res.all$pct.bkg.seq.hits <- 100 * (res.all$bkg.seq.hits / res.all$bkg.seq.count)

  res.all$motif.consensus <- to_df(motifs)$consensus[res.all$motif.i]
  res.all$target.enrichment <- res.all$pct.target.seq.hits / res.all$pct.bkg.seq.hits

  if (nrow(res.all) < 1 || is.null(res.all)) {
    message(" ! No enriched motifs")
    return(res.all)
  }

  res.all

}

enrich_mots2 <- function(motifs, sequences, bkg.sequences, threshold,
  verbose, RC, use.freq, threshold.type, motcount, return.scan.results,
  nthreads, args, use.gaps, allow.nonfinite, warn.NA, no.overlaps,
  no.overlaps.by.strand, no.overlaps.strat, respect.strand,
  motif_pvalue.method, scan_sequences.qvals.method, mode, pseudocount) {

  seq.names <- names(sequences)
  if (is.null(seq.names)) seq.names <- seq_len(length(sequences))
  bkg.names <- names(bkg.sequences)
  if (is.null(seq.names)) bkg.names <- seq_len(length(bkg.sequences))

  seq.widths <- width(sequences)
  bkg.widths <- width(bkg.sequences)

  if (verbose > 0) message(" > Scanning input sequences")
  results <- scan_sequences(motifs, sequences, threshold, threshold.type,
    RC, use.freq, verbose = verbose - 1, nthreads = nthreads,
    use.gaps = use.gaps, allow.nonfinite = allow.nonfinite, warn.NA = warn.NA,
    no.overlaps = no.overlaps, no.overlaps.by.strand = no.overlaps.by.strand,
    no.overlaps.strat = no.overlaps.strat, respect.strand = respect.strand,
    motif_pvalue.method = motif_pvalue.method,
    calc.qvals.method = scan_sequences.qvals.method)

  if (verbose > 0) message(" > Scanning background sequences")
  results.bkg <- scan_sequences(motifs, bkg.sequences, threshold,
    threshold.type, RC, use.freq, verbose = verbose - 1, nthreads = nthreads,
    use.gaps = use.gaps, allow.nonfinite = allow.nonfinite, warn.NA = warn.NA,
    no.overlaps = no.overlaps, no.overlaps.by.strand = no.overlaps.by.strand,
    no.overlaps.strat = no.overlaps.strat, respect.strand = respect.strand,
    motif_pvalue.method = motif_pvalue.method,
    calc.qvals.method = scan_sequences.qvals.method)

  if (verbose > 0) message(" > Processing output")

  results2 <- split_by_motif_enrich(motifs, results)
  results.bkg2 <- split_by_motif_enrich(motifs, results.bkg)

  if (length(results2) == 0) {
    out <- DataFrame()
    if (return.scan.results) {
      out@metadata <- list(scan.target = results, scan.bkg = results.bkg,
                           args = args)
    } else {
      out@metadata <- list(args = args)
    }
    return(out)
  }

  if (length(results.bkg2) == 0) results.bkg2 <- DataFrame()

  if (verbose > 0) message(" > Testing motifs for enrichment")

  if (verbose > 3) tmp_pb <- FALSE
  results.all <- mapply(function(x, y, z)
                          enrich_mots2_subworker(x, y, z, seq.widths,
                                                  bkg.widths, sequences,
                                                  RC, bkg.sequences, verbose,
                                                  use.gaps, pseudocount, mode),
                        results2, results.bkg2, motifs,
                        SIMPLIFY = FALSE)

  results.all <- do.call(rbind, results.all)
  results.all$motif.i <- as.integer(seq_along(motifs))

  if (return.scan.results) {
    results.all@metadata <- list(scan.target = results, scan.bkg = results.bkg,
                                 args = args)
  } else {
    results.all@metadata <- list(args = args)
  }

  results.all

}

split_by_motif_enrich <- function(motifs, results) {

  mot.names <- vapply(motifs, function(x) x@name, character(1))
  results.sep <- lapply(mot.names, function(x) results[results$motif == x, ])

  results.sep

}

enrich_mots2_subworker <- function(results, results.bkg, motifs,
                                   seq.widths, bkg.widths, sequences, RC,
                                   bkg.sequences, verbose, use.gaps,
                                   pseudocount, mode) {

  seq.hits <- as.vector(results$start)
  bkg.hits <- as.vector(results.bkg$start)

  if (length(seq.hits) > 0 && length(bkg.hits) == 0 && pseudocount == 0) {
    warning(wmsg("Found hits for motif '", motifs@name,
                 "' in target sequences but none in bkg, ",
                 "significance will not be calculated and instead a ",
                 "P-value of 0 will be assigned. Set a pseudocount ",
                 "greater than 0 to calculate a P-value."), immediate. = TRUE,
            call. = FALSE)
    skip.calc <- TRUE
  } else {
    skip.calc <- FALSE
  }

  if (length(seq.hits) == 0 && length(bkg.hits) == 0)
    return(DataFrame(
        motif = motifs@name,
        motif.i = NA_integer_,
        motif.consensus = NA_character_,
        target.hits = 0L,
        target.seq.hits = 0L,
        target.seq.count = length(sequences),
        bkg.hits = 0L,
        bkg.seq.hits = 0L,
        bkg.seq.count = length(bkg.sequences),
        Pval = 1.0,
        Qval = NA_real_,
        Eval = NA_real_
    ))

  if (mode == "seq.hits") {
    seq.total <- length(seq.widths)
    bkg.total <- length(bkg.widths)
  } else {
    seq.total <- (mean(seq.widths) - ncol(motifs@motif) + 1) * length(sequences)
    if (RC) seq.total <- seq.total * 2
    bkg.total <- (mean(bkg.widths) - ncol(motifs@motif) + 1) * length(bkg.sequences)
    if (RC) bkg.total <- bkg.total * 2
  }

  if (use.gaps && motifs@gapinfo@isgapped) {
    x <- motifs@gapinfo@maxgap - motifs@gapinfo@mingap
    x <- prod(x + 1)
    if (mode == "total.hits") {
      seq.total <- seq.total * x
      bkg.total <- bkg.total * x
    }
  }

  seq.hits.n.n <- length(seq.hits)
  bkg.hits.n.n <- length(bkg.hits)

  if (mode == "total.hits") {
    seq.hits.n <- length(seq.hits)
    bkg.hits.n <- length(bkg.hits)
  } else {
    seq.hits.n <- sum(!duplicated(results$sequence.i))
    bkg.hits.n <- sum(!duplicated(results.bkg$sequence.i))
  }

  seq.no <- seq.total - seq.hits.n
  bkg.norm <- seq.total / bkg.total
  bkg.no <- bkg.total - bkg.hits.n

  results.table <- matrix(
    c(seq.hits.n, seq.no,
      as.integer(bkg.hits.n * bkg.norm), as.integer(bkg.no * bkg.norm)
     ) + pseudocount, nrow = 2, byrow = TRUE)

  # if (verbose > 3) message(" > Testing motif for enrichment: ", motifs@name)

  if (!skip.calc) {
    hits.p <- fisher.test(results.table, alternative = "greater")$p.value
  } else {
    hits.p <- 0
  }

  if (verbose > 3) {
    message("       ", motifs@name, " occurrences p-value: ", hits.p)
  }

  results <- DataFrame(
    motif = motifs@name,
    motif.i = NA_integer_,
    motif.consensus = NA_character_,
    target.hits = seq.hits.n.n,
    target.seq.hits = length(unique(as.vector(results$sequence.i))),
    target.seq.count = length(sequences),
    bkg.hits = bkg.hits.n.n,
    bkg.seq.hits = length(unique(as.vector(results.bkg$sequence.i))),
    bkg.seq.count = length(bkg.sequences),
    Pval = hits.p,
    Qval = NA_real_,
    Eval = NA_real_
  )

  results

}
