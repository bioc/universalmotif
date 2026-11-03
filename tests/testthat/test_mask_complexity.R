context("mask_complexity()")

DUST_A <- paste0("ACGTACGTAAAAAAAAAAAAAAAAAAAAAAAAACGTTGCAGTCGATCGATGCATGCATCAGTCAGC",
  "TAGCTAGCATCGACTAGCTAGCATGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGTTGACCAGT")
DUST_B <- paste0("GGTTGTACTGTACAGGAACTTATAGAATAATTCTGGATGTATATATATATATATATATATATAT",
  "ATATATTCCAAACGAGCCTCCCGATAGTATCGAAATGGGGGGGGGGGGGATTTCCGTAGTTATACGTCTGCCAA")

## Regions as start-end pairs, counted from 1 with the end included.
dust_regions <- function(seq, ...) {
  found <- dust_intervals_cpp(seq, ...)[[1]]
  paste(found$start, found$end, sep = "-")
}

test_that("regions match dustmasker", {
  ## Reference regions from dustmasker 1.0.0 (BLAST 2.16.0), shifted from
  ## its 0-based positions.
  expect_equal(dust_regions(DUST_A, 64, 20, FALSE, 1), c("9-33", "51-127"))
  expect_equal(dust_regions(DUST_B, 64, 20, FALSE, 1), c("40-70", "101-113"))
  expect_equal(dust_regions(DUST_A, 32, 10, FALSE, 1), c("9-33", "64-127"))
  rna <- gsub("t", "u", tolower(DUST_B))
  expect_equal(dust_regions(rna, 64, 20, FALSE, 1), c("40-70", "101-113"))
})

test_that("other letters split sequences and lowercase can be masked", {
  split <- paste0(substr(DUST_B, 1, 55), "NNN", substr(DUST_B, 56, nchar(DUST_B)))
  expect_equal(dust_regions(split, 64, 20, FALSE, 1),
    c("40-55", "59-73", "104-116"))
  soft <- paste0("acgt", substr(DUST_B, 5, 120), "gattaca")
  expect_equal(dust_regions(soft, 64, 20, TRUE, 1),
    c("1-4", "40-70", "101-113", "121-127"))
  masked <- mask_complexity(Biostrings::BStringSet(c(s = soft)),
    lowercase = TRUE)
  expect_s4_class(masked, "BStringSet")
  expect_equal(substr(as.character(masked), 1, 6), c(s = "----GT"))
})

test_that("sequences keep their class and names", {
  seqs <- Biostrings::DNAStringSet(c(a = DUST_A, b = DUST_B, c = "ACGT"))
  masked <- mask_complexity(seqs, letter = "N", nthreads = 2)
  expect_s4_class(masked, "DNAStringSet")
  expect_equal(names(masked), c("a", "b", "c"))
  expect_equal(substr(as.character(masked[["b"]]), 37, 73),
    paste0("ATG", strrep("N", 31), "TCC"))
  expect_equal(as.character(masked[["c"]]), "ACGT")
  rna <- Biostrings::RNAStringSet(gsub("T", "U", DUST_B))
  expect_s4_class(mask_complexity(rna), "RNAStringSet")
  expect_equal(sum(Biostrings::letterFrequency(mask_complexity(seqs), "-")),
    25 + 77 + 31 + 13)
})

test_that("ranges work with mask_ranges()", {
  skip_if_not_installed("GenomicRanges")
  seqs <- Biostrings::DNAStringSet(c(a = DUST_A, b = DUST_B))
  r <- mask_complexity(seqs, return.ranges = TRUE)
  expect_s4_class(r, "GRanges")
  expect_equal(as.character(GenomicRanges::seqnames(r)), c("a", "a", "b", "b"))
  expect_equal(as.character(mask_ranges(seqs, r)),
    as.character(mask_complexity(seqs)))
  r <- mask_complexity(unname(seqs), return.ranges = TRUE)
  expect_equal(unique(as.character(GenomicRanges::seqnames(r))), c("1", "2"))
})

test_that("bad input is refused", {
  seqs <- Biostrings::DNAStringSet(DUST_A)
  expect_error(mask_complexity(DUST_A), "XStringSet")
  expect_error(mask_complexity(Biostrings::AAStringSet("MKV")), "DNA and RNA")
  expect_error(mask_complexity(seqs, window.size = 3), "above 3")
  expect_error(mask_complexity(seqs, threshold = 0), "above 0")
  expect_error(mask_complexity(seqs, threshold = 2.5), "whole number")
  expect_error(mask_complexity(seqs, letter = "--"), "single character")
  expect_error(mask_complexity(seqs, lowercase = NA), "TRUE or FALSE")
})
