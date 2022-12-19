## usethis namespace: start
#' @useDynLib Polly, .registration = TRUE
## usethis namespace: end
NULL

## usethis namespace: start
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL

#' PollySI() outputs a matrix of SNP and indel loci at which both the major and minor allele frequency are greater than or equal to an adjustable threshold.
#' @export
PollySI <- function(freq_name, desired_scaffolds, threshold) {
  .Call(`_Polly_PollySI`, freq_name, desired_scaffolds, threshold)
}

#' MicroGenotyper() finds the genotypes of each specified individual at each of the microsatellite loci indicated by the lookup table.
#' @export
MicroGenotyper <- function(bam_vec, lookup_file, desired_scaffolds, output_names) {
  invisible(.Call(`_Polly_MicroGenotyper`, bam_vec, lookup_file, desired_scaffolds, output_names))
}

#' Identifies user-defined polymorphic microsatellites from the output of MicroGenotyper.
#' @export
PollyMicros <- function(CSV_names, desired_scaffolds, threshold, output_name) {
  .Call(`_Polly_PollyMicros`, CSV_names, desired_scaffolds, threshold, output_name)
}

#' Combines the polymorphic microsatellites detected by PollyMicros() with the the polymorphic SNPs and indels detected by PollySI() to group polymorphic microsatellites, SNPs, and indels into regions of a user-defined size.
#' @export
Polly <- function(region_size, desired_scaffolds, PollySIInputName, PollyMicrosInputName, FinalOutputName) {
  invisible(.Call(`_Polly_Polly`, region_size, desired_scaffolds, PollySIInputName, PollyMicrosInputName, FinalOutputName))
}