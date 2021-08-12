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

#' Performs 3 operations: 1) Genotypes microsatellites. 2). Identifies user-defined polymorphic SNPs, indels, and microsatellites. 3). Groups polymorphic microsatellites, SNPs, and indels into regions of a user-defined size.
#' @export
Polly <- function(region_size, desired_scaffolds, SIthreshold, freq_name, bam_vec, lookup_file, MicroThreshold, SIOutputName, MicroOutputNames, PollyMicroOutputName, FinalOutputName) {
  invisible(.Call(`_Polly_Polly`, region_size, desired_scaffolds, SIthreshold, freq_name, bam_vec, lookup_file, MicroThreshold, SIOutputName, MicroOutputNames, PollyMicroOutputName, FinalOutputName))
}
