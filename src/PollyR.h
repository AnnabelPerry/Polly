#ifndef POLLY_SI_H
#define POLLY_SI_H

//' PollySI() outputs a matrix of SNP and indel loci at which both the major and minor allele frequency are greater than or equal to an adjustable threshold.
//' @param freq_name A string indicating the VCFtools .frq file. Be sure to include the ".frq" suffix. See Supplemental Materials for instructions on creation.
//' @param desired_scaffolds A vector of scaffolds. Only SNPs and indels on the scaffolds in this vector will be included in the final output.
//' @param threshold Numeric. The minimum frequency an allele must have to be considered highly polymorphic.
// [[Rcpp::export]]
Rcpp::CharacterMatrix PollySI(std::string freq_name, std::vector<std::string> desired_scaffolds, double threshold);

#endif // POLLY_SI_H

#ifndef MICROGENOTYPER_H
#define MICROGENOTYPER_H

//' MicroGenotyper() finds the genotypes of each specified individual at each of the microsatellite loci indicated by the lookup table.
//' @param bam_vec A vector of .bam files, each from a separate individual, to be genotyped.
//' @param lookup_file A string denoting the .csv in which the known microsatellite loci for the species are stored. See Supplemental Materials for instructions on creation.
//' @param desired_scaffolds A vector of scaffolds. Only microsatellites on these scaffolds will be genotyped.
//' @param output_names A vector of CSV files into which the microsatellite genotypes of each .bam file will be outputted. Must be in the same order as the input names specified in bam_vec. Be sure to include ".csv" at the end of each name.
// [[Rcpp::export]]
void MicroGenotyper(std::vector<std::string> bam_vec, std::string lookup_file, std::vector<std::string> desired_scaffolds, std::vector<std::string> output_names);

#endif // MICROGENOTYPER_H

#ifndef POLLYMICROS_H
#define POLLYMICROS_H

//' Identifies user-defined polymorphic microsatellites from the output of MicroGenotyper.
//' @param CSV_names A vector of CSV output files from MicroGenotyper. Be sured to include ".csv" to the end of each name.
//' @param desired_scaffolds A vector of scaffolds. Only microsatellites on these scaffolds will be checked for polymorphism.
//' @param threshold An integer describing the minimum number of alleles a microsatellite must have in the population to be considered polymorphic.
//' @param output_name A string naming the CSV file into which polymorphic microsatellites ought to be listed. Be sure to include ".csv" suffix.
// [[Rcpp::export]]
void PollyMicros(std::vector<std::string> CSV_names, std::vector<std::string> desired_scaffolds, unsigned int threshold, std::string output_name);

#endif // POLLYMICROS_H

#ifndef POLLY_H
#define POLLY_H

//' Performs 3 operations: 1) Genotypes microsatellites. 2). Identifies user-defined polymorphic SNPs, indels, and microsatellites. 3). Groups polymorphic microsatellites, SNPs, and indels into regions of a user-defined size.
//' @param region_size Size (in base pairs) of regions into which polymorphic microsatellites, SNPs, and indels ought to be grouped.
//' @param desired_scaffolds A vector of scaffolds. Only microsatellites, SNPs, and indels on these scaffolds will be genotyped and checked for polymorphism.
//' @param SIthreshold Numeric. The minimum frequency an allele must have to be considered highly polymorphic.
//' @param freq_name A string indicating the VCFtools .frq file. Be sure to include the ".frq" suffix. See Supplemental Materials for instructions on creation.
//' @param bam_vec A vector of .bam files, each from a separate individual, to be genotyped.
//' @param lookup_file A string denoting the .csv in which the known microsatellite loci for the species are stored. See Supplemental Materials for instructions on creation.
//' @param MicroThreshold An integer describing the minimum number of alleles a microsatellite must have in the population to be considered polymorphic.
//' @param SIOutputName A string naming the CSV file into which polymorphic SNPs and indels should be outputted. Be sure to include ".csv" suffix.
//' @param MicroOutputNames A vector of CSV files into which the microsatellite genotypes of each .bam file will be outputted. Must be in the same order as the input names specified in bam_vec. Be sure to include ".csv" at the end of each name.
//' @param PollyMicroOutputName A string naming the CSV file into which polymorphic microsatellites ought to be listed. Be sure to include the ".csv" suffix.
//' @param FinalOutputName A string naming the CSV file into which the highly-polymorphic SNP, indel, and microsatellite clusters ought to be outputted. Be sure to include the ".csv" suffix.
// [[Rcpp::export]]
void Polly(int region_size, std::vector<std::string> desired_scaffolds, double SIthreshold, std::string freq_name, std::vector<std::string> bam_vec, std::string lookup_file, int MicroThreshold, std::string SIOutputName, std::vector<std::string> MicroOutputNames, std::string PollyMicroOutputName, std::string FinalOutputName);

#endif // POLLY_H
