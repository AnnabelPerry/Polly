#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <iterator>
#include <regex>
#include <limits>
#include <map>
#include <stdlib.h>
#include <bits/stdc++.h>
#include "htslib/sam.h"
#include <Rcpp.h>
#include "Polly.h"
#include <sys/resource.h>
#include <stdio.h>
#include <errno.h>
#include <unistd.h>
#include <algorithm>

// Create constants to enable opening of tab and comma separated files
// NOTE: dos2unix CSV file prior to opening
const std::regex comma_it(",");
const std::regex tab("\t");

/*
Function III: Whittle VCF frequency file down to only rows with at least two
alleles of frequency >= threshold
*/
Rcpp::CharacterMatrix PollySI(std::string freq_name, std::vector<std::string> desired_scaffolds, double threshold){
  // Prepare objects for for loop
	std::string colon = ":";
	std::string scaff;
	std::string allele1;
	std::string allele2;

	// Open frequency file as a VoV
	std::vector<std::vector<std::string>> freq_table;
	std::ifstream file(freq_name);
	std::string line{};
	while (file && getline(file, line)) {
	  std::vector<std::string> row{ std::sregex_token_iterator(line.begin(),line.end(),tab,-1), std::sregex_token_iterator() };
	  freq_table.push_back(row);
	}
	file.close();
	// Erase the row with the column names
	freq_table.erase(freq_table.begin());

	// Create distinct vectors for each column
	std::vector<std::string> col0;
	std::vector<std::string> col1;
	std::vector<std::string> col2;
	std::vector<std::string> col3;
	std::vector<std::string> col4;
	std::vector<std::string> col5;

	/*
	Check each known SNP/Indel locus. If the locus has at least two alleles with
	frequency threshold or greater, add the locus to the output 2D vector.
	*/
	for(unsigned int row = 0; row < freq_table.size(); row++){
		/*
		First, make sure the scaffold of the current row is one of the desired
		scaffolds
		*/
		scaff = freq_table[row][0];
		// Check if the row occurs on a desired scaffold
        if(std::count(desired_scaffolds.begin(), desired_scaffolds.end(), scaff) != 0){
    		allele1 = freq_table[row][4];
    		allele1.erase(0, allele1.find(colon) + colon.size());
    		allele2 = freq_table[row][5];
    		allele2.erase(0, allele2.find(colon) + colon.size());

    		// Convert the frequency of each allele to a double.
    		std::stringstream str_freq1(allele1);
    		double freq1 = 0.0;
    		str_freq1 >> freq1;

    		std::stringstream str_freq2(allele2);
    		double freq2 = 0.0;
    		str_freq2 >> freq2;

    		/*
    		 *  Check if both alleles have frequencies greater than threshold.
    		 *   If so, count the row as a polymorphic SNP/indel locus and output
    		 *   its columns to the appropriate vector.
    		 */
    		if((freq1 >= threshold)&&(freq2 >= threshold)){
            col0.push_back(freq_table[row][0]);
            col1.push_back(freq_table[row][1]);
            col2.push_back(freq_table[row][2]);
            col3.push_back(freq_table[row][3]);
            col4.push_back(freq_table[row][4]);
            col5.push_back(freq_table[row][5]);
    			allele1 = "";
    			allele2 = "";
    			continue;
    		}else {
    			allele1 = "";
    			allele2 = "";
    			continue;
    		}
    		allele1 = "";
    		allele2 = "";
    		continue;
		}else{
		    // If the current row's scaffold is not one of the desired scaffolds, don't add it to the final frequencies.
		    continue;
		}
	}
	// Add each column to a master vector
	std::vector<std::string> master_vec = col0;
	for(unsigned int i = 0; i < col1.size(); i++){
	  master_vec.push_back(col1[i]);
	}
	for(unsigned int i = 0; i < col2.size(); i++){
	  master_vec.push_back(col2[i]);
	}
	for(unsigned int i = 0; i < col3.size(); i++){
	  master_vec.push_back(col3[i]);
	}
	for(unsigned int i = 0; i < col4.size(); i++){
	  master_vec.push_back(col4[i]);
	}
	for(unsigned int i = 0; i < col5.size(); i++){
	  master_vec.push_back(col5[i]);
	}
	Rcpp::CharacterMatrix output(col0.size(), 6, master_vec.begin());
	return output;
}

// Function IV: Genotypes microsats in an individual
void MicroGenotyper(std::vector<std::string> bam_vec, std::string lookup_file, std::vector<std::string> desired_scaffolds, std::vector<std::string> output_names){
  // Open lookup table as VoV
  std::vector<std::vector<std::string>> lookup_table;
  std::ifstream file(lookup_file);
  std::string csv_line{};

  while (file && getline(file, csv_line)) {
    std::vector<std::string> row{ std::sregex_token_iterator(csv_line.begin(),csv_line.end(),comma_it,-1), std::sregex_token_iterator() };
    lookup_table.push_back(row);
  }

  // Erases the row with column names
  lookup_table.erase(lookup_table.begin());
  file.close();
  // Objects needed in following loop
  std::string ind_bam;
  const char * char_bam;
  std::string chr;
  std::string lookup_locus;
  std::string allele_pair;
  std::string forward_slash = "/";
  // sam_bp and ID_bp_counter must be intitialized as unsigned integer to avoid downstream compilation errors
  unsigned int sam_bp;
  int motif_bp;
  unsigned int ID_bp_counter;
  int frequent_variant;
  int current_variant;
  int current_copies;
  // Must be set to 0 so you can check if it has been found later
  int second_variant = 0;
  int second_copies = 0;
  int largest_copies = 0;
  int repeats_count = 0;
  bool break_both = false;
  std::string seq_string = "";
  std::vector<std::string> seq_vec;
  std::vector<std::vector<std::string>> mapped_reads;
  int first_pos;
  std::string first_pos_str;
  std::string locus_chr;
  const char * char_locus_chr;


  // Repeat this function for each bam file in the bam vector
  for(unsigned int bam = 0; bam < bam_vec.size(); bam++){

    // EDIT: Delete before final product
    int running_loci_count = 0;
    int read_depth = 0;
    std::ofstream RD_LOCI;
    RD_LOCI.open("RD_LociCount.txt");
//////////////////////////////////////////////////////

    // Store name of the current bam in a single string
    ind_bam = bam_vec[bam];
    // Convert string describing individual bam file to a character constant
    char_bam = ind_bam.c_str();
    // Open a bam file
    samFile *bam_file = sam_open(char_bam, "r");
    // Load the header of the bam file
    bam_hdr_t *header = sam_hdr_read(bam_file);
    // Initialize the alignment of the bam file
    bam1_t *b = bam_init1();

    // File into which genotypes for each locus will be reported
    std::ofstream output_file;
    std::string current_output = output_names[bam];
    output_file.open(current_output);
    output_file << "Scaffold,Locus,Motif,Genotype\n";

    for(unsigned int row = 0; row < lookup_table.size(); row++){
      /*
       * This map will TEMPORARILY house the count of all genotypes with a
       * certain number of repeats in a single sam file.
       * NOTE: Don't get confused! The KEY and VALUE are BOTH integers, but
       * the KEY is the GENOTYPE while the VALUE is the NUMBER OF READS with
       * that genotype.
       */
       // UNTESTED: Clear mapped reads each time you enter a new locus
        // Re-initalize objects needed in for-loop
        ID_bp_counter = 0;
        repeats_count = 0;
        frequent_variant;
        largest_copies = 0;
        second_variant = 0;
        second_copies = 0;
        break_both = false;
        seq_string = "";
        seq_vec.clear();
        mapped_reads.clear();
        first_pos_str = "";
      std::map<int, int> scratch_genotypes;

      chr = lookup_table[row][3];

      /*
       * Check if the current scaffold is one of the desired scaffolds
       */
      if(std::count(desired_scaffolds.begin(), desired_scaffolds.end(), chr) != 0){
        lookup_locus = lookup_table[row][0];
        ///////////////////////////////////////////////////////////////////////
        /////////////// Memory Consuming Step (~5800 KB/locus) ///////////////
        ///////////////////////////////////////////////////////////////////////
        /*
         * Obtain all reads which could possibly house the current look-
         * up locus.
         */
        // Load the current read index (0) in the bam file
        hts_idx_t *idx = sam_index_load(bam_file, char_bam);
        // Create an iterator to move through the bam file
        locus_chr = chr;
        locus_chr.append(":");
        locus_chr.append(lookup_locus);
        locus_chr.append("-");
        locus_chr.append(lookup_locus);
        char_locus_chr = locus_chr.c_str();
        hts_itr_t *iterator = sam_itr_querys(idx, header, char_locus_chr);
        /*
         * So long as there is another iterator after the current one, check that the
         * read id is valid. If so, store the start and end positions in a VoV
         */
        while (sam_itr_next(bam_file, iterator, b) >= 0) {
          if (b->core.tid < 0){
            continue;
          }else if(b->core.qual >= 10){
            first_pos = b->core.pos;
            first_pos_str = std::__cxx11::to_string(b->core.pos);
            seq_vec.push_back(first_pos_str);
            for(int i = 0; i < b->core.l_qseq; i++){
              if(bam_seqi(bam_get_seq(b), i) == 1){
                seq_string.append("A");
              }else if(bam_seqi(bam_get_seq(b), i) == 2){
                seq_string.append("C");
              }else if(bam_seqi(bam_get_seq(b), i) == 4){
                seq_string.append("G");
              }else if(bam_seqi(bam_get_seq(b), i) == 8){
                seq_string.append("T");
              }
            }
            seq_vec.push_back(seq_string);
            seq_string.clear();
            mapped_reads.push_back(seq_vec);
            seq_vec.clear();
          }
        }
        hts_itr_destroy(iterator);
        ///////////////////////////////////////////////////////////////////////

        /*
         * If the vector of mapped reads is empty, then the individual does
         * not have any reads mapped to this locus, so this individual's
         * reads at this locus should be skipped.
         */
        if (mapped_reads.empty()){
          continue;
        }else{
            // EDIT: Delete before finalizing package
            running_loci_count++;
    
            /*
             * Genotype each read, adding to the INDIVIDUAL's genotype
             * table.
             */
            for(unsigned int i = 0; i < mapped_reads.size(); i++){
              // EDIT: Delete before final product
              read_depth += mapped_reads.size();
              
              /*
               * Set the first locus of the microsat equal to the the point
               * at which the first motif SHOULD occur if the current sam
               * sequence has a motif.
               * Some microsat loci are present, but don't begin at exactly
               * the same locus as the reference microsat. So, you want to
               * also check the loci -/+2 bp from the known microsat locus
               * to see if the microsat starts there.
			   * HTSlib read indexing is lowered by 1, while you must search 2 bp prior to
			   * start to account for search window, so subtract by 1 to get appropriate start
			   * locus
               */
              sam_bp = std::__cxx11::stoi(lookup_table[row][0]) - std::__cxx11::stoi(mapped_reads[i][0]) - 1;
              /*
               * Set the first base pair in the motif equal to 0 so you can
               * check if the motif is found in the current sam sequence.
               */
              motif_bp = 0;
              /*
               * Check if the current read has a microsatellite starting
               * anywhere within 5bp of the microsat locus indicated in
               * the lookup table.
               */
              for(int pos = 0; pos < 5; pos++){
                if(mapped_reads[i][1][sam_bp] != lookup_table[row][2][motif_bp]){
                  sam_bp++;
                  continue;
                }else if(mapped_reads[i][1][sam_bp] == lookup_table[row][2][motif_bp]){
                  /*
                   * Begin counting the number of base pairs which are
                   * identical between the lookup motif and the current
                   * sam sequence.
                   */
                  while(sam_bp < (mapped_reads[i][1].size() + 2)){
                    /*
                     * If the current number of identical base pairs
                     * is less than the size of the motif and the
                     * current motif and sam_bp match each other,
                     * count this as an identical base pair and move
                     * to the next base pair.
                     */
                    if((ID_bp_counter < lookup_table[row][2].size()) && ((mapped_reads[i][1][sam_bp] == lookup_table[row][2][motif_bp]))){
                      ID_bp_counter++;
                      motif_bp++;
                      sam_bp++;
                      continue;
                      /*
                       * If a full motif is detected, increment the
                       * count of repeats by one and start at the
                       * beginning of the motif to see if there are any
                       * other repeats in the current sequence.
                       */
                    }else if(ID_bp_counter == lookup_table[row][2].size()){
                      repeats_count++;
                      motif_bp = 0;
                      ID_bp_counter = 0;
                      continue;
                      /*
                       * If repeats have already been detected but the
                       * current sam bp does not align with the
                       * corresponding motif bp OR the read has ended,
                       * the microsat has ended. Thus, you must break
                       * out of JUST the current loop to add the genotype
                       * information to the output 2D vector.
                       */
                    }else if((repeats_count > 0) && ((mapped_reads[i][1][sam_bp] != lookup_table[row][2][motif_bp]) || (sam_bp > mapped_reads[i][1].length()))){
                      motif_bp = 0;
                      break;
                      /*
                       * If repeats have not been detected but the
                       * current sam bp no longer aligns with the motif,
                       * don't output the lookup row's information.
                       * Simply skip to the next read, as this one
                       * lacks a full motif.
                       */
                    }else if((repeats_count == 0)&&((mapped_reads[i][1][sam_bp] != lookup_table[row][2][motif_bp]))){
                      break_both = true;
                      motif_bp = 0;
                      ID_bp_counter = 0;
                      break;
                    }
                  }
                  /*
                   * If the for loop was broken because the sam sequence
                   * did not contain any repeats, move to the next base
                   * pair in position.
                   */
                  if(break_both){
                    break_both = false;
                    continue;
                  }
                  /*
                   * Once you are finished genotyping a read, check if
                   * it has been counted in this individual at this
                   * locus before. If not, add a new entry to the
                   * scratch genotypes and move to the next read in the
                   * file.
                   */
                  std::map<int, int>::iterator rep;
                  rep = scratch_genotypes.find(repeats_count); // Map find searches only for KEYS, so you need not specify ->first
                  if(rep == scratch_genotypes.end()){
                    scratch_genotypes.insert(std::pair<int, int>(repeats_count,1));
                    repeats_count = 0;
                    sam_bp = 0;
                    motif_bp = 0;
                    ID_bp_counter = 0;
                    break;
                    /*
                     * If the genotype of the current read HAS been
                     * recorded previously, increment its count by one,
                     * clear variables, and move to the next read.
                     */
                  }else{
                    scratch_genotypes.find(repeats_count)->second = scratch_genotypes.find(repeats_count)->second + 1;
                    repeats_count = 0;
                    sam_bp = 0;
                    motif_bp = 0;
                    ID_bp_counter = 0;
                    break;
                  }
                  /*
                   If the current sam read doesn't have a microsat motif at the indicated position, skip it.
                   */
                }else{
                  continue;
                }
                continue;
              }
            }
            repeats_count = 0;
            sam_bp = 0;
            motif_bp = 0;
            ID_bp_counter = 0;
            /*
             Once you've genotyped every read (including SE reads) for the current locus of the current individual, identify the two alleles which were
             most common and output them in a single string as "allele1/allele2".
             */
    
            // First, make sure reads were actually genotyped (it is possible that all reads for this locus had very low map quality scores)
            if(!scratch_genotypes.empty()){
              // Iterate through each member of the map and find the value of "second" which is largest
              std::map<int, int>::iterator variant;
              for(variant = scratch_genotypes.begin(); variant != scratch_genotypes.end(); variant++){
                current_variant = variant -> first;
                current_copies = variant -> second;
                if(current_copies > largest_copies){
                  frequent_variant = current_variant;
    
                  largest_copies = current_copies;
                  continue;
                }else{
                  continue;
                }
              }
              // Check if the current most frequent variant has 2x as many reads as any other variant.
              std::map<int, int>::iterator variant2;
              bool two_times_larger = true;
              for(variant2 = scratch_genotypes.begin(); variant2 != scratch_genotypes.end(); variant2++){
                current_variant = variant2 ->first;
                current_copies = variant2 ->second;
                // First, make sure the current number of repeats does not equal the most common number of repeats to ensure you don't count the most common variant's number of copies as NOT two times larger than it's own number of copies.
                if((current_variant != frequent_variant) && (largest_copies >= 2*current_copies)){
                  continue;
                }else if((current_variant != frequent_variant) && (largest_copies < 2*current_copies)){
                  two_times_larger = false;
                  current_copies = 0;
                  break;
                }
              }
              /*
               Homozygous Conditions:
               If the variant which appears on the greatest number of reads has 2x as many reads as any other variant, assume the individual has only
               ONE true genotype and is thus HOMOZYGOUS.
               Output the genotype for this locus and move to the next locus.
               */
              if(two_times_larger){
                allele_pair = std::to_string(frequent_variant);
                allele_pair.append(forward_slash);
                allele_pair.append(std::to_string(frequent_variant));
    
                output_file << chr << "," << lookup_locus << "," << lookup_table[row][2] << "," << allele_pair << std::endl;
                // Clean all relevant objects for the next locus.
                mapped_reads.clear();
                second_copies = 0;
                largest_copies = 0;
                current_copies = 0;
                continue;
                /*
                 Heterozygous Conditions:
                 If the allele with the greatest number of copies does NOT have 2x as many copies as any other allele, find the allele found in the second
                 greatest number of reads and add both to the output.
                 */
              }else{
                // Find the second most frequent variant in the individual's reads for this locus
                std::map<int, int>::iterator variant3;
                for(variant3 = scratch_genotypes.begin(); variant3 != scratch_genotypes.end(); variant3++){
                  current_variant = variant3 -> first;
                  current_copies = variant3 -> second;
                  if((current_variant != frequent_variant) && (current_copies <= largest_copies) && (current_copies >= second_copies)){
                    second_variant = current_variant;
                    second_copies = current_copies;
                    continue;
                  }else{
                    continue;
                  }
                }
                // If a second largest variant was found, output it as an official allele and move to the next locus.
                if(second_variant != 0){
                  allele_pair = std::to_string(frequent_variant);
                  allele_pair.append(forward_slash);
                  allele_pair.append(std::to_string(second_variant));
    
                  output_file << chr << "," << lookup_locus << "," << lookup_table[row][2] << "," << allele_pair << std::endl;
                  /*
                   * Clean all relevant objects for the next locus.
                   */
                  mapped_reads.clear();
                  second_copies = 0;
                  largest_copies = 0;
                  current_copies = 0;
                  continue;
                }
              }
              // If the map quality of ALL reads in this individual were too low to be genotyped, skip to the next locus
            }else{
              continue;
            }
        }
        /*
        * If the current scaffold is not one of the desired scaffolds, don't
        * find genotypes for this locus.
        */
        }else{
            continue;
        }
    }
    bam_destroy1(b);
    bam_hdr_destroy(header);
    sam_close(bam_file);
    //EDIT: Remove before final product
    RD_LOCI << "Total number of microsat loci found in " << ind_bam << ": " << running_loci_count << std::endl << "Total Read Depth of " << ind_bam << ": " << read_depth << std::endl;
    }
}

// Function V: Extracts highly-polymorphic microsats from a set of genotype files

void PollyMicros(std::vector<std::string> CSV_names, std::vector<std::string> desired_scaffolds, unsigned int threshold, std::string output_name){
    // Containers for temporary and final outputs
    std::string temp_string;
    std::string allele_str;
    std::string frq_str;
    std::ofstream output;
    output.open(output_name);
    output << "Scaffold,Locus,Motif,Alleles,Allele Frequencies\n";
    // Doubles used in calculating and storing allele frequencies
    double current_allele;
    double total_alleles;
    double temp_frq;

    std::string colon = ":";
    std::string slash = "/";
    // The genotype and two alleles reported in a row
    std::string geno;
    std::string allele1;
    std::string allele2;
    // Vector into which new alleles are temporarily stored
    std::vector<std::string> temp_vec;
    // Vector into which unique alleles will temporarily be stored
    std::vector<std::string> unique_alleles;

    /*
    * Map describing all alleles at a locus as well as their number of occurrences
    * in the sample.
    * Key: Scaffold:Locus:Motif
    * Value: Vector of alleles
    */
    std::map<std::string, std::vector<std::string>> Alleles;
    std::string key;
    std::string map_val;

    // Iterate through each CSV file
    for(unsigned int i = 0; i < CSV_names.size(); i++){
      // Open CSV file
      std::vector<std::vector<std::string>> indv;
      std::ifstream file(CSV_names[i]);
      std::string csv_line{};
      while (file && getline(file, csv_line)) {
        std::vector<std::string> row{ std::sregex_token_iterator(csv_line.begin(),csv_line.end(),comma_it,-1), std::sregex_token_iterator() };
        indv.push_back(row);
      }
      // Erases the row with column names
      indv.erase(indv.begin());
      file.close();

      // Check each row in the current CSV
        for(unsigned int row = 0; row < indv.size(); row++){
            // Check if the row occurs on a desired scaffold
            if(std::count(desired_scaffolds.begin(), desired_scaffolds.end(), indv[row][0]) != 0){
                /*
                * Check if this is the first CSV.
                * If this is the first CSV, initialize a new map key:value for
                * each new row.
                */
                if(i == 0){
                    // Build key
                    key = indv[row][0];
                    key.append(colon);
                    key.append(indv[row][1]);
                    key.append(colon);
                    key.append(indv[row][2]);

                    // Build value entries
                    geno = indv[row][3];
                    allele1 = geno.substr(0, geno.find(slash));
                    geno.erase(0, geno.find(slash) + slash.size());
                    allele2 = geno;
                    temp_vec.push_back(allele1);
                    temp_vec.push_back(allele2);

                    // Initialize map with key and values and move to next row
                    Alleles.insert({key, temp_vec});
                    key.clear();
                    temp_vec.clear();
                    continue;
                /*
                * If this is NOT the first CSV, only initialize a new key:value
                * if the current key has yet to be recorded
                */
                }else{
                    // Build key
                    key = indv[row][0];
                    key.append(colon);
                    key.append(indv[row][1]);
                    key.append(colon);
                    key.append(indv[row][2]);

                    // Build value entries
                    geno = indv[row][3];
                    allele1 = geno.substr(0, geno.find(slash));
                    geno.erase(0, geno.find(slash) + slash.size());
                    allele2 = geno;

                    // If key is present, append alleles to corresponding vector
                    std::map<std::string,std::vector<std::string>>::iterator it;
                    it = Alleles.find(key);
                    if(it != Alleles.end()){
                        it->second.push_back(allele1);
                        it->second.push_back(allele2);
                        continue;
                    // If key is NOT present, initialize new map entry
                    }else{
                        temp_vec.push_back(allele1);
                        temp_vec.push_back(allele2);
                        Alleles.insert({key, temp_vec});
                        key.clear();
                        temp_vec.clear();
                        continue;
                    }
                }
            // If the row does not occur on a desired scaffold, skip it
            }else{
                continue;
            }
        }
    }
    // Iterate through each entry in the map
    std::map<std::string,std::vector<std::string>>::iterator it2;
    for(it2 = Alleles.begin(); it2 != Alleles.end(); it2++){
        temp_vec = it2->second;
        unique_alleles = temp_vec;
        std::vector<std::string>::iterator uni_it;
        // Count number of unique alleles at this locus in the sample
        std::sort(unique_alleles.begin(), unique_alleles.end());
        // Resizing the vector so to remove the undefined terms
        auto iter = std::unique(unique_alleles.begin(), unique_alleles.end());
        unique_alleles.erase(iter, unique_alleles.end());
        // Check if number of unique alleles is above threshold for polymorphism
        if(unique_alleles.size() >= threshold){
            /*
            * Begin outputting polymorphic locus
            * Store key in temporary string.
            */
            temp_string = it2->first;
            // Add scaffold to output file
            output << temp_string.substr(0, temp_string.find(colon)) << ",";
            temp_string.erase(0, temp_string.find(colon) + colon.size());
            // Add locus to output vector
            output << temp_string.substr(0, temp_string.find(colon)) << ",";
            temp_string.erase(0, temp_string.find(colon) + colon.size());
            // Add motif to output vector
            output << temp_string << ",";
            temp_string.clear();

            // Store total number of alleles
            total_alleles = temp_vec.size();
            // Count number of occurrences of each unique allele
            for(unsigned int a = 0; a < unique_alleles.size(); a++){
                /*
                * Find number of occurrences of the unique allele in the WHOLE
                * sample.
                * Multiply number of occurrences by 1.0 to convert to a double.
                */
                current_allele = 1.0*std::count(temp_vec.begin(), temp_vec.end(), unique_alleles[a]);
                // Calculate allele frequency for current allele
                temp_frq = current_allele/total_alleles;
                /*
                * If you are at the last allele, finalize the running strings
                * and output the vector to the matrix.
                */
                if(a == unique_alleles.size()-1){
                    allele_str.append(unique_alleles[a]);
                    frq_str.append(std::to_string(temp_frq));
                    output << allele_str << ",";
                    output << frq_str << std::endl;
                    frq_str.clear();
                    allele_str.clear();
                    break;
                /*
                * If you are NOT at the last allele, add the allele frequency and
                * allele to the running strings of alleles and frequencies.
                */
                }else{
                    allele_str.append(unique_alleles[a]);
                    allele_str.append(slash);
                    frq_str.append(std::to_string(temp_frq));
                    frq_str.append(slash);
                    continue;
                }
            }
        // If not, this locus is not polymorphic and should be skipped
        }else{
            unique_alleles.clear();
            temp_vec.clear();
            continue;
        }
    }
    output.close();
}

/*
* Function VI: Groups polymorphic SNPs, indels, and microsats into regions, then
* performs a posthoc to filter out overlapping regions and to summarize data.
*/
void Polly(int region_size, std::vector<std::string> desired_scaffolds, std::string PollySIInputName, std::string PollyMicrosInputName, std::string FinalOutputName){
    ////////////////////////////// Define Objects //////////////////////////////
    std::string colon = ":";
    std::string open_bracket = "[";
    std::string comma = ", ";
    std::string slash = "/";
    std::string scaff;
    std::string SI_element0;
    std::string SI_element1;
    std::string SI_element4;
    std::string SI_element5;
	std::string allele_freqs;
    std::string loci_string;
    std::string motif_repeat;
    std::string motif;
    std::string temp_alleles;
    std::string temp_freqs;
    std::string sub;
    
    bool break_both;
    bool first_cent_allele;
    bool first_noncent_allele;
    
    int central_locus;
    int temp_locus;
	int loci_count = 0;
	int allele_count_int;
	int SNP_count = 0;
	int indel_count = 0;
	int di = 0;
	int tri = 0;
	int tetra = 0;
	int penta_plus = 0;
	int repeat = 0;
	int cent_longest_allele = 0;
	int longest_allele = 0;
	int min_diff = 0;
	
	std::vector<int> repeats;
	std::vector<int> temp_allele_pairs;
	std::vector<int> abs_diffs;
	
	std::vector<std::vector<int>> allele_pairs;
	
	// Unsigned integers to enable comparison to size integers in RCpp
	unsigned int temp_UI = 0;
	unsigned int temp_UI2 = 0;
	
    ////////////////////////// Input/Output Prep     ///////////////////////////
    // Open the table of highly polymorphic microsatellites 
    std::vector<std::vector<std::string>> PM_table;
    std::ifstream file1(PollyMicrosInputName);
    std::string csv_line1{};
    while (file1 && getline(file1, csv_line1)) {
        std::vector<std::string> row1{ std::sregex_token_iterator(csv_line1.begin(),csv_line1.end(),comma_it,-1), std::sregex_token_iterator() };
        PM_table.push_back(row1);
    }
    // Erases the row with column names
    PM_table.erase(PM_table.begin());
    file1.close();
    
    // Open the table of highly polymorphic SNPs and Indels.
    std::vector<std::vector<std::string>> SI_table;
    std::ifstream file2(PollySIInputName);
    std::string csv_line2{};
    while (file2 && getline(file2, csv_line2)) {
        std::vector<std::string> row2{ std::sregex_token_iterator(csv_line2.begin(),csv_line2.end(),comma_it,-1), std::sregex_token_iterator() };
        SI_table.push_back(row2);
    }
    // Erases the row with column names
    SI_table.erase(SI_table.begin());
    file2.close();
    
    // Open the file into which the final output will be read.
	std::ofstream FinalOutput;
    FinalOutput.open(FinalOutputName);
    FinalOutput << "Scaffold,Total Loci (Count),Total Alleles (Count),SNP Count,Indel Count,Dinucleotide Allele Count,Trinucleotide Allele Count,Tetranucleotide Allele Count,Pentanucleotide (and Above) Allele Count,Loci (Positions),Allele:Frequency,Length of Longest Allele (Central Microsatellite),Minimum Difference in Repeats (Central Microsatellite)\n";

    
    //////////////////////// Group Polymorphic Regions /////////////////////////

    for(int microsat1 = 0; microsat1 < PM_table.size(); microsat1++){
        scaff = PM_table[microsat1][0];
    	/*
    	* Check if the current scaffold is a desired scaffold
    	*/
    	if(std::find(desired_scaffolds.begin(), desired_scaffolds.end(), scaff) != desired_scaffolds.end()){
    	    // Wipe all counters to start fresh for the new CENTRAL microsatellite
            SI_element0.clear();
            SI_element1.clear();
            SI_element4.clear();
            SI_element5.clear();
            allele_pairs.clear();
        	allele_freqs = open_bracket;
            loci_string = PM_table[microsat1][1];
            motif_repeat.clear();
            motif = PM_table[microsat1][2];
            temp_alleles = PM_table[microsat1][3];
            temp_freqs = PM_table[microsat1][4];
            central_locus = stoi(PM_table[microsat1][1]);
            
            loci_count = 1;
        	SNP_count = 0;
        	indel_count = 0;
        	di = 0;
        	tri = 0;
        	tetra = 0;
        	penta_plus = 0;
        	repeat = 0;
        	cent_longest_allele = 0;
        	longest_allele = 0;
        	min_diff = 0;
        	abs_diffs.clear();
        	repeats.clear();
        	break_both = false;
        	allele_count_int = 0;
        	
            /*
            * Upon entering a central microsat, create a new map of ranges of longest
            * microsats within <region_size> of the central microsat so you can
            * later remove any SNPs/Indels which overlap with one or more of the
            * microsats in the reigon.
            */
            std::map<int, int> loci_range;
            first_cent_allele = true;

    		/*
    		* Start the running string of alleles and frequencies
    		* The end goal is a string like this:
    		* [allele1.1:freq1.1/allele1.2:freq1.2], [allele2.1:freq2.1/allele2.2:freq2.2]
    		*/
    		while(!temp_alleles.empty()){
    		    if(temp_alleles.find(slash) != std::string::npos){
    		        motif_repeat = motif;
        			motif_repeat.append(temp_alleles.substr(0,temp_alleles.find(slash)));
        			motif_repeat.append(colon);
        			motif_repeat.append(temp_freqs.substr(0,temp_freqs.find(slash)));
        			motif_repeat.append(slash);
        			allele_freqs.append(motif_repeat);
        			
        			repeat = stoi(temp_alleles.substr(0,temp_alleles.find(slash)));
                    repeats.push_back(repeat);
                    
                    
    		        /*
					* If this is the first allele, store the current microsat
					* locus as a new entry in the map of loci ranges and store
					* the current repeat*motif as the longest repeat
					*/
					if(first_cent_allele){
					    first_cent_allele = false;
					    cent_longest_allele = repeat*motif.size();
						loci_range.insert(std::pair<int, int>(central_locus, cent_longest_allele + central_locus));
					/*
					* If the current microsat locus HAS been stored in the map
					* already, check if end locus of the current repeat number
					* is greater than the end locus of the previous longest repeat.
					* If so, store the current repeat's range as the largest
					* range for the microsat locus.
					*/
					}else{
    					// Create unsigned integer so later comparison is not invalid in RCpp
    					temp_UI = loci_range.find(central_locus)->second;
					    if(((repeat*motif.size()) + central_locus) > temp_UI){
					        cent_longest_allele = repeat*motif.size();
    						loci_range.find(central_locus)->second = cent_longest_allele + central_locus;
					    }
					}
        			temp_alleles.erase(0,temp_alleles.find(slash) + slash.size());
        			temp_freqs.erase(0,temp_freqs.find(slash) + slash.size());
        			
        			/*
        			* Increment the appropriate motif size count and continue
        			* to the next allele
        			*/
					if(motif.size() == 2){
						di++;
						allele_count_int++;
						continue;
					}else if(motif.size() == 3){
						tri++;
						allele_count_int++;
						continue;
					}else if(motif.size() == 4){
						tetra++;
						allele_count_int++;
						continue;
					}else if(motif.size() >= 5){
						penta_plus++;
						allele_count_int++;
						continue;
					}
    		    }else{
    		        motif_repeat = motif;
        			motif_repeat.append(temp_alleles);
        			motif_repeat.append(colon);
        			motif_repeat.append(temp_freqs);
        			motif_repeat.append("]");
                	allele_freqs.append(motif_repeat);
    		        temp_freqs.clear();
    		        
    		        repeat = stoi(temp_alleles);
                    repeats.push_back(repeat);
                    
    		        temp_alleles.clear();
    		        /*
					* If this is the first allele, store the current microsat
					* locus as a new entry in the map of loci ranges and store
					* the current repeat*motif as the longest repeat
					*/
					if(first_cent_allele){
					    first_cent_allele = false;
					    cent_longest_allele = repeat*motif.size();
						loci_range.insert(std::pair<int, int>(central_locus, cent_longest_allele + central_locus));
					/*
					* If the current microsat locus HAS been stored in the map
					* already, check if end locus of the current repeat number
					* is greater than the end locus of the previous longest repeat.
					* If so, store the current repeat's range as the largest
					* range for the microsat locus.
					*/
					}else{
    					// Create unsigned integer so later comparison is not invalid in RCpp
    					temp_UI = loci_range.find(central_locus)->second;
					    if(((repeat*motif.size()) + central_locus) > temp_UI){
					        cent_longest_allele = repeat*motif.size();
    						loci_range.find(central_locus)->second = cent_longest_allele + central_locus;
					    }
					}
					
    		        /*
    		        * Increment the appropriate motif size count and stop
    		        * searching for alleles at this locus
    		        */
					if(motif.size() == 2){
						di++;
						allele_count_int++;
    		            break;
					}else if(motif.size() == 3){
						tri++;
						allele_count_int++;
    		            break;
					}else if(motif.size() == 4){
						tetra++;
						allele_count_int++;
    		            break;
					}else if(motif.size() >= 5){
						penta_plus++;
						allele_count_int++;
    		            break;
					}
    		    }
    		}
    		
    		/*
			* Find and store the minimum difference in repeat length between
			* microsat alleles of the central microsat.
			* Collect each possible combination of repeats.
			*/
			for(int i = 0; i < repeats.size(); i++){
        	    for(int j = 0; j < repeats.size(); j++){
        	        if(j == i){
        	            continue;
        	        }else{
            	        temp_allele_pairs.push_back(repeats[i]);
            		    temp_allele_pairs.push_back(repeats[j]);
            		    allele_pairs.push_back(temp_allele_pairs);
        	            temp_allele_pairs.clear();
        	        }
        		}
        	continue;
        	}
        	
        	
        	// Store absolute value of differences in each allele length in vector
			for(int p = 0; p < allele_pairs.size(); p++){
				abs_diffs.push_back(abs(allele_pairs[p][0] - allele_pairs[p][1]));
				continue;
			}
			
			// Find smallest difference and store as smallest difference
			min_diff = abs_diffs[0];
			for(int r = 1; r < abs_diffs.size(); r++){
			    if((abs_diffs[r] < min_diff) & (abs_diffs[r] != 0)){
			        min_diff = abs_diffs[r];
			    }
			    continue;
			}
			
			/*
    		* Now find all microsats occurring within <region_size> of the
    		* current microsat.
    		*/
    		for(int microsat2 = 0; microsat2 < PM_table.size(); microsat2++){
    			/*
    			* Check that microsat2...
    			*	1. occurs on the same scaffold
    			*   2. does NOT match microsat1 and
    			*	3. occurs within <region_size>bp of microsat1.
    			* If all 3 conditions are met, record the microsat's alleles,
    			* frequencies, and locus.
    			*/
    			temp_locus = stoi(PM_table[microsat2][1]);
    			first_noncent_allele = true;
    			if((PM_table[microsat1][0] == PM_table[microsat2][0]) && (microsat1 != microsat2) && (abs(central_locus - temp_locus) < region_size/2)){
    			    motif = PM_table[microsat2][2];
    				// Record the microsat2 locus.
    				loci_string.append(comma);
    				loci_string.append(PM_table[microsat2][1]);
    				loci_count++;
    				// Start the new box of alleles
    				allele_freqs.append(comma);
    				allele_freqs.append(open_bracket);
                    temp_alleles = PM_table[microsat2][3];
                    temp_freqs = PM_table[microsat2][4];
    				
        			/*
        			* Add the motif, repeats, and frequencies of the current
        			* microsat locus to the allele_freqs string
        			*/
                    while(!temp_alleles.empty()){
            		    if(temp_alleles.find(slash) != std::string::npos){
            		        motif_repeat = motif;
                			motif_repeat.append(temp_alleles.substr(0,temp_alleles.find(slash)));
                			motif_repeat.append(colon);
                			motif_repeat.append(temp_freqs.substr(0,temp_freqs.find(slash)));
                			motif_repeat.append(slash);
                			allele_freqs.append(motif_repeat);
                			
                			repeat = stoi(temp_alleles.substr(0,temp_alleles.find(slash)));
            		        /*
        					* If this is the first allele, store the current microsat
        					* locus as a new entry in the map of loci ranges
        					*/
        					sub = allele_freqs.substr(allele_freqs.size() - 1, allele_freqs.size()-1);
        					if(first_noncent_allele){
        					    first_noncent_allele = false;
        					    longest_allele = repeat*motif.size();
        						loci_range.insert(std::pair<int, int>(temp_locus, longest_allele + temp_locus));
        					/*
        					* If the current microsat locus HAS been stored in the map
        					* already, check if end locus of the current repeat number
        					* is greater than the end locus of the previous longest repeat.
        					* If so, store the current repeat's range as the largest
        					* range for the microsat locus.
        					*/
        					}else{
            					// Create unsigned integer so later comparison is not invalid in RCpp
            					temp_UI2 = loci_range.find(temp_locus)->second;
        					    if(((repeat*motif.size()) + temp_locus) > temp_UI2){
        					        longest_allele = repeat*motif.size();
            						loci_range.find(temp_locus)->second = longest_allele + temp_locus;
        					    }
        					}
                			temp_alleles.erase(0,temp_alleles.find(slash) + slash.size());
                			temp_freqs.erase(0,temp_freqs.find(slash) + slash.size());
                			
                			/*
                			* Increment the appropriate motif size count and continue
                			* to the next allele
                			*/
        					if(motif.size() == 2){
        						di++;
        						allele_count_int++;
        						continue;
        					}else if(motif.size() == 3){
        						tri++;
        						allele_count_int++;
        						continue;
        					}else if(motif.size() == 4){
        						tetra++;
        						allele_count_int++;
        						continue;
        					}else if(motif.size() >= 5){
        						penta_plus++;
        						allele_count_int++;
        						continue;
        					}
            		    }else{
            		        motif_repeat = motif;
                			motif_repeat.append(temp_alleles);
                			motif_repeat.append(colon);
                			motif_repeat.append(temp_freqs);
                			motif_repeat.append("]");
                			allele_freqs.append(motif_repeat);
            		        temp_freqs.clear();
            		        
            		        repeat = stoi(temp_alleles);
                            repeats.push_back(repeat);
                            
            		        temp_alleles.clear();
    		        
            		        /*
        					* If this is the first allele, store the current microsat
        					* locus as a new entry in the map of loci ranges and store
        					* the current repeat*motif as the longest repeat
        					*/
        					if(allele_freqs.find(comma) == std::string::npos){
        					    longest_allele = repeat*motif.size();
        						loci_range.insert(std::pair<int, int>(temp_locus, longest_allele + temp_locus));
            					/*
            					* If the current microsat locus HAS been stored in the map
            					* already, check if end locus of the current repeat number
            					* is greater than the end locus of the previous longest repeat.
            					* If so, store the current repeat's range as the largest
            					* range for the microsat locus.
            					*/
        					}else{
            					// Create unsigned integer so later comparison is not invalid in RCpp
            					temp_UI = loci_range.find(temp_locus)->second;
        					    if(((repeat*motif.size()) + temp_locus) > temp_UI){
        					        longest_allele = repeat*motif.size();
            						loci_range.find(temp_locus)->second = longest_allele + temp_locus;
        					    }
        					}
            		        /*
            		        * Increment the appropriate motif size count and stop
            		        * searching for alleles at this locus
            		        */
        					if(motif.size() == 2){
        						di++;
        						allele_count_int++;
            		            break;
        					}else if(motif.size() == 3){
        						tri++;
        						allele_count_int++;
            		            break;
        					}else if(motif.size() == 4){
        						tetra++;
        						allele_count_int++;
            		            break;
        					}else if(motif.size() >= 5){
        						penta_plus++;
        						allele_count_int++;
            		            break;
        					}
            		    }
            		}
    			}else{
    			    continue;
    			}
    		}
        		
    		/*
    		* Once you've found all microsats within <region_size>bp of
    		* microsat1, find the SNPs and Indels occurring within
    		* <region_size>bp of microsat1
    		*/
    		for(int row = 0; row < SI_table.size(); row++){
    			/*
    			* If the SNP/Indel occurs within <region_size>bp of microsat1, 
    			* record it, up the loci count by 1, and up the allele count by 2.
    			* Check the identity of the SNP/Indel and up the appropriate
    			* count by 1.
				*/
    			SI_element0 = SI_table[row][0];
    			SI_element1 = SI_table[row][1];
    			if((scaff == SI_element0) && (abs(central_locus - stoi(SI_element1)) < region_size/2)){
    				/*
    				* Before recording an SNP or indel, ensure it does not 
    				* overlap with any pre-recorded microsats in the region.
    				*/
    				std::map<int, int>::iterator locus_range;
    				for(locus_range = loci_range.begin(); locus_range != loci_range.end(); locus_range++){
    				    SI_element1 = SI_table[row][1];
    				    if((stoi(SI_element1) >= locus_range->first) && (stoi(SI_element1) <= locus_range->second)){
    					    break_both = true;
    						break;
    					}else{
    					    continue;
    					}
    				}
    				/*
    				* If the SNP/Indel locus overlaps with a pre-recorded
    				* microsatellite, don't record it.
    				*/
    				if(break_both){
    				    break_both = false;
    				    continue;
    				/*
    				* If the SNP/Indel locus does NOT overlap with a pre-
    				* recorded microsatellite, record it.
    				*/
    				}else{
    					allele_freqs.append(", [");
        				allele_freqs.append(SI_table[row][4]);
        				allele_freqs.append(slash);
        				allele_freqs.append(SI_table[row][5]);
        				allele_freqs.append("]");
        				loci_string.append(comma);
        				loci_string.append(SI_table[row][1]);
        				loci_count++;
        				allele_count_int += 2;
        				// If both alleles are size 1, the locus is a SNP
        				SI_element4 = SI_table[row][4];
        				SI_element5 = SI_table[row][5];
    					if((SI_element4.substr(0,SI_element4.find(colon)).size() == 1) && (SI_element5.substr(0,SI_element5.find(colon)).size() == 1)){
    						SNP_count++;
    						continue;
    					// If one or more alleles has size > 1, it is an indel
    					}else{
    						indel_count++;
    						continue;
    					}
					}
				/*
				* If the SNP/Indel does not occur within <region_size>bp of
				* microsat1, simply move to the next row.
				*/
				}else{
					continue;
				}
			}
    	/*
    	* Once you've found every variant within <region_size>bp of microsat1,
    	* output the row to the output file and move to the next microsat.
    	*/

    	FinalOutput << scaff << "," << std::__cxx11::to_string(loci_count) << "," << std::__cxx11::to_string(allele_count_int) << "," << std::__cxx11::to_string(SNP_count) << "," << std::__cxx11::to_string(indel_count) << "," << std::__cxx11::to_string(di) << "," << std::__cxx11::to_string(tri) << "," << std::__cxx11::to_string(tetra) << "," << std::__cxx11::to_string(penta_plus) << "," << loci_string << "," << allele_freqs << "," << std::__cxx11::to_string(cent_longest_allele) << "," << std::__cxx11::to_string(min_diff) << std::endl;
    	continue;
    	}
	}
    FinalOutput.close();
}
