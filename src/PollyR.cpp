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
#include "PollyR.h"

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
        // Store name of the current bam in a single string
        ind_bam = bam_vec[bam];
        // Convert string describing individual bam file to a character constant
        char_bam = ind_bam.c_str();

    	// File into which genotypes for each locus will be reported
    	std::ofstream output_file;
    	std::string current_output = output_names[bam];
	    output_file.open(current_output);
	    output_file << "Scaffold,Locus,Motif,Genotype\n";

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

    	for(unsigned int row = 0; row < lookup_table.size(); row++){
          /*
        	* This map will TEMPORARILY house the count of all genotypes with a
        	* certain number of repeats in a single sam file.
        	* NOTE: Don't get confused! The KEY and VALUE are BOTH integers, but
        	* the KEY is the GENOTYPE while the VALUE is the NUMBER OF READS with
        	* that genotype.
        	*/
        	std::map<int, int> scratch_genotypes;

    		chr = lookup_table[row][3];

    		/*
    		* Check if the current scaffold is one of the desired scaffolds
    		*/
    		if(std::count(desired_scaffolds.begin(), desired_scaffolds.end(), chr) != 0){
        		lookup_locus = lookup_table[row][0];
                /*
                * Obtain all reads which could possibly house the current look-
                * up locus.
                */
                // Open a bam file
                samFile *bam_file = sam_open(char_bam, "r");
                // Load the current read index (0) in the bam file
                hts_idx_t *idx = sam_index_load(bam_file, char_bam);
                // Load the header of the bam file
                bam_hdr_t *header = sam_hdr_read(bam_file);
                // Initialize the alignment of the bam file
                bam1_t *b = bam_init1();
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
                bam_destroy1(b);
                hts_itr_destroy(iterator);
                bam_hdr_destroy(header);
                sam_close(bam_file);

        	    /*
        		* If the vector of mapped reads is empty, then the individual does
        		* not have any reads mapped to this locus, so this individual's
        		* reads at this locus should be skipped.
        		*/
        		if (mapped_reads.empty()){
        			continue;
        		}

        		/*
        		* Genotype each read, adding to the INDIVIDUAL's genotype
        		* table.
        		*/
        		for(unsigned int i = 0; i < mapped_reads.size(); i++){
        			/*
        			* Set the first locus of the microsat equal to the the point
        			* at which the first motif SHOULD occur if the current sam
        			* sequence has a motif.
        			* Some microsat loci are present, but don't begin at exactly
        			* the same locus as the reference microsat. So, you want to
        			* also check the loci -/+2 bp from the known microsat locus
        			* to see if the microsat starts there.
        			*/
        			sam_bp = std::__cxx11::stoi(lookup_table[row][0]) - std::__cxx11::stoi(mapped_reads[i][0]);
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
        			    allele_pair.append("/");
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
            			    allele_pair.append("/");
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
        	/*
        	* If the current scaffold is not one of the desired scaffolds, don't
        	* find genotypes for this locus.
        	*/
        	}else{
        	    continue;
        	}
    	}
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
* Function VI: Extracts and groups polymorphic SNPs, indels, and microsats into
* regions, then performs a posthoc to filter out overlapping regions and to summarize
* data.
*/
void Polly(int region_size, std::vector<std::string> desired_scaffolds, double SIthreshold, std::string freq_name, std::vector<std::string> bam_vec, std::string lookup_file, int MicroThreshold, std::string SIOutputName, std::vector<std::string> MicroOutputNames, std::string PollyMicroOutputName, std::string FinalOutputName){
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

    // Objects to be used in the for loop
    std::string period = ".";
    std::string percent = "%";
    std::string colon = ":";
    std::string open_bracket = "[";
    std::string comma = ", ";
    std::string space = " ";
    std::string chr;
    std::string lookup_locus;
    std::string allele_identifier;
    std::string locus_identifier;
    std::string str_row;
    std::string allele_pair;
    std::string SI_element0;
    std::string SI_element1;
    std::string SI_element4;
    std::string SI_element5;
	// sam_bp and ID_bp_counter must be defined as unsigned integer to avoit future compilation errors
    unsigned int sam_bp;
    int motif_bp;
    unsigned int ID_bp_counter;
    int repeats_count;
	int frequent_variant;
	int largest_copies;
	int current_variant;
	int current_copies;
    int second_variant;
    int second_copies;
    bool break_both;
    bool allele_found;
    bool locus_found;
    std::string ind_bam;
    const char * char_bam;
    std::string seq_string = "";
    std::vector<std::string> seq_vec;
    std::vector<std::vector<std::string>> mapped_reads;
    int first_pos;
    std::string first_pos_str;
    std::string locus_chr;
    const char * char_locus_chr;

    /*
    * For each Scaffold%Locus pair, this map stores every possible allele
    * across the sample space as well as the number of times that allele
    * was found.
    * Outer Key: Scaffold%Locus
    * Inner Key: MotifRepeat
    * Inner Value: Number of observations
    */
    std::map<std::string, std::map<std::string, int>> allele_counts;

    // Repeat this function for each bam file in the bam vector
    for(unsigned int bam = 0; bam < bam_vec.size(); bam++){
        // Objects which must be reset to with each new individual
    	repeats_count = 0;
	    largest_copies = 0;
        second_variant = 0;
        second_copies = 0;
        break_both = false;
        ind_bam = bam_vec[bam];
        char_bam = ind_bam.c_str();

    	// File into which genotypes for each locus will be reported
    	std::ofstream MicroOutput;
    	std::string current_output = MicroOutputNames[bam];
	    MicroOutput.open(current_output);
	    MicroOutput << "Scaffold,Locus,Motif,Genotype\n";
    	for(unsigned int row = 0; row < lookup_table.size(); row++){
    	    /*
        	* This map will TEMPORARILY house the count of all genotypes with a
        	* certain number of repeats in a single sam file. Key is the number of
        	* repeats, value is the count of reads with that genotype.
        	*/
        	std::map<int, int> scratch_genotypes;

    		chr = lookup_table[row][3];
    		mapped_reads.clear();
    		/*
    		Check if the current scaffold is one of the desired scaffolds
    		*/
    		if(std::count(desired_scaffolds.begin(), desired_scaffolds.end(), chr) != 0){
        		lookup_locus = lookup_table[row][0];
        		/*
        		* The locus identifier is the string which uniquely identifies
        		* the present locus.
        		* Locus identifiers look like this: ScyDAA6_1_HRSCAF_23%2874
        		*/
        		locus_identifier = chr;
        		locus_identifier.append(percent);
        		locus_identifier.append(lookup_locus);
        		/*
                * Obtain all reads which could possibly house the current look-
                * up locus.
                */
                // Open a bam file
                samFile *bam_file = sam_open(char_bam, "r");
                // Load the current read index (0) in the bam file
                hts_idx_t *idx = sam_index_load(bam_file, char_bam);
                // Load the header of the bam file
                bam_hdr_t *header = sam_hdr_read(bam_file);
                // Initialize the alignment of the bam file
                bam1_t *b = bam_init1();
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
                bam_destroy1(b);
                hts_itr_destroy(iterator);
                bam_hdr_destroy(header);
                sam_close(bam_file);

        	    /*
        		* If the vector of mapped reads is empty, then the individual does
        		* not have any reads mapped to this locus, so this individual's
        		* reads at this locus should be skipped.
        		*/
        		if (mapped_reads.empty()){
        			continue;
        		}

        		/*
        		* Genotype each read, adding to the INDIVIDUAL's genotype
        		* table.
        		*/
        		for(unsigned int i = 0; i < mapped_reads.size(); i++){
        			/*
        			* Set the first locus of the microsat equal to the the point
        			* at which the first motif SHOULD occur if the current sam
        			* sequence has a motif.
        			* Some microsat loci are present, but don't begin at exactly
        			* the same locus as the reference microsat. So, you want to
        			* also check the loci -/+2 bp from the known microsat locus
        			* to see if the microsat starts there.
        			*/
        			sam_bp = std::__cxx11::stoi(lookup_table[row][0]) - std::__cxx11::stoi(mapped_reads[i][0]);
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
        			* Homozygous Conditions:
        			* If the variant which appears on the greatest number of reads
        			* has 2x as many reads as any other variant, assume the
        			* individual has only ONE true genotype and is thus HOMOZYGOUS.
        			* Output the genotype for this locus, record the allele count,
        			* and move to the next locus.
        			*/
        			if(two_times_larger){

        			    allele_pair = std::to_string(frequent_variant);
        			    allele_pair.append("/");
        			    allele_pair.append(std::to_string(frequent_variant));
        				MicroOutput << chr << "," << lookup_locus << "," << lookup_table[row][2] << "," << allele_pair << std::endl;
        				/*
        				* First, build the allele identifier - a motif followed
        				* by a number of repeats (for example, TT7)
        				*/
    					allele_identifier = lookup_table[row][2];
    					allele_identifier.append(std::__cxx11::to_string(frequent_variant));

    					// Next, create iterators for the locus map and the allele map
    					std::map<std::string, std::map<std::string, int>>::iterator outer_iterator;
                        std::map<std::string, int>::iterator inner_iterator;

                        /*
                        * If the map is not empty, check if the allele-to-be
                        * -recorded has been "found" before.
                        */
                        for(outer_iterator = allele_counts.begin(); outer_iterator != allele_counts.end(); outer_iterator++){
                            // If the current locus has been recorded before, check if the allele for this locus has been recorded as well
                            if(outer_iterator->first == locus_identifier){
    	                        locus_found = true;
                                for(inner_iterator = outer_iterator->second.begin(); inner_iterator != outer_iterator->second.end(); inner_iterator++){
                                    /*
                                    If the current allele HAS been found before,
                                    * increment its count by 2 to record the
                                    * new homozygote
                                    */
                                    if(inner_iterator->first == allele_identifier){
                                        inner_iterator->second = inner_iterator->second + 2;
                                        allele_identifier.clear();
                                        frequent_variant = 0;
                        				two_times_larger = false;
                        				scratch_genotypes.clear();
                        				allele_found = true;
                                        break;
                                    }else{
                                        continue;
                                    }
                                }
                            // If the locus identifier has yet to be found but there still remains other rows, check those rows
                            }else{
                                locus_found = false;
                                continue;
                            }
                            /*
                            * If you have reached the last allele under the
                            * current locus and have NOT found the allele,
                            * add a new entry of the allele to the map and
                            * stop searching.
                            */
                            if((outer_iterator->first == locus_identifier) && (!allele_found)){
                                outer_iterator->second.insert(std::pair<std::string, int>(allele_identifier, 2));
                                allele_identifier.clear();
                                frequent_variant = 0;
                        		two_times_larger = false;
                        		scratch_genotypes.clear();
                                break;
                            }else if(allele_found){
                                allele_found = false;
								allele_identifier.clear();
                                frequent_variant = 0;
                        		two_times_larger = false;
                        		scratch_genotypes.clear();
                                break;
                            }
                        }
                        // UNTESTED:
                        // If the current locus was not found after iterating through every entry of allele_counts, add the current locus and allele as a new entry and move to the next locus
                        if(!locus_found){
                            std::map<std::string, int> temp_map{{allele_identifier,2}};
                            allele_counts.insert(std::pair<std::string, std::map<std::string, int>>(locus_identifier, temp_map));
                            frequent_variant = 0;
                        	two_times_larger = false;
                        	scratch_genotypes.clear();
                            allele_identifier.clear();
                            /*
                			* Clean all relevant objects for the next locus.
                			*/
                			second_copies = 0;
                			largest_copies = 0;
                			current_copies = 0;
                            continue;
                        }else{
                            locus_found = false;
							/*
                			* Clean all relevant objects for the next locus.
                			*/
                			second_copies = 0;
                			largest_copies = 0;
                			current_copies = 0;
                            continue;
                        }
                        /*
                		* Clean all relevant objects for the next locus.
                		*/
            			second_copies = 0;
            			largest_copies = 0;
            			current_copies = 0;
                        continue;
        			/*
        			* Heterozygous Conditions:
        			* If the allele with the greatest number of copies does NOT
        			* have 2x as many copies as any other allele, find the allele
        			* found in the second greatest number of reads and add both
        			* to the output.
        			*/
        			}else{
    					allele_identifier = lookup_table[row][2];
    					allele_identifier.append(std::__cxx11::to_string(frequent_variant));
    					std::map<std::string, std::map<std::string, int>>::iterator outer_iterator;
                        std::map<std::string, int>::iterator inner_iterator;

        				/*
        				* Search the outer map until you find the current
        				* locus. Check if the allele-to-be-recorded has been
        				* "found" before.
        				*/
                        for(outer_iterator = allele_counts.begin(); outer_iterator != allele_counts.end(); outer_iterator++){
                            if(outer_iterator->first == locus_identifier){
                                locus_found = true;
                                for(inner_iterator = outer_iterator->second.begin(); inner_iterator != outer_iterator->second.end(); inner_iterator++){
                                    /*
                                    * If the current allele HAS been found
                                    * before, increment its count by 1 to
                                    * record one of the heterozygote's alleles.
                                    * Do NOT clear the objects, as you still
                                    * must record the second most common allele
                                    */
                                    if(inner_iterator->first == allele_identifier){
                                        inner_iterator->second = inner_iterator->second + 1;
                                        allele_identifier.clear();
                        				allele_found = true;
                                        break;
                                    }else{
                                        continue;
                                    }
                                }
                        // If the locus identifier has yet to be found but there still remains other rows, check those rows
                            }else{
                                locus_found = false;
                                continue;
                            }

                        // If you have reached the last allele under the current locus and have NOT found the allele, add a new entry of the allele to the map and stop searching.
                            if((outer_iterator->first == locus_identifier) && (!allele_found)){
                                outer_iterator->second.insert(std::pair<std::string, int>(allele_identifier, 1));
								allele_identifier.clear();
                                break;
                            }else if(allele_found){
								allele_found = false;
								allele_identifier.clear();
								break;
							}
                        }
                        // UNTESTED:
                        // If the current locus was not found after iterating through every entry of allele_counts, add the current locus and allele as a new entry and move to the next locus
                        if(!locus_found){
                            std::map<std::string, int> temp_map{{allele_identifier,1}};
                            allele_counts.insert(std::pair<std::string, std::map<std::string, int>>(locus_identifier, temp_map));
                            allele_identifier.clear();
                        }else{
                            locus_found = false;
                        }
        				/*
        				* Find the second most frequent variant in the individual's
        				* reads for this locus
        				*/
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
            			    allele_pair.append("/");
            			    allele_pair.append(std::to_string(second_variant));
        					MicroOutput << chr << "," << lookup_locus << "," << lookup_table[row][2] << "," << allele_pair << std::endl;

        					allele_identifier = lookup_table[row][2];
    						allele_identifier.append(std::__cxx11::to_string(second_variant));
    						std::map<std::string, std::map<std::string, int>>::iterator outer_iterator;
                            std::map<std::string, int>::iterator inner_iterator;

            				/*
            				* Search the outer map until you find the current
            				* locus. Check if the allele-to-be-recorded has been
            				* "found" before.
            				*/
                            for(outer_iterator = allele_counts.begin(); outer_iterator != allele_counts.end(); outer_iterator++){
                                if(outer_iterator->first == locus_identifier){
                                    locus_found = true;
                                    for(inner_iterator = outer_iterator->second.begin(); inner_iterator != outer_iterator->second.end(); inner_iterator++){
                                        /*
                                        * If the current allele HAS been found before,
                                        * increment its count by 1 to record one
                                        * of the heterozygote's alleles.
                                        */
                                        if(inner_iterator->first == allele_identifier){
                                            inner_iterator->second = inner_iterator->second + 1;
                                            allele_identifier.clear();
                                            scratch_genotypes.clear();
                            				allele_found = true;
                                            break;
                                        }else{
                                            continue;
                                        }
                                    }
                                }else{
                                    locus_found = false;
                                    continue;
                                }
								// If you have reached the last allele under the current locus and have NOT found the allele, add a new entry of the allele to the map and stop searching.
								if((outer_iterator->first == locus_identifier) && (!allele_found)){
									outer_iterator->second.insert(std::pair<std::string, int>(allele_identifier, 1));
									allele_identifier.clear();
									frequent_variant = 0;
									second_variant = 0;
									break;
								}else if(allele_found){
									allele_found = false;
									allele_identifier.clear();
									frequent_variant = 0;
									second_variant = 0;
									break;
								}
                            }
                            // UNTESTED:
                            // If the current locus was not found after iterating through every entry of allele_counts, add the current locus and allele as a new entry and move to the next locus
                            if(!locus_found){
                                std::map<std::string, int> temp_map{{allele_identifier,1}};
                                allele_counts.insert(std::pair<std::string, std::map<std::string, int>>(locus_identifier, temp_map));
                                second_variant = 0;
                            	scratch_genotypes.clear();
                                allele_identifier.clear();
                                /*
                    			* Clean all relevant objects for the next locus.
                    			*/
                    			second_copies = 0;
                    			largest_copies = 0;
                    			current_copies = 0;
                                continue;
                            }else{
                                locus_found = false;
								/*
                    			* Clean all relevant objects for the next locus.
                    			*/
                    			second_copies = 0;
                    			largest_copies = 0;
                    			current_copies = 0;
                                continue;
                            }
        					/*
                    		* Clean all relevant objects for the next locus.
                    		*/
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
    		/*
    		* If the current scaffold is not one of the desired scaffolds, skip
    		* the locus
    		*/
    		}else{
        	    continue;
        	}
    	}
    	MicroOutput.close();
    }

    // You now have a map with all alleles and the associated counts for each locus in the lookup table.

    // Open the table of highly polymorphic SNPs and Indels.
    Rcpp::CharacterMatrix SI_table = PollySI(freq_name, desired_scaffolds, SIthreshold);
    std::ofstream SIOutput;
    SIOutput.open(SIOutputName);
    SIOutput << "CHROM,POS,N_ALLELES,N_CHR,{ALLELE:FREQ}\n";
	  for(int i = 0; i < SI_table.nrow(); i++){
	    for(int j = 0; j < SI_table.ncol(); j++){
	        /*
	        * If you've reached the last entry in the row, output a line end
	        * instead of a comma
	        */
	        if(j == (SI_table.ncol() - 1)){
	            SIOutput << SI_table(i, j) << std::endl;
	        }else{
	            SIOutput << SI_table(i, j) << ",";
	        }
	    }
	}
	SIOutput.close();
    /*
    * Sort through the map of official allele counts and retain only loci with
    * >= <MicroThreshold> alleles, then find high-polymorphism SNPs/Indels
    * occurring within 100bp of these loci.
    */

    // These iterators peg each microsat LOCUS.
    std::map<std::string, std::map<std::string, int>>::iterator microsat1;
    std::map<std::string, std::map<std::string, int>>::iterator microsat2;

    // These iterators peg each microsat ALLELE
    std::map<std::string, int>::iterator allele1;
    std::map<std::string, int>::iterator allele2;

	int num_unique_alleles1; // number of unique alleles for a locus
	int num_unique_alleles2;
	std::string locus1;
	std::string locus2;

	/*
	* This double will hold the total number of alleles which were detected
	* across the sample for a single region.
	*/
	double total_alleles;

	// This string will temporarily house all alleles and their associated frequencies for export to the final file.
	std::string allele_freqs;
	// This string will temporarily house all loci occurring within 100bp of each other, separated by commas.
	std::string loci;

	// Open the file into which the final output will be read.
	std::ofstream FinalOutput;
    FinalOutput.open(FinalOutputName);
    FinalOutput << "Scaffold,Total Loci (Count),Total Alleles (Count),SNP Count,Indel Count,Dinucleotide Allele Count,Trinucleotide Allele Count,Tetranucleotide Allele Count,Pentanucleotide (and Above) Allele Count,Loci (Positions),Allele:Frequency\n";

	// Open file into which polymorphic micros will be read
	std::ofstream PollyMicroOutput;
    PollyMicroOutput.open(PollyMicroOutputName);
    PollyMicroOutput << "Scaffold,Locus,Motif,Allele:Frequency\n";

	// Objects to be used repeatedly in following loops
    std::string motif_repeat;
    std::string motif;
    int repeat;

	// Counters for summary stats
	int di = 0;
	int tri = 0;
	int tetra = 0;
	int penta_plus = 0;
	int SNP_count = 0;
	int indel_count = 0;
	int loci_count;
	int allele_count = 0;

	// Unsigned integers to enable comparison to size integers in RCpp
	unsigned int temp_UI = 0;
	unsigned int temp_UI2 = 0;

    for(microsat1 = allele_counts.begin(); microsat1 != allele_counts.end(); microsat1++){
        num_unique_alleles1 = 0;
        /*
        * Upon entering a central microsat, create a new map of ranges of longest
        * microsats within <region_size> of the central microsat so you can
        * later remove any SNPs/Indels which overlap with one or more of the
        * microsats in the reigon.
        */
        std::map<int, int> loci_range;
        num_unique_alleles1 = microsat1->second.size();
    	// Find the NUMBER of alleles at the current microsat locus.

    	/*
    	* If the allele count of the current microsat is greater than
    	* <MicroThreshold>, then it is highly polymorphic. Calculate the
    	* frequency of each of its alleles and output the locus and its
    	* associated alleles.
    	*/
    	if(num_unique_alleles1 >= MicroThreshold){
    		/*
    		* First, isolate the current locus' chromosome. This chromosome
    		* will be outputted in the chromosome column of the final output
    		* file.
    		*/
    		chr = microsat1->first.substr(0, microsat1->first.find(percent));
    		/*
    		* Next, isolate the current locus and add it to the running
    		* string of loci to be outputted. Set the total loci count to 1.
    		*/
    		locus1 = microsat1->first.substr((microsat1->first.find(percent) + percent.length()), (microsat1->first.length() - (microsat1->first.find(percent) + percent.length())));
    		loci = locus1;
    		loci_count = 1;
    		/*
    		* Add the chromosome and locus to the polymorphic micro output.
    		*/
    		PollyMicroOutput << chr << "," << locus1 << ",";
        	/*
        	* Zero out the old total_alleles double and iterate through every
        	* allele in the inner map of allele_counts to add the count of alleles
        	* to total_alleles.
        	* This adjusts allele frequencies at this locus for the fact that
        	* some reads of the locus will not be genotyped.
        	*/
        	total_alleles = 0.0;
        	std::map<std::string, int>::iterator a;
        	for(a = microsat1->second.begin(); a != microsat1->second.end(); a++){
        	    // Convert integer allele count to double
        	    total_alleles += 1.0*(a->second);
        	    continue;
        	}

    		/*
    		* Start the running string of alleles and frequencies
    		* The end goal is a string like this:
    		* [allele1.1:freq1.1/allele1.2:freq1.2], [allele2.1:freq2.1/allele2.2:freq2.2]
    		*/
    		allele_freqs = open_bracket;
    		for(allele1 = microsat1->second.begin(); allele1 != microsat1->second.end(); allele1++){
    			motif_repeat = allele1->first;
    			/*
    			* Extract the motif so you can increment the appropriate allele
    			* count.
    			* If the allele has a two-digit repeat number, output motif as
    			* all units from the first base pair to 2 characters before the
    			* end. Furthermore, check if the current allele is LONGER than
    			* the preceding longest allele FOR THIS MICROSAT. If so, store
    			* it as the longest allele for this locus.
				*/
				if((motif_repeat.substr((motif_repeat.length() - 2), 1) != "N") && (motif_repeat.substr((motif_repeat.length() - 2), 1) != "A") && (motif_repeat.substr((motif_repeat.length() - 2), 1) != "G") && (motif_repeat.substr((motif_repeat.length() - 2), 1) != "C") && (motif_repeat.substr((motif_repeat.length() - 2), 1) != "T")){
				    motif = motif_repeat.substr(0, (motif_repeat.length() - 2));
					repeat = stoi(motif_repeat.substr(motif_repeat.length() - 2, 2));
					/*
					* If this is the first allele, store the current microsat
					* locus as a new entry in the map of loci ranges AND add the
					* motif to the motif column of the polymorphic output
					*/
					if(allele_freqs == open_bracket){
						loci_range.insert(std::pair<int, int>(stoi(locus1), repeat*motif.size() + stoi(locus1)));
						PollyMicroOutput << motif << "," << repeat;
					/*
					* If the current microsat locus HAS been stored in the map
					* already, check if end locus of the current repeat number
					* is greater than the end locus of the previous longest repeat.
					* If so, store the current repeat's range as the largest
					* range for the microsat locus.
					* Also, this means the motif has alreadt been stored in the
					* polymorphic microsat output. Store the new repeat preceded
					* by a SPACE to separate it from other repeat:freqs.
					*/
					// Create unsigned integer so later comparison is not invalid in RCpp
					temp_UI = loci_range.find(stoi(locus1))->second;
					}else if(((repeat*motif.size()) + stoi(locus1)) > temp_UI){
						loci_range.find(stoi(locus1))->second = ((repeat*motif.size()) + stoi(locus1));
						PollyMicroOutput << space << repeat;
					}
					// Increment the appropriate motif size count
					if(motif.size() == 2){
						di++;
						allele_count++;
					}else if(motif.size() == 3){
						tri++;
						allele_count++;
					}else if(motif.size() == 4){
						tetra++;
						allele_count++;
					}else if(motif.size() >= 5){
						penta_plus++;
						allele_count++;
					}
				/*
				* If instead the repeat number has one digit, output the motif
				* as the characters from the first character to 1 character
				* before the end.
				* Also, if the current repeat number is longer than the reigning
				* longest allele, output it as the longest allele
				*/
				}else if(((motif_repeat.substr((motif_repeat.length() - 1), 1) != "N") && (motif_repeat.substr((motif_repeat.length() - 1), 1) != "A") && (motif_repeat.substr((motif_repeat.length() - 1), 1) != "G") && (motif_repeat.substr((motif_repeat.length() - 1), 1) != "C") && (motif_repeat.substr((motif_repeat.length() - 1), 1) != "T")) && ((motif_repeat.substr((motif_repeat.length() - 2), 1) == "N") || (motif_repeat.substr((motif_repeat.length() - 2), 1) == "A") || (motif_repeat.substr((motif_repeat.length() - 2), 1) == "G") || (motif_repeat.substr((motif_repeat.length() - 2), 1) == "C") || (motif_repeat.substr((motif_repeat.length() - 2), 1) == "T"))){
        			motif = motif_repeat.substr(0, (motif_repeat.length() - 1));
					repeat = stoi(motif_repeat.substr((motif_repeat.length() - 1), 1));
					/*
					* If this is the first allele, store the current microsat locus as a new entry in the map of loci ranges
					*/
					if(allele_freqs ==open_bracket){
						loci_range.insert(std::pair<int, int>(stoi(locus1), repeat*motif.size() + stoi(locus1)));
					/*
					* If the current microsat locus HAS been stored in the map
					* already, check if end locus of the current repeat number
					* is greater than the end locus of the previous longest repeat
					*/
					}else if(((repeat*motif.size()) + stoi(locus1)) > temp_UI){
						loci_range.find(stoi(locus1))->second = ((repeat*motif.size()) + stoi(locus1));
					}
					// Increment the appropriate motif size count
					if(motif.size() == 2){
						di++;
						allele_count++;
					}else if(motif.size() == 3){
						tri++;
						allele_count++;
					}else if(motif.size() == 4){
						tetra++;
						allele_count++;
					}else if(motif.size() >= 5){
						penta_plus++;
						allele_count++;
					}
				}
				// Now add the motif and repeat to the allele_freqs string
    			allele_freqs.append(motif_repeat);
    			allele_freqs.append(colon);
    			/*
    			* Now, calculate the allele frequency for the current allele and
    			* add it to the string AND to the polymorphic micro output
    			*/
    			allele_freqs.append(std::__cxx11::to_string((allele1->second/total_alleles)));
    			allele_freqs.append("/");
    			PollyMicroOutput << colon << (allele1->second/total_alleles);
    		}

    		/*
    		* Once you've searched each allele, remove the extra / and cap the
    		* allele_freq string with a closing bracket and a comma.
    		* Then, clear the longest allele.
    		* Add an endline to the Polly micro output
    		*/
    		allele_freqs.pop_back();
    		allele_freqs.append("], ");
    		PollyMicroOutput << std::endl;

    		/*
    		* Now find all microsats occurring within <region_size> of the
    		* current microsat.
    		*/
    		for(microsat2 = allele_counts.begin(); microsat2 != allele_counts.end(); microsat2++){
    			num_unique_alleles2 = microsat2 -> second.size();
    			/*
    			* Check that microsat2...
    			*	1. does NOT match microsat1
    			*	2. has at least MicroThreshold unique alleles and
    			*	3. occurs within <region_size>bp of microsat1.
    			* If all three conditions are met, record the microsat's alleles,
    			* frequencies, and locus.
    			* All polymorphic micros have already been recorded in the first
    			* for loop, so don't worry about recording the here.
    			*/
    			locus2 = microsat2->first.substr((microsat2->first.find(percent) + percent.length()), (microsat2->first.length() - (microsat2->first.find(percent) + percent.length())));
    			if((microsat1->first != microsat2->first) && (num_unique_alleles2 >= MicroThreshold) && (abs(stoi(locus1) - stoi(locus2)) < region_size/2)){
    				// Record the microsat2 locus.
    				loci.append(", ");
    				loci.append(locus2);
    				loci_count++;
    				// Start the new box of alleles
    				allele_freqs.append(open_bracket);
    				for(allele2 = microsat2->second.begin(); allele2 != microsat2->second.end(); allele2++){
    				    motif_repeat = allele2->first;
    					/*
            			* Extract the motif so you can increment the appropriate allele
            			* count.
            			* If the allele has a two-digit repeat number, output motif as
            			* all units from the first base pair to 2 characters before the
            			* end. Furthermore, check if the current allele is LONGER than
            			* the preceding longest allele FOR THIS MICROSAT. If so, store
            			* it as the longest allele for this locus.
        				*/
        				if((motif_repeat.substr((motif_repeat.length() - 2), 1) != "N") && (motif_repeat.substr((motif_repeat.length() - 2), 1) != "A") && (motif_repeat.substr((motif_repeat.length() - 2), 1) != "G") && (motif_repeat.substr((motif_repeat.length() - 2), 1) != "C") && (motif_repeat.substr((motif_repeat.length() - 2), 1) != "T")){
        				    motif = motif_repeat.substr(0, (motif_repeat.length() - 2));
        					repeat = stoi(motif_repeat.substr(motif_repeat.length() - 2, 2));
        					/*
        					* If this is the first allele of microsat2, store
        					* the current microsat locus as a new entry in the
        					* map of loci ranges
        					*/
        					if(std::__cxx11::to_string(allele_freqs[allele_freqs.length() - 1])== open_bracket){
        						loci_range.insert(std::pair<int, int>(stoi(locus2), repeat*motif.size() + stoi(locus2)));
        					/*
        					* If the current microsat locus HAS been stored in the map
        					* already, check if end locus of the current repeat number
        					* is greater than the end locus of the previous longest repeat.
        					* If so, store the current repeat's range as the largest
        					* range for the microsat locus.
        					*/
							    temp_UI2 = loci_range.find(stoi(locus2))->second;
        					}else if(((repeat*motif.size()) + stoi(locus2)) > temp_UI2){
        						loci_range.find(stoi(locus2))->second = ((repeat*motif.size()) + stoi(locus2));
        					}
        					// Increment the appropriate motif size count
        					if(motif.size() == 2){
        						di++;
        						allele_count++;
        					}else if(motif.size() == 3){
        						tri++;
        						allele_count++;
        					}else if(motif.size() == 4){
        						tetra++;
        						allele_count++;
        					}else if(motif.size() >= 5){
        						penta_plus++;
        						allele_count++;
        					}
        				/*
        				* If instead the repeat number has one digit, output the motif
        				* as the characters from the first character to 1 character
        				* before the end.
        				* Also, if the current repeat number is longer than the reigning
        				* longest allele, output it as the longest allele
        				*/
        				}else if(((motif_repeat.substr((motif_repeat.length() - 1), 1) != "N") && (motif_repeat.substr((motif_repeat.length() - 1), 1) != "A") && (motif_repeat.substr((motif_repeat.length() - 1), 1) != "G") && (motif_repeat.substr((motif_repeat.length() - 1), 1) != "C") && (motif_repeat.substr((motif_repeat.length() - 1), 1) != "T")) && ((motif_repeat.substr((motif_repeat.length() - 2), 1) == "N") || (motif_repeat.substr((motif_repeat.length() - 2), 1) == "A") || (motif_repeat.substr((motif_repeat.length() - 2), 1) == "G") || (motif_repeat.substr((motif_repeat.length() - 2), 1) == "C") || (motif_repeat.substr((motif_repeat.length() - 2), 1) == "T"))){
        				    motif = motif_repeat.substr(0, (motif_repeat.length() - 1));
        					repeat = stoi(motif_repeat.substr((motif_repeat.length() - 1), 1));
        					/*
        					* If this is the first allele, store the current
        					* microsat locus as a new entry in the map of loci ranges
        					*/
        					if(std::__cxx11::to_string(allele_freqs[allele_freqs.length() - 1]) == open_bracket){
        						loci_range.insert(std::pair<int, int>(stoi(locus2), ((repeat*motif.size()) + stoi(locus2))));
        					/*
        					* If the current microsat locus HAS been stored in the map
        					* already, check if end locus of the current repeat number
        					* is greater than the end locus of the previous longest repeat
        					*/
        					}else if(((repeat*motif.size()) + stoi(locus2)) > temp_UI2){
        						loci_range.find(stoi(locus2))->second = ((repeat*motif.size()) + stoi(locus2));
        					}
        					// Increment the appropriate motif size count
        					if(motif.size() == 2){
        						di++;
        						allele_count++;
        					}else if(motif.size() == 3){
        						tri++;
        						allele_count++;
        					}else if(motif.size() == 4){
        						tetra++;
        						allele_count++;
        					}else if(motif.size() >= 5){
        						penta_plus++;
        						allele_count++;
        					}
        				}
        				// Now add the motif and repeat to the allele_freqs string
    					allele_freqs.append(motif_repeat);
    					allele_freqs.append(colon);
    					// Now, calculate the allele frequency for the current allele and add it to the string.
    					allele_freqs.append(std::__cxx11::to_string((allele2->second/total_alleles)));
    					allele_freqs.append("/");
    					// Increment total allele count by 1
    			        allele_count += 1;
    				}
    				/*
    				* Once you've searched each allele, remove the extra / and
    				* cap the allele_freq string with a closing bracket and a comma.
    				*/
    				allele_freqs.pop_back();
    				allele_freqs.append("], ");

    				/*
    				* If the microsat2 does not occur within <region_size>bp of
    				* microsat1, simply move to the next microsat.
    				*/
    				}else{
    					locus2.clear();
    					continue;
    				}
    			    /*
    			    * Once you've recorded the current microsat2, clear the
    			    * relevant objects and move to the next microsat2.
    			    */
    				locus2.clear();
    				continue;
    			}

    	        PollyMicroOutput.close();
    			/*
    			* Once you've found all microsats within <region_size>bp of
    			* microsat1, find the SNPs and Indels occurring within
    			* <region_size>bp of microsat1
    			*/
    			for(int row = 0; row < SI_table.nrow(); row++){
    				/*
    				* If the SNP/Indel occurs within <region_size>bp of
    				* microsat1, record it, up the loci count by 1, and up the
    				* allele count by 2.
    				* Check the identity of the SNP/Indel and up the appropriate
    				* count by 1.
    				*/
    				/*
    				 *  Store matrix element in string prior to converting to integer or
    				 *  making a comparison.
    				*/
    				SI_element0 = SI_table(row, 0);
    			  SI_element1 = SI_table(row,1);
    				if((chr == SI_element0) && (abs(stoi(locus1) - stoi(SI_element1)) < region_size/2)){
    				    /*
    				    * Before recording an SNP or indel, ensure it does not
    				    * overlap with any pre-recorded microsats in the region.
    				    */
    				    std::map<int, int>::iterator locus_range;
        				for(locus_range = loci_range.begin(); locus_range != loci_range.end(); locus_range++){
    						SI_element1 = SI_table(row, 1);
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
        					allele_freqs.append("[");
        					allele_freqs.append(SI_table(row, 4));
        					allele_freqs.append("/");
        					allele_freqs.append(SI_table(row, 5));
        					allele_freqs.append("], ");
        					loci.append(", ");
        					loci.append(SI_table(row, 1));
        					loci_count++;
        					allele_count += 2;
        					// If both alleles are size 1, the locus is a SNP
        					SI_element4 = SI_table(row, 4);
        					SI_element5 = SI_table(row, 5);
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
    		* If microsat1 doesn't have <MicroThreshold> unique alleles, move to
    		* the next microsat.
    		*/
    		} else{
    			continue;
    		}
    	/*
    	* Once you've found every variant within <region_size>bp of microsat1,
    	* remove the trailing space and comma from the allele frequencies,
    	* output the row to the output file, clear the relevant objects, and
    	* move to the next microsat.
    	*/
    	// Remove space
    	allele_freqs.pop_back();
    	// Remove comma
    	allele_freqs.pop_back();

    	FinalOutput << chr << "," << loci_count << "," << allele_count << "," << SNP_count << "," << indel_count << "," << di << "," << tri << "," << tetra << "," << penta_plus << "," << loci << "," << allele_freqs << std::endl;
    	chr.clear();
    	allele_freqs.clear();
    	loci.clear();
    	locus1.clear();
    	di = 0;
	    tri = 0;
	    tetra = 0;
    	penta_plus = 0;
	    SNP_count = 0;
	    indel_count = 0;
	    loci_count = 0;
	    allele_count = 0;
    	continue;
    	}
    	FinalOutput.close();
	}
