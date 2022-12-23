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

#' Converts MicroGenotyper() outputs to GenePop format
#' @export
PollyPop <- function(MicroGenotyperInputNames, PollyMicrosInputNames, 
                     output_name, GenePop_header = ""){
  # Function 1: Take a name in "pop/csv" format...
  GenePopIndv <- function(pop_indv){
    print(paste(c("------------ GenePop Conversion: Pop = ", pop_indv[1],
                  " Individual = ", pop_indv[2], " ------------"), collapse = ""))
    #  1. Define Function 2
    #  Function 2: Take a MicroID as input.
    AlleleConverter <- function(id){
      print(paste(c("ID: ", which(colnames(MegaDF) == id), "/", ncol(MegaDF)), 
                  collapse = ""))
      # Use the alleles in MicroGenotypes to convert the alleles
      # in "df" to GenePop format and record the individual's genotype in temp
      if(id %in% df$MicroID){
        cur_alleles <- strsplit(df$Genotype[df$MicroID == id], "/")[[1]]
        all_alleles <- sort(levels(as.factor(
          strsplit(MicroGenotypes$Alleles[MicroGenotypes$MicroID == id], "/")[[1]]
        )))
        genepop_alleles <- sprintf("%02d", 1:length(all_alleles))
        cur_alleles <- paste(
          c(genepop_alleles[all_alleles == cur_alleles[1]],
            genepop_alleles[all_alleles == cur_alleles[2]]), collapse = "")
      }else{
        # If the current id is NOT present in the current individual, record this 
        # individual's genotype as 00
        cur_alleles <- sprintf("%04d", 0000)
      }
      return(cur_alleles)
    }
    
    # 2. Read CSV (Don't read in the motif column)
    df <- read.csv(pop_indv[2])[,-3]
    
    # 3. Convert each scaffold to the numeric value
    for(i in 1:nrow(scaff_key)){
      df$Scaffold[df$Scaffold == scaff_key$original[i]] <- scaff_key$code[i]
    }
    
    # 4. Convert each scaffold + locus to MicroID
    df$MicroID <- apply(df[,1:2], MARGIN = 1, FUN = paste, collapse = "-")
    df <- df[,-c(1,2)]
    
    # 5. Apply function 2 to each MicroID and store output as a row of a 
    # dataframe with the same column names as "MegaDF"
    temp <- vapply(sort(MicroGenotypes$MicroID), FUN = AlleleConverter, 
                   FUN.VALUE = vector(mode = "character", length = 1))
    
    return(c(pop_indv[2], temp))
  }
  
  # Check if the number of populations is the same for MicroGenotyper and 
  # PollyMicros. If not, throw an error.
  if(length(MicroGenotyperInputNames) != length(PollyMicrosInputNames)){
    return("ERROR: The number of populations described in the MicroGenotyper() 
           outputs is not the same as the number of populations described in the 
           PollyMicros() outputs. Please make sure the list of vectors inputted 
           for MicroGenotyperInputNames is the same length as the vector 
           inputted for PollyMicrosInputNames.")
  }
  
  
  # Step 1: Find all unique microsatellite IDs in the polymorphism files
  #       1.1) Initialize a dataframe of each polymorphic MicroID and all alleles
  MicroGenotypes <- data.frame(matrix(ncol = 3, nrow = 0))
  colnames(MicroGenotypes) <- c("Scaffold", "Locus", "Alleles")
  
  #       1.2) Create a key of all scaffolds and numeric values
  scaffs <- c()
  for(pop in 1:length(PollyMicrosInputNames)){
    scaffs <- append(scaffs, levels(as.factor(read.csv(PollyMicrosInputNames[[pop]])[,1])))
  }
  scaff_key <- data.frame(original = levels(as.factor(scaffs)), 
                          code = 1:length(levels(as.factor(scaffs))))
  
  #       1.3) Convert each scaffold name in the polymorphism files to numerics
  for(pop in 1:length(PollyMicrosInputNames)){
    MicroGenotypes <- rbind(MicroGenotypes, 
                            read.csv(PollyMicrosInputNames[[pop]])[,c(1,2,4)])
  }
  # Replace each scaffold name with key
  for(i in 1:nrow(scaff_key)){
    MicroGenotypes$Scaffold[MicroGenotypes$Scaffold == scaff_key$original[i]] <- 
      scaff_key$code[i]
  }
  #       1.4) Merge each scaffold name and corresponding locus into a MicroID
  MicroGenotypes$MicroID <- apply(MicroGenotypes[,1:2], MARGIN = 1, FUN = paste, 
                                  collapse = "-")
  MicroGenotypes <- MicroGenotypes[,-c(1,2)]
  #       1.5) Merge alleles of the same microgenotype but from different pops
  #            into single columns
  # 1.5.1) Create vector describing the names of the duplicated MicroIDs (dup_ids)
  dup_ids <- levels(as.factor(MicroGenotypes$MicroID[duplicated(MicroGenotypes$MicroID)]))
  
  # 1.5.2) Create temporary dataframe retaining all MicroIDs which are unique
  #        to a single population
  temp <- MicroGenotypes[!(MicroGenotypes$MicroID %in% dup_ids),]
  
  # 1.5.3) Combine alleles from different populations of a duplicated MicroID 
  #        and add to a single row describing that MicroID, then append the
  #        row to temp
  for(id in 1:length(dup_ids)){
    print("----------------- Collecting Alleles of IDs Found Across Multiple Populations -----------------")
    print(paste(c(id, "/", length(dup_ids)), collapse = ""))
    temp <- rbind(temp, data.frame(
      Alleles = paste(c(MicroGenotypes$Alleles[MicroGenotypes$MicroID == dup_ids[id]]), 
                      collapse = "/"),
      MicroID = dup_ids[id]))
  }
  MicroGenotypes <- temp
  print("")
  
  
  
  # Step 2: Initialize a dataframe (MegaDF) storing...
  #       1. Each population-individual as a row name
  #       2. Each microsatellite ID as a column name (SORTED)
  MegaDF <- data.frame(matrix(
    nrow = sum(unlist(lapply(MicroGenotyperInputNames, length)) + 1), 
    ncol = nrow(MicroGenotypes) + 1))
  colnames(MegaDF) <- c("", sort(MicroGenotypes$MicroID))
  
  # Step 3: Convert each individual's genotype to GenePop format
  #      3.1) For each individual MicroGenotyper() file, convert name to 
  #           pop/indv format
  #      3.2) Feed name into "GenePopIndv" and append output as a row to MegaDF
  cur_row <- 1
  for(p in 1:length(MicroGenotyperInputNames)){
    # Add a "pop" prior to each population
    MegaDF[cur_row, ] <- c("pop", rep("", ncol(MegaDF) - 1))
    cur_row <- cur_row + 1
    for(i in 1:length(MicroGenotyperInputNames[[p]])){
      cur_popindv <- c(p, MicroGenotyperInputNames[[p]][i])
      MegaDF[cur_row, ] <- GenePopIndv(cur_popindv)
      cur_row <- cur_row + 1
    }
  }
  
  # Step 4: Convert MegaDF to GenePop format
  OutputDF <- data.frame(matrix(ncol = ncol(MegaDF), 
                                nrow = nrow(MegaDF) + 2))
  OutputDF[1,] <- c(GenePop_header, rep("", ncol(MegaDF)-1))
  OutputDF[2,] <- c(colnames(MegaDF)[1], unlist(lapply(colnames(MegaDF)[-1], paste, ",")))
  OutputDF[3:nrow(OutputDF), ] <- MegaDF
  
  # Step 5: Output GenePop-formatted file AND scaffold key
  write.table(OutputDF, file = output_name, 
              sep = "\t", row.names = F, quote = F, col.names = F)
  write.csv(scaff_key, file = "Scaffold Key.csv", row.names = F)
  return(NULL)
}