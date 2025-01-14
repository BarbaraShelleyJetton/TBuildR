# here is the code so far, it will get fasta files for you from a species list / gene list and then smashes
# mitochondrial genes together to make a frankenstein mitochondrial sequence. If you have questions or anything
# let me know. I also use some packages for smaller alignments / trimming / functions for random tree building shit. 
# so let me know if you want those or if there is anything else I could do to help. <3

# anyways, here are the libraries used for this

library(tidyverse)
library(janitor)
library(Biostrings)
library(DECIPHER)
library(foreach)
library(doParallel)
library(rentrez)

# first obtain a csv of species you want to use, then clean that file to have only unique species or however else 
# you want that shit sorted. I did that somewhere else and lost it, but it was V easy with janitor

species_list <- read_csv("unique_diptera.csv")

# the gene list is used for genes you want to search for.
# you can search for more than one at a time, I just forgot to test it with mito_smash...Im sure it'll work

gene_list <- c("coi")

# now we search genbank for gene hits for your species

generate_and_perform_search <- function(species_list, gene_list, retmax_value = 10) {
  search_results <- list()
  for (species in species_list$species) {
    for (gene in gene_list) {
      search_term <- paste0(species, "[Organism] AND ", gene, "[ALL]")
      search_result <- entrez_search(db = "nucleotide", term = search_term, retmax = retmax_value, use_history = TRUE)
      
      if (search_result$count > 0) {
        search_results[[paste(species, gene)]] <- search_result
      } else {
        cat("No hits found for: ", search_term, "\n")
      }
    }
  }
  
  return(search_results)
}

search_results <- generate_and_perform_search(species_list, gene_list)

# then we get a summary for the hits we got

get_sum <- function(search_results){
  summaries <- list()
  for (key in names(search_results)){
    ids <- search_results[[key]]$ids
    sum_result <- entrez_summary(db = "nucleotide", id = ids)
    summaries[[key]] <- sum_result
  }
  return(summaries)
}

species_summaries <- get_sum(search_results)

# using the summaries, we can filter the accession numbers. there is a max length to ensure we are not 
# obtaining entire genomes (can be changed to whatever number you want)
# this function also ensures that the gene you are searching for is in the title of the accession (hopefully)

get_filtered_accessions <- function(species_summaries, gene_list, desired_length = NULL) {
  filtered_accessions <- list()
  
  for (key in names(species_summaries)) {
    summaries <- species_summaries[[key]]
    filtered_accession <- NULL
    
    if (is.list(summaries) && length(summaries) == 33) {
      if (!is.null(summaries$slen) && !is.null(summaries$accession) && !is.null(summaries$title)) {
        contains_gene <- any(sapply(gene_list, function(gene) grepl(gene, summaries$title, ignore.case = TRUE)))
        
        if (contains_gene && (is.null(desired_length) || summaries$slen <= desired_length)) {
          filtered_accession <- summaries$accession
        }
      }
    } else if (is.list(summaries)) {
      for (i in seq_along(summaries)) {
        current_hit <- summaries[[i]]
        current_length <- current_hit$slen
        current_title <- current_hit$title
        
        contains_gene <- any(sapply(gene_list, function(gene) grepl(gene, current_title, ignore.case = TRUE)))
        
        if (contains_gene && (is.null(desired_length) || current_length <= desired_length)) {
          filtered_accession <- current_hit$accession
          break
        }
      }
    }
    filtered_accessions[[key]] <- filtered_accession
  }
  
  return(filtered_accessions)
}

filtered_accessions <- get_filtered_accessions(species_summaries, gene_list, desired_length = 18000)

# after filtering, we can now get a fasta file containing the gene info we found 

get_fasta <- function(accession_numbers, output_file) {
  all_recs <- character()
  for (accession_list in accession_numbers) {
    for (accession in accession_list) {
      rec <- entrez_fetch(db = "nuccore", id = accession, rettype = "fasta")
      all_recs <- c(all_recs, rec)
    }
  }
  writeLines(all_recs, con = output_file)
}

get_fasta(filtered_accessions, "coi.fasta")

# I plan on combining these functions so you dont have to go through all those steps eventually, bear with me.

# now is the funky part, I pulled all my genes separate because I hate myself apparently. 
# this function goes through each of my fasta files and creates a list that stores my species with its various genes
# I have not tried it on a fasta file with multiple genes, so if you want to try that go ahead. I am pretty
# sure it will work? Anyways, this is the list mito_smash iterates through to do its thang

# again I need to condense my functions into one, just havent done it. you need to make sure extract_species_name
# and merge_all_seqs are in your functions when running process_fasta_files

# btw this uses parallel processing because it seemed faster to me with the amount of data I am working with. 
# it runs on all cores except 1, if you want that to be lower change the num_cores detectCores() -1 to however cores
# you dont want used. I have a non parallel version as well but I broke it. oops.

# run this

extract_species_name <- function(header, species_names) {
  matches <- character(length(species_names))
  for (i in seq_along(species_names)) {
    pattern <- paste0("\\b", species_names[i], "\\b")
    match <- regmatches(header, regexpr(pattern, header, ignore.case = TRUE))
    if (length(match) > 0) {
      matches[i] <- match
    }
  }
  matches <- matches[matches != ""]
  if (length(matches) > 0) {
    return(matches)
  } else {
    return(NA)
  }
}

# and this

merge_all_seqs <- function(x, y) {
  merged <- x
  for (key in names(y)) {
    if (key %in% names(x)) {
      merged[[key]] <- c(merged[[key]], y[[key]])
    } else {
      merged[[key]] <- y[[key]]
    }
  }
  return(merged)
}

# then this

process_fasta_files <- function(file_paths, species_names) {
  num_cores <- detectCores() - 1 # change this if you want
  
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  clusterExport(cl, c("readDNAStringSet", "extract_species_name", "species_names"))
  
  all_seqs <- foreach(file_path = file_paths, .combine = merge_all_seqs) %dopar% {
    dat <- readDNAStringSet(file_path)
    
    
    species <- character(length(dat))
    for (i in 1:length(dat)) {
      species[i] <- extract_species_name(names(dat)[i], species_names)
    }
    
    
    species <- species[!is.na(species)]
    all_seqs_per_file <- list()
    for (i in seq_along(species)) {
      if (!(species[i] %in% names(all_seqs_per_file))) {
        all_seqs_per_file[[species[i]]] <- list()
      }
      all_seqs_per_file[[species[i]]][[basename(file_path)]] <- dat[[i]]
    }
    return(all_seqs_per_file)
  }
  
  stopCluster(cl)
  
  
  all_seqs_filtered <- all_seqs[!sapply(all_seqs, function(x) all(sapply(x, is.null)))]
  
  return(all_seqs_filtered)
}

file_paths <- c("coi.fasta", "coii.fasta", "cytb.fasta") # put your fasta files that you want to mito_smash in here

species_names <- species_list$species

all_seqs_filtered <- process_fasta_files(file_paths, species_names)

# sometimes the list from process_fasta_files comes back with a lower number, aka lower species count (I think this 
# is becuase there are a. repeats or b. some specie names are not being matched from the accessions)

# at this point you can look at all_seqs_filtered and determine if there are any species you want removed, or dont, you do you. 

# time to mito_smash!!!!!!!! for this part you need a reference sequence (I am using a complete mitochondrial genome 
# from Drosophila melanogaster)

# this also uses parallel processing. change the detect cores if you like

mito_smash <- function(all_seqs_filtered, ref_mito, output_file) {
  
  file_conn <- file(output_file, "w")
  
  for (species_name in names(all_seqs_filtered)) {
    species_data <- all_seqs_filtered[[species_name]]
    
    gene_seqs <- list()
    
    for (gene_name in names(species_data)) {
      if (!is.null(species_data[[gene_name]])) {
        gene_seqs[[gene_name]] <- species_data[[gene_name]]
      }
    }
    
    seqs <- c(DNAStringSet(unlist(gene_seqs)), ref_mito)
    
    aln <- AlignSeqs(seqs, processors = (parallel::detectCores() - 1)) # change here if you want
    
    cons_seq <- (ConsensusSequence(aln[-length(aln)]))
    
    total_chars <- nchar(cons_seq)
    
    
    full_chunks <- floor(total_chars / 80)
    
    
    cons_seq_chunks <- substring(cons_seq, seq(1, by = 80, length.out = full_chunks),
                                 seq(80, by = 80, length.out = full_chunks))
    
    
    if (total_chars %% 80 != 0) {
      remaining_chars <- total_chars %% 80
      last_chunk <- substring(cons_seq, first = total_chars - remaining_chars + 1, last = total_chars)
      cons_seq_chunks <- c(cons_seq_chunks, last_chunk)
    }
    
    cons_seq_chunks_new <- paste(cons_seq_chunks, collapse = "\n")
    line <- paste0(">", species_name, "\n", cons_seq_chunks_new, "\n")
    writeLines(line, file_conn)
  }
  
  
  close(file_conn)
}

ref_mito <- readDNAStringSet("reference_mitochondrion_sequence.fasta")

mito_smash(all_seqs_filtered, ref_mito, "mito_smash.fasta")

# once you have your mito_smash fasta, you want to run it through MAFFT. the previous alignment was aligning the individual
# species genes to the reference, so MAFFT will align all species together.

