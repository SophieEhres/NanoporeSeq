## list of functions for consensus annotation

`%notin%` <- function(x, table) {
  !(x %in% table)
}

# annotate stringtie transcripts with either gene_id/ref_gene_id if it is ENSG
concat_gene_and_ref_gene <- function(annotation_gr) {
  annotation_gr$unique_gene_id <- ifelse( # make new unique_gene_id column
    grepl("ENSG", annotation_gr$gene_id), # if it has ENSG, take gene_id
    annotation_gr$gene_id,
    ifelse(
      grepl("MCF10a", annotation_gr$gene_id) & grepl("ENSG", annotation_gr$ref_gene_id), # if gene_id has MCF10a and ref_gene_id has ENSG, take ref_gene_id
      annotation_gr$ref_gene_id,
      annotation_gr$gene_id # else take gene_id (novel genes)
    )
  )

  n_genes <- length(unique(annotation_gr$gene_id))
  n_gene_id_ENSG <- length(unique(annotation_gr$gene_id[grepl("ENSG", annotation_gr$gene_id)]))
  n_gene_id_MCF10a <- length(unique(annotation_gr$gene_id[grepl("MCF10a", annotation_gr$gene_id) &
    grepl("ENSG", annotation_gr$ref_gene_id)]))
  n_gene_id_none <- length(unique(annotation_gr$gene_id[!grepl("ENSG", x = annotation_gr$gene_id) &
    !grepl("ENSG", x = annotation_gr$ref_gene_id)]))

  cat("Total genes:", n_genes, "\n")
  cat("ENSG gene_ids: ", n_gene_id_ENSG, "\n")
  cat("MCF10a gene_ids with ENSG ref_gene_id: ", n_gene_id_MCF10a, "\n")
  cat("MCF10a gene_ids with no ENSG: ", n_gene_id_none, "\n")
  return(annotation_gr)
}

# flatten exons into genes from all transcripts
make_gene_from_exons <- function(exon_gr) {
  require(GenomicRanges)
  require(GenomicFeatures)

  # split exons by unique_gene_id
  gene_ranges_list <- split(exon_gr, exon_gr$unique_gene_id)

  # reduce to get gene features
  gene_ranges <- GenomicRanges::reduce(gene_ranges_list)

  # Flatten into one range per gene (min start, max end)
  gene_flat <- lapply(gene_ranges, function(gr) {
    GRanges(
      seqnames = seqnames(gr)[1],
      ranges = IRanges(start = min(start(gr)), end = max(end(gr))),
      strand = strand(gr)[1]
    )
  })


  gene_flat_combined <- do.call(c, unname(gene_flat))
  gene_flat_combined$unique_gene_id <- names(gene_flat)
  gene_flat_combined$type <- "gene"

  n_genes <- length(gene_flat_combined)

  cat("Gene features generated: ", n_genes, "\n")

  return(gene_flat_combined)
}

# add metadata from gencode to a gene_gr
add_metadata_from_gff_genes <- function(gene_gr, gff_gene_gr, columns_to_add) {
  # match IDs (gff) to unique_gene_id (gtf)
  gene_gr$unique_gene_id <- as.character(gene_gr$unique_gene_id)
  gff_gene_gr$ID <- as.character(gff_gene_gr$ID)
  match_idx <- match(gene_gr$unique_gene_id, gff_gene_gr$ID)

  # add metadata by match ID
  for (col in columns_to_add) {
    mcols(gene_gr)[[col]] <- mcols(gff_gene_gr)[[col]][match_idx]
  }

  mcols(gene_gr)$type <- "gene"

  gene_gr$biotype[is.na(gene_gr$biotype)] <- "novel_gene"
  gene_gr <- sort(gene_gr)

  mcols(gene_gr)$type <- "gene"

  n_duplicates <- length(gene_gr[duplicated(gene_gr$unique_gene_id), ]) ## check only one per gene

  if (n_duplicates > 0) {
    warning(paste("Duplicates made:", n_duplicates, "\n"))
  } else {
    cat("No duplicates made", "\n")
  }
  
  gene_number <- length(gene_gr)
  gene_number_with_metadata_from_gff <- length(gene_gr[gene_gr$biotype != "novel_gene"])
  gene_number_novel <- length(gene_gr[gene_gr$biotype == "novel_gene"])
  cat(paste(gene_number_with_metadata_from_gff, "genes out of", gene_number, "have data from gff", "\n"))
  cat(paste(gene_number_novel, "genes out of", gene_number, "are not found in gff", "\n"))

  return(gene_gr)
}

# add metadata from gencode to known and novel transcripts
add_metadata_from_gff_transcripts <- function(gtf, gencode_transcript_gff, columns_to_adds) {
  # filter for transcript_ids in gencode_gff
  gtf_known <- gtf[gtf$transcript_id %in% gencode_transcript_gff$transcript_id]
  gtf_novel <- gtf[gtf$transcript_id %notin% gencode_transcript_gff$transcript_id]

  # Match transcript_ids
  match_idx <- match(gtf_known$transcript_id, gencode_transcript_gff$transcript_id)


  # Loop to add metadata to all rows based on transcript_id match
  for (col in columns_to_add) {
    mcols(gtf_known)[[col]] <- mcols(gencode_transcript_gff)[[col]][match_idx]
  }

  # add the same data for novel transcripts

  gtf_novel$biotype <- "novel_RNA"
  gtf_novel$description <- NA
  gtf_novel$Name <- NA
  gtf_novel$version <- NA
  gtf_novel$Parent <- ifelse(
    test = grepl(pattern = "ENSG", x = gtf_novel$unique_gene_id),
    yes = gtf_novel$unique_gene_id,
    no = paste0("gene:", gtf_novel$unique_gene_id)
  )

  gtf_final <- c(gtf_known, gtf_novel)


  known <- length(gtf_known[gtf_known$type == "transcript"])
  annotated_known <- length(gtf_known$biotype[complete.cases(gtf_known$biotype) & gtf_known$type == "transcript"])
  novel <- length(gtf_novel[gtf_novel$type == "transcript"])

  cat("Known transcripts: ", known, "\n")
  cat("annotated known transcripts: ", annotated_known, "\n")
  cat("Novel transcripts: ", novel, "\n")
  return(gtf_final)
}

# get start and stop features from cds genomic ranges
get_start_and_stop <- function(cds_gr) {
  cds_split <- split(cds_gr, cds_gr$transcript_id) # split by transcript_id

  cds_flat <- lapply(cds_split, function(gr) { # flatten to get one feature per CDS
    GRanges(
      seqnames = seqnames(gr)[1],
      ranges = IRanges(start = min(start(gr)), end = max(end(gr))),
      strand = strand(gr)[1],
      mcols = mcols(gr[1])
    )
  })

  cds_flat <- do.call(c, unname(cds_flat)) # get a granges object from a grangeslist

  number_CDS <- length(cds_flat)
  number_gencode_transcripts <- length(unique(cds_gr$transcript_id))

  cat(paste("There are", number_CDS, "CDS annotations out of", number_gencode_transcripts, "transcripts"))

  starts_gr <- GRanges(
    seqnames = seqnames(cds_flat),
    ranges = IRanges(start = start(cds_flat), width = 1),
    strand = strand(cds_flat)
  )
  
  # reassign correctly the mcols()
  mcols(starts_gr) <- NULL
  mcols(starts_gr) <- mcols(cds_flat)
  
  starts_gr$type <- "start_codon"


  # same for stop codon

  stops_gr <- GRanges(
    seqnames = seqnames(cds_flat),
    ranges = IRanges(start = end(cds_flat), width = 1),
    strand = strand(cds_flat)
  )
  
  mcols(stops_gr) <- NULL
  mcols(stops_gr) <- mcols(cds_flat)
  stops_gr$type <- "stop_codon"

  start_stop_gr <- c(starts_gr, stops_gr)

  return(start_stop_gr)
}

# get fasta of transcript granges
get_fasta_of_transcripts <- function(transcript_gr, genome) {
  if (all(grepl(seqlevels(genome), pattern = "chr")) & all(grepl(seqlevels(transcript_gr), pattern = "chr"))) {
    cat("all seqlevels are the same - proceed \n")

    # if chr in genome // no chr in gtf
  } else if (all(grepl(seqlevels(genome), pattern = "chr")) & all(!grepl(seqlevels(transcript_gr), pattern = "chr"))) {
    seqlevels(genome) <- sapply(str_split(seqlevels(genome), pattern = "r"), "[", 2)
    cat("changing seqlevels to fit gtf: chr -> no chr \n")

    # if no chr in genome // chr in gtf
  } else if (all(!grepl(seqlevels(genome), pattern = "chr")) & all(grepl(seqlevels(transcript_gr), pattern = "chr"))) {
    seqlevels(genome) <- paste("chr", seqlevels(genome), sep = "")
    cat("changing seqlevels to fit gtf: no chr -> chr \n")
  }

  # group exons by transcript
  exons <- transcript_gr[transcript_gr$type == "exon"]
  exons_by_tx <- split(exons, exons$transcript_id)

  # get sequences
  transcripts_fasta <- DNAStringSet(lapply(exons_by_tx, function(tx) {
    seqs <- getSeq(genome, tx)
    if (as.character(unique(strand(tx))) == "-") {
      reverseComplement(do.call(xscat, rev(seqs)))
    } else {
      do.call(xscat, seqs)
    }
  }))

  names(transcripts_fasta) <- names(exons_by_tx)
  return(transcripts_fasta)
}

# validate size of predicted ORFs
validate_and_fix_cds <- function(cpc2_coding) {
  # Verify that ORFs are not longer than transcripts
  invalid_orf <- which(cpc2_coding$ORF_start + cpc2_coding$peptide_length * 3 > cpc2_coding$transcript_length)
  if (length(invalid_orf) > 0) {
    cat(paste("adjusting lengths for", length(invalid_orf), "ORFs longer than transcript \n"))
    # Adjust length of peptide
    cpc2_coding$peptide_length[invalid_orf] <- floor((cpc2_coding$transcript_length[invalid_orf] - cpc2_coding$ORF_start[invalid_orf]) / 3)
  }
  
  # recalculate end
  cpc2_coding$ORF_end <- cpc2_coding$ORF_start + cpc2_coding$peptide_length * 3 - 1
  
  return(cpc2_coding)
}

# get cds ranges of predicted ORFs
get_cds_from_cpc2 <- function(cpc2_result) {
  # set colnames
  colnames(cpc2_result) <- c(
    "ID", "transcript_length", "peptide_length",
    "Fickett_score", "PI", "ORF_integrity", "ORF_start",
    "coding_probability", "label"
  )

  # set more stringent threshold at probability 0.7
  cpc2_result$label[cpc2_result$coding_probability < 0.7] <- "noncoding"

  # keep only predicted coding
  cpc2_coding <- cpc2_result[cpc2_result$coding_probability >= 0.7, ]

  # add ORF_end to match to transcript
  cpc2_coding$ORF_end <- c(cpc2_coding$ORF_start + cpc2_coding$peptide_length * 3)

  cpc2_coding_validated <- validate_and_fix_cds(cpc2_coding)

  #### get cds granges ####
  cds_novel <- GRanges(
    seqnames = Rle(cpc2_coding_validated$ID),
    ranges = IRanges(start = cpc2_coding_validated$ORF_start, end = cpc2_coding_validated$ORF_end),
    strand = "*" # Strand is derived from exons later
  )

  names(cds_novel) <- cpc2_coding_validated$ID

  ntranscripts <- nrow(cpc2_result)
  ncoding <- nrow(cpc2_coding_validated)

  cat("CPC2 predicts ", ncoding, " ORFs from ", ntranscripts, " transcripts (>0.7 prob) \n")

  return(cds_novel)
}

# calculate phase for every cds exon
calculate_cds_phase <- function(cds_gr, exons_by_tx) {
  
  # Group by tx
  tx_ids <- unique(cds_gr$transcript_id)
  
  cds_with_phase_list <- lapply(tx_ids, function(tx_id) {
    # select tx
    cds_tx <- cds_gr[cds_gr$transcript_id == tx_id]
    
    # Get exons of tx
    if (!tx_id %in% names(exons_by_tx)) {
      return(NULL)
    }
    exons_tx <- exons_by_tx[[tx_id]]
    
    # Verify that CDS are in exons
    overlaps <- findOverlaps(cds_tx, exons_tx, type = "within")
    if (length(overlaps) == 0) {
      cat("CDS of", tx_id, "not in exons\n")
      return(NULL)
    }
    
    # Get strand from exons
    strand_tx <- as.character(strand(exons_tx)[1])
    strand(cds_tx) <- strand_tx
    
    # Sort by strand
    if (strand_tx == "+") {
      cds_tx <- sort(cds_tx)
    } else {
      cds_tx <- sort(cds_tx, decreasing = TRUE)
    }
    
    # Calculate phases
    cum_len <- 0
    phase_vec <- integer(length(cds_tx))
    
    for (i in seq_along(cds_tx)) {
      phase_vec[i] <- cum_len %% 3
      cum_len <- cum_len + width(cds_tx[i])
    }
    
    cds_tx$phase <- phase_vec
    cds_tx$type <- "CDS"
    
    return(cds_tx)
  })
  
  # Remove nulls and make GRanges
  cds_with_phase_list <- cds_with_phase_list[!sapply(cds_with_phase_list, is.null)]
  
  if (length(cds_with_phase_list) > 0) {
    return(do.call(c, cds_with_phase_list))
  } else {
    return(GRanges())
  }
}

## map cds to transcript structure
map_cds_to_transcript <- function(cds_gr, novel_gtf) {
  # take only exons of gtf
  novel_exons <- novel_gtf[novel_gtf$type == "exon"]

  # split exons by transcript_id
  novel_exons_by_tx <- split(novel_exons, novel_exons$transcript_id)

  # take only coding
  novel_exons_by_tx_coding <- novel_exons_by_tx[names(novel_exons_by_tx) %in% names(cds_gr)]

  # Map to genome
  cds_mapped <- pmapFromTranscripts(cds_gr, novel_exons_by_tx_coding)
  
  # Flatten mapped CDS
  cds_flat <- unlist(cds_mapped)
  cds_flat$transcript_id <- rep(names(cds_mapped), lengths(cds_mapped))
  
  # change Feature type
  cds_flat$type <- "CDS"
  
  # Validation of mapping
  valid_mappings <- !is.na(cds_mapped@partitioning@end)
  cat("Valid mappings:", sum(valid_mappings), "of", length(valid_mappings), "\n")

  # add phase information
  cds_phased <- calculate_cds_phase(cds_gr = cds_flat, exons_by_tx = novel_exons_by_tx_coding)
  
  return(cds_phased)

}

## generate UTRs from cds granges
get_utrs <- function(gtf, cds_phased) {
  # take only exons of gtf
  exons <-gtf[gtf$type == "exon"]
  
  # split exons by transcript_id
  exons_by_tx <- split(exons, exons$transcript_id)
  
  utr_granges <- GRanges()
  tx_ids <- unique(names(cds_phased))
  
  for (tx in tx_ids) {
    cds_tx <- cds_phased[cds_phased$transcript_id == tx] # get the right cds
    tx_model <- unlist(exons_by_tx[tx]) # get the right tx model
    
    if (as.character(strand(tx_model[1])) == "+") { # for + strand only
      
      cds_start <- min(start(cds_tx)) # get the start of the CDS
      
      utrs <- unlist(subtract(tx_model, cds_tx, ignore.strand = FALSE)) # subtract CDS from exons
      
      mcols(utrs) <- mcols(cds_tx[1]) # readd metadata
      mcols(utrs)$type <- "five_prime_UTR" # change type
      mcols(utrs)$type[start(utrs) > cds_start] <- "three_prime_UTR"
    } else if (as.character(strand(tx_model[1])) == "-") {
      cds_stop <- max(end(cds_tx))
      utrs <- unlist(subtract(tx_model, cds_tx, ignore.strand = FALSE)) # subtract CDS from exons
      
      mcols(utrs) <- mcols(cds_tx[1]) # readd metadata
      mcols(utrs)$type <- "three_prime_UTR" # change type
      mcols(utrs)$type[end(utrs) < cds_stop] <- "five_prime_UTR"
    }
    
    utr_granges <- c(utr_granges, utrs)
  }
  
  utr_granges$phase<-NA
  
  
  return(utr_granges)
}


