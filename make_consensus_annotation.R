library(GenomicRanges)
library(rtracklayer)

setwd("D:/Dropbox/SauvageauLab Dropbox/SauvageauLab/Members/EHRESMANN_Sophie/Projects/EMT_Timecourse/sequencing/nanopore/dorado_minimap2_stringtie_sqanti_merge/")

setwd("/Users/ehresms/SauvageauLab Dropbox/SauvageauLab/Members/EHRESMANN_Sophie/Projects/EMT_Timecourse/sequencing/nanopore/dorado_minimap2_stringtie_sqanti_merge/")

source("make_consensus_annotations_functions.R")
# import stringtie gtf as granges
stringtie_gtf <- import("SE_MCF10a_Nanopore_Dorado_Restrander_Minimap2_Stringtie3_Sqanti_merge.gtf")

# import gencode gff as granges
gencode_gff <- import("Homo_sapiens-GCA_009914755.4-2022_07-genes_chr.gff3")
# get only genes
gencode_genes <- gencode_gff[grep(mcols(gencode_gff)$type, pattern = "gene")]

# annotate stringtie transcripts with either gene_id/ref_gene_id if it is ENSG
stringtie_gtf_unique <- concat_gene_and_ref_gene(stringtie_gtf)
stringtie_gtf <- stringtie_gtf_unique

#### flatten exons of stringtie transcripts to create a gene feature ####

# get exons only
stringtie_exons <- stringtie_gtf[stringtie_gtf$type == "exon"]

stringtie_gene_flat_combined <- make_gene_from_exons(exon_gr = stringtie_exons)

### add metadata from gff ####
# only add columns with data - can adjust
columns_to_add <- c("biotype", "description", "Name", "version")

stringtie_gene_flat_combined_gff <- add_metadata_from_gff_genes(
  gene_gr = stringtie_gene_flat_combined,
  gff_gene_gr = gencode_genes,
  columns_to_add = columns_to_add
)
## make tx2gene with unique_gene_id
library(dplyr)

stringtie_transcripts <- stringtie_gtf[stringtie_gtf$type == "transcript"]
tx2gene <- as.data.frame(unique(mcols(stringtie_transcripts)[, c("transcript_id", "unique_gene_id")]))

write.csv(file = "SE_MCF10a_tx2gene.csv", tx2gene)

### add metadata to transcripts ####

# filter gencode_gff for transcripts to add data per transcript
gencode_transcripts <- gencode_gff[gencode_gff$type %in% c(
  "lnc_RNA", "pseudogenic_transcript", "mRNA", "snRNA", "miRNA",
  "ncRNA", "unconfirmed_transcript", "snoRNA", "scRNA", "rRNA"
)]

## annotate using same columns to add as for genes
stringtie_gtf_annotated <- add_metadata_from_gff_transcripts(
  gtf = stringtie_gtf,
  gencode_transcript_gff = gencode_transcripts,
  columns_to_add = columns_to_add
)

#### make CDS and UTR features ####
library(stringr)

# take the features already made
gencode_UTR_CDS <- gencode_gff[gencode_gff$type %in% c("CDS", "five_prime_UTR", "three_prime_UTR")]
gencode_UTR_CDS_transcript_id <- sapply(
  str_split(as.character(
    gencode_UTR_CDS$Parent
  ), pattern = ":"),
  "[", 2
)
# add transcript_id

gencode_UTR_CDS$transcript_id <- gencode_UTR_CDS_transcript_id

# get UTR/CDS features only for transcript_ids in stringtie gtf
gencode_UTR_CDS_select <- gencode_UTR_CDS[gencode_UTR_CDS_transcript_id %in% stringtie_gtf_known$transcript_id] ### use these annotations for final file

# get start/stop features using gencode CDS

gencode_CDS_select <- gencode_UTR_CDS_select[gencode_UTR_CDS_select$type == "CDS"]

gencode_start_stop <- get_start_and_stop(gencode_CDS_select)

gencode_UTR_CDS_start_stop <- c(gencode_UTR_CDS_select, gencode_start_stop)

#### generate features for novel transcripts ####
library(BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0)
library(Biostrings)

# get fasta sequences of novel transcripts
genome <- BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0

stringtie_gtf_novel <- stringtie_gtf_annotated[stringtie_gtf$transcript_id %notin% gencode_gff$transcript_id]

novel_transcripts_fasta <- get_fasta_of_transcripts(transcript_gr = stringtie_gtf_novel, genome = genome)

writeXStringSet(novel_transcript_fasta, filepath = "SE_MCF10a_noveltranscripts.fa")

# run CPC2 locally
system2("python", args = c(
  "/Users/ehresms/computational/tools/CPC2_standalone-1.0.1/bin/CPC2.py",
  "-i", "SE_MCF10a_noveltranscripts.fa",
  "-o", "SE_MCF10a_noveltranscripts_cpc2_results",
  "--ORF"
))
# read CPC2 results
cpc2_result <- read.table("SE_MCF10a_noveltranscripts_cpc2_results.txt", sep = "\t", header = FALSE)

# get cds ranges for novel transcripts (transcript-relative mapping)
cds_novel <- get_cds_from_cpc2(cpc2_result = cpc2_result)

# map back to transcript structure
library(GenomicFeatures)

# get granges of cds per transcript with phase information
cds_phased<-map_cds_to_transcript(cds_gr = cds_novel,
                                  novel_gtf = stringtie_gtf_novel)
  

#### generate start and stop gr for novel ####

novel_start_stop <- get_start_and_stop(cds_phased)


## generate utrs form CDS features and transcripts
utr_granges <- get_utrs(gtf = stringtie_gtf_novel,
                        cds_phased = cds_phased)


novel_cds_utr_start_stop<-c(cds_phased,
                            utr_granges,
                            novel_start_stop)

novel_cds_utr_start_stop<-sort(novel_cds_utr_start_stop)
novel_cds_utr_start_stop<-unname(novel_cds_utr_start_stop)
mcols(novel_cds_utr_start_stop)<-mcols(novel_cds_utr_start_stop)[,!grepl(colnames(mcols(novel_cds_utr_start_stop)),
                                                                        pattern = "mcols")]


novel_test<-novel_cds_utr_start_stop[mcols(novel_cds_utr_start_stop)$unique_gene_id == "MCF10a.21"]

########


rtracklayer::export(novel_cds_utr_start_stop, format = "gff3", con = "SE_MCF10a_novel_cds_utr_start_stop")

#### merge everything together ####

# gene_gr is stringtie_gene_flat_combined
# tx_gr is stringtie_gtf

## merge gencode other features
gencode_cds_utr_start_stop <- c(gencode_UTR_CDS_select, gencode_starts, gencode_stops)
rtracklayer::export(gencode_cds_utr_start_stop, format = "gff3", con = "SE_MCF10a_gencode_cds_utr_start_stop")

## harmonize mcols
columns_needed <- c(
  "source", "type", "score", "phase", "ID", "biotype", "description", "unique_gene_id",
  "gene_id", "version", "Parent", "transcript_id", "exon_id", "exon_number", "Name", "protein_id",
  "gene_name", "ref_gene_id"
)

gr_list <- list(
  stringtie_gene_flat_combined,
  stringtie_gtf,
  novel_cds_utr_start_stop,
  gencode_cds_utr_start_stop
)


gr_list_filled <- lapply(gr_list, function(gr) {
  missing <- setdiff(columns_needed, colnames(mcols(gr)))
  for (col in missing) {
    mcols(gr)[[col]] <- NA
  }
  mcols(gr) <- mcols(gr)[, columns_needed]
  gr
})

## make final gff

final_gff <- do.call(c, gr_list_filled)

final_gff <- sort(final_gff)

final_gff$original_gene_id <- final_gff$gene_id
final_gff$gene_id <- final_gff$unique_gene_id

final_gff$source <- "StringTie"

rtracklayer::export(final_gff, format = "gff3", con = "SE_MCF10a_stringtie_annotated_fromGencode.gff3")
