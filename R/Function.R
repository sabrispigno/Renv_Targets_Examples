#-------------------------------
# 1. Load FASTQ files
#-------------------------------
load_fastq_files <- function(path) {
  fnFs <- sort(list.files(path, pattern = "_R1_001.fastq", full.names = TRUE))
  fnRs <- sort(list.files(path, pattern = "_R2_001.fastq", full.names = TRUE))
  sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
  list(fnFs = fnFs, fnRs = fnRs, sample.names = sample.names)
}

#-------------------------------
# 2. Filter and Trim
#-------------------------------
filter_and_trim <- function(fnFs, fnRs, sample.names, path, truncLen = c(240, 160)) {
  filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
  filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
  names(filtFs) <- sample.names
  names(filtRs) <- sample.names
  
  out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                       truncLen = truncLen,
                       maxN = 0, maxEE = c(2, 2),
                       truncQ = 2, rm.phix = TRUE,
                       compress = TRUE, multithread = TRUE)
  list(filtFs = filtFs, filtRs = filtRs, out = out)
}

#-------------------------------
# 3. Learn Error Rates
#-------------------------------
learn_error_rates <- function(filtFs, filtRs) {
  errF <- learnErrors(filtFs, multithread = TRUE)
  errR <- learnErrors(filtRs, multithread = TRUE)
  list(errF = errF, errR = errR)
}

#-------------------------------
# 4. Sample Inference
#-------------------------------
infer_samples <- function(filtFs, filtRs, errF, errR) {
  dadaFs <- dada(filtFs, err = errF, multithread = TRUE)
  dadaRs <- dada(filtRs, err = errR, multithread = TRUE)
  list(dadaFs = dadaFs, dadaRs = dadaRs)
}

#-------------------------------
# 5. Merge Pairs
#-------------------------------
merge_pairs <- function(dadaFs, filtFs, dadaRs, filtRs) {
  mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)
}

#-------------------------------
# 6. Construct Sequence Table & Remove Chimeras
#-------------------------------
make_sequence_table <- function(mergers) {
  seqtab <- makeSequenceTable(mergers)
  seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus",
                                      multithread = TRUE, verbose = TRUE)
  list(seqtab = seqtab, seqtab.nochim = seqtab.nochim)
}

#-------------------------------
# 7. Track Reads Through Pipeline
#-------------------------------
track_reads <- function(out, dadaFs, dadaRs, mergers, seqtab.nochim, sample.names) {
  getN <- function(x) sum(getUniques(x))
  track <- cbind(out,
                 sapply(dadaFs, getN),
                 sapply(dadaRs, getN),
                 sapply(mergers, getN),
                 rowSums(seqtab.nochim))
  colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
  rownames(track) <- sample.names
  track
}

#-------------------------------
# 8. Assign Taxonomy
#-------------------------------
assign_taxonomy <- function(seqtab.nochim, ref_db) {
  taxa <- assignTaxonomy(seqtab.nochim, ref_db, multithread = TRUE)
  taxa
}

#-------------------------------
# 9. Evaluate Accuracy (for mock data)
#-------------------------------
evaluate_accuracy <- function(seqtab.nochim, path) {
  unqs.mock <- seqtab.nochim["Mock", ]
  unqs.mock <- sort(unqs.mock[unqs.mock > 0], decreasing = TRUE)
  cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n")
  
  mock.ref <- getSequences(file.path(path, "HMP_MOCK.v35.fasta"))
  match.ref <- sum(sapply(names(unqs.mock), function(x) any(grepl(x, mock.ref))))
  cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n")
}

#-----------------------------------
# 10. Create Sample Metadata
#-----------------------------------
make_sample_metadata <- function(seqtab.nochim) {
  samples.out <- rownames(seqtab.nochim)
  subject <- sapply(strsplit(samples.out, "D"), `[`, 1)
  gender <- substr(subject, 1, 1)
  subject <- substr(subject, 2, 999)
  day <- as.integer(sapply(strsplit(samples.out, "D"), `[`, 2))
  
  samdf <- data.frame(Subject = subject, Gender = gender, Day = day)
  samdf$When <- ifelse(samdf$Day > 100, "Late", "Early")
  rownames(samdf) <- samples.out
  samdf
}

#-----------------------------------
# 11. Make Phyloseq Object
#-----------------------------------
make_phyloseq_object <- function(seqtab.nochim, taxa, samdf) {
  phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE),
           sample_data(samdf),
           tax_table(taxa))
}

#-----------------------------------
# 12. Prepare Phyloseq (cleanup, rename, add DNA)
#-----------------------------------
prep_phyloseq <- function(ps) {
  ps <- prune_samples(sample_names(ps) != "Mock", ps)  # remove mock
  dna <- Biostrings::DNAStringSet(taxa_names(ps))
  names(dna) <- taxa_names(ps)
  ps <- merge_phyloseq(ps, dna)
  taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
  ps
}

#-----------------------------------
# 13. Plot Alpha Diversity
#-----------------------------------
plot_alpha_diversity <- function(ps) {
  p <- plot_richness(ps, x = "Day", measures = c("Shannon", "Simpson"), color = "When")
  print(p)
  p
}

#-----------------------------------
# 14. Ordination (NMDS Bray-Curtis)
#-----------------------------------
plot_ordination_nmds <- function(ps) {
  ps.prop <- transform_sample_counts(ps, function(otu) otu / sum(otu))
  ord.nmds.bray <- ordinate(ps.prop, method = "NMDS", distance = "bray")
  p <- plot_ordination(ps.prop, ord.nmds.bray, color = "When", title = "Bray NMDS")
  print(p)
  p
}

#-----------------------------------
# 15. Plot Top 20 Taxa Barplot
#-----------------------------------
plot_top20_taxa <- function(ps) {
  top20 <- names(sort(taxa_sums(ps), decreasing = TRUE))[1:20]
  ps.top20 <- transform_sample_counts(ps, function(OTU) OTU / sum(OTU))
  ps.top20 <- prune_taxa(top20, ps.top20)
  p <- plot_bar(ps.top20, x = "Day", fill = "Family") + facet_wrap(~When, scales = "free_x")
  print(p)
  p
}



