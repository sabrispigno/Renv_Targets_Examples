# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline

# Load packages required to define the pipeline:
library(targets)
# library(tarchetypes) # Load other packages as needed.

# Set target options:
tar_option_set(
  packages = c("dada2", "phyloseq", "Biostrings", "ggplot2"), # Packages that your targets need for their tasks.
  format = "rds", # Optionally set the default storage format. qs is fast.
  #
  # Pipelines that take a long time to run may benefit from
  # optional distributed computing. To use this capability
  # in tar_make(), supply a {crew} controller
  # as discussed at https://books.ropensci.org/targets/crew.html.
  # Choose a controller that suits your needs. For example, the following
  # sets a controller that scales up to a maximum of two workers
  # which run as local R processes. Each worker launches when there is work
  # to do and exits if 60 seconds pass with no tasks to run.
  #
  #   controller = crew::crew_controller_local(workers = 2, seconds_idle = 60)
  #
  # Alternatively, if you want workers to run on a high-performance computing
  # cluster, select a controller from the {crew.cluster} package.
  # For the cloud, see plugin packages like {crew.aws.batch}.
  # The following example is a controller for Sun Grid Engine (SGE).
  #
  #   controller = crew.cluster::crew_controller_sge(
  #     # Number of workers that the pipeline can scale up to:
  #     workers = 10,
  #     # It is recommended to set an idle time so workers can shut themselves
  #     # down if they are not running tasks.
  #     seconds_idle = 120,
  #     # Many clusters install R as an environment module, and you can load it
  #     # with the script_lines argument. To select a specific verison of R,
  #     # you may need to include a version string, e.g. "module load R/4.3.2".
  #     # Check with your system administrator if you are unsure.
  #     script_lines = "module load R"
  #   )
  #
  # Set other options as needed.
)

# Run the R scripts in the R/ folder with your custom functions:
tar_source()
# tar_source("other_functions.R") # Source other scripts as needed.
# Define pipeline
list(
  
  # ---- 1. Input path ----
  tar_target(
    dada_path,
    "MiSeq_SOP",
    format = "file"
  ),
  
  # ---- 2. Load FASTQ files ----
  tar_target(
    dada_fastq_files,
    load_fastq_files(dada_path)
  ),
  
  # ---- 3. Filter and trim ----
  tar_target(
    dada_filt,
    filter_and_trim(
      dada_fastq_files$fnFs,
      dada_fastq_files$fnRs,
      dada_fastq_files$sample.names,
      dada_path
    )
  ),
  
  # ---- 4. Learn error rates ----
  tar_target(
    dada_errs,
    learn_error_rates(dada_filt$filtFs, dada_filt$filtRs)
  ),
  
  # ---- 5. Denoise ----
  tar_target(
    dada_inf,
    infer_samples(
      dada_filt$filtFs,
      dada_filt$filtRs,
      dada_errs$errF,
      dada_errs$errR
    )
  ),
  
  # ---- 6. Merge paired reads ----
  tar_target(
    dada_mergers,
    merge_pairs(
      dada_inf$dadaFs,
      dada_filt$filtFs,
      dada_inf$dadaRs,
      dada_filt$filtRs
    )
  ),
  
  # ---- 7. Sequence table & chimera removal ----
  tar_target(
    dada_seqtab,
    make_sequence_table(dada_mergers)$seqtab
  ),
  
  tar_target(
    dada_seqtab_nochim,
    make_sequence_table(dada_mergers)$seqtab.nochim
  ),
  
  # ---- 8. Track reads ----
  tar_target(
    dada_track,
    track_reads(
      dada_filt$out,
      dada_inf$dadaFs,
      dada_inf$dadaRs,
      dada_mergers,
      dada_seqtab_nochim,
      dada_fastq_files$sample.names
    )
  ),
  
  # ---- 9. Assign taxonomy ----
  tar_target(
    dada_taxa,
    assign_taxonomy(
      dada_seqtab_nochim,
      file.path(dada_path, "silva_nr_v132_train_set.fa.gz")
    )
  ),
  
  # ---- 10. Evaluate accuracy ----
  tar_target(
    dada_eval,
    evaluate_accuracy(dada_seqtab_nochim, dada_path)
  ),
  
  # ---- 11. Sample metadata ----
  tar_target(
    dada_samdf,
    make_sample_metadata(dada_seqtab_nochim)
  ),
  
  # ---- 12. Build & prepare phyloseq object ----
  tar_target(
    dada_ps,
    {
      ps <- make_phyloseq_object(dada_seqtab_nochim, dada_taxa, dada_samdf)
      prep_phyloseq(ps)
    }
  ),
  
  # ---- 13. Alpha diversity plot ----
  tar_target(
    dada_plot_alpha,
    {
      p <- plot_alpha_diversity(dada_ps)
      dir.create("results", showWarnings = FALSE)
      ggsave("results/alpha_diversity.png", p, width = 6, height = 4)
      "results/alpha_diversity.png"
    },
    format = "file"
  ),
  
  # ---- 14. NMDS ordination plot ----
  tar_target(
    dada_plot_ord,
    {
      p <- plot_ordination_nmds(dada_ps)
      dir.create("results", showWarnings = FALSE)
      ggsave("results/ordination_nmds.png", p, width = 6, height = 4)
      "results/ordination_nmds.png"
    },
    format = "file"
  ),
  
  # ---- 15. Top 20 taxa barplot ----
  tar_target(
    dada_plot_bar,
    {
      p <- plot_top20_taxa(dada_ps)
      dir.create("results", showWarnings = FALSE)
      ggsave("results/bar_top20_taxa.png", p, width = 6, height = 4)
      "results/bar_top20_taxa.png"
    },
    format = "file"
  )
)