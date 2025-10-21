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

# Define the pipeline
list(
  # Input path
  tar_target(
    path, 
    "MiSeq_SOP", 
    format = "file"
  ),
  
  # Load FASTQ files
  tar_target(
    fastq_files,
    load_fastq_files(path)
  ),
  
  # Filter and trim
  tar_target(
    filt,
    filter_and_trim(fastq_files$fnFs, fastq_files$fnRs, fastq_files$sample.names, path)
  ),
  
  # Learn errors
  tar_target(
    errs,
    learn_error_rates(filt$filtFs, filt$filtRs)
  ),
  
  # Denoise
  tar_target(
    inf,
    infer_samples(filt$filtFs, filt$filtRs, errs$errF, errs$errR)
  ),
  
  # Merge
  tar_target(
    mergers,
    merge_pairs(inf$dadaFs, filt$filtFs, inf$dadaRs, filt$filtRs)
  ),
  
  # Sequence table and chimera removal (split!)
  tar_target(
    seqtab,
    make_sequence_table(mergers)$seqtab
  ),
  
  tar_target(
    seqtab_nochim,
    make_sequence_table(mergers)$seqtab.nochim
  ),
  
  # Tracking reads
  tar_target(
    track,
    track_reads(filt$out, inf$dadaFs, inf$dadaRs, mergers, seqtab_nochim, fastq_files$sample.names)
  ),
  
  # Assign taxonomy
  tar_target(
    taxa,
    assign_taxonomy(seqtab_nochim, file.path(path, "silva_nr_v132_train_set.fa.gz"))
  ),
  
  # Evaluate accuracy (optional)
  tar_target(
    evaluation,
    evaluate_accuracy(seqtab_nochim, path)
  ),
  
  # Sample metadata
  tar_target(
    dada_samdf,
    make_sample_metadata(seqtab_nochim)
  ),
  
  # Build and prepare phyloseq object
  tar_target(
    dada_ps,
    {
      ps <- make_phyloseq_object(seqtab_nochim, taxa, dada_samdf)
      prep_phyloseq(ps)
    }
  ),
  
  # Alpha diversity plot
  tar_target(
    plot_alpha_file,
    {
      p <- plot_alpha_diversity(dada_ps)
      ggsave("results/alpha_diversity.png", p, width = 6, height = 4)
      "results/alpha_diversity.png"
    },
    format = "file"
  ),
  
  # NMDS ordination plot
  tar_target(
    plot_ord,
    {
      p <- plot_ordination_nmds(dada_ps)
      ggsave("results/ordination_nmds.png", p, width = 6, height = 4)
      "results/ordination_nmds.png"
    },
    format = "file"
  ),
  
  # Top 20 taxa barplot
  tar_target(
    plot_bar,
    {
      p <- plot_top20_taxa(dada_ps)
      ggsave("results/Bar_top20_taxa.png", p, width = 6, height = 4)
      "results/Bar_top20_taxa.png"
    },
    format = "file"
  )
)
