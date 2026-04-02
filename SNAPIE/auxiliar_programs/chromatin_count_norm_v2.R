#!/usr/bin/env Rscript

#' =============================================================================
#' Chromatin Fragment Counter - Quantification Pipeline
#' =============================================================================
#' 
#' Date: 2025-08-31
#' Description: Pipeline for counting chromatin fragments at user-defined sites
#'              Generates raw counts, CPM, and optionally normalized matrices
#' 
#' Usage: Rscript chromatin_counter.R [options]
#'        ./chromatin_counter.R [options]
#' 
#' Dependencies: See install_required_packages() function below
#' =============================================================================

#' Display comprehensive help information
show_help <- function() {
    cat("=============================================================================\n")
    cat("Chromatin Fragment Counter - Quantification Pipeline\n")
    cat("=============================================================================\n\n")
    
    cat("DESCRIPTION:\n")
    cat("  Pipeline for counting chromatin fragments at user-defined genomic sites.\n")
    cat("  Supports both batch mode (multiple samples via samplesheet) and single sample mode.\n")
    cat("  Generates three types of output matrices:\n")
    cat("  1. Raw fragment counts at target sites\n")
    cat("  2. CPM (Counts Per Million mapped fragments) at target sites\n")
    cat("  3. Reference-normalized counts (optional)\n\n")
    
    cat("USAGE:\n")
    cat("  # Batch mode (multiple samples)\n")
    cat("  Rscript chromatin_counter.R --samplesheet PATH --frags-dir PATH --target-sites PATH [OPTIONS]\n")
    cat("  \n")
    cat("  # Single sample mode\n")
    cat("  Rscript chromatin_counter.R --sample-name NAME --fragment-file PATH --target-sites PATH [OPTIONS]\n\n")
    
    cat("BATCH MODE (REQUIRED):\n")
    cat("  --samplesheet PATH      Samplesheet file (TSV format)\n")
    cat("  --frags-dir PATH        Directory containing fragment files (batch mode, default: ./frags)\n")
    cat("  --target-sites PATH     Target sites BED file\n\n")
    
    cat("SINGLE SAMPLE MODE (REQUIRED):\n")
    cat("  --sample-name NAME      Name for the single sample\n")
    cat("  --fragment-file PATH    Fragment BED file for the sample\n")
    cat("  --target-sites PATH     Target sites BED file\n\n")
    
    cat("OPTIONAL PARAMETERS:\n")
    cat("  -h, --help              Show this help message and exit\n")
    cat("  --output-dir PATH       Output directory (default: ./output)\n")
    cat("  --reference-sites PATH  Reference sites for normalization (optional)\n")
    cat("  --include-all-chr       Include all chromosomes (default: autosomes only)\n")
    cat("  --regenerate-counts     Force regeneration of count matrices\n")
    cat("  --save-intermediate     Save intermediate matrices as text files\n")
    cat("  --verbose              Enable verbose output\n\n")
    
    cat("EXAMPLES:\n")
    cat("  # Batch mode - basic usage (raw counts and CPM only)\n")
    cat("  ./chromatin_counter.R \\\n")
    cat("    --samplesheet samples.tsv \\\n")
    cat("    --frags-dir ./frags \\\n")
    cat("    --target-sites peaks_of_interest.bed\n\n")
    
    cat("  # Batch mode - with reference normalization\n")
    cat("  ./chromatin_counter.R \\\n")
    cat("    --samplesheet samples.tsv \\\n")
    cat("    --frags-dir ./frags \\\n")
    cat("    --target-sites peaks_of_interest.bed \\\n")
    cat("    --reference-sites housekeeping_regions.bed\n\n")
    
    cat("  # Single sample mode - basic usage\n")
    cat("  ./chromatin_counter.R \\\n")
    cat("    --sample-name sample_001 \\\n")
    cat("    --fragment-file ./frags/sample_001.bed \\\n")
    cat("    --target-sites peaks_of_interest.bed\n\n")
    
    cat("  # Single sample mode - with reference normalization\n")
    cat("  ./chromatin_counter.R \\\n")
    cat("    --sample-name sample_001 \\\n")
    cat("    --fragment-file ./frags/sample_001.bed \\\n")
    cat("    --target-sites peaks_of_interest.bed \\\n")
    cat("    --reference-sites dhs_sites.bed\n\n")
    
    cat("  # Include all chromosomes (including sex chromosomes and contigs)\n")
    cat("  ./chromatin_counter.R \\\n")
    cat("    --samplesheet samples.tsv \\\n")
    cat("    --target-sites peaks_of_interest.bed \\\n")
    cat("    --include-all-chr \\\n")
    cat("    --verbose\n\n")
    
    cat("=============================================================================\n")
    cat("INPUT FILE FORMATS\n")
    cat("=============================================================================\n\n")
    
    cat("1. SAMPLESHEET FILE (TSV format, batch mode only):\n")
    cat("   sample_name\n")
    cat("   sample_001\n")
    cat("   sample_002\n")
    cat("   sample_003\n")
    cat("   sample_004\n")
    cat("   \n")
    cat("   Requirements:\n")
    cat("   - Tab-separated format (or single column for sample names only)\n")
    cat("   - Header line required\n")
    cat("   - sample_name column: unique identifier for each sample\n")
    cat("   - Additional columns ignored (conditions not needed)\n\n")
    
    cat("2. TARGET SITES BED FILE (required):\n")
    cat("   chr1	1000	2000	peak_1\n")
    cat("   chr1	5000	6000	peak_2\n")
    cat("   chr2	10000	11000	peak_3\n")
    cat("   chr3	15000	16000	peak_4\n")
    cat("   \n")
    cat("   Requirements:\n")
    cat("   - Standard BED format (minimum 3 columns)\n")
    cat("   - Column 1: Chromosome (chr1, chr2, etc.)\n")
    cat("   - Column 2: Start position (0-based)\n")
    cat("   - Column 3: End position (1-based)\n")
    cat("   - Optional Column 4: Peak name/ID (auto-generated if missing)\n")
    cat("   - Tab-separated format\n")
    cat("   - No header line\n\n")
    
    cat("3. REFERENCE SITES BED FILE (optional, for normalization):\n")
    cat("   chr1	50000	51000	ref_1\n")
    cat("   chr2	100000	101000	ref_2\n")
    cat("   chr3	150000	151000	ref_3\n")
    cat("   \n")
    cat("   Same format as target sites BED file.\n")
    cat("   Used as denominator for normalization calculations.\n")
    cat("   Examples: housekeeping genes, DHS sites, or other reference regions.\n\n")
    
    cat("4. FRAGMENT FILES (required):\n")
    cat("   Location: --frags-dir (default: ./frags)\n")
    cat("   Pattern: {sample_name}*.bed\n")
    cat("   Format: Standard BED format containing fragment coordinates\n")
    cat("   chr1	1500	1600\n")
    cat("   chr1	2000	2100\n")
    cat("   chr2	5500	5650\n\n")
    
    cat("=============================================================================\n")
}

# Check for help request before loading any packages
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0 && (args[1] %in% c("-h", "--help", "help"))) {
    show_help()
    quit(save = "no", status = 0)
}

#' Parse command line arguments
parse_arguments <- function() {
    # Default parameters
    params <- list(
        samplesheet = NULL,
        sample_name = NULL,
        fragment_file = NULL,
        target_sites = NULL,
        reference_sites = NULL,
        frags_dir = "./frags", 
        output_dir = "./output",
        include_all_chr = FALSE,
        regenerate_counts = FALSE,
        save_intermediate = FALSE,
        verbose = FALSE
    )
    
    args <- commandArgs(trailingOnly = TRUE)
    
    i <- 1
    while (i <= length(args)) {
        arg <- args[i]
        
        if (arg == "--samplesheet" && i < length(args)) {
            params$samplesheet <- args[i + 1]
            i <- i + 2
        } else if (arg == "--sample-name" && i < length(args)) {
            params$sample_name <- args[i + 1]
            i <- i + 2
        } else if (arg == "--fragment-file" && i < length(args)) {
            params$fragment_file <- args[i + 1]
            i <- i + 2
        } else if (arg == "--target-sites" && i < length(args)) {
            params$target_sites <- args[i + 1]
            i <- i + 2
        } else if (arg == "--reference-sites" && i < length(args)) {
            params$reference_sites <- args[i + 1]
            i <- i + 2
        } else if (arg == "--frags-dir" && i < length(args)) {
            params$frags_dir <- args[i + 1]
            i <- i + 2
        } else if (arg == "--output-dir" && i < length(args)) {
            params$output_dir <- args[i + 1]
            i <- i + 2
        } else if (arg == "--include-all-chr") {
            params$include_all_chr <- TRUE
            i <- i + 1
        } else if (arg == "--regenerate-counts") {
            params$regenerate_counts <- TRUE
            i <- i + 1
        } else if (arg == "--save-intermediate") {
            params$save_intermediate <- TRUE
            i <- i + 1
        } else if (arg == "--verbose") {
            params$verbose <- TRUE
            i <- i + 1
        } else {
            cat("Warning: Unknown argument:", arg, "\n")
            i <- i + 1
        }
    }
    
    # Validate required parameters based on mode
    if (is.null(params$target_sites)) {
        cat("Error: --target-sites is required\n")
        cat("Use --help for usage information\n")
        quit(save = "no", status = 1)
    }
    
    # Check mode validity
    batch_mode <- !is.null(params$samplesheet)
    single_mode <- !is.null(params$sample_name) && !is.null(params$fragment_file)
    
    if (!batch_mode && !single_mode) {
        cat("Error: Must specify either:\n")
        cat("  Batch mode: --samplesheet\n")
        cat("  Single mode: --sample-name AND --fragment-file\n")
        cat("Use --help for usage information\n")
        quit(save = "no", status = 1)
    }
    
    if (batch_mode && single_mode) {
        cat("Error: Cannot use both batch mode and single sample mode simultaneously\n")
        cat("Choose either --samplesheet OR (--sample-name + --fragment-file)\n")
        quit(save = "no", status = 1)
    }
    
    # Validate single sample mode requirements
    if (single_mode) {
        if (is.null(params$sample_name)) {
            cat("Error: --sample-name is required for single sample mode\n")
            quit(save = "no", status = 1)
        }
        if (is.null(params$fragment_file)) {
            cat("Error: --fragment-file is required for single sample mode\n")
            quit(save = "no", status = 1)
        }
        if (!file.exists(params$fragment_file)) {
            cat("Error: Fragment file does not exist:", params$fragment_file, "\n")
            quit(save = "no", status = 1)
        }
    }
    
    # Set analysis mode
    params$analysis_mode <- ifelse(batch_mode, "batch", "single")
    
    return(params)
}

#' Install required packages if not already installed
install_required_packages <- function() {
    cat("Checking and installing required packages...\n")
    
    # CRAN packages
    cran_packages <- c("dplyr", "readr", "stringr")
    
    # Bioconductor packages
    bioc_packages <- c("rtracklayer", "GenomicRanges")
    
    # Install CRAN packages
    for (pkg in cran_packages) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
            cat("Installing", pkg, "from CRAN...\n")
            install.packages(pkg, dependencies = TRUE, quiet = TRUE)
        }
    }
    
    # Install Bioconductor packages
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager", quiet = TRUE)
    }
    
    for (pkg in bioc_packages) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
            cat("Installing", pkg, "from Bioconductor...\n")
            BiocManager::install(pkg, dependencies = TRUE, quiet = TRUE, ask = FALSE)
        }
    }
    
    cat("Package installation completed.\n")
}

#' Load required libraries with error handling
load_required_libraries <- function() {
    cat("Loading required libraries...\n")
    
    required_libraries <- c("rtracklayer", "GenomicRanges", "dplyr", "readr", "stringr")
    
    for (lib in required_libraries) {
        tryCatch({
            #library(lib, character.only = TRUE, quietly = TRUE)
            suppressPackageStartupMessages(library(lib, character.only = TRUE, quietly = TRUE))
        }, error = function(e) {
            cat("Error loading library:", lib, "\n")
            cat("Error message:", e$message, "\n")
            quit(save = "no", status = 1)
        })
    }
    
    cat("All libraries loaded successfully.\n")
}

#' Log message with timestamp
log_message <- function(message, verbose = TRUE) {
    if (verbose) {
        cat("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", message, "\n", sep = "")
    }
}

#' Create organized output directory structure
create_output_structure <- function(output_dir, verbose = FALSE) {
    log_message("=== Creating Output Directory Structure ===", verbose)
    
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
    }
    
    subdirs <- c("matrices", "logs")
    
    for (subdir in subdirs) {
        subdir_path <- file.path(output_dir, subdir)
        if (!dir.exists(subdir_path)) {
            dir.create(subdir_path, recursive = TRUE)
            log_message(paste("Created directory:", subdir_path), verbose)
        }
    }
}

#' Read and validate samplesheet
read_and_validate_samplesheet <- function(samplesheet_path, verbose = FALSE) {
    log_message("=== Reading and Validating Samplesheet ===", verbose)
    
    if (!file.exists(samplesheet_path)) {
        stop("Samplesheet file does not exist: ", samplesheet_path)
    }
    
    tryCatch({
        samplesheet <- read.table(samplesheet_path, 
                                 header = TRUE, 
                                 sep = "\t", 
                                 stringsAsFactors = FALSE)
    }, error = function(e) {
        stop("Error reading samplesheet: ", e$message, 
             "\nEnsure file is tab-separated with header containing sample_name column")
    })
    
    log_message(paste("Read samplesheet with", nrow(samplesheet), "samples"), verbose)
    
    # Validate required sample_name column
    if (!"sample_name" %in% colnames(samplesheet)) {
        stop("Samplesheet missing required column: sample_name")
    }
    
    # Check for duplicate samples
    if (any(duplicated(samplesheet$sample_name))) {
        duplicates <- samplesheet$sample_name[duplicated(samplesheet$sample_name)]
        stop("Duplicate sample names found: ", paste(duplicates, collapse = ", "))
    }
    
    return(samplesheet)
}

#' Find fragment files for each sample
find_fragment_files <- function(sample_names, frags_dir, verbose = FALSE) {
    log_message("=== Finding Fragment Files ===", verbose)
    log_message(paste("Search directory:", frags_dir), verbose)
    
    if (!dir.exists(frags_dir)) {
        stop("Fragment directory does not exist: ", frags_dir)
    }
    
    matched_files <- list()
    missing_samples <- character()
    
    for (sample_name in sample_names) {
        # Find files matching the pattern
        found_files <- list.files(frags_dir, 
                                 pattern = paste0(sample_name, ".*\\.bed"), 
                                 full.names = TRUE, 
                                 recursive = TRUE)
        
        if (length(found_files) == 0) {
            missing_samples <- c(missing_samples, sample_name)
            log_message(paste("WARNING: No fragment files found for sample:", sample_name), verbose)
        } else if (length(found_files) == 1) {
            matched_files[[sample_name]] <- found_files[1]
            log_message(paste("Found fragment file for", sample_name, ":", basename(found_files[1])), verbose)
        } else {
            matched_files[[sample_name]] <- found_files[1]
            log_message(paste("WARNING: Multiple fragment files found for", sample_name, ". Using:", basename(found_files[1])), verbose)
        }
    }
    
    log_message(paste("Successfully matched", length(matched_files), "out of", length(sample_names), "fragment files"), verbose)
    
    if (length(missing_samples) > 0) {
        stop("Missing fragment files for samples: ", paste(missing_samples, collapse = ", "))
    }
    
    return(matched_files)
}

#' Read BED file and create GRanges object
read_bed_file <- function(bed_file, file_type = "sites", verbose = FALSE) {
    log_message(paste("=== Reading", file_type, "BED file ==="), verbose)
    log_message(paste("File:", bed_file), verbose)
    
    if (!file.exists(bed_file)) {
        stop(paste(file_type, "BED file does not exist:", bed_file))
    }
    
    tryCatch({
        bed_data <- read.table(bed_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    }, error = function(e) {
        stop("Error reading BED file: ", e$message)
    })
    
    if (ncol(bed_data) < 3) {
        stop("BED file must have at least 3 columns (chr, start, end)")
    }
    
    # Assign column names
    colnames(bed_data)[1:3] <- c("chromosome", "start_position", "end_position")
    
    # Generate site names if not provided (4th column)
    if (ncol(bed_data) >= 4) {
        bed_data$site_name <- bed_data[, 4]
    } else {
        bed_data$site_name <- paste0("site_", 1:nrow(bed_data))
    }
    
    # Create site IDs: chr_start_end
    bed_data$site_id <- paste(bed_data$chromosome, bed_data$start_position, bed_data$end_position, sep = "_")
    
    log_message(paste("Loaded", nrow(bed_data), file_type, "from BED file"), verbose)
    
    # Convert to GRanges object
    site_ranges <- makeGRangesFromDataFrame(bed_data, 
                                          seqnames.field = "chromosome",
                                          start.field = "start_position", 
                                          end.field = "end_position",
                                          keep.extra.columns = TRUE)
    
    # Set proper names using site_id
    names(site_ranges) <- bed_data$site_id
    
    # Apply chromosome filtering to target/reference sites as well
    # This ensures consistency with fragment filtering
    if (file_type == "target sites" || file_type == "reference sites") {
        # Keep only standard chromosomes for consistency
        standard_chr <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")
        site_ranges <- subset(site_ranges, seqnames %in% standard_chr)
        log_message(paste("  Filtered to", length(site_ranges), "sites on standard chromosomes"), verbose)
    }
    
    return(site_ranges)
}

#' Save matrix in both RDS and text format
save_matrix_formats <- function(matrix_data, base_filename, output_dir, save_intermediate, verbose = FALSE) {
    # Always save RDS
    rds_file <- file.path(output_dir, "matrices", paste0(base_filename, ".RDS"))
    saveRDS(matrix_data, rds_file)
    log_message(paste("Saved RDS format:", rds_file), verbose)
    
    # Save as text (always for main outputs, optionally for intermediate)
    if (save_intermediate || base_filename %in% c("raw_counts_matrix", "cpm_matrix", "normalized_matrix")) {
        txt_file <- file.path(output_dir, "matrices", paste0(base_filename, ".txt"))
        
        # Convert matrix to dataframe with row names as first column
        if (is.matrix(matrix_data)) {
            df_to_save <- data.frame(
                peak_coordinates = rownames(matrix_data),
                matrix_data,
                stringsAsFactors = FALSE
            )
        } else {
            df_to_save <- matrix_data
        }
        
        write.table(df_to_save, txt_file, 
                   sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
        log_message(paste("Saved text format:", txt_file), verbose)
    }
}

#' Import BED file with proper chromosome filtering and warning suppression
import_bed_with_filtering <- function(bed_file, include_all_chr = FALSE, verbose = FALSE) {
    # Import with warning suppression
    suppressWarnings({
        ranges <- import(bed_file, format = 'bed')
    })
    
    # Apply chromosome filtering
    if (include_all_chr) {
        # Standard chromosomes only (chr1-chr22, chrX, chrY, chrM)
        standard_chr <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")
        ranges <- subset(ranges, seqnames %in% standard_chr)
    } else {
        # Autosomal chromosomes only (chr1-chr22)
        ranges <- subset(ranges, seqnames %in% paste0("chr", 1:22))
    }
    
    return(ranges)
}

#' Count fragments at specified sites
count_fragments_at_sites <- function(site_ranges, fragment_files, output_dir, 
                                    matrix_name, include_all_chr = FALSE, regenerate = FALSE, 
                                    save_intermediate = FALSE, verbose = FALSE) {
    log_message(paste("=== Counting Fragments at", matrix_name, "==="), verbose)
    
    count_matrix_file <- file.path(output_dir, 'matrices', paste0(matrix_name, '.RDS'))
    
    # Check if count matrix already exists and regeneration is not forced
    if (file.exists(count_matrix_file) && !regenerate) {
        log_message("Loading existing count matrix...", verbose)
        return(readRDS(count_matrix_file))
    }
    
    sample_names <- names(fragment_files)
    
    # Initialize count matrix
    fragment_count_matrix <- matrix(0, 
                                   nrow = length(site_ranges), 
                                   ncol = length(fragment_files))
    
    colnames(fragment_count_matrix) <- sample_names
    rownames(fragment_count_matrix) <- names(site_ranges)
    
    # Log chromosome filtering mode
    if (include_all_chr) {
        log_message("Using standard chromosomes (chr1-chr22, chrX, chrY, chrM)", verbose)
    } else {
        log_message("Using autosomal chromosomes only (chr1-chr22)", verbose)
    }
    
    # Process each fragment file
    for (i in seq_along(sample_names)) {
        sample_name <- sample_names[i]
        current_file <- fragment_files[[sample_name]]
        
        log_message(paste('Processing sample', i, 'of', length(sample_names), ':', sample_name), verbose)
        
        tryCatch({
            # Import fragment file with proper filtering and warning suppression
            fragment_ranges <- import_bed_with_filtering(current_file, include_all_chr, verbose)
            
            log_message(paste("  Retained", length(fragment_ranges), "fragments after chromosome filtering"), verbose)
            
            # Count overlaps with sites (no fragment collapsing)
            suppressWarnings({
                overlap_counts <- countOverlaps(site_ranges, fragment_ranges)
            })
            
            # Store counts in matrix
            fragment_count_matrix[, i] <- overlap_counts
            
            log_message(paste("  Total overlaps:", sum(overlap_counts)), verbose)
            
        }, error = function(e) {
            stop("Error processing fragment file for sample ", sample_name, ": ", e$message)
        })
    }
    
    # Save count matrix
    save_matrix_formats(fragment_count_matrix, matrix_name, output_dir, save_intermediate, verbose)
    
    # Print summary
    total_counts_per_sample <- colSums(fragment_count_matrix)
    log_message("Fragment count summary per sample:", verbose)
    if (verbose) {
        for (i in seq_along(sample_names)) {
            log_message(paste("  ", sample_names[i], ":", total_counts_per_sample[i], "total fragments"), verbose)
        }
    }
    
    return(fragment_count_matrix)
}

#' Generate CPM matrix from raw counts
generate_cpm_matrix <- function(raw_counts, output_dir, save_intermediate = FALSE, verbose = FALSE) {
    log_message("=== Generating CPM Matrix ===", verbose)
    
    # Calculate total reads per sample
    total_reads_per_sample <- colSums(raw_counts)
    
    # Calculate CPM (Counts Per Million)
    cpm_matrix <- sweep(raw_counts, MARGIN = 2, total_reads_per_sample / 1e6, '/')
    
    # Save CPM matrix
    save_matrix_formats(cpm_matrix, "cpm_matrix", output_dir, save_intermediate, verbose)
    
    log_message("CPM matrix generated successfully", verbose)
    
    return(cpm_matrix)
}

#' Generate reference-normalized matrix
generate_normalized_matrix <- function(raw_counts, reference_counts, output_dir, save_intermediate = FALSE, verbose = FALSE) {
    log_message("=== Generating Reference-Normalized Matrix ===", verbose)
    
    # Calculate total reference counts per sample
    total_reference_per_sample <- colSums(reference_counts)
    
    log_message("Reference totals per sample:", verbose)
    if (verbose) {
        for (i in seq_along(total_reference_per_sample)) {
            log_message(paste("  ", names(total_reference_per_sample)[i], ":", total_reference_per_sample[i]), verbose)
        }
    }
    
    # Check for zero reference counts
    if (any(total_reference_per_sample == 0)) {
        zero_samples <- names(total_reference_per_sample)[total_reference_per_sample == 0]
        stop("Zero reference counts found for samples: ", paste(zero_samples, collapse = ", "))
    }
    
    # Normalize: divide each sample's target counts by its total reference counts
    normalized_matrix <- sweep(raw_counts, MARGIN = 2, total_reference_per_sample, '/')
    
    # Save normalized matrix
    save_matrix_formats(normalized_matrix, "normalized_matrix", output_dir, save_intermediate, verbose)
    
    log_message("Reference-normalized matrix generated successfully", verbose)
    
    return(normalized_matrix)
}

#' Generate analysis summary
generate_summary <- function(target_sites, reference_sites, sample_names, output_dir, 
                            analysis_mode, include_all_chr, verbose = FALSE) {
    log_message("=== Generating Analysis Summary ===", verbose)
    
    summary_lines <- c(
        "=============================================================================",
        "CHROMATIN FRAGMENT COUNTER - ANALYSIS SUMMARY",
        "=============================================================================",
        "",
        paste("Analysis completed:", Sys.time()),
        paste("Output directory:", output_dir),
        paste("Analysis mode:", toupper(analysis_mode)),
        paste("Chromosome filtering:", ifelse(include_all_chr, "All chromosomes", "Autosomes only (chr1-chr22)")),
        "",
        "INPUT DATA SUMMARY:",
        paste("- Target sites:", length(target_sites)),
        paste("- Reference sites:", ifelse(!is.null(reference_sites), length(reference_sites), "Not provided")),
        paste("- Samples analyzed:", length(sample_names)),
        "",
        "SAMPLE LIST:",
        paste("  ", 1:length(sample_names), ". ", sample_names, sep = "", collapse = "\n"),
        "",
        "OUTPUT MATRICES GENERATED:",
        "1. Raw counts matrix (raw_counts_matrix.txt/.RDS)",
        "2. CPM matrix (cpm_matrix.txt/.RDS)", 
        ifelse(!is.null(reference_sites), "3. Reference-normalized matrix (normalized_matrix.txt/.RDS)", ""),
        "",
        "MATRIX DIMENSIONS:",
        paste("- Rows (sites):", length(target_sites)),
        paste("- Columns (samples):", length(sample_names)),
        "",
        "=============================================================================",
        "Analysis completed successfully!",
        "============================================================================="
    )
    
    # Write summary to file
    summary_file <- file.path(output_dir, "logs", "analysis_summary.txt")
    writeLines(summary_lines, summary_file)
    
    if (verbose) {
        cat(paste(summary_lines, collapse = "\n"), "\n")
    }
    
    log_message(paste("Analysis summary saved to:", summary_file), verbose)
}

#' Main analysis pipeline function
main_analysis_pipeline <- function() {
    # Parse command line arguments
    params <- parse_arguments()
    verbose_logging <- params$verbose
    
    log_message("=============================================================================", verbose_logging)
    log_message("CHROMATIN FRAGMENT COUNTER - SIMPLIFIED QUANTIFICATION PIPELINE", verbose_logging)
    log_message("=============================================================================", verbose_logging)
    log_message(paste("Analysis mode:", toupper(params$analysis_mode)), verbose_logging)
    log_message(paste("Output directory:", params$output_dir), verbose_logging)
    
    # Create output directory structure
    create_output_structure(params$output_dir, verbose_logging)
    
    # Install and load required packages
    install_required_packages()
    load_required_libraries()
    
    # Handle different analysis modes
    if (params$analysis_mode == "batch") {
        # Batch mode: use samplesheet
        samplesheet <- read_and_validate_samplesheet(params$samplesheet, verbose_logging)
        sample_names <- samplesheet$sample_name
        
        # Find fragment files
        fragment_files <- find_fragment_files(sample_names, params$frags_dir, verbose_logging)
        
    } else {
        # Single sample mode: use direct file inputs
        log_message("=== Single Sample Mode Configuration ===", verbose_logging)
        sample_names <- params$sample_name
        fragment_files <- setNames(list(params$fragment_file), params$sample_name)
        
        log_message(paste("Sample name:", params$sample_name), verbose_logging)
        log_message(paste("Fragment file:", params$fragment_file), verbose_logging)
        
        # Create a minimal samplesheet for consistency
        samplesheet <- data.frame(sample_name = params$sample_name, stringsAsFactors = FALSE)
    }
    
    # Read target sites BED file
    target_sites <- read_bed_file(params$target_sites, "target sites", verbose_logging)
    
    # Read reference sites BED file (if provided)
    reference_sites <- NULL
    if (!is.null(params$reference_sites)) {
        reference_sites <- read_bed_file(params$reference_sites, "reference sites", verbose_logging)
    }
    
    cat("\n")
    log_message("=== STARTING QUANTIFICATION WORKFLOW ===", verbose_logging)
    
    # Step 1: Count fragments at target sites (raw counts matrix)
    raw_counts <- count_fragments_at_sites(
        site_ranges = target_sites,
        fragment_files = fragment_files,
        output_dir = params$output_dir,
        matrix_name = "raw_counts_matrix",
        include_all_chr = params$include_all_chr,
        regenerate = params$regenerate_counts,
        save_intermediate = params$save_intermediate,
        verbose = verbose_logging
    )
    
    # Step 2: Generate CPM matrix
    cpm_matrix <- generate_cpm_matrix(
        raw_counts = raw_counts,
        output_dir = params$output_dir,
        save_intermediate = params$save_intermediate,
        verbose = verbose_logging
    )
    
    # Step 3: Generate normalized matrix (if reference sites provided)
    if (!is.null(reference_sites)) {
        reference_counts <- count_fragments_at_sites(
            site_ranges = reference_sites,
            fragment_files = fragment_files,
            output_dir = params$output_dir,
            matrix_name = "reference_counts_matrix",
            include_all_chr = params$include_all_chr,
            regenerate = params$regenerate_counts,
            save_intermediate = params$save_intermediate,
            verbose = verbose_logging
        )
        
        normalized_matrix <- generate_normalized_matrix(
            raw_counts = raw_counts,
            reference_counts = reference_counts,
            output_dir = params$output_dir,
            save_intermediate = params$save_intermediate,
            verbose = verbose_logging
        )
    }
    
    # Step 4: Generate analysis summary
    generate_summary(
        target_sites = target_sites,
        reference_sites = reference_sites,
        sample_names = sample_names,
        output_dir = params$output_dir,
        analysis_mode = params$analysis_mode,
        include_all_chr = params$include_all_chr,
        verbose = verbose_logging
    )
    
    # Final summary
    cat("\n")
    log_message("=============================================================================", verbose_logging)
    log_message("ANALYSIS COMPLETED SUCCESSFULLY", verbose_logging)
    log_message("=============================================================================", verbose_logging)
    log_message("SUMMARY:", verbose_logging)
    log_message(paste("- Analysis mode:", toupper(params$analysis_mode)), verbose_logging)
    log_message(paste("- Samples processed:", length(sample_names)), verbose_logging)
    log_message(paste("- Target sites quantified:", length(target_sites)), verbose_logging)
    log_message(paste("- Chromosome filtering:", ifelse(params$include_all_chr, "All chromosomes", "Autosomes only")), verbose_logging)
    if (!is.null(reference_sites)) {
        log_message(paste("- Reference sites used:", length(reference_sites)), verbose_logging)
    }
    log_message("", verbose_logging)
    log_message("OUTPUT FILES:", verbose_logging)
    log_message(paste("- Raw counts matrix:", file.path(params$output_dir, "matrices", "raw_counts_matrix.txt")), verbose_logging)
    log_message(paste("- CPM matrix:", file.path(params$output_dir, "matrices", "cpm_matrix.txt")), verbose_logging)
    if (!is.null(reference_sites)) {
        log_message(paste("- Normalized matrix:", file.path(params$output_dir, "matrices", "normalized_matrix.txt")), verbose_logging)
    }
    log_message("=============================================================================", verbose_logging)
}

# Execute main pipeline
if (!interactive()) {
    main_analysis_pipeline()
}
