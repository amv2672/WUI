library(Seurat)
library(Sierra)
library(dplyr)
library(stringr)
library(purrr)
library(tidyr)
library(ggplot2)
library(tibble)


# CONFIGURATION
PAS_FILE   <- "hg38.PAS.main.tsv"   # path to PolyASite TEZ reference

# Cell types to process: names = celltype_harmonyclusters values in Seurat
CELL_TYPE_MAP <- c(
  "Epithelial_cells"  = "epithelial",
  "Fibroblast"        = "Fibroblast",
  "B_plasma_cells"    = "B_plasma_cells",
  "Macrophage"        = "Macrophage",
  "Endothelial_cells" = "Endothelial_cells",
  "Myocytes"          = "Myocytes",
  "Dendritic_cells"   = "Dendritic_cells",
  "T_cells"           = "T_cells",
  "Mast_cells"        = "Mast_cells"
)

MIN_TOTAL_COUNTS   <- 0
MIN_CELLS_PER_PEAK <- 10


#Step 1: Subset overlapping cells

common_cells <- intersect(colnames(genes.seurat), colnames(peaks.seurat))
genes.subset <- genes.seurat[, common_cells]
peaks.seurat.full <- peaks.seurat[, common_cells]
genes.subset <- JoinLayers(genes.subset)


#  Step 2: Filter genes by expression
gene_counts_mat <- genes.subset[["RNA"]]$counts
total_gene_reads <- rowSums(gene_counts_mat)
cells_per_gene   <- rowSums(gene_counts_mat > 0)

active_genes <- names(total_gene_reads[total_gene_reads >= 100 & cells_per_gene >= 50])
message("Active genes retained: ", length(active_genes))


#Step 3: Select UTR3 peaks 

utr3_peak_list <- map(active_genes, function(gene) {
  tryCatch({
    peaks <- Sierra::SelectGenePeaks(peaks.seurat.full, gene = gene, feature.type = "UTR3")
    if (length(peaks) < 2) return(NULL)
    data.frame(Gene = gene, Peak = as.character(peaks), stringsAsFactors = FALSE)
  }, error = function(e) NULL)
})

utr3_peak_df <- bind_rows(compact(utr3_peak_list))

utr3_peak_df <- utr3_peak_df %>%
  mutate(
    Coord  = str_extract(Peak, "chr[^:]+:\\d+-\\d+"),
    Chr    = str_extract(Coord, "chr[^:]+"),
    Start  = as.numeric(str_extract(Coord, "(?<=:)\\d+")),
    End    = as.numeric(str_extract(Coord, "\\d+$")),
    Strand = str_extract(Peak, "[^:]+$")
  )

message("Genes with >=2 UTR3 peaks (before TEZ filter): ",
        n_distinct(utr3_peak_df$Gene))


# Step 3b: Filter peaks to TEZ using hg38.PAS.main.tsv 
# The TEZ (Terminal Exon Zone) is the 3'-most genomic region of the terminal
# exon. PAS sites within the TEZ represent true tandem APA events within the
# same 3'UTR. This comes from PolyA_DB v4.

message("Loading TEZ reference: ", PAS_FILE)
pas <- read.table(PAS_FILE, sep = "\t", header = TRUE,
                  stringsAsFactors = FALSE, quote = "")

# Build one-row-per-gene TEZ boundary lookup
# TEZ_span format: "chr16:+:12568506-12574287"
tez_lookup <- pas %>%
  dplyr::filter(!is.na(TEZ_span), TEZ_span != "") %>%
  dplyr::select(GeneSymbol, TEZ_span) %>%
  dplyr::distinct(GeneSymbol, .keep_all = TRUE) %>%
  dplyr::mutate(
    tez_chr   = str_extract(TEZ_span, "^chr[^:]+"),
    tez_start = as.numeric(str_extract(TEZ_span, "(?<=:)\\d+(?=-)")),
    tez_end   = as.numeric(str_extract(TEZ_span, "\\d+$"))
  ) %>%
  dplyr::rename(Gene = GeneSymbol) %>%
  dplyr::select(Gene, tez_chr, tez_start, tez_end)

message("Genes with TEZ annotation: ", nrow(tez_lookup))

tez_filter <- function(gene_df, tez_tbl) {
  g   <- gene_df$Gene[1]
  tez <- tez_tbl %>% dplyr::filter(Gene == g)

  if (nrow(tez) == 0) return(NULL)   # gene not in PAS database → exclude

  # Keep only peaks whose coordinates fall within the TEZ span
  inside <- gene_df %>%
    dplyr::filter(
      Chr   == tez$tez_chr,
      Start >= tez$tez_start,
      End   <= tez$tez_end
    )

  if (nrow(inside) >= 2) return(inside)
  return(NULL)
}

utr3_filtered_list <- utr3_peak_df %>%
  group_by(Gene) %>%
  group_split() %>%
  map(~ tez_filter(.x, tez_lookup)) %>%
  compact()

utr3_peak_df <- bind_rows(utr3_filtered_list)

message("Genes with >=2 peaks inside TEZ: ",
        n_distinct(utr3_peak_df$Gene))


#  Step 4: Assign weights

master_weights <- utr3_peak_df %>%
  mutate(Peak_Sort = if_else(Strand %in% c("+", "1"), Start, -End)) %>%
  group_by(Gene) %>%
  arrange(Peak_Sort, .by_group = TRUE) %>%
  mutate(
    Num_Peaks = n(),
    Weight    = (row_number() - 1) / (Num_Peaks - 1)
  ) %>%
  ungroup()

write.csv(master_weights, "Master_Peak_Weights_UTR3_tandem_only_Choi_HPVneg_all.csv",
          row.names = FALSE)
message("Master weights saved: ", nrow(master_weights), " peak-gene entries")


# Steps 5–7: Loop through each cell type 

genes_to_run <- unique(master_weights$Gene)
n_genes      <- length(genes_to_run)
message("Genes to process per cell type: ", n_genes)

for (ct_key in names(CELL_TYPE_MAP)) {

  ct_label <- CELL_TYPE_MAP[[ct_key]]
  message("\n", strrep("=", 60))
  message("Processing cell type: ", ct_key, "  (label: ", ct_label, ")")
  message(strrep("=", 60))

  # Step 5: Subset to this cell type 

  peaks.subset <- tryCatch(
    subset(peaks.seurat.full, subset = celltype_harmonyclusters == ct_key),
    error = function(e) {
      message("  Could not subset '", ct_key, "': ", e$message)
      return(NULL)
    }
  )

  if (is.null(peaks.subset) || ncol(peaks.subset) < 50) {
    message("  Skipping — too few cells (",
            if (is.null(peaks.subset)) 0 else ncol(peaks.subset), ")")
    next
  }

  message("  Cells: ", ncol(peaks.subset))

  peaks.subset$Identity <- peaks.subset$Condition
  Idents(peaks.subset) <- "Identity"

  skip_plot     <- 0
  skip_mismatch <- 0
  final_results <- list()

  pb <- txtProgressBar(min = 0, max = n_genes, style = 3)

  for (i in seq_along(genes_to_run)) {

    g <- genes_to_run[i]
    setTxtProgressBar(pb, i)

    g_info  <- master_weights %>% dplyr::filter(Gene == g)
    g_peaks <- g_info$Peak

    plot2 <- tryCatch({
      suppressWarnings(suppressMessages(
        PlotRelativeExpressionBox(peaks.subset, peaks.to.plot = g_peaks)
      ))
    }, error = function(e) NULL)

    if (is.null(plot2) || is.null(plot2$data) || nrow(plot2$data) == 0) {
      skip_plot <- skip_plot + 1
      next
    }

    rel_expr_df <- plot2$data

    # Aggregate per condition
    cond_expr <- rel_expr_df %>%
      group_by(Peak, Identity) %>%
      summarise(Expression = sum(Expression), .groups = "drop") %>%
      pivot_wider(names_from = Identity, values_from = Expression,
                  values_fill = 0)

    expr_mat <- cond_expr %>%
      tibble::column_to_rownames("Peak") %>%
      as.matrix()

    # Match weights
    w <- g_info$Weight[match(rownames(expr_mat), g_info$Peak)]
    if (any(is.na(w))) {
      skip_mismatch <- skip_mismatch + 1
      next
    }

    # Filter conditions by minimum cell 
    peak_cells_ok <- sapply(colnames(expr_mat), function(cond) {
      counts <- rel_expr_df %>%
        dplyr::filter(Peak %in% g_peaks & Identity == cond) %>%
        dplyr::pull(Expression)
      sum(counts > 0) >= MIN_CELLS_PER_PEAK
    })

    if (sum(peak_cells_ok) == 0) next

    expr_mat <- expr_mat[, peak_cells_ok, drop = FALSE]

    # Proportions & WUI
    peak_props <- sweep(expr_mat, 2, colSums(expr_mat), "/")
    wui        <- colSums(peak_props * w)

    for (cond in colnames(expr_mat)) {
      peak_df <- data.frame(
        Gene                      = g,
        Peak                      = rownames(expr_mat),
        Identity                  = cond,
        Expression                = expr_mat[, cond],
        Total_Expression_Identity = sum(expr_mat[, cond]),
        Proportion                = peak_props[, cond],
        Weight                    = w,
        WUI                       = wui[cond],
        stringsAsFactors          = FALSE
      )
      final_results[[length(final_results) + 1]] <- peak_df
    }
  }

  close(pb)
  message("\n  Skipped (no plot data): ", skip_plot)
  message("  Skipped (peak-weight mismatch): ", skip_mismatch)

  if (length(final_results) == 0) {
    message("  No results for ", ct_key, ", skipping output.")
    next
  }

  # Step 6: Build final detailed table 

  final_df <- bind_rows(final_results)

  final_df <- final_df %>%
    dplyr::left_join(
      master_weights %>% dplyr::select(Gene, Peak, Peak_Sort, Weight) %>% dplyr::distinct(),
      by = c("Gene", "Peak")
    ) %>%
    dplyr::mutate(
      Coord  = str_extract(Peak, "chr[^:]+:\\d+-\\d+"),
      Start  = as.numeric(str_extract(Coord, "(?<=:)\\d+")),
      End    = as.numeric(str_extract(Coord, "\\d+$")),
      Strand = str_extract(Peak, "[^:]+$")
    )

  if ("Weight.x" %in% colnames(final_df)) {
    final_df <- final_df %>%
      dplyr::rename(Weight = Weight.x) %>%
      dplyr::select(-Weight.y)
  }

  peak_numbers <- final_df %>%
    dplyr::select(Gene, Peak, Peak_Sort) %>%
    dplyr::distinct() %>%
    dplyr::arrange(Gene, Peak_Sort) %>%
    dplyr::group_by(Gene) %>%
    dplyr::mutate(Peak_Number = row_number()) %>%
    dplyr::ungroup()

  final_df <- final_df %>%
    dplyr::left_join(peak_numbers, by = c("Gene", "Peak", "Peak_Sort"))

  #  Step 7: Save outputs 

  detailed_file <- paste0("dWUI_detailed_peaks_", ct_label, "_HPVneg_tandem_only.csv")
  summary_file  <- paste0("dWUI_summary_",        ct_label, "_HPVneg_tandem_only.csv")
  wide_file     <- paste0("WUI_per_condition_wide_", ct_label, "_HPVneg_tandem_only.csv")

  write.csv(final_df, detailed_file, row.names = FALSE)
  message("  Detailed peaks saved: ", detailed_file)

  dWUI_summary <- final_df %>%
    dplyr::select(Gene, Identity, WUI) %>%
    dplyr::distinct() %>%
    pivot_wider(names_from = Identity, values_from = WUI) %>%
    dplyr::mutate(dWUI = Cancer - Normal)

  write.csv(dWUI_summary, summary_file, row.names = FALSE)
  message("  Summary saved: ", summary_file)
  message("  Genes: ", nrow(dWUI_summary),
          "  Shortening: ", sum(dWUI_summary$dWUI < 0, na.rm = TRUE),
          "  Lengthening: ", sum(dWUI_summary$dWUI > 0, na.rm = TRUE))

  wui_wide <- final_df %>%
    dplyr::select(Gene, Identity, WUI) %>%
    dplyr::distinct() %>%
    pivot_wider(names_from = Identity, values_from = WUI)

  write.csv(wui_wide, wide_file, row.names = FALSE)
  message("  WUI wide saved: ", wide_file)
}

message("\n", strrep("=", 60))
message("All cell types complete.")
message(strrep("=", 60))
