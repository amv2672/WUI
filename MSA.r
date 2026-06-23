# 3'UTR Multi-species MSA for TargetScan — Human / Mouse / Chimpanzee
#
# Filter: human sequence required + at least one of mouse or chimp.
# Output: tab-delimited TargetScan MSA file (gene, taxon_id, aligned_seq).
#

# ── Libraries ──────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(biomaRt)
  library(Biostrings)
  library(rtracklayer)
  library(msa)
  library(seqinr)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(Rsamtools)
})

# ── Parameters ─────────────────────────────────────────────────────────────────
MAX_UTR_LEN <- 10000   # skip genes whose longest species sequence exceeds this

# ── Input / output paths ───────────────────────────────────────────────────────
HUMAN_dPAS_BED <- "~/Downloads/New_Analysis_dWUI/Choi_epi_target_genes_long_utr_pas.bed"
HUMAN_pPAS_BED <- "~/Downloads/New_Analysis_dWUI/Choi_epi_target_genes_short_utr_pas.bed"
OUTPUT_MSA_FILE <- "~/Downloads/New_Analysis_dWUI/3UTR_epi_choihpvneg_MSA_for_TargetScan.txt"

MOUSE_ORTHOLOG_BED <- "/Volumes/LaCie/lab/CellLines_mRNA_Expression/human_to_mouse_exons_hglft_genome_2eaecd_250c00.bed"
CHIMP_ORTHOLOG_BED <- "/Volumes/LaCie/lab/CellLines_mRNA_Expression/human_to_chimp_exons_hglft_genome_e9acb_250a00.bed"
UTR_GTF_PATH       <- "/Volumes/LaCie/lab/scUTRquant/extdata/targets/utrome_hg38_v1/utrome.e30.t5.gc39.pas3.f0.9999.w500.gtf"

HUMAN_FASTA_PATH   <- "/Volumes/LaCie/lab/Downloads/grch38_ncbi_dataset/ncbi_dataset/data/GCA_000001405.15/GCA_000001405.15_GRCh38_genomic.fna"
MOUSE_FASTA_PATH   <- "/Volumes/LaCie/lab/Downloads/mm10_ncbi_dataset/ncbi_dataset/data/GCF_000001635.20/GCF_000001635.20_GRCm38_genomic.fna"
CHIMP_FASTA_PATH   <- "~/Downloads/miRNA/chimp_ncbi_dataset/ncbi_dataset/data/GCF_002880755.1/GCF_002880755.1_Clint_PTRv2_genomic.fna"

# ── 1. Read BED files and build APA 3'UTR ranges ──────────────────────────────
read_pas_bed <- function(path, pas_type) {
  readr::read_delim(path,
    col_names = c("chr", "start", "end", "name", "score", "strand"),
    delim = "\t", show_col_types = FALSE) %>%
    mutate(
      Gene = sub(paste0("_", pas_type, "$"), "", name),
      !!paste0("end.", pas_type) := end
    ) %>%
    dplyr::select(chr, Gene, !!paste0("end.", pas_type), strand)
}

df_dPAS    <- read_pas_bed(HUMAN_dPAS_BED, "dPAS")
df_pPAS    <- read_pas_bed(HUMAN_pPAS_BED, "pPAS")

# Collapse multiple PAS peaks per gene to one representative per gene:
# use the most distal dPAS and most proximal pPAS to define the APA window.
# Without this, inner_join creates all pairwise combos and the max span
# across combos can be 10-100x larger than the true APA region.
df_dPAS <- df_dPAS %>%
  group_by(chr, Gene, strand) %>%
  summarise(end.dPAS = if (strand[1] == "+") max(end.dPAS) else min(end.dPAS),
            .groups = "drop")
df_pPAS <- df_pPAS %>%
  group_by(chr, Gene, strand) %>%
  summarise(end.pPAS = if (strand[1] == "+") min(end.pPAS) else max(end.pPAS),
            .groups = "drop")

df_combined <- inner_join(df_dPAS, df_pPAS, by = c("chr", "Gene", "strand"))

# ── 2. UTRome: get spliced exons within each APA window ───────────────────────
# The pPAS→dPAS genomic span often contains introns (e.g. CROT: 39kb genomic,
# ~3.8kb exonic). The UTRome GTF has the exon structure; we collect all UTRome
# exons overlapping [pPAS+1, dPAS], merge them, and concatenate for mRNA sequence.

message("Importing UTRome GTF...")
utrome_exons <- import(UTR_GTF_PATH)
utrome_exons <- utrome_exons[utrome_exons$type == "exon"]
message("UTRome exons loaded: ", length(utrome_exons))

get_apa_exons_utrome <- function(gene, chr_name, strand_val, pPAS_end, dPAS_end) {
  g_ex <- utrome_exons[!is.na(utrome_exons$gene_name) & utrome_exons$gene_name == gene &
                        as.character(seqnames(utrome_exons)) == chr_name]
  if (length(g_ex) == 0) return(NULL)

  # APA genomic window (same orientation logic as before)
  if (strand_val == "+") {
    apa_start <- pPAS_end + 1; apa_end <- dPAS_end
  } else {
    apa_start <- min(dPAS_end, pPAS_end - 1)
    apa_end   <- max(dPAS_end, pPAS_end - 1)
  }
  if (apa_start >= apa_end) return(NULL)

  # Keep exons overlapping the window, trim to window boundaries
  g_ex <- g_ex[end(g_ex) >= apa_start & start(g_ex) <= apa_end]
  if (length(g_ex) == 0) return(NULL)
  start(g_ex) <- pmax(start(g_ex), apa_start)
  end(g_ex)   <- pmin(end(g_ex),   apa_end)
  g_ex <- g_ex[width(g_ex) > 0]
  if (length(g_ex) == 0) return(NULL)

  # Merge overlapping exons from different isoforms
  g_red <- reduce(g_ex)
  strand(g_red) <- strand_val
  mcols(g_red)$Gene <- gene
  mcols(g_red)$name <- gene
  g_red
}

apa_exon_list <- mapply(
  get_apa_exons_utrome,
  gene       = df_combined$Gene,
  chr_name   = df_combined$chr,
  strand_val = df_combined$strand,
  pPAS_end   = df_combined$end.pPAS,
  dPAS_end   = df_combined$end.dPAS,
  SIMPLIFY   = FALSE
)
apa_exon_list <- Filter(Negate(is.null), apa_exon_list)
message("Genes with UTRome APA exons: ", length(apa_exon_list))

# Export exon-level BED — use this for mouse/chimp LiftOver for best accuracy
apa_exon_gr <- unlist(GRangesList(apa_exon_list), use.names = FALSE)
export(apa_exon_gr, "human_apa_exons_for_liftover.bed")
message("Exon-level BED saved: human_apa_exons_for_liftover.bed")
message("Tip: re-run LiftOver on this BED (not gene-level) for exact mouse/chimp exon coords")

# ── 3. BioMart: orthologs for mouse and chimp ─────────────────────────────────
target_genes <- unique(names(apa_exon_list))
target_genes <- target_genes[nzchar(target_genes) & !is.na(target_genes)]

human_mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

orthologs_raw <- getBM(
  attributes = c(
    "external_gene_name",
    "mmusculus_homolog_associated_gene_name",
    "ptroglodytes_homolog_associated_gene_name",
    "strand"
  ),
  filters = "external_gene_name",
  values  = target_genes,
  mart    = human_mart
)
if (nrow(orthologs_raw) == 0) stop("BioMart returned no orthologs — check gene symbols.")

orthologs <- orthologs_raw %>%
  rename(
    Human      = external_gene_name,
    Mouse      = mmusculus_homolog_associated_gene_name,
    Chimp      = ptroglodytes_homolog_associated_gene_name,
    Strand_raw = strand
  ) %>%
  mutate(Human_Strand = ifelse(Strand_raw == 1, "+", "-")) %>%
  dplyr::select(-Strand_raw) %>%
  filter(nzchar(Mouse) | nzchar(Chimp))   # keep if at least one ortholog

message("Genes with ≥1 ortholog (mouse or chimp): ", length(unique(orthologs$Human)))

# Filter apa_exon_list to genes with ≥1 ortholog
genes_keep    <- unique(orthologs$Human)
apa_exon_list <- apa_exon_list[names(apa_exon_list) %in% genes_keep]
message("Genes after ortholog filter: ", length(apa_exon_list))

# ── 4. Helpers: open indexed FASTA and harmonize seqlevels ────────────────────
open_fasta <- function(path) {
  fa <- FaFile(path); open(fa); fa
}

# scanFaIndex returns only bare accessions (CM000663.2, NC_000067.6) — no
# descriptions. For human GCA assembly the accessions are GenBank CM series,
# which UCSC maps to RefSeq NC series — completely different. Hard-code the
# GRCh38 GCA map. Mouse/chimp are GCF (NC_ series) so UCSC lookup works.
get_chr_accession_map <- function(ucsc_genome) {
  switch(ucsc_genome,
    # GCA_000001405.15 — GenBank CM accessions (not NC)
    hg38 = c(
      chr1  = "CM000663.2", chr2  = "CM000664.1", chr3  = "CM000665.1",
      chr4  = "CM000666.1", chr5  = "CM000667.1", chr6  = "CM000668.1",
      chr7  = "CM000669.1", chr8  = "CM000670.1", chr9  = "CM000671.1",
      chr10 = "CM000672.1", chr11 = "CM000673.1", chr12 = "CM000674.1",
      chr13 = "CM000675.1", chr14 = "CM000676.1", chr15 = "CM000677.1",
      chr16 = "CM000678.1", chr17 = "CM000679.1", chr18 = "CM000680.1",
      chr19 = "CM000681.1", chr20 = "CM000682.1", chr21 = "CM000683.1",
      chr22 = "CM000684.1", chrX  = "CM000685.1", chrY  = "CM000686.1",
      chrM  = "J01415.2"
    ),
    # GCF_000001635.20 — mm10 / GRCm38
    mm10 = c(
      chr1  = "NC_000067.6", chr2  = "NC_000068.7", chr3  = "NC_000069.6",
      chr4  = "NC_000070.6", chr5  = "NC_000071.6", chr6  = "NC_000072.6",
      chr7  = "NC_000073.6", chr8  = "NC_000074.6", chr9  = "NC_000075.6",
      chr10 = "NC_000076.6", chr11 = "NC_000077.6", chr12 = "NC_000078.6",
      chr13 = "NC_000079.6", chr14 = "NC_000080.6", chr15 = "NC_000081.6",
      chr16 = "NC_000082.6", chr17 = "NC_000083.6", chr18 = "NC_000084.6",
      chr19 = "NC_000085.6", chrX  = "NC_000086.7", chrY  = "NC_000087.7",
      chrM  = "NC_005089.1"
    ),
    # GCF_002880755.1 — Clint_PTRv2 / panTro5
    panTro5 = c(
      chr1  = "NC_036879.1", chr2A = "NC_036880.1", chr2B = "NC_036881.1",
      chr3  = "NC_036882.1", chr4  = "NC_036883.1", chr5  = "NC_036884.1",
      chr6  = "NC_036885.1", chr7  = "NC_036886.1", chr8  = "NC_036887.1",
      chr9  = "NC_036888.1", chr10 = "NC_036889.1", chr11 = "NC_036890.1",
      chr12 = "NC_036891.1", chr13 = "NC_036892.1", chr14 = "NC_036893.1",
      chr15 = "NC_036894.1", chr16 = "NC_036895.1", chr17 = "NC_036896.1",
      chr18 = "NC_036897.1", chr19 = "NC_036898.1", chr20 = "NC_036899.1",
      chr21 = "NC_036900.1", chr22 = "NC_036901.1", chrX  = "NC_036902.1",
      chrY  = "NC_036903.1", chrM  = "NC_001643.1"
    ),
    stop("Unknown genome: '", ucsc_genome, "'. Add it to get_chr_accession_map().")
  )
}

harmonize_bed <- function(bed_gr, fa, species, ucsc_genome) {
  fa_accessions       <- sub(" .*", "", seqlevels(scanFaIndex(fa)))
  fa_acc_unver        <- sub("\\..*", "", fa_accessions)
  chr_map             <- get_chr_accession_map(ucsc_genome)
  chr_map_unver       <- sub("\\..*", "", chr_map)

  valid   <- seqlevels(bed_gr) %in% names(chr_map) &
             chr_map_unver[seqlevels(bed_gr)] %in% fa_acc_unver
  present <- seqlevels(bed_gr)[valid]
  if (length(present) == 0)
    stop(species, ": no chr names matched FASTA accessions. ",
         "Run head(seqlevels(scanFaIndex(fa)), 5) to inspect FASTA accessions.")

  # Use the versioned accession that's actually in the FASTA
  target_acc <- fa_accessions[match(chr_map_unver[present], fa_acc_unver)]
  bed_gr     <- keepSeqlevels(bed_gr, present, pruning.mode = "coarse")
  seqlevels(bed_gr) <- target_acc
  message(species, " ranges after harmonization: ", length(bed_gr))
  bed_gr
}

extract_seqs <- function(bed_gr, fa, species) {
  seqs <- getSeq(fa, bed_gr)
  names(seqs) <- mcols(bed_gr)$name
  seqs <- seqs[width(seqs) > 0]
  message(species, " sequences extracted: ", length(seqs))
  seqs
}

# ── 5. Human sequences — UTRome exon extraction ───────────────────────────────
human_fa <- open_fasta(HUMAN_FASTA_PATH)

# Pre-build the chr → versioned FASTA accession lookup once
chr_to_acc   <- get_chr_accession_map("hg38")
fa_acc_all   <- sub(" .*", "", seqlevels(scanFaIndex(human_fa)))
fa_acc_unver <- sub("\\..*", "", fa_acc_all)
chr_unver    <- sub("\\..*", "", chr_to_acc)

human_utr <- DNAStringSet(sapply(names(apa_exon_list), function(g) {
  exs <- apa_exon_list[[g]]

  # Map UTRome chr names → FASTA accessions
  lvls    <- seqlevels(exs)
  new_acc <- fa_acc_all[match(chr_unver[lvls], fa_acc_unver)]
  if (any(is.na(new_acc))) {
    lvls    <- lvls[!is.na(new_acc)]
    new_acc <- new_acc[!is.na(new_acc)]
    exs     <- keepSeqlevels(exs, lvls, pruning.mode = "coarse")
  }
  if (length(exs) == 0) return(NA_character_)
  seqlevels(exs) <- new_acc

  # Sort exons in mRNA order: ascending (+), descending (-)
  exs <- sort(exs, decreasing = as.character(strand(exs[1])) == "-")

  paste(as.character(getSeq(human_fa, exs)), collapse = "")
}))

human_utr <- human_utr[!is.na(as.character(human_utr)) & nchar(as.character(human_utr)) > 0]
message("Human sequences extracted: ", length(human_utr))

# ── 6. Mouse sequences ────────────────────────────────────────────────────────
# Using gene-level liftover BED. For genes with intronic APA regions (e.g. CROT),
# re-run LiftOver on human_apa_exons_for_liftover.bed and update MOUSE_ORTHOLOG_BED.
mouse_fa     <- open_fasta(MOUSE_FASTA_PATH)
mouse_bed_gr <- harmonize_bed(import(MOUSE_ORTHOLOG_BED), mouse_fa, "Mouse", "mm10")
mouse_utr    <- extract_seqs(mouse_bed_gr, mouse_fa, "Mouse")

# ── 7. Chimp sequences ────────────────────────────────────────────────────────
chimp_fa     <- open_fasta(CHIMP_FASTA_PATH)
chimp_bed_gr <- harmonize_bed(import(CHIMP_ORTHOLOG_BED), chimp_fa, "Chimp", "panTro5")
chimp_utr    <- extract_seqs(chimp_bed_gr, chimp_fa, "Chimp")

# ── 8. Filter: keep genes with human + ≥1 other species ──────────────────────
human_genes <- names(human_utr)
genes_pass  <- human_genes[human_genes %in% names(mouse_utr) |
                            human_genes %in% names(chimp_utr)]
message("Genes passing human + ≥1 species filter: ", length(genes_pass))

human_utr <- human_utr[genes_pass]
mouse_utr <- mouse_utr[intersect(genes_pass, names(mouse_utr))]
chimp_utr <- chimp_utr[intersect(genes_pass, names(chimp_utr))]

# ── 9. MSA per gene ───────────────────────────────────────────────────────────
species_id_map <- c(Human = 9606L, Mouse = 10090L, Chimp = 9598L)

result_list <- lapply(genes_pass, function(g) {
  human_seq <- as.character(human_utr[[g]])
  human_len <- nchar(human_seq)

  # Filter on human length only — mouse/chimp syntenic regions can be larger
  # due to species-specific insertions; the aligner handles length differences.
  if (human_len > MAX_UTR_LEN) {
    message("Skipping ", g, " (human ", human_len, " bp > MAX_UTR_LEN)")
    return(NULL)
  }

  seqs <- c(Human = human_seq)
  if (g %in% names(mouse_utr)) seqs["Mouse"] <- as.character(mouse_utr[[g]])
  if (g %in% names(chimp_utr)) seqs["Chimp"] <- as.character(chimp_utr[[g]])

  if (length(seqs) == 1) return(NULL)   # should not happen after filter

  aln      <- msa(DNAStringSet(seqs), method = "ClustalW", order = "input")
  aln_seqs <- as.character(msaConvert(aln, type = "seqinr::alignment")$seq)
  names(aln_seqs) <- names(seqs)

  row <- data.frame(Gene = g, stringsAsFactors = FALSE)
  for (sp in names(seqs)) row[[sp]] <- aln_seqs[[sp]]
  row
})

result_list <- Filter(Negate(is.null), result_list)
final_df    <- bind_rows(result_list)
message("MSA complete: ", nrow(final_df), " genes")

# ── 10. Reshape to TargetScan format ─────────────────────────────────────────
species_cols <- intersect(c("Human", "Mouse", "Chimp"), colnames(final_df))

targetscan_long <- final_df %>%
  pivot_longer(cols = all_of(species_cols),
               names_to  = "Species",
               values_to = "Aligned_Sequence") %>%
  filter(!is.na(Aligned_Sequence)) %>%
  mutate(Species_ID = species_id_map[Species]) %>%
  dplyr::select(Gene_Symbol = Gene, Species_ID, Aligned_Sequence)

# Sanity check: no all-gap sequences
all_gap <- targetscan_long %>% filter(grepl("^-+$", Aligned_Sequence))
if (nrow(all_gap) > 0)
  warning(nrow(all_gap), " all-gap sequences found — review alignment quality.")

write_tsv(targetscan_long, OUTPUT_MSA_FILE, col_names = FALSE)
message("Saved: ", OUTPUT_MSA_FILE)
message("Rows: ", nrow(targetscan_long),
        "  (", nrow(final_df), " genes × up to 3 species)")
head(targetscan_long)
