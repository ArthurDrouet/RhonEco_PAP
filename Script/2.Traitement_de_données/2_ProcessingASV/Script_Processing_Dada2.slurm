#!/usr/bin/env Rscript
 
#### PROCESSING ASV with Dada2 PIPELINE ####
# Date :   11/2026
# Lieu :   Lyon, France
# Projet : Piège à sédiments métabaroding
# Caractèrisation de diverisité multitaxonomique (Diatomées, Insectes, Poissons) dans les sédiments capturé par des pièges à sédiments dans le Rhône.

# Auteur : Arthur DROUET
# Equipe : Johan PANSU (LEHNA, CNRS), Arthur DROUET (RIVERLY, INRAE)
# Titre :  Preprocess Script

# Objectif :
# Utiliser Dada2 pour nettoyer les séquences et obetenur des ASVs.

# Steps of Processing :
# 1. Paths and working directory
# 2. Sample names
# 3. Filter & trim
# 4. Learning errors + plots
# 5. Denoising (dada)
# 6. Merge pairs + sequence table
# 7. Remove chimeras
# 8. Track reads
# 9. Export
# 10. Save Image


### 1. Paths and working directory ####
cat("1. Paths and working directory\n")

# Load directory
setwd("/beegfs/data/adrouet/Stage_Lehna/MicroLones_Git")

WORKDIR   <- "Results/ASV"
SCRIPTDIR <- "/beegfs/home/adrouet/Stage_Lehna/MicroLones_Git/Scripts/ASV/2_ProcessingASV"
path      <- file.path(WORKDIR, "Demultiplex", "1_Preprocessing", "RemovePrimer")
outdir    <- file.path(WORKDIR, "Demultiplex", "2_ProcessingASV_250180_24", "Dada2")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "Filt"),   recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "Obj"),    recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "Figure"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "Output"),   recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "Log"),   recursive = TRUE, showWarnings = FALSE)

# LOG 
LOGFILE <- file.path(outdir, "Log", "Dada2.log")
sink(LOGFILE, split = TRUE)
on.exit({
  sink()
}, add = TRUE)

# Optional: save/load the full R session (useful as a checkpoint)
# save.image(file = file.path(SCRIPTDIR, "Script_Dada2.RData"))
# RDATA <- file.path(SCRIPTDIR, "Script_Dada2.RData")

# Load Package Dada2
library(dada2)
cat("dada2 version:", as.character(packageVersion("dada2")), "\n")
library(ggplot2)
cat("ggplot2 version:", as.character(packageVersion("ggplot2")), "\n")
library(ShortRead)
cat("ShortRead version:", as.character(packageVersion("ShortRead")), "\n")

# Fonction
############################################################################################
############################################################################################
plot_EE <- function(fq_paths, out_pdf, title) {
  EE_all <- numeric()
  cat("Computing EE:", title, "\n")
  for (p in fq_paths) {
    fq <- readFastq(p)
    q  <- as(quality(fq), "matrix")
    EE <- rowSums(10^(-q/10))
    EE_all <- c(EE_all, EE)
    rm(fq, q, EE)
    gc()
  }

  pdf(out_pdf, width = 8, height = 6)
  h <- hist(EE_all, breaks = 200, plot = FALSE)
  plot(h,
       main = paste0(title, "\nN = ", length(EE_all)),
       xlab = "Expected Errors",
       ylab = "Reads",
       xlim = c(0, 10),
       col = "grey80",
       border = "grey30")
  abline(v = 1:5, lty = 2, col = "red")
  dev.off()
  cat("Saved:", out_pdf, "\n")
}
############################################################################################
############################################################################################


### 2. Sample names & Errors diagnostic ####
cat("2. Sample names & Errors diagnostic\n")

fnFs <- sort(list.files(path, pattern = "_R1_RemovePrimer\\.fastq(\\.gz)?$", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2_RemovePrimer\\.fastq(\\.gz)?$", full.names = TRUE))

cat("Forward:", length(fnFs), " | Reverse:", length(fnRs), "\n")

sample.names <- sub("_R1_RemovePrimer\\.fastq(\\.gz)?$", "", basename(fnFs))
cat("Samples detected:", length(sample.names), "\n")

RUN_FULL <- FALSE
if (RUN_FULL) {

RUN_FULL <- FALSE
if (RUN_FULL) {

#Quality profile
# Forward
pdf(file.path(outdir, "Figure", "QualityProfile_Forward_files.pdf"),
    width = 12, height = 8)
plotQualityProfile(fnFs[1:6])
dev.off()
cat("Saved: QualityProfile_Forward_files.pdf\n")

# Reverse
pdf(file.path(outdir, "Figure", "QualityProfile_Reverse_files.pdf"),
    width = 12, height = 8)
plotQualityProfile(fnRs[1:6])
dev.off()
cat("Saved: QualityProfile_Reverse_files.pdf\n")

# EE distribution GLOBAL — Forward (all samples pooled)
plot_EE(fnFs,
        file.path(outdir, "Figure", "EE_Forward_PREFILT.pdf"),
        "Forward EE PRE FILTER")

# EE distribution GLOBAL — Reverse (all samples pooled)
plot_EE(fnRs,
        file.path(outdir, "Figure", "EE_Reverse_PREFILT.pdf"),
        "Reverse EE PRE FILTER")

}

### 3. Filter & trim ####
cat("3. Filter & trim\n")

filtFs <- file.path(outdir, "Filt", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(outdir, "Filt", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(
  fnFs, filtFs,
  fnRs, filtRs,
  truncLen = c(250, 180),
  maxN = 0,
  maxEE = c(2, 4),
  truncQ = 2,
  rm.phix = TRUE,
  compress = TRUE,
  multithread = TRUE
)

print(head(out))
saveRDS(out, file.path(outdir, "Obj", "filterAndTrim_out.rds"))
write.csv(out, file.path(outdir, "Output", "post.filt.csv"), quote = FALSE)

# Keep only samples with reads after filtering (2nd column) AND existing files
keep_idx <- out[, 2] > 0
filtFs <- filtFs[keep_idx]
filtRs <- filtRs[keep_idx]
sample.names <- sample.names[keep_idx]

keep_files <- file.exists(filtFs) & file.exists(filtRs)
filtFs <- filtFs[keep_files]
filtRs <- filtRs[keep_files]
sample.names <- sample.names[keep_files]

cat("Samples kept:", length(sample.names), "\n")


### 4. Learn errors ####
cat("4. Learn errors\n")

errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

saveRDS(errF, file.path(outdir, "Obj", "errF.rds"))
saveRDS(errR, file.path(outdir, "Obj", "errR.rds"))

pdf(file.path(outdir, "Figure", "DADA2_ErrorModels.pdf"), width = 8, height = 6)
plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)
dev.off()


### 5. Denoise ####
cat("5. Denoise\n")

dadaFs <- dada(filtFs, err = errF, multithread = TRUE)
dadaRs <- dada(filtRs, err = errR, multithread = TRUE)

print(dadaFs[[1]])

### 6. Merge + seqtab + fasta ####
cat("6. Merge + seqtab + fasta\n")

mergers <- mergePairs(
  dadaFs, filtFs, 
  dadaRs, filtRs, 
  minOverlap  = 12,   # 12 : default
  maxMismatch = 0, # 0 : default
  verbose = TRUE)

print(head(mergers[[1]]))

seqtab <- makeSequenceTable(mergers)

saveRDS(seqtab, file.path(outdir, "Obj", "seqtab.rds"))
write.csv(seqtab, file.path(outdir, "Output", "SeqTable.csv"), quote = FALSE)

# Follow ASV
asv_seqs <- colnames(seqtab)
asv_ids <- paste0("ASV", seq_len(ncol(seqtab)))

# Build mapping 
asv_map <- data.frame(ASV_ID = asv_ids, Sequence = asv_seqs, stringsAsFactors = FALSE)
write.table(asv_map, file.path(outdir, "Output", "ASV_map.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

# Build fasta
fasta_lines <- c(rbind(paste0(">", asv_ids), asv_seqs))
writeLines(fasta_lines, file.path(outdir, "Output", "ASV.fasta"))
seqtab_ID <- seqtab

# Build table with ASV_ID
colnames(seqtab_ID) <- asv_ids
saveRDS(seqtab_ID, file.path(outdir, "Obj", "seqtab_ASVID.rds"))
write.table(seqtab_ID,
            file = file.path(outdir, "Output", "SeqTable_ASVID.tsv"),
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)



}

########################  ################## a supr 
# ---- 1) Reload seqtab ----
seqtab_path <- file.path(outdir, "Obj", "seqtab.rds")
stopifnot(file.exists(seqtab_path))

seqtab <- readRDS(seqtab_path)

cat(
  "Reloaded seqtab:",
  nrow(seqtab), "samples x",
  ncol(seqtab), "ASVs\n"
)

# ---- 2) Rebuild ASV sequences / IDs ----
asv_seqs <- colnames(seqtab)
asv_ids  <- paste0("ASV", seq_len(ncol(seqtab)))
########################  ####################### jusque la 

### Statistics & plots for SeqTab ####
cat("### Statistics & plots for SeqTab ###\n")

# Total reads + number of ASVs
cat("Total reads: ", sum(seqtab), "\n", sep = "")
cat("Number of ASVs: ", ncol(seqtab), "\n", sep = "")

# Reads per sample
cat("Reads per sample (summary):\n")
reads_per_sample <- rowSums(seqtab)
print(summary(reads_per_sample))

df_sample <- data.frame(
  sample = rownames(seqtab),
  reads  = as.numeric(reads_per_sample),
  stringsAsFactors = FALSE
)

p_sample <- ggplot(df_sample, aes(x = "", y = reads)) +
  geom_violin(fill = "grey80") +
  geom_boxplot(width = 0.15, outlier.size = 0.6) +
  scale_y_log10() +
  labs(
    title = "Distribution of read counts per sample",
    y = "Read count (log10)",
    x = ""
  ) +
  theme_bw()

print(p_sample)
ggsave(
  filename = file.path(outdir, "Figure", "SeqTab_ReadsPerSample_ViolinBox.png"),
  plot = p_sample, width = 7, height = 5, dpi = 300
)

# Reads per ASV
cat("Reads per ASV (summary):\n")
reads_per_asv <- colSums(seqtab)
print(summary(reads_per_asv))

df_asv <- data.frame(
  ASV   = asv_ids,
  reads = as.numeric(reads_per_asv),
  stringsAsFactors = FALSE
)

p_asv_violin <- ggplot(df_asv, aes(x = "", y = reads)) +
  geom_violin(fill = "grey80") +
  geom_boxplot(width = 0.15, outlier.size = 0.6) +
  scale_y_log10() +
  labs(
    title = "Distribution of ASV abundances",
    y = "Reads per ASV (log10)",
    x = ""
  ) +
  theme_bw()

print(p_asv_violin)
ggsave(
  filename = file.path(outdir, "Figure", "SeqTab_ReadsPerASV_ViolinBox.png"),
  plot = p_asv_violin, width = 7, height = 5, dpi = 300
)

p_asv_hist <- ggplot(df_asv, aes(x = reads)) +
  geom_histogram(bins = 100, fill = "grey70", color = "black") +
  scale_x_log10() +
  labs(
    title = "Histogram of reads per ASV",
    x = "Reads per ASV (log10)",
    y = "Number of ASVs"
  ) +
  theme_bw()

print(p_asv_hist)
ggsave(
  filename = file.path(outdir, "Figure", "SeqTab_ReadsPerASV_Histogram.png"),
  plot = p_asv_hist, width = 7, height = 5, dpi = 300
)

# Singletons
cat("Singleton ASVs (reads == 1): ", sum(reads_per_asv == 1), "\n", sep = "")

# ASV length
cat("ASV sequence length (bp) summary:\n")
asv_length <- nchar(asv_seqs)
print(summary(asv_length))

df_len <- data.frame(
  ASV    = asv_ids,
  length = as.numeric(asv_length),
  stringsAsFactors = FALSE
)

p_len_hist <- ggplot(df_len, aes(x = length)) +
  geom_histogram(bins = 50, fill = "grey70", color = "black") +
  labs(
    title = "Distribution of ASV sequence lengths",
    x = "Sequence length (bp)",
    y = "Number of ASVs"
  ) +
  theme_bw()

print(p_len_hist)
ggsave(
  filename = file.path(outdir, "Figure", "SeqTab_ASV_Length_Histogram.png"),
  plot = p_len_hist, width = 7, height = 5, dpi = 300
)


### 7. Remove chimeras ####
cat("7. Remove chimeras\n")

seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE)

saveRDS(seqtab.nochim, file.path(outdir, "Obj", "seqtab.nochim.rds"))
write.csv(seqtab.nochim, file.path(outdir, "Output", "SeqTable_NoChim.csv"), quote = FALSE)

# Build fasta NoChim
asv_seqs_nochim <- colnames(seqtab.nochim)
asv_ids_nochim <- asv_ids[match(asv_seqs_nochim, asv_seqs)] # Reprendre les ASV_IDs de seqtab et filtré

# Build mapping NoChim
asv_map_nochim <- data.frame(ASV_ID = asv_ids_nochim, Sequence = asv_seqs_nochim, stringsAsFactors = FALSE)
write.table(asv_map_nochim, file.path(outdir, "Output", "ASV_map_nochim.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

# Build fasta NoChim
fasta_lines_nochim <- c(rbind(paste0(">", asv_ids_nochim), asv_seqs_nochim))
writeLines(fasta_lines_nochim, file.path(outdir, "Output", "ASV_nochim.fasta"))

# Build table NoChim with ASV_ID (same to SeqTab)
seqtab.nochim_ID <- seqtab.nochim
colnames(seqtab.nochim_ID) <- asv_ids_nochim
saveRDS(seqtab.nochim_ID, file.path(outdir, "Obj", "seqtab.nochim_ASVID.rds"))
write.table(seqtab.nochim_ID,
            file = file.path(outdir, "Output", "SeqTable_NoChim_ASVID.tsv"),
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

### Statistics & plots for SeqTab_NoChim ####
cat("### Statistics & plots for SeqTab_NoChim ###\n")

# Total reads + number of ASVs
cat("Total reads: ", sum(seqtab.nochim), "\n", sep = "")
cat("Number of ASVs: ", ncol(seqtab.nochim), "\n", sep = "")

# Reads per sample
cat("Reads per sample (summary):\n")
reads_per_sample <- rowSums(seqtab.nochim)
print(summary(reads_per_sample))

# Reads per ASV
cat("Reads per ASV (summary):\n")
reads_per_asv <- colSums(seqtab.nochim)
print(summary(reads_per_asv))

# Singletons
cat("Singleton ASVs (reads == 1): ", sum(reads_per_asv == 1), "\n", sep = "")

# Why is removed
asv_removed_seqs <- setdiff(asv_seqs, asv_seqs_nochim)
asv_removed_ids  <- asv_ids[match(asv_removed_seqs, asv_seqs)]
cat("ASV lost (chimeras removed): ", length(asv_removed_ids), "\n", sep = "")


### 8. Track reads through the pipeline ####
cat("8. Track reads through the pipeline\n")
getN <- function(x) sum(getUniques(x))

# Align 'out' rows to the kept samples (out rownames are input file paths)
out_kept <- out[keep_idx, , drop = FALSE]
out_kept <- out_kept[keep_files, , drop = FALSE]
rownames(out_kept) <- sample.names


track <- cbind(
  input     = out_kept[, 1],
  filtered  = out_kept[, 2],
  denoisedF = sapply(dadaFs, getN),
  denoisedR = sapply(dadaRs, getN),
  merged    = sapply(mergers, getN),
  nonchim   = rowSums(seqtab.nochim)
)

write.csv(track, file.path(outdir, "Output", "track.csv"), quote = FALSE)

# Resume Track
col_sums <- colSums(track)
cat(paste(names(col_sums), col_sums, sep = "\t", collapse = "\n"), "\n")


### 9. Save Image (optional) ####
cat("9. Save Image\n")

save.image(file = file.path(SCRIPTDIR, "Script_ASV.RData"))





























RUN_FULL <- FALSE
if (RUN_FULL) {



### 3. Test Filter & trim — benchmark de paramètres ####
cat("3. Test Filter & trim — parameter benchmarking\n")

TESTROOT <- file.path(outdir, "Filt_param_tests")
dir.create(TESTROOT, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "Obj"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "Output"), recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# TEST01 — full parameters (reference)
# ----------------------------
cat("\nTEST01 — full parameters (reference)\n")

TESTDIR <- file.path(TESTROOT, "T01_full")
dir.create(TESTDIR, recursive = TRUE, showWarnings = FALSE)

filtFs <- file.path(TESTDIR, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(TESTDIR, paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(
  fnFs, filtFs,
  fnRs, filtRs,
  truncLen = c(250, 200),
  maxN     = 0,
  maxEE    = c(2, 2),
  truncQ   = 2,
  rm.phix  = TRUE,
  compress = TRUE,
  multithread = TRUE
)

saveRDS(out, file.path(outdir, "Obj", "out_T01_full.rds"))
write.csv(out, file.path(outdir, "Output", "post.filt_T01_full.csv"), quote = FALSE)
cat("T01 total in/out:", sum(out[, 1]), sum(out[, 2]), "\n")

# ----------------------------
# TEST02 — default (no extra params)
# ----------------------------
cat("\nTEST02 — default (no extra params)\n")

TESTDIR <- file.path(TESTROOT, "T02_default")
dir.create(TESTDIR, recursive = TRUE, showWarnings = FALSE)

filtFs <- file.path(TESTDIR, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(TESTDIR, paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(
  fnFs, filtFs,
  fnRs, filtRs,
  compress = TRUE,
  multithread = TRUE
)

saveRDS(out, file.path(outdir, "Obj", "out_T02_default.rds"))
write.csv(out, file.path(outdir, "Output", "post.filt_T02_default.csv"), quote = FALSE)
cat("T02 total in/out:", sum(out[, 1]), sum(out[, 2]), "\n")

# ----------------------------
# TEST03 — rm.phix
# ----------------------------
cat("\nTEST03 — rm.phix=TRUE\n")

TESTDIR <- file.path(TESTROOT, "T03_rmPhix_TRUE")
dir.create(TESTDIR, recursive = TRUE, showWarnings = FALSE)

filtFs <- file.path(TESTDIR, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(TESTDIR, paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(
  fnFs, filtFs,
  fnRs, filtRs,
  rm.phix = TRUE,
  compress = TRUE,
  multithread = TRUE
)

saveRDS(out, file.path(outdir, "Obj", "out_T03_rmPhix_TRUE.rds"))
write.csv(out, file.path(outdir, "Output", "post.filt_T03_rmPhix_TRUE.csv"), quote = FALSE)
cat("T03 total in/out:", sum(out[, 1]), sum(out[, 2]), "\n")

# ----------------------------
# TEST04 — truncQ=2
# ----------------------------
cat("\nTEST04 — truncQ=2\n")

TESTDIR <- file.path(TESTROOT, "T04_truncQ_2")
dir.create(TESTDIR, recursive = TRUE, showWarnings = FALSE)

filtFs <- file.path(TESTDIR, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(TESTDIR, paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(
  fnFs, filtFs,
  fnRs, filtRs,
  truncQ = 2,
  compress = TRUE,
  multithread = TRUE
)

saveRDS(out, file.path(outdir, "Obj", "out_T04_truncQ_2.rds"))
write.csv(out, file.path(outdir, "Output", "post.filt_T04_truncQ_2.csv"), quote = FALSE)
cat("T04 total in/out:", sum(out[, 1]), sum(out[, 2]), "\n")

# ----------------------------
# TEST05 — maxN=0
# ----------------------------
cat("\nTEST05 — maxN=0\n")

TESTDIR <- file.path(TESTROOT, "T05_maxN_0")
dir.create(TESTDIR, recursive = TRUE, showWarnings = FALSE)

filtFs <- file.path(TESTDIR, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(TESTDIR, paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(
  fnFs, filtFs,
  fnRs, filtRs,
  maxN = 0,
  compress = TRUE,
  multithread = TRUE
)

saveRDS(out, file.path(outdir, "Obj", "out_T05_maxN_0.rds"))
write.csv(out, file.path(outdir, "Output", "post.filt_T05_maxN_0.csv"), quote = FALSE)
cat("T05 total in/out:", sum(out[, 1]), sum(out[, 2]), "\n")

# ----------------------------
# TEST06 — rm.phix + maxN + truncQ
# ----------------------------
cat("\nTEST06 — rm.phix + maxN + truncQ\n")

TESTDIR <- file.path(TESTROOT, "T06_rmPhix_maxN_truncQ")
dir.create(TESTDIR, recursive = TRUE, showWarnings = FALSE)

filtFs <- file.path(TESTDIR, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(TESTDIR, paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(
  fnFs, filtFs,
  fnRs, filtRs,
  maxN    = 0,
  truncQ  = 2,
  rm.phix = TRUE,
  compress = TRUE,
  multithread = TRUE
)

saveRDS(out, file.path(outdir, "Obj", "out_T06_rmPhix_maxN_truncQ.rds"))
write.csv(out, file.path(outdir, "Output", "post.filt_T06_rmPhix_maxN_truncQ.csv"), quote = FALSE)
cat("T06 total in/out:", sum(out[, 1]), sum(out[, 2]), "\n")

# ----------------------------
# TEST07 — maxEE=2,2
# ----------------------------
cat("\nTEST07 — maxEE=c(2,2)\n")

TESTDIR <- file.path(TESTROOT, "T07_maxEE_2_2")
dir.create(TESTDIR, recursive = TRUE, showWarnings = FALSE)

filtFs <- file.path(TESTDIR, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(TESTDIR, paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(
  fnFs, filtFs,
  fnRs, filtRs,
  maxEE = c(2, 2),
  compress = TRUE,
  multithread = TRUE
)

saveRDS(out, file.path(outdir, "Obj", "out_T07_maxEE_2_2.rds"))
write.csv(out, file.path(outdir, "Output", "post.filt_T07_maxEE_2_2.csv"), quote = FALSE)
cat("T07 total in/out:", sum(out[, 1]), sum(out[, 2]), "\n")

# ----------------------------
# TEST08 — truncLen=250,200
# ----------------------------
cat("\nTEST08 — truncLen=c(250,200)\n")

TESTDIR <- file.path(TESTROOT, "T08_truncLen_250_200")
dir.create(TESTDIR, recursive = TRUE, showWarnings = FALSE)

filtFs <- file.path(TESTDIR, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(TESTDIR, paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(
  fnFs, filtFs,
  fnRs, filtRs,
  truncLen = c(250, 200),
  compress = TRUE,
  multithread = TRUE
)

saveRDS(out, file.path(outdir, "Obj", "out_T08_truncLen_250_200.rds"))
write.csv(out, file.path(outdir, "Output", "post.filt_T08_truncLen_250_200.csv"), quote = FALSE)
cat("T08 total in/out:", sum(out[, 1]), sum(out[, 2]), "\n")

# ----------------------------
# TEST09 — (rm.phix + maxN + truncQ) + truncLen
# ----------------------------
cat("\nTEST09 — (rm.phix + maxN + truncQ) + truncLen\n")

TESTDIR <- file.path(TESTROOT, "T09_rmPhix_maxN_truncQ_truncLen")
dir.create(TESTDIR, recursive = TRUE, showWarnings = FALSE)

filtFs <- file.path(TESTDIR, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(TESTDIR, paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(
  fnFs, filtFs,
  fnRs, filtRs,
  truncLen = c(250, 200),
  maxN    = 0,
  truncQ  = 2,
  rm.phix = TRUE,
  compress = TRUE,
  multithread = TRUE
)

saveRDS(out, file.path(outdir, "Obj", "out_T09_rmPhix_maxN_truncQ_truncLen.rds"))
write.csv(out, file.path(outdir, "Output", "post.filt_T09_rmPhix_maxN_truncQ_truncLen.csv"), quote = FALSE)
cat("T09 total in/out:", sum(out[, 1]), sum(out[, 2]), "\n")

# ----------------------------
# TEST10 — (rm.phix + maxN + truncQ) + maxEE
# ----------------------------
cat("\nTEST10 — (rm.phix + maxN + truncQ) + maxEE\n")

TESTDIR <- file.path(TESTROOT, "T10_rmPhix_maxN_truncQ_maxEE")
dir.create(TESTDIR, recursive = TRUE, showWarnings = FALSE)

filtFs <- file.path(TESTDIR, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(TESTDIR, paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(
  fnFs, filtFs,
  fnRs, filtRs,
  maxEE   = c(2, 2),
  maxN    = 0,
  truncQ  = 2,
  rm.phix = TRUE,
  compress = TRUE,
  multithread = TRUE
)

saveRDS(out, file.path(outdir, "Obj", "out_T10_rmPhix_maxN_truncQ_maxEE.rds"))
write.csv(out, file.path(outdir, "Output", "post.filt_T10_rmPhix_maxN_truncQ_maxEE.csv"), quote = FALSE)
cat("T10 total in/out:", sum(out[, 1]), sum(out[, 2]), "\n")

# ----------------------------
# TEST11 — truncLen + maxEE
# ----------------------------
cat("\nTEST11 — truncLen + maxEE\n")

TESTDIR <- file.path(TESTROOT, "T11_truncLen_maxEE")
dir.create(TESTDIR, recursive = TRUE, showWarnings = FALSE)

filtFs <- file.path(TESTDIR, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(TESTDIR, paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(
  fnFs, filtFs,
  fnRs, filtRs,
  maxEE    = c(2, 2),
  truncLen = c(250, 200),
  compress = TRUE,
  multithread = TRUE
)

saveRDS(out, file.path(outdir, "Obj", "out_T11_truncLen_maxEE.rds"))
write.csv(out, file.path(outdir, "Output", "post.filt_T11_truncLen_maxEE.csv"), quote = FALSE)
cat("T11 total in/out:", sum(out[, 1]), sum(out[, 2]), "\n")

# ----------------------------
# TEST12 — full parameters (modif)
# ----------------------------
cat("\nTEST12 — full parameters maxEE33\n")

TESTDIR <- file.path(TESTROOT, "T12_full_MaxEE33")
dir.create(TESTDIR, recursive = TRUE, showWarnings = FALSE)

filtFs <- file.path(TESTDIR, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(TESTDIR, paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(
  fnFs, filtFs,
  fnRs, filtRs,
  truncLen = c(250, 200),
  maxN     = 0,
  maxEE    = c(3, 3),
  truncQ   = 2,
  rm.phix  = TRUE,
  compress = TRUE,
  multithread = TRUE
)

saveRDS(out, file.path(outdir, "Obj", "out_T12_full_MaxEE33.rds"))
write.csv(out, file.path(outdir, "Output", "post.filt_T12_full_MaxEE33.csv"), quote = FALSE)
cat("T12 total in/out:", sum(out[, 1]), sum(out[, 2]), "\n")

# ----------------------------
# TEST13 — full parameters (modif)
# ----------------------------
cat("\nTEST13 — full parameters maxEE44\n")

TESTDIR <- file.path(TESTROOT, "T13_full_MaxEE44")
dir.create(TESTDIR, recursive = TRUE, showWarnings = FALSE)

filtFs <- file.path(TESTDIR, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(TESTDIR, paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(
  fnFs, filtFs,
  fnRs, filtRs,
  truncLen = c(250, 200),
  maxN     = 0,
  maxEE    = c(4, 4),
  truncQ   = 2,
  rm.phix  = TRUE,
  compress = TRUE,
  multithread = TRUE
)

saveRDS(out, file.path(outdir, "Obj", "out_T13_full_MaxEE44.rds"))
write.csv(out, file.path(outdir, "Output", "post.filt_T13_full_MaxEE44.csv"), quote = FALSE)
cat("T13 total in/out:", sum(out[, 1]), sum(out[, 2]), "\n")

# ----------------------------
# TEST14 — full parameters (modif)
# ----------------------------
cat("\nTEST14 — full parameters maxEE34\n")

TESTDIR <- file.path(TESTROOT, "T14_full_MaxEE34")
dir.create(TESTDIR, recursive = TRUE, showWarnings = FALSE)

filtFs <- file.path(TESTDIR, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(TESTDIR, paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(
  fnFs, filtFs,
  fnRs, filtRs,
  truncLen = c(250, 200),
  maxN     = 0,
  maxEE    = c(3, 4),
  truncQ   = 2,
  rm.phix  = TRUE,
  compress = TRUE,
  multithread = TRUE
)

saveRDS(out, file.path(outdir, "Obj", "out_T14_full_MaxEE34.rds"))
write.csv(out, file.path(outdir, "Output", "post.filt_T14_full_MaxEE34.csv"), quote = FALSE)
cat("T14 total in/out:", sum(out[, 1]), sum(out[, 2]), "\n")

# ----------------------------
# TEST15 — full parameters (modif)
# ----------------------------
cat("\nTEST15 — full parameters maxEE24\n")

TESTDIR <- file.path(TESTROOT, "T15_full_MaxEE24")
dir.create(TESTDIR, recursive = TRUE, showWarnings = FALSE)

filtFs <- file.path(TESTDIR, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(TESTDIR, paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(
  fnFs, filtFs,
  fnRs, filtRs,
  truncLen = c(250, 200),
  maxN     = 0,
  maxEE    = c(2, 4),
  truncQ   = 2,
  rm.phix  = TRUE,
  compress = TRUE,
  multithread = TRUE
)

saveRDS(out, file.path(outdir, "Obj", "out_T15_full_MaxEE24.rds"))
write.csv(out, file.path(outdir, "Output", "post.filt_T15_full_MaxEE24.csv"), quote = FALSE)
cat("T15 total in/out:", sum(out[, 1]), sum(out[, 2]), "\n")

# ----------------------------
# TEST16 — full parameters (modif)
# ----------------------------
cat("\nTEST16 — full parameters maxEE23\n")

TESTDIR <- file.path(TESTROOT, "T16_full_MaxEE23")
dir.create(TESTDIR, recursive = TRUE, showWarnings = FALSE)

filtFs <- file.path(TESTDIR, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(TESTDIR, paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(
  fnFs, filtFs,
  fnRs, filtRs,
  truncLen = c(250, 200),
  maxN     = 0,
  maxEE    = c(2, 3),
  truncQ   = 2,
  rm.phix  = TRUE,
  compress = TRUE,
  multithread = TRUE
)

saveRDS(out, file.path(outdir, "Obj", "out_T16_full_MaxEE23.rds"))
write.csv(out, file.path(outdir, "Output", "post.filt_T16_full_MaxEE23.csv"), quote = FALSE)
cat("T16 total in/out:", sum(out[, 1]), sum(out[, 2]), "\n")

# ----------------------------
# TEST17 — full parameters (modif)
# ----------------------------
cat("\nTEST17 — full parameters maxEE24\n")

TESTDIR <- file.path(TESTROOT, "T17_full_MaxEE43")
dir.create(TESTDIR, recursive = TRUE, showWarnings = FALSE)

filtFs <- file.path(TESTDIR, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(TESTDIR, paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(
  fnFs, filtFs,
  fnRs, filtRs,
  truncLen = c(250, 200),
  maxN     = 0,
  maxEE    = c(4, 3),
  truncQ   = 2,
  rm.phix  = TRUE,
  compress = TRUE,
  multithread = TRUE
)

saveRDS(out, file.path(outdir, "Obj", "out_T17_full_MaxEE43.rds"))
write.csv(out, file.path(outdir, "Output", "post.filt_T17_full_MaxEE43.csv"), quote = FALSE)
cat("T17 total in/out:", sum(out[, 1]), sum(out[, 2]), "\n")

# ----------------------------
# TEST18 — full parameters (modif)
# ----------------------------
cat("\nTEST18 — full parameters maxEE23_trunc250170\n")

TESTDIR <- file.path(TESTROOT, "T18_full_maxEE23_trunc250170")
dir.create(TESTDIR, recursive = TRUE, showWarnings = FALSE)

filtFs <- file.path(TESTDIR, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(TESTDIR, paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(
  fnFs, filtFs,
  fnRs, filtRs,
  truncLen = c(250, 170),
  maxN     = 0,
  maxEE    = c(2, 3),
  truncQ   = 2,
  rm.phix  = TRUE,
  compress = TRUE,
  multithread = TRUE
)

saveRDS(out, file.path(outdir, "Obj", "out_T18_full_maxEE23_trunc250170.rds"))
write.csv(out, file.path(outdir, "Output", "post.filt_T18_full_maxEE23_trunc250170.csv"), quote = FALSE)
cat("T18 total in/out:", sum(out[, 1]), sum(out[, 2]), "\n")

# ----------------------------
# TEST19 — full parameters (modif)
# ----------------------------
cat("\nTEST19 — full parameters maxEE33_trunc250170\n")

TESTDIR <- file.path(TESTROOT, "T19_full_maxEE33_trunc250170")
dir.create(TESTDIR, recursive = TRUE, showWarnings = FALSE)

filtFs <- file.path(TESTDIR, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(TESTDIR, paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(
  fnFs, filtFs,
  fnRs, filtRs,
  truncLen = c(250, 170),
  maxN     = 0,
  maxEE    = c(3, 3),
  truncQ   = 2,
  rm.phix  = TRUE,
  compress = TRUE,
  multithread = TRUE
)

saveRDS(out, file.path(outdir, "Obj", "out_T19_full_maxEE33_trunc250170.rds"))
write.csv(out, file.path(outdir, "Output", "post.filt_T19_full_maxEE33_trunc250170.csv"), quote = FALSE)
cat("T19 total in/out:", sum(out[, 1]), sum(out[, 2]), "\n")

# ----------------------------
# TEST20 — full parameters (modif)
# ----------------------------
cat("\nTEST20 — full parameters maxEE34_trunc250170\n")

TESTDIR <- file.path(TESTROOT, "T20_full_maxEE34_trunc250170")
dir.create(TESTDIR, recursive = TRUE, showWarnings = FALSE)

filtFs <- file.path(TESTDIR, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(TESTDIR, paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(
  fnFs, filtFs,
  fnRs, filtRs,
  truncLen = c(250, 170),
  maxN     = 0,
  maxEE    = c(3, 4),
  truncQ   = 2,
  rm.phix  = TRUE,
  compress = TRUE,
  multithread = TRUE
)

saveRDS(out, file.path(outdir, "Obj", "out_T20_full_maxEE34_trunc250170.rds"))
write.csv(out, file.path(outdir, "Output", "post.filt_T20_full_maxEE34_trunc250170.csv"), quote = FALSE)
cat("T20 total in/out:", sum(out[, 1]), sum(out[, 2]), "\n")

# ----------------------------
# TEST21 — full parameters (modif)
# ----------------------------
cat("\nTEST21 — full parameters maxEE34_trunc250170\n")

TESTDIR <- file.path(TESTROOT, "T21_full_maxEE34_trunc250180")
dir.create(TESTDIR, recursive = TRUE, showWarnings = FALSE)

filtFs <- file.path(TESTDIR, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(TESTDIR, paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(
  fnFs, filtFs,
  fnRs, filtRs,
  truncLen = c(250, 180),
  maxN     = 0,
  maxEE    = c(3, 4),
  truncQ   = 2,
  rm.phix  = TRUE,
  compress = TRUE,
  multithread = TRUE
)

saveRDS(out, file.path(outdir, "Obj", "out_T21_full_maxEE34_trunc250180.rds"))
write.csv(out, file.path(outdir, "Output", "post.filt_T21_full_maxEE34_trunc250180.csv"), quote = FALSE)
cat("T21 total in/out:", sum(out[, 1]), sum(out[, 2]), "\n")

# ----------------------------
# TEST22 — full parameters (modif)
# ----------------------------
cat("\nTEST22 — Trunc250200\n")

TESTDIR <- file.path(TESTROOT, "T22_Trunc250200")
dir.create(TESTDIR, recursive = TRUE, showWarnings = FALSE)

filtFs <- file.path(TESTDIR, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(TESTDIR, paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(
  fnFs, filtFs,
  fnRs, filtRs,
  truncLen = c(250, 200),
  maxN     = 0,
  truncQ   = 2,
  rm.phix  = TRUE,
  compress = TRUE,
  multithread = TRUE
)

saveRDS(out, file.path(outdir, "Obj", "out_T22_Trunc250200.rds"))
write.csv(out, file.path(outdir, "Output", "post.filt_T22_Trunc250200.csv"), quote = FALSE)
cat("T22 total in/out:", sum(out[, 1]), sum(out[, 2]), "\n")

# ----------------------------
# TEST23 — full parameters (modif)
# ----------------------------
cat("\nTEST23 — Trunc250190\n")

TESTDIR <- file.path(TESTROOT, "T23_Trunc250190")
dir.create(TESTDIR, recursive = TRUE, showWarnings = FALSE)

filtFs <- file.path(TESTDIR, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(TESTDIR, paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(
  fnFs, filtFs,
  fnRs, filtRs,
  truncLen = c(250, 190),
  maxN     = 0,
  truncQ   = 2,
  rm.phix  = TRUE,
  compress = TRUE,
  multithread = TRUE
)

saveRDS(out, file.path(outdir, "Obj", "out_T23_Trunc250190.rds"))
write.csv(out, file.path(outdir, "Output", "post.filt_T23_Trunc250190.csv"), quote = FALSE)
cat("T23 total in/out:", sum(out[, 1]), sum(out[, 2]), "\n")

# ----------------------------
# TEST24 — full parameters (modif)
# ----------------------------
cat("\nTEST24 — Trunc250180\n")

TESTDIR <- file.path(TESTROOT, "T24_Trunc250180")
dir.create(TESTDIR, recursive = TRUE, showWarnings = FALSE)

filtFs <- file.path(TESTDIR, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(TESTDIR, paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(
  fnFs, filtFs,
  fnRs, filtRs,
  truncLen = c(250, 180),
  maxN     = 0,
  truncQ   = 2,
  rm.phix  = TRUE,
  compress = TRUE,
  multithread = TRUE
)

saveRDS(out, file.path(outdir, "Obj", "out_T24_Trunc250180.rds"))
write.csv(out, file.path(outdir, "Output", "post.filt_T24_Trunc250180.csv"), quote = FALSE)
cat("T24 total in/out:", sum(out[, 1]), sum(out[, 2]), "\n")

# ----------------------------
# TEST25 — full parameters (modif)
# ----------------------------
cat("\nTEST25 — Trunc250170\n")

TESTDIR <- file.path(TESTROOT, "T25_Trunc250170")
dir.create(TESTDIR, recursive = TRUE, showWarnings = FALSE)

filtFs <- file.path(TESTDIR, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(TESTDIR, paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(
  fnFs, filtFs,
  fnRs, filtRs,
  truncLen = c(250, 170),
  maxN     = 0,
  truncQ   = 2,
  rm.phix  = TRUE,
  compress = TRUE,
  multithread = TRUE
)

saveRDS(out, file.path(outdir, "Obj", "out_T25_Trunc250170.rds"))
write.csv(out, file.path(outdir, "Output", "post.filt_T25_Trunc250170.csv"), quote = FALSE)
cat("T25 total in/out:", sum(out[, 1]), sum(out[, 2]), "\n")

# ----------------------------
# Summary (cohérent avec les fichiers réellement créés)
# ----------------------------
summary_tests <- data.frame(
  test = c(
    "T01_full",
    "T02_default",
    "T03_rmPhix_TRUE",
    "T04_truncQ_2",
    "T05_maxN_0",
    "T06_rmPhix_maxN_truncQ",
    "T07_maxEE_2_2",
    "T08_truncLen_250_200",
    "T09_rmPhix_maxN_truncQ_truncLen",
    "T10_rmPhix_maxN_truncQ_maxEE",
    "T11_truncLen_maxEE",
    "T12_full_MaxEE33",
    "T13_full_MaxEE44",
    "T14_full_MaxEE34",
    "T15_full_MaxEE24",
    "T16_full_MaxEE23",
    "T17_full_MaxEE43",
    "T18_full_maxEE23_trunc250170",
    "T19_full_maxEE33_trunc250170",
    "T20_full_maxEE34_trunc250170",
    "T21_full_maxEE34_trunc250180",
    "T22_Trunc250200",
    "T23_Trunc250190",
    "T24_Trunc250180",
    "T25_Trunc250170"
  ),
  total_in = c(
    sum(readRDS(file.path(outdir, "Obj", "out_T01_full.rds"))[, 1]),
    sum(readRDS(file.path(outdir, "Obj", "out_T02_default.rds"))[, 1]),
    sum(readRDS(file.path(outdir, "Obj", "out_T03_rmPhix_TRUE.rds"))[, 1]),
    sum(readRDS(file.path(outdir, "Obj", "out_T04_truncQ_2.rds"))[, 1]),
    sum(readRDS(file.path(outdir, "Obj", "out_T05_maxN_0.rds"))[, 1]),
    sum(readRDS(file.path(outdir, "Obj", "out_T06_rmPhix_maxN_truncQ.rds"))[, 1]),
    sum(readRDS(file.path(outdir, "Obj", "out_T07_maxEE_2_2.rds"))[, 1]),
    sum(readRDS(file.path(outdir, "Obj", "out_T08_truncLen_250_200.rds"))[, 1]),
    sum(readRDS(file.path(outdir, "Obj", "out_T09_rmPhix_maxN_truncQ_truncLen.rds"))[, 1]),
    sum(readRDS(file.path(outdir, "Obj", "out_T10_rmPhix_maxN_truncQ_maxEE.rds"))[, 1]),
    sum(readRDS(file.path(outdir, "Obj", "out_T11_truncLen_maxEE.rds"))[, 1]),
    sum(readRDS(file.path(outdir, "Obj", "out_T12_full_MaxEE33.rds"))[, 1]),
    sum(readRDS(file.path(outdir, "Obj", "out_T13_full_MaxEE44.rds"))[, 1]),
    sum(readRDS(file.path(outdir, "Obj", "out_T14_full_MaxEE34.rds"))[, 1]),
    sum(readRDS(file.path(outdir, "Obj", "out_T15_full_MaxEE24.rds"))[, 1]),
    sum(readRDS(file.path(outdir, "Obj", "out_T16_full_MaxEE23.rds"))[, 1]),
    sum(readRDS(file.path(outdir, "Obj", "out_T17_full_MaxEE43.rds"))[, 1]),
    sum(readRDS(file.path(outdir, "Obj", "out_T18_full_maxEE23_trunc250170.rds"))[, 1]),
    sum(readRDS(file.path(outdir, "Obj", "out_T19_full_maxEE33_trunc250170.rds"))[, 1]),
    sum(readRDS(file.path(outdir, "Obj", "out_T20_full_maxEE34_trunc250170.rds"))[, 1]),
    sum(readRDS(file.path(outdir, "Obj", "out_T21_full_maxEE34_trunc250180.rds"))[, 1]),
    sum(readRDS(file.path(outdir, "Obj", "out_T22_Trunc250200.rds"))[, 1]),
    sum(readRDS(file.path(outdir, "Obj", "out_T23_Trunc250190.rds"))[, 1]),
    sum(readRDS(file.path(outdir, "Obj", "out_T24_Trunc250180.rds"))[, 1]),
    sum(readRDS(file.path(outdir, "Obj", "out_T25_Trunc250170.rds"))[, 1])
  ),
  total_out = c(
    sum(readRDS(file.path(outdir, "Obj", "out_T01_full.rds"))[, 2]),
    sum(readRDS(file.path(outdir, "Obj", "out_T02_default.rds"))[, 2]),
    sum(readRDS(file.path(outdir, "Obj", "out_T03_rmPhix_TRUE.rds"))[, 2]),
    sum(readRDS(file.path(outdir, "Obj", "out_T04_truncQ_2.rds"))[, 2]),
    sum(readRDS(file.path(outdir, "Obj", "out_T05_maxN_0.rds"))[, 2]),
    sum(readRDS(file.path(outdir, "Obj", "out_T06_rmPhix_maxN_truncQ.rds"))[, 2]),
    sum(readRDS(file.path(outdir, "Obj", "out_T07_maxEE_2_2.rds"))[, 2]),
    sum(readRDS(file.path(outdir, "Obj", "out_T08_truncLen_250_200.rds"))[, 2]),
    sum(readRDS(file.path(outdir, "Obj", "out_T09_rmPhix_maxN_truncQ_truncLen.rds"))[, 2]),
    sum(readRDS(file.path(outdir, "Obj", "out_T10_rmPhix_maxN_truncQ_maxEE.rds"))[, 2]),
    sum(readRDS(file.path(outdir, "Obj", "out_T11_truncLen_maxEE.rds"))[, 2]),
    sum(readRDS(file.path(outdir, "Obj", "out_T12_full_MaxEE33.rds"))[, 2]),
    sum(readRDS(file.path(outdir, "Obj", "out_T13_full_MaxEE44.rds"))[, 2]),
    sum(readRDS(file.path(outdir, "Obj", "out_T14_full_MaxEE34.rds"))[, 2]),
    sum(readRDS(file.path(outdir, "Obj", "out_T15_full_MaxEE24.rds"))[, 2]),
    sum(readRDS(file.path(outdir, "Obj", "out_T16_full_MaxEE23.rds"))[, 2]),
    sum(readRDS(file.path(outdir, "Obj", "out_T17_full_MaxEE43.rds"))[, 2]),
    sum(readRDS(file.path(outdir, "Obj", "out_T18_full_maxEE23_trunc250170.rds"))[, 2]),
    sum(readRDS(file.path(outdir, "Obj", "out_T19_full_maxEE33_trunc250170.rds"))[, 2]),
    sum(readRDS(file.path(outdir, "Obj", "out_T20_full_maxEE34_trunc250170.rds"))[, 2]),
    sum(readRDS(file.path(outdir, "Obj", "out_T21_full_maxEE34_trunc250180.rds"))[, 2]),
    sum(readRDS(file.path(outdir, "Obj", "out_T22_Trunc250200.rds"))[, 2]),
    sum(readRDS(file.path(outdir, "Obj", "out_T23_Trunc250190.rds"))[, 2]),
    sum(readRDS(file.path(outdir, "Obj", "out_T24_Trunc250180.rds"))[, 2]),
    sum(readRDS(file.path(outdir, "Obj", "out_T25_Trunc250170.rds"))[, 2])
  )
)

summary_tests$kept_ratio <- with(summary_tests, ifelse(total_in > 0, total_out / total_in, NA_real_))

write.csv(summary_tests,
          file.path(outdir, "Output", "filter_tests_summary.csv"),
          row.names = FALSE,
          quote = FALSE)

print(summary_tests)



suppressPackageStartupMessages({
  library(ShortRead)
})

# ------------------------------------------------------------
# Helper: compute pooled EE from a vector of fastq paths
# (memory-safe-ish: accumulates in numeric, ok for moderate size;
# for very large datasets, see note below)
# ------------------------------------------------------------
get_EE_pooled <- function(fq_paths, verbose = TRUE) {
  EE_all <- numeric()

  for (p in fq_paths) {
    if (!file.exists(p)) {
      if (verbose) message("Missing: ", basename(p))
      next
    }
    if (verbose) message("Reading: ", basename(p))

    fq <- readFastq(p)
    q  <- as(quality(fq), "matrix")
    EE <- rowSums(10^(-q/10))

    EE_all <- c(EE_all, EE)

    rm(fq, q, EE)
    gc()
  }

  EE_all
}

# ------------------------------------------------------------
# Tests to plot (after filterAndTrim outputs)
# ------------------------------------------------------------
tests <- list(
  T22 = list(dir = file.path(TESTROOT, "T22_Trunc250200"),
             label = "T22 truncLen=250/200"),
  T23 = list(dir = file.path(TESTROOT, "T23_Trunc250190"),
             label = "T23 truncLen=250/190"),
  T24 = list(dir = file.path(TESTROOT, "T24_Trunc250180"),
             label = "T24 truncLen=250/180"),
  T25 = list(dir = file.path(TESTROOT, "T25_Trunc250170"),
             label = "T25 truncLen=250/170")
)

# Make sure output dir exists
dir.create(file.path(outdir, "Figure"), recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 1) FORWARD pooled EE — one panel per test
# ------------------------------------------------------------
pdf(file.path(outdir, "Figure", "EE_Distribution_Forward_ALLsamples_PostFilt_T22-25.pdf"),
    width = 12, height = 9)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

for (tn in names(tests)) {
  TESTDIR <- tests[[tn]]$dir
  lab     <- tests[[tn]]$label

  filtFs <- file.path(TESTDIR, paste0(sample.names, "_F_filt.fastq.gz"))

  EE_all <- get_EE_pooled(filtFs, verbose = FALSE)

  if (length(EE_all) == 0) {
    plot.new()
    title(main = paste0(lab, "\n(no reads / files missing)"))
    next
  }

  h <- hist(EE_all, breaks = 200, plot = FALSE)

  plot(h,
       main = paste0(lab, "\nForward pooled (N=", length(EE_all), ")"),
       xlab = "Expected Errors",
       ylab = "Reads",
       xlim = c(0, 10),
       col = "grey80",
       border = "grey30")

  abline(v = c(1, 2, 3), lty = 2, col = "red")
  mtext(paste0("max bin = ", max(h$counts)), side = 3, line = -1, cex = 0.8)
}

dev.off()
cat("Saved: EE_Distribution_Forward_ALLsamples_PostFilt_T22-25.pdf\n")

# ------------------------------------------------------------
# 2) REVERSE pooled EE — one panel per test
# ------------------------------------------------------------
pdf(file.path(outdir, "Figure", "EE_Distribution_Reverse_ALLsamples_PostFilt_T22-25.pdf"),
    width = 12, height = 9)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

for (tn in names(tests)) {
  TESTDIR <- tests[[tn]]$dir
  lab     <- tests[[tn]]$label

  filtRs <- file.path(TESTDIR, paste0(sample.names, "_R_filt.fastq.gz"))

  EE_all <- get_EE_pooled(filtRs, verbose = FALSE)

  if (length(EE_all) == 0) {
    plot.new()
    title(main = paste0(lab, "\n(no reads / files missing)"))
    next
  }

  h <- hist(EE_all, breaks = 200, plot = FALSE)

  plot(h,
       main = paste0(lab, "\nReverse pooled (N=", length(EE_all), ")"),
       xlab = "Expected Errors",
       ylab = "Reads",
       xlim = c(0, 10),
       col = "grey80",
       border = "grey30")

  abline(v = c(1, 2, 3), lty = 2, col = "red")
  mtext(paste0("max bin = ", max(h$counts)), side = 3, line = -1, cex = 0.8)
}

dev.off()
cat("Saved: EE_Distribution_Reverse_ALLsamples_PostFilt_T22-25.pdf\n")

# ============================================================
# MERGING tests (subset): T14, T20, T21  — PATH FIX
# ============================================================

# paramètres mergePairs (modifiables)
minOv  <- 12
maxMis <- 1

MERGE_TAG <- paste0("minOv", minOv, "_maxMis", maxMis)

# base dir REEL des fastq filtrés
MERGE_TESTROOT <- file.path(outdir, "Filt_param_tests")

# ----------------------------
# TEST14 — MERGING
# ----------------------------
cat("\nTEST14 — MERGING T14_full_MaxEE34\n")

TESTID  <- "T14_full_MaxEE34"
TESTDIR <- file.path(MERGE_TESTROOT, TESTID)

out <- readRDS(file.path(outdir, "Obj", paste0("out_", TESTID, ".rds")))
sample_test <- rownames(out)

filtFs <- file.path(TESTDIR, paste0(sample_test, "_F_filt.fastq.gz"))
filtRs <- file.path(TESTDIR, paste0(sample_test, "_R_filt.fastq.gz"))
names(filtFs) <- sample_test
names(filtRs) <- sample_test

keep_idx <- out[, 2] > 0
out2   <- out[keep_idx, , drop = FALSE]
filtFs <- filtFs[keep_idx]
filtRs <- filtRs[keep_idx]

keep_files <- file.exists(filtFs) & file.exists(filtRs)
out2   <- out2[keep_files, , drop = FALSE]
filtFs <- filtFs[keep_files]
filtRs <- filtRs[keep_files]

cat("Samples kept:", nrow(out2), "\n")
if (nrow(out2) == 0) stop("No sample kept for merging in ", TESTID,
                          " (check MERGE_TESTROOT/TESTID contains *_filt.fastq.gz).")

cat("4. Learn errors\n")
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

saveRDS(errF, file.path(outdir, "Obj", paste0("errF_", TESTID, ".rds")))
saveRDS(errR, file.path(outdir, "Obj", paste0("errR_", TESTID, ".rds")))

pdf(file.path(outdir, "Figure", paste0("DADA2_ErrorModels_", TESTID, ".pdf")), 8, 6)
plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)
dev.off()

cat("5. Denoise\n")
dadaFs <- dada(filtFs, err = errF, multithread = TRUE)
dadaRs <- dada(filtRs, err = errR, multithread = TRUE)

denoised_F <- sum(sapply(dadaFs, function(d) sum(getUniques(d), na.rm = TRUE)), na.rm = TRUE)
denoised_R <- sum(sapply(dadaRs, function(d) sum(getUniques(d), na.rm = TRUE)), na.rm = TRUE)

cat("6. Merge + seqtab\n")
mergers <- mergePairs(
  dadaFs, filtFs,
  dadaRs, filtRs,
  minOverlap  = minOv,
  maxMismatch = maxMis,
  verbose     = TRUE
)

merged_reads <- sum(sapply(mergers, function(m) {
  if (is.null(m) || nrow(m) == 0) return(0)
  sum(m[, "abundance"], na.rm = TRUE)
}), na.rm = TRUE)

seqtab <- makeSequenceTable(mergers)
seqtab_reads <- sum(rowSums(seqtab), na.rm = TRUE)
nASV <- ncol(seqtab)

saveRDS(mergers, file.path(outdir, "Obj",
        paste0("mergers_", TESTID, "_", MERGE_TAG, ".rds")))
saveRDS(seqtab, file.path(outdir, "Obj",
        paste0("seqtab_", TESTID, "_", MERGE_TAG, ".rds")))
write.csv(seqtab, file.path(outdir, "Output",
        paste0("SeqTable_", TESTID, "_", MERGE_TAG, ".csv")), quote = FALSE)


merging_summary_T14 <- data.frame(
  test = TESTID,
  n_samples_kept = nrow(out2),
  total_in = sum(out[, 1], na.rm = TRUE),
  total_filtered = sum(out[, 2], na.rm = TRUE),
  denoised_F = denoised_F,
  denoised_R = denoised_R,
  merged_reads = merged_reads,
  seqtab_reads = seqtab_reads,
  nASV = nASV,
  minOverlap = minOv,
  maxMismatch = maxMis,
  stringsAsFactors = FALSE
)

saveRDS(merging_summary_T14, file.path(outdir, "Obj", paste0("merging_summary_", TESTID, "_", MERGE_TAG, ".rds")))
write.csv(merging_summary_T14,
          file.path(outdir, "Output", paste0("merging_summary_", TESTID, "_", MERGE_TAG, ".csv")),
          row.names = FALSE, quote = FALSE)

print(merging_summary_T14)


# ----------------------------
# TEST20 — MERGING
# ----------------------------
cat("\nTEST20 — MERGING T20_full_maxEE34_trunc250170\n")

TESTID  <- "T20_full_maxEE34_trunc250170"
TESTDIR <- file.path(MERGE_TESTROOT, TESTID)

out <- readRDS(file.path(outdir, "Obj", paste0("out_", TESTID, ".rds")))
sample_test <- rownames(out)

filtFs <- file.path(TESTDIR, paste0(sample_test, "_F_filt.fastq.gz"))
filtRs <- file.path(TESTDIR, paste0(sample_test, "_R_filt.fastq.gz"))
names(filtFs) <- sample_test
names(filtRs) <- sample_test

keep_idx <- out[, 2] > 0
out2   <- out[keep_idx, , drop = FALSE]
filtFs <- filtFs[keep_idx]
filtRs <- filtRs[keep_idx]

keep_files <- file.exists(filtFs) & file.exists(filtRs)
out2   <- out2[keep_files, , drop = FALSE]
filtFs <- filtFs[keep_files]
filtRs <- filtRs[keep_files]

cat("Samples kept:", nrow(out2), "\n")
if (nrow(out2) == 0) stop("No sample kept for merging in ", TESTID)

cat("4. Learn errors\n")
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

saveRDS(errF, file.path(outdir, "Obj", paste0("errF_", TESTID, ".rds")))
saveRDS(errR, file.path(outdir, "Obj", paste0("errR_", TESTID, ".rds")))

pdf(file.path(outdir, "Figure", paste0("DADA2_ErrorModels_", TESTID, ".pdf")), 8, 6)
plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)
dev.off()

cat("5. Denoise\n")
dadaFs <- dada(filtFs, err = errF, multithread = TRUE)
dadaRs <- dada(filtRs, err = errR, multithread = TRUE)

denoised_F <- sum(sapply(dadaFs, function(d) sum(getUniques(d), na.rm = TRUE)), na.rm = TRUE)
denoised_R <- sum(sapply(dadaRs, function(d) sum(getUniques(d), na.rm = TRUE)), na.rm = TRUE)

cat("6. Merge + seqtab\n")
mergers <- mergePairs(
  dadaFs, filtFs,
  dadaRs, filtRs,
  minOverlap  = minOv,
  maxMismatch = maxMis,
  verbose     = TRUE
)

merged_reads <- sum(sapply(mergers, function(m) {
  if (is.null(m) || nrow(m) == 0) return(0)
  sum(m[, "abundance"], na.rm = TRUE)
}), na.rm = TRUE)

seqtab <- makeSequenceTable(mergers)
seqtab_reads <- sum(rowSums(seqtab), na.rm = TRUE)
nASV <- ncol(seqtab)

saveRDS(mergers, file.path(outdir, "Obj",
        paste0("mergers_", TESTID, "_", MERGE_TAG, ".rds")))
saveRDS(seqtab, file.path(outdir, "Obj",
        paste0("seqtab_", TESTID, "_", MERGE_TAG, ".rds")))
write.csv(seqtab, file.path(outdir, "Output",
        paste0("SeqTable_", TESTID, "_", MERGE_TAG, ".csv")), quote = FALSE)

merging_summary_T20 <- data.frame(
  test = TESTID,
  n_samples_kept = nrow(out2),
  total_in = sum(out[, 1], na.rm = TRUE),
  total_filtered = sum(out[, 2], na.rm = TRUE),
  denoised_F = denoised_F,
  denoised_R = denoised_R,
  merged_reads = merged_reads,
  seqtab_reads = seqtab_reads,
  nASV = nASV,
  minOverlap = minOv,
  maxMismatch = maxMis,
  stringsAsFactors = FALSE
)

saveRDS(merging_summary_T20, file.path(outdir, "Obj", paste0("merging_summary_", TESTID, "_", MERGE_TAG, ".rds")))
write.csv(merging_summary_T20,
          file.path(outdir, "Output", paste0("merging_summary_", TESTID, "_", MERGE_TAG, ".csv")),
          row.names = FALSE, quote = FALSE)

print(merging_summary_T20)


# ----------------------------
# TEST21 — MERGING
# ----------------------------
cat("\nTEST21 — MERGING T21_full_maxEE34_trunc250180\n")

TESTID  <- "T21_full_maxEE34_trunc250180"
TESTDIR <- file.path(MERGE_TESTROOT, TESTID)

out <- readRDS(file.path(outdir, "Obj", paste0("out_", TESTID, ".rds")))
sample_test <- rownames(out)

filtFs <- file.path(TESTDIR, paste0(sample_test, "_F_filt.fastq.gz"))
filtRs <- file.path(TESTDIR, paste0(sample_test, "_R_filt.fastq.gz"))
names(filtFs) <- sample_test
names(filtRs) <- sample_test

keep_idx <- out[, 2] > 0
out2   <- out[keep_idx, , drop = FALSE]
filtFs <- filtFs[keep_idx]
filtRs <- filtRs[keep_idx]

keep_files <- file.exists(filtFs) & file.exists(filtRs)
out2   <- out2[keep_files, , drop = FALSE]
filtFs <- filtFs[keep_files]
filtRs <- filtRs[keep_files]

cat("Samples kept:", nrow(out2), "\n")
if (nrow(out2) == 0) stop("No sample kept for merging in ", TESTID)

cat("4. Learn errors\n")
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

saveRDS(errF, file.path(outdir, "Obj", paste0("errF_", TESTID, ".rds")))
saveRDS(errR, file.path(outdir, "Obj", paste0("errR_", TESTID, ".rds")))

pdf(file.path(outdir, "Figure", paste0("DADA2_ErrorModels_", TESTID, ".pdf")), 8, 6)
plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)
dev.off()

cat("5. Denoise\n")
dadaFs <- dada(filtFs, err = errF, multithread = TRUE)
dadaRs <- dada(filtRs, err = errR, multithread = TRUE)

denoised_F <- sum(sapply(dadaFs, function(d) sum(getUniques(d), na.rm = TRUE)), na.rm = TRUE)
denoised_R <- sum(sapply(dadaRs, function(d) sum(getUniques(d), na.rm = TRUE)), na.rm = TRUE)

cat("6. Merge + seqtab\n")
mergers <- mergePairs(
  dadaFs, filtFs,
  dadaRs, filtRs,
  minOverlap  = minOv,
  maxMismatch = maxMis,
  verbose     = TRUE
)

merged_reads <- sum(sapply(mergers, function(m) {
  if (is.null(m) || nrow(m) == 0) return(0)
  sum(m[, "abundance"], na.rm = TRUE)
}), na.rm = TRUE)

seqtab <- makeSequenceTable(mergers)
seqtab_reads <- sum(rowSums(seqtab), na.rm = TRUE)
nASV <- ncol(seqtab)

saveRDS(mergers, file.path(outdir, "Obj",
        paste0("mergers_", TESTID, "_", MERGE_TAG, ".rds")))
saveRDS(seqtab, file.path(outdir, "Obj",
        paste0("seqtab_", TESTID, "_", MERGE_TAG, ".rds")))
write.csv(seqtab, file.path(outdir, "Output",
        paste0("SeqTable_", TESTID, "_", MERGE_TAG, ".csv")), quote = FALSE)

merging_summary_T21 <- data.frame(
  test = TESTID,
  n_samples_kept = nrow(out2),
  total_in = sum(out[, 1], na.rm = TRUE),
  total_filtered = sum(out[, 2], na.rm = TRUE),
  denoised_F = denoised_F,
  denoised_R = denoised_R,
  merged_reads = merged_reads,
  seqtab_reads = seqtab_reads,
  nASV = nASV,
  minOverlap = minOv,
  maxMismatch = maxMis,
  stringsAsFactors = FALSE
)

saveRDS(merging_summary_T21, file.path(outdir, "Obj", paste0("merging_summary_", TESTID, "_", MERGE_TAG, ".rds")))
write.csv(merging_summary_T21,
          file.path(outdir, "Output", paste0("merging_summary_", TESTID, "_", MERGE_TAG, ".csv")),
          row.names = FALSE, quote = FALSE)


print(merging_summary_T21)


# ============================================================
# Summary MERGING (T14, T20, T21)
# ============================================================
summary_merging <- rbind(merging_summary_T14, merging_summary_T20, merging_summary_T21)

summary_merging$kept_ratio_filter <- with(summary_merging,
                                          ifelse(total_in > 0, total_filtered / total_in, NA_real_))
summary_merging$kept_ratio_merge <- with(summary_merging,
                                         ifelse(total_filtered > 0, merged_reads / total_filtered, NA_real_))

write.csv(summary_merging,
          file.path(outdir, "Output", "merging_tests_summary_T14_T20_T21.csv"),
          row.names = FALSE, quote = FALSE)

print(summary_merging)

# ============================================================
# MERGING tests (subset): T14, T20, T21  — ROBUST VERSION
# ============================================================

suppressPackageStartupMessages({
  library(dada2)
})

# paramètres mergePairs (modifiables)
minOv  <- 12
maxMis <- 1
MERGE_TAG <- paste0("minOv", minOv, "_maxMis", maxMis)

# base dir REEL des fastq filtrés
MERGE_TESTROOT <- file.path(outdir, "Filt_param_tests")

# ----------------------------
# Helper: list paired filt fastqs in a TESTDIR
# ----------------------------
get_paired_filt_fastqs <- function(testdir) {
  Fs <- list.files(testdir, pattern = "_F_filt\\.fastq\\.gz$", full.names = TRUE)
  Rs <- list.files(testdir, pattern = "_R_filt\\.fastq\\.gz$", full.names = TRUE)

  if (length(Fs) == 0) Fs <- list.files(testdir, pattern = "_F_filt\\.fastq$", full.names = TRUE)
  if (length(Rs) == 0) Rs <- list.files(testdir, pattern = "_R_filt\\.fastq$", full.names = TRUE)

  Fs <- sort(Fs); Rs <- sort(Rs)

  keyF <- sub("_F_filt\\.fastq(\\.gz)?$", "", basename(Fs))
  keyR <- sub("_R_filt\\.fastq(\\.gz)?$", "", basename(Rs))
  common <- intersect(keyF, keyR)

  if (length(common) == 0) {
    return(list(filtFs = character(0), filtRs = character(0)))
  }

  Fs2 <- Fs[match(common, keyF)]
  Rs2 <- Rs[match(common, keyR)]

  names(Fs2) <- common
  names(Rs2) <- common

  list(filtFs = Fs2, filtRs = Rs2)
}

# ----------------------------
# Helper: run one test (T14/T20/T21)
# ----------------------------
run_merge_test <- function(TESTID) {

  cat("\n====================================\n")
  cat("MERGING —", TESTID, "\n")

  TESTDIR <- file.path(MERGE_TESTROOT, TESTID)
  if (!dir.exists(TESTDIR)) stop("Missing TESTDIR: ", TESTDIR)

  # list actual filt fastqs present in folder
  pair <- get_paired_filt_fastqs(TESTDIR)
  filtFs <- pair$filtFs
  filtRs <- pair$filtRs

  if (length(filtFs) == 0) {
    stop("No *_F_filt.fastq(.gz)/*_R_filt.fastq(.gz) pairs found in: ", TESTDIR)
  }

  # optional: use out (reads.out) to drop samples with 0 reads after filtering
  out_path <- file.path(outdir, "Obj", paste0("out_", TESTID, ".rds"))
  if (file.exists(out_path)) {
    out <- readRDS(out_path)

    # out rows can be basenames of input fastq (often *_R1*.fastq.gz)
    # -> try to build a key comparable to our filt keys:
    #    remove "_R1..." suffix if present, and/or any extension
    out_key <- rownames(out)
    out_key2 <- sub("_R1.*$", "", out_key)       # if naming contains _R1...
    out_key2 <- sub("\\.fastq\\.gz$", "", out_key2)
    out_key2 <- sub("\\.fastq$", "", out_key2)

    # our filt keys are names(filtFs) which look like:
    # RH24-..._L001   (because file is ..._L001_F_filt.fastq.gz)
    filt_key <- names(filtFs)

    m <- match(filt_key, out_key2)

    # if many NA, we don't filter on out (safer than dropping everything)
    na_ratio <- mean(is.na(m))
    if (!is.nan(na_ratio) && na_ratio < 0.5) {
      keep_idx <- !is.na(m) & (out[m, "reads.out"] > 0)
      filtFs <- filtFs[keep_idx]
      filtRs <- filtRs[keep_idx]
    } else {
      cat("Warning: could not reliably map out rows to filt fastqs (NA ratio =", na_ratio, "). Skipping out-based filtering.\n")
    }
  } else {
    cat("Note: out_*.rds not found for", TESTID, "-> skipping out-based filtering.\n")
  }

  # keep only existing files
  keep_files <- file.exists(filtFs) & file.exists(filtRs)
  filtFs <- filtFs[keep_files]
  filtRs <- filtRs[keep_files]

  cat("Samples kept:", length(filtFs), "\n")
  if (length(filtFs) == 0) stop("No sample kept for merging in ", TESTID)

  # Learn errors
  cat("4. Learn errors\n")
  errF <- learnErrors(filtFs, multithread = TRUE)
  errR <- learnErrors(filtRs, multithread = TRUE)

  saveRDS(errF, file.path(outdir, "Obj", paste0("errF_", TESTID, ".rds")))
  saveRDS(errR, file.path(outdir, "Obj", paste0("errR_", TESTID, ".rds")))

  pdf(file.path(outdir, "Figure", paste0("DADA2_ErrorModels_", TESTID, ".pdf")), 8, 6)
  plotErrors(errF, nominalQ = TRUE)
  plotErrors(errR, nominalQ = TRUE)
  dev.off()

  # Denoise
  cat("5. Denoise\n")
  dadaFs <- dada(filtFs, err = errF, multithread = TRUE)
  dadaRs <- dada(filtRs, err = errR, multithread = TRUE)

  denoised_F <- sum(sapply(dadaFs, function(d) sum(getUniques(d), na.rm = TRUE)), na.rm = TRUE)
  denoised_R <- sum(sapply(dadaRs, function(d) sum(getUniques(d), na.rm = TRUE)), na.rm = TRUE)

  # Merge + seqtab
  cat("6. Merge + seqtab\n")
  mergers <- mergePairs(
    dadaFs, filtFs,
    dadaRs, filtRs,
    minOverlap  = minOv,
    maxMismatch = maxMis,
    verbose     = TRUE
  )

  merged_reads <- sum(sapply(mergers, function(m) {
    if (is.null(m) || nrow(m) == 0) return(0)
    sum(m[, "abundance"], na.rm = TRUE)
  }), na.rm = TRUE)

  seqtab <- makeSequenceTable(mergers)
  seqtab_reads <- sum(rowSums(seqtab), na.rm = TRUE)
  nASV <- ncol(seqtab)

  saveRDS(mergers, file.path(outdir, "Obj", paste0("mergers_", TESTID, "_", MERGE_TAG, ".rds")))
  saveRDS(seqtab,  file.path(outdir, "Obj", paste0("seqtab_",  TESTID, "_", MERGE_TAG, ".rds")))

  write.csv(seqtab,
            file.path(outdir, "Output", paste0("SeqTable_", TESTID, "_", MERGE_TAG, ".csv")),
            quote = FALSE)

  # summary
  merging_summary <- data.frame(
    test = TESTID,
    n_samples_kept = length(filtFs),
    denoised_F = denoised_F,
    denoised_R = denoised_R,
    merged_reads = merged_reads,
    seqtab_reads = seqtab_reads,
    nASV = nASV,
    minOverlap = minOv,
    maxMismatch = maxMis,
    stringsAsFactors = FALSE
  )

  saveRDS(merging_summary,
          file.path(outdir, "Obj", paste0("merging_summary_", TESTID, "_", MERGE_TAG, ".rds")))
  write.csv(merging_summary,
            file.path(outdir, "Output", paste0("merging_summary_", TESTID, "_", MERGE_TAG, ".csv")),
            row.names = FALSE, quote = FALSE)

  print(merging_summary)
  merging_summary
}

# ----------------------------
# Run tests
# ----------------------------
merging_summary_T14 <- run_merge_test("T14_full_MaxEE34")
merging_summary_T20 <- run_merge_test("T20_full_maxEE34_trunc250170")
merging_summary_T21 <- run_merge_test("T21_full_maxEE34_trunc250180")

# ----------------------------
# Summary MERGING (T14, T20, T21)
# ----------------------------
summary_merging <- rbind(merging_summary_T14, merging_summary_T20, merging_summary_T21)

# ratios computed only if you also want them:
# (ici on n'a pas forcément total_in / total_filtered si mapping out incertain)
summary_merging$kept_ratio_merge <- with(summary_merging,
                                        ifelse(is.finite(seqtab_reads) & seqtab_reads > 0,
                                               merged_reads / seqtab_reads, NA_real_))

write.csv(summary_merging,
          file.path(outdir, "Output", "merging_tests_summary_T14_T20_T21.csv"),
          row.names = FALSE, quote = FALSE)

print(summary_merging)

}