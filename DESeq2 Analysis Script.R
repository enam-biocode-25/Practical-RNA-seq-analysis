rna-seq-gene-expression-workflow/
  ├── data/
  │   ├── fastq/
  │   ├── reference/
  │   │   ├── transcripts.fa
│   │   └── transcript_to_gene.tsv
│   └── metadata/
  │       └── sample_metadata.tsv
├── scripts/
  │   ├── 01_fastqc.sh
│   ├── 02_trim_galore.sh
│   ├── 03_salmon_index.sh
│   ├── 04_salmon_quant.sh
│   ├── 05_deseq2_analysis.R
│   ├── 06_pca_plot.R
│   ├── 07_volcano_plot.R
│   └── 08_heatmap.R
├── results/
  │   ├── fastqc/
  │   ├── trimmed/
  │   ├── salmon_index/
  │   ├── salmon_quant/
  │   ├── deseq2/
  │   └── figures/
  └── run_pipeline.sh



## DESEQ2 ANALYSIS
library(tximport)
library(DESeq2)
library(readr)
library(dplyr)
library(tibble)

dir.create("results/deseq2", recursive = TRUE, showWarnings = FALSE)

samples <- read.delim("data/metadata/sample_metadata.tsv", stringsAsFactors = FALSE)

files <- file.path("results/salmon_quant", samples$sample_id, "quant.sf")
names(files) <- samples$sample_id



# Check that all quant files exist
missing_files <- files[!file.exists(files)]
if (length(missing_files) > 0) {
  stop("Missing quant.sf files:\n", paste(missing_files, collapse = "\n"))
}

tx2gene <- read.delim(
  "data/reference/transcript_to_gene.tsv",
  header = FALSE,
  stringsAsFactors = FALSE
)
colnames(tx2gene) <- c("TXNAME", "GENEID")

txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

samples$condition <- factor(samples$condition, levels = c("control", "treated"))
rownames(samples) <- samples$sample_id

dds <- DESeqDataSetFromTximport(
  txi = txi,
  colData = samples,
  design = ~ condition
)



# Optional light prefiltering
dds <- dds[rowSums(counts(dds)) >= 10, ]

dds <- DESeq(dds)

res <- results(dds, contrast = c("condition", "treated", "control"))
res_df <- as.data.frame(res) %>%
  rownames_to_column("gene") %>%
  arrange(padj)

write.csv(res_df, "results/deseq2/deseq2_results_all.csv", row.names = FALSE)

sig_df <- res_df %>%
  filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) >= 1)

write.csv(sig_df, "results/deseq2/deseq2_results_significant.csv", row.names = FALSE)

# Save normalized counts
norm_counts <- counts(dds, normalized = TRUE) %>%
  as.data.frame() %>%
  rownames_to_column("gene")

write.csv(norm_counts, "results/deseq2/normalized_counts.csv", row.names = FALSE)



# Save transformed object for plotting scripts
vsd <- vst(dds, blind = FALSE)
saveRDS(dds, "results/deseq2/dds.rds")
saveRDS(vsd, "results/deseq2/vsd.rds")

### PCA PLOT

library(DESeq2)
library(ggplot2)

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

vsd <- readRDS("results/deseq2/vsd.rds")

p <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(p, "percentVar"))

ggplot(p, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw()

ggsave("results/figures/pca_plot.png", width = 7, height = 5, dpi = 300)




### VOLCANO PLOT

library(readr)
library(dplyr)
library(ggplot2)

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

res <- read_csv("results/deseq2/deseq2_results_all.csv", show_col_types = FALSE)

res <- res %>%
  mutate(
    status = case_when(
      !is.na(padj) & padj < 0.05 & log2FoldChange >= 1 ~ "Up",
      !is.na(padj) & padj < 0.05 & log2FoldChange <= -1 ~ "Down",
      TRUE ~ "NS"
    ),
    neglog10padj = -log10(padj)
  )

ggplot(res, aes(x = log2FoldChange, y = neglog10padj, shape = status)) +
  geom_point(alpha = 0.7, size = 2) +
  theme_bw() +
  labs(
    title = "Volcano plot",
    x = "log2 Fold Change",
    y = "-log10 adjusted p-value"
  )

ggsave("results/figures/volcano_plot.png", width = 7, height = 5, dpi = 300)




### HEATMAP

library(DESeq2)
library(pheatmap)
library(readr)
library(dplyr)

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

vsd <- readRDS("results/deseq2/vsd.rds")
res <- read_csv("results/deseq2/deseq2_results_all.csv", show_col_types = FALSE)

top_genes <- res %>%
  filter(!is.na(padj)) %>%
  arrange(padj) %>%
  slice_head(n = 30) %>%
  pull(gene)

mat <- assay(vsd)[top_genes, , drop = FALSE]
mat <- mat - rowMeans(mat)

anno <- as.data.frame(colData(vsd)[, "condition", drop = FALSE])

pheatmap(
  mat,
  annotation_col = anno,
  show_rownames = TRUE,
  fontsize_row = 8,
  filename = "results/figures/heatmap_top30_genes.png",
  width = 8,
  height = 10
)



