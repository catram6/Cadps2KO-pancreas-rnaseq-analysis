# ============================================================
# DESeq2 analysis
# WT: 6-month-old vs 14-month-old
# ============================================================

# Load package
library(DESeq2)


# ------------------------------------------------------------
# 1. Load sample information
# ------------------------------------------------------------

samples <- read.csv("WT/sample_table.csv")

samples$condition <- factor(
    samples$condition,
    levels = c("6month", "14month")
)


# ------------------------------------------------------------
# 2. Load count matrix
# ------------------------------------------------------------

counts_df <- read.csv(
    "WT/WT_count.csv",
    row.names = 1
)

count_matrix <- as.matrix(counts_df)

storage.mode(count_matrix) <- "numeric"


# ------------------------------------------------------------
# 3. Remove rows containing NA
# ------------------------------------------------------------

count_matrix <- na.omit(count_matrix)


# ------------------------------------------------------------
# 4. Create DESeq2 dataset
# ------------------------------------------------------------

dds <- DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = samples,
    design = ~ condition
)


# ------------------------------------------------------------
# 5. Filter low-count genes
# ------------------------------------------------------------

keep <- rowSums(counts(dds)) >= 10

dds <- dds[keep, ]


# ------------------------------------------------------------
# 6. Run DESeq2
# ------------------------------------------------------------

dds <- DESeq(dds)


# ------------------------------------------------------------
# 7. Extract differential expression results
# ------------------------------------------------------------

res <- results(dds)


# ------------------------------------------------------------
# 8. Identify DEGs
# ------------------------------------------------------------

# log2FC threshold was set to 2 because
# the number of DEGs exceeded 3,000 with a threshold of 1,
# which prevented downstream Metascape analysis.

# Upregulated genes
up <- subset(
    as.data.frame(res),
    padj < 0.05 & log2FoldChange > 2
)

# Downregulated genes
down <- subset(
    as.data.frame(res),
    padj < 0.05 & log2FoldChange < -2
)


# ------------------------------------------------------------
# 9. Add gene IDs
# ------------------------------------------------------------

res$gene <- rownames(res)

up$gene <- rownames(up)

down$gene <- rownames(down)


# ------------------------------------------------------------
# 10. Save results
# ------------------------------------------------------------

write.csv(
    res,
    "WT/WT_deg.csv",
    row.names = FALSE
)

write.csv(
    up,
    "WT/WT_updeg.csv",
    row.names = FALSE
)

write.csv(
    down,
    "WT/WT_downdeg.csv",
    row.names = FALSE
)
