# ============================================================
# DESeq2 analysis
# 6-month-old WT vs KO
# ============================================================

# Load package
library(DESeq2)


# ------------------------------------------------------------
# 1. Load sample information
# ------------------------------------------------------------

samples <- read.csv("2023/sample_table.csv")

samples$condition <- factor(
    samples$condition,
    levels = c("WT", "KO")
)


# ------------------------------------------------------------
# 2. Load count matrix
# ------------------------------------------------------------

counts_df <- read.csv(
    "2023/6month_count.csv",
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

# Upregulated genes
up <- subset(
    as.data.frame(res),
    padj < 0.05 & log2FoldChange > 1
)

# Downregulated genes
down <- subset(
    as.data.frame(res),
    padj < 0.05 & log2FoldChange < -1
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
    "2023/6month_deg.csv",
    row.names = FALSE
)

write.csv(
    up,
    "2023/6month_updeg.csv",
    row.names = FALSE
)

write.csv(
    down,
    "2023/6month_downdeg.csv",
    row.names = FALSE
)
