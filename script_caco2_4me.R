#---- 4-ME, TRANSCRIPTOMICS CODE ----
#Packages 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("tximport")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("DOSE") 

#Read files
library(readr)
RawDat <- read_tsv("/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/Laura_JB_caco2_RNA-seq.tsv", col_names = TRUE)
counts <- RawDat[,-c(1,2)]
gene_info <- RawDat[,c(1:2)]
counts <- as.data.frame(counts)
coldata <-read_tsv("/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/Laura_JB_caco2_RNA-seq_meta.tsv", col_names = TRUE)
coldata <- as.data.frame(coldata)
rownames(coldata) <- coldata$sample
coldata$group <- as.factor(coldata$group)
coldata$replicate <- as.factor(coldata$replicate)

#Samples 
files <- c("I223" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I223/quant.sf",
           "I224" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I224/quant.sf",
           "I225" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I225/quant.sf",
           "I226" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I226/quant.sf",
           "I227" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I227/quant.sf",
           "I228" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I228/quant.sf",
           "I229" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I229/quant.sf",
           "I230" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I230/quant.sf",
           "I231" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I231/quant.sf",
           "I232" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I232/quant.sf",
           "I233" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I233/quant.sf",
           "I234" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I234/quant.sf",
           "I235" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I235/quant.sf",
           "I236" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I236/quant.sf",
           "I237" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I237/quant.sf",
           "I238" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I238/quant.sf",
           "I239" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I239/quant.sf",
           "I240" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I240/quant.sf",
           "I241" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I241/quant.sf",
           "I242" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I242/quant.sf",
           "I243" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I243/quant.sf",
           "I244" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I244/quant.sf",
           "I245" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I245/quant.sf",
           "I246" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I246/quant.sf")

library(tximport)

#Data import, transcript level quantification using Salmon
txi <- tximport(files,
                type = "salmon",
                txOut = TRUE,
                countsFromAbundance = "lengthScaledTPM")


#Round to integer for DESeq2 analysis
counts_int <- round(txi$counts)

library(DESeq2)

#DESeq2 using groups as experimental design
dds2 <- DESeqDataSetFromMatrix(countData = counts_int,
                               colData = coldata,
                               design = ~ group)

#Keeps only genes with sample counts over >= 10
keep <- apply(counts(dds2), 1, function(x) all(x >= 10))
dds2_filtered <- dds2[keep, ]

#DESeq2, size factor estimation, dispersion estimation and GLM fitting.
dds2 <- DESeq(dds2_filtered)

#Log transformation
rld <- rlog(dds2, blind = FALSE)

#PCA plot using Group factor
pcaPlot <- plotPCA(rld, intgroup = "group")
pcaPlot

#3 bad replicates with extreme outliers are observed, these files need to be removed
tx2gene <- RawDat[,c(1,2)]

files <- c("I223" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I223/quant.sf",
           "I224" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I224/quant.sf",
           "I225" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I225/quant.sf",
           "I226" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I226/quant.sf",
           "I227" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I227/quant.sf",
           "I229" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I229/quant.sf",
           "I230" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I230/quant.sf",
           "I231" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I231/quant.sf",
           "I232" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I232/quant.sf",
           "I233" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I233/quant.sf",
           "I234" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I234/quant.sf",
           "I235" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I235/quant.sf",
           "I237" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I237/quant.sf",
           "I238" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I238/quant.sf",
           "I239" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I239/quant.sf",
           "I240" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I240/quant.sf",
           "I241" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I241/quant.sf",
           "I243" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I243/quant.sf",
           "I244" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I244/quant.sf",
           "I245" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I245/quant.sf",
           "I246" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I246/quant.sf")

#Import
txi <- tximport(files,
                type = "salmon",
                tx2gene = tx2gene,
                countsFromAbundance = "lengthScaledTPM")

#Removal of faulty replicates from the meta data (per. the PCA analysis)
coldata_mod <- coldata[-c(6,14,20),]
raw_counts_df <- as.data.frame(txi$counts)

#Add a column with gene identifiers
raw_counts_df$gene <- rownames(raw_counts_df)

library(reshape2)

#Reshape the data from wide to long format
raw_counts_long <- melt(raw_counts_df, id.vars = "gene", 
                        variable.name = "sample", 
                        value.name = "counts")
library(ggplot2)

#Boxplot with a log2 transformation
ggplot(raw_counts_long, aes(x = sample, y = log2(counts))) +
  geom_boxplot(fill = "lightgreen") +
  labs(title = "Boxplot of Raw Counts (log2-transformed)",
       x = "Samples",
       y = "log2(raw counts)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#As this boxplot shows a need for normalization, the DESeq2 Pipeline is continued
counts_int <- round(txi$counts)
dds <- DESeqDataSetFromMatrix(countData = counts_int,
                              colData = coldata_mod,
                              design = ~ group)
keep <- apply(counts(dds), 1, function(x) all(x >= 10))
dds_filtered <- dds[keep, ]

#Running DESeq2
dds <- DESeq(dds_filtered)

normalized_counts <- counts(dds,normalized=TRUE)
norm_counts_df <- as.data.frame(normalized_counts)

#Adding Gene Identifiers
norm_counts_df$gene <- rownames(normalized_counts)

#Reshape the data from wide to long format
norm_counts_long <- melt(norm_counts_df, id.vars = "gene", 
                         variable.name = "sample", 
                         value.name = "counts")

#Create the boxplot with a log2 transformation
ggplot(norm_counts_long, aes(x = sample, y = log2(counts))) +
  geom_boxplot(fill = "lightblue") +
  labs(title = "Boxplot of Normalized Counts (log2-transformed)",
       x = "Samples",
       y = "log2(Normalized counts)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#The data is now normalized

#PCA plot without the faulty replicates
rld <- rlog(dds, blind = FALSE)
rld_mat <- assay(rld)
pca <- prcomp(t(rld_mat))
variance_explained <- (pca$sdev^2 / sum(pca$sdev^2)) * 100

#Dataframe containging two principal components and metadata
pca_df <- data.frame(PC1 = pca$x[, 1],
                     PC2 = pca$x[, 2],
                     sample = colnames(rld),
                     group = rld$group,       # assumes 'group' exists in colData(rld)
                     replicate = rld$replicate)  # likewise for 'replicate'

library(ggrepel)

ggplot(pca_df, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = replicate), size = 4, box.padding = 0.35, point.padding = 0.5) +
  xlab(paste0("PC1: ", round(variance_explained[1], 1), "% variance")) +
  ylab(paste0("PC2: ", round(variance_explained[2], 1), "% variance")) +
  ggtitle("PCA of RNA-seq Samples") +
  theme_bw()

#Log2 Fold Change calculation
res_1 <- results(dds, contrast = c("group", "Caco2-4ME-short", "Caco2-ctrl"))
res_2 <- results(dds, contrast = c("group", "Caco2-4ME-long", "Caco2-ctrl"))
res_3 <- results(dds, contrast = c("group", "Caco2-LPS", "Caco2-ctrl"))
res_4 <- results(dds, contrast = c("group", "Caco2-4ME+LPS", "Caco2-ctrl"))
res_5 <- results(dds, contrast = c("group", "Caco2-4ME+LPS+4ME", "Caco2-ctrl"))

#LOG2 FC BY LPS
res_4_lc <- results(dds, contrast = c("group", "Caco2-4ME+LPS", "Caco2-LPS"))
res_5_lc <- results(dds, contrast = c("group", "Caco2-4ME+LPS+4ME", "Caco2-LPS"))

# ---- TIGHT JUNCTION GENES ----

library(tibble)
library(dplyr)
library(tidyr)
library(forcats)

#Using paste to find all claudins and magi
cldn <- paste0("CLDN", 1:27)
magi <- paste0("MAGI", 1:3) 

#Using lists to orgnaize protein-coding genes
genes_tj <- list(
  `Transmembrane proteins` = c(
    cldn,              
    "OCLN",           
    "F11R",           
    "MARVELD2",        
    "MARVELD3"       
  ),
  `Adapter- and scaffold proteins` = c(
    "TJP1","TJP2","TJP3",  
    "CGN","CGNL1"         
  ),
  `Other associated proteins` = c(
    "INADL",              
    "MPDZ",           
    magi                  
  )
)

#Building a list to look up specific genes
canonicalize <- function(x) {
  x <- gsub("[[:space:]-]+", "", x)
  x <- ifelse(x == "PATJ", "INADL", x)
  x
}
tj_df <- enframe(genes_tj, name = "category", value = "gene") %>%
  unnest(gene) %>%
  mutate(gene = canonicalize(gene)) %>%
  distinct(category, gene)
genes_of_interest_tj <- unique(tj_df$gene)


#DE-results 
df1 <- as.data.frame(res_1) %>% rownames_to_column("gene") %>%
  mutate(gene = canonicalize(gene), dataset = "2H 4-ME")
df2 <- as.data.frame(res_2) %>% rownames_to_column("gene") %>%
  mutate(gene = canonicalize(gene), dataset = "22H 4-ME")
df3 <- as.data.frame(res_3) %>% rownames_to_column("gene") %>%
  mutate(gene = canonicalize(gene), dataset = "20H LPS")
df4 <- as.data.frame(res_4) %>% rownames_to_column("gene") %>%
  mutate(gene = canonicalize(gene), dataset = "Combined, Short")
df5 <- as.data.frame(res_5) %>% rownames_to_column("gene") %>%
  mutate(gene = canonicalize(gene), dataset = "Combined, Long")

#DE-resultsm changes log2 (LPS as baseline)
df4_lc <- as.data.frame(res_4_lc) %>% rownames_to_column("gene") %>%
  mutate(gene = canonicalize(gene), dataset = "Combined, Short")
df5_lc <- as.data.frame(res_5_lc) %>% rownames_to_column("gene") %>%
  mutate(gene = canonicalize(gene), dataset = "Combined, Long")

#Binding rows
tj_all <- dplyr::bind_rows(df1, df2, df3, df4, df5) %>%
  dplyr::inner_join(tj_df, by = "gene") %>%
  dplyr::select(category, gene, dataset, log2FoldChange, padj)

tj_changed <- dplyr::bind_rows(df4_lc, df5_lc) %>%
  dplyr::inner_join(tj_df, by = "gene") %>%
  dplyr::select(category, gene, dataset, log2FoldChange, padj)

#HEATMAPS

#ALL
plot_df_all_tj <- tj_all %>%
  mutate(
    star = ifelse(!is.na(padj) & padj < 0.05, "*", ""),
    dataset = factor(dataset, levels = c("2H 4-ME", "22H 4-ME", "20H LPS", "Combined, Short", "Combined, Long"))
  ) %>%
  group_by(category) %>%
  mutate(gene_f = fct_reorder(gene, abs(log2FoldChange), .fun = max, .desc = TRUE)) %>%
  ungroup()
maxval <- max(abs(tj_all$log2FoldChange), na.rm = TRUE)
ggplot(plot_df_all_tj, aes(x = dataset, y = gene_f, fill = log2FoldChange)) +
  geom_tile(color = "grey90") +
  geom_text(aes(label = star), color = "black", size = 3.5) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 0, limits = c(-maxval, maxval),
                       na.value = "grey90") +
  facet_wrap(~ category, scales = "free_y", ncol = 1) +
  labs(title = "Tight Junction-Related Genes, All Conditions",
       x = "Dataset", y = "Gene", fill = "log2 Fold Change") +
  theme_minimal(base_size = 16) +
  theme(axis.text.x = element_text(size = 14, hjust = 0.5),
        strip.text = element_text(face = "bold"))


#EXTRA PLOT
plot_df_changed <- tj_changed %>%
  mutate(
    star = ifelse(!is.na(padj) & padj < 0.05, "*", ""),
    dataset = factor(dataset, levels = c("Combined, Short", "Combined, Long"))
  ) %>%
  group_by(category) %>%
  mutate(gene_f = fct_reorder(gene, abs(log2FoldChange), .fun = max, .desc = TRUE)) %>%
  ungroup()
maxval <- max(abs(plot_df_changed$log2FoldChange), na.rm = TRUE)
ggplot(plot_df_changed, aes(x = dataset, y = gene_f, fill = log2FoldChange)) +
  geom_tile(color = "grey90") +
  geom_text(aes(label = star), color = "black", size = 3.5) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 0, limits = c(-maxval, maxval),
                       na.value = "grey90") +
  facet_wrap(~ category, scales = "free_y", ncol = 1) +
  labs(title = "Tight Junction-Related Genes, Log2 Ratio Relative to LPS",
       x = "Dataset", y = "Gene", fill = "log2 Fold Change") +
  theme_minimal(base_size = 16) +
  theme(axis.text.x = element_text(size = 14, hjust = 0.5),
        strip.text = element_text(face = "bold"))


# ---- MAPK ----

library(org.Hs.eg.db)

#Combined as a list
dfs <- list(
  "2H 4-ME"       = df1,
  "22H 4-ME"      = df2,
  "20H LPS"       = df3,
  "Combined, Short" = df4,
  "Combined, Long" = df5
)

dfs_lc <- list(
  "Combined, Short"       = df4_lc,
  "Combined, Long"      = df5_lc
)

#ALL
get_kegg_symbols <- function(kegg_id = "hsa04010") {
  kg <- KEGGREST::keggGet(kegg_id)[[1]]$GENE
  if (is.null(kg) || !length(kg)) return(character())
  entrez <- sub("hsa:", "", kg[seq(1, length(kg), by = 2)])
  AnnotationDbi::select(org.Hs.eg.db, keys = entrez,
                        keytype = "ENTREZID", columns = "SYMBOL") %>%
    dplyr::pull(SYMBOL) %>% unique() %>% toupper()
}
mapk_syms <- get_kegg_symbols("hsa04010")

#ALL
plot_df_mapk_all <- bind_rows(dfs, .id = "dataset") %>%
  mutate(
    dataset    = factor(dataset, levels = names(dfs)),
    gene_up    = toupper(gene),
    star       = ifelse(!is.na(padj) & padj < 0.05, "*", "")
  ) %>%
  filter(gene_up %in% mapk_syms)

#Filtering the top genes
top_genes <- plot_df_mapk_all %>%
  group_by(gene) %>%
  summarise(maxabs = max(abs(log2FoldChange), na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(maxabs)) %>%
  slice_head(n = 30) %>%
  pull(gene)

plot_df_mapk_all <- plot_df_mapk_all %>%
  filter(gene %in% top_genes) %>%
  mutate(gene = factor(gene, levels = rev(top_genes)))
clim <- max(abs(plot_df_mapk_all$log2FoldChange), na.rm = TRUE); if (!is.finite(clim)) clim <- NA
ggplot(plot_df_mapk_all, aes(x = dataset, y = gene, fill = log2FoldChange)) +
  geom_tile(color = "grey90") +
  geom_text(aes(label = star), size = 3) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 0, limits = c(-clim, clim), oob = scales::squish) +
  labs(title = "Differential Expression of MAPK Genes, All Conditions",
       x = NULL, y = "Gene", fill = "log2FC") +
  theme_minimal(base_size = 16) +
  theme(axis.text.x = element_text(size = 14, hjust = 0.5),
        panel.grid.minor = element_blank())

#EXTRA PLOT, LOG2FC RELATIVE TO LPS
plot_df_lc <- bind_rows(dfs_lc, .id = "dataset") %>%
  mutate(
    dataset    = factor(dataset, levels = names(dfs_lc)),
    gene_up    = toupper(gene),
    star       = ifelse(!is.na(padj) & padj < 0.05, "*", "")
  ) %>%
  filter(gene_up %in% mapk_syms)

#Filtering the top genes
top_genes_lc <- plot_df_lc %>%
  group_by(gene) %>%
  summarise(maxabs = max(abs(log2FoldChange), na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(maxabs)) %>%
  slice_head(n = 30) %>%
  pull(gene)

plot_df_lc <- plot_df_lc %>%
  filter(gene %in% top_genes_lc) %>%
  mutate(gene = factor(gene, levels = rev(top_genes_lc)))
clim <- max(abs(plot_df_lc$log2FoldChange), na.rm = TRUE); if (!is.finite(clim)) clim <- NA
ggplot(plot_df_lc, aes(x = dataset, y = gene, fill = log2FoldChange)) +
  geom_tile(color = "grey90") +
  geom_text(aes(label = star), size = 3) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 0, limits = c(-clim, clim), oob = scales::squish) +
  labs(title = "Differential Expression of MAPK Genes, Log2 Ratio Relative to LPS",
       x = NULL, y = "Gene", fill = "log2FC") +
  theme_minimal(base_size = 16) +
  theme(axis.text.x = element_text(size = 14, hjust = 0.5),
        panel.grid.minor = element_blank())

#---- NF-kB ----

kegg_top_df <- function(dfs, kegg_id, top_n = 30, padj_cut = 0.05, sig_only = FALSE) {
  # KEGG -> ENTREZ -> SYMBOL (uppercase for robust matching)
  kg <- KEGGREST::keggGet(kegg_id)[[1]]$GENE
  if (is.null(kg) || !length(kg)) return(tibble())
  entrez <- sub("^hsa:", "", kg[seq(1, length(kg), by = 2)])
  syms <- AnnotationDbi::select(org.Hs.eg.db, keys = entrez,
                                keytype = "ENTREZID", columns = "SYMBOL") |>
    dplyr::pull(SYMBOL) |> unique() |> toupper()
  out <- dplyr::bind_rows(dfs, .id = "dataset") |>
    dplyr::mutate(
      dataset = factor(dataset, levels = names(dfs)),
      gene_up = toupper(.data$gene),
      sig     = !is.na(.data$padj) & .data$padj < padj_cut,
      star    = ifelse(sig, "*", "")
    ) |>
    dplyr::filter(.data$gene_up %in% syms)
  if (sig_only) out <- dplyr::filter(out, .data$sig)
  if (!nrow(out)) return(out)
  ord <- out |>
    dplyr::group_by(.data$gene) |>
    dplyr::summarise(maxabs = max(abs(.data$log2FoldChange), na.rm = TRUE), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(.data$maxabs)) |>
    dplyr::slice_head(n = top_n) |>
    dplyr::pull(.data$gene)
  
  out |>
    dplyr::filter(.data$gene %in% ord) |>
    dplyr::mutate(gene = factor(.data$gene, levels = rev(ord)))
}

#ALL
plt_nfkb_all <- kegg_top_df(dfs, kegg_id = "hsa04064", top_n = 30, padj_cut = 0.05, sig_only = FALSE)
if (nrow(plt_nfkb_all)) {
  clim <- max(abs(plt_nfkb_all$log2FoldChange), na.rm = TRUE)
  ggplot(plt_nfkb_all, aes(x = dataset, y = gene, fill = log2FoldChange)) +
    geom_tile(color = "grey90") +
    geom_text(aes(label = star), size = 3) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                         midpoint = 0, limits = c(-clim, clim), oob = scales::squish) +
    labs(title = "Differential Expression of NF-κB Genes, All Conditions",
         x = NULL, y = "Gene", fill = "log2FC") +
    theme_minimal(base_size = 16) +
    theme(axis.text.x = element_text(size = 14, hjust = 0.5),
          panel.grid.minor = element_blank())
}

#2H 4-ME, 20H LPS & Combined
plt_nfkb_lc <- kegg_top_df(dfs_lc, kegg_id = "hsa04064", top_n = 30, padj_cut = 0.05, sig_only = FALSE)
if (nrow(plt_nfkb_lc)) {
  clim <- max(abs(plt_nfkb_lc$log2FoldChange), na.rm = TRUE)
  ggplot(plt_nfkb_lc, aes(x = dataset, y = gene, fill = log2FoldChange)) +
    geom_tile(color = "grey90") +
    geom_text(aes(label = star), size = 3) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                         midpoint = 0, limits = c(-clim, clim), oob = scales::squish) +
    labs(title = "Differential Expression of NF-κB Genes, Log2 Ratio Relative to LPS",
         x = NULL, y = "Gene", fill = "log2FC") +
    theme_minimal(base_size = 16) +
    theme(axis.text.x = element_text(size = 14, hjust = 0.5),
          panel.grid.minor = element_blank())
}

#---- MAPK & NF-kB OVERLAPPING GENES ----

#Overlapping genes by KEGG
mapk_symbols <- toupper(get_kegg_symbols("hsa04010"))
nfkb_symbols <- toupper(get_kegg_symbols("hsa04064"))
overlap_syms <- intersect(mapk_symbols, nfkb_symbols)

ibm_cols <- c(
  "2H 4-ME"      = "#0F62FE", 
  "20H LPS"      = "#8A3FFC",  
  "Combined, Short"     = "#EE5396",  
  "22H 4-ME"     = "#FF832B",  
  "Combined, Long"= "#F1C21B"   
)


#ALL
plot_df_overlap_all <- bind_rows(dfs, .id = "dataset") %>%
  mutate(
    dataset = factor(dataset, levels = names(dfs)),
    gene_up = toupper(gene),
    star    = ifelse(!is.na(padj) & padj < 0.05, "*", "")
  ) %>%
  filter(gene_up %in% overlap_syms)

if (nrow(plot_df_overlap_all) == 0) {
  message("Ingen overlap-gener fundet i dine datasets.")
} else {
  ord <- plot_df_overlap_all %>%
    group_by(gene) %>%
    summarise(maxabs = max(abs(log2FoldChange), na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(maxabs)) %>%
    pull(gene)
  
  plot_df_overlap_all <- plot_df_overlap_all %>%
    mutate(
      gene    = factor(gene, levels = rev(ord)),
      dataset = fct_drop(dataset) 
    )
  
  pal <- ibm_cols[levels(plot_df_overlap_all$dataset)]
  maxval <- max(abs(plot_df_overlap_all$log2FoldChange), na.rm = TRUE)
  off    <- 0.04 * maxval
  ggplot(plot_df_overlap_all, aes(x = gene, y = log2FoldChange, fill = dataset)) +
    geom_hline(yintercept = 0, colour = "grey75") +
    geom_col(position = position_dodge(width = 0.8), width = 0.72) +
    geom_text(aes(label = star,
                  y = ifelse(log2FoldChange >= 0, log2FoldChange + off, log2FoldChange - off)),
              position = position_dodge(width = 0.8), size = 3.2) +
    coord_flip() +
    scale_fill_manual(values = pal, drop = TRUE) +
    scale_y_continuous(limits = c(-maxval, maxval),
                       expand = expansion(mult = c(0.02, 0.08))) +
    labs(title = "MAPK ∩ NF-κB Overlapping Genes",
         x = "Gene", y = "log2 Fold Change", fill = "Dataset") +
    theme_minimal(base_size = 16)
}

#HEATMAP
heat_df_all <- plot_df_overlap_all %>%
  mutate(dataset = fct_drop(dataset))  
clim <- max(abs(heat_df_all$log2FoldChange), na.rm = TRUE)
ggplot(heat_df_all, aes(x = dataset, y = gene, fill = log2FoldChange)) +
  geom_tile(color = "grey90") +
  geom_text(aes(label = star), size = 3) +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0,
    limits = c(-clim, clim), oob = scales::squish, na.value = "grey90"
  ) +
  scale_x_discrete(drop = TRUE) +
  labs(title = "MAPK ∩ NF-κB Overlapping Genes",
       x = NULL, y = "Gene", fill = "log2FC") +
  theme_minimal(base_size = 16) +
  theme(axis.text.x = element_text(size = 14, hjust = 0.5))

#CHANGED LOG
plot_df_overlap_lc <- bind_rows(dfs_lc, .id = "dataset") %>%
  mutate(
    dataset = factor(dataset, levels = names(dfs_lc)),
    gene_up = toupper(gene),
    star    = ifelse(!is.na(padj) & padj < 0.05, "*", "")
  ) %>%
  filter(gene_up %in% overlap_syms)

if (nrow(plot_df_overlap_lc) == 0) {
  message("Ingen overlap-gener fundet i dine datasets.")
} else {
  ord <- plot_df_overlap_lc %>%
    group_by(gene) %>%
    summarise(maxabs = max(abs(log2FoldChange), na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(maxabs)) %>%
    pull(gene)
  
  plot_df_overlap_lc <- plot_df_overlap_lc %>%
    mutate(
      gene    = factor(gene, levels = rev(ord)),
      dataset = fct_drop(dataset) 
    )
  
  pal <- ibm_cols[levels(plot_df_overlap_lc$dataset)]
  maxval <- max(abs(plot_df_overlap_lc$log2FoldChange), na.rm = TRUE)
  off    <- 0.04 * maxval
  ggplot(plot_df_overlap_lc, aes(x = gene, y = log2FoldChange, fill = dataset)) +
    geom_hline(yintercept = 0, colour = "grey75") +
    geom_col(position = position_dodge(width = 0.8), width = 0.72) +
    geom_text(aes(label = star,
                  y = ifelse(log2FoldChange >= 0, log2FoldChange + off, log2FoldChange - off)),
              position = position_dodge(width = 0.8), size = 3.2) +
    coord_flip() +
    scale_fill_manual(values = pal, drop = TRUE) +
    scale_y_continuous(limits = c(-maxval, maxval),
                       expand = expansion(mult = c(0.02, 0.08))) +
    labs(title = "MAPK ∩ NF-κB Overlapping Genes, Log2 Ratio Relative to LPS",
         x = "Gene", y = "log2 Fold Change", fill = "Dataset") +
    theme_minimal(base_size = 16)
}

#HEATMAP
heat_df_lc <- plot_df_overlap_lc %>%
  mutate(dataset = fct_drop(dataset)) 
clim <- max(abs(heat_df_lc$log2FoldChange), na.rm = TRUE)
ggplot(heat_df_lc, aes(x = dataset, y = gene, fill = log2FoldChange)) +
  geom_tile(color = "grey90") +
  geom_text(aes(label = star), size = 3) +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0,
    limits = c(-clim, clim), oob = scales::squish, na.value = "grey90"
  ) +
  scale_x_discrete(drop = TRUE) +
  labs(title = "MAPK ∩ NF-κB Overlapping Genes, Log2 Ratio Relative to LPS",
       x = NULL, y = "Gene", fill = "log2FC") +
  theme_minimal(base_size = 16) +
  theme(axis.text.x = element_text(size = 14, hjust = 0.5))


#---- ROS & GSH (MITOCARTA) ----
library(readxl)
library(tidyverse)
library(purrr)
mito_genes <- readxl::read_excel("~/Desktop/Human.MitoCarta3.0.xls", sheet = "A Human MitoCarta3.0")
mito_symbols <- mito_genes$Symbol %>% unique()
mito_pathways <- readxl::read_excel("~/Desktop/Human.MitoCarta3.0.xls", sheet = "C MitoPathways")
library(stringr)
ros_gsh_syms <- mito_pathways %>%
  filter(!is.na(MitoPathway),
         str_detect(MitoPathway,
                    regex("reactive\\s*oxygen|\\bROS\\b|glutathione", ignore_case = TRUE)),
         !is.na(Genes)) %>%
  pull(Genes) %>%
  paste(collapse = ",") %>%                   
  str_split("[,;|]") %>% `[[`(1) %>%
  str_trim() %>%
  str_replace_all("[^A-Za-z0-9_-]", "") %>%   
  toupper() %>%
  { .[. != ""] } %>%
  unique()
if (!length(ros_gsh_syms)) stop("Ingen ROS/Glutathione-symboler fundet i MitoPathways.")

long_all <- bind_rows(dfs, .id = "dataset") %>%
  mutate(
    dataset    = factor(dataset, levels = names(dfs)),
    gene_upper = toupper(gene),
    in_set     = gene_upper %in% ros_gsh_syms,
    sig        = !is.na(padj) & padj < 0.05
  )

long_lc <- bind_rows(dfs_lc, .id = "dataset") %>%
  mutate(
    dataset    = factor(dataset, levels = names(dfs_lc)),
    gene_upper = toupper(gene),
    in_set     = gene_upper %in% ros_gsh_syms,
    sig        = !is.na(padj) & padj < 0.05
  )

plot_df_all <- long_all %>%
  filter(in_set) %>%
  mutate(star = ifelse(sig, "*", ""))
if (nrow(plot_df_all)) {
  ord <- plot_df_all %>%
    group_by(gene) %>%
    summarise(maxabs = max(abs(log2FoldChange), na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(maxabs)) %>% pull(gene)
  plot_df_all <- plot_df_all %>% mutate(gene = factor(gene, levels = rev(ord)))
  clim <- max(abs(plot_df_all$log2FoldChange), na.rm = TRUE); if (!is.finite(clim)) clim <- NA
  ggplot(plot_df_all, aes(x = dataset, y = gene, fill = log2FoldChange)) +
    geom_tile(color = "grey85") +
    geom_text(aes(label = star), size = 3.3) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                         midpoint = 0, limits = c(-clim, clim),
                         oob = scales::squish) +
    scale_x_discrete(position = "bottom") +
    labs(title = "ROS & Glutathione Metabolism",
         x = NULL, y = "Gene", fill = "log2FC") +
    theme_minimal(base_size = 16) +
    theme(axis.text.x = element_text(size = 14, hjust = 0.5))
} else {
  message("Ingen ROS/GSH pathway-gener matcher dine res-tabeller.")
}

#CHANGED LOG
plot_df_lc <- long_lc %>%
  filter(in_set) %>%
  mutate(star = ifelse(sig, "*", ""))
if (nrow(plot_df_lc)) {
  ord <- plot_df_lc %>%
    group_by(gene) %>%
    summarise(maxabs = max(abs(log2FoldChange), na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(maxabs)) %>% pull(gene)
  plot_df_lc <- plot_df_lc %>% mutate(gene = factor(gene, levels = rev(ord)))
  clim <- max(abs(plot_df_lc$log2FoldChange), na.rm = TRUE); if (!is.finite(clim)) clim <- NA
  ggplot(plot_df_lc, aes(x = dataset, y = gene, fill = log2FoldChange)) +
    geom_tile(color = "grey85") +
    geom_text(aes(label = star), size = 3.3) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                         midpoint = 0, limits = c(-clim, clim),
                         oob = scales::squish) +
    scale_x_discrete(position = "bottom") +
    labs(title = "ROS & Glutathione Metabolism, Log2 Ratio Relative to LPS",
         x = NULL, y = "Gene", fill = "log2FC") +
    theme_minimal(base_size = 16) +
    theme(axis.text.x = element_text(size = 14, hjust = 0.5))
} else {
  message("Ingen ROS/GSH pathway-gener matcher dine res-tabeller.")
}

#---- RESPONSE TO OXIDATIVE STRESS ----

go_sig_df <- function(dfs, go_id = "GO:0006979", padj_cut = 0.05) {
  goset <- AnnotationDbi::select(org.Hs.eg.db, keys = go_id, keytype = "GO", columns = "ENTREZID") |>
    dplyr::filter(!is.na(ENTREZID)) |>
    dplyr::distinct(ENTREZID) |>
    dplyr::pull(ENTREZID)
  gosym <- AnnotationDbi::select(org.Hs.eg.db, keys = goset, keytype = "ENTREZID", columns = "SYMBOL") |>
    dplyr::filter(!is.na(SYMBOL)) |>
    dplyr::distinct(SYMBOL) |>
    dplyr::pull(SYMBOL) |>
    toupper()
  out <- dplyr::bind_rows(dfs, .id = "dataset") |>
    dplyr::mutate(
      dataset = factor(dataset, levels = names(dfs)),
      gene_up = toupper(gene),
      sig     = !is.na(padj) & padj < padj_cut,
      star    = ifelse(sig, "*", "")
    ) |>
    dplyr::filter(gene_up %in% gosym, sig)
  if (nrow(out) == 0) return(out)
  ord <- out |>
    dplyr::group_by(gene) |>
    dplyr::summarise(mx = max(abs(log2FoldChange), na.rm = TRUE), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(mx)) |>
    dplyr::pull(gene)
  
  out |>
    dplyr::mutate(gene = factor(gene, levels = rev(ord)))
}

plt_all_os <- go_sig_df(dfs, go_id = "GO:0006979", padj_cut = 0.05)

if (nrow(plt_all_os) == 0) {
  message("Ingen signifikante oxidative-stress gener fundet.")
} else {
  clim <- max(abs(plt_all_os$log2FoldChange), na.rm = TRUE); if (!is.finite(clim)) clim <- NA
  ggplot2::ggplot(plt_all_os, ggplot2::aes(x = dataset, y = gene, fill = log2FoldChange)) +
    ggplot2::geom_tile(color = "grey80") +
    ggplot2::geom_text(ggplot2::aes(label = star), size = 3.4) +
    ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                                  midpoint = 0, limits = c(-clim, clim),
                                  oob = scales::squish) +
    ggplot2::scale_x_discrete(position = "bottom") +
    ggplot2::labs(title = "Response to Oxidative Stress, All Conditions",
                  x = NULL, y = "Gene", fill = "log2FC") +
    ggplot2::theme_minimal(base_size = 15) +
    theme(axis.text.x = element_text(size = 14, hjust = 0.5))
}

#CHANGED LOG2
go_sig_df_lc <- function(dfs_lc, go_id = "GO:0006979", padj_cut = 0.05) {
  goset <- AnnotationDbi::select(org.Hs.eg.db, keys = go_id, keytype = "GO", columns = "ENTREZID") |>
    dplyr::filter(!is.na(ENTREZID)) |>
    dplyr::distinct(ENTREZID) |>
    dplyr::pull(ENTREZID)
  gosym <- AnnotationDbi::select(org.Hs.eg.db, keys = goset, keytype = "ENTREZID", columns = "SYMBOL") |>
    dplyr::filter(!is.na(SYMBOL)) |>
    dplyr::distinct(SYMBOL) |>
    dplyr::pull(SYMBOL) |>
    toupper()
  out <- dplyr::bind_rows(dfs_lc, .id = "dataset") |>
    dplyr::mutate(
      dataset = factor(dataset, levels = names(dfs_lc)),
      gene_up = toupper(gene),
      sig     = !is.na(padj) & padj < padj_cut,
      star    = ifelse(sig, "*", "")
    ) |>
    dplyr::filter(gene_up %in% gosym, sig)
  if (nrow(out) == 0) return(out)
  ord <- out |>
    dplyr::group_by(gene) |>
    dplyr::summarise(mx = max(abs(log2FoldChange), na.rm = TRUE), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(mx)) |>
    dplyr::pull(gene)
  
  out |>
    dplyr::mutate(gene = factor(gene, levels = rev(ord)))
}

plt_all_lc <- go_sig_df_lc(dfs_lc, go_id = "GO:0006979", padj_cut = 0.05)

if (nrow(plt_all_lc) == 0) {
  message("Ingen signifikante oxidative-stress gener fundet.")
} else {
  clim <- max(abs(plt_all_lc$log2FoldChange), na.rm = TRUE); if (!is.finite(clim)) clim <- NA
  ggplot2::ggplot(plt_all_lc, ggplot2::aes(x = dataset, y = gene, fill = log2FoldChange)) +
    ggplot2::geom_tile(color = "grey80") +
    ggplot2::geom_text(ggplot2::aes(label = star), size = 3.4) +
    ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                                  midpoint = 0, limits = c(-clim, clim),
                                  oob = scales::squish) +
    ggplot2::scale_x_discrete(position = "bottom") +
    ggplot2::labs(title = "Response to Oxidative Stress, Log2 Ratio Relative to LPS",
                  x = NULL, y = "Gene", fill = "log2FC") +
    ggplot2::theme_minimal(base_size = 16) +
    theme(axis.text.x = element_text(size = 14, hjust = 0.5))
}

#---- VOLCANO PLOTS ----
#2H 4-ME
volcano_one <- function(res_1, title = "Volcano", padj_cut = 0.05, lfc_cut = 1, top_n = 10) {
  df <- as.data.frame(res_1) %>%
    rownames_to_column("gene") %>%
    filter(!is.na(padj)) %>%
    mutate(
      padj_safe = ifelse(padj <= 0, .Machine$double.xmin, padj),
      neglog10  = -log10(padj_safe),
      sig       = case_when(
        padj < padj_cut & log2FoldChange >  lfc_cut ~ "Up",
        padj < padj_cut & log2FoldChange < -lfc_cut ~ "Down",
        TRUE ~ "NS"
      )
    )
  top_genes <- df %>% arrange(padj) %>% slice_head(n = top_n)
  
  ggplot(df, aes(x = log2FoldChange, y = neglog10, color = sig)) +
    geom_point(alpha = 0.7, size = 1.6, na.rm = TRUE) +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed") +
    geom_hline(yintercept = -log10(padj_cut), linetype = "dashed") +
    ggrepel::geom_text_repel(data = top_genes, aes(label = gene), size = 5,
                             max.overlaps = 20, show.legend = FALSE) +
    scale_color_manual(values = c(Up = "#D62728", Down = "#1F77B4", NS = "grey70")) +
    labs(title = title, x = "log2 Fold Change", y = expression(-log[10]("padj")), color = NULL) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 18),   
      legend.title = element_text(size = 14),                
      legend.text  = element_text(size = 14)                 
    )
}

volcano_one(res_1, title = "Differential Expression: 2H 4-ME vs Ctrl")

#22H 4-ME
volcano_two <- function(res_2, title = "Volcano", padj_cut = 0.05, lfc_cut = 1, top_n = 10) {
  df <- as.data.frame(res_2) %>%
    rownames_to_column("gene") %>%
    filter(!is.na(padj)) %>%
    mutate(
      padj_safe = ifelse(padj <= 0, .Machine$double.xmin, padj),
      neglog10  = -log10(padj_safe),
      sig       = case_when(
        padj < padj_cut & log2FoldChange >  lfc_cut ~ "Up",
        padj < padj_cut & log2FoldChange < -lfc_cut ~ "Down",
        TRUE ~ "NS"
      )
    )
  top_genes <- df %>% arrange(padj) %>% slice_head(n = top_n)
  
  ggplot(df, aes(x = log2FoldChange, y = neglog10, color = sig)) +
    geom_point(alpha = 0.7, size = 1.6, na.rm = TRUE) +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed") +
    geom_hline(yintercept = -log10(padj_cut), linetype = "dashed") +
    ggrepel::geom_text_repel(data = top_genes, aes(label = gene), size = 5,
                             max.overlaps = 20, show.legend = FALSE) +
    scale_color_manual(values = c(Up = "#D62728", Down = "#1F77B4", NS = "grey70")) +
    labs(title = title, x = "log2 Fold Change", y = expression(-log[10]("padj")), color = NULL) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 18),   
      legend.title = element_text(size = 14),                
      legend.text  = element_text(size = 14)                 
    )
}

volcano_two(res_2, title = "Differential Expression 22H 4-ME vs Ctrl")

#LPS
volcano_three <- function(res_3, title = "Volcano", padj_cut = 0.05, lfc_cut = 1, top_n = 10) {
  df <- as.data.frame(res_3) %>%
    rownames_to_column("gene") %>%
    filter(!is.na(padj)) %>%
    mutate(
      padj_safe = ifelse(padj <= 0, .Machine$double.xmin, padj),
      neglog10  = -log10(padj_safe),
      sig       = case_when(
        padj < padj_cut & log2FoldChange >  lfc_cut ~ "Up",
        padj < padj_cut & log2FoldChange < -lfc_cut ~ "Down",
        TRUE ~ "NS"
      )
    )
  top_genes <- df %>% arrange(padj) %>% slice_head(n = top_n)
  
  ggplot(df, aes(x = log2FoldChange, y = neglog10, color = sig)) +
    geom_point(alpha = 0.7, size = 1.6, na.rm = TRUE) +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed") +
    geom_hline(yintercept = -log10(padj_cut), linetype = "dashed") +
    ggrepel::geom_text_repel(data = top_genes, aes(label = gene), size = 5,
                             max.overlaps = 20, show.legend = FALSE) +
    scale_color_manual(values = c(Up = "#D62728", Down = "#1F77B4", NS = "grey70")) +
    labs(title = title, x = "log2 Fold Change", y = expression(-log[10]("padj")), color = NULL) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 18),   
      legend.title = element_text(size = 14),                
      legend.text  = element_text(size = 14)                 
    )
}

volcano_three(res_3, title = "Differential Expression LPS vs Ctrl")

#COMBINED, SHORT
volcano_four <- function(res_4, title = "Volcano", padj_cut = 0.05, lfc_cut = 1, top_n = 10) {
  df <- as.data.frame(res_4) %>%
    rownames_to_column("gene") %>%
    filter(!is.na(padj)) %>%
    mutate(
      padj_safe = ifelse(padj <= 0, .Machine$double.xmin, padj),
      neglog10  = -log10(padj_safe),
      sig       = case_when(
        padj < padj_cut & log2FoldChange >  lfc_cut ~ "Up",
        padj < padj_cut & log2FoldChange < -lfc_cut ~ "Down",
        TRUE ~ "NS"
      )
    )
  top_genes <- df %>% arrange(padj) %>% slice_head(n = top_n)
  
  ggplot(df, aes(x = log2FoldChange, y = neglog10, color = sig)) +
    geom_point(alpha = 0.7, size = 1.6, na.rm = TRUE) +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed") +
    geom_hline(yintercept = -log10(padj_cut), linetype = "dashed") +
    ggrepel::geom_text_repel(data = top_genes, aes(label = gene), size = 5,
                             max.overlaps = 20, show.legend = FALSE) +
    scale_color_manual(values = c(Up = "#D62728", Down = "#1F77B4", NS = "grey70")) +
    labs(title = title, x = "log2 Fold Change", y = expression(-log[10]("padj")), color = NULL) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 18),   
      legend.title = element_text(size = 14),                
      legend.text  = element_text(size = 14)                 
    )
}

volcano_four(res_4, title = "Differential Expression Combined, Short vs Ctrl")

#Combined, Long 
volcano_five <- function(res_5, title = "Volcano", padj_cut = 0.05, lfc_cut = 1, top_n = 10) {
  df <- as.data.frame(res_5) %>%
    rownames_to_column("gene") %>%
    filter(!is.na(padj)) %>%
    mutate(
      padj_safe = ifelse(padj <= 0, .Machine$double.xmin, padj),
      neglog10  = -log10(padj_safe),
      sig       = case_when(
        padj < padj_cut & log2FoldChange >  lfc_cut ~ "Up",
        padj < padj_cut & log2FoldChange < -lfc_cut ~ "Down",
        TRUE ~ "NS"
      )
    )
  top_genes <- df %>% arrange(padj) %>% slice_head(n = top_n)
  
  ggplot(df, aes(x = log2FoldChange, y = neglog10, color = sig)) +
    geom_point(alpha = 0.7, size = 1.6, na.rm = TRUE) +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed") +
    geom_hline(yintercept = -log10(padj_cut), linetype = "dashed") +
    ggrepel::geom_text_repel(data = top_genes, aes(label = gene), size = 5,
                             max.overlaps = 20, show.legend = FALSE) +
    scale_color_manual(values = c(Up = "#D62728", Down = "#1F77B4", NS = "grey70")) +
    labs(title = title, x = "log2 Fold Change", y = expression(-log[10]("padj")), color = NULL) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 18),   
      legend.title = element_text(size = 14),                
      legend.text  = element_text(size = 14)                 
    )
}

volcano_five(res_5, title = "Differential Expression Combined, Long vs Ctrl")

#LOG CHANGED

#COMBINED, SHORT
volcano_four <- function(res_4_lc, title = "Volcano", padj_cut = 0.05, lfc_cut = 0.5, top_n = 10) {
  df <- as.data.frame(res_4_lc) %>%
    rownames_to_column("gene") %>%
    filter(!is.na(padj)) %>%
    mutate(
      padj_safe = ifelse(padj <= 0, .Machine$double.xmin, padj),
      neglog10  = -log10(padj_safe),
      sig       = case_when(
        padj < padj_cut & log2FoldChange >  lfc_cut ~ "Up",
        padj < padj_cut & log2FoldChange < -lfc_cut ~ "Down",
        TRUE ~ "NS"
      )
    )
  top_genes <- df %>% arrange(padj) %>% slice_head(n = top_n)
  
  ggplot(df, aes(x = log2FoldChange, y = neglog10, color = sig)) +
    geom_point(alpha = 0.7, size = 1.6, na.rm = TRUE) +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed") +
    geom_hline(yintercept = -log10(padj_cut), linetype = "dashed") +
    ggrepel::geom_text_repel(data = top_genes, aes(label = gene), size = 5,
                             max.overlaps = 20, show.legend = FALSE) +
    scale_color_manual(values = c(Up = "#D62728", Down = "#1F77B4", NS = "grey70")) +
    labs(title = title, x = "log2 Fold Change", y = expression(-log[10]("padj")), color = NULL) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 18),   
      legend.title = element_text(size = 14),                
      legend.text  = element_text(size = 14)                 
    )
}

volcano_four(res_4_lc, title = "Differential Expression Combined, Short vs LPS")

#Combined, Long 
volcano_five <- function(res_5_lc, title = "Volcano", padj_cut = 0.05, lfc_cut = 0.5, top_n = 10) {
  df <- as.data.frame(res_5_lc) %>%
    rownames_to_column("gene") %>%
    filter(!is.na(padj)) %>%
    mutate(
      padj_safe = ifelse(padj <= 0, .Machine$double.xmin, padj),
      neglog10  = -log10(padj_safe),
      sig       = case_when(
        padj < padj_cut & log2FoldChange >  lfc_cut ~ "Up",
        padj < padj_cut & log2FoldChange < -lfc_cut ~ "Down",
        TRUE ~ "NS"
      )
    )
  top_genes <- df %>% arrange(padj) %>% slice_head(n = top_n)
  
  ggplot(df, aes(x = log2FoldChange, y = neglog10, color = sig)) +
    geom_point(alpha = 0.7, size = 1.6, na.rm = TRUE) +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed") +
    geom_hline(yintercept = -log10(padj_cut), linetype = "dashed") +
    ggrepel::geom_text_repel(data = top_genes, aes(label = gene), size = 5,
                             max.overlaps = 20, show.legend = FALSE) +
    scale_color_manual(values = c(Up = "#D62728", Down = "#1F77B4", NS = "grey70")) +
    labs(title = title, x = "log2 Fold Change", y = expression(-log[10]("padj")), color = NULL) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 18),   
      legend.title = element_text(size = 14),                
      legend.text  = element_text(size = 14)                 
    )
}

volcano_five(res_5_lc, title = "Differential Expression Combined, Long vs LPS")

#---- SAMPLE-SAMPLE CORRELATION HEATMAP ----
library(pheatmap)

#Combination in long format
all_long <- bind_rows(
  dplyr::select(df1, gene, dataset, log2FoldChange),
  dplyr::select(df2, gene, dataset, log2FoldChange),
  dplyr::select(df3, gene, dataset, log2FoldChange),
  dplyr::select(df4, gene, dataset, log2FoldChange),
  dplyr::select(df5, gene, dataset, log2FoldChange)
)

#Filtering for common genes in the contrast data
keep_genes <- all_long %>%
  dplyr::distinct(gene, dataset) %>%
  dplyr::count(gene, name = "n_ds") %>%
  dplyr::filter(n_ds == dplyr::n_distinct(all_long$dataset)) %>%
  dplyr::pull(gene)

#Wide matrix
wide_mat <- all_long %>%
  dplyr::filter(gene %in% keep_genes) %>%
  dplyr::mutate(dataset = factor(dataset, levels = unique(c(df1$dataset, df2$dataset, df3$dataset, df4$dataset, df5$dataset)))) %>%
  tidyr::pivot_wider(names_from = dataset, values_from = log2FoldChange) %>%
  dplyr::arrange(gene)
mat <- as.matrix(dplyr::select(wide_mat, -gene))

#Pearson-Correlation
cor_mat <- cor(mat, method = "pearson", use = "pairwise.complete.obs")

#Heatmap
pheatmap::pheatmap(
  cor_mat,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  display_numbers = TRUE,
  number_format = "%.2f",
  color = colorRampPalette(c("blue","white","red"))(50),
  main = "Pearson correlation between contrasts",
  fontsize = 12,      
  fontsize_col = 14,     
  fontsize_row = 12,     
  border_color = NA
)

#---- PCA (Prettier) ----

pca_df <- pca_df %>%
  dplyr::mutate(
    Group = dplyr::case_when(
      group == "Caco2-ctrl"        ~ "UT",
      group == "Caco2-LPS"         ~ "20H LPS",
      group == "Caco2-4ME-short"   ~ "2H 4-ME",
      group == "Caco2-4ME-long"    ~ "22H 4-ME",
      group == "Caco2-4ME+LPS"     ~ "Combined, Short",
      group == "Caco2-4ME+LPS+4ME" ~ "Combined, Long",
      TRUE                         ~ group  
    ),
    group = factor(group,
                   levels = c("UT","2H 4-ME","22H 4-ME","20H LPS","Combined, Short","Combined, Long")
    )
  )

pca_df <- pca_df %>%
  dplyr::mutate(
    label = replicate,                             
    label = str_replace(label, "^ctrl",               "UT"),
    label = str_replace(label, "^LPS",                "20H LPS"),
    label = str_replace(label, "^4ME-short",          "2H 4-ME"),
    label = str_replace(label, "^4ME-long",           "22H 4-ME"),
    label = str_replace(label, "^4ME\\+LPS\\+4ME",    "Combined, Long"),
    label = str_replace(label, "^4ME\\+LPS",          "Combined, Short")
  )

cols <- c(
  "UT"              = "#4C78A8",
  "2H 4-ME"         = "#F58518",
  "22H 4-ME"        = "#E45756",
  "20H LPS"         = "#72B7B2",
  "Combined, Short" = "#54A24B",
  "Combined, Long"  = "#B279A2"
)

ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 5) +
  ggrepel::geom_text_repel(aes(label = label), size = 5,
                           box.padding = 0.35, point.padding = 0.5) +
  xlab(paste0("PC1: ", round(variance_explained[1], 1), "% variance")) +
  ylab(paste0("PC2: ", round(variance_explained[2], 1), "% variance")) +
  ggtitle("PCA of RNA-seq Samples") +
  scale_color_manual(values = cols, drop = FALSE) +
  theme_bw() +
  theme(
    plot.title   = element_text(size = 20, hjust = 0.5), 
    axis.title.x = element_text(size = 16),             
    axis.title.y = element_text(size = 16),          
    axis.text    = element_text(size = 14),                          
    legend.title = element_text(size = 16),       
    legend.text  = element_text(size = 14)                            
  )

#---- PEARSON CORRELATION HEATMAP ----
library(matrixStats)

rld <- vst(dds, blind = TRUE) 
mat <- assay(rld)
rv <- matrixStats::rowVars(mat)
mat <- mat[rv > 0, ]
top_n <- min(2000, nrow(mat))
top_idx <- order(rv[rv > 0], decreasing = TRUE)[seq_len(top_n)]
mat_top <- mat[top_idx, , drop = FALSE]

#Making sure labels & everything is correct
meta <- data.frame(
  sample    = colnames(mat_top),
  group_raw = as.character(colData(rld)$group),
  repl_raw  = as.character(colData(rld)$replicate),
  stringsAsFactors = FALSE
) %>%
  mutate(
    Group = case_when(
      group_raw == "Caco2-ctrl"             ~ "UT",
      group_raw == "Caco2-LPS"              ~ "20H LPS",
      group_raw == "Caco2-4ME-short"        ~ "2H 4-ME",
      group_raw == "Caco2-4ME-long"         ~ "22H 4-ME",
      group_raw == "Caco2-4ME+LPS"          ~ "Combined, Short",
      group_raw == "Caco2-4ME+LPS+4ME"      ~ "Combined, Long",
      TRUE ~ group_raw
    ),
    label = repl_raw %>%
      str_replace("^ctrl",              "UT") %>%
      str_replace("^LPS",               "20H LPS") %>%
      str_replace("^4ME-short",         "2H 4-ME") %>%
      str_replace("^4ME-long",          "22H 4-ME") %>%
      str_replace("^4ME\\+LPS\\+4ME",   "Combined, Long") %>%
      str_replace("^4ME\\+LPS",         "Combined, Short")
  )

#Pearson-Correlation
cor_mat <- cor(mat_top, method = "pearson", use = "pairwise.complete.obs")
stopifnot(all(colnames(cor_mat) == meta$sample))
annot <- data.frame(group = meta$group)
rownames(annot) <- meta$sample
labels <- setNames(meta$label, meta$sample)

ann_cols <- list(
  Group = c(
    "UT"              = "#4C78A8",
    "2H 4-ME"         = "#F58518",
    "22H 4-ME"        = "#E45756",
    "20H LPS"         = "#72B7B2",
    "Combined, Short" = "#54A24B",
    "Combined, Long"  = "#B279A2"
  )
)

#Heatmap
dist_fun <- function(x) as.dist(1 - x)
cols   <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(200)
pheatmap(
  cor_mat,
  clustering_distance_rows = dist_fun(cor_mat),
  clustering_distance_cols = dist_fun(cor_mat),
  clustering_method = "complete",
  display_numbers = TRUE,
  number_format = "%.2f",
  number_color = "grey20",
  color = cols,
  main  = "Pearson correlation between samples",
  annotation_col = annot,
  annotation_colors = ann_cols,
  labels_col = labels[colnames(cor_mat)],
  labels_row = labels[rownames(cor_mat)],
  fontsize_number = 8,
  border_color = NA,
  legend = TRUE,
  annotation_legend = TRUE
)

#---- KEGG ----
library(clusterProfiler)
library(enrichplot)

#RNA-seq 2H 4ME 
rna_up1 <- df1 %>%
  filter(padj < 0.05) %>%
  pull(gene)
rna_entrez1 <- bitr(rna_up1, fromType = "SYMBOL",
                    toType = "ENTREZID",
                    OrgDb = org.Hs.eg.db)
rna_kegg1 <- enrichKEGG(gene = rna_entrez1$ENTREZID,
                        organism = "hsa",  
                        pvalueCutoff = 0.05)

p1 <- suppressWarnings(barplot(rna_kegg1, showCategory = 10, title = "2H 4-ME"))

#RNA-seq 22H 4ME
rna_up2 <- df2 %>%
  filter(padj < 0.05) %>%
  pull(gene)
rna_entrez2 <- bitr(rna_up2, fromType = "SYMBOL",
                    toType = "ENTREZID",
                    OrgDb = org.Hs.eg.db)
rna_kegg2 <- enrichKEGG(gene = rna_entrez2$ENTREZID,
                        organism = "hsa",  
                        pvalueCutoff = 0.05)

p2 <- suppressWarnings(barplot(rna_kegg2, showCategory = 10, title = "22H 4-ME"))

#20H LPS
rna_up3 <- df3 %>%
  filter(padj < 0.05) %>%
  pull(gene)
rna_entrez3 <- bitr(rna_up3, fromType = "SYMBOL",
                    toType = "ENTREZID",
                    OrgDb = org.Hs.eg.db)
rna_kegg3 <- enrichKEGG(gene = rna_entrez3$ENTREZID,
                        organism = "hsa",  
                        pvalueCutoff = 0.05)

p3 <- suppressWarnings(barplot(rna_kegg3, showCategory = 10, title = "20H LPS"))

#Combined, Short 
rna_up4 <- df4 %>%
  filter(padj < 0.05) %>%
  pull(gene)
rna_entrez4 <- bitr(rna_up4, fromType = "SYMBOL",
                    toType = "ENTREZID",
                    OrgDb = org.Hs.eg.db)
rna_kegg4 <- enrichKEGG(gene = rna_entrez4$ENTREZID,
                        organism = "hsa",  
                        pvalueCutoff = 0.05)

p4 <- suppressWarnings(barplot(rna_kegg4, showCategory = 10, title = "Combined, Short"))

#Combined, Long
rna_up5 <- df5 %>%
  filter(padj < 0.05) %>%
  pull(gene)
rna_entrez5 <- bitr(rna_up5, fromType = "SYMBOL",
                    toType = "ENTREZID",
                    OrgDb = org.Hs.eg.db)
rna_kegg5 <- enrichKEGG(gene = rna_entrez5$ENTREZID,
                        organism = "hsa",  
                        pvalueCutoff = 0.05)

p5 <- suppressWarnings(barplot(rna_kegg5, showCategory = 10, title = "Combined, Long"))

#Combined, short - LC
rna_up6 <- df4_lc %>%
  filter(padj < 0.05 & log2FoldChange > 1) %>%
  pull(gene)
rna_entrez6 <- bitr(rna_up6, fromType = "SYMBOL",
                    toType = "ENTREZID",
                    OrgDb = org.Hs.eg.db)
rna_kegg6 <- enrichKEGG(gene = rna_entrez6$ENTREZID,
                        organism = "hsa",  
                        pvalueCutoff = 0.05)

p6 <- suppressWarnings(barplot(rna_kegg6, showCategory = 10, title = "Combined, Short - Log2 Relative to LPS"))

#Combined, Long - LC
rna_up7 <- df5_lc %>%
  filter(padj < 0.05 & log2FoldChange > 1) %>%
  pull(gene)
rna_entrez7 <- bitr(rna_up7, fromType = "SYMBOL",
                    toType = "ENTREZID",
                    OrgDb = org.Hs.eg.db)
rna_kegg7 <- enrichKEGG(gene = rna_entrez7$ENTREZID,
                        organism = "hsa",  
                        pvalueCutoff = 0.05)

p7 <- suppressWarnings(barplot(rna_kegg7, showCategory = 10, title = "Combined, Long - Log2 Relative to LPS"))

library(patchwork)

p1+p2+p3+p4+p5

p6+p7

p7

#---- GOBP ----


#2H 4ME
rna_up_genes1 <- df1 %>%
  filter(padj < 0.05) %>%
  pull(gene) %>%
  unique()
rna_gene_df1 <- bitr(rna_up_genes1, fromType = "SYMBOL",
                     toType = "ENTREZID",
                     OrgDb = org.Hs.eg.db)
rna_go1 <- enrichGO(gene         = rna_gene_df1$ENTREZID,
                    OrgDb        = org.Hs.eg.db,
                    keyType      = "ENTREZID",
                    ont          = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05)
p8 <- suppressWarnings(barplot(rna_go1, showCategory = 10, title = "2H 4-ME"))

#22H 4ME
rna_up_genes2 <- df2 %>%
  filter(padj < 0.05) %>%
  pull(gene) %>%
  unique()
rna_gene_df2 <- bitr(rna_up_genes2, fromType = "SYMBOL",
                     toType = "ENTREZID",
                     OrgDb = org.Hs.eg.db)
rna_go2 <- enrichGO(gene         = rna_gene_df2$ENTREZID,
                    OrgDb        = org.Hs.eg.db,
                    keyType      = "ENTREZID",
                    ont          = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05)
p9 <- suppressWarnings(barplot(rna_go2, showCategory = 10, title = "22H 4-ME"))

#20H LPS
rna_up_genes3 <- df3 %>%
  filter(padj < 0.05) %>%
  pull(gene) %>%
  unique()
rna_gene_df3 <- bitr(rna_up_genes3, fromType = "SYMBOL",
                     toType = "ENTREZID",
                     OrgDb = org.Hs.eg.db)
rna_go3 <- enrichGO(gene         = rna_gene_df3$ENTREZID,
                    OrgDb        = org.Hs.eg.db,
                    keyType      = "ENTREZID",
                    ont          = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05)
p10 <- suppressWarnings(barplot(rna_go3, showCategory = 10, title = "20H LPS"))

#Combined, Short 
rna_up_genes4 <- df4 %>%
  filter(padj < 0.05) %>%
  pull(gene) %>%
  unique()
rna_gene_df4 <- bitr(rna_up_genes4, fromType = "SYMBOL",
                     toType = "ENTREZID",
                     OrgDb = org.Hs.eg.db)
rna_go4 <- enrichGO(gene         = rna_gene_df4$ENTREZID,
                    OrgDb        = org.Hs.eg.db,
                    keyType      = "ENTREZID",
                    ont          = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05)
p11 <- suppressWarnings(barplot(rna_go4, showCategory = 10, title = "Combined, Short"))

#Combined, Long
rna_up_genes5 <- df5 %>%
  filter(padj < 0.05) %>%
  pull(gene) %>%
  unique()
rna_gene_df5 <- bitr(rna_up_genes5, fromType = "SYMBOL",
                     toType = "ENTREZID",
                     OrgDb = org.Hs.eg.db)
rna_go5 <- enrichGO(gene         = rna_gene_df5$ENTREZID,
                    OrgDb        = org.Hs.eg.db,
                    keyType      = "ENTREZID",
                    ont          = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05)
p12 <- suppressWarnings(barplot(rna_go5, showCategory = 10, title = "Combined, Long"))

#Combined, Short, LOG2 changed
rna_up_genes6 <- df4_lc %>%
  filter(padj < 0.05) %>%
  pull(gene) %>%
  unique()
rna_gene_df6 <- bitr(rna_up_genes6, fromType = "SYMBOL",
                     toType = "ENTREZID",
                     OrgDb = org.Hs.eg.db)
rna_go6 <- enrichGO(gene         = rna_gene_df6$ENTREZID,
                    OrgDb        = org.Hs.eg.db,
                    keyType      = "ENTREZID",
                    ont          = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05)
p13 <- suppressWarnings(barplot(rna_go6, showCategory = 10, title = "Combined, Short - Log2 Relative to LPS"))

#Combined, Long, Log2 changed 
rna_up_genes7 <- df5_lc %>%
  filter(padj < 0.05) %>%
  pull(gene) %>%
  unique()
rna_gene_df7 <- bitr(rna_up_genes7, fromType = "SYMBOL",
                     toType = "ENTREZID",
                     OrgDb = org.Hs.eg.db)
rna_go7 <- enrichGO(gene         = rna_gene_df6$ENTREZID,
                    OrgDb        = org.Hs.eg.db,
                    keyType      = "ENTREZID",
                    ont          = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05)
p14 <- suppressWarnings(barplot(rna_go7, showCategory = 10, title = "Combined, Long - Log2 Relative to LPS"))

#Plot all
p8+p9+p10+p11+p12

#Plot LC
p13+p14

