# Install and load required packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "ReactomePA", "DOSE"))
install.packages(c("tidyverse", "rWikiPathways"))

library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(rWikiPathways)
library(tidyverse)
library(DOSE)

# Input genes
genes <- c("FZD7", "ZNF248", "TLR2", "JAG1", "TCF7L2", "SULF2", "SLMAP", "LY9", 
           "DMXL2", "RBM47", "UBASH3B", "ATG7", "LILRB2", "AIF1", "NPL", 
           "GRK3", "GLUL", "GJB6")

# Convert gene symbols to ENTREZ IDs
gene_ids <- bitr(genes, fromType = "SYMBOL",
                 toType = "ENTREZID",
                 OrgDb = org.Hs.eg.db)

# KEGG pathway analysis
kegg_enrich <- enrichKEGG(
  gene = gene_ids$ENTREZID,
  organism = "hsa",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

# Create and save KEGG plot
kegg_data <- data.frame(
  Pathway = head(kegg_enrich@result$Description, 15),
  GeneRatio = sapply(head(kegg_enrich@result$GeneRatio, 15), 
                     function(x) eval(parse(text = x))),
  Pvalue = head(kegg_enrich@result$pvalue, 15)
)

kegg_plot <- ggplot(kegg_data, 
                    aes(x = GeneRatio, y = reorder(Pathway, GeneRatio))) +
  geom_point(aes(size = -log10(Pvalue)), color = "#FF6B6B") +
  scale_size_continuous(range = c(3, 10)) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "right"
  ) +
  labs(
    title = "KEGG Pathway Enrichment Analysis",
    x = "Gene Ratio",
    y = "Pathway",
    size = "-log10(p-value)"
  )

print(kegg_plot)
ggsave("kegg_enrichment_plot.pdf", kegg_plot, width = 12, height = 8)

# WikiPathways analysis
wiki_enrich <- enrichWP(
  gene = gene_ids$ENTREZID,
  organism = "Homo sapiens",
  pvalueCutoff = 0.05
)

# Create and save WikiPathways plot
wiki_data <- data.frame(
  Pathway = head(wiki_enrich@result$Description, 15),
  GeneRatio = sapply(head(wiki_enrich@result$GeneRatio, 15), 
                     function(x) eval(parse(text = x))),
  Pvalue = head(wiki_enrich@result$pvalue, 15)
)

wiki_plot <- ggplot(wiki_data, 
                    aes(x = GeneRatio, y = reorder(Pathway, GeneRatio))) +
  geom_point(aes(size = -log10(Pvalue)), color = "#4ECDC4") +
  scale_size_continuous(range = c(3, 10)) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "right"
  ) +
  labs(
    title = "WikiPathways Enrichment Analysis",
    x = "Gene Ratio",
    y = "Pathway",
    size = "-log10(p-value)"
  )

print(wiki_plot)
ggsave("wiki_enrichment_plot.pdf", wiki_plot, width = 12, height = 8)

# Reactome pathway analysis
reactome_enrich <- enrichPathway(
  gene = gene_ids$ENTREZID,
  organism = "human",
  pvalueCutoff = 0.05
)

# Create and save Reactome plot
reactome_data <- data.frame(
  Pathway = head(reactome_enrich@result$Description, 15),
  GeneRatio = sapply(head(reactome_enrich@result$GeneRatio, 15), 
                     function(x) eval(parse(text = x))),
  Pvalue = head(reactome_enrich@result$pvalue, 15)
)

reactome_plot <- ggplot(reactome_data, 
                        aes(x = GeneRatio, y = reorder(Pathway, GeneRatio))) +
  geom_point(aes(size = -log10(Pvalue)), color = "#45B7D1") +
  scale_size_continuous(range = c(3, 10)) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "right"
  ) +
  labs(
    title = "Reactome Pathway Enrichment Analysis",
    x = "Gene Ratio",
    y = "Pathway",
    size = "-log10(p-value)"
  )

print(reactome_plot)
ggsave("reactome_enrichment_plot.pdf", reactome_plot, width = 12, height = 8)

# Print status message
cat("\nAll plots have been created and saved:\n",
    "1. kegg_enrichment_plot.pdf\n",
    "2. wiki_enrichment_plot.pdf\n",
    "3. reactome_enrichment_plot.pdf\n")