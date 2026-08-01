#!/usr/bin/env Rscript
# Generate synthetic salmon.merged.{gene,transcript}_counts.tsv + tx2gene.tsv that
# match nf-core's format and use REAL ENSMUSG ids (so annotation is exercised too).
# Used only for development before the real pipeline output exists.
#   export TUY_DEV_DATA=/hpf/.../rnaseq_TUY35595/dev_data
#   Rscript analysis/make_synthetic_data.R
source("/hpf/projects/msalter/sghazis/rnaseq_TUY35595/analysis/_common.R")
suppressPackageStartupMessages(library(org.Mm.eg.db))

devdir <- Sys.getenv("TUY_DEV_DATA"); stopifnot(nzchar(devdir))
dir.create(devdir, showWarnings = FALSE, recursive = TRUE)
set.seed(42)

meta <- load_metadata(); samples <- rownames(meta)
NG <- 1500
genes <- head(keys(org.Mm.eg.db, keytype = "ENSEMBL"), 4000)
genes <- sample(genes, NG)
gsym  <- mapIds(org.Mm.eg.db, genes, "SYMBOL", "ENSEMBL")

# baseline + condition/sex effects on a subset to create real signal
base <- rnbinom(NG, mu = 200, size = 2) + 5
mk <- function(s) {
  cond <- meta[s, "condition"]; sx <- meta[s, "sex"]
  mu <- base
  up <- sample(NG, 120); mu[up] <- mu[up] * if (cond == "NP") 3 else 1
  dn <- sample(NG, 80);  mu[dn] <- mu[dn] * if (sx == "M") 1.8 else 1
  rnbinom(NG, mu = mu, size = 3)
}
gmat <- sapply(samples, mk); rownames(gmat) <- genes

gene_df <- data.frame(gene_id = genes, gene_name = gsym, gmat, check.names = FALSE)
write.table(gene_df, file.path(devdir, "salmon.merged.gene_counts.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# transcripts: 1-3 per gene, unique tx ids, counts split from the gene
ntx <- sample(1:3, NG, replace = TRUE)
tx  <- data.frame(tx = sprintf("ENSMUST%011d.1", seq_len(sum(ntx))),
                  gene_id = rep(genes, times = ntx),
                  stringsAsFactors = FALSE)
frac <- runif(nrow(tx), .2, 1)
txmat <- t(sapply(seq_len(nrow(tx)), function(j)
  rpois(length(samples), gmat[match(tx$gene_id[j], genes), ] * frac[j] / 2)))
colnames(txmat) <- samples
tx_df <- data.frame(tx = tx$tx, gene_id = tx$gene_id, txmat, check.names = FALSE)
write.table(tx_df, file.path(devdir, "salmon.merged.transcript_counts.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

t2g <- data.frame(tx$tx, tx$gene_id, gsym[match(tx$gene_id, genes)])
write.table(t2g, file.path(devdir, "tx2gene.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

cat("Synthetic dev data written to", devdir, "\n",
    NG, "genes,", nrow(tx), "transcripts,", length(samples), "samples\n")
