#!/usr/bin/env Rscript
###############################################################################
# nanquplus_from_tob.R
#
# Build quartet table (parallel) + p-values from gene trees, then resolve level-1
# networks using a user-supplied ToB (tree-of-blobs) file.
#
# Usage:
#   Rscript nanquplus_from_tob.R \
#     <gene_trees.tre> <tob.newick> <out_prefix> \
#     [num_cores] [alpha] [beta] [t3_model]
#
# Examples:
#   Rscript nanquplus_from_tob.R gene_trees.tre my_tob.nwk out/run1 16 0.05 0.95 T3
#
# Outputs:
#   <out_prefix>.ptable.rds            # p-value-augmented quartet table (pT)
#   <out_prefix>.tob_labeled.nwk       # ToB with internal nodes labeled
#   <out_prefix>.level1_networks.rds   # list of inferred level-1 networks (MSCquartets objects)
#   <out_prefix>.level1_networks.txt   # printed summary of networks
###############################################################################

suppressPackageStartupMessages({
  library(ape)
  library(MSCquartets)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  cat("Usage:\n",
      "  Rscript nanquplus_from_tob.R <gene_trees.tre> <tob.newick> <out_prefix> [num_cores] [alpha] [beta] [t3_model]\n\n",
      "Example:\n",
      "  Rscript nanquplus_from_tob.R gene_trees.tre my_tob.nwk out/run1 16 0.05 0.95 T3\n",
      sep = "")
  quit(status = 2)
}

gene_trees_file <- args[1]
tob_file        <- args[2]
out_prefix      <- args[3]

num_cores <- if (length(args) >= 4) as.integer(args[4]) else 4L
alpha     <- if (length(args) >= 5) as.numeric(args[5]) else 0.05
beta      <- if (length(args) >= 6) as.numeric(args[6]) else 0.95
t3_model  <- if (length(args) >= 7) as.character(args[7]) else "T3"

# ---- basic checks ----
if (!file.exists(gene_trees_file)) stop("Gene trees file not found: ", gene_trees_file)
if (!file.exists(tob_file))        stop("ToB file not found: ", tob_file)

if (is.na(num_cores) || num_cores < 1) stop("num_cores must be a positive integer.")
if (is.na(alpha) || alpha <= 0 || alpha >= 1) stop("alpha must be in (0,1).")
if (is.na(beta)  || beta  <= 0 || beta  >= 1) stop("beta must be in (0,1).")

message("== Inputs ==")
message("Gene trees: ", gene_trees_file)
message("ToB:        ", tob_file)
message("Out prefix: ", out_prefix)
message("Cores:      ", num_cores)
message("alpha:      ", alpha)
message("beta:       ", beta)
message("T3 model:   ", t3_model)

# ---- 1) Read gene trees ----
message("\n[1/6] Reading gene trees...")
gtrees <- read.tree(gene_trees_file)
if (is.null(gtrees) || length(gtrees) < 1) stop("No trees read from: ", gene_trees_file)
class(gtrees) <- "multiPhylo"

# ---- 2) Build quartet counts in parallel ----
message("[2/6] Computing quartet count table in parallel...")
if (num_cores >= 2) {
  QT <- quartetTableParallel(gtrees, numCores = num_cores)
} else {
  QT <- quartetTable(gtrees)
}


# ---- 3) Resolve quartets (required before tests) ----
message("[3/6] Resolving quartet table...")
RQT <- quartetTableResolved(QT)

# ---- 4) Add p-values (tree test + star test) ----
message("[4/6] Computing p-values...")
pT <- quartetTreeTestInd(RQT, model = t3_model)
pT <- quartetStarTestInd(pT)

ptable_out <- paste0(out_prefix, ".ptable.rds")
saveRDS(pT, ptable_out)
message("Saved pT (p-value table) to: ", ptable_out)

# ---- 5) Read and label ToB ----
message("[5/6] Reading and labeling ToB internal nodes...")
# ToB <- read.tree(tob_file)
ToB <- read.tree(text = write.tree(read.tree(tob_file)))

ToB <- unroot(ToB)
ToB$root.edge <- NULL
ToB$node.label <- NULL

ToB$edge.length <- NULL
if (is.null(ToB) || length(ToB$tip.label) < 4) stop("ToB must have at least 4 tips.")
ToB <- unroot(ToB)

ToB <- labelIntNodes(ToB, plot = FALSE)

cat("is.rooted =", is.rooted(ToB), "\n")
cat("root.edge =", ToB$root.edge, "\n")

tob_labeled_out <- paste0(out_prefix, ".tob_labeled.nwk")
write.tree(ToB, file = tob_labeled_out)
message("Saved labeled ToB to: ", tob_labeled_out)

# ---- sanity: taxa match ----
tob_taxa <- sort(ToB$tip.label)
# pT can be different types depending on MSCquartets version; try to extract taxa set robustly
get_ptable_taxa <- function(pTobj) {
  if (is.matrix(pTobj) || is.data.frame(pTobj)) {
    # qcCF tables often have taxa indicator columns; try colnames that look like taxa
    cn <- colnames(pTobj)
    if (!is.null(cn)) return(sort(cn))
  }
  return(NULL)
}

# ---- 6) Resolve level-1 networks using your ToB ----
message("[6/6] Resolving level-1 network(s) from your ToB...")


rm(QT, RQT); gc()

# Disable plotting safely (no namespace hacking)
plot.phylo <- function(...) invisible(NULL)
plot <- function(...) invisible(NULL)

ensure_edge_lengths <- function(tr, len = 1) {
  if (is.null(tr$edge.length) || any(!is.finite(tr$edge.length))) {
    tr$edge.length <- rep(len, nrow(tr$edge))
  }
  tr
}

# if (is.rooted(ToB)) {
#   ToB <- unroot(ToB)
# }



ToB <- ensure_edge_lengths(ToB, len = 1)

nets <- resolveLevel1(ToB, pT, alpha = alpha, beta = beta, plot = FALSE)

nets_txt_out <- paste0(out_prefix, ".level1_networks.txt")
capture.output({
  cat("resolveLevel1 output object:\n")
  print(nets)
  cat("\n---\nSession info:\n")
  print(sessionInfo())
}, file = nets_txt_out)
message("Saved networks summary to: ", nets_txt_out)


final_nwk <- nets[[1]][[1]]
writeLines(final_nwk, paste0(out_prefix, ".level1_network.nwk"))
message("Saved final level-1 network to: ",
        paste0(out_prefix, ".level1_network.nwk"))


message("\nDone.")
