library(SeqArray)
library(pegas)
library(ape)
library(igraph)

## ----------------------------
## Inputs
## ----------------------------
gds_path <- "/scratch/rjp5nc/UK2022_2024/allsites_mito/usdobtusa_mito_allsites.annot2.ALL.gds"
out_dir  <- "/scratch/rjp5nc/UK2022_2024/allsites_mito/pegas_outputs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

out_png <- file.path(out_dir, "minimum_spanning_haplotype_network.png")
out_rds <- file.path(out_dir, "minimum_spanning_haplotype_network_objects.rds")

## ----------------------------
## 1) Read GDS (haploid mito) -> 0/1/NA matrix
## ----------------------------
g <- seqOpen(gds_path)

samples <- seqGetData(g, "sample.id")
pos     <- seqGetData(g, "position")
gt      <- seqGetData(g, "genotype")     # [ploidy, sample, variant]
gt1     <- gt[1, , ]                     # haploid: [sample, variant]

x <- gt1
x[x == 0] <- NA                          # missing
x[!is.na(x)] <- as.integer(x[!is.na(x)] != 1)  # ref->0, alt->1

seqClose(g)

## ----------------------------
## 2) Keep variable sites only
## ----------------------------
keep_non_na <- colSums(!is.na(x)) > 0
var_sites <- vapply(seq_len(ncol(x)), function(j) {
  v <- x[, j]; v <- v[!is.na(v)]
  length(unique(v)) >= 2
}, logical(1))

keep <- keep_non_na & var_sites
x <- x[, keep, drop = FALSE]
pos_use <- pos[keep]

## ----------------------------
## 3) Build DNAbin (placeholder bases; missing = n)
## ----------------------------
dna_mat <- matrix(
  "n",
  nrow = nrow(x),
  ncol = ncol(x),
  dimnames = list(samples, paste0("pos_", pos_use))
)
dna_mat[x == 0] <- "a"
dna_mat[x == 1] <- "g"

dna <- as.DNAbin(dna_mat)

## ----------------------------
## 4) Haplotypes + distances + MST
## ----------------------------
haps <- haplotype(dna)

# haplotype frequencies (robust; avoids igraph masking)
hap_freq <- lengths(attr(haps, "index"))
names(hap_freq) <- rownames(haps)   # label nodes by hap rownames (I, II, III, ...)

# use haplotypes directly as DNAbin via matrix roundtrip (stable)
hap_dna <- as.DNAbin(as.matrix(haps))

D <- dist.dna(
  hap_dna,
  model = "N",
  pairwise.deletion = TRUE,
  as.matrix = TRUE
)

mst_adj <- ape::mst(D)

g_net <- graph_from_adjacency_matrix(mst_adj, mode = "undirected", diag = FALSE)

## Edge weights = mutational steps; labels are dashes
E(g_net)$weight <- apply(ends(g_net, E(g_net)), 1, function(e) D[e[1], e[2]])
E(g_net)$label  <- vapply(E(g_net)$weight, function(w) paste(rep("-", w), collapse = ""), character(1))

## Node sizes proportional to haplotype frequency
V(g_net)$size <- 6 + 2 * as.numeric(hap_freq[V(g_net)$name])

## ----------------------------
## 5) Plot
## ----------------------------
png(out_png, width = 2000, height = 1600, res = 300)

plot(
  g_net,
  layout = layout_with_kk(g_net),
  vertex.color = "white",
  vertex.frame.color = "grey30",
  vertex.label = V(g_net)$name,
  edge.color = "grey40",
  edge.width = 1.2,
  edge.label = E(g_net)$label,
  edge.label.cex = 0.9
)

dev.off()

## ----------------------------
## 6) Save reusable objects
## ----------------------------
saveRDS(
  list(
    dna = dna,
    haps = haps,
    hap_freq = hap_freq,
    hap_dna = hap_dna,
    D = D,
    mst_adj = mst_adj,
    g_net = g_net,
    pos_use = pos_use,
    samples = samples
  ),
  out_rds
)

cat("Saved plot: ", out_png, "\n", sep = "")
cat("Saved objects: ", out_rds, "\n", sep = "")
