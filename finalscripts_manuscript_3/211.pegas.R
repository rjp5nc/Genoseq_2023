library(SeqArray)
library(pegas)
library(ape)
library(data.table)

gds_path <- "/scratch/rjp5nc/UK2022_2024/allsites_mito/usdobtusa_mito_allsites.annot2.ALL.gds"
out_dir  <- "/scratch/rjp5nc/UK2022_2024/allsites_mito/pegas_outputs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## -----------------------
## Open GDS and extract data
## -----------------------
g <- seqOpen(gds_path)

pos     <- seqGetData(g, "position")
samples <- seqGetData(g, "sample.id")
gt      <- seqGetData(g, "genotype")   # [ploidy, sample, variant]

# haploid mito
gt1 <- gt[1, , ]   # [sample, variant]

# recode to 0/1/NA
x <- gt1
x[x == 0] <- NA
x[!is.na(x)] <- as.integer(x[!is.na(x)] != 1)

## -----------------------
## Keep variable sites only
## -----------------------
keep_non_na <- colSums(!is.na(x)) > 0
var_sites <- vapply(seq_len(ncol(x)), function(j) {
  v <- x[, j]
  v <- v[!is.na(v)]
  length(unique(v)) >= 2
}, logical(1))

keep <- keep_non_na & var_sites

x <- x[, keep, drop = FALSE]
pos_use <- pos[keep]

seqClose(g)

## -----------------------
## Build DNAbin for pegas
## -----------------------
# placeholder bases: ref="a", alt="g", missing="n"
dna_mat <- matrix(
  "n",
  nrow = nrow(x),
  ncol = ncol(x),
  dimnames = list(samples, paste0("pos_", pos_use))
)

dna_mat[x == 0] <- "a"
dna_mat[x == 1] <- "g"

dna <- as.DNAbin(dna_mat)

## -----------------------
## pegas haplotypes + network
## -----------------------
haps <- haplotype(dna)
net  <- haploNet(haps)

## -----------------------
## Save reusable R objects
## -----------------------
saveRDS(dna,  file.path(out_dir, "mito_dna_dnabin.rds"))
saveRDS(haps, file.path(out_dir, "mito_haplotypes.rds"))
saveRDS(net,  file.path(out_dir, "mito_haplotype_network.rds"))
saveRDS(pos_use, file.path(out_dir, "mito_variable_positions.rds"))

## -----------------------
## Save tabular summaries
## -----------------------

# sample → haplotype assignment
hap_assign <- data.table(
  sample = unlist(attr(haps, "index")),
  hap_id = rep(seq_along(attr(haps, "index")),
               lengths(attr(haps, "index")))
)

fwrite(hap_assign,
       file.path(out_dir, "sample_to_haplotype.csv"))

# haplotype frequencies
hap_freq <- hap_assign[, .N, by = hap_id][order(-N)]
fwrite(hap_freq,
       file.path(out_dir, "haplotype_frequencies.csv"))

message("Saved pegas haplotype outputs to: ", out_dir)

library(pegas)

out_dir <- "/scratch/rjp5nc/UK2022_2024/allsites_mito/pegas_outputs"

haps <- readRDS(file.path(out_dir, "mito_haplotypes.rds"))
net  <- readRDS(file.path(out_dir, "mito_haplotype_network.rds"))

hap_freq <- read.csv(file.path(out_dir, "haplotype_frequencies.csv"))
sizes <- hap_freq$N

png("/scratch/rjp5nc/UK2022_2024/allsites_mito/mito_haplotype_network_pegas.png",
    width = 2000, height = 2000, res = 300)

plot(net,
     size = sizes,
     scale.ratio = 0.6,
     show.mutation = 1)

dev.off()
