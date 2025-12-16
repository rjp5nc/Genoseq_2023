library(SeqArray)
library(pegas)
library(ape)

gds_path <- "/scratch/rjp5nc/UK2022_2024/allsites_mito/usdobtusa_mito_allsites.annot2.ALL.gds"

g <- seqOpen(gds_path)

# pull data
pos     <- seqGetData(g, "position")
samples <- seqGetData(g, "sample.id")
gt      <- seqGetData(g, "genotype")     # [ploidy, sample, variant]
gt1     <- gt[1, , ]                     # haploid mito: [sample, variant]

# recode SeqArray coding -> 0/1/NA (0=ref, 1=alt)
x <- gt1
x[x == 0] <- NA                          # missing
x[!is.na(x)] <- as.integer(x[!is.na(x)] != 1)  # 1 (ref) -> 0, >=2 (alt) -> 1

# keep variable sites only (recommended)
keep_non_na <- colSums(!is.na(x)) > 0
var_sites <- vapply(seq_len(ncol(x)), function(j) {
  v <- x[, j]; v <- v[!is.na(v)]
  length(unique(v)) >= 2
}, logical(1))
keep <- keep_non_na & var_sites

x <- x[, keep, drop = FALSE]
pos_use <- pos[keep]

seqClose(g)

# build a DNAbin alignment for pegas:
# 0 -> "A" (placeholder), 1 -> "G" (placeholder), NA -> "n"
# (For haplotypes, the actual nucleotide identity is not required; consistency is.)
mat <- matrix("n", nrow = nrow(x), ncol = ncol(x),
              dimnames = list(samples, paste0("p", pos_use)))

mat[x == 0] <- "a"
mat[x == 1] <- "g"

dna <- as.DNAbin(mat)

# pegas haplotypes
haps <- haplotype(dna)          # unique haplotypes
hap_assign <- with(stack(setNames(attr(haps, "index"), rownames(haps))), data.frame(sample = values, hap = ind))

# haplotype network
net <- haploNet(haps)

# plot network (sizes proportional to hap counts)
hap_freq <- table(attr(haps, "index"))
sizes <- as.numeric(hap_freq)

png("/scratch/rjp5nc/UK2022_2024/allsites_mito/mito_haplotype_network_pegas.png",
    width = 2000, height = 2000, res = 300)

plot(net, size = sizes, scale.ratio = 0.6, show.mutation = 1)

dev.off()
