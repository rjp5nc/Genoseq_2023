#ijob -A berglandlab -c10 -p standard --mem=40G

#module load gcc/11.4.0  openmpi/4.1.4 icu R/4.3.1

#R

#install.packages(c("usethis","pkgdown","gh","rcmdcheck"), repos="https://cloud.r-project.org")


#install.packages("devtools")
#library(devtools)

#devtools::install_github("pievos101/PopGenome")

library(PopGenome)

#bcftools view -h /scratch/rjp5nc/UK2022_2024/daphnia_phylo/trimmed_10bp_repeatmasked_vcf/trimmed10bp_masked_usobtusa.bgz.vcf.gz | \
#grep "^##contig" | cut -d'=' -f3 | cut -d',' -f1 > contigs_list.txt

# 1. Load the VCF
vcf_file <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/trimmed_10bp_repeatmasked_vcf/trimmed10bp_masked_usobtusa.bgz.vcf.gz"
contigs_file <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/trimmed_10bp_repeatmasked_vcf/contigs_list.txt"

# Read into R as a vector
contigs <- readLines(contigs_file)
print(contigs)

# Loop over all contigs and read each one
genomes <- lapply(contigs, function(ctg) {
  readVCF(
    vcf_file,
    numcols = 10000,
    tid = ctg,
    frompos = 1,
    topos = 999999999,
    include.unknown = TRUE
  )
})

# Merge them into one PopGenome object
genome <- do.call(c, genomes)

genome_win <- sliding.window.transform(
  genome,
  width = 50000,   # window size in bp
  jump = 50000     # step size = window size → non-overlapping
)

pops <- read.table("/scratch/rjp5nc/UK2022_2024/mito_vcf/pops_fixed_filtered.txt",
                   header = FALSE, stringsAsFactors = FALSE)
colnames(pops) <- c("sample", "population")

to_remove <- pops$population == "PBO66"

# Keep only non-PBO66 samples and their population labels
pops$sample <- pops$sample[!to_remove]
pops$population <- pops$population[!to_remove]

# Split into list of populations
pop_list <- split(pops$sample, pops$population)

# Set populations
genome_win <- set.populations(genome_win, pop_list)

# Set Gilmer as outgroup
genome_win <- set.outgroup(genome_win, pop_list$Gilmer, diploid = FALSE)

# Calculate π (within-pop diversity)
genome_win <- diversity.stats(genome_win)
pi_vals <- genome_win@nuc.diversity.within

# Calculate Dxy (between populations)
genome_win <- diversity.stats.between(genome_win)
dxy_matrix <- genome_win@nuc.diversity.between

# Show results
pi_vals
dxy_matrix

colnames(pi_vals) <- names(pop_list)

# Rename Dxy results
colnames(dxy_matrix) <- combn(names(pop_list), 2, FUN = paste, collapse = "/")



genome_win <- neutrality.stats(genome_win, detail = TRUE)  # This computes theta_Watterson
theta_values <- genome_win@theta_Watterson

colnames(theta_values) <- names(pop_list)

theta_values


tajima.D <- genome_win@Tajima.D 
colnames(tajima.D) <- names(pop_list)

tajima.D

write.csv(pi_vals, file="/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/pi_vals.csv", row.names=FALSE)
write.csv(dxy_matrix, file="/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/dxy_matrix.csv", row.names=FALSE)
write.csv(theta_values, file="/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/theta_values.csv", row.names=FALSE)
write.csv(tajima.D, file="/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/tajimaD.csv", row.names=FALSE)
