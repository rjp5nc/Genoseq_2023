
library(stringr)
library(tidyr)
library(tibble)
library(dplyr)
library(ggplot2)
library(data.table)
library(ggtree)
library(ape)
library(dplyr)

locationandclone <- read.csv("C:/Users/rjpor/Downloads/2022_2024seqmetadata20250811.csv")
multiqc <- read.csv("C:/Users/rjpor/Downloads/multiqc_fastqcfull.csv")
bamout <- read.csv("C:/Users/rjpor/Downloads/bamout.csv")
multiqcdup <- read.csv("C:/Users/rjpor/Downloads/multiqc_fastqc2.csv")
mappedreads <- read.csv("C:/Users/rjpor/Downloads/unmappedinsortedbam.csv")
#need to open mapped reads and save for proper format


mappedreads$BAM_File <- gsub(".sorted.bam", "", mappedreads$BAM_File)

bamout2 <- subset(bamout, Average_Depth >= 2)

bamout$Filename <- gsub(".sort.dedup.bam", "", bamout$Filename)

bamout2$Filename <- gsub(".sort.dedup.bam", "", bamout2$Filename)
#write.csv(bamout2, "/Users/rjpor/Downloads/bamout2.csv")


plates <- multiqc[c('plate', 'well', "else")] <- str_split_fixed(multiqc$Filename, '_', 3)
plates <- as.data.frame(plates)
colnames(plates) <- c("Plate","well", "other")

Comb <- cbind(multiqc, plates)
Combshort <- Comb[,c(1,5,9,24,25)]

Combshort$deduptotal <- Combshort$Total.Sequences*Combshort$total_deduplicated_percentage/100
merged_data <- aggregate(deduptotal ~ plate + well, data = Combshort, sum)
merged_data2 <- aggregate(Total.Sequences ~ plate + well, data = Combshort, sum)

setDT(merged_data)
setDT(merged_data2)
setkey(merged_data, plate, well)
setkey(merged_data2,  plate, well)
fordups <- merge(merged_data, merged_data2, all.x=TRUE, all.y=TRUE)

fordups$plate_well <- paste(fordups$plate, "_", fordups$well)
fordups$plate_well <- gsub(" ", "", fordups$plate_well)







merged_data <- merged_data %>%
  separate(well, into = c("Row", "Column"), sep = "(?<=[A-Za-z])", remove = FALSE)
# 

# merged_data <- merged_data %>%
#   mutate(
#     RowGroup = case_when(
#       Row %in% c("A", "B") ~ "Set1",
#       Row %in% c("C", "D") ~ "Set2",
#       Row %in% c("E", "F") ~ "Set3",
#       Row %in% c("G", "H") ~ "Set4"
#     )
#   )

merged_data <- merged_data %>%
  mutate(
    RowGroup = case_when(
      Row %in% c("A") ~ "Set1",
      Row %in% c("B") ~ "Set2",
      Row %in% c("C") ~ "Set3",
      Row %in% c("D") ~ "Set4",
      Row %in% c("E") ~ "Set5",
      Row %in% c("F") ~ "Set6",
      Row %in% c("G") ~ "Set7",
      Row %in% c("H") ~ "Set8"
    )
  )

merged_data <- merged_data %>%
  mutate(concentration = case_when(
    plate == "Rockpool1" ~ 27.5,
    plate == "Rockpool2" ~ 28.2,
    plate == "Rockpool3" ~ 29.4,
    plate == "Rockpool4" ~ 50,
    plate == "Gilmer5" ~ 39.2,
    plate == "RobertUK" ~ 20.8
  ))

merged_data$plateconc <- 1/(merged_data$concentration/50)
merged_data$depth <- merged_data$deduptotal*150/150000000
merged_data$plate_well <- paste(merged_data$plate,merged_data$well)

ggplot(subset(merged_data), aes(x = as.numeric(Column), y = RowGroup, fill = depth)) +
  geom_tile(color = "white") +
  facet_grid(RowGroup ~ plate, scales = "free_y", space = "free_y") +
  scale_fill_viridis_c(option = "plasma", name = "Value") +
  theme_minimal() + scale_x_continuous(breaks = seq(1, 12, by = 1)) +
  geom_text(aes(label = round(depth)), color = "white", size = 3,angle = 90) +
  theme(
    strip.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.spacing = unit(1, "lines")
  ) +
  labs(
    title = "96-Well Plate Visualization",
    x = "Column",
    y = "Row"
  )



merged_data$plate_well <- gsub(" ", "_", merged_data$plate_well)

locationandclone2 <- locationandclone
colnames(locationandclone2)[1] ="plate_well"



setDT(merged_data)
setDT(locationandclone2)
setkey(merged_data, plate_well)
setkey(locationandclone2,  plate_well)
merged2 <- merge(merged_data, locationandclone2, all.x=TRUE, all.y=TRUE)



pdf("C:\\Users\\rjpor\\Downloads\\clones.pdf", width=15, height=10)

wellxclone <- merged2[,c(1, 11, 12,14,15)]


write.csv(wellxclone, "C:\\Users\\rjpor\\Downloads\\wellxclone.csv")

ggplot(subset(merged2), aes(x = as.numeric(Column), y = RowGroup, fill = deduptotal)) +
  geom_tile(color = "white") +
  facet_grid(RowGroup ~ plate, scales = "free_y", space = "free_y") +
  scale_fill_viridis_c(option = "plasma", name = "Value") +
  theme_minimal() + scale_x_continuous(breaks = seq(1, 12, by = 1)) +
  geom_text(aes(label = clone), color = "white", size = 3,angle = 90) +
  theme(
    strip.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.spacing = unit(1, "lines")
  ) +
  labs(
    title = "96-Well Plate Visualization",
    x = "Column",
    y = "Row"
  )

dev.off()



pca_data <- read.csv("C:/Users/rjpor/Downloads/2022seqPCA.csv")



ggplot(subset(pca_data, Plate =="UK"|Plate=="us"), aes(x = PC1, y = PC2, color=Plate)) +
  geom_point() +
  labs(title = "PCA of Genomic Data", x = "PC1", y = "PC2") +
  theme_minimal()

ggplot(subset(pca_data, Plate =="UK"|Plate=="us"|Plate=="Gilmer"), aes(x = PC1, y = PC2, color=Plate)) +
  geom_point() +
  labs(title = "PCA of Genomic Data", x = "PC1", y = "PC2") +
  theme_minimal()


ggplot(pca_data, aes(x = PC1, y = PC2, color=Well)) +
  geom_point() + ylim(-.1,.1)+
  labs(title = "PCA of Genomic Data", x = "PC1", y = "PC2") +
  theme_minimal()+ guides(color="none")




load("C:/Users/rjpor/Downloads/nj_tree_gilmer.RData")



#tree$tip.label <- pca_data$Well[match(tree$tip.label, pca_data$Sample)]

tree_plot <- ggtree(tree, options(ignore.negative.edge=TRUE))



# Extract tree data
tree_data <- tree_plot$data
tree_data$label <- gsub("\\.sorted", "", tree_data$label)




pca_data
setDT(pca_data)
setDT(locationandclone)
setkey(pca_data, Well)
setkey(locationandclone, Well)
pca_datawithclone <- merge(pca_data, locationandclone, all.x=TRUE, all.y=TRUE)

write.csv(pca_datawithclone, "C:\\Users\\rjpor\\Downloads\\pca_withclone.csv")

# Perform the merge using left_join (same as merge but more readable)
merged_data2 <- left_join(tree_data, pca_datawithclone, by = c("label" = "Well"))

merged_data2 <- merged_data2 %>%
  separate(label, into = c("plate", "column"), sep = "_", remove = FALSE)


merged_data2$Plate

#tree_plot <- ggtree(merged_data2, layout="circular")
tree_plot <- ggtree(merged_data2, layout="rectangular")

pdf("C:\\Users\\rjpor\\Downloads\\treewgsgilmer.pdf", width=20, height=20)

tree_plot +
  geom_tiplab(aes(label = clone, color = date)) +  # Use new labels and color by group
  #scale_color_manual(values = c("us" = "blue", "UK" = "red", "Gilmer" = "green")) +  # Custom colors
  theme_tree2()  # Optional tree theme

dev.off()

aaaa <- tree_plot +
  geom_tiplab(aes(label = clone, color = date)) +  # Use new labels and color by group
  #scale_color_manual(values = c("us" = "blue", "UK" = "red", "Gilmer" = "green")) +  # Custom colors
  theme_tree2()  

ggsave("C:/Users/rjpor/Downloads/treecirdate.png", plot = aaaa, width = 20, height = 20, dpi = 300)

aaab <- tree_plot +
  geom_tiplab(aes(label = clone, color = accuratelocation)) +  # Use new labels and color by group
  #scale_color_manual(values = c("us" = "blue", "UK" = "red", "Gilmer" = "green")) +  # Custom colors
  theme_tree2()  

ggsave("C:/Users/rjpor/Downloads/treecir.png", plot = aaab, width = 20, height = 20, dpi = 300)















tree10000 <- tree

tree10000$tip.label <- gsub("\\.sort\\.dedup", "", tree10000$tip.label)


#tree$tip.label <- pca_data$Well[match(tree$tip.label, pca_data$Sample)]
tree_plot <- ggtree(tree10000, layout="circular",options(ignore.negative.edge=TRUE))

# Extract tree data
tree_data <- tree_plot$data



library(dplyr)

# Perform the merge using left_join (same as merge but more readable)
merged_data2 <- left_join(tree_data, pca_datawithclone, by = c("label" = "Well"))


# Update the tree plot data with merged data
tree_plot$data <- merged_data2


stuff <- as.data.table(tree_plot$data)

outgroup_sample <- "Sample4"  # Replace with the correct sample ID

#pdf("C:\\Users\\rjpor\\Downloads\\treemitoallfull.pdf", width=30, height=30)

tree_plot + 
  geom_tiplab(aes(label = clone)) +  # Use new labels and color by group
#  scale_color_manual(values = c("us" = "blue", "UK" = "red", "Gilmer" = "green")) +  # Custom colors
  theme_tree2()  # Optional tree theme

dev.off()

#sim_matrix <- read.csv("C:\\Users\\rjpor\\Downloads\\similarity_matrixGilmer.csv", row.names = 1)
sim_matrix <- read.csv("C:\\Users\\rjpor\\Downloads\\ibs_matrix_gilmer.csv", row.names = 1)
sim_matrix <- as.matrix(sim_matrix)  # Coerce it into a matrix

# Create a named vector for mapping
Sample<-pca_datawithclone$Well
d<-as.data.table(Sample)
#d$clone<- pca_datawithclone$clone
d$clone<- pca_datawithclone$date

name_map <-d

name_dict <- setNames(name_map$clone, name_map$Sample)



# Check the structure of sim_matrix
str(sim_matrix)

# Ensure it's a numeric matrix
sim_matrix <- as.matrix(sim_matrix)  # Coerce it into a matrix
mode(sim_matrix) <- "numeric"        # Ensure it's numeric

# Check for NAs and replace them (if necessary)

sim_matrix[is.na(sim_matrix)] <- 1  # Replace NAs with 0

# Check the dimensions to ensure it's square
if (nrow(sim_matrix) != ncol(sim_matrix)) {
  stop("The similarity matrix is not square!")
}



# Step 2: Replace row names in the similarity matrix
rownames(sim_matrix) <- name_dict[rownames(sim_matrix)]

# Step 3: Replace column names in the similarity matrix
colnames(sim_matrix) <- name_dict[colnames(sim_matrix)]

# subset_matrix <- sim_matrix["Elvis_03", "Elvis_03"]
# print(subset_matrix)

row_indices <- which(rownames(sim_matrix) == "D8.4A")
col_indices <- which(colnames(sim_matrix) == "Islands_02")

# Remove all occurrences of "B"
subset_matrix <- sim_matrix[row_indices, col_indices]
print(subset_matrix)


pdf("C:\\Users\\rjpor\\Downloads\\heatmapusobtusagilmer.pdf", width=30, height=30)


# Plot the heatmap
heatmap(sim_matrix, 
        main = "Genotype Similarity Heatmap", 
        col = colorRampPalette(c("blue", "red"))(50), 
        scale = "none", 
        margins = c(8, 8))

dev.off()








depth_data <- data.frame(Sample = names(avg_depth_per_sample), 
                         Depth = gsub("\\.sort\\.dedup$", "", avg_depth_per_sample), 
                         stringsAsFactors = FALSE)
depth_data$Sample <- str_remove(depth_data$Sample, ".sort.dedup")
depth_data$Depth <- as.numeric(depth_data$Depth)

p <- ggplot(bamout, aes(x = Average_Depth)) +
  geom_histogram(binwidth = 1, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "Histogram of Read Depth", x = "Average Read Depth", y = "Frequency") +
  theme_bw()

# Save plot as PNG
ggsave("C:/Users/rjpor/Downloads/avg_depth_histogramwholegenome.png", plot = p, width = 6, height = 4, dpi = 300)





totaldedupseqplot <- ggplot(fordups, aes(x = deduptotal)) +
  geom_histogram(binwidth = 1000000, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "Histogram of Deduplicated Reads", x = "Number of Reads", y = "Frequency") +
  theme_bw()+xlim(0,150000000) + ylim(0,25)

ggsave("C:/Users/rjpor/Downloads/totaldeduphist.png", plot = totaldedupseqplot, width = 6, height = 4, dpi = 300)

totalseqplot <- ggplot(fordups, aes(x = Total.Sequences)) +
  geom_histogram(binwidth = 1000000, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "Histogram of Total Reads", x = "Number of Reads", y = "Frequency") +
  theme_bw()+xlim(0,150000000) + ylim(0,25)

ggsave("C:/Users/rjpor/Downloads/totalseqhist.png", plot = totalseqplot, width = 6, height = 4, dpi = 300)


fordups$duppercent <- (1 - fordups$deduptotal/fordups$Total.Sequences)
fordups$dups <- fordups$Total.Sequences - fordups$deduptotal

dupreadsplot <- ggplot(fordups, aes(x = duppercent)) +
  geom_histogram(binwidth = 0.01, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "Histogram of Duplicated Reads", x = "PCR Duplication %", y = "Frequency") +
  theme_bw()+xlim(0,1)

ggsave("C:/Users/rjpor/Downloads/dupreadshist.png", plot = dupreadsplot, width = 6, height = 4, dpi = 300)





mappedreads
colnames(mappedreads) <- c("plate_well","Mapped", "Unmapped")


fordups

setDT(mappedreads)
setDT(fordups)
setkey(fordups, plate_well)
setkey(mappedreads,  plate_well)
fordupsmapped <- merge(mappedreads, fordups)



readsxdups <- ggplot(fordupsmapped, aes(x = Total.Sequences, y = Mapped, group= plate_well)) +
  geom_point() + scale_y_log10(limits = c(1, 1000000000)) + scale_x_log10(limits = c(1, 1000000000)) + 
  labs(title = "Total Reads x Mapped Reads", x = "Total Reads", y = "Mapped Reads") 
readsxdups

ggsave("C:/Users/rjpor/Downloads/readsxdups.png", plot = readsxdups, width = 6, height = 4, dpi = 300)


readsxduppercent <- ggplot(fordupsmapped, aes(x = Total.Sequences, y = duppercent, group= plate_well)) +
  geom_point() + scale_x_log10(limits = c(1, 1000000000)) + 
  labs(title = "Total Reads x Duplicate percent", x = "Total Reads", y = "Duplicate percent") 
readsxduppercent

ggsave("C:/Users/rjpor/Downloads/readsxduppercent.png", plot = readsxduppercent, width = 6, height = 4, dpi = 300)

Mappedxduppercent <- ggplot(fordupsmapped, aes(x = Mapped, y = duppercent, group= plate_well)) +
  geom_point() + scale_x_log10(limits = c(1, 1000000000)) + 
  labs(title = "Mapped Reads x Duplicate percent", x = "Mapped Reads", y = "Duplicate percent") 
Mappedxduppercent

ggsave("C:/Users/rjpor/Downloads/Mappedxduppercent.png", plot = Mappedxduppercent, width = 6, height = 4, dpi = 300)


fordupsmapped$mappedovertotal <- fordupsmapped$Mapped/fordupsmapped$Total.Sequences

Mappedovertotal <- ggplot(fordupsmapped, aes(x = mappedovertotal)) +
  geom_histogram(binwidth = 0.01, fill = "blue", color = "black", alpha = 0.7) + xlim(0,1)+
  labs(title = "Histogram of Mapped reads over Total reads", x = "Mapped/TotalReads", y = "Frequency") +
  theme_bw()

ggsave("C:/Users/rjpor/Downloads/Mappedovertotal.png", plot = Mappedovertotal, width = 6, height = 4, dpi = 300)





















#For species

counts <- read.csv("C:/Users/rjpor/Downloads/reference_counts.csv")
species <- read.csv("C:/Users/rjpor/Downloads/species.csv")

df_summary <- counts %>%
  group_by(Filename) %>%
  mutate(Total_Count = sum(Count)) %>%  # Calculate total count per file
  mutate(Normalized_Count = Count / Total_Count)  # Divide each count by total count

df_summary <- as.data.table(df_summary)

df_summaryeuprep <- subset(df_summary, Total_Count >= 40)

df_summarysubeteupulobtusa <- subset(df_summaryeuprep, Reference == "ElvisCOI" | Reference == "mtdna_D8_119")

df_summarysubeteupulobtusa$Normalized_Count <- as.numeric(df_summarysubeteupulobtusa$Normalized_Count)

df_summarysubeteupulobtusawide <- dcast(
  df_summarysubeteupulobtusa[Reference %in% c("ElvisCOI", "mtdna_D8_119")],
  Filename ~ Reference,
  value.var = "Normalized_Count",
  fun.aggregate = mean   # if there is only one row per Filename/Reference, mean == that value
)

setDT(df_summarysubeteupulobtusawide)  # convert in-place



df_summarysubeteupulobtusawide[is.na(df_summarysubeteupulobtusawide)] <- 0

df_summarysubeteupulobtusawidenoothers <- subset(df_summarysubeteupulobtusawide, mtdna_D8_119 >= 0.3 | ElvisCOI >= 0.3)

df_summarysubeteupulobtusawidenoothers[, Species :=
                                 ifelse(ElvisCOI > mtdna_D8_119, "D. obtusa",
                                        ifelse(mtdna_D8_119 > ElvisCOI, "D. pulex", "tie"))
]


onlyeuplot <- ggplot(subset(df_summarysubeteupulobtusawidenoothers), aes(x = mtdna_D8_119, y = ElvisCOI, col = Species)) +
  geom_point() +
  labs(title = "Reference Counts", x = "D. pulex counts", y = "D. obtusa counts")+ theme_bw() 
onlyeuplotlegend <- onlyeuplot+ 
  theme(strip.text.x = element_blank(),
        strip.background = element_rect(colour="white", fill="white"),
        legend.position=c(.8,.75)
  )





ggsave("C:/Users/rjpor/Downloads/onlyeuplotpercent_allratios.pdf", plot = onlyeuplotlegend, width = 4, height = 4, dpi = 300)







filtered_plot80_onlyeuro <- ggplot(subset(filtered0.80, Reference == "ElvisCOI" | Reference == "mtdna_D8_119"), aes(x = Total_Count, y = Normalized_Count, col= Reference)) +
  geom_point() + scale_x_log10(limits = c(10, 100000)) + ylim(0.8,1)+
  labs(title = "Reference counts", x = "Total counts", y = "Normalized Counts") 
filtered_plot80_onlyeuro
ggsave("C:/Users/rjpor/Downloads/filtered_plot80_onlyeuro.png", plot = filtered_plot80_onlyeuro, width = 6, height = 4, dpi = 300)






filtered0.40 <- subset(df_summary, Total_Count >= 20 & Normalized_Count >= 0.40)

Filteredhist40 <- ggplot(filtered0.40, aes(x = Count)) +
  geom_histogram(binwidth = 0.1, fill = "blue", color = "black", alpha = 0.7) + scale_x_log10(limits = c(1, 1000000)) +
  labs(title = "Histogram read counts", x = "Counts", y = "Frequency") +
  theme_bw()
Filteredhist40


filtered_plot40 <- ggplot(filtered0.40, aes(x = Total_Count, y = Normalized_Count, col= Reference)) +
  geom_point() + scale_x_log10(limits = c(10, 100000)) +
  labs(title = "Reference counts", x = "Total counts", y = "Normalized Counts") 
filtered_plot40
ggsave("C:/Users/rjpor/Downloads/filtered_plot40.png", plot = filtered_plot40, width = 6, height = 4, dpi = 300)



filtered0.80 <- subset(df_summary, Total_Count >= 40 & Normalized_Count >= 0.80)

filtered_plot80 <- ggplot(filtered0.80, aes(x = Total_Count, y = Normalized_Count, col= Reference)) +
  geom_point() + scale_x_log10(limits = c(10, 100000)) +
  labs(title = "Reference counts", x = "Total counts", y = "Normalized Counts") 
filtered_plot80
ggsave("C:/Users/rjpor/Downloads/filtered_plot80.png", plot = filtered_plot80, width = 6, height = 4, dpi = 300)

unfiltered <- ggplot(df_summary, aes(x = Total_Count, y = Normalized_Count, col= Reference)) +
  geom_point() + scale_x_log10(limits = c(10, 100000)) +
  labs(title = "Reference counts", x = "Total counts", y = "Normalized Counts") 
unfiltered
ggsave("C:/Users/rjpor/Downloads/unfiltered.png", plot = unfiltered, width = 6, height = 4, dpi = 300)




df_summaryonlyobtusa <- subset(df_summary, Reference == "ElvisCOI")
df_summaryonlyobtusa <- df_summaryonlyobtusa[,c(1,2,3)]

df_summaryonlypulex <- subset(df_summary, Reference == "mtdna_D8_119")
df_summaryonlypulex <- df_summaryonlypulex[,c(1,2,3)]

df_summaryonlyeu <- full_join(df_summaryonlypulex,df_summaryonlyobtusa, by = "Filename")

df_summaryonlyeu$Count.x[is.na(df_summaryonlyeu$Count.x)] <- 0
df_summaryonlyeu$Count.y[is.na(df_summaryonlyeu$Count.y)] <- 0
df_summaryonlyeu <- df_summaryonlyeu[startsWith(df_summaryonlyeu$Filename, "/scratch/rjp5nc/UK2022_2024/allshortreads/counts/RobertUK"), ]
df_summaryonlyeu <- subset(df_summaryonlyeu, Count.x >= 40 | Count.y >= 40)

df_summaryonlyeu <- df_summaryonlyeu %>%
  mutate(Species = case_when(
    Count.x > Count.y ~ "D. pulex",
    Count.x < Count.y ~ "D. obtusa",
    TRUE ~ NA_character_  # for ties or NAs
  ))

onlyeuplot <- ggplot(subset(df_summaryonlyeu, Count.x >= 40 | Count.y >= 40), aes(x = Count.x, y = Count.y, col= Species)) +
  geom_point() + scale_x_log10(limits = c(-5, 100000))+ scale_y_log10(limits = c(-1, 100000))+
  labs(title = "Reference Counts", x = "D. pulex counts", y = "D. obtusa counts")+ theme_bw() 
  onlyeuplot
  
  
  
  
  sum(df_summaryonlyeu$Count.x == 0 & (df_summaryonlyeu$Count.y >= 40))
  sum(df_summaryonlyeu$Count.x >= 40 & (df_summaryonlyeu$Count.y == 0))
  
  subset(df_summaryonlyeu, Species == "D. pulex")
  subset(df_summaryonlyeu, Species == "D. obtusa")
  
  
  
  subset(df_summaryonlyeu, Count.x >= 500 & Count.y >= 500)
  smallcounts <- subset(df_summaryonlyeu, Count.x <= 1 | Count.y <= 1)
  smallcounts <- subset(df_summaryonlyeu, Count.x <= 0)
  smallcounts <- subset(smallcounts, Count.x >= 40 | Count.y >= 40)
  
ggsave("C:/Users/rjpor/Downloads/onlyeuplot.png", plot = onlyeuplot, width = 6, height = 4, dpi = 300)



df_summaryonlyeu$totalpuob <- df_summaryonlyeu$Count.x + df_summaryonlyeu$Count.y

df_summaryonlyeu$percentpulex <- df_summaryonlyeu$Count.x / df_summaryonlyeu$totalpuob
df_summaryonlyeu$percentobtusa <- df_summaryonlyeu$Count.y / df_summaryonlyeu$totalpuob



onlyeuplot <- ggplot(subset(df_summaryonlyeu, Count.x >= 40 | Count.y >= 40), aes(x = percentpulex, y = percentobtusa, col= Species)) +
  geom_point() +
  labs(title = "Reference Counts", x = "D. pulex counts", y = "D. obtusa counts")+ theme_bw() 
onlyeuplotlegend <- onlyeuplot+ 
  theme(strip.text.x = element_blank(),
        strip.background = element_rect(colour="white", fill="white"),
        legend.position=c(.8,.75)
  )

ggsave("C:/Users/rjpor/Downloads/onlyeuplotpercent.pdf", plot = onlyeuplotlegend, width = 4, height = 4, dpi = 300)




filtered_plot80_onlyeuro <- ggplot(subset(filtered0.80, Reference == "ElvisCOI" | Reference == "mtdna_D8_119"), aes(x = Total_Count, y = Normalized_Count, col= Reference)) +
  geom_point() + scale_x_log10(limits = c(10, 100000)) + ylim(0.8,1)+
  labs(title = "Reference counts", x = "Total counts", y = "Normalized Counts") 
filtered_plot80_onlyeuro
ggsave("C:/Users/rjpor/Downloads/filtered_plot80_onlyeuro.png", plot = filtered_plot80_onlyeuro, width = 6, height = 4, dpi = 300)






filtered0.40_sep <- filtered0.40 %>%
  separate(Filename, into = c("a", "b","c","d","e","f", "g"), sep = "/", remove = FALSE)
filtered0.40_sep <- filtered0.40_sep %>%
  separate(g, into = c("well", "h"), sep = ".counts.", remove = FALSE)
filtered0.40_sep <- filtered0.40_sep[,c(9,11,12,13,14)]

filtered0.80_sep <- filtered0.80 %>%
  separate(Filename, into = c("a", "b","c","d","e","f", "g"), sep = "/", remove = FALSE)
filtered0.80_sep <- filtered0.80_sep %>%
  separate(g, into = c("well", "h"), sep = ".counts.", remove = FALSE)
filtered0.80_sep <- filtered0.80_sep[,c(9,11,12,13,14)]


filtered0.70 <- subset(df_summary, Total_Count >= 20 & Normalized_Count >= 0.70)
filtered0.70_sep <- filtered0.70 %>%
  separate(Filename, into = c("a", "b","c","d","e","f", "g"), sep = "/", remove = FALSE)
filtered0.70_sep <- filtered0.70_sep %>%
  separate(g, into = c("well", "h"), sep = ".counts.", remove = FALSE)
filtered0.70_sep <- filtered0.70_sep[,c(9,11,12,13,14)]



merged_species <- left_join(filtered0.70_sep, species, by = c("Reference" = "coi_accession"))

write.csv(merged_species, "C:/Users/rjpor/Downloads/filtered_species0.70.csv")







#species for 502 wells properly assigned in 0.80 - 5 species. Below was for 0.40, where there were more

unique(merged_species$species)

ambiguaparvula <- subset(merged_species, species == "Daphnia ambigua" | species == "Daphnia parvula")

others <- subset(merged_species, species == "Daphnia pileata" | species == "Daphnia catawba" | species == "Daphnia longispina")

obtusaNA2 <- subset(merged_species, species == "Daphnia obtusaNA2")

Eudappul <- subset(merged_species, species == "EU Daphnia pulex")

Eudapobtusa <- subset(merged_datawithclone, Reference == "ElvisCOI")
Eudappul <- subset(merged_datawithclone, Reference == "mtdna_D8_119")
ambigua <- subset(merged_datawithclone, Reference == "AF523699.1")

catawba <- subset(merged_species, species == "Daphnia catawba")





filtered_plot40line <- ggplot(filtered0.40, aes(x = Reference, y = Normalized_Count, group=Filename)) +
  geom_point() + geom_line() +
  labs(title = "Reference counts", x = "Total counts", y = "Normalized Counts") 
filtered_plot40line
ggsave("C:/Users/rjpor/Downloads/filtered_plot40line.png", plot = filtered_plot40line, width = 6, height = 4, dpi = 300)


filtered_plot80line <- ggplot(filtered0.80, aes(x = Reference, y = Normalized_Count, group=Filename)) +
  geom_point() + geom_line() +
  labs(title = "Reference counts", x = "Total counts", y = "Normalized Counts") 
filtered_plot80line
ggsave("C:/Users/rjpor/Downloads/filtered_plot80line.png", plot = filtered_plot80line, width = 6, height = 4, dpi = 300)


filtered_plotnofiltline <- ggplot(df_summary, aes(x = Reference, y = Normalized_Count, group=Filename)) +
  geom_point() + geom_line() +
  labs(title = "Reference counts", x = "Total counts", y = "Normalized Counts") +
  theme(
    strip.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.spacing = unit(1, "lines")
  )
filtered_plotnofiltline
ggsave("C:/Users/rjpor/Downloads/filtered_plotnofiltline.png", plot = filtered_plotnofiltline, width = 6, height = 4, dpi = 300)


merged_species <- read.csv("C:/Users/rjpor/Downloads/filtered_species0.20.csv")

wellxclone <- read.csv("C:\\Users\\rjpor\\Downloads\\wellxclone.csv")


merged_datawithclone <- left_join(filtered0.80_sep, wellxclone, by = c("well" = "plate_well"))



Rockpool <- subset(merged_datawithclone, location == "Rockpool")

Rockpool <- as.data.table(Rockpool)

ggplot(Rockpool, aes(x = Total_Count)) +
  geom_histogram(binwidth = 10, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "Histogram of Read Depth", x = "Average Read Depth", y = "Frequency") +
  theme_bw()+ xlim(0,1000)


Rockpool100 <- subset(Rockpool, Reference == "AY380446.1")

Rockpoolinds <- Rockpool[, .N, by=list(clone, date, Reference)]

usethis <- Rockpool100$well

usethis <- as.data.frame(usethis)

write.csv(usethis,"C:\\Users\\rjpor\\Downloads\\Rockpoolobtusainds.csv")





















metrics <- read.csv("C:/Users/rjpor/Downloads/metricsoutput.csv")


ggplot(metrics, aes(x = PERCENT_DUPLICATION)) +
  geom_histogram(binwidth = .01, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "Histogram of Duplication Rate", x = "Dups", y = "Frequency") +
  theme_bw() + xlim(0,1)

readdepth <- read.csv("C:/Users/rjpor/Downloads/readdepth.csv")

ggplot(readdepth, aes(x = Average_Coverage)) +
  geom_histogram(binwidth = 1, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "Histogram of Duplication Rate", x = "Dups", y = "Frequency") +
  theme_bw() + xlim(0,40)

221523/150000






relatedness <- read.csv("C:/Users/rjpor/Downloads/estimated_relatedness_not_gilmer.csv")

cloneshort <- locationandclone[,c(1,2,3)]


merged_data6 <- left_join(relatedness, cloneshort, by = c("ID1" = "Well"))
merged_data7 <- left_join(merged_data6, cloneshort, by = c("ID2" = "Well"))


#ibd <- snpgdsIBDMoM(genofile, sample.id=gilmer_good_samples, snp.id=snp.id,
#                    maf=0.05, missing.rate=0.05, num.thread=2, verbose = TRUE, autosome.only = FALSE)

#relatedness <- snpgdsIBDSelection(ibd)

ggplot(merged_data7, aes(x = kinship)) +
  geom_histogram(binwidth = 0.01, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "kinship", x = "kinship", y = "Frequency") + xlim(0.0,0.55)+
  theme_bw()

ggplot(merged_data7, aes(x = k1)) +
  geom_histogram(binwidth = 0.01, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "kinship", x = "kinship", y = "Frequency") + xlim(0.0,1.05)+
  theme_bw()

ggplot(merged_data7, aes(x = k0, y=kinship)) +
  geom_point() +
  labs(title = "kinship", x = "k0", y = "kinship") +
  theme_bw()

ggplot(merged_data7, aes(x = k0, y=k1)) +
  geom_point() +
  labs(title = "kinship", x = "k0", y = "k1") +
  theme_bw()


k1of1 <- subset(merged_data7, k1 =="1")





ggplot(k1of1, aes(x = k0, y=kinship)) +
  geom_point() +
  labs(title = "kinship", x = "k0", y = "kinship") +
  theme_bw()

highkinship<-subset(merged_data7, kinship >= 0.4)






paramfile <- read.csv("C:/Users/rjpor/Downloads/robert_paramfile.txt", header = TRUE)



pulexparam <- subset(paramfile, Species == "Daphnia pulex" & Continent == "Europe")



write.csv(pulexparam, "C:/Users/rjpor/Downloads/robert_pulexparamfile.txt")









pca_data <- read.csv("C:/Users/rjpor/Downloads/seqarray_pca_not_gilmer.csv")


merged_data5 <- left_join(pca_data, locationandclone, by = c("Sample" = "Well"))



PCA_no_gilmer <- ggplot(subset(merged_data5), aes(x = PC1, y = PC2, shape=accuratelocation, color=date)) +
  geom_point() +
  labs(title = "PCA of Genomic Data", x = "PC1", y = "PC2") +
  theme_bw()


ggsave("C:/Users/rjpor/Downloads/PCA_no_gilmer.png", plot = PCA_no_gilmer, width = 9, height = 6, dpi = 300)






pca_data <- read.csv("C:/Users/rjpor/Downloads/seqarray_pca_gilmer.csv")


merged_data5 <- left_join(pca_data, locationandclone, by = c("Sample" = "Well"))



PCA_no_gilmer <- ggplot(subset(merged_data5), aes(x = PC1, y = PC2, shape=accuratelocation, color=date)) +
  geom_point() +
  labs(title = "PCA of Genomic Data", x = "PC1", y = "PC2") +
  theme_bw()


ggsave("C:/Users/rjpor/Downloads/PCA_gilmer.png", plot = PCA_no_gilmer, width = 9, height = 6, dpi = 300)






df_summaryobtusapulex <- subset(df_summary, Reference == "mtdna_D8_119" | Reference == "ElvisCOI")   

df_summaryobtusapulexsep <- df_summaryobtusapulex %>%
  separate(Filename, into = c("a", "b","c","d","e","f", "g"), sep = "/", remove = FALSE)
df_summaryobtusapulexsep <- df_summaryobtusapulexsep %>%
  separate(g, into = c("well", "h"), sep = ".counts.", remove = FALSE)
df_summaryobtusapulexsep <- df_summaryobtusapulexsep[,c(9,11,12,13,14)]



df_summarypulexsep <- subset(df_summaryobtusapulexsep, Reference == "mtdna_D8_119")
df_summaryobtusasep <- subset(df_summaryobtusapulexsep, Reference == "ElvisCOI")


df_summary_wide <- full_join(df_summarypulexsep,df_summaryobtusasep, by = c("well" = "well"))

df_summary_wide <- df_summary_wide %>%
  separate(well, into = c("plate", "well"), sep = "_", remove = FALSE)


df_summary_wide <- subset(df_summary_wide, plate == "RobertUK" & well != "B11" & well != "D9")

df_summary_wideplot <- ggplot(subset(df_summary_wide), aes(x = Normalized_Count.x, y = Normalized_Count.y, col=Total_Count.x)) +
  geom_point() +
  labs(title = "obtusa x pulex", x = "percent pulex", y = "percent obtusa") +
  theme_bw()
df_summary_wideplot

setDT(df_summary_wide)
df_summary_wide$Normalized_Count.y[is.na(df_summary_wide$Normalized_Count.y)] <- 0

df_summary_wideplot <- ggplot(subset(df_summary_wide, Total_Count.x >= 40 & Normalized_Count.x >= 0.7 | Total_Count.y >= 40 & Normalized_Count.y >=0.7), aes(x = Normalized_Count.x*100, y = Normalized_Count.y*100, col=species)) +
  geom_point() +
  labs(title = "Mapping to COI, total counts > 100", x = "D. pulex %", y = "D. obtusa %") +
  theme_bw()
df_summary_wideplot

ggsave("C:/Users/rjpor/Downloads/df_summary_wideplot.png", plot = df_summary_wideplot, width = 9, height = 6, dpi = 300)

subset(df_summary_wide, Total_Count.x >= 40 & Normalized_Count.y <= 0.7 & Normalized_Count.x <=0.7)


df_summary_wide[, species :=
                  fifelse(Normalized_Count.y > 0.80, "D. obtusa",
                          fifelse(Normalized_Count.x > 0.80, "D. pulex", NA_character_))
]


df_summary_wideplot_colored <- ggplot(subset(df_summary_wide, Total_Count.x >= 40 & Normalized_Count.x >= 0.8 | Total_Count.y >= 40 & Normalized_Count.y >=0.8), aes(x = Normalized_Count.x*100, y = Normalized_Count.y*100, col=species)) +
  geom_point() +
  labs(title = "Mapping to COI, total counts > 100", x = "D. pulex %", y = "D. obtusa %") +
  theme_bw()
df_summary_wideplot_colored



df_summary_wideplot_coloredlegend <- df_summary_wideplot_colored+ 
  theme(strip.text.x = element_blank(),
        strip.background = element_rect(colour="white", fill="white"),
        legend.position=c(.8,.75)
  )


ggsave("C:/Users/rjpor/Downloads/Fig2A.png", plot = df_summary_wideplot_coloredlegend, width = 9, height = 6, dpi = 300)


df_summary_wideplot_coloredlegend + p


#.nwk tree redone



library(treeio)
library(ggtree)
library(data.table)

pca_data <- read.csv("C:/Users/rjpor/Downloads/2022seqPCA.csv")
treefile <- "C:/Users/rjpor/Downloads/mito_tree.nwk"
tree <- read.tree(treefile)


library(ape)





treefile <- "C:/Users/rjpor/Downloads/dna.tree"
tree <- read.nexus(treefile)

# View the tree
plot(tree)


#tree$tip.label <- pca_data$Well[match(tree$tip.label, pca_data$Sample)]

tree_plot <- ggtree(tree, options(ignore.negative.edge=TRUE))

# Extract tree data
tree_data <- tree_plot$data

#write.csv(pca_datawithclone, "C:\\Users\\rjpor\\Downloads\\pca_withclone.csv")

# Perform the merge using left_join (same as merge but more readable)
merged_data2 <- left_join(tree_data, XXXXXXXXXXX, by = c("label" = "Well"))

merged_data2 <- merged_data2 %>%
  separate(label, into = c("plate", "column"), sep = ".", remove = FALSE)


merged_data2$Plate
unique(merged_data2$Plate)


#tree_plot <- ggtree(merged_data2, layout="circular")
tree_plot <- ggtree(tree_data, layout="rectangular")

pdf("C:\\Users\\rjpor\\Downloads\\treemito_nexus_rect.pdf", width=20, height=20)

tree_plot +
  geom_tiplab() +           #Use new labels and color by group
  #scale_color_manual(values = c("us" = "blue", "UK" = "red", "Gilmer" = "green")) +  # Custom colors
  theme_tree2()  # Optional tree theme

dev.off()


pdf("C:\\Users\\rjpor\\Downloads\\treemito_nwk_rect.pdf", width=20, height=20)

tree_plot +
  geom_tiplab(aes(label = clone, color = date)) +  # Use new labels and color by group
  #scale_color_manual(values = c("us" = "blue", "UK" = "red", "Gilmer" = "green")) +  # Custom colors
  theme_tree2()  # Optional tree theme

dev.off()


pattern <- paste(obtusalist$Filename, collapse = "|")

tips_to_keep <- setdiff(
  tree$tip.label[grepl(pattern, tree$tip.label)],
  "mtdna.RobertUK_B11"
)


tree_rooted <- root(tree, outgroup = "mtdna.Daphnia_magna", resolve.root = TRUE)
tips_to_keep <- tree$tip.label[grepl("RobertUK", tree$tip.label) | tree$tip.label == "mtdna.Daphnia_magna"]





# Plot and save to PDF
pdf("C:\\Users\\rjpor\\Downloads\\treemito_nexus_rect_rooted_magna.pdf", width = 20, height = 60)

ggtree(tree_rooted) +
  geom_tiplab() +
  theme_tree2()

dev.off()


tree_rooted <- root(tree, outgroup = "mtdna_D8_119.index_8nt_N7_484-index_8nt_N5_700_SL378987_Spring_2018_Pond22_21", resolve.root = TRUE)
tree_rooted <- root(tree, outgroup = "mtdna.RobertUK_E12", resolve.root = TRUE)



tips_to_keep <- setdiff(
  tree$tip.label[grepl("RobertUK", tree$tip.label) | grepl("Pond22_21", tree$tip.label)],
  "mtdna.RobertUK_B11"
)







tips_to_keep <- setdiff(
  tree$tip.label[grepl("RobertUK", tree$tip.label) | grepl("Rockpool3_A11", tree$tip.label)],
  "mtdna.RobertUK_B11"
)

obtusalist <- subset(df_summary, Normalized_Count >= 0.8 & Reference == "ElvisCOI")
obtusalist$Filename <- gsub("/scratch/rjp5nc/UK2022_2024/allshortreads/counts/", "", obtusalist$Filename)
obtusalist$Filename <- gsub(".counts.txt", "", obtusalist$Filename)
obtusalist$Reference <- gsub("ElvisCOI", "D. obtusa", obtusalist$Reference)


pulexlist <- subset(df_summary, Normalized_Count >= 0.8 & Reference == "mtdna_D8_119")
pulexlist$Filename <- gsub("/scratch/rjp5nc/UK2022_2024/allshortreads/counts/", "", pulexlist$Filename)
pulexlist$Filename <- gsub(".counts.txt", "", pulexlist$Filename)
pulexlist$Reference <- gsub("mtdna_D8_119", "D. pulex", pulexlist$Reference)

ambigualist <- df_summary

ambigualist <- subset(df_summary, Normalized_Count >= 0.8 & Reference == "AF523699.1")
ambigualist$Filename <- gsub("/scratch/rjp5nc/UK2022_2024/allshortreads/counts/", "", ambigualist$Filename)
ambigualist$Filename <- gsub(".counts.txt", "", ambigualist$Filename)
ambigualist$Reference <- gsub("AF523699.1", "D. ambigua", ambigualist$Reference)

alllist <-rbind(pulexlist, obtusalist,ambigualist)
alllist$Name <- paste("mtdna.",alllist$Filename)
alllist$Name <- gsub(" ", "", alllist$Name)
alllist$Name <- factor(alllist$Name)

alllist <- alllist[,c(1,2,6)]

cloneidpluswell <- pca_datawithclone[,c(1,7,8,9,10,11,12)]

cloneidpluswell2 <- left_join(alllist, cloneidpluswell, by = c("Filename" = "Well"))

treefile <- "C:/Users/rjpor/Downloads/dna.tree"

tree <- read.nexus(treefile)

tree_rooted <- root(tree, outgroup = "mtdna.Rockpool3_A11", resolve.root = TRUE)

pruned_tree <- drop.tip(tree_rooted, setdiff(tree$tip.label, tips_to_keep))

pruned_tree_df2 <- left_join(pruned_tree, alllist, by = c("label" = "Name"))
pruned_tree_df3 <- left_join(pruned_tree_df2, cloneidpluswell2, by = c("Filename","Reference"))

pruned_tree_df3data <- pruned_tree_df3$data

tips_to_drop <- cloneidpluswell2$Filename[cloneidpluswell2$clone == "BLANK"] 
tips_to_drop_prefixed <- paste0("mtdna.", tips_to_drop)
pruned_tree_no_blank <- drop.tip(pruned_tree_df3, tips_to_drop_prefixed)


ggtree(pruned_tree_no_blank, layout = "circular") +
  geom_tiplab(aes(label = clone, col=Reference)) +
theme_tree2()


pdf("C:\\Users\\rjpor\\Downloads\\treemito_nexus_Uk_rooted_ambigua.pdf", width = 10, height = 10)


ggtree(pruned_tree_no_blank, layout = "circular") +
  geom_tiplab(aes(label = clone, col = Reference)) +
  scale_color_manual(values = c(
    "D. ambigua" = "green",
    "D. pulex" = "#00BFC4",
    "D. obtusa" = "#F8766D"
  )) +
  theme_tree2()

dev.off()







ggtree(pruned_tree_data3, layout = "rectangular") +
  geom_tiplab(aes(label = clone, col=Reference)) +
  theme_tree2()





tree_rooted <- root(tree, outgroup = "mtdna.RobertUK_A6", resolve.root = TRUE)

tips_to_drop <- cloneidpluswell2$Filename[cloneidpluswell2$clone == "BLANK" | cloneidpluswell2$Reference == "D. pulex"| cloneidpluswell2$accuratelocation == "P759"]
tips_to_drop_prefixed <- paste0("mtdna.", tips_to_drop)

pruned_tree <- drop.tip(tree_rooted, setdiff(tree$tip.label, tips_to_keep))

pruned_tree_df2 <- left_join(pruned_tree, alllist, by = c("label" = "Name"))
pruned_tree_df3 <- left_join(pruned_tree_df2, cloneidpluswell2, by = c("Filename","Reference"))

pruned_tree_no_blank <- drop.tip(pruned_tree_df3, tips_to_drop_prefixed)

pdf("C:\\Users\\rjpor\\Downloads\\treemito_nexus_obtusa_only.pdf", width = 10, height = 10)

ggtree(pruned_tree_no_blank, layout = "circular") +
  geom_tiplab(aes(label = clone, col=accuratelocation)) +
  theme_tree2()

dev.off()






#EU pulex

treefile <- "C:/Users/rjpor/Downloads/dna.tree"

tree <- read.nexus(treefile)

tree_rooted <- root(tree, outgroup = "mtdna.RobertUK_D6", resolve.root = TRUE)

tips_to_drop <- cloneidpluswell2$Filename[cloneidpluswell2$clone == "BLANK" | cloneidpluswell2$Reference == "D. obtusa"| cloneidpluswell2$accuratelocation == "P759"]
tips_to_drop_prefixed <- paste0("mtdna.", tips_to_drop)

pruned_tree <- drop.tip(tree_rooted, setdiff(tree$tip.label, tips_to_keep))

pruned_tree_df2 <- left_join(pruned_tree, alllist, by = c("label" = "Name"))
pruned_tree_df3 <- left_join(pruned_tree_df2, cloneidpluswell2, by = c("Filename","Reference"))

pruned_tree_no_blank <- drop.tip(pruned_tree_df3, tips_to_drop_prefixed)

pdf("C:\\Users\\rjpor\\Downloads\\treemito_nexus_pulex_only.pdf", width = 10, height = 10)

ggtree(pruned_tree_no_blank, layout = "circular") +
  geom_tiplab(aes(label = clone, col=accuratelocation)) +
  theme_tree2()

dev.off()






samples_to_remove <- c(
  "mtdna.SRR14370488", "mtdna.SRR14370491", "mtdna.SRR14370483",
  "mtdna.SRR14370492",
  "mtdna.SRR14370486",
  "mtdna.SRR14370489",
  "mtdna.SRR14370484",
  "mtdna.SRR14370485",
  "mtdna.SRR14370482",
  "mtdna.SRR24967773",
  "mtdna.SRR14370487",
  "mtdna.SRR14370481",
  "mtdna.SRR14370490")

tree_pruned <- drop.tip(tree_rooted, samples_to_remove)


tree_plot_circ <- ggtree(tree_rooted, layout="circular")

pdf("C:\\Users\\rjpor\\Downloads\\treemito_nwk_circ.pdf", width=200, height=200)

tree_plot_circ  +
  geom_tiplab() +           #Use new labels and color by group
  #scale_color_manual(values = c("us" = "blue", "UK" = "red", "Gilmer" = "green")) +  # Custom colors
  theme_tree2()  # Optional tree theme

dev.off()







aaaa <- tree_plot +
  geom_tiplab(aes(label = clone, color = date)) +  # Use new labels and color by group
  #scale_color_manual(values = c("us" = "blue", "UK" = "red", "Gilmer" = "green")) +  # Custom colors
  theme_tree2()  

ggsave("C:/Users/rjpor/Downloads/treecirdate.png", plot = aaaa, width = 20, height = 20, dpi = 300)

aaab <- tree_plot +
  geom_tiplab(aes(label = clone, color = accuratelocation)) +  # Use new labels and color by group
  #scale_color_manual(values = c("us" = "blue", "UK" = "red", "Gilmer" = "green")) +  # Custom colors
  theme_tree2()  

ggsave("C:/Users/rjpor/Downloads/treecir.png", plot = aaab, width = 20, height = 20, dpi = 300)


