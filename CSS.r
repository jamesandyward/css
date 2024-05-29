##### Script to generate inputs for CSS and then run CSS using the MINOTAUR package, and finally to generate manhattan plots of raw and smoothed CSS values #####

# Load libraries 
library(data.table)
library(dplyr)
library(tidyr)
library(MINOTAUR)

# Fst 
# Fst values can be generated using VCFtools and the --weir-fst-pop function, for example:
# vcftools --vcf ../data/wgs_data.vcf --weir-fst-pop ../phased_data/popmaps/test_pop.txt --weir-fst-pop ../phased_data/popmaps/ref_pop.txt  --out est_ref_Fst_ind_SNP
# This will generate a file with Fst values for each SNP

#Read in the Fst data and format
data_fst <- fread("FST.weir.fst")
colnames(data_fst) <- c("chromosome", "position", "fst")
data_fst <- data_fst[order(data_fst$chromosome, data_fst$position),] # Order data by chromosome and position
data_fst$SNPid <- paste0(data_fst$chromosome,":",data_fst$position) # Append a SNPid (chromosome:position) to the dataframe

# Write the dataframe to a file in the event of a crash so data can be reloaded quickly
write.table(data_fst, "CSS_results/Fst.txt", sep = "\t", row.names = FALSE)
data_fst <- fread("CSS_results/Fst.txt", sep = "\t")

# ∆SAF
# ∆SAF values can be calculated using the VCFtools --freq option and some manipulation in R
# Read in freq data generated using VCFtools --freq option
# Generate the files for the selected population and the non-selected population
selected_frq <- fread("SAF.frq", header = FALSE, skip = 1)
non_selected_frq <- fread("SAF.frq", header = FALSE, skip = 1)

# Set column names
colnames(selected_frq) <- c("CHROM", "POS", "N_ALLELES", "N_CHR", "Allele_Freq1", "Allele_Freq2")
colnames(non_selected_frq) <- c("CHROM", "POS", "N_ALLELES", "N_CHR", "Allele_Freq1", "Allele_Freq2")

# Separate the last two columns into allele and frequency columns
# This will take some time to complete
# Could split into chromosomes to speed up, but I have not done that here
selected_frq <- selected_frq %>%
  separate(Allele_Freq1, into = c("Allele1", "Freq1"), sep = ":", convert = TRUE) %>%
  separate(Allele_Freq2, into = c("Allele2", "Freq2"), sep = ":", convert = TRUE)

non_selected_frq <- non_selected_frq %>%
  separate(Allele_Freq1, into = c("Allele1", "Freq1"), sep = ":", convert = TRUE) %>%
  separate(Allele_Freq2, into = c("Allele2", "Freq2"), sep = ":", convert = TRUE)

# Identify major allele and its frequency in the selected population (i.e. target population)
selected_frq <- selected_frq %>%
  mutate(MajorAllele = ifelse(Freq1 > Freq2, Allele1, Allele2),
         MajorAlleleFreq = ifelse(Freq1 > Freq2, Freq1, Freq2))

# Find the frequency of the same allele in the non-selected population
data_saf <- selected_frq %>%
  left_join(non_selected_frq, by = c("CHROM", "POS"), suffix = c("_selected", "_non_selected")) %>%
  mutate(NonSelectedAlleleFreq = ifelse(MajorAllele == Allele1_non_selected, Freq1_non_selected,
                                        ifelse(MajorAllele == Allele2_non_selected, Freq2_non_selected, NA)))

# Select only relevant columns
data_saf <- data_saf[, c("CHROM", "POS", "MajorAllele", "MajorAlleleFreq", "NonSelectedAlleleFreq")]

# Calculate the directional change in the selected allele frequency (∆SAF)
data_saf <- data_saf %>%
  mutate(DeltaSAF = MajorAlleleFreq - NonSelectedAlleleFreq)

# Standardize ∆SAF values to a normal distribution with mean 0 and standard deviation 1 as per Randhawa et al. 2014
data_saf$DeltaSAF_Z <- scale(data_saf$DeltaSAF)

# Order data by chromosome and position. Remove unnecessary columns and add SNPid
data_saf <- data_saf[order(data_saf$CHROM, data_saf$POS),]
data_saf <- data_saf[,c(1,2,7)]
colnames(data_saf) <- c("chromosome", "position", "deltaSAFz")
data_saf$SNPid <- paste0(data_saf$chromosome,":",data_saf$position)

# Write the dateframe to a file in the event of a crash so data can be reloaded quickly
write.table(data_saf, "CSS_results/SAF.txt", sep = "\t", row.names = FALSE)
data_saf <- fread("CSS_results/SAF.txt", sep = "\t")

# Remove unneccessary dataframes from the environment
rm(selected_frq, non_selected_frq)

# XPEHH
# XPEHH values are calculated using selscan 
# This is done on a chromosome by chromosome basis, for example:
# selscan --xpehh --vcf test_pop_phased_chr1.vcf.gz --vcf-ref ref_pop_phased_chr1.vcf.gz --map chr1_selscan.map --out chr1_selscan_test-ref --threads 20
# Each chromosome output is then merged into a single file using simple awk command:
# awk 'FNR==1 && NR!=1{next;}{print}' *_selscan_test-ref.xpehh.out | sort -k 1n > all_selscan_test-ref.xpehh.out

# Read in the XPEHH data and format
data_xpehh<-fread("selscan.win50.xpehh.out")
colnames(data_xpehh)[1] <- "chrpos"
data_xpehh<- separate(data_xpehh, chrpos, into = c("chromosome", "position"), sep = "-")
data_xpehh$chromosome <- as.numeric(data_xpehh$chromosome)
data_xpehh$position <- as.integer(data_xpehh$position) 
data_xpehh$XPEHHz <- scale(data_xpehh$xpehh)
data_xpehh$SNPid <- paste0(data_xpehh$chromosome,":",data_xpehh$position)
data_xpehh <- data_xpehh[,c(1,2,12,13)]
colnames(data_xpehh) <- c("chromosome", "position", "XPEHHz", "SNPid")

# Write the dateframe to a file in the event of a crash so data can be reloaded quickly
write.table(data_xpehh, "CSS_results/XPEHH.txt", sep = "\t", row.names = FALSE)
data_xpehh <- fread("CSS_results/XPEHH.txt", sep = "\t")

# Merge data frames by SNPid, select chromosome, position, SNPid and selection stat columns
merged_df <- data_fst %>%
  full_join(data_saf, by = "SNPid") %>%
  full_join(data_xpehh, by = "SNPid") %>% drop_na()
df<- merged_df[,c(1,2,4,3,7,10)]
colnames(df) <- c("chromosome", "position", "SNPid", "fst", "deltaSAFz", "XPEHHz")

# Remove dataframes no longer in use
rm(data_fst, data_saf, data_xpehh, merged_df)

# Write the dataframe to a file.
write.table(df, "CSS_results/merged_statistics.txt", sep = "\t", row.names = FALSE)

# Calculate CSS values for the data.
# Run CSS on the dataframe using the css() function from the MINOTAUR package

# Specify the columns and test with the selection statistics
columns <- 4:6

# CSS does not take covariance into account (i.e., no S argument)
df$css <- CSS(df, columns, right.tailed = rep(TRUE, length(columns)))

# Plot histogram of the statistics 
hist(df$css, breaks=25, xlab="CSS value", main="Histogram of CSS statistic")

# Write the dateframe containing the css values to a file
df <- df[,c(1,2,3,7)]
write.table(df, "CSS_results/CSS.txt", sep = "\t", row.names = FALSE)

# Run python smoothing script on the CSS.txt

# Smooth CSS statistics using the sliding_window.py script and merge the resulting output using the merge_csv.py script
smoothed_results <- fread("CSS_results/CSS_windows.csv", header = T)

# Add window midpoint and SNPid
smoothed_results$position <- with(smoothed_results, (start + end) / 2)
smoothed_results$SNPid <- paste0(smoothed_results$chromosome,":",smoothed_results$position)

# Order data by chromosome and position
smoothed_results <- smoothed_results[order(smoothed_results$chromosome, smoothed_results$start),]

# Find cluster regions
threshold_0_1_percent <- quantile(smoothed_results$mean_css, probs = 0.999, na.rm = TRUE)
threshold_1_percent <- quantile(smoothed_results$mean_css, probs = 0.99, na.rm = TRUE)
significant_snps <- smoothed_results[smoothed_results$mean_css >= threshold_0_1_percent,]
top_1_percent <- smoothed_results[smoothed_results$mean_css >= threshold_1_percent,]

# Define the function to check for at least five flanking SNPs from the top 1%
has_five_flanking_snps <- function(snp, top_1_percent) {
  start_pos <- snp$position - 100000  # 100kb window
  end_pos <- snp$position + 100000 # 100kb window
  flanking_snps <- top_1_percent[top_1_percent$chromosome == snp$chromosome & 
                                   top_1_percent$position >= start_pos & 
                                   top_1_percent$position <= end_pos, ]
  return(nrow(flanking_snps) >= 5)
}

# Identify the significant SNPs
significant_snps <- significant_snps[sapply(1:nrow(significant_snps), function(i) {
  has_five_flanking_snps(significant_snps[i,], top_1_percent)
}), ]

# Merge clusters that are less than 1Mb apart
significant_snps[, cluster := cumsum(c(1, diff(position) > 1000000)), by = chromosome]
cluster_regions <- significant_snps[, .(start = min(position), end = max(position)), by = .(chromosome, cluster)]

# Significant SNPs 
highlighted_snps <- significant_snps$SNPid

# Write cluster regions to a txt file
write.table(cluster_regions, "CSS_results/CSS_Cluster_Regions.txt", sep = "\t", row.names = FALSE)

# Manhattan Plots

# Load required libraries 
library(ggplot2)
library(dplyr)

# Load in data. 1. Raw CSS values and 2. Smoothed CSS values if you have not done so already
df <- fread("CSS_results/CSS.txt", sep = "\t", header = T)
df$chromosome <- as.numeric(df$chromosome)
smoothed_results <- fread("CSS_results/CSS_windows.csv", header = T)

# Create a new variable for alternating colors
df$color <- with(df, as.numeric(chromosome) %% 2)
smoothed_results$color <- with(smoothed_results, as.numeric(chromosome) %% 2)

# Add window midpoint and SNPid
smoothed_results$position <- with(smoothed_results, (start + end) / 2)
smoothed_results$SNPid <- paste0(smoothed_results$chromosome,":",smoothed_results$position)

# Order data by chromosome and position
smoothed_results <- smoothed_results[order(smoothed_results$chromosome, smoothed_results$position),]

# Generate cumulative max position for each chromosome to offset the position
cumulative_max_position <- df %>% 
  group_by(chromosome) %>% 
  summarise(max_position = max(position)) %>% 
  mutate(cumulative_max_position = lag(cumsum(as.numeric(max_position)), default = 0)) 

# Join cumulative_max_position back to df and smoothed_data
df <- df %>% left_join(cumulative_max_position, by = "chromosome")
smoothed_results <- smoothed_results %>% left_join(cumulative_max_position, by = "chromosome")
significant_snps <- significant_snps %>% left_join(cumulative_max_position, by = "chromosome")

# Create new position (continuous genomic coordinate) based on position and cumulative_max_position
df$continuous_position <- df$position + df$cumulative_max_position
smoothed_results$continuous_position <- smoothed_results$position + smoothed_results$cumulative_max_position
significant_snps$continuous_position <- significant_snps$position + significant_snps$cumulative_max_position

# Compute the maximum CSS for text annotation
max_p_value <- max(c(max((df$css)), max((smoothed_results$css))))

# Calculate the middle point for each chromosome
axis_set <- smoothed_results %>%
  group_by(chromosome) %>%
  summarise(center = (min(continuous_position) + max(continuous_position)) / 2)

# Write image to file
png(file="CSS_results/CSS_final.png", width=10, height=5, units="in", res=300)

# Manhattan plot
ggplot() +
  # Add raw css data points
  geom_point(data = df, 
             aes(x = continuous_position, 
                 y = css), 
             color = "grey",
             shape = 4,
             size = 0.5, 
             alpha = 0.2) +
  # Add smoothed css data points
  geom_point(data = smoothed_results, 
             aes(x = continuous_position, 
                 y = mean_css, 
                 color = as.factor(color)), 
             size = 1) +
  # Add points to highlight
  geom_point(data = significant_snps, 
             aes(x = continuous_position, 
                 y = mean_css), 
             color = "red", size = 1, alpha = 1) +
  # Axis cosmetics
  scale_x_continuous(label = axis_set$chromosome, breaks = axis_set$center) +
  scale_y_continuous(expand = c(0,1)) +
  scale_color_manual(values = rep(c("#999999", "#E69F00"), unique(length(axis_set$chromosome)))) +
  labs(x = "Chromosome", 
       y = "CSS", 
       color = "Chromosome") +
  # Add horizontal line for significance level
  geom_hline(yintercept = threshold_0_1_percent, linetype="dashed", color = "red", size = 0.5) +
  scale_x_continuous(label = axis_set$chromosome, breaks = axis_set$center) +
  scale_y_continuous(expand = c(0,0.5)) +
  scale_color_manual(values = rep(c("#1D1D74", "#eb9b05"), unique(length(axis_set$chromosome)))) +
  scale_size_continuous(range = c(0.5,3)) +
  labs(x = "Chromosome", 
       y = "CSS") +
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5, hjust = 1)
  )

dev.off()
