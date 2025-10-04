# bio4
output:
  pdf_document: default
  html_document: default
---

```{r}
# Create a folder for  assignment
dir.create("bioassignment")
setwd("bioassignment")  

# Download the files
download.file("https://github.com/ghazkha/Assessment4/raw/main/gene_expression.tsv", 
              "gene_expression.tsv")
download.file("https://github.com/ghazkha/Assessment4/raw/main/growth_data.csv", 
              "growth_data.csv")
```
```{r}
# Read the file
gene_data <- read.table("gene_expression.tsv", header = TRUE, row.names = 1)

# 1. Show first six genes
head(gene_data)

# 2. Add mean column
gene_data$mean_expression <- rowMeans(gene_data)
head(gene_data)

# 3. Top 10 genes with highest mean
top_10 <- gene_data[order(-gene_data$mean_expression), ][1:10, ]
top_10

# 4. Count genes with mean < 10
sum(gene_data$mean_expression < 10)

# 5. Make histogram
hist(gene_data$mean_expression, main = "Gene Expression Means", 
     xlab = "Mean Expression", col = "lightblue")
```
```{r}
# Read the file
growth_data <- read.csv("growth_data.csv")

# 6. Show column names 
colnames(growth_data)

# 7. Calculate mean and standard deviation
# For start (2005) at control site
start_control <- growth_data$Circumf_2005_cm[growth_data$Site == "control"]
mean(start_control)
sd(start_control)

# For end (2020) at control site
end_control <- growth_data$Circumf_2020_cm[growth_data$Site == "control"]
mean(end_control)
sd(end_control)

# For start (2005) at treatment site
start_treatment <- growth_data$Circumf_2005_cm[growth_data$Site == "treatment"]
mean(start_treatment)
sd(start_treatment)

# For end (2020) at treatment site
end_treatment <- growth_data$Circumf_2020_cm[growth_data$Site == "treatment"]
mean(end_treatment)
sd(end_treatment)

# 8. Make boxplot 
library(tidyr)
library(ggplot2)

# Convert data from wide to long format
long_data <- growth_data %>%
  pivot_longer(cols = starts_with("Circumf"),
               names_to = "Year",
               values_to = "Circumference")

# Clean up year names
long_data$Year <- gsub("Circumf_", "", long_data$Year)
long_data$Year <- gsub("_cm", "", long_data$Year)

# Create boxplot
ggplot(long_data, aes(x = Year, y = Circumference, fill = Site)) +
  geom_boxplot() +
  labs(title = "Tree Circumference by Site and Year")

# 9. Calculate mean growth over last 10 years (2010-2020)
# For control site
control_growth <- growth_data$Circumf_2020_cm[growth_data$Site == "control"] - 
                  growth_data$Circumf_2010_cm[growth_data$Site == "control"]
mean(control_growth)

# For treatment site
treatment_growth <- growth_data$Circumf_2020_cm[growth_data$Site == "treatment"] - 
                    growth_data$Circumf_2010_cm[growth_data$Site == "treatment"]
mean(treatment_growth)

# 10. t-test to compare growth 


print("Number of control trees:")
print(sum(growth_data$Site == "control"))

print("Number of treatment trees:")
print(sum(growth_data$Site == "treatment"))

# Check 
if(sum(growth_data$Site == "control") >= 2 && sum(growth_data$Site == "treatment") >= 2) {
  # Calculate growth for control group
  control_growth <- growth_data$Circumf_2020_cm[growth_data$Site == "control"] - 
                    growth_data$Circumf_2010_cm[growth_data$Site == "control"]
  
  # Calculate growth for treatment group
  treatment_growth <- growth_data$Circumf_2020_cm[growth_data$Site == "treatment"] - 
                      growth_data$Circumf_2010_cm[growth_data$Site == "treatment"]
  
  # Remove any NA values
  control_growth <- control_growth[!is.na(control_growth)]
  treatment_growth <- treatment_growth[!is.na(treatment_growth)]
  
  # Check if we still have enough data after removing NAs
  if(length(control_growth) >= 2 && length(treatment_growth) >= 2) {
    t_test_result <- t.test(control_growth, treatment_growth)
    print(t_test_result)
  } else {
    print("Not enough data after removing NA values")
    print(paste("Control growth observations:", length(control_growth)))
    print(paste("Treatment growth observations:", length(treatment_growth)))
  }
} else {
  print("Not enough trees in one or both groups for t-test")
  print(paste("Control trees:", sum(growth_data$Site == "control")))
  print(paste("Treatment trees:", sum(growth_data$Site == "treatment")))
}
```
title: "Untitled"
author: "Ann"
date: "2025-09-18"
output:
  pdf_document: default
  html_document: default
---

```{r}
library("R.utils")

URL="https://ftp.ensemblgenomes.ebi.ac.uk/pub/bacteria/release-62/fasta/bacteria_40_collection/mesomycoplasma_hyopneumoniae_gca_004768725/cds/Mesomycoplasma_hyopneumoniae_gca_004768725.ASM476872v1.cds.all.fa.gz"
#download.file(URL,destfile="Mesomycoplasma_cds.fa.gz")
#gunzip("Mesomycoplasma_cds.fa.gz")
list.files()
```
```{r}
# ------------------------------
# R Markdown code to check file location
# ------------------------------

# 1. Check R's current working directory
getwd()  # This shows the folder R is currently using

# 2. List all files in the working directory
list.files()  # Shows files in the current folder

# 3Files specify the path
folder_path <- "C:/Users/Ann/Documents/myrepo/bioassignment" 

# 4. List all files in that folder
list.files(folder_path)

# 5. Checking specific FASTA files exist
file.exists(paste0(folder_path, "/e_coli_cds.fa"))
file.exists(paste0(folder_path, "/mesomycoplasma_cds.fa"))

```

```{r}
library(seqinr)

# Read the files with exact names
ecoli_fasta <- read.fasta("ecoli_cds.fa")
Mesomycoplasma_fasta <- read.fasta("Mesomycoplasma_cds.fa")

# Calculate total CDS lengths
total_ecoli_length <- sum(sapply(ecoli_fasta, length))
total_Mesomycoplasma_length <- sum(sapply(Mesomycoplasma_fasta, length))

# Print results
cat("Total E.coli CDS length:", total_ecoli_length, "\n")
cat("Total Mesomycoplasma CDS length:", total_Mesomycoplasma_length, "\n")
library(seqinr)
# Count the sequences
num_ecoli <- length(ecoli_fasta)
num_myco <- length(Mesomycoplasma_fasta)


# Calculate total lengths
total_ecoli_length <- sum(sapply(ecoli_fasta, length))
total_Mesomycoplasma_length <- sum(sapply(Mesomycoplasma_fasta, length))

# Calculate number of genes
num_ecoli_genes <- length(ecoli_fasta)
num_Mesomycoplasma_genes <- length(Mesomycoplasma_fasta)

# Calculate average gene lengths
avg_ecoli_length <- total_ecoli_length / num_ecoli_genes
avg_Mesomycoplasma_length <- total_Mesomycoplasma_length / num_Mesomycoplasma_genes

# Create a data frame
gene_lengths_df <- data.frame(
  Organism = c("E.coli", "Mesomycoplasma"),
  Total_CDS_Length = c(total_ecoli_length, total_Mesomycoplasma_length),
  Number_of_Genes = c(num_ecoli_genes, num_Mesomycoplasma_genes),
  Average_Gene_Length = c(avg_ecoli_length, avg_Mesomycoplasma_length)
)

# Print the table
gene_lengths_df


```


```{r}
# Count the sequences
num_ecoli <- length(ecoli_fasta)
num_myco <- length(Mesomycoplasma_fasta)

# Get the total length of every sequence
total_length_ecoli <- sum(sapply(ecoli_fasta, length))
total_length_myco <- sum(sapply(Mesomycoplasma_fasta, length))

# Calculate average CDS length
avg_length_ecoli <- round(total_length_ecoli / num_ecoli, 1)
avg_length_myco <- round(total_length_myco / num_myco, 1)

# Make a table (showing total length in bp and kbp, plus average CDS length)
length_table <- data.frame(
  Organism = c("E. coli", "Mesomycoplasma"),
  Number_of_CDS = c(num_ecoli, num_myco),
  Total_Length_bp = c(total_length_ecoli, total_length_myco),
  Total_Length_kbp = round(c(total_length_ecoli/1000, total_length_myco/1000), 1),
  Average_CDS_Length = c(avg_length_ecoli, avg_length_myco)
)

# Print the table
print(length_table)

```
```{r}
# Get the length of each individual sequence
ecoli_lengths <- sapply(ecoli_fasta, length)
myco_lengths <- sapply(Mesomycoplasma_fasta, length)

# Combine for plotting
all_lengths <- c(ecoli_lengths, myco_lengths)
all_organisms <- c(rep("E. coli", length(ecoli_lengths)), rep("Mesomycoplasma", length(myco_lengths)))
plot_data <- data.frame(Length = all_lengths, Organism = all_organisms)

# Make the boxplot
boxplot(Length ~ Organism, data = plot_data,
        main = "Coding Sequence Length",
        ylab = "Length (bp)",
        xlab = "Organism",
        col = c("lightgreen", "lightcoral"))

# Calculate summary statistics (Mean and Median)
summary_table <- data.frame(
  Organism = c("E. coli", "Mesomycoplasma"),
  Mean_Length = round(c(mean(ecoli_lengths), mean(myco_lengths)), 1),
  Median_Length = c(median(ecoli_lengths), median(myco_lengths))
)

# Print the summary table
print(summary_table)

```
```{r}
# Count all bases in all sequences for each organism
all_ecoli_bases <- unlist(ecoli_fasta)
all_myco_bases <- unlist(Mesomycoplasma_fasta)

# Count frequency of A, T, G, C
ecoli_base_freq <- table(all_ecoli_bases)
myco_base_freq <- table(all_myco_bases)

# Combine frequencies into a single data frame for easier comparison
base_freq_table <- data.frame(
  Base = c("A", "T", "G", "C"),
  E_coli = as.integer(ecoli_base_freq[c("a","t","g","c")]),
  Mesomycoplasma = as.integer(myco_base_freq[c("a","t","g","c")])
)

# Print the table
print(base_freq_table)

# Make side-by-side barplots
par(mfrow = c(1, 2), mar=c(5,4,4,2)) # Two plots side-by-side, adjust margins

barplot(ecoli_base_freq[c("a","t","g","c")],
        main="E. coli Base Frequency",
        col=rainbow(4),
        ylab="Count")

barplot(myco_base_freq[c("a","t","g","c")],
        main="Mesomycoplasma Base Frequency",
        col=rainbow(4),
        ylab="Count")

```
```{r}
library(seqinr)
ecoli_cds <- read.fasta(file = "ecoli_cds.fa", seqtype = "DNA")
meso_cds <- read.fasta(file = "Mesomycoplasma_cds.fa", seqtype = "DNA")
str(head(ecoli_cds))    # Should show a list of sequences
str(head(meso_cds))
ecoli_cds_lengths <- sapply(ecoli_cds, length)
meso_cds_lengths <- sapply(meso_cds, length)

boxplot(list(Ecoli = ecoli_cds_lengths, Mesomycoplasma = meso_cds_lengths),
        main = "CDS Length Comparison",
        ylab = "CDS length (bp)",
        col = c("skyblue","tomato"))


# Mean and Median
mean_ecoli <- mean(ecoli_cds_lengths)
median_ecoli <- median(ecoli_cds_lengths)
mean_meso <- mean(meso_cds_lengths)
median_meso <- median(meso_cds_lengths)

cat("E. coli: mean =", mean_ecoli, ", median =", median_ecoli, "\n")
cat("Mesomycoplasma: mean =", mean_meso, ", median =", median_meso, "\n")

```
```{r}
# Nucleotide frequency:
ecoli_concat <- toupper(paste(sapply(ecoli_cds, function(x) paste(getSequence(x), collapse = "")), collapse = ""))
meso_concat <- toupper(paste(sapply(meso_cds, function(x) paste(getSequence(x), collapse = "")), collapse = ""))

ecoli_nt_freq <- table(strsplit(ecoli_concat, "")[[1]])
meso_nt_freq <- table(strsplit(meso_concat, "")[[1]])

par(mfrow=c(1,2))
barplot(ecoli_nt_freq, main="E. coli Nucleotide Frequency", col="skyblue")
barplot(meso_nt_freq, main="Mesomycoplasma Nucleotide Frequency", col="tomato")

# Amino acid frequency (protein translation)
library(Biostrings)
ecoli_AA <- lapply(ecoli_cds, function(x) as.character(translate(DNAString(paste(getSequence(x), collapse = "")))))
ecoli_AA_concat <- paste(unlist(ecoli_AA), collapse = "")
ecoli_aa_freq <- table(strsplit(ecoli_AA_concat, "")[[1]])

Mesomycoplasma_AA <- lapply(meso_cds, function(x) as.character(translate(DNAString(paste(getSequence(x), collapse = "")))))
Mesomycoplasma_AA_concat <- paste(unlist(Mesomycoplasma_AA), collapse = "")
Mesomycoplasma_aa_freq <- table(strsplit(Mesomycoplasma_AA_concat, "")[[1]])

par(mfrow=c(1,2))
barplot(ecoli_aa_freq, main="E. coli Amino Acid Frequency", col="skyblue")
barplot(Mesomycoplasma_aa_freq, main="Mesomycoplasma Amino Acid Frequency", col="tomato")

```
```{r}
# Codon frequency function
codon_table <- function(seq_list) {
  codons <- unlist(lapply(seq_list, function(x) {
    s <- toupper(paste(getSequence(x), collapse = ""))
    groups <- substring(s, seq(1,nchar(s)-2,3), seq(3, nchar(s), 3))
    groups
  }))
  table(codons)
}

ecoli_codon_freq <- codon_table(ecoli_cds)
meso_cds <- read.fasta(file = "Mesomycoplasma_cds.fa", seqtype = "DNA")
meso_cds_lengths <- sapply(meso_cds, length)
meso_codon_freq <- codon_table(meso_cds)

# Convert to proportion:
ecoli_codon_prop <- ecoli_codon_freq / sum(ecoli_codon_freq)
meso_codon_prop <- meso_codon_freq / sum(meso_codon_freq)

# Compare using barplots
par(mfrow=c(2,1))
barplot(sort(ecoli_codon_prop), las=2, main="E. coli Codon Usage (%)", col="skyblue")
barplot(sort(meso_codon_prop), las=2, main="Mesomycoplasma Codon Usage (%)", col="tomato")

```
```{r}
# K-mer counting function
kmer_counter <- function(seq_vec, k=3) {
  kmers <- unlist(lapply(seq_vec, function(s) {
    n <- nchar(s)
    if (n < k) return(character(0))
    sapply(1:(n-k+1), function(i) substr(s, i, i+k-1))
  }))
  table(kmers)
}
# For E. coli
ecoli_aa_strings <- sapply(ecoli_AA, paste, collapse = "")
k=3
ecoli_kmer <- kmer_counter(ecoli_aa_strings, k)
Mesomycoplasma_aa_strings <- sapply(Mesomycoplasma_AA, paste, collapse = "")
Mesomycoplasma_kmer <- kmer_counter(Mesomycoplasma_aa_strings, k)

ecoli_kmer_prop <- ecoli_kmer / sum(ecoli_kmer)
Mesomycoplasma_kmer_prop <- Mesomycoplasma_kmer / sum(Mesomycoplasma_kmer)

# Most over- and under-represented in Mesomycoplasma:
Mesomycoplasma_kmer_sorted <- sort(Mesomycoplasma_kmer_prop, decreasing=TRUE)
over_Mesomycoplasma <- head(Mesomycoplasma_kmer_sorted, 10)
under_Mesomycoplasma <- tail(Mesomycoplasma_kmer_sorted, 10)

print(over_Mesomycoplasma)
print(under_Mesomycoplasma)

# Compare to E. coli for the same k-mers:
ecoli_over_match <- ecoli_kmer_prop[names(over_Mesomycoplasma)]
ecoli_under_match <- ecoli_kmer_prop[names(under_Mesomycoplasma)]

# Plot comparisons:
barplot(rbind(as.numeric(over_Mesomycoplasma), as.numeric(ecoli_over_match)),
        beside=TRUE, names.arg=names(over_Mesomycoplasma),
        legend.text=c("Mesomycoplasma", "E. coli"),
        main="Most over-represented 3-mers", col=c("tomato","skyblue"))

barplot(rbind(as.numeric(under_Mesomycoplasma), as.numeric(ecoli_under_match)),
        beside=TRUE, names.arg=names(under_Mesomycoplasma),
        legend.text=c("Mesomycoplasma", "E. coli"),
        main="Most under-represented 3-mers", col=c("tomato","skyblue"))

```
