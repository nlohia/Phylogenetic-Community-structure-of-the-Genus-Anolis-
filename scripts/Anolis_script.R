# Phylogenetic Community Structure of Anolis Genus
# Author: Nandini Lohia
# Purpose: Comprehensive analysis of Anolis species distribution and phylogenetic relationships
#Research Questions:
#1.What is the geographical distribution of Anolis species?
#2.What are the evolutionary relationships among Anolis species?
#3.How are Anolis species distributed across different communities?

# 1. PACKAGE INSTALLATION AND DEPENDENCY MANAGEMENT
# Install and load necessary packages for data manipulation, biological analysis, phylogenetic analysis, spatial analysis, and visualization
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  # Data manipulation and visualization
  tidyverse, dplyr, tidyr, 
  
  # Biological sequence analysis
  Biostrings, DECIPHER, seqinr, rentrez, msa, muscle,
  
  # Phylogenetic analysis
  ape, phangorn, treeio, picante,
  
  # Spatial and ecological analysis
  rgbif, sf, vegan, 
  
  # Statistical and computational tools
  cluster, factoextra, future, furrr, parallel,
  
  # Advanced visualization
  plotly, RColorBrewer, ggtree, tidytree, ggplot2
)

# 2. SET RANDOM SEEDS FOR REPRODUCIBILITY
# Setting random seeds for reproducibility of analysis at different stages
set.seed(123)  # For community analysis
set.seed(42)   # For phylogenetic analysis

# 3. OCCURRENCE DATA RETRIEVAL
# Research Question: What is the geographical distribution of Anolis species?
# Retrieve occurrence data from Global Biodiversity Information Facility (GBIF)
# Filtering and cleaning the data for relevant columns: species, latitude, longitude, country, state, and IUCN status
occ_search_results <- occ_search(scientificName = "Anolis", hasCoordinate = TRUE, limit = 500)
occurrence_data <- occ_search_results$data %>%
  filter(!is.na(decimalLatitude), !is.na(decimalLongitude), !is.na(scientificName)) %>%
  distinct(species, decimalLatitude, decimalLongitude, .keep_all = TRUE) %>%
  dplyr::select(species, decimalLatitude, decimalLongitude, country, stateProvince, iucnRedListCategory)

# View the cleaned occurrence data
head(occurrence_data)

# 4. SEQUENCE DATA RETRIEVAL AND PROCESSING
# Research Questions: Quality control of genetic sequences from GenBank
# Example commented-out code to search for specific gene sequences (e.g., 16S rRNA) for Anolis species
#Anolis_16S_search <- entrez_search(db = "nuccore", term = "Anolis AND 16S", retmax = 300)
#Anolis_16S_search1 <- entrez_search(db = "nuccore", term = "Anolis AND 16S AND 400:650[SLEN]", retmax = 300)

# Fetch sequences (this part is commented out as it's hypothetical)
#Anolis_16S_sequences <- entrez_fetch(db = "nuccore", id = Anolis_16S_search1$ids, rettype = "fasta", retmode = "text")

# Load the sequences into R
Anolis_16S_summ <- readDNAStringSet("./data/Anolis_16S_sequences.fasta")

# Quality Control and Filtering of Sequences
# Remove sequences with 'N' (missing bases), sequences with excessive gaps (>50%), and duplicates
Anolis_16S_summ_cleaned <- Anolis_16S_summ[!grepl("N", Anolis_16S_summ)]  # Remove 'N'
max_gap_percentage <- 0.5  # Allow only 50% gaps
Anolis_16S_summ_cleaned <- Anolis_16S_summ_cleaned[sapply(Anolis_16S_summ_cleaned, function(seq) {
  sum(strsplit(as.character(seq), NULL)[[1]] == "-") / length(seq) < max_gap_percentage
})]

# Remove duplicate sequences and sequences shorter than 100 bases
Anolis_16S_summ_cleaned <- unique(Anolis_16S_summ_cleaned)
min_length <- 100
Anolis_16S_summ_cleaned <- Anolis_16S_summ_cleaned[nchar(Anolis_16S_summ_cleaned) > min_length]

# Create a dataframe for cleaned sequences
dfAnolis_16S <- data.frame(Anolis16S_Title = names(Anolis_16S_summ_cleaned), 
                           Anolis_16S_Sequence = paste(Anolis_16S_summ_cleaned))
dfAnolis_16S$Species_Name <- word(dfAnolis_16S$Anolis16S_Title, 2L, 3L)

# Clean species names and filter valid species
dfAnolis_16S <- dfAnolis_16S %>%
  mutate(Species_Name = gsub("^A\\.\\s*", "Anolis ", Species_Name)) %>%
  mutate(Species_Name = word(Species_Name, 1, 2)) %>%
  filter(grepl("^Anolis\\s\\w+", Species_Name)) %>%
  select(Anolis16S_Title, Species_Name, Anolis_16S_Sequence)

# View the cleaned sequence data
head(dfAnolis_16S)

# 5. PHYLOGENETIC TREE CONSTRUCTION
# Research Questions: Evolutionary relationships of Anolis species
# Create the DNAStringSet object for the cleaned sequence data
dna_sequences <- DNAStringSet(dfAnolis_16S$Anolis_16S_Sequence)
names(dna_sequences) <- dfAnolis_16S$Species_Name

# Remove duplicate species names
dfAnolis_16S_unique <- dfAnolis_16S[!duplicated(dfAnolis_16S$Species_Name), ]
dna_sequences_unique <- DNAStringSet(dfAnolis_16S_unique$Anolis_16S_Sequence)
names(dna_sequences_unique) <- dfAnolis_16S_unique$Species_Name

# Perform multiple sequence alignment using the MUSCLE algorithm
aligned_sequences_unique <- muscle::muscle(dna_sequences_unique)

# Convert aligned sequences to phyDat format for phylogenetic analysis
aligned_phyDat_unique <- as.phyDat(aligned_sequences_unique)

# Compute a distance matrix for phylogenetic tree construction
dist_matrix <- dist.ml(aligned_phyDat_unique, model = "JC69")

# Build a Neighbor-Joining tree from the distance matrix
nj_tree <- nj(dist_matrix)

# Plot the NJ tree
plot(nj_tree, main = "Neighbor-Joining Phylogenetic Tree of Anolis species", 
     cex = 0.4, edge.width = 0.3, no.margin = TRUE)

# Maximum Likelihood tree optimization
ml_tree <- pml(nj_tree, aligned_phyDat_unique)
ml_tree_optimized <- optim.pml(ml_tree, model = "GTR", optInv = TRUE, optGamma = TRUE)

# Print summary of the optimized Maximum Likelihood tree
summary(ml_tree_optimized)

# Plot the optimized tree with adjusted graphical parameters
par(mar = c(1, 1, 1, 1))  # Adjust margins
plot(ml_tree_optimized$tree, cex = 0.4)  # Scale tree label size

# 6. COMMUNITY STRUCTURE ANALYSIS
# Research Questions: Species distribution across communities
# Find species common to both occurrence and sequence data
species_occurrence <- unique(occurrence_data$species)
species_sequences <- unique(dfAnolis_16S_unique$Species_Name)
common_species <- intersect(species_occurrence, species_sequences)

# Filter occurrence and sequence data to only include common species
occurrence_common <- occurrence_data %>% 
  filter(species %in% common_species)

sequence_common <- dfAnolis_16S_unique %>% 
  filter(Species_Name %in% common_species)

# Check results for common species
print(common_species)
head(occurrence_common)
head(sequence_common)

# Visualize species occurrence on a map using Plotly
unique_countries <- unique(occurrence_common$country)
colors <- rainbow(length(unique_countries))
country_colors <- setNames(colors, unique_countries)
occurrence_common$color <- country_colors[occurrence_common$country]

# Create the interactive map
plot_ly(occurrence_common, 
        type = 'scattermapbox',  
        lat = ~decimalLatitude,  
        lon = ~decimalLongitude,  
        mode = 'markers',  
        marker = list(color = ~color, size = 8, opacity = 0.7),
        text = ~paste0("Species: ", species, "<br>",  
                       "Country: ", country, "<br>",  
                       "State/Province: ", stateProvince, "<br>",  
                       "Latitude: ", decimalLatitude, "<br>",  
                       "Longitude: ", decimalLongitude),  
        hoverinfo = 'text') %>%  
  layout(mapbox = list(style = 'open-street-map', 
                       zoom = 3, 
                       center = list(lat = 0, lon = -90), 
                       bearing = 0, 
                       pitch = 0),  
         margin = list(l = 0, r = 0, t = 0, b = 0))  # Remove margins

# 6. COMMUNITY STRUCTURE ANALYSIS
# Research Questions: Species distribution across communities
# Find common species between occurrence and sequence data
# Combine occurrence data with phylogenetic information

# Step 1: Create the community matrix
community_matrix <- occurrence_common %>%
  group_by(country, species) %>%
  summarise(n = n(), .groups = "drop") %>%  # Summarize the number of occurrences per species per country
  pivot_wider(names_from = species, values_from = n, values_fill = 0) %>% # Pivot data into a matrix format
  column_to_rownames(var = "country")  # Set country names as row names

# Check the resulting community matrix
head(community_matrix)  # View the first few rows of the matrix

# Step 2: Prepare data for phylogenetic community structure
# Combine occurrence data with phylogenetic information, ensuring tree tip labels match community matrix species names
pruned_tree <- keep.tip(ml_tree_optimized$tree, colnames(community_matrix))  # Prune tree to match species in the community matrix

# Step 3: Calculate phylogenetic community structure metrics
phylo_community_structure <- ses.mpd(community_matrix, cophenetic(pruned_tree))  # Calculate mean pairwise distance

# Remove rows with missing or invalid 'mpd.obs' values (e.g., negative distances)
phylo_community_structure_clean <- phylo_community_structure %>%
  filter(!is.na(mpd.obs) & mpd.obs >= 0)  # Filter out invalid data

# Step 4: Visualize the phylogenetic community structure with the clean data
ggplot(phylo_community_structure_clean, aes(x = reorder(row.names(phylo_community_structure_clean), mpd.obs), y = mpd.obs)) +
  geom_bar(stat = "identity", fill = "skyblue") +  # Create bar plot for mean pairwise distance
  labs(title = "Mean Pairwise Phylogenetic Distance by Country",
       x = "Country", y = "MPD Observed") +  # Label axes
  theme_minimal() +  # Apply minimal theme
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels for readability
  geom_text(aes(label = round(mpd.obs, 2)), vjust = -0.3, size = 3)  # Display MPD values on top of bars

# Step 5: Midpoint rooting of the phylogenetic tree
rooted_tree <- midpoint(pruned_tree)  # Root tree at its midpoint

# Verify the rooted tree (ensure it's properly rooted)
is.rooted(rooted_tree)

# Step 6: Calculate Phylogenetic Diversity (PD)
pd_results <- pd(community_matrix, rooted_tree)  # Calculate phylogenetic diversity using the community matrix

# Print results of phylogenetic diversity
print(pd_results)

# Step 7: Visualize Phylogenetic Diversity by Country
# Prepare data for plotting
pd_plot_data <- data.frame(
  Country = rownames(pd_results),
  PD = pd_results$PD,  # Phylogenetic diversity
  SR = pd_results$SR   # Species richness
)

# Create a bar plot to visualize phylogenetic diversity
ggplot(pd_plot_data, aes(x = Country, y = PD)) +
  geom_bar(stat = "identity", fill = "forestgreen") +  # Bar plot for PD values
  labs(title = "Phylogenetic Diversity by Country",
       x = "Country", 
       y = "Phylogenetic Diversity") +
  theme_minimal() +  # Minimal theme for readability
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels
  geom_text(aes(label = round(PD, 2)), vjust = -0.3, size = 3)  # Display PD values on top of bars

# Step 8: Convert tree to treedata object for annotation
tree_data <- as_tibble(rooted_tree)

# Create community annotation data: assign each species to its most frequent community
community_annotations <- data.frame(
  species = colnames(community_matrix),
  community = apply(community_matrix, 2, function(x) names(which.max(x)))  # Assign each species to the community with the highest frequency
)

# Create a color palette for the communities
n_communities <- length(unique(community_annotations$community))  # Get number of unique communities
community_colors <- setNames(
  brewer.pal(n = max(3, min(n_communities, 9)), "Set3"),  # Choose a color palette with a suitable number of colors
  unique(community_annotations$community)  # Map each community to a color
)

# Step 9: Plot the enhanced phylogenetic tree with community distribution
p <- ggtree(rooted_tree) %<+% community_annotations +
  geom_tippoint(aes(color = community), size = 3) +  # Plot tips with community colors
  scale_color_manual(values = community_colors) +  # Apply custom colors
  theme_tree2() +  # Enhance tree visualization
  ggplot2::labs(title = "Phylogenetic Tree with Community Distribution",
                color = "Community")  # Label the plot

# Add a scale bar to the tree plot
p <- p + geom_treescale()
plot(p)  # Display the plot

# Step 10: Compute ses.mpd with null model and runs for significance
ses_mpd_results <- ses.mpd(community_matrix, 
                           cophenetic(rooted_tree),
                           null.model = "taxa.labels",  # Null model based on species labels
                           abundance.weighted = FALSE,  # Do not weight by abundance
                           runs = 999)  # Number of randomizations

# Clean up and filter results (remove missing values)
df_clean <- na.omit(ses_mpd_results)

# Create a data frame for analysis with observed MPD, Z-score, and p-value
mpd_clean <- data.frame(
  Community = rownames(df_clean),
  MPD_Observed = df_clean$mpd.obs,  # Observed MPD
  MPD_Z_Score = df_clean$mpd.obs.z,  # Z-score for observed MPD
  MPD_P_Value = df_clean$mpd.obs.p  # p-value for significance
)

# Interpretation based on Z-scores: Phylogenetic clustering, overdispersion, or random structure
mpd_clean$Interpretation <- case_when(
  mpd_clean$MPD_P_Value < 0.05 & mpd_clean$MPD_Z_Score < -1.96 ~ "Phylogenetic Clustering",  # Clustered if p < 0.05 and Z < -1.96
  mpd_clean$MPD_P_Value < 0.05 & mpd_clean$MPD_Z_Score > 1.96 ~ "Phylogenetic Overdispersion",  # Overdispersed if p < 0.05 and Z > 1.96
  TRUE ~ "Random Phylogenetic Structure"  # Random if p >= 0.05
)

# Step 11: Plotting the Z-scores for interpretation of phylogenetic community structure
ggplot(mpd_clean, aes(x = reorder(Community, MPD_Z_Score), y = MPD_Z_Score, fill = Interpretation)) +
  geom_bar(stat = "identity") +  # Bar plot for Z-scores
  scale_fill_manual(values = c("blue", "gray", "gray")) +  # Custom fill colors for each interpretation
  labs(title = "Phylogenetic Community Structure",
       x = "Community", 
       y = "MPD Standardized Effect Size") +  # Axis labels
  theme_minimal() +  # Minimal theme for clarity
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels
  geom_hline(yintercept = 0, linetype = "dashed") +  # Add a dashed line at y=0
  geom_hline(yintercept = c(-1.96, 1.96), linetype = "dotted", color = "red")  # Add red lines for Z-scores ±1.96

# Step 12: Additional insights: Calculate summary statistics grouped by interpretation
summary_stats <- mpd_clean %>%
  group_by(Interpretation) %>%
  summarise(
    Count = n(),  # Number of communities per interpretation
    Mean_Z_Score = mean(MPD_Z_Score),  # Mean Z-score per interpretation
    Mean_P_Value = mean(MPD_P_Value)  # Mean p-value per interpretation
  )

# Print summary statistics
cat("Summary of Phylogenetic Community Structure:\n")
print(summary_stats)  # Display the summary statistics for interpretation