#Phyloseq Object Creation

#### load libraries ####
library(phyloseq)
library(ape) 
library(tidyverse)
library(vegan)

#### Load data #### 
metafp <- "team9_project/qimme2_exported/parkinsons_metadata.txt"
meta <- read_delim(metafp, delim="\t")

otufp <- "team9_project/qimme2_exported/feature-table.txt"
otu <- read_delim(file = otufp, delim="\t", skip=1)

taxfp <- "team9_project/qimme2_exported/taxonomy.tsv"
tax <- read_delim(taxfp, delim="\t")

phylotreefp <- "team9_project/qimme2_exported/tree.nwk"
phylotree <- read.tree(phylotreefp)

#### Format OTU table ####
#save everything except first column (OTU ID) into a matrix
otu_mat <- as.matrix(otu[,-1])

#Make first column (#OTU ID) the rownames of the new matrix
rownames(otu_mat) <- otu$`#OTU ID`

#Use the "otu_table" function to make an OTU table
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE) 

#### Format sample metadata ####
#Save everything except sampleid as new data frame
samp_df <- as.data.frame(meta[,-1])

#Make sampleids the rownames
rownames(samp_df) <- meta[["#SampleID"]]

#Make phyloseq sample data with sample_data() function
SAMP <- sample_data(samp_df)

#### Formatting taxonomy ####
#Convert taxon strings to a table with separate taxa rank columns
tax_mat <- tax %>% 
  separate(col=Taxon, sep="; "
           , into = c("Domain","Phylum","Class","Order","Family","Genus","Species")) %>%
  as.matrix() # Saving as a matrix

#Save everything except feature IDs 
tax_mat <- tax_mat[,-1]

# Make sampleids the rownames
rownames(tax_mat) <- tax$`Feature ID`

# Make taxa table
TAX <- tax_table(tax_mat)

#### Create phyloseq object ####
# Merge all into a phyloseq object
parkinsons_phyloseq <- phyloseq(OTU, SAMP, TAX, phylotree)

#save phyloseq object
save(parkinsons_phyloseq, file = "parkinsons_phyloseq.RData")