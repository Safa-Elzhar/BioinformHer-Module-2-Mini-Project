
#Step 1: Install and Load Required R Packages
# Install from CRAN/Bioconductor
install.packages("rentrez")
# Load the libraries
library(rentrez)
library(Biostrings)


# Step 2: Retrieve Human HBB Gene Sequence from NCBI
# We can search for the HBB gene using its gene ID (3043) or accession number NM_000518 (mRNA) or NP_000509 (protein).
## Retrieve the nucleotide (mRNA) sequence
hbb_nuc <- entrez_fetch(db = "nuccore", id = "NM_000518", rettype = "fasta")
# print sequence
cat(hbb_nuc)
## Retrieve the protein sequence
# Human
hbb_human <- entrez_fetch(db = "protein", id = "NP_000509.1", rettype = "fasta")

# Chimpanzee
hbb_chimp <- entrez_fetch(db = "protein", id = "XP_508242.1", rettype = "fasta")

# Cow
hbb_cow <- entrez_fetch(db = "protein", id = "NP_776342.1", rettype = "fasta")

# Mouse
hbb_mouse <- entrez_fetch(db = "protein", id = "ADD52650.1", rettype = "fasta")

# Chicken
hbb_chicken <- entrez_fetch(db = "protein", id = "AAA48996.1", rettype = "fasta")

# Zebrafish
hbb_zebrafish <- entrez_fetch(db = "protein", id = "NP_001013045.1", rettype = "fasta")

# print sequence
cat(hbb_human)
## Step 3: Save the Sequence to a FASTA File 
writeLines(hbb_nuc, "HBB_human_mRNA.fasta")
writeLines(hbb_human, "HBB_human.fasta")
writeLines(hbb_chimp, "HBB_chimp_protein.fasta")
writeLines(hbb_cow, "HBB_cow_protein.fasta")
writeLines(hbb_mouse, "HBB_mouse_protein.fasta")
writeLines(hbb_chicken, "HBB_chicken_protein.fasta")
writeLines(hbb_zebrafish, "HBB_zebrafish_protein.fasta")

##  Step 4: Read the sequences 
human <- readAAStringSet("HBB_human.fasta")
chimp <- readAAStringSet("HBB_chimp_protein.fasta")
zfish <- readAAStringSet("HBB_zebrafish_protein.fasta")

##  Step 5: Perform pairwise alignments

# Human vs Chimpanzee
alignment1 <- pairwiseAlignment(human[[1]], chimp[[1]], substitutionMatrix = "BLOSUM62", gapOpening = -10, gapExtension = -0.5)

# Human vs Zebrafish
alignment2 <- pairwiseAlignment(human[[1]], zfish[[1]], substitutionMatrix = "BLOSUM62", gapOpening = -10, gapExtension = -0.5)

##  Step 6: View results
# Show alignment summaries
alignment1
alignment2

# Show percent identity
pid1 <- pid(alignment1)
pid2 <- pid(alignment2)

cat("Human vs Chimpanzee % Identity:", pid1, "%\n")
cat("Human vs Zebrafish % Identity:", pid2, "%\n")
### Interpretation
#The differences in alignment scores and sequence similarities between the two alignments highlight 
#the varying degrees of evolutionary conservation of the HBB gene among different species.
#The high similarity between human and chimpanzee HBB sequences underscores their recent common ancestry,
#while the greater divergence observed in zebrafish reflects a more ancient evolutionary split.

BiocManager::install("msa") 
library(msa)

human <- readAAStringSet("HBB_human.fasta")
chimp <- readAAStringSet("HBB_chimp_protein.fasta")
cow <- readAAStringSet("HBB_cow_protein.fasta")
mouse <- readAAStringSet("HBB_mouse_protein.fasta")
chicken <- readAAStringSet("HBB_chicken_protein.fasta")
zebrafish <- readAAStringSet("HBB_zebrafish_protein.fasta")

sequences <- c(human, chimp, cow, mouse, chicken, zebrafish)
alignment <- msa(sequences, method = "ClustalOmega")
print(alignment, show = "alignment")
# Save the Alignment to a File
writeXStringSet(unmasked(alignment), filepath = "HBB_alignment.fasta")
#Visualize the Alignment as pdf
install.packages("tinytex")
tinytex::install_tinytex()
library(tinytex)

msaPrettyPrint(alignment,
               output = "pdf",
               showNames = "left",
               showLogo = "top",
               askForOverwrite = FALSE,
               verbose = FALSE,
               file = "HBB_alignment.pdf")


#Constructing a Phylogenetic Tree in R
install.packages("ape")
install.packages("phangorn")
install.packages("ggtree")    # For visualization (requires Bioconductor)
install.packages("ggplot2")   # For plotting

library(ape)
library(phangorn)
library(ggtree)
library(ggplot2)

## step 7 Convert Alignment for Phylogenetic Analysis
##Convert the alignment to a format suitable for distance calculation:


# Convert alignment to a matrix
alignment_matrix <- as.matrix(alignment)

# Convert to DNAbin object
alignment_dnabin <- as.DNAbin(alignment_matrix)

##Compute Distance Matrix
#Calculate the pairwise distances between sequences:
dist_matrix <- dist.dna(alignment_dnabin, model = "raw")

# Build the tree
tree <- nj(dist_matrix)

# Plot the tree
plot(tree, main = "Phylogenetic Tree of HBB Sequences", type = "phylogram")
##The phylogenetic tree based on hemoglobin beta (HBB) gene sequences reflects
#evolutionary relationships among vertebrates. 
# Human (NP\_000509.1) and chimpanzee (XP\_508242.1) sequences cluster closely, indicating their recent common ancestry and high genetic similarity.
#Mouse (ADD52650.1) and cow (NP\_776342.1), as fellow mammals, form a separate clade, sharing a more distant common ancestor with primates. 
#Chicken (AAA48996.1), representing birds, diverges earlier from the mammalian lineage.
# while zebrafish, a fish species, branches off even earlier, highlighting its more distant relationship to the other vertebrates.
#This arrangement underscores the progressive divergence of species from a common ancestor, consistent with evolutionary theory.

