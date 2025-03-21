# Script for processing the traits data from Ronquist et al 2020

# Read in the data
D  <- read.delim("ronquist_2020_SE_traits.csv",sep=";")

# Remove empty rows
D <- D[!is.na(D$Sorting.Number),]

# Split combined taxon names into family and subfamily names
x <- character()
for (i in 1:nrow(D)) { word <- D$Combined.Taxon[i]; if (grepl("-",word)) { x[i] <- strsplit(word,"-")[[1]][1] } else x[i] <- word; }
D$Family <- x

for (i in 1:nrow(D)) { word <- D$Combined:Taxon[i]; if (grepl("-",word)) { x[i] <- strsplit(word,"-")[[1]][2] } else x[i] <- ""; }
D$Subfamily <- x

# Match against NCBI classification

