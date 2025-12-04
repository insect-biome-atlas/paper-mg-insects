# Analyze results from lca info in blastn searches of expanded chesters sequence dataset
source("lca_fxns.R")

# Read in blastn results
B <- read.delim("blastn.tsv")

# Read in taxonomic data for original expanded sequence dataset
C <- read.delim("chesters_expanded_taxonomy_raw.tsv")

# Read in lca results
D <- read.delim("blastn_lca_100_hits.tsv")
D1 <- read.delim("blastn_lca_0_95_idty.tsv")

# Read in NCBI taxonomy data
T <- read.delim("../ncbi/ncbi_name_ranked_classification.tsv")

# Define function for updating taxonomy
update_taxonomy <- function(D) {

    # Adjust some taxonomy mismatches
    # using updated taxonomy

    # Neuropteroid orders grouped in expanded chesters
    Neuropterida <- c("Neuroptera","Megaloptera","Raphidioptera")
    D$order[D$order %in% Neuropterida] <- "Neuropterida"

    # Psocoptera and Phthiraptera now grouped in Psocodea
    Psocodea <- c("Psocoptera","Phthiraptera")
    D$order[D$order %in% Psocodea] <- "Psocodea"

    # Reclinomonas classification updated
    idx <- which(D$genus=="Reclinomonas")
    D$kingdom[idx] <- "Diphoda"
    D$phylum[idx] <- "Discoba"
    D$class[idx] <- "Jakobea"

    # Aphanomyces classification updated
    idx <- which(D$genus=="Aphanomyces")
    D$class[idx] <- "Peronosporea"

    # Peronospora classification updated
    idx <- which(D$genus=="Peronospora")
    D$class[idx] <- "Peronosporea"

    # Diplura classification updated
    idx <- which(D$family %in% c("Campodeidae","Projapygidae","Octostigmatidae","Anajapygidae"))
    D$order[idx] <- "Rhabdura"
    idx <- which(D$family %in% c("Japygidae","Heterojapygidae","Parajapygidae"))
    D$order[idx] <- "Dicellurata"
    idx <- which(D$order %in% c("Dicellurata","Rhabdura"))
    D$class[idx] <- "Diplura"

    # Protura classification updated
    idx <- which(D$family %in% c("Sinentomidae","Fujientomidae"))
    D$order[idx] <- "Sinentomata"
    idx <- which(D$family %in% c("Eosentomidae"))
    D$order[idx] <- "Eosentomata"
    idx <- which(D$family %in% c("Hesperentomidae","Acerentomidae","Protentomidae"))
    D$order[idx] <- "Acerentomata"

    # Get annotation from expanded chesters data
    idx <- match(D$query, C$TipLabel)
    D$Chesters_phylum <- C$Phylum[idx]
    D$Chesters_class <- C$Class[idx]
    D$Chesters_order <- C$Order[idx]
    D$Chesters_family <- C$Family[idx]
    D$Chesters_genus <- C$Genus[idx]
    D$Chesters_species <- C$Species[idx]

    D
}

# Update the taxonomy
D <- update_taxonomy(D)
D1 <- update_taxonomy(D1)

# Identify class misses that are likely errors

# Cycle over the problematic tip labels (class misses)
# Tentatively identified as the tip labels with ambiguous or conflicting class hits
class_miss_stats <- data.frame()
W <- D[D$Chesters_class != D$class & D$class!="None",]
cat("Checking for class level misses\n")
for (i in 1:nrow(W)) {

    next

    # For each of them, sort out the query hits
    # and add taxonomy info
    tip_label <- W$query[i]
    cat("Processing tip_label: ",tip_label,"\n")
    B2 <- B[B$query==tip_label,]
    if (nrow(B2)==0)
        next
    B2 <- add_lca_info(B2, T)

    # Remove all "unclassified" and "" records at class level
    B2 <- B2[!B2$Class=="" & !grepl("unclassified",B2$Class),]

    # Remove the top hit to correct class and all identical accns
    chesters_class <- C$Class[match(tip_label,C$TipLabel)]
    idx <- match(chesters_class, B2$Class)
    if (!is.na(idx)) {
        accn <- B2$accession[idx[1]]
        B2 <- B2[!B2$accession==accn,]
    }

    # Identify as potential error if the class list is not clean
    # (or if it does not contain the right class)
    hit_class <- character()
    prop_corr <- 0.0
    n_hits <- 0
    hit_classes <- ""
    hit_freqs <- ""
    hit_class95 <- character()
    prop_corr95 <- 0.0
    n_hits95 <- 0
    hit_classes95 <- ""
    hit_freqs95 <- ""
    if (nrow(B2) > 0) {
        hit_class <- unique(B2$Class)
        a <- table(B2$Class)
        if (chesters_class %in% names(a))
            prop_corr <- as.numeric(a[chesters_class] / sum(a))
        hit_classes <- paste(names(a),collapse=";")
        hit_freqs <- paste(a,collapse=";")
        write.table(B2,paste0("class_miss_taxonomies/",tip_label,"_hit_taxonomy.tsv"), sep="\t", row.names=FALSE)
        B2 <- B2[B2$per.ident>=95.0,]
        if (nrow(B2) > 0) {
            hit_class95 <- unique(B2$Class)
            a <- table(B2$Class)
            if (chesters_class %in% names(a))
                prop_corr95 <- as.numeric(a[chesters_class] / sum(a))
            hit_orders95 <- paste(names(a),collapse=";")
            hit_freqs95 <- paste(a,collapse=";")
        }
    }
    class_miss_stats <- rbind(class_miss_stats,
                              list(TipLabel=tip_label,
                                   num_classes=length(hit_class),
                                   prop_corr=prop_corr,
                                   hit_classes=hit_classes,
                                   hit_freqs=hit_freqs,
                                   num_classes95=length(hit_class95),
                                   prop_corr95=prop_corr95,
                                   hit_classes95=hit_classes95,
                                   hit_freqs95=hit_freqs95))
}

write.table(class_miss_stats, "class_miss_stats.tsv", row.names=FALSE, sep="\t")


# Identify order misses that are likely errors

# Cycle over the problematic tip labels (order misses)
# Tentatively identified as the tip labels with ambiguous or conflicting class hits
order_miss_stats <- data.frame()
V <- D[D$Chesters_order != D$order & D$order!="None",]
cat("Checking for order level misses\n")
for (i in 1:nrow(V)) {

    # For each of them, sort out the query hits
    # and add taxonomy info
    tip_label <- V$query[i]
    cat("Processing tip_label: ",tip_label,"\n")
    B2 <- B[B$query==tip_label,]
    if (nrow(B2)==0)
        next
    B2 <- add_lca_info(B2, T)

    # Remove all "unclassified" and "" records at class level
    B2 <- B2[!B2$Order=="" & !grepl("unclassified",B2$Order),]

    # Remove the top hit to correct class and all identical accns
    chesters_order <- C$Order[match(tip_label,C$TipLabel)]
    idx <- match(chesters_order, B2$Order)
    if (!is.na(idx)) {
        accn <- B2$accession[idx[1]]
        B2 <- B2[!B2$accession==accn,]
    }

    # Identify as potential error if the order list is not clean
    # (or if it does not contain the right order)
    hit_order <- character()
    prop_corr <- 0.0
    n_hits <- 0
    hit_orders <- ""
    hit_freqs <- ""
    hit_order95 <- character()
    prop_corr95 <- 0.0
    n_hits95 <- 0
    hit_orders95 <- ""
    hit_freqs95 <- ""
    if (nrow(B2) > 0) {
        hit_order <- unique(B2$Order)
        a <- table(B2$Order)
        if (chesters_order %in% names(a))
            prop_corr <- as.numeric(a[chesters_order] / sum(a))
        hit_orders <- paste(names(a),collapse=";")
        hit_freqs <- paste(a,collapse=";")
        write.table(B2,paste0("order_miss_taxonomies/",tip_label,"_hit_taxonomy.tsv"), sep="\t", row.names=FALSE)
        B2 <- B2[B2$per.ident>=95.0,]
        if (nrow(B2) > 0) {
            hit_order95 <- unique(B2$Order)
            a <- table(B2$Order)
            if (chesters_order %in% names(a))
                prop_corr95 <- as.numeric(a[chesters_order] / sum(a))
            hit_orders95 <- paste(names(a),collapse=";")
            hit_freqs95 <- paste(a,collapse=";")
        }
    }

    order_miss_stats <- rbind(order_miss_stats, list(TipLabel=tip_label,
                                            num_orders=length(hit_order),
                                            prop_corr=prop_corr,
                                            hit_orders=hit_orders,
                                            hit_freqs=hit_freqs,
                                            num_orders95=length(hit_order95),
                                            prop_corr95=prop_corr95,
                                            hit_orders95=hit_orders95,
                                            hit_freqs95=hit_freqs95))

}

write.table(order_miss_stats, "order_miss_stats.tsv", row.names=FALSE, sep="\t")
