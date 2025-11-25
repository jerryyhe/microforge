# metagenomics_functions.R

# FUNCTIONS

#' Split the SQM Output Taxonomy Column into Individual Columns
#'
#'
#'
#' @param data An input dataframe containing a columne called `Taxon` for which the function will split.
#' @param taxa_list A string vector of taxa names that will be used to rename columns.
#' @return The input dataframe with the taxon column separated into the values of the `taxa_list`.
#' @export
split_taxonomy <- function(data, taxa_list = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";\\w_", taxa_col = "Taxon") {

  taxonomy <- taxa_list
  col <- data[[taxa_col]]

  taxon_df <- data %>%
    mutate(col = gsub("n_[^;]+;", "", col)) %>%
    separate(col, into = taxonomy, sep = sep)

  return(taxon_df)
}

#' Calculate Relative Abundance from
#'
#'
#'
#' @param data An input dataframe containing columns titles by taxon (allowed values are: Kingdom, Phylum, Class, Order, Family, Genus, Species). The dataframe must also contain columns for each sample suffixed with "_bases" that contain coverage metrics to calculate relative abundance with.
#' @param taxa_rank A string indicating the taxonomic level to calculated. Allowed values are "k", "p", "c", "o", "f", "g", "s", for the first letter of each taxonomic rank.
#' @return A dataframe with each taxa and the associated relative abundance.
#' @export
calculate_relative_abundance <- function(data, taxa_rank = "f") {

  rank_map <- c(
    k = "Kingdom",
    p = "Phylum",
    c = "Class",
    o = "Order",
    f = "Family",
    g = "Genus",
    s = "Species"
  )

  grouping_taxon <- rank_map[[taxa_rank]]

  rel_abundance <- data %>%
    select(Rank, Taxon, dplyr::ends_with("bases")) %>%
    split_taxonomy() %>%
    filter(Rank == taxa_rank) %>%
    filter(Kingdom == "k_Bacteria" | Kingdom == "Unknown") %>%
    pivot_longer(cols = dplyr::ends_with("bases"),
                 names_to = "sample",
                 values_to = "bases") %>%
    mutate(sample = gsub(" bases", "", sample)) %>%
    group_by(sample) %>%
    mutate(total_bases = sum(bases, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(relative_abundance = round(100 * (bases/total_bases)), 4) %>%
    arrange(desc(relative_abundance)) %>%
    mutate(Taxon = fct_reorder(factor(.data[[grouping_taxon]]), relative_abundance, .na_rm = TRUE))

  return(rel_abundance)
}
