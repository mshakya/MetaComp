#' @importFrom data.table fread
NULL

#' Extracts taxonomic information from the otu_table and ncbi taxonomy file to
#' to create phyloseq formatted tax_table. ncbi taxonomy file is provided with the package
#' taxonomy.tsv.gz
#'
#'
#' @param otu table produced by convert_to_otu_table
#'
#' @param taxonomic level (First letter must be caps). It should be the same taxonomic level
#' that its corresponding otu table is at.
#'
#' @return phyloseq formatted taxonomic table
#'
#' @export

convert_to_taxa_table <- function(otu_table, TAXON){

  # get the path of the taxonomy table concatenated to command to unzip it
  ncbi_taxa_filepath <- base::system.file("extdata", "taxonomy.tsv.gz", package="MetaComp")
  ncbi_taxa_gz_filepath <- base::paste('gunzip -c ', ncbi_taxa_filepath, sep="")

  # read in the taxonomy file.
  ncbi_taxa <- data.table::fread(ncbi_taxa_gz_filepath, header=T)
  base::colnames(ncbi_taxa) <- c("taxid", "dont_know", "parent_taxid", "LEVEL", "NAME")

  # get taxa name from otu table
  taxa_name <- base::row.names(otu_table)


  #create an empty taxa table
  taxa_table <- base::data.frame(superkingdom = character(),
                           phylum = character(),
                           class = character(),
                           order = character(),
                           family = character(),
                           genus = character(),
                           species = character(),
                           stringsAsFactors = F)



  # loop through name of taxa and build a taxa_table
  for (taxa in taxa_name){
    taxon_level <- base::as.character(subset(ncbi_taxa, NAME==taxa)$LEVEL[1])
    parent_taxID <- base::as.character(subset(ncbi_taxa, NAME==taxa)$parent_taxid[1])
    taxID <- base::as.character(subset(ncbi_taxa, NAME==taxa)$taxid[1])
    one_row <- base::data.frame(taxon_level = taxa)
    colnames(one_row) <- taxon_level
    if (parent_taxID != "1") {
      while (parent_taxID != "1") {
        # get taxon LEVEL of parent taxa
        taxon_level <- base::as.character(subset(ncbi_taxa, taxid==parent_taxID)$LEVEL[1])
        # get taxon NAME of parent taxa
        taxon_name <- base::as.character(subset(ncbi_taxa, taxid==parent_taxID)$NAME[1])
        # get taxid
        taxID <- base::as.character(subset(ncbi_taxa, taxid==parent_taxID)$taxid[1])
        # add that to the one_row data frame
        one_row[taxon_level] <- taxon_name
        # get parent_taxID and loop through until parent_taxID doesnt equal to 131567 (highest classification)
        parent_taxID <- base::as.character(subset(ncbi_taxa, taxid==taxID)$parent_taxid[1])
        taxID <- base::as.character(subset(ncbi_taxa, taxid==taxid)$taxid[1])
      }
      taxa_table <- dplyr::bind_rows(taxa_table, one_row)
    }
    else {
      taxa_table <- dplyr::bind_rows(taxa_table, one_row)
    }


  }

  base::colnames(taxa_table) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  base::subset(taxa_table, select=c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
  base::rownames(taxa_table) <- taxa_table[, TAXON]
  phyloseq::tax_table(base::as.matrix(taxa_table))
}
