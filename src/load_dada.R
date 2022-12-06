# Functions to load DADA2-output (based on the standard output)

########################
## `create_tax_names` ##
########################

create_tax_names <- function(tax_table) {
  tax_names <- apply(tax_table, 1, function(x) {
    if(any(names(na.omit(x)) %in% "species")) {
      name <- str_c(tail(na.omit(x), 2), collapse = "_")
    } else { 
      name <- tail(na.omit(x), 1)
    }
    return(name)
  })
  
  tax_names %<>% 
    str_replace_all("-| ", "_") %>%
    str_remove_all("\\[|\\]|\\(|\\)")
  
  return(str_c(tax_names, 1:nrow(tax_table), sep = "_"))
}

##########################
## `create_taxa_filter` ##
##########################

create_taxa_filter <- function(taxa) {
  #remove Chloroplast/Mitochondria/Archae/Eukaryota/no kingdom-level annotation
  #see Park et al. (Raes). Microbiome. 2019. / Winkler et al. Cell. 2020.
  taxa_filter <- (str_detect(replace_na(taxa[, "family"], ""), "Mitochondria") | 
                    str_detect(replace_na(taxa[, "order"], ""), "Chloroplast") |
                    str_detect(replace_na(taxa[, "domain"], ""), "Archae|Eukaryota") |
                    is.na(taxa[, "domain"]))
  return(taxa_filter)
}

####################
## `create_phylo` ##
####################

create_phylo <- function(seqtab, taxa, meta) {
  ps <- phyloseq(otu_table(t(seqtab), taxa_are_rows=T), 
                 sample_data(meta), # TODO add meta
                 tax_table(taxa))
  dna <- Biostrings::DNAStringSet(taxa_names(ps))
  names(dna) <- taxa_names(ps)
  ps <- merge_phyloseq(ps, dna)
  taxa_names(ps) <- create_tax_names(taxa)
  tax_table(ps) <- cbind(tax_table(ps), ASV = taxa_names(ps))
  
  ps_filtered <- ps %>% 
    prune_taxa(!create_taxa_filter(taxa), .) %>%
    prune_samples(!sample_sums(.) == 0, .)
  
  n_taxa_excl <- ntaxa(ps) - ntaxa(ps_filtered)
  
  print(glue::glue("ASVs beloning to the phylogenetic groups 'Mitochondria', 'Chloroplast', 'Archae' or 'Eukaryota' were
                   excluded (n = {n_taxa_excl} ASVs), resulting in a total of {scales::comma(nrow(taxa))} ASVs after exclusion."))
  
  return(ps)
}

#fix the function: final return object as well as the printing line to say how many filtered out and how many remaining.
create_phylo_fixed <- function(seqtab, taxa, meta) {
  ps <- phyloseq(otu_table(t(seqtab), taxa_are_rows=T), 
                 sample_data(meta), # TODO add meta
                 tax_table(taxa))
  dna <- Biostrings::DNAStringSet(taxa_names(ps))
  names(dna) <- taxa_names(ps)
  ps <- merge_phyloseq(ps, dna)
  taxa_names(ps) <- create_tax_names(taxa)
  tax_table(ps) <- cbind(tax_table(ps), ASV = taxa_names(ps))
  
  ps_filtered <- ps %>% 
    prune_taxa(!create_taxa_filter(taxa), .) #%>%
    #prune_samples(!sample_sums(.) == 0, .)
  
  n_taxa_excl <- ntaxa(ps) - ntaxa(ps_filtered)
  
  print(glue::glue("ASVs beloning to the phylogenetic groups 'Mitochondria', 'Chloroplast', 'Archae' or 'Eukaryota' were
                   excluded (n = {n_taxa_excl} ASVs), resulting in a total of {scales::comma(ntaxa(ps_filtered))} ASVs after exclusion."))
  
  return(ps_filtered)
}
