# Note, many of the functions below are implemented in the microbiomer-package.
# See here: https://github.com/wsteenhu/microbiomer

# Save/load function for 'heavy' objects
sl <- function(name, ..., overwrite = FALSE, dir_path = here::here("results", "RData", subdir_name)) {
  # Possibility to add name as name or literal character string
  name <- as.character(substitute(name))
  assign(name, 
    if(file.exists(glue::glue("{dir_path}/{name}.Rds")) && !overwrite) {
     readRDS(glue::glue("{dir_path}/{name}.Rds"))
    }
    else { 
     dir.create(dir_path, showWarnings = F, recursive = T)
     saveRDS(..., file=glue::glue("{dir_path}/{name}.Rds"))
     readRDS(glue::glue("{dir_path}/{name}.Rds"))
    }, envir=.GlobalEnv)
}

# Format printing numbers with comma separating thousands
f <- function(x, ...) {
  format(x = x, big.mark=",", ...)
}

# Round
r <- function(x, ...) {
  round(x = x, digits = 3, ...)
}

# Keep trailing zero's
# Significant figures: # http://stackoverflow.com/questions/3245862/
s <- function(vec, digits = 3, format = "f", flag = ""){
  return(gsub("\\.$", "", 
              formatC(vec, 
                      digits = digits, 
                      # use "fg" for significant digits
                      format = format, 
                      flag = flag)))
}


# Conversion functions

to_RA <- function(ps) {
  phyloseq::transform_sample_counts(ps, function(OTU) OTU / sum(OTU))
}

pres_abund_filter <- function(ps, pres = 2, abund = 0.001, verbose = TRUE) { # Subramanian filter
  is_raw <- max(phyloseq::otu_table(ps)) > 1 # detect is ps has raw reads
  
  if(is_raw) { # if raw reads; convert to RA for filtering
    ps_raw <- ps
    ps <- ps_raw %>% to_RA()
  }
  
  ps_filt <- phyloseq::filter_taxa(ps, function(x) sum(x > abund) >= pres, TRUE)
  
  if(verbose) {
    message(glue::glue("A total of {phyloseq::ntaxa(ps_filt)} ASVs were found to be present at or above a level of confident detection ({(abund * 100)}% relative abundance) in at least {pres} samples (n = {phyloseq::ntaxa(ps) - phyloseq::ntaxa(ps_filt)} ASVs excluded)."))
  }
  
  if(!is_raw) {
    return(ps_filt)
  } else {
    return(phyloseq::prune_taxa(phyloseq::taxa_names(ps_filt), ps_raw))
  }
}

mean_abund_filter <- function(ps, mean_abund = 0.001) {
  filter_taxa(ps, function(x) mean(x) > mean_abund, TRUE) } 

asinsqrt <- function(otu_table) {
  asin(sqrt(otu_table))
}

get_topn <- function(ps, n = 15, residuals = TRUE) {
  otu_tab <- phyloseq::otu_table(ps)
  otu_tab_n <- otu_tab[order(rowMeans(otu_tab), decreasing = TRUE)[1:n], ]
  otu_tab_res <- otu_tab[-order(rowMeans(otu_tab), decreasing = TRUE)[1:n], ]
  
  if(residuals) {
    otu_tab_n <- otu_tab_n %>%
      t %>%
      data.frame(residuals = colSums(otu_tab_res), check.names = FALSE) %>%
      t
  }
  return(phyloseq::phyloseq(phyloseq::otu_table(otu_tab_n, taxa_are_rows = TRUE),
                            phyloseq::sample_data(ps)))
}

prep_bar <- function(ps, n, residuals = TRUE) {
  excl_cols <- c("sample_id", colnames(phyloseq::sample_data(ps)))
  
  df_topn <- ps %>%
    get_topn(n = n, residuals = residuals) %>%
    ps_to_df(sample_name = "sample_id") %>%
    tidyr::pivot_longer(-dplyr::all_of(excl_cols), names_to = "OTU", values_to = "value") %>%
    dplyr::mutate(OTU = format_OTU(.data$OTU) %>% forcats::fct_inorder() %>% forcats::fct_rev()) %>%
    dplyr::arrange(.data$sample_id) %>%
    dplyr::mutate(sample_id = forcats::fct_inorder(.data$sample_id))
  
  return(df_topn)
}

otu_tab_to_df <- function(ps, sample_name = "sample_id") {
  otu_tab <- phyloseq::otu_table(ps)
  df <- otu_tab %>%
    t %>%
    data.frame() %>%
    tibble::rownames_to_column(sample_name)
  return(df)
}

meta_to_df <- function(ps, sample_name = "sample_id") {
  meta <- phyloseq::sample_data(ps)
  df <- meta %>%
    data.frame() %>%
    tibble::rownames_to_column(sample_name)
  return(df)
}

ps_to_df <- function(ps, sample_name = "sample_id") {
  df_meta <- meta_to_df(ps, sample_name)
  df_otu <- otu_tab_to_df(ps, sample_name)
  
  dplyr::left_join(df_meta, df_otu, by = sample_name) %>% tibble::as_tibble()
}

log10_px <- function(ps, pseudocount = 1) {
  ps %>% transform_sample_counts(., function(x) log10(x + pseudocount)) }

format_OTU <- function(OTU_names, short = F, parse = F) {
  
  OTU_names_tb <- tibble::tibble(OTU_names)
  
  OTU_names_num <- OTU_names_tb %>%
    dplyr::mutate(num = purrr::map_chr(OTU_names, ~stringr::str_extract(., "(?<=_)[0-9]+$")))
  
  OTU_names_tb_format <- OTU_names_num %>%
    dplyr::filter(!is.na(.data$num) & !duplicated(OTU_names)) %>%
    dplyr::mutate(
      name1 = unlist(purrr::map(stringr::str_split(OTU_names, "__|_"), ~ head(.x, -1) %>% glue::glue_collapse(., " "))), #everything except last
      name2 = unlist(purrr::map(stringr::str_split(OTU_names, "__|_"), ~ head(.x, 1))), #last
      long_name = as.character(glue::glue("*{name1}* ({num})")),
      short_name = as.character(glue::glue("*{name2}* ({num})")),
      parse_name = as.character(glue::glue("italic('{name2}')~({num})")))
  
  OTU_names_nonum <- OTU_names_num %>%
    dplyr::filter(is.na(.data$num) & !duplicated(OTU_names))
  
  OTU_names_final <- dplyr::left_join(
    OTU_names_tb,
    dplyr::bind_rows(OTU_names_tb_format, OTU_names_nonum) %>%
      dplyr::mutate(dplyr::across(dplyr::ends_with("_name"), ~dplyr::if_else(is.na(.data$num), stringr::str_to_title(OTU_names), .x)))
    , by = "OTU_names")
  
  if(!short & !parse) {
    return(OTU_names_final %>% dplyr::pull(.data$long_name))
  } else if(parse) {
    return(OTU_names_final %>% dplyr::pull(.data$parse_name))
  } else {
    return(OTU_names_final %>% dplyr::pull(.data$short_name)) }
}


# Plotting functions

prep_bar <- function(ps, n, residuals = TRUE) {
  excl_cols <- c("sample_id", colnames(phyloseq::sample_data(ps)))
  
  df_topn <- ps %>%
    get_topn(n = n, residuals = residuals) %>%
    ps_to_df(sample_name = "sample_id") %>%
    tidyr::pivot_longer(-dplyr::all_of(excl_cols), names_to = "OTU", values_to = "value") %>%
    dplyr::mutate(OTU = format_OTU(OTU) %>% forcats::fct_inorder() %>% forcats::fct_rev()) %>%
    dplyr::arrange(sample_id) %>%
    dplyr::mutate(sample_id = forcats::fct_inorder(sample_id))
  
  return(df_topn)
}

make_color_scheme <- function(name, n) {
  max_n <-  RColorBrewer::brewer.pal.info %>%
    tibble::rownames_to_column("name_pal") %>%
    dplyr::filter(.data$name_pal == name) %>%
    dplyr::pull(.data$maxcolors)
  
  grDevices::colorRampPalette(RColorBrewer::brewer.pal(
    dplyr::case_when(n < 3 ~ 3, n < max_n ~ n, TRUE ~ as.numeric(max_n)), name))(n)
}

create_bar <- function(ps = NULL, df_topn = NULL, id = "sample_id", y = "value", n = 15, ncol_legend = 3, name_legend = "ASV", RA = TRUE, colour = "white") {
  
  # accepts either ps (running prep_bar) or a dataframe already prepared with prep_bar
  if(any(class(ps)=="phyloseq")) { df_topn <- prep_bar(ps = ps, n = n) }
  
  bar <- df_topn %>%
    ggplot2::ggplot(ggplot2::aes_string(x = id, y = y, fill = "OTU")) +
    ggplot2::geom_bar(stat="identity", colour=colour) +
    ggplot2::scale_fill_manual(name = name_legend, values = c("grey90", rev(make_color_scheme("Paired", n)))) +
    ggplot2::labs(x = "" , y = "Number of reads") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
                   legend.position = "bottom", legend.text = ggtext::element_markdown(),
                   panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank()) +
    ggplot2::guides(fill = ggplot2::guide_legend(ncol = ncol_legend))
  
  if(RA) {
    bar <- bar +
      ggplot2::scale_y_continuous(expand = c(0.01, 0.01), labels = scales::percent) +
      ggplot2::ylab("Relative abundance")
  }
  return(bar)
  #https://bookdown.org/rdpeng/RProgDA/non-standard-evaluation.html
  
  #TODO: add y = value to plot mean values 
}

create_bar_revleg <- function(ps = NULL, df_topn = NULL, id = "sample_id", y = "value", n = 15, ncol_legend = 3, name_legend = "ASV", RA = TRUE, colour = "white") {
  
  # accepts either ps (running prep_bar) or a dataframe already prepared with prep_bar
  if(any(class(ps)=="phyloseq")) { df_topn <- prep_bar(ps = ps, n = n) }
  
  bar <- df_topn %>%
    ggplot2::ggplot(ggplot2::aes_string(x = id, y = y, fill = "OTU")) +
    ggplot2::geom_bar(stat="identity", colour=colour) +
    ggplot2::scale_fill_manual(name = name_legend, values = c("grey90", rev(make_color_scheme("Paired", n)))) +
    ggplot2::labs(x = "" , y = "Number of reads") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
                   legend.position = "bottom", legend.text = ggtext::element_markdown(),
                   panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank()) +
    ggplot2::guides(fill = ggplot2::guide_legend(ncol = ncol_legend, reverse = T))
  
  if(RA) {
    bar <- bar +
      ggplot2::scale_y_continuous(expand = c(0.01, 0.01), labels = scales::percent) +
      ggplot2::ylab("Relative abundance")
  }
  return(bar)
  #https://bookdown.org/rdpeng/RProgDA/non-standard-evaluation.html
  
  #TODO: add y = value to plot mean values 
}

create_bar_revleg2 <- function(ps = NULL, df_topn = NULL, id = "sample_id", y = "value", n = 15, ncol_legend = 3, name_legend = "ASV", RA = TRUE, colour = "white") {
  
  # accepts either ps (running prep_bar) or a dataframe already prepared with prep_bar
  if(any(class(ps)=="phyloseq")) { df_topn <- prep_bar(ps = ps, n = n) }
  
  bar <- df_topn %>%
    ggplot2::ggplot(ggplot2::aes_string(x = id, y = y, fill = "OTU")) +
    ggplot2::geom_bar(stat="identity", colour=colour) +
    ggplot2::scale_fill_manual(name = name_legend, values = c("grey90", rev(make_color_scheme("Paired", n)))) +
    ggplot2::labs(x = "" , y = "Number of reads") + theme_bw(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
                   legend.position = "bottom", legend.text = ggtext::element_markdown(),
                   panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank()) +
    ggplot2::guides(fill = ggplot2::guide_legend(ncol = ncol_legend, reverse = T))
  
  if(RA) {
    bar <- bar +
      ggplot2::scale_y_continuous(expand = c(0.01, 0.01), labels = scales::percent) +
      ggplot2::ylab("Relative abundance")
  }
  return(bar)
  #https://bookdown.org/rdpeng/RProgDA/non-standard-evaluation.html
  
  #TODO: add y = value to plot mean values 
}


create_bar2 <- function(ps = NULL, df_topn = NULL, id = "sample_id", y = "value", n = 15, ncol_legend = 3, name_legend = "ASV", RA = TRUE, colour = "white") {
  
  # accepts either ps (running prep_bar) or a dataframe already prepared with prep_bar
  if(any(class(ps)=="phyloseq")) { df_topn <- prep_bar(ps = ps, n = n) }
  
  bar <- df_topn %>%
    ggplot2::ggplot(ggplot2::aes_string(x = id, y = y, fill = "OTU")) +
    ggplot2::geom_bar(stat="identity", colour=colour) +
    ggplot2::scale_fill_manual(name = name_legend, values = c("grey90", rev(make_color_scheme("Spectral", n)))) +
    ggplot2::labs(x = "" , y = "Number of reads") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
                   legend.position = "bottom", legend.text = ggtext::element_markdown(),
                   panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank()) +
    ggplot2::guides(fill = ggplot2::guide_legend(ncol = ncol_legend))
  
  if(RA) {
    bar <- bar +
      ggplot2::scale_y_continuous(expand = c(0.01, 0.01), labels = scales::percent) +
      ggplot2::ylab("Relative abundance")
  }
  return(bar)
  #https://bookdown.org/rdpeng/RProgDA/non-standard-evaluation.html
  
  #TODO: add y = value to plot mean values 
}

create_dendro_gg <- function(hc, hang_height=0.05) {
  
  library(dendextend)
  library(ggdendro)
  
  dendro_data_bl <- hc %>% as.dendrogram %>% hang.dendrogram(hang_height=hang_height) %>% dendro_data
  dendro_data_bl$segments$yend[dendro_data_bl$segments$yend<0] <- 0
  
  dendro_data_gr <- hc %>% as.dendrogram %>% dendro_data
  hc_order <- dendro_data_gr$labels$label
  
  plot <- ggplot(segment(dendro_data_gr)) +
    geom_segment(aes(x=x, y=y, xend=xend, yend=yend),colour="grey75") + 
    geom_segment(data=segment(dendro_data_bl),aes(x=x, y=y, xend=xend, yend=yend)) +
    scale_x_continuous(expand = rep(1/length(hc_order)/2, 2)) + 
    scale_y_continuous(expand=c(0,0.02)) +
    #theme(plot.margin=unit(c(0,0,0,0),"lines")) +
    theme_void()
  return(list(plot=plot, hc_order=hc_order)) 
}

meta_bar <- function(data, var, name = NULL, color_scale = NULL, ...) {
  data <- data %>% mutate(index=1:nrow(.))
  
  p <- ggplot(data, aes_string(x = "index", y = 1, fill = var)) + 
    geom_tile() + 
    scale_y_continuous(expand = c(0,0), breaks = 1, labels = name) + 
    scale_x_continuous(expand = c(0,0)) + 
    theme(axis.title = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(), legend.key.size = unit(0.2, "in"), legend.spacing = unit(0.2, "line"), plot.margin=unit(c(0, 0, 0, 0), "lines"))
  
  if(!is.null(color_scale)) { 
    p <- p + color_scale
  } else {
    n_col <- length(levels(meta_hm_order[[var]])) %>% as.numeric()
    p <- p + scale_fill_manual(values = c(make_color_scheme("BrBG", n_col -1), "grey70"))
  }
  return(p)
}

######################
#additions when doing clustering analysis

vegan_otu <-  function(physeq){   #source: https://rdrr.io/github/taowenmicro/ggClusterNet/src/R/utls.R
  OTU <-  otu_table(physeq)
  if(taxa_are_rows(OTU)){
    OTU <-  t(OTU)
  }
  return(as(OTU,"matrix"))
}


stacked.bar=function(ps,n_otus,hc_order=T,method="average",cols=NULL,var1=NULL,var2=NULL, plot=TRUE,by_genus=F, ra=F){
  
  if(by_genus==F){
     if(ra==F){
    otu_RA <- ps %>%
      transform_sample_counts(., function(x) x/sum(x)) 
     }else{otu_RA=ps}
    
    otu_RA= otu_RA%>%
      otu_table
  }else{
    otu_RA= t(genus_table(ps=ps))
  }
  
  otu_RA_ord <- otu_RA[order(rowMeans(otu_RA), decreasing=T), ]
  
  if(hc_order) {
    bc <- vegdist(t(otu_RA_ord), method="bray")
    hc <- hclust(bc, method=method) 
    otu_RA_ord <- otu_RA_ord[, hc$order]
  } 
  
  stack.df <- t(otu_RA_ord[1:n_otus, ]) %>% 
    data.frame(., Residuals=1-rowSums(.), check.names = F) %>%
    rownames_to_column("sample_id") %>%
    gather(key=OTU, value=RA, -sample_id) %>% 
    mutate(OTU = fct_rev(fct_inorder(OTU)),
           sample_id=fct_inorder(sample_id))
  
  stack.df_meta <- left_join(stack.df, sample_data(ps) %>%
                                 data.frame(check.names=F) %>%
                                 rownames_to_column("sample_id") %>%
                                 #select(sample_id, {{var1}},{{var2}}) %>%
                                 mutate(sample_id=as.factor(sample_id)), by="sample_id") %>%
      mutate(OTU = format_OTU(OTU) %>% forcats::fct_inorder() %>% forcats::fct_rev())
    
  # define your color palette (see https://medialab.github.io/iwanthue/) or use Rcolorbrewer
  #cols<-c("gray50", "#00a3b8","#f17f09","#bd83ff","#e3ba00","#6e3e8a","#51e082","#ff73c3","#00952f","#ff6180","#00b281","#ae3500","#6bd7d7", "#ff9e75","#00897f","#2c5d34") 
  set.seed(19)
  
  if(is.null(cols)){
    cols <- c("grey90", rev(make_color_scheme("Paired", n=n_otus)))} else{cols=cols}
  #
  
  bar1 <- stack.df_meta %>%
    ggplot(aes(x = sample_id, y = RA, fill = OTU)) + 
    geom_bar(stat = "identity", position = "stack", colour = "white") + 
    scale_fill_manual(name = "ASV", values = cols) +
    theme(legend.key.size =  unit(0.2, "in"), legend.spacing = unit(0.2, "line"), plot.margin=unit(c(0,0,0,0),"lines"), 
          axis.title.x = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position = "bottom",
          legend.text = ggtext::element_markdown()) +
    guides(fill = guide_legend(ncol = 4, reverse = T)) + 
    scale_x_discrete(name="Sample") + 
    scale_y_continuous(labels=scales::percent) +
    ylab("Relative abundance") 
  
  #if(!is.null(var)){
  #  return(bar1+var_grid(~var, scales = "free_x", space="free_x") )
  #}else{ 
  
  if(plot){return(bar1)
  }else{return(stack.df_meta)
  }}


meta_plot <- function(data, var, name = NULL, fillScale = NULL, ...) {
  data <- data %>% mutate(index=1:nrow(.))
  p <- ggplot(data, aes_string(x="index", y=1, fill=var)) + 
    geom_tile() + 
    scale_y_continuous(expand=c(0,0), breaks=1, labels=name) + 
    scale_x_continuous(expand=c(0,0)) + 
    theme(axis.title=element_blank(), axis.ticks=element_blank(), axis.text.x=element_blank(), legend.key.size =  unit(0.2, "in"), legend.spacing = unit(0.2, "line"), plot.margin=unit(c(0,0,0,0),"lines"))
  if(!is.null(fillScale)) p <- p + fillScale
  return(p)
}



#Heatmap- cluster plot - copied from Justyna and adjusted to work for me
cluster_plot_hm =  function(ps,k=NULL, n_otus,min_size=3,clust.var=NULL, plot=T,method="average",by_genus=F,clust_colours=F){
  #1.Create distance metric and dendogram
  if(by_genus==F){ 
    bc=vegdist(vegan_otu(ps),method = "bray")
  }else{
    bc=vegdist(genus_table(ps),method = "bray")
  }
  hc= hclust(bc, method={{method}})
  dendro_full <- create_dendro_gg(hc)
  hc_order <- dendro_full$hc_order
  dendro <- dendro_full$plot
  #2.Create heatmap
  
  hm_data <- stacked.bar(ps, n_otus = n_otus, hc_order = T, method = method, plot=F) %>%
    mutate(RA=ifelse(RA==0, NA, RA))
  
  hm <- hm_data %>%
    mutate(sample_id=factor(sample_id, levels=hc_order)) %>%
    ggplot(aes(x=sample_id, y=OTU, fill=RA)) +
    geom_tile(colour="white") +
    scale_fill_gradientn(name="Relative abundance", colors=brewer.pal(9, "Reds"), na.value = "grey90") +
    theme_bw() + 
    theme(axis.text.x = element_blank(),
          axis.text.y = ggtext::element_markdown(size=12),
          axis.title =element_text(size=12), legend.position = "bottom") +
    xlab("Sample") + ylab("ASV")
  
  #3.Add metadata
  if(is.null(clust.var)){
    meta_clust <- hm_data %>%
      dplyr:: slice(match(hc_order, sample_id)) %>%
      mutate(clust=as.character(cutree(hc, k)[hc$order])) #k is the number of clusters
  }else{
    meta_clust <- hm_data %>%
      dplyr:: slice(match(hc_order, sample_id))  
    meta_clust$clust=meta_clust[,{{clust.var}}]   #rename cluster variable if provided
    
    #meta_clust <-hm_data
    #meta_clust$clust=meta_clust[,{{clust.var}}]  ## ${{clust.var}} doesn't work 
  }
  
  if(is.null(clust.var)){
    min_size = {{min_size}} # recode small clusters to NA
    meta_clust_NA <- meta_clust %>% 
      group_by(clust) %>%
      mutate(n=n()) %>% ungroup %>%
      mutate(clust=ifelse(n<=min_size, NA, clust) %>% as.factor() %>% fct_infreq())
  }else{
    meta_clust_NA <- meta_clust
  }
  
  cols_clust <- brewer.pal(length(levels(meta_clust_NA$clust)),"BrBG")
  clust_plot <- meta_clust_NA %>%
    meta_plot("clust", "Cluster", fillScale=scale_fill_manual(name="Cluster", 
                                                              values=cols_clust, na.value="grey"))+
    theme(axis.text.y = element_text(size=12), legend.position = "none")
  
  #5. Combine
  comb_plot <- cowplot::plot_grid(dendro, clust_plot, hm , ncol=1, rel_heights = c(1.5,0.3,5), align="v", axis="lr")

  if(plot){print(comb_plot)}else{return(meta_clust_NA)}
}

#same as above but input is RA ps object
cluster_plot_hm_RA =  function(ps,k=NULL, n_otus,min_size=3,clust.var=NULL, plot=T,method="average",by_genus=F,clust_colours=F){
  #1.Create distance metric and dendogram
  if(by_genus==F){ 
    bc=vegdist(vegan_otu(ps),method = "bray")
  }else{
    bc=vegdist(genus_table(ps),method = "bray")
  }
  hc= hclust(bc, method={{method}})
  dendro_full <- create_dendro_gg(hc)
  hc_order <- dendro_full$hc_order
  dendro <- dendro_full$plot
  #2.Create heatmap
  
  hm_data <- stacked.bar(ps, n_otus = n_otus, hc_order = T, method = method, plot=F, ra=T) %>%
    mutate(RA=ifelse(RA==0, NA, RA))
  
  hm <- hm_data %>%
    mutate(sample_id=factor(sample_id, levels=hc_order)) %>%
    ggplot(aes(x=sample_id, y=OTU, fill=RA)) +
    geom_tile(colour="white") +
    scale_fill_gradientn(name="Relative abundance", colors=brewer.pal(9, "Reds"), na.value = "grey90") +
    theme_bw() + 
    theme(axis.text.x = element_blank(),
          axis.text.y = ggtext::element_markdown(size=12),
          axis.title =element_text(size=12), legend.position = "bottom") +
    xlab("Sample") + ylab("ASV")
  
  #3.Add metadata
  if(is.null(clust.var)){
    meta_clust <- hm_data %>%
      dplyr:: slice(match(hc_order, sample_id)) %>%
      mutate(clust=as.character(cutree(hc, k)[hc$order])) #k is the number of clusters
  }else{
    meta_clust <- hm_data %>%
      dplyr:: slice(match(hc_order, sample_id))  
    meta_clust$clust=meta_clust[,{{clust.var}}]   #rename cluster variable if provided
    
    #meta_clust <-hm_data
    #meta_clust$clust=meta_clust[,{{clust.var}}]  ## ${{clust.var}} doesn't work 
  }
  
  if(is.null(clust.var)){
    min_size = {{min_size}} # recode small clusters to NA
    meta_clust_NA <- meta_clust %>% 
      group_by(clust) %>%
      mutate(n=n()) %>% ungroup %>%
      mutate(clust=ifelse(n<=min_size, NA, clust) %>% as.factor() %>% fct_infreq())
  }else{
    meta_clust_NA <- meta_clust
  }
  
  cols_clust <- brewer.pal(length(levels(meta_clust_NA$clust)),"BrBG")
  #replace the colour that is pretty much white
  cols_clust <- replace(cols_clust, cols_clust=="#F5F5F5", "#D07D38")
  clust_plot <- meta_clust_NA %>%
    meta_plot("clust", "Cluster", fillScale=scale_fill_manual(name="Cluster", 
                                                              values=cols_clust, na.value="grey"))+
    theme(axis.text.y = element_text(size=12), legend.position = "none")
  
  #5. Combine
  comb_plot <- cowplot::plot_grid(dendro, clust_plot, hm , ncol=1, rel_heights = c(1.5,0.3,5), align="v", axis="lr")
  
  if(plot){print(comb_plot)}else{return(meta_clust_NA)}
}

#as above but with an additional sample type bar added
cluster_plot_hm_sampletype =  function(ps,k=NULL, n_otus,min_size=3,clust.var=NULL, plot=T,method="average",by_genus=F,clust_colours=F){
  #1.Create distance metric and dendogram
  if(by_genus==F){ 
    bc=vegdist(vegan_otu(ps),method = "bray")
  }else{
    bc=vegdist(genus_table(ps),method = "bray")
  }
  hc= hclust(bc, method={{method}})
  dendro_full <- create_dendro_gg(hc)
  hc_order <- dendro_full$hc_order
  dendro <- dendro_full$plot
  #2.Create heatmap
  
  hm_data <- stacked.bar(ps, n_otus = n_otus, hc_order = T, method = method, plot=F) %>%
    mutate(RA=ifelse(RA==0, NA, RA))
  
  hm <- hm_data %>%
    mutate(sample_id=factor(sample_id, levels=hc_order)) %>%
    ggplot(aes(x=sample_id, y=OTU, fill=RA)) +
    geom_tile(colour="white") +
    scale_fill_gradientn(name="Relative abundance", colors=brewer.pal(9, "Reds"), na.value = "grey90") +
    theme_bw() + 
    theme(axis.text.x = element_blank(),
          axis.text.y = ggtext::element_markdown(size=12),
          axis.title =element_text(size=12), legend.position = "bottom") +
    xlab("Sample") + ylab("ASV")
  
  #3.Add metadata
  if(is.null(clust.var)){
    meta_clust <- hm_data %>%
      dplyr:: slice(match(hc_order, sample_id)) %>%
      mutate(clust=as.character(cutree(hc, k)[hc$order])) #k is the number of clusters
  }else{
    meta_clust <- hm_data %>%
      dplyr:: slice(match(hc_order, sample_id))  
    meta_clust$clust=meta_clust[,{{clust.var}}]   #rename cluster variable if provided
    
    #meta_clust <-hm_data
    #meta_clust$clust=meta_clust[,{{clust.var}}]  ## ${{clust.var}} doesn't work 
  }
  
  if(is.null(clust.var)){
    min_size = {{min_size}} # recode small clusters to NA
    meta_clust_NA <- meta_clust %>% 
      group_by(clust) %>%
      mutate(n=n()) %>% ungroup %>%
      mutate(clust=ifelse(n<=min_size, NA, clust) %>% as.factor() %>% fct_infreq())
  }else{
    meta_clust_NA <- meta_clust
  }
  
  cols_clust <- brewer.pal(length(levels(meta_clust_NA$clust)),"BrBG")
  cols_clust <- replace(cols_clust, cols_clust=="#F5F5F5", "#D07D38")
  clust_plot <- meta_clust_NA %>%
    meta_plot("clust", "Cluster", fillScale=scale_fill_manual(name="Cluster", 
                                                              values=cols_clust, na.value="grey"))+
    theme(axis.text.y = element_text(size=12), legend.position = "none")
  
  #add the sample type metadata
  cols_sample <- brewer.pal(length(levels(meta_clust$type_timepoint)),"Pastel1")
  
  sampletype <- meta_clust %>%
    meta_plot("type_timepoint", "Sample timepoint", fillScale=scale_fill_manual(name="Sample timepoint", values=cols_sample, na.value="grey")) +
    theme(axis.text.y = element_text(size=12), legend.position = "right") + guides(fill = guide_legend(ncol = 3))
  
  
  #5. Combine
  comb_plot <- cowplot::plot_grid(dendro, clust_plot, sampletype, hm , ncol=1, rel_heights = c(1.5,0.3,0.3, 4), align="v", axis="lr")
  
  if(plot){print(comb_plot)}else{return(meta_clust_NA)}
}

#as above but updated colours and based on relative abundances
cluster_plot_hm_sampletype2_RA =  function(ps,k=NULL, n_otus,min_size=3,clust.var=NULL, plot=T,method="average",by_genus=F,clust_colours=F){
  #1.Create distance metric and dendogram
  if(by_genus==F){ 
    bc=vegdist(vegan_otu(ps),method = "bray")
  }else{
    bc=vegdist(genus_table(ps),method = "bray")
  }
  hc= hclust(bc, method={{method}})
  dendro_full <- create_dendro_gg(hc)
  hc_order <- dendro_full$hc_order
  dendro <- dendro_full$plot
  #2.Create heatmap
  
  hm_data <- stacked.bar(ps, n_otus = n_otus, hc_order = T, method = method, plot=F, ra=T) %>%
    mutate(RA=ifelse(RA==0, NA, RA))
  
  hm <- hm_data %>%
    mutate(sample_id=factor(sample_id, levels=hc_order)) %>%
    ggplot(aes(x=sample_id, y=OTU, fill=RA)) +
    geom_tile(colour="white") +
    scale_fill_gradientn(name="Relative abundance", colors=brewer.pal(9, "Reds"), na.value = "grey90") +
    theme_bw() + 
    theme(axis.text.x = element_blank(),
          axis.text.y = ggtext::element_markdown(size=12),
          axis.title =element_text(size=12), legend.position = "bottom") +
    xlab("Sample") + ylab("ASV")
  
  #3.Add metadata
  if(is.null(clust.var)){
    meta_clust <- hm_data %>%
      dplyr:: slice(match(hc_order, sample_id)) %>%
      mutate(clust=as.character(cutree(hc, k)[hc$order])) #k is the number of clusters
  }else{
    meta_clust <- hm_data %>%
      dplyr:: slice(match(hc_order, sample_id))  
    meta_clust$clust=meta_clust[,{{clust.var}}]   #rename cluster variable if provided
    
    #meta_clust <-hm_data
    #meta_clust$clust=meta_clust[,{{clust.var}}]  ## ${{clust.var}} doesn't work 
  }
  
  if(is.null(clust.var)){
    min_size = {{min_size}} # recode small clusters to NA
    meta_clust_NA <- meta_clust %>% 
      group_by(clust) %>%
      mutate(n=n()) %>% ungroup %>%
      mutate(clust=ifelse(n<=min_size, NA, clust) %>% as.factor() %>% fct_infreq())
  }else{
    meta_clust_NA <- meta_clust
  }
  
  cols_clust <- brewer.pal(length(levels(meta_clust_NA$clust)),"BrBG")
  cols_clust <- replace(cols_clust, cols_clust=="#F5F5F5", "#D07D38")
  clust_plot <- meta_clust_NA %>%
    meta_plot("clust", "Cluster", fillScale=scale_fill_manual(name="Cluster", 
                                                              values=cols_clust, na.value="grey"))+
    theme(axis.text.y = element_text(size=12), legend.position = "none")
  
  #add the sample type metadata
  cols_sample <- c("#1B9E77", "#FDC46B", "#7570B3")
  
  sampletype <- meta_clust %>%
    meta_plot("type_timepoint", "Sample type", fillScale=scale_fill_manual(name="", values=cols_sample, na.value="grey", 
                                                                           labels=c("0_V1" = "Term\nmeconium", "1_V1" = "Preterm\nmeconium",
                                                                                    "1_V2" = "Preterm\npre-discharge"))) +
    theme(axis.text.y = element_text(size=12), legend.position = "right", legend.title=element_blank()) + guides(fill = guide_legend(ncol = 3))
  
  
  #5. Combine
  comb_plot <- cowplot::plot_grid(dendro, clust_plot, sampletype, hm , ncol=1, rel_heights = c(1.5,0.3,0.3, 4), align="v", axis="lr")
  
  if(plot){print(comb_plot)}else{return(meta_clust_NA)}
}


cluster_plot_bar <- function(ps,k=NULL, n_otus,min_size=3, method="complete") {
  
  #calculate bray curtis matrix
  # extracting data from phyloseq object to dataframes
  ASVm <- as(otu_table(ps), "matrix")
  
  # transpose if necessary
  tASVm <- t(ASVm)
  
  # coerce to data.frame
  ASV_raw <- as.data.frame(tASVm)
  
  # calculate bray-curtis dissimilarity matrix
  bc <- vegdist(ASV_raw, method = "bray") # ASV_raw needs to be raw counts
  
  #do clustering
  hc <- hclust(bc, method=method)
  
  #do dendrogram
  dendro_full <- create_dendro_gg(hc)
  hc_order <- dendro_full$hc_order
  dendro <- dendro_full$plot
  
  #now do the barplot
  #make ps object into RA
  ps_RA <- ps %>% to_RA()
  
  # Order OTU table by decreasing abundance
  otu_RA_ord <- otu_table(ps_RA)[order(rowSums(otu_table(ps_RA)), decreasing=T),]
  
  # melt OTU table in long format, i.e. change the format of the OTU table such that all OTUs are listed in one column and their relative abundances in another column (note the use of the pipe operator %>%). The below code will list the relative abundances for the top 15 most abundant OTUs and will sum up the remaining OTUs as "Residuals" for ease of plotting
  stack.df <- t(otu_RA_ord[1:n_otus, ]) %>% 
    data.frame(., Residuals=1-rowSums(.), check.names = F) %>%
    rownames_to_column("sample.id") %>%
    gather(key=OTU, value=RA, -sample.id) %>% 
    mutate(OTU = fct_rev(fct_inorder(OTU)))
  
  # append relevant metadata (time)
  stack.df_meta <- left_join(sample_data(ps_RA) %>%
                               data.frame(check.names=F) %>%
                               rownames_to_column("sample.id") %>%
                               #select(sample.id, time) %>%
                               mutate(sample.id=as.character(sample.id)), stack.df, by="sample.id") %>%
    mutate(OTU = format_OTU(OTU) %>% forcats::fct_inorder() %>% forcats::fct_rev())
  
  cols <- c("grey90", rev(make_color_scheme("Paired", n=n_otus)))
  
  bar <- stack.df_meta  %>%
    mutate(sample.id=factor(sample.id, levels=hc_order)) %>%
    ggplot(aes(x=sample.id, y=RA, fill=OTU)) +
    geom_bar(stat="identity", colour="white") +
    xlab("Sample") + 
    ylab("Relative abundance") +
    scale_fill_manual(name = "ASV", values=cols) +
    scale_y_continuous(labels=scales::percent, expand=c(0.01, 0.01)) + 
    theme(legend.key.size =  unit(0.2, "in"), legend.spacing = unit(0.2, "line"), plot.margin=unit(c(0,0,0,0),"lines"), 
          axis.title.x = element_blank(), panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), panel.background = element_blank(), 
          axis.text.x=element_blank(), axis.ticks.x=element_blank(),
          legend.text = ggtext::element_markdown()) +
    ggplot2::guides(fill = ggplot2::guide_legend(reverse = T))
  
  
  #add metadata
  meta_plot <- function(data, var, name = NULL, fillScale = NULL, ...) {
    data <- data %>% mutate(index=1:nrow(.))
    p <- ggplot(data, aes_string(x="index", y=1, fill=var)) + 
      geom_tile() + 
      scale_y_continuous(expand=c(0,0), breaks=1, labels=name) + 
      scale_x_continuous(expand=c(0,0)) + 
      theme(axis.title=element_blank(), axis.ticks=element_blank(), axis.text.x=element_blank(), legend.key.size =  unit(0.2, "in"), legend.spacing = unit(0.2, "line"), plot.margin=unit(c(0,0,0,0),"lines"), legend.position="none")
    if(!is.null(fillScale)) p <- p + fillScale
    return(p)
  }
  
  meta_clust <- stack.df_meta %>%
    dplyr::slice(match(hc_order, sample.id)) %>%
    mutate(clust=as.character(cutree(hc, k)[hc$order]))
  
  # recode small clusters to NA
  meta_clust_NA <- meta_clust %>% 
    group_by(clust) %>%
    mutate(n=n()) %>% ungroup %>%
    mutate(clust=ifelse(n<=min_size, NA, clust) %>% as.factor() %>% fct_infreq())
  
  cols_clust <- brewer.pal(length(levels(meta_clust_NA$clust)),"BrBG")
  cols_clust <- replace(cols_clust, cols_clust="#F5F5F5", "#D07D38")
  
  clust <- meta_clust_NA %>%
    meta_plot("clust", "Cluster", fillScale=scale_fill_manual(name="Cluster", 
                                                              values=cols_clust, na.value="grey"))+
    theme(axis.text.y = element_text(size=12), legend.position = "none")
  
  #combine plots
  cowplot::plot_grid(dendro, clust, bar + theme(legend.position = "bottom"), ncol=1, rel_heights = c(1.5,0.3,5), align="v", axis="lr")
  
}

#as above but based on RA ps object
cluster_plot_bar_RA <- function(ps,k=NULL, n_otus,min_size=3, method="complete") {
  
  #calculate bray curtis matrix
  # extracting data from phyloseq object to dataframes
  ASVm <- as(otu_table(ps), "matrix")
  
  # transpose if necessary
  tASVm <- t(ASVm)
  
  # coerce to data.frame
  ASV_raw <- as.data.frame(tASVm)
  
  # calculate bray-curtis dissimilarity matrix
  bc <- vegdist(ASV_raw, method = "bray") # ASV_raw needs to be raw counts
  
  #do clustering
  hc <- hclust(bc, method=method)
  
  #do dendrogram
  dendro_full <- create_dendro_gg(hc)
  hc_order <- dendro_full$hc_order
  dendro <- dendro_full$plot
  
  #now do the barplot
  #make ps object into RA - not needed as input already RA table
  ps_RA <- ps #%>% to_RA()
  
  # Order OTU table by decreasing abundance
  otu_RA_ord <- otu_table(ps_RA)[order(rowSums(otu_table(ps_RA)), decreasing=T),]
  
  # melt OTU table in long format, i.e. change the format of the OTU table such that all OTUs are listed in one column and their relative abundances in another column (note the use of the pipe operator %>%). The below code will list the relative abundances for the top 15 most abundant OTUs and will sum up the remaining OTUs as "Residuals" for ease of plotting
  stack.df <- t(otu_RA_ord[1:n_otus, ]) %>% 
    data.frame(., Residuals=1-rowSums(.), check.names = F) %>%
    rownames_to_column("sample.id") %>%
    gather(key=OTU, value=RA, -sample.id) %>% 
    mutate(OTU = fct_rev(fct_inorder(OTU)))
  
  # append relevant metadata (time)
  stack.df_meta <- left_join(sample_data(ps_RA) %>%
                               data.frame(check.names=F) %>%
                               rownames_to_column("sample.id") %>%
                               #select(sample.id, time) %>%
                               mutate(sample.id=as.character(sample.id)), stack.df, by="sample.id") %>%
    mutate(OTU = format_OTU(OTU) %>% forcats::fct_inorder() %>% forcats::fct_rev())
  
  cols <- c("grey90", rev(make_color_scheme("Paired", n=n_otus)))
  
  bar <- stack.df_meta  %>%
    mutate(sample.id=factor(sample.id, levels=hc_order)) %>%
    ggplot(aes(x=sample.id, y=RA, fill=OTU)) +
    geom_bar(stat="identity", colour="white") +
    xlab("Sample") + 
    ylab("Relative abundance") +
    scale_fill_manual(name = "ASV", values=cols) +
    scale_y_continuous(labels=scales::percent, expand=c(0.01, 0.01)) + 
    theme(legend.key.size =  unit(0.2, "in"), legend.spacing = unit(0.2, "line"), plot.margin=unit(c(0,0,0,0),"lines"), 
          axis.title.x = element_blank(), panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), panel.background = element_blank(), 
          axis.text.x=element_blank(), axis.ticks.x=element_blank(),
          legend.text = ggtext::element_markdown()) +
    ggplot2::guides(fill = ggplot2::guide_legend(reverse = T))
  
  
  #add metadata
  meta_plot <- function(data, var, name = NULL, fillScale = NULL, ...) {
    data <- data %>% mutate(index=1:nrow(.))
    p <- ggplot(data, aes_string(x="index", y=1, fill=var)) + 
      geom_tile() + 
      scale_y_continuous(expand=c(0,0), breaks=1, labels=name) + 
      scale_x_continuous(expand=c(0,0)) + 
      theme(axis.title=element_blank(), axis.ticks=element_blank(), axis.text.x=element_blank(), legend.key.size =  unit(0.2, "in"), legend.spacing = unit(0.2, "line"), plot.margin=unit(c(0,0,0,0),"lines"), legend.position="none")
    if(!is.null(fillScale)) p <- p + fillScale
    return(p)
  }
  
  meta_clust <- stack.df_meta %>%
    dplyr::slice(match(hc_order, sample.id)) %>%
    mutate(clust=as.character(cutree(hc, k)[hc$order]))
  
  # recode small clusters to NA
  meta_clust_NA <- meta_clust %>% 
    group_by(clust) %>%
    mutate(n=n()) %>% ungroup %>%
    mutate(clust=ifelse(n<=min_size, NA, clust) %>% as.factor() %>% fct_infreq())
  
  cols_clust <- brewer.pal(length(levels(meta_clust_NA$clust)),"BrBG")
  cols_clust <- replace(cols_clust, cols_clust=="#F5F5F5", "#D07D38")
  
  clust <- meta_clust_NA %>%
    meta_plot("clust", "Cluster", fillScale=scale_fill_manual(name="Cluster", 
                                                              values=cols_clust, na.value="grey"))+
    theme(axis.text.y = element_text(size=12), legend.position = "none")
  
  #combine plots
  cowplot::plot_grid(dendro, clust, bar + theme(legend.position = "bottom"), ncol=1, rel_heights = c(1.5,0.3,5), align="v", axis="lr")
  
}


#as above but adding another bar that specifies which sample type is in a cluster
cluster_plot_bar_sampletype <- function(ps,k=NULL, n_otus,min_size=3, method="complete") {
  
  #calculate bray curtis matrix
  # extracting data from phyloseq object to dataframes
  ASVm <- as(otu_table(ps), "matrix")
  
  # transpose if necessary
  tASVm <- t(ASVm)
  
  # coerce to data.frame
  ASV_raw <- as.data.frame(tASVm)
  
  # calculate bray-curtis dissimilarity matrix
  bc <- vegdist(ASV_raw, method = "bray") # ASV_raw needs to be raw counts
  
  #do clustering
  hc <- hclust(bc, method=method)
  
  #do dendrogram
  dendro_full <- create_dendro_gg(hc)
  hc_order <- dendro_full$hc_order
  dendro <- dendro_full$plot
  
  #now do the barplot
  #make ps object into RA
  ps_RA <- ps %>% to_RA()
  
  # Order OTU table by decreasing abundance
  otu_RA_ord <- otu_table(ps_RA)[order(rowSums(otu_table(ps_RA)), decreasing=T),]
  
  # melt OTU table in long format, i.e. change the format of the OTU table such that all OTUs are listed in one column and their relative abundances in another column (note the use of the pipe operator %>%). The below code will list the relative abundances for the top 15 most abundant OTUs and will sum up the remaining OTUs as "Residuals" for ease of plotting
  stack.df <- t(otu_RA_ord[1:n_otus, ]) %>% 
    data.frame(., Residuals=1-rowSums(.), check.names = F) %>%
    rownames_to_column("sample.id") %>%
    gather(key=OTU, value=RA, -sample.id) %>% 
    mutate(OTU = fct_rev(fct_inorder(OTU)))
  
  # append relevant metadata (time)
  stack.df_meta <- left_join(sample_data(ps_RA) %>%
                               data.frame(check.names=F) %>%
                               rownames_to_column("sample.id") %>%
                               #select(sample.id, time) %>%
                               mutate(sample.id=as.character(sample.id)), stack.df, by="sample.id") %>%
    mutate(OTU = format_OTU(OTU) %>% forcats::fct_inorder() %>% forcats::fct_rev())
  
  cols <- c("grey90", rev(make_color_scheme("Paired", n=n_otus)))
  
  bar <- stack.df_meta  %>%
    mutate(sample.id=factor(sample.id, levels=hc_order)) %>%
    ggplot(aes(x=sample.id, y=RA, fill=OTU)) +
    geom_bar(stat="identity", colour="white") +
    xlab("Sample") + 
    ylab("Relative abundance") +
    scale_fill_manual(name = "ASV", values=cols) +
    scale_y_continuous(labels=scales::percent, expand=c(0.01, 0.01)) + 
    theme(legend.key.size =  unit(0.2, "in"), legend.spacing = unit(0.2, "line"), plot.margin=unit(c(0,0,0,0),"lines"), 
          axis.title.x = element_blank(), panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), panel.background = element_blank(), 
          axis.text.x=element_blank(), axis.ticks.x=element_blank(),
          legend.text = ggtext::element_markdown()) +
    ggplot2::guides(fill = ggplot2::guide_legend(reverse = T))
  
  
  #add metadata
  meta_plot <- function(data, var, name = NULL, fillScale = NULL, ...) {
    data <- data %>% mutate(index=1:nrow(.))
    p <- ggplot(data, aes_string(x="index", y=1, fill=var)) + 
      geom_tile() + 
      scale_y_continuous(expand=c(0,0), breaks=1, labels=name) + 
      scale_x_continuous(expand=c(0,0)) + 
      theme(axis.title=element_blank(), axis.ticks=element_blank(), axis.text.x=element_blank(), legend.key.size =  unit(0.2, "in"), 
            legend.spacing = unit(0.2, "line"), plot.margin=unit(c(0,0,0,0),"lines"), legend.position="right") +
    guides(fill = guide_legend(ncol = 5))
    if(!is.null(fillScale)) p <- p + fillScale
    return(p)
  }
  
  meta_clust <- stack.df_meta %>%
    dplyr::slice(match(hc_order, sample.id)) %>%
    mutate(clust=as.character(cutree(hc, k)[hc$order]))
  
  # recode small clusters to NA
  meta_clust_NA <- meta_clust %>% 
    group_by(clust) %>%
    mutate(n=n()) %>% ungroup %>%
    mutate(clust=ifelse(n<=min_size, NA, clust) %>% as.factor() %>% fct_infreq())
  
  cols_clust <- brewer.pal(length(levels(meta_clust_NA$clust)),"BrBG")
  cols_clust <- replace(cols_clust, cols_clust=="#F5F5F5", "#D07D38")
  
  clust <- meta_clust_NA %>%
    meta_plot("clust", "Cluster", fillScale=scale_fill_manual(name="Cluster", 
                                                              values=cols_clust, na.value="grey"))+
    theme(axis.text.y = element_text(size=12), legend.position = "none")
  
  #add the sample type metadata
  cols_sample <- brewer.pal(length(levels(meta_clust$type_timepoint)),"Pastel1")
  
  sampletype <- meta_clust %>%
    meta_plot("type_timepoint", "Sample timepoint", fillScale=scale_fill_manual(name="Sample timepoint", values=cols_sample, na.value="grey"))
  
  #combine plots
  cowplot::plot_grid(dendro, clust, sampletype, bar + theme(legend.position = "bottom"), ncol=1, rel_heights = c(1.5,0.4,0.4,4), align="v", axis="lr")
  
}

#as above but adding another bar that specifies which sample type is in a cluster - UPDATED COLOURS
cluster_plot_bar_sampletype2 <- function(ps,k=NULL, n_otus,min_size=3, method="complete") {
  
  #calculate bray curtis matrix
  # extracting data from phyloseq object to dataframes
  ASVm <- as(otu_table(ps), "matrix")
  
  # transpose if necessary
  tASVm <- t(ASVm)
  
  # coerce to data.frame
  ASV_raw <- as.data.frame(tASVm)
  
  # calculate bray-curtis dissimilarity matrix
  bc <- vegdist(ASV_raw, method = "bray") # ASV_raw needs to be raw counts
  
  #do clustering
  hc <- hclust(bc, method=method)
  
  #do dendrogram
  dendro_full <- create_dendro_gg(hc)
  hc_order <- dendro_full$hc_order
  dendro <- dendro_full$plot
  
  #now do the barplot
  #make ps object into RA
  ps_RA <- ps %>% to_RA()
  
  # Order OTU table by decreasing abundance
  otu_RA_ord <- otu_table(ps_RA)[order(rowSums(otu_table(ps_RA)), decreasing=T),]
  
  # melt OTU table in long format, i.e. change the format of the OTU table such that all OTUs are listed in one column and their relative abundances in another column (note the use of the pipe operator %>%). The below code will list the relative abundances for the top 15 most abundant OTUs and will sum up the remaining OTUs as "Residuals" for ease of plotting
  stack.df <- t(otu_RA_ord[1:n_otus, ]) %>% 
    data.frame(., Residuals=1-rowSums(.), check.names = F) %>%
    rownames_to_column("sample.id") %>%
    gather(key=OTU, value=RA, -sample.id) %>% 
    mutate(OTU = fct_rev(fct_inorder(OTU)))
  
  # append relevant metadata
  stack.df_meta <- left_join(sample_data(ps_RA) %>%
                               data.frame(check.names=F) %>%
                               rownames_to_column("sample.id") %>%
                               #select(sample.id, time) %>%
                               mutate(sample.id=as.character(sample.id)), stack.df, by="sample.id") %>%
    mutate(OTU = format_OTU(OTU) %>% forcats::fct_inorder() %>% forcats::fct_rev())
  
  cols <- c("grey90", rev(make_color_scheme("Paired", n=n_otus)))
  
  bar <- stack.df_meta  %>%
    mutate(sample.id=factor(sample.id, levels=hc_order)) %>%
    ggplot(aes(x=sample.id, y=RA, fill=OTU)) +
    geom_bar(stat="identity", colour="white") +
    xlab("Sample") + 
    ylab("Relative abundance") +
    scale_fill_manual(name = "ASV", values=cols) +
    scale_y_continuous(labels=scales::percent, expand=c(0.01, 0.01)) + 
    theme(legend.key.size =  unit(0.2, "in"), legend.spacing = unit(0.2, "line"), plot.margin=unit(c(0,0,0,0),"lines"), 
          axis.title.x = element_blank(), panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), panel.background = element_blank(), 
          axis.text.x=element_blank(), axis.ticks.x=element_blank(),
          legend.text = ggtext::element_markdown()) +
    ggplot2::guides(fill = ggplot2::guide_legend(reverse = T))
  
  
  #add metadata
  meta_plot <- function(data, var, name = NULL, fillScale = NULL, ...) {
    data <- data %>% mutate(index=1:nrow(.))
    p <- ggplot(data, aes_string(x="index", y=1, fill=var)) + 
      geom_tile() + 
      scale_y_continuous(expand=c(0,0), breaks=1, labels=name) + 
      scale_x_continuous(expand=c(0,0)) + 
      theme(axis.title=element_blank(), axis.ticks=element_blank(), axis.text.x=element_blank(), legend.key.size =  unit(0.2, "in"), 
            legend.spacing = unit(0.2, "line"), plot.margin=unit(c(0,0,0,0),"lines"), legend.position="right") +
      guides(fill = guide_legend(ncol = 5))
    if(!is.null(fillScale)) p <- p + fillScale
    return(p)
  }
  
  meta_clust <- stack.df_meta %>%
    dplyr::slice(match(hc_order, sample.id)) %>%
    mutate(clust=as.character(cutree(hc, k)[hc$order]))
  
  # recode small clusters to NA
  meta_clust_NA <- meta_clust %>% 
    group_by(clust) %>%
    mutate(n=n()) %>% ungroup %>%
    mutate(clust=ifelse(n<=min_size, NA, clust) %>% as.factor() %>% fct_infreq())
  
  cols_clust <- brewer.pal(length(levels(meta_clust_NA$clust)),"BrBG")
  cols_clust <- replace(cols_clust, cols_clust=="#F5F5F5", "#D07D38")
  
  clust <- meta_clust_NA %>%
    meta_plot("clust", "Cluster", fillScale=scale_fill_manual(name="Cluster", 
                                                              values=cols_clust, na.value="grey"))+
    theme(axis.text.y = element_text(size=12), legend.position = "none")
  
  #add the sample type metadata
  cols_sample <- c("#1B9E77", "#FDC46B", "#7570B3")
  
  sampletype <- meta_clust %>%
    meta_plot("type_timepoint", "Sample type", fillScale=scale_fill_manual(name="", values=cols_sample, na.value="grey", 
                                                                           labels=c("0_V1" = "Term\nmeconium", "1_V1" = "Preterm\nmeconium",
                                                                                    "1_V2" = "Preterm\npre-discharge"))) +
    theme(axis.text.y = element_text(size=12), legend.position = "right", legend.title=element_blank()) + guides(fill = guide_legend(ncol = 3))
  
  #combine plots
  cowplot::plot_grid(dendro, clust, sampletype, bar + theme(legend.position = "bottom"), ncol=1, rel_heights = c(1.5,0.4,0.4,4), align="v", axis="lr")
  
}

#as above but giving relative abundance ps object as the input
cluster_plot_bar_sampletype2_RA <- function(ps,k=NULL, n_otus,min_size=3, method="complete") {
  
  #calculate bray curtis matrix
  # extracting data from phyloseq object to dataframes
  ASVm <- as(otu_table(ps), "matrix")
  
  # transpose if necessary
  tASVm <- t(ASVm)
  
  # coerce to data.frame
  ASV_raw <- as.data.frame(tASVm)
  
  # calculate bray-curtis dissimilarity matrix
  bc <- vegdist(ASV_raw, method = "bray") # ASV_raw needs to be raw counts
  
  #do clustering
  hc <- hclust(bc, method=method)
  
  #do dendrogram
  dendro_full <- create_dendro_gg(hc)
  hc_order <- dendro_full$hc_order
  dendro <- dendro_full$plot
  
  #now do the barplot
  #make ps object into RA - not needed because here already giving the function relative abundance ps object
  ps_RA <- ps #%>% to_RA()
  
  # Order OTU table by decreasing abundance
  otu_RA_ord <- otu_table(ps_RA)[order(rowSums(otu_table(ps_RA)), decreasing=T),]
  
  # melt OTU table in long format, i.e. change the format of the OTU table such that all OTUs are listed in one column and their relative abundances in another column (note the use of the pipe operator %>%). The below code will list the relative abundances for the top 15 most abundant OTUs and will sum up the remaining OTUs as "Residuals" for ease of plotting
  stack.df <- t(otu_RA_ord[1:n_otus, ]) %>% 
    data.frame(., Residuals=1-rowSums(.), check.names = F) %>%
    rownames_to_column("sample.id") %>%
    gather(key=OTU, value=RA, -sample.id) %>% 
    mutate(OTU = fct_rev(fct_inorder(OTU)))
  
  # append relevant metadata
  stack.df_meta <- left_join(sample_data(ps_RA) %>%
                               data.frame(check.names=F) %>%
                               rownames_to_column("sample.id") %>%
                               #select(sample.id, time) %>%
                               mutate(sample.id=as.character(sample.id)), stack.df, by="sample.id") %>%
    mutate(OTU = format_OTU(OTU) %>% forcats::fct_inorder() %>% forcats::fct_rev())
  
  cols <- c("grey90", rev(make_color_scheme("Paired", n=n_otus)))
  
  bar <- stack.df_meta  %>%
    mutate(sample.id=factor(sample.id, levels=hc_order)) %>%
    ggplot(aes(x=sample.id, y=RA, fill=OTU)) +
    geom_bar(stat="identity", colour="white") +
    xlab("Sample") + 
    ylab("Relative abundance") +
    scale_fill_manual(name = "ASV", values=cols) +
    scale_y_continuous(labels=scales::percent, expand=c(0.01, 0.01)) + 
    theme(legend.key.size =  unit(0.2, "in"), legend.spacing = unit(0.2, "line"), plot.margin=unit(c(0,0,0,0),"lines"), 
          axis.title.x = element_blank(), panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), panel.background = element_blank(), 
          axis.text.x=element_blank(), axis.ticks.x=element_blank(),
          legend.text = ggtext::element_markdown()) +
    ggplot2::guides(fill = ggplot2::guide_legend(reverse = T))
  
  
  #add metadata
  meta_plot <- function(data, var, name = NULL, fillScale = NULL, ...) {
    data <- data %>% mutate(index=1:nrow(.))
    p <- ggplot(data, aes_string(x="index", y=1, fill=var)) + 
      geom_tile() + 
      scale_y_continuous(expand=c(0,0), breaks=1, labels=name) + 
      scale_x_continuous(expand=c(0,0)) + 
      theme(axis.title=element_blank(), axis.ticks=element_blank(), axis.text.x=element_blank(), legend.key.size =  unit(0.2, "in"), 
            legend.spacing = unit(0.2, "line"), plot.margin=unit(c(0,0,0,0),"lines"), legend.position="right") +
      guides(fill = guide_legend(ncol = 5))
    if(!is.null(fillScale)) p <- p + fillScale
    return(p)
  }
  
  meta_clust <- stack.df_meta %>%
    dplyr::slice(match(hc_order, sample.id)) %>%
    mutate(clust=as.character(cutree(hc, k)[hc$order]))
  
  # recode small clusters to NA
  meta_clust_NA <- meta_clust %>% 
    group_by(clust) %>%
    mutate(n=n()) %>% ungroup %>%
    mutate(clust=ifelse(n<=min_size, NA, clust) %>% as.factor() %>% fct_infreq())
  
  cols_clust <- brewer.pal(length(levels(meta_clust_NA$clust)),"BrBG")
  cols_clust <- replace(cols_clust, cols_clust=="#F5F5F5", "#D07D38")
  
  clust <- meta_clust_NA %>%
    meta_plot("clust", "Cluster", fillScale=scale_fill_manual(name="Cluster", 
                                                              values=cols_clust, na.value="grey"))+
    theme(axis.text.y = element_text(size=12), legend.position = "none")
  
  #add the sample type metadata
  cols_sample <- c("#1B9E77", "#FDC46B", "#7570B3")
  
  sampletype <- meta_clust %>%
    meta_plot("type_timepoint", "Sample type", fillScale=scale_fill_manual(name="", values=cols_sample, na.value="grey", 
                                                                           labels=c("0_V1" = "Term\nmeconium", "1_V1" = "Preterm\nmeconium",
                                                                                    "1_V2" = "Preterm\npre-discharge"))) +
    theme(axis.text.y = element_text(size=12), legend.position = "right", legend.title=element_blank()) + guides(fill = guide_legend(ncol = 3))
  
  #combine plots
  cowplot::plot_grid(dendro, clust, sampletype, bar + theme(legend.position = "bottom"), ncol=1, rel_heights = c(1.5,0.4,0.4,4), align="v", axis="lr")
  
}


grouped_stacked_bar <- function(ps, n_otus, grouping_var=NULL) {
  
  # select OTUs that you want to use in your plotting
  top <- names(sort(taxa_sums(ps), decreasing = TRUE)[1:n_otus])
  
  # aggregate otu RA's according to group, here type_timepoint, and transform the counts to percentages
  phy_m <- merge_samples(ps, group = grouping_var) %>%
    transform_sample_counts(., function(OTU) OTU/sum(OTU))
  # merge_samples: sums OTUs, in this case relative abundances, by group
  # transform_sample_counts: calculates counts back to relative abundances, but you can also specify a different function
  otu_m <- data.frame(t(otu_table(phy_m)))
  
  # Aggregate RA's of OTUs outside your selection 'top' into 'Residuals' and bind this to the OTU selection  
  other_top <- summarise_each(otu_m[-which(rownames(otu_m) %in% top),], list(~ sum(.)))
  otu_m_top <- rbind(otu_m[which(rownames(otu_m) %in% top),], Residuals = other_top)
  
  # Convert to long format so that you can use it for ggplot
  otu_m_top <- otu_m_top %>% rownames_to_column(var="ASV") #first rownames to column to be able to merge and refer to it
  
  #create a dataframe that is ordered by the abundance of the bacteria - from the top vector as this is ordered
  top_res <- c(top, "Residuals")
  top_res <- factor(top_res, levels=fct_inorder(top_res))
  
  top_ordered <- data.frame(ASV=top_res)
  
  top_df <- left_join(top_ordered, otu_m_top, by="ASV")
  
  library(reshape)
  otum_long_top <- melt(top_df, id.vars = c("ASV"), measure.vars = 2:ncol(top_df), variable_name=grouping_var)
  detach("package:reshape", unload=TRUE)
  otum_long_top <- otum_long_top %>% #rename(RA=value) %>%
    mutate(ASV2=format_OTU(ASV) %>% forcats::fct_inorder() %>% forcats::fct_rev())

  # plot
  top_abu_bar <- ggplot(data = otum_long_top,
                        aes(x = otum_long_top[,2],
                            y = value,
                            fill = ASV2)) +
    geom_bar(stat="identity", color = "white", width = 0.8, position = position_stack()) +   # stat = "identity" makes sure you get stacked bars
    scale_fill_manual(values = c("grey90", rev(make_color_scheme("Paired", n_otus))), 
                      guide = guide_legend(ncol = 2, reverse = T), name="ASV") + 
    scale_y_continuous(expand = c(0, 0), labels = scales::percent) +
    theme_bw(base_size = 13) + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line.x = element_line(color = "grey80"),
          axis.line.y = element_line(color = "grey80"), 
          #axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "right",
          #legend.direction = "horizontal",
          #legend.text = element_text(size = 10),
          #legend.title = element_blank(),
          #legend.key.width = unit(0.4, "cm"),
          #legend.key.height = unit(0.3, "cm"),
          axis.title.x = element_blank(), legend.text = ggtext::element_markdown(),
          legend.key.width=unit(0.6,"cm")) +
    labs(y="Mean relative abundance (%)")
  
  top_abu_bar
}
