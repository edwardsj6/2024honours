#functions

# script read in ----------------------------------------------------------





#read in files, generate cv data group wise (gene level) and count ids
generate_report <- function(files, is_diann=F, map_concs =T) {
  
  ###global report for workspace
  Global_Report <- data.frame()
  print('reading files')
  #stich files together
  for (file in files) {
    data <- read.delim(file, sep = "\t")
    
    Global_Report <- rbind(Global_Report, data)
  }
  
  
  useful_columns_sn <- c('R.FileName', "PG.Genes", "PG.ProteinAccessions", 
                         "PG.Quantity", 'PG.IsSingleHit', 'PEP.IsProteotypic',
                         'EG.PrecursorId','PG.QValue..Run.Wise.', 'PG.PEP..Run.Wise.',
                         'EG.DatapointsPerPeak', 'EG.DatapointsPerPeak..MS1.')
  useful_columns_diann <- c('Run', "Genes", "Protein.Group", 
                            "PG.MaxLFQ", 'Proteotypic',
                            'Stripped.Sequence','PG.Q.Value', 'PEP'
  ) #run specific pep and qval
  print('selecting columns')
  
  if (is_diann==T) {
    Report <- Global_Report[,useful_columns_diann]
    Report <- Report %>%
      rename(
        R.FileName = Run, 
        PG.Genes = Genes, 
        PG.ProteinAccessions = Protein.Group, 
        PG.Quantity = PG.MaxLFQ, 
        PEP.IsProteotypic = Proteotypic,
        EG.PrecursorId = Stripped.Sequence,
        PG.QValue..Run.Wise. = PG.Q.Value,
        PG.PEP..Run.Wise. = PEP
      ) %>%
      mutate(PEP.IsProteotypic = case_when(
        PEP.IsProteotypic == 1 ~ 'True',
        TRUE ~ 'False'
      ))
  } else {
    Report <- Global_Report[,useful_columns_sn]
  }
  
  print('mapping samples')
  #map instruments
  
  #add columns
  Report <- Report %>%
    separate(R.FileName, c('proj', 'item', 'no'), remove =F)
  
  
  
  return(Report) 
}

generate_cvs <- function(data, group_column) {
  # Create a unique ID based on the specified group column and gene
  data <- data %>%
    mutate(ID_gene = paste0(!!sym(group_column), "_", PG.Genes))
  
  # Initialize vectors to store CVs, ID_genes, completeness, and group information
  cv_vector <- c()
  id_gene_vector <- c()
  completeness_vector <- c()
  group_vector <- c()
  
  # Get unique groups from the data
  groups <- unique(data[[group_column]])
  
  # Loop over each group
  for (gp in groups) {
    # Filter data for the current group
    tmp_report <- data %>%
      filter(!!sym(group_column) == gp)
    
    # Calculate the number of unique samples in the current group
    group_length <- length(unique(tmp_report$sample_no))
    
    # Get unique ID_genes for the current group
    ID_genes <- unique(tmp_report$ID_gene)
    
    # Loop over each ID_gene in the current group
    for (ID_gene in ID_genes) {
      # Filter data for the current ID_gene
      data_tmp <- tmp_report %>%
        filter(ID_gene == !!ID_gene)
      
      # Calculate completeness
      completeness <- nrow(data_tmp) / group_length
      
      # Calculate coefficient of variation (CV)
      prot_cv <- sd(data_tmp$log2, na.rm = TRUE) / abs(mean(data_tmp$log2, na.rm = TRUE))
      
      # Append the results to the vectors
      cv_vector <- c(cv_vector, prot_cv)
      id_gene_vector <- c(id_gene_vector, ID_gene)
      completeness_vector <- c(completeness_vector, completeness)
      group_vector <- c(group_vector, gp)
    }
  }
  
  # Create a data frame from the vectors
  cv_map <- data.frame(
    'cv' = cv_vector, 
    'ID_gene' = id_gene_vector, 
    'completeness' = completeness_vector
    
  )
  
  # Merge the CV mapping with the original data
  final_data <- merge(data, cv_map, by = c('ID_gene'), all.x = TRUE)
  
  return(final_data)
}





# script detections -------------------------------------------------------

#create a report summary of unique genes, unfilterde and list(list(filter1), list(filter2)) etc
generate_summary_report <- function(data_in, filter_levels) {
  # Filter proteotypic data
  
  Report_filt <- data_in %>%
    filter(PEP.IsProteotypic == 'True')
  
  # Count total IDs for a method, filtered for proteotypic
  
  Report_filt <- Report_filt
  methods <- unique(Report_filt$group)
  
  methods_vec <- c()
  replicate_vec <- c()
  
  IDstot_vec <- c()
  IDsfilt_vec <- vector("list", length(filter_levels))
  
  avgIDstot_vec <- c()
  avgIDsfilt_vec <- vector("list", length(filter_levels))
  
  sdIDstot_vec <- c()
  sdIDsfilt_vec <- vector("list", length(filter_levels))
  
  avgprec_vec <- c()
  
  for (meth in methods) {
    # Filter data for the specific method and sample group
    data <- Report_filt[Report_filt$group == meth,] 
    
    # Calculate unique genes with replicate information
    dIDs <- calculate_unique_genes(data)
    IDstot_vec <- c(IDstot_vec, dIDs$n_unique_genes)
    avgIDstot_vec <- c(avgIDstot_vec, dIDs$avg_unique_genes)
    sdIDstot_vec <- c(sdIDstot_vec, dIDs$sd_genes)
    replicate_vec <- c(replicate_vec, dIDs$Replicate)
    methods_vec <- c(methods_vec, rep(meth, nrow(dIDs)))
    
    # Loop through filter levels and calculate unique genes for each
    for (i in seq_along(filter_levels)) {
      filters <- filter_levels[[i]]
      
      dIDsfilt <- do.call(calculate_unique_genes, c(list(data), filters))
      
      IDsfilt_vec[[i]] <- c(IDsfilt_vec[[i]], dIDsfilt$n_unique_genes)
      avgIDsfilt_vec[[i]] <- c(avgIDsfilt_vec[[i]], dIDsfilt$avg_unique_genes)
      sdIDsfilt_vec[[i]] <- c(sdIDsfilt_vec[[i]], dIDsfilt$sd_genes)
    }
    
    # Calculate average precursors per replicate
    davgprec <- data %>%
      group_by(R.FileName) %>%
      summarise(n_unique_genes = n_distinct(EG.PrecursorId)) %>%
      mutate(avg_unique_genes = mean(n_unique_genes), 
             sd_genes = sd(n_unique_genes)) %>%
      select(n_unique_genes)
    avgprec_vec <- c(avgprec_vec, davgprec$n_unique_genes)
  }
  
  # Combine the results into a data frame
  plasma_summary <- data.frame(
    method = methods_vec,
    IDstot = IDstot_vec,
    avgIDstot = avgIDstot_vec,
    sdIDstot = sdIDstot_vec,
    avgprec = avgprec_vec,
    replicate = replicate_vec
  )
  
  # Add filtered results for each filter level
  for (i in seq_along(filter_levels)) {
    plasma_summary[[paste0("IDsfilt_level_", i)]] <- IDsfilt_vec[[i]]
    plasma_summary[[paste0("avgIDsfilt_level_", i)]] <- avgIDsfilt_vec[[i]]
    plasma_summary[[paste0("sdIDsfilt_level_", i)]] <- sdIDsfilt_vec[[i]]
  }
  
  return(plasma_summary)
}


calculate_unique_genes <- function(data,
                                   cv_threshold = NULL, 
                                   pep_threshold = NULL,
                                   pg_qval = NULL,
                                   single_hits = NULL) {
  if (!is.null(cv_threshold)) {
    data <- data %>% filter(cv <= cv_threshold)
  }
  if (!is.null(pep_threshold)) {
    data <- data %>% filter(PG.PEP..Run.Wise. <= pep_threshold)
  }
  if (!is.null(pg_qval)) {
    data <- data %>% filter(PG.QValue..Run.Wise. <= pg_qval)
  }
  if (!is.null(single_hits)) {
    data <- data %>% filter(PG.IsSingleHit == single_hits)
  }
  
  result <- data %>%
    group_by(R.FileName) %>%
    summarise(n_unique_genes = n_distinct(PG.Genes)) %>%
    mutate(avg_unique_genes = mean(n_unique_genes), 
           sd_genes = sd(n_unique_genes)) %>%
    select(R.FileName, n_unique_genes, avg_unique_genes, sd_genes)
  result <- result %>%
    mutate(Replicate = sub(".*_(R\\d+).*", "\\1", R.FileName))
  
  return(result)
}


# quantitation analysis ---------------------------------------------------

'''
notes for jack 
orhoganl test for quant would be mixing study wiht different specis (or sample matrixes like hela and plama)
though lcinical proteomics this is as close to the real thing

'''

impute_normal <- function(g1_vec, g1_vec_len, g2_vec_len) {
  set.seed(123)
  
  sd_LOD = lod_sd  # log2
  mu_LOD = lod_mu  # log2
  mu = mu_LOD
  sigma = sd_LOD * 0.3  
  
  imp_g2 <- rnorm(g2_vec_len, mu, sigma)
  
  imp_g1 <- g1_vec 
  g1_vec_missing <- g1_vec_len - length(g1_vec)
  if (g1_vec_missing > 0) {
    imputed_g1 <- rnorm(g1_vec_missing, mean = mu_LOD, sd = sigma)
    imp_g1 <- c(g1_vec, imputed_g1)
  }
  
  return(list(g1_vec = imp_g1, g2_vec = imp_g2))
}




#t-test for two sample groups for a MS-method adjusted p-vals
differential_expression_analysis <- function(data, sample_group1, sample_group2, imp_comp=.7, completeness = .5, test_table = test_lookup) {
  # Filter data based on FDR thresholds for PEP and QVAL
  
  data <- data %>%
    filter(sex %in% c('Female'), #since pool and tech are labelled in sex
           group %in% c(sample_group1, sample_group2),
           
           PEP.IsProteotypic=='True') %>% 
    select(R.FileName, log2, PG.Genes, group, sample_no) %>% 
    distinct()
  
  
  # Ensure sample completeness based on a threshold percentage
  group1_sample_size <- length(unique(data %>% filter(group == sample_group1) %>% pull(sample_no)))
  group2_sample_size <- length(unique(data %>% filter(group == sample_group2) %>% pull(sample_no)))
  
  group1_required <- ceiling(group1_sample_size * completeness)  # Calculate required number of samples for group 1
  group2_required <- ceiling(group2_sample_size * completeness)  # Calculate required number of samples for group 2
  
  genes <- unique(data$PG.Genes)
  
  FCvec <- c()
  pvec <- c()
  genevec <- c()
  impvec <- c()
  
  
  # Find which test to apply for the given sample groups
  test_row <- test_table %>%
    filter((group_1 == sample_group1 & group_2 == sample_group2) | 
             (group_1 == sample_group2 & group_2 == sample_group1)) %>%
    select(test_type)
  test_type <- ifelse(nrow(test_row) > 0, test_row$test_type[1], "Welch's test")
  
  
  for (gene in genes) {
    
    gene_data <- data %>%
      
      filter(PG.Genes == gene)
    
    # Get data for group 1
    gene_data_group1 <- gene_data %>%
      filter(group == sample_group1) %>%
      select(R.FileName, log2, PG.Genes, group) %>% 
      distinct() %>%
      pull(log2)
    
    
    # Get data for group 2
    gene_data_group2 <- gene_data %>%
      filter(group == sample_group2) %>%
      select(R.FileName, log2, PG.Genes, group) %>% 
      distinct() %>%
      pull(log2)
    
    # Perform imputation if data in one group is less than required and completely missing in the other
    imp <- '
    sd_LOD = lod_sd  # log2
    mu_LOD = lod_mu  # log2
    mu = mu_LOD
    sigma = sd_LOD * 0.3  
    
    imp1 <- rnorm((group1_sample_size - length(gene_data_group1)), 
                              mu, sigma)
    gene_data_group1 <- c(gene_data_group1, imp1)
    
    imp2 <- rnorm((group2_sample_size - length(gene_data_group2)), 
                  mu, sigma)
    gene_data_group2 <- c(gene_data_group2, imp2)
    '
    is_imp = F
    if (length(gene_data_group2) == 0 &length(gene_data_group1) >=imp_comp *group1_sample_size) {
      imputed_values <- impute_normal(gene_data_group1, group1_sample_size, group2_sample_size)  # Impute missing values for group 2
      gene_data_group1 <- imputed_values$g1_vec
      gene_data_group2 <- imputed_values$g2_vec
      is_imp = T
    }
    
    if (length(gene_data_group1) == 0 &length(gene_data_group2) >= imp_comp *group2_sample_size  ) {
      imputed_values <- impute_normal(gene_data_group2, group2_sample_size, group1_sample_size)  # Impute missing values for group 1
      gene_data_group1 <- imputed_values$g2_vec
      gene_data_group2 <- imputed_values$g1_vec
      is_imp = T
    }
    
    #If neither group has sufficient samples after imputation, skip
    if (length(gene_data_group1) <= group1_required || length(gene_data_group2) <= group2_required) {
      next
    }
    
    # Conduct appropriate test based on test_lookup
    t <- tryCatch({
      if (test_type == "t-test") {
        t.test(gene_data_group1, gene_data_group2, var.equal = TRUE)  # Perform standard t-test
      } else {
        t.test(gene_data_group1, gene_data_group2, var.equal = FALSE)  # Perform Welch's t-test
      }
    }, error = function(e) {
      return(NULL)
    })
    
    if (is.null(t) || is.na(t$p.value)) {
      next
    }
    
    # Calculate fold change and store results
    FC <- t$estimate[2] - t$estimate[1]
    pvec <- c(pvec, t$p.value)
    FCvec <- c(FCvec, -FC)
    genevec <- c(genevec, gene)
    impvec <- c(impvec, is_imp)
  }
  
  padj <- p.adjust(pvec, method = 'fdr')
  df <- data.frame('PG.Genes' = genevec, 'fc' = FCvec, 'adjp' = padj, 'pval' = pvec, 'imp'= impvec)
  
  return(df)
}


#correlation calculations for sample 101,102,103
correlation_analysis <- function(Report, ms_method, proteins = NULL, samples = c(101, 102, 103)) {
  if (is.null(proteins)) {
    proteins <- unique(Report %>% filter(MS_method == ms_method) %>% pull(PG.Genes))
  }
  
  results <- data.frame()
  expected_data_points <- length(samples) * 3
  
  for (prot in proteins) {
    data <- Report %>%
      filter(MS_method == ms_method, PG.Genes == prot, sample %in% samples) %>%
      select(PG.Genes, R.FileName, log2, conc_PRP, MS_method) %>%
      distinct() %>%
      mutate(conc_PRP_transformed = log(conc_PRP + 0.1, 10))
    
    if (nrow(data) == expected_data_points) {
      x <- data$conc_PRP_transformed
      y <- data$log2
      
      cor_test <- tryCatch({
        cor.test(x, y, method = 'pearson')
      }, error = function(e) {
        return(NULL)
      })
      if (is.null(cor_test)) {
        next
      }
      
      cor_r <- as.numeric(cor_test$estimate)
      cor_pval <- as.numeric(cor_test$p.value)
      
      results <- rbind(results, data.frame(
        PG.Genes = prot,
        correlation = cor_r,
        p_value = cor_pval
      ))
    }
  }
  
  results$p_adjusted <- p.adjust(results$p_value, method = "fdr")
  return(results)
}

#generate cormats so each methods cormat can be quickly called
generate_correlation_matrix <- function(methods) {
  for (meth in methods) {
    cormat <- correlation_analysis(Report, meth)
    
    
    cormat_name <- paste0(meth, '_cormat')
    assign(cormat_name, cormat, envir = .GlobalEnv)
  }
}

#count correlated proteins at threshold
count_corr <- function(tmp_data, fdr_levels) {
  
  methods <- unique(tmp_data$MS_method)
  
  # Initialize an empty list to store results for each FDR level
  results_list <- list()
  
  for (fdr in fdr_levels) {
    vec_props <- c()
    vec_meths <- c()
    
    for (meth in methods) {
      
      var <- paste0(meth, '_cormat')
      data <- get(var, envir = .GlobalEnv)  # Get the corresponding data frame
      
      count <- data %>%
        highlight_qc_proteins() %>%
        filter(is_qc_protein == TRUE) %>%
        nrow()
      
      count_signif <- data %>%
        highlight_qc_proteins() %>%
        filter(is_qc_protein == TRUE, p_adjusted < fdr) %>%
        nrow()
      
      prop_signif <- if (count > 0) count_signif / count else 0
      
      vec_props <- c(vec_props, prop_signif)
      vec_meths <- c(vec_meths, meth)
    }
    
    # Store the results for the current FDR level in a data frame
    result <- data.frame(
      method = vec_meths,
      prop_signif = vec_props,
      fdr_level = rep(fdr, length(vec_meths))
    )
    
    # Append the results to the list
    results_list[[as.character(fdr)]] <- result
  }
  
  # Combine all results into a single data frame
  final_result <- do.call(rbind, results_list)
  
  return(final_result)
  
}

#r-squared analysis
rsquared_analysis <- function(Report, ms_method, proteins = NULL) {
  if (is.null(proteins)) {
    proteins <- unique(Report %>% filter(MS_method == ms_method) %>% pull(PG.Genes))
  }
  
  results <- data.frame()
  
  # Iterate through each protein
  for (prot in proteins) {
    data <- Report %>%
      filter(MS_method == ms_method, PG.Genes == prot) %>%
      select(PG.Genes, R.FileName, log2, conc_PRP, MS_method) %>%
      distinct() %>%
      mutate(conc_PRP_transformed = log(conc_PRP + 0.1, 10))
    
    # Only analyze if there are 12 data points
    if (nrow(data) == 12) {
      x <- data$conc_PRP_transformed
      y <- data$log2
      
      model <- tryCatch({
        lm(y ~ x)
      }, error = function(e) {
        return(NULL)
      })
      if (is.null(model)) {
        next
      }
      
      summary_model <- summary(model)
      
      rsquared <- summary_model$r.squared
      fstat <- summary_model$fstatistic
      p_value <- pf(fstat[1], fstat[2], fstat[3], lower.tail = FALSE)
      
      results <- rbind(results, data.frame(
        PG.Genes = prot,
        R_squared = rsquared,
        p_value = p_value
      ))
    }
  }
  
  # Plotting
  #plot <- ggplot(data, aes(x = conc_PRP_transformed, y = log2)) +
  #geom_point() +
  #geom_smooth(method = "lm", se = FALSE, color = "blue") +
  #ggtitle(paste("Rsquared analysis for", prot, round(rsquared,3), '-- pval', -log(p_value,10))) +
  #xlab(expression("log"[10]*"(conc_PRP + 0.1) ")) +
  #ylab("Log2 Intensity") +
  #theme_minimal()
  
  #print(plot)
  results$p_adjusted <- p.adjust(results$p_value, method = "fdr")
  return(results)
}


# visualisations ----------------------------------------------------------

#directional log log plot for corr, method is vec len=2
log_log_corr <- function(data, meth) {
  meth1 =meth[1] 
  meth2 =meth[2]
  
  cormat1_name <- paste0(meth1, '_cormat')
  corr_df1 <- get(cormat1_name, envir = .GlobalEnv)
  cormat2_name <- paste0(meth2, '_cormat')
  corr_df2 <- get(cormat2_name, envir = .GlobalEnv)
  
  combined_df   <- merge(corr_df1, corr_df2, by='PG.Genes', all=T, suffix = c("_meth1", "_meth2"))
  combined_df <- combined_df %>%
    mutate(meth1 = sign(correlation_meth1) * -log(p_adjusted_meth1,10), 
           meth2 = sign(correlation_meth2) * -log(p_adjusted_meth2,10))
  
  plot <- combined_df %>%
    highlight_qc_proteins() %>%
    arrange(is_qc_protein) %>%
    ggplot(aes(x=meth1, y=meth2)) + 
    geom_point(aes(colour=is_qc_protein)) + 
    scale_colour_manual(values =c( 'grey','red')) +
    scale_x_continuous(breaks=seq(-17, 17, by= 1)) + 
    scale_y_continuous(breaks = seq(-17,17, by=1)) + 
    geom_vline(xintercept = c(-1.3, +1.3), linetype = 2) +
    geom_vline(xintercept = 0)+ 
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = c(-1.3, +1.3), linetype = 2)+
    geom_abline(slope =1) +
    labs(x=meth1, y=meth2, title = paste0(meth1, ' vs ', meth2, " -- Corr")) +
    #geom_text_repel(aes(
    #  label = ifelse(PG.Genes=='GSN', PG.Genes, "") #ifelse(condition, True, False)
    #), size =4, color = "red", max.overlaps = Inf, min.segment.length = 0, box.padding = .5) +
    theme 
  print(plot)
  
  
}


#direction log log plot for FC - samples MUST be specifiec, method is vec len=2
log_log_de <- function(data, sample_u, sample_f) {
  de1 <- differential_expression_analysis(Report, sample_u,sample_f) #
  de2 <- differential_expression_analysis(Report, sample_u,sample_f) #
  combined_df   <- merge(de1, de2, by='PG.Genes', all=T, suffix = c("meth1", "meth2"))
  
  #add corr data
  combined_df <- merge(combined_df, corr_df, all.x=T)
  combined_df <- combined_df %>%
    highlight_qc_proteins() %>%
    arrange(is_qc_protein)
  
  combined_df <- combined_df %>%
    mutate(meth1 = sign(fcmeth1) * -log(adjpmeth1,10), 
           meth2 = sign(fcmeth2) * -log(adjpmeth2,10))
  
  #plotting
  plot <- combined_df %>%
    ggplot(aes(x=meth1, y=meth2)) + 
    geom_point(aes(colour=is_qc_protein)) + 
    scale_colour_manual(values =c('grey', 'red')) +
    scale_x_continuous(breaks=seq(-17, 17, by= 1)) + 
    scale_y_continuous(breaks = seq(-17,17, by=1)) + 
    geom_vline(xintercept = c(-1.3, +1.3), linetype = 2) +
    geom_vline(xintercept = 0)+ 
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = c(-1.3, +1.3), linetype = 2)+
    geom_abline(slope =1) +
    labs(x=meth1, y=meth2, title = paste0(meth1, ' vs ', meth2, " -- FC" )) + 
    theme 
  print(plot)
}


plot <- log_log_corr(Report, c('Eclipse_15_2', 'Astral_30_2'))
"heterogenous smaple (101-104) serached togehte introduce falsae ids for simple samples <-- notes for jack"
# simulations  ------------------------------------------------------------

#genearte correlation plot  
generate_correlation_plot <- function(n = 100, correlation = 0.7) {
  set.seed(42)  # For reproducibility
  
  # Generate random x values
  x <- runif(n, min = 0, max = 1)
  
  # Generate y values based on x and the desired correlation
  y <- correlation * x + sqrt(1 - correlation^2) * runif(n, min = 0, max = 1)
  
  # Create a data frame
  data <- data.frame(x = x, y = y)
  
  # Perform a correlation test
  cor_test <- cor.test(data$x, data$y)
  calculated_correlation <- cor_test$estimate
  
  # Plot
  plot <- ggplot(data, aes(x = x, y = y)) +
    geom_point(color = "blue") +
    geom_smooth(method = "lm", color = "red", se = FALSE) +  # Add a regression line
    labs(title = paste("Desired Correlation =", correlation, 
                       "\nCalculated Correlation =", round(calculated_correlation, 2)),
         x = "X",
         y = "Y") +
    theme_minimal() +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1))
  
  print(plot) 
  
} 
