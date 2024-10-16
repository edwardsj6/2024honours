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
  
  print('mapping instruments')
  #map instruments
  Report$Proj <- substring(Report$R.FileName, 1,5)
  MS_map <- data.frame("Proj" = c('P4776', 'P4777', 'P4778'), "MS" = c('Astral', 'TimsTOF', 'Eclipse'))
  Report <- merge(Report, MS_map, by = 'Proj', all.x = TRUE)                    
  
  #add columns
  Report$log2 <- log(Report$PG.Quantity,2)
  Report$sample <- substring(Report$R.FileName, 7,9)
  Report$grad <- substring(Report$R.FileName, 11,12)
  Report$dia <- substring(Report$R.FileName, 14,14)
  Report$method <- substring(Report$R.FileName, 11,14)
  Report$ID <- paste0(Report$MS, '_', substring(Report$R.FileName, 7,14))
  
  
  #generate cvs 
  print('generating cvs')
  Report <- generate_cvs(Report)
  
  
  #conc map 
  if (map_concs == T) {
    print('mapping concs')
    conc_map <- data.frame(sample = c('101', '102','103', '104'), conc_PRP = c(1, .5, .1,0))
    Report <- merge(Report, conc_map, by = 'sample', all.x = TRUE)
  }
  
  
  
  return(Report) 
}

generate_cvs <- function(data) {
  data <- data %>%
    mutate(ID_gene = paste0(ID, "_", PG.Genes))
  
  cv_report <- data         

  #remove dupes from fragements
  cv_report <- cv_report %>% select(R.FileName, PG.Genes, log2, method, ID_gene)
  cv_report$unique_id <- paste(cv_report$R.FileName, cv_report$PG.Genes, cv_report$log2, cv_report$methods)
  cv_report <- cv_report[!duplicated(cv_report$unique_id),]
  
  #init vectors
  t1 <- c()
  t2 <- c()
  t3 <- c()
  
  #calculate CVs
  ID_genes <- unique(cv_report$ID_gene)
  for (ID_gene in ID_genes) {
    #filter data and calc values 
    
    data_tmp  <- cv_report[cv_report$ID_gene == ID_gene,]
    data_tmp <- na.omit(data_tmp)
    completeness <- nrow(data_tmp)/3
    prot_cv <- sd(data_tmp$log2)/mean(data_tmp$log2)
    
    #update vecs 
    t1 <- c(t1, prot_cv)
    t2 <- c(t2, ID_gene)
    t3 <- c(t3, completeness)
  }
  
  #create cv mapping and merge datasets
  cv_map <- data.frame('cv' = t1, 'ID_gene' = t2, 'completeness' = t3)
  data <- merge(data, cv_map, by = 'ID_gene', all.x = TRUE)
  
  return(data) 
}


generate_report_quant <- function(files, is_diann=F, map_concs =T) {
  
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
                         'EG.DatapointsPerPeak', 'EG.DatapointsPerPeak..MS1.',
                         'FG.Quantity', 'FG.IntMID')
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
  
  print('mapping instruments')
  #map instruments
  Report$Proj <- substring(Report$R.FileName, 1,5)
  MS_map <- data.frame("Proj" = c('P4776', 'P4777', 'P4778'), "MS" = c('Astral', 'TimsTOF', 'Eclipse'))
  Report <- merge(Report, MS_map, by = 'Proj', all.x = TRUE)                    
  
  #add columns
  Report$log2 <- log(Report$PG.Quantity,2)
  Report$sample <- substring(Report$R.FileName, 7,9)
  Report$grad <- substring(Report$R.FileName, 11,12)
  Report$dia <- substring(Report$R.FileName, 14,14)
  Report$method <- substring(Report$R.FileName, 11,14)
  Report$ID <- paste0(Report$MS, '_', substring(Report$R.FileName, 7,14))
  
  
 
  
  #conc map 
  if (map_concs == T) {
    print('mapping concs')
    conc_map <- data.frame(sample = c('101', '102','103', '104'), conc_PRP = c(1, .5, .1,0))
    Report <- merge(Report, conc_map, by = 'sample', all.x = TRUE)
  }
  
  
  
  return(Report) 
}


# script detections -------------------------------------------------------

#create a report summary of unique genes, unfilterde and list(list(filter1), list(filter2)) etc
generate_summary_report <- function(data_in, filter_levels) {
  # Filter proteotypic data
  #data_in <- Report_test
  Report_filt <- data_in %>%
    filter(PEP.IsProteotypic == 'True')
  
  # Count total IDs for a method, filtered for proteotypic
  Report_filt$MS_method <- paste0(Report_filt$MS, '_', Report_filt$method)
  Report_filt <- Report_filt
  methods <- unique(Report_filt$MS_method)
  
  methods_vec <- c()
  replicate_vec <- c()
  
  IDstot_vec <- c()
  IDsfilt_vec <- vector("list", length(filter_levels))
  
  avgIDstot_vec <- c()
  avgIDsfilt_vec <- vector("list", length(filter_levels))
  
  sdIDstot_vec <- c()
  sdIDsfilt_vec <- vector("list", length(filter_levels))
  
  avgprec_vec <- c()
  prec_vec <- c()
  
  
  for (meth in methods) {
    # Filter data for the specific method and sample group
    data <- Report_filt[Report_filt$MS_method == meth & Report_filt$sample == 104,] 
    
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
      select(avg_unique_genes, sd_genes, n_unique_genes)
    avgprec_vec <- c(avgprec_vec, davgprec$avg_unique_genes)
    prec_vec <- c(prec_vec, davgprec$n_unique_genes )
  }
  
  # Combine the results into a data frame
  plasma_summary <- data.frame(
    method = methods_vec,
    IDstot = IDstot_vec,
    avgIDstot = avgIDstot_vec,
    sdIDstot = sdIDstot_vec,
    avgprec = avgprec_vec,
    replicate = replicate_vec,
    precs = prec_vec
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
#t-test for two sample groups for a MS-method adjusted p-vals
differential_expression_analysis <- function(data, sample_group1, sample_group2, ms_method, fdr_threshold = 0.01) {
  # Filter and reshape the data
  sample_data <- data %>%
    filter(sample %in% c(sample_group1, sample_group2), MS_method == ms_method) %>%
    select(PG.Genes, sample, log2, R.FileName) %>%
    distinct() %>% 
    select(PG.Genes, sample, log2)
  
  wide_data <- sample_data %>%
    pivot_wider(names_from = sample, values_from = log2, names_prefix = "n")
  
  genes <- unique(sample_data$PG.Genes)
  
  FCvec <- c()
  pvec <- c()
  genevec <- c()
  
  for (gene in genes) {
    de_group <- wide_data %>% 
      filter(PG.Genes == gene)
    mu_before = unlist(de_group[[paste0("n", sample_group1)]]) # x 
    mu_after = unlist(de_group[[paste0("n", sample_group2)]]) # y
    if (is.null(mu_before) || is.null(mu_after) || 
        length(mu_before) < 2 || length(mu_after) < 2) {
      next 
    }
    
    #issues with for loop/ next not working, t.test spits out error. take that
    #as trigger to skip iternation
    t <- tryCatch({
      t.test(mu_before, mu_after, var.equal = TRUE)
    }, error = function(e) {
      # Skip to the next iteration without printing an error message
      return(NULL)
    })
    if (is.null(t)) {
      next
    }
    
    #t <- t.test(mu_before, mu_after, var.equal = TRUE)
    FC <- t$estimate[2] - t$estimate[1]
    pvec <- c(pvec, t$p.value)
    FCvec <- c(FCvec, FC)
    genevec <- c(genevec, gene)
  }
  
  padj <- p.adjust(pvec, method = 'fdr')
  df <- data.frame('PG.Genes' = genevec, 'fc' = FCvec, 'adjp' = padj)
  
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
      mutate(conc_PRP_transformed = log(conc_PRP, 10))
    
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

plot_protein_precursor_errorbars <- function(summary_d, ms = "Astral", coeff = 9.3, fill_colors, fill_labels, x_labels, 
                                             show_protein_error = TRUE, show_precursor_error = TRUE) {
  
  # Generate summary stats for protein data
  summary_stats <- summary_d %>%
    pivot_longer(cols = c(IDstot, IDsfilt_level_1, IDsfilt_level_2, precs),
                 names_to = 'Metric',
                 values_to = 'Values') %>%
    group_by(method, Metric) %>%
    summarise(mean = mean(Values),
              sd = sd(Values),
              se = sd(Values) / sqrt(n()),  # Calculate standard error
              .groups = 'drop') %>%
    separate(method, c('ms', 'grad', 'dia'), remove = F) %>%
    filter(ms == !!ms) %>%
    mutate(Metric = factor(Metric, c('IDstot', 'IDsfilt_level_1', 'IDsfilt_level_2', 'precs')))
  
  # Generate data for plotting precursors
  line <- summary_d %>%
    select(method, avgprec, precs) %>%
    group_by(method) %>%
    summarise(avgprec = mean(precs),
              se_prec = sd(precs) / sqrt(n()),  # Calculate standard error
              .groups = 'drop') %>%
    separate(method, c('ms', 'grad', 'dia'), sep = "_", remove = F) %>%
    filter(ms == !!ms)
  
  # Main plot
  t <- summary_d %>%
    select(method, IDstot, IDsfilt_level_1, IDsfilt_level_2, replicate) %>%
    pivot_longer(cols = c(IDstot, IDsfilt_level_1, IDsfilt_level_2),
                 names_to = 'Metric',
                 values_to = 'Values') %>%
    mutate(Metric = factor(Metric, c('IDstot', 'IDsfilt_level_1', 'IDsfilt_level_2'))) %>%
    separate(method, c('ms', 'grad', 'dia'), sep = "_", remove = F) %>%
    mutate(mm = paste0(method, Metric)) %>%
    group_by(mm) %>% 
    mutate(mean = mean(Values)) %>% 
    filter(ms == !!ms) %>%
    
    ggplot(aes(x = method, y = mean, fill = grad)) + 
    geom_bar(aes(alpha = Metric), colour = 'black', size = 0.5, stat = 'identity', position = 'dodge') + 
    
    geom_point(aes(x = method, y = Values, fill = grad, group = Metric), shape = 21,
               size = 1.5, position = position_jitterdodge(jitter.width = 0.06, dodge.width = 0.9), show.legend = TRUE) +
    
    scale_fill_manual(values = fill_colors,
                      labels = fill_labels,
                      name = "Gradient") + 
    
    scale_alpha_manual(values = c(0.1, 0.25, 0.5),
                       labels = c('Total', "Filtered, FDR < 5%", "Filtered, FDR < 1%"),
                       name = "Protein ID's") +
    
    scale_x_discrete(labels = x_labels, name = NULL) +
    
    labs(title = ms) +
    
    coord_cartesian(expand = FALSE, xlim = c(0.3, 9.7), ylim = c(0, 1200)) +
    
    theme +
    
    geom_line(data = line, aes(x = method, y = avgprec / coeff, group = grad), linewidth = 0.5) + 
    geom_point(data = line, aes(x = method, y = avgprec / coeff), size = 2.5, colour = 'red', show.legend = FALSE) +
    
    scale_y_continuous(name = "Number of Proteins",
                       sec.axis = sec_axis(~.*coeff, name = 'Number of Precursors', breaks = seq(0, 12000, by = 1500)),
                       breaks = seq(0, 1200, by = 100)) + 
    labs() +
    
    theme(axis.text.y.right = element_text(colour = 'red'),
          axis.title.y.right = element_text(colour = 'red'),
          axis.line.y.right = element_line(colour = 'red'),
          axis.ticks.y.right = element_line(colour = 'red'))
  
  # Add protein error bars if requested
  if (show_protein_error) {
    t <- t + 
      geom_errorbar(data = summary_stats, 
                    aes(x = as.numeric(as.factor(method)) + (as.numeric(Metric) - 2) * 0.3, 
                        group = Metric, ymin = mean - se, ymax = mean + se),
                    width = 0.13, size = 0.5)
  }
  
  # Add precursor error bars if requested
  if (show_precursor_error) {
    t <- t + 
      geom_errorbar(data = line, 
                    aes(x = method, 
                        y = avgprec / coeff, 
                        ymin = (avgprec - se_prec) / coeff, 
                        ymax = (avgprec + se_prec) / coeff), 
                    width = 0.08, size = 0.5, show.legend = FALSE)
  }
  
  # Print the plot
  print(t)
}

# Example usage:
plot_protein_precursor_errorbars(summary_d, ms = "Astral", coeff = 9.3, 
                                 fill_colors = c('#D95F02', "#E7298A", "#66A61E"), 
                                 fill_labels = c("15 mins", "30 mins", "45 mins"), 
                                 x_labels = c('A1', 'A2', 'A3', 'A1', 'A2', 'A3', 'A1', 'A2','A3'),
                                 show_protein_error = T, 
                                 show_precursor_error = T)


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
log_log_de <- function(data, meth, sample_u, sample_f) {
  meth1 = meth[1]
  meth2 = meth[2]
  de1 <- differential_expression_analysis(Report, sample_u,sample_f, meth1) #
  de2 <- differential_expression_analysis(Report, sample_u,sample_f, meth2) #
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
