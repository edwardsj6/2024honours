#IMPUTATION



quant_report <- generate_report_quant(files)
write.csv(quant_report, "quant_report.csv", row.names = FALSE)

# Step 1: Group by EG.PrecursorId and calculate the mean of FG.Quantity
new_report <- quant_report %>%
  group_by(EG.PrecursorId) %>%
  summarise(mean_FG_Quantity = mean(FG.Quantity, na.rm = TRUE)) %>%
  ungroup()

# Step 2: Remove the last two characters of EG.PrecursorId
new_report <- new_report %>%
  mutate(EG.PrecursorId_trimmed = substr(EG.PrecursorId, 1, nchar(EG.PrecursorId) - 2))

# Step 3: Group by the modified EG.PrecursorId and sum the mean FG.Quantity
final_report <- new_report %>%
  group_by(EG.PrecursorId_trimmed) %>%
  summarise(total_FG_Quantity = sum(mean_FG_Quantity, na.rm = TRUE)) %>%
  ungroup()

# View the final report
print(final_report)

log_log_de <- function(data, meth, sample_u, sample_f) {
  
  de1 <- differential_expression_analysis(Report, sample_u[1],sample_u[2], meth) #
  de2 <- differential_expression_analysis(Report, sample_f[1],sample_f[2], meth) #
  combined_df   <- merge(de1, de2, by='EG.PrecursorId', all=T, suffix = c("meth1", "meth2"))
  
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
    labs(x=sample_u, y=sample_f, title = paste0(meth1, ' vs ', meth2, " -- FC" )) + 
    theme 
  #print(plot)
}

differential_expression_analysis <- function(data, sample_group1, sample_group2, ms_method, fdr_threshold = 0.01) {
  # Filter and reshape the data
  sample_data <- data %>%
    filter(sample %in% c(sample_group1, sample_group2), MS_method == ms_method) %>%
    select(EG.PrecursorId, sample, log2, R.FileName) %>%
    distinct() %>% 
    select(EG.PrecursorId, sample, log2)
  
  wide_data <- sample_data %>%
    pivot_wider(names_from = sample, values_from = log2, names_prefix = "n")
  
  genes <- unique(sample_data$EG.PrecursorId)
  
  FCvec <- c()
  pvec <- c()
  genevec <- c()
  
  for (gene in genes) {
    de_group <- wide_data %>% 
      filter(EG.PrecursorId == gene)
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
  df <- data.frame('EG.PrecursorId' = genevec, 'fc' = FCvec, 'adjp' = padj)
  
  return(df)
}

log_log_de(quant_report, 'Astral_20_1', c(104, 101), c(103,101))

# missing value assesment -------------------------------------------------

#check intensityts across astral
Report %>%
  filter(MS == 'Astral', dia == 2) %>%
  select(R.FileName, log2, PG.Genes, sample, method) %>% 
  distinct() %>%
  group_by(PG.Genes) %>%
  mutate(RLE = log2 - median(log2, na.rm = TRUE)) %>%  # Calculate RLE by subtracting the median
  ggplot(aes(x = R.FileName, y = RLE, fill = sample)) + 
  geom_boxplot() + 
  labs(x = "Sample", y = "Relative Log Expression (RLE)") +
  theme

gene_list <- unique(Report$PG.Genes)
methods <- unique(Report$MS_method)

df_missing <- data.frame(genes = gene_list)

for (meth in methods) {
  method_data <- Report %>%
    filter(MS_method == meth, sample %in% c(101))
  missing_vec <-c()
  
  c1 <- length(unique(method_data$sample)) *3
  for (gene in gene_list) {
    samples_for_gene <- method_data %>%
      filter(PG.Genes ==  gene) %>%
      select(PG.Quantity) %>%
      distinct() %>%
      nrow()
      
    missingness = 1 - (samples_for_gene/c1) ####chanege for samples
    missing_vec <- c(missing_vec, missingness)
  }
  
  df_missing <- data.frame(df_missing, 'name' = missing_vec)
  
  
}
colnames(df_missing) <- c('genes', methods)

df_long <- df_missing %>%
  pivot_longer(
    cols = -genes,      # Exclude the 'genes' column from pivoting
    names_to = "method", # Name for the new column that will hold the method names
    values_to = "missingness" # Name for the new column that will hold the missingness values
  )

# Create the box plot
ggplot(df_long, aes(x = method, y = missingness)) +
  geom_boxplot() +
  labs(
    title = "Gene Missingness Across Methods",
    x = "MS Method",
    y = "Proportion of Missingness"
  ) +
  theme_minimal()


# diffrerntauil expression  -----------------------------------------------

#data, sample groups as vector n=2, 
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


# correlations analyasi  --------------------------------------------------


prott = rsquared_analysis(Report, 'Astral_30_2', proteins ='PKM')
print(prott)
Report %>%
  filter(MS_method == 'Astral_30_2', PG.Genes == 'PKM' , sample %in% c(101,102,103,104)) %>%
  select(PG.Genes, R.FileName, log2, conc_PRP, MS_method) %>%
  distinct() %>%
  mutate(xaxis = (log(conc_PRP+.1, 10))) %>%
  ggplot(aes(x= xaxis, y=log2)) + 
  geom_point() + 
  xlim(-1,0.05) +
  labs(x=expression("log"[10]*"(conc_PRP + 0.1) "), 
       title=paste0("Protein = ",prot))+
  theme


#generate correlation function 

correlation_analysis <- function(Report, ms_method, proteins = NULL) {
  if (is.null(proteins)) {
    proteins <- unique(Report %>% filter(MS_method == ms_method) %>% pull(PG.Genes))
  }
  
  results <- data.frame()
  
  #correlation across 100%, 50%, 10% since 0% doesnt contian them
  for (prot in proteins) {
    data <- Report %>%
      filter(MS_method == ms_method, PG.Genes == prot, sample %in% c(101,102,103)) %>%
      select(PG.Genes, R.FileName, log2, conc_PRP, MS_method) %>%
      distinct() %>%
      mutate(conc_PRP_transformed = (log(conc_PRP+.1, 10)+1))
    
    #not complete data -- no corr calc
    if (nrow(data) == 9) {
      x <- data$conc_PRP_transformed
      y <- data$log2

      cor_test <- tryCatch({
        cor_test <- cor.test(x, y, method = 'pearson')
      }, error = function(e) {
        # Skip to the next iteration without printing an error message
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
      # Plotting
      #plot <- ggplot(data, aes(x = conc_PRP, y = log2)) +
        #geom_point() +
        #geom_smooth(method = "lm", se = FALSE, color = "blue") +
        #ggtitle(paste("Correlation for", prot, round(cor_r,2))) +
        #xlab("Concentration (PRP)") +
        #ylab("Log2 Intensity") +
        #theme_minimal()
      
      #print(plot)
    }
  }
  
  results$p_adjusted <- p.adjust(results$p_value, method = "fdr")
  return(results)
}

################### EHERHREHHHERHERHR VVVVVVVVVVVVVVVVVVVVVVVVVVVVV HERHER VVVVVVVVVVVVVVVV
#it owuld be unfair to assess a method on proteins it cant see, as that is a refelction of detection, not quantitation



generate_correlation_matrix <- function(methods) {
  for (meth in methods) {
    cormat <- correlation_analysis(Report, meth, proteins = qc_proteins_all)
    
    
    cormat_name <- paste0(meth, '_cormat_qc')
    assign(cormat_name, cormat, envir = .GlobalEnv)
  }
}

methods = unique(Report$MS_method)

generate_correlation_matrix(methods)

Astral_20_1_cormat %>%
  ggplot(aes(x=correlation, y=p_adjusted)) +
  geom_point()

Report %>% 
  filter(MS_method=='Astral_30_3', PG.Genes %in% c('lacZ')) %>%
  separate(R.FileName, c('proj', 'samp', 'grad', 'dia', 'rep')) %>%
  select(rep, PG.Quantity, log2, conc_PRP) %>%
  ggplot(aes(x=conc_PRP, y = PG.Quantity))+ theme  +geom_point(aes(colour = rep), size =2) 
  
filtered_proteins <- vp %>%
  filter(fc > 0.95 & fc < 1.05 & adjp < 0.05) %>%
  select(PG.Genes)

# Convert to a list if needed
protein_list <- filtered_proteins$PG.Genes

# Print or return the list of proteins
print(protein_list)

# combind it all together -------------------------------------------------

#### %%%%%%%%%%%%
#### Volcano plots
#### %%%%%%%%%%%%

methods <- unique(Report$MS_method)
instruments <- unique(Report$MS)
gradients <- unique(Report$grad) # Assuming you have a Gradient column

# Initialize an empty list to store plots
plot_list <- list()

# Initialize an empty list to store plots
plot_list <- list()

for (ms in instruments) {
  plot_list <- list()
  test <- Report %>%
    filter(MS == ms)
  methods <- unique(test$MS_method)
  for (meth in methods) {
    # Calculate differential expression
    vp <- differential_expression_analysis(Report, 104, 102, meth)
    
    # Highlight QC proteins and calculate the number of QC proteins identified
    vp <- vp %>% highlight_qc_proteins(platelet_proteins)
    
    # Count the number of identified QC proteins
    num_qc_proteins <- sum(vp$is_qc_protein, na.rm = TRUE)
    
    # Add column for lacZ gene identification
    vp <- vp %>% mutate(is_lacz = ifelse(PG.Genes == "lacZ", TRUE, FALSE))
    
    # Create plot for this method and add to list
    plot <- vp %>%
      arrange(is_qc_protein, is_lacz) %>%
      ggplot(aes(x = fc, y = -log10(adjp))) +
      geom_point(aes(color = interaction(is_qc_protein)), show.legend = FALSE) + 
      scale_colour_manual(values = c('grey', 'blue')) +
      xlim(-3, 4) + ylim(0, 6) + 
      geom_vline(xintercept = c(-1), linetype = 'dashed') +
      labs(
        x = expression("log"[2]*"(Fold Change)"), 
        y = expression("-log"[10]*"(adjusted p-value)"),
        title = paste0(meth, " - ", num_qc_proteins, " QC Proteins Identified")
      ) +
      geom_text_repel(aes(
        label = ifelse(is_lacz, paste0(PG.Genes, "-", round(fc,2)), "") #ifelse(condition, True, False)
      ), size = 4, color = "red", max.overlaps = Inf, min.segment.length = 0, box.padding = .5) +
      geom_hline(yintercept = 2, colour = 'black') +
      theme
    
    plot_list[[meth]] <- plot
  }
  t <- grid.arrange(grobs = plot_list, ncol = 3)
  print(t)
}

# Combine plots into a single grid
grid.arrange(grobs = plot_list, ncol = 1)
# Combine plots into a single grid
library(gridExtra)
grid.arrange(grobs = plot_list, ncol = 1)

#### %%%%%%%%%%%%
#### Log log plots
#### %%%%%%%%%%%%

method_pairs = list(
  c('Astral_30_1', 'Astral_30_2'),
  c('Astral_30_1', 'Astral_30_3'),
  c('Astral_30_3', 'Astral_30_2')
  
)


for (meth_pair in method_pairs) {
  print(meth_pair)
  log_log_corr(Report, meth_pair)
  log_log_de(Report, meth_pair, 104,101)
}



# method correlation assesmnt ---------------------------------------------


Astral_30_2_cormat %>%
  highlight_qc_proteins() %>%
  arrange(is_qc_protein) %>%
  ggplot(aes(x=correlation, y=-log(p_adjusted,10))) + 
  geom_point(aes(colour=is_qc_protein)) + 
  scale_colour_manual(values = c('grey','red'))

log_log_de(Report, c('Astral_30_3', 'Astral_30_2'), 101, 103)





# correlation counting prot -----------------------------------------------------------

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

test <- count_corr(Report, c(.05, .01, .001))

#geom point
test %>% #ASTRAL
  separate(method, into = c('ms', 'grad', 'dia'), sep = "_", remove = FALSE) %>%
  mutate(ms_method = paste0(grad, '_', dia),
         fdr_level = factor(fdr_level, c(0.05, 0.01, 0.001))) %>%
  filter(ms == 'Astral') %>%
  
  ggplot(aes(x = ms_method, y = prop_signif, colour = as.factor(fdr_level))) +
  geom_line(aes( group=fdr_level), size=1) +
  scale_x_discrete(labels = c(1,2,3,1,2,3,1,2,3)) +
  geom_point(aes(shape=grad),size=3) + 
  ylim(0,1) +
  theme
  
test %>% #eCLIPSE
  separate(method, into = c('ms', 'grad', 'dia'), sep = "_", remove = FALSE) %>%
  mutate(ms_method = paste0(grad, '_', dia),
         fdr_level = factor(fdr_level, c(0.05, 0.01, 0.001))) %>%
  filter(ms == 'Eclipse') %>%
  
  ggplot(aes(x = ms_method, y = prop_signif, colour = as.factor(fdr_level))) +
  geom_line(aes( group=fdr_level), size=1) +
  scale_x_discrete(labels = c(1,2,3,1,2,3,1,2,3)) +
  geom_point(aes(shape=grad),size=3) + 
  ylim(0,1) +
  theme

test %>% #TimsTOF
  separate(method, into = c('ms', 'grad', 'dia'), sep = "_", remove = FALSE) %>%
  mutate(ms_method = paste0(grad, '_', dia),
         fdr_level = factor(fdr_level, c(0.05, 0.01, 0.001))) %>%
  filter(ms == 'TimsTOF') %>%
  
  ggplot(aes(x = ms_method, y = prop_signif, colour = as.factor(fdr_level))) +
  geom_line(aes( group=fdr_level), size=1) +
  scale_x_discrete(labels = c(1,2,3,1,2,3,1,2,3)) +
  geom_point(aes(shape=grad),size=3) + 
  ylim(0,1) +
  theme

#### bar plot
test %>%
  separate(method, into = c('ms', 'grad', 'dia'), sep = "_", remove = FALSE) %>%
  mutate(ms_method = paste0(grad, '_', dia), 
  fdr_level = factor(fdr_level, c(0.01, 0.001, 0.0001))) %>%
  
  ggplot(aes(x = method, y = prop_signif, fill = ms, alpha = as.factor(fdr_level))) +
  geom_bar(stat = 'identity', position = 'dodge') +
  scale_fill_manual(
    values = c('#D95F02', "#E7298A", "#66A61E" ),
    labels = c("Astral", "Eclipse", "TimsTOF"),
    name = "MS"
  ) + 
  scale_alpha_manual(
    values  = c(0.1, 0.25, 0.5),
    labels = c('5%', "1%", "0.1%"),
    name = "FDR Level"
  ) +
  scale_alpha_manual(values = c(0.25, .5, .9)) + 
  labs(x = "Method", y = "Proportion of Significant QC Proteins", 
       fill = "MS Method", alpha = "FDR Level") +
  theme +
  coord_cartesian(expand =F) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# SIMULATION correlation plot  -------------------------------------------------------


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
generate_correlation_plot(n=100, correlation = 0.990)

