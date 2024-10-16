#quant analysis
files = c('astral/dia1.tsv',
          'astral/dia2.tsv',
          'astral/dia3.tsv',
          'astral/dia1-15.tsv',
          'astral/dia2-15.tsv',
          'astral/dia3-15.tsv',
          'astral/dia1-20.tsv',
          'astral/dia2-20.tsv',
          'astral/dia3-20.tsv', 
          "P4778/15mins/Report.tsv",
          "P4778/30mins/Report.tsv",
          "P4778/45mins/Report.tsv",
          
          "P4777/15mins/Report.tsv",
          "P4777/30mins/Report.tsv",
          "P4777/45mins/Report.tsv")

quant_report <- generate_report_quant(files)

report_q <- quant_report %>%
  filter(sample %in% c(103,102,101)) %>%#sample 10% 50% 100%
  mutate(MS_method = paste0(MS, '_', method))
quant_report #quant report data
prots <- qc_prot$Gene.names #qc proteins


# graphs for explanation  -------------------------------------------------



data <- report_q %>%
  filter(MS_method == 'Astral_30_2', PG.Genes == 'GAPDH') %>%
  select(PG.Genes, R.FileName, log2, conc_PRP) %>%
  distinct() %>%
  mutate(conc_PRP_transformed = log(conc_PRP, 10))

# Perform the Pearson correlation test
cor_test <- cor.test(data$conc_PRP_transformed, data$log2, method = "pearson")

# Extract the correlation coefficient and p-value
cor_value <- round(cor_test$estimate, 3)
p_value <- formatC(cor_test$p.value, format = "e", digits = 2)

# Plot the data with the correlation and p-value in the title
ggplot(data, aes(x = conc_PRP_transformed, y = log2)) +
  geom_point(size=2, colour ='red') +
  geom_smooth(method = "lm", se = FALSE, color = "black", size=.5) +
  labs(
    title = paste0(gene, " - Pearson: ", cor_value, ", P-value: ", p_value),
    x = expression('log'[10]*'(Proportion PRP)'),
    y = expression('Protein Intensity (log'[2]*')')
  ) +
  scale_y_continuous(breaks = seq(0,30, by=.5)) + 
  theme  -> p

# Print the plot
print(p)
test <- correlation_analysis(report_q, 'Eclipse_15_3', proteins = prots)

# script for total --------------------------------------------------------

methodss = unique(report_q$MS_method)
generate_correlation_proportion <- function(methods, Report, proteins, shared_genes, fdrs = c(0.05, 0.001, 0.0001)) {
  
  final_result <- data.frame()
  
  for (meth in methods) {
    # Generate correlation analysis for the current method
    tmp <- correlation_analysis(Report, ms_method = meth, proteins = proteins)
    
    # Get the total number of proteins (ids) identified in the correlation analysis
    ids <- nrow(tmp)
    
    # Initialize a vector to store the proportions for each FDR level
    counts_vec <- c()
    
    # Loop over each FDR threshold
    for (fdr_level in fdrs) {
      # Count the number of proteins where the adjusted p-value is less than the FDR threshold
      count_above_threshold <- tmp %>%
        filter(p_adjusted < fdr_level) %>%
        nrow()
      
      # Calculate the proportion of proteins passing the threshold
      prop <- count_above_threshold / ids
      counts_vec <- c(counts_vec, prop)
    }
    
    # Calculate the average datapoints per peak for shared genes in the current method
    avg_datapoints_per_peak <- Report %>%
      filter(MS_method == meth, PG.Genes %in% shared_genes) %>% 
      select(PG.Genes, EG.DatapointsPerPeak) %>%
      summarise(Average_Datapoints_Per_Peak = mean(EG.DatapointsPerPeak, na.rm = TRUE)) %>%
      pull(Average_Datapoints_Per_Peak)
    
    # Combine method name with the proportions for all FDR levels and the average datapoints per peak
    method_result <- data.frame(
      method = meth,
      level1 = counts_vec[1],
      level2 = counts_vec[2],
      level3 = counts_vec[3],
      avg_datapoints_per_peak = avg_datapoints_per_peak  # Add the average datapoints per peak column
    )
    
    # Append the result for the current method to the final result data frame
    final_result <- rbind(final_result, method_result)
  }
  
  return(final_result)
}


prop_report <- generate_correlation_proportion(methodss, report_q, proteins = prots, shared_genes=shared_genes)

scale = 8
prop_report %>%
  separate(method, c('ms', 'grad', 'dia'), remove = FALSE) %>%
  mutate(dia, factor(dia, c(1,2,3))) %>%
  filter(ms == 'Astral') %>%
  
  # Plot with points and lines for level2
  ggplot(aes(x = method, y = avg_datapoints_per_peak)) + 
  geom_line(aes(group = grad, colour = grad)) +
  geom_point(aes(colour = grad), size = 3) +
  
  # Add color scale for gradient
  scale_colour_manual(values = c('#C9AA81', "#808000", "#8A360F"), 
                      labels = c('15 mins', '22 mins', '30 mins'),
                      name = 'Gradient') +
  
  # Add points for avg_datapoints_per_peak with a secondary y-axis
  geom_line(aes(x = method, y = avg_datapoints_per_peak/scale, group=grad)) +  # Dashed line for avg_datapoints_per_peak
  geom_point(aes(x = method, y = avg_datapoints_per_peak/scale), color = "red", size = 2) +
  
  
  # Adding secondary y-axis for avg_datapoints_per_peak
  scale_y_continuous(name = "Proportion",
                     sec.axis = sec_axis(~ .*scale, 
                                         name = "Average Datapoints Per Peak",
                                         breaks=seq(0,10, by=1)), limits=c(0,10)) +  # Secondary axis with actual values
  

  # X-axis labels for methods
  scale_x_discrete(labels = c('A1', 'A2', 'A3', 'A1', 'A2', 'A3', 'A1', 'A2', 'A3'), name = NULL) +
  theme + 
  theme(axis.text.y.right = element_text(colour='red'),
        axis.title.y.right = element_text(colour='red'),
        axis.line.y.right = element_line(colour='red'),
        axis.ticks.y.right = element_line(colour='red'))

prop_report %>%
  separate(method, c('ms', 'grad', 'dia'), remove = FALSE) %>%
  mutate(dia, factor(dia, c(1,2,3))) %>%
  filter(ms == 'TimsTOF') %>%
  
  # Plot with points and lines for level2
  ggplot(aes(x = method, y = avg_datapoints_per_peak)) + 
  geom_line(aes(group = grad, colour = grad), show.legend = F) +
  geom_point(aes(colour = grad), size = 3, show.legend = F) +
  
  # Add color scale for gradient
  scale_colour_manual(values = c('#C9AA81', "#8A360F",'#4B0082'), 
                      labels = c('15 mins', '22 mins', '30 mins'),
                      name = 'Gradient') +
  
  # Adding secondary y-axis for avg_datapoints_per_peak
  scale_y_continuous(name = 'Average Datapoints Per Peak',
                                         
                                         breaks=seq(0,10, by=1), limits=c(0,8)) +  # Secondary axis with actual values
  
  
  # X-axis labels for methods
  scale_x_discrete(labels = c('T1', 'T2', 'T3', 'T1', 'T2', 'T3', 'T1', 'T2', 'T3'), name = NULL) +
  theme 


















# Example usage:
plot_protein_precursor_errorbars(summary_d, ms = "Astral", coeff = 9.3, 
                                 fill_colors = c('#C9AA81', "#808000", "#8A360F"), 
                                 fill_labels = c("22 mins", "30 mins", "30 mins"), 
                                 x_labels = c('A1', 'A2', 'A3', 'A1', 'A2', 'A3', 'A1', 'A2', 'A3'),
                                 show_protein_error = F, 
                                 show_precursor_error = TRUE)

plot_protein_precursor_errorbars(summary, ms = "Eclipse", coeff = 9.3, 
                                 fill_colors = c('#C9AA81', "#8A360F",'#4B0082'), 
                                 fill_labels = c("22 mins", "30 mins", "30 mins"), 
                                 x_labels = c('E1', 'E2', 'E3', 'E1', 'E2', 'E3', 'E1', 'E2', 'E3'),
                                 show_protein_error = F, 
                                 show_precursor_error = TRUE)



# shared genes datapoitns per peak  ---------------------------------------

# Create a list of genes for each method
genes_per_method <- report_q %>%
  filter(PEP.IsProteotypic=='True') %>%
  select(MS_method, PG.Genes) %>%
  distinct() %>%
  group_by(MS_method) %>%
  summarize(Genes = list(PG.Genes))  # Create a list of genes per method

# Find the intersection of genes shared by all methods
shared_genes <- Reduce(intersect, genes_per_method$Genes)
write.csv(shared_genes,'test.csv')
shared_genes %>% filter(shared_genes %in% default_qc_proteins)
avg_datapoints_per_peak <- report_q %>%
  filter(MS_method == 'Astral_30_2', 
         PG.Genes %in% shared_genes) %>% 
  select(PG.Genes, EG.DatapointsPerPeak) %>%
  summarise(Average_Datapoints_Per_Peak = mean(EG.DatapointsPerPeak, na.rm = TRUE)) 

print(avg_datapoints_per_peak)









