
# formal analysis ---------------------------------------------------------
sn_report_group_sex <- read.delim('sn_report_group_sex.tsv', sep ='\t')
sn_report_group <- read.delim('sn_report_group.tsv', sep ='\t')

theme <- theme_classic() + theme(
  text = element_text((family = "Arial"), size = 12, colour ="black"),  # Set all text to Arial, size 12
  plot.title = element_text((family = "Arial"), size = 14),  # Adjust title size (optional)
  #strip.text.x = element_text((family = "Arial"), size = 12),  # Adjust strip text size (optional)
  axis.text = element_text((family = "Arial"), size = 12, colour ="black"),
  legend.text = element_text((family = "Arial"), size = 12, colour ="black")
  # ... adjust other text elements as needed
)

# read in data ------------------------------------------------------------
report %>% 
  filter(sample_no!=1) %>%
  group_by(sample_no) %>%
  summarise(count = n_distinct(PG.Genes)) %>% 
  ungroup() %>%
  select(sample_no, count) %>%
  distinct() %>% 
  summarise(t = mean(count))
  
file = 'sn_report_groupsNEW.tsv'
group_search <- generate_report(file, map_concs=F)
map <- read.csv('P4954_samplemap.csv')
map <- map %>% 
  rename(no = sample_no)
group_report <-  merge(group_search, map, by='no', all.x=T)

report <- group_report %>%
  rename(sex = Gender, sample_no = no) %>% 
  mutate(log2 = log(PG.Quantity, 2)) %>%
  rename(age = Age.at.time.of.sample.collection ,
         group = Age.at.time.of.sample.collection.1)

# check qcs ---------------------------------------------------------------

colours = c('#a2f2c9', '#f7e07b', '#87ceeb', '#EAB2EA', '#d8e2dc', '#EA8484')
groups = c(0,3,6,9,12,15)
colour_map = data.frame(colours, groups)
#rle

report %>%
  select(sample_no, log2, PG.Genes, group, sex) %>% 
  filter() %>% 
  distinct() %>%
  group_by(PG.Genes) %>%
  mutate(RLE = log2 - median(log2, na.rm = TRUE)) %>%  # Calculate RLE by subtracting the median
  ungroup() %>%
  arrange(group) %>%
  # Reorder sample_no by group
 # mutate(sample_no = factor(sample_no, c(seq(0,96,by=1)))) %>% # fct_reorder(as.factor(as.numeric(sample_no)), sample_no)) %>%
  ggplot(aes(x = sample_no, y = RLE, fill = as.factor(group))) + 
  geom_boxplot(outlier.shape = NA) + 
  labs(x = "Sample Number", y = "Relative Log Expression", fill = "Age Group") +  # Set legend title
  scale_fill_manual(values = c('#a2f2c9', '#f7e07b', '#87ceeb', '#EAB2EA', '#d8e2dc', '#EA8484')) +  # Use a colorblind-friendly palette
  theme +
  scale_x_discrete(breaks = c(1,seq(0,96, by=8)))+
  scale_y_continuous(breaks=seq(-10,10, by=2), limits = c(-4,6)) +
  theme(axis.text.x = element_text( size=12, vjust=.5))


#arrange by gorup 
report %>%
  select(sample_no, log2, PG.Genes, group, sex) %>% 
  filter(sex %in% c('Male', 'Female')) %>%
  distinct() %>%
  group_by(PG.Genes) %>%
  mutate(RLE = log2 - median(log2, na.rm = TRUE)) %>%  # Calculate RLE by subtracting the median
  ungroup() %>%
  arrange(group) %>%
  # Reorder sample_no by group
  mutate(sample_no = factor(sample_no, levels = unique(sample_no[order(group, sample_no)]))) %>%  # Reorder sample_no by group
  ggplot(aes(x = sample_no, y = RLE, fill = as.factor(group))) + 
  geom_boxplot(outlier.shape = NA) + 
  labs(x = "", y = "Relative Log Expression", fill = "Age Group") +  # Set legend title
  scale_fill_manual(values = c('#a2f2c9', '#f7e07b', '#87ceeb', '#EAB2EA', '#d8e2dc', '#EA8484')) +  # Use a colorblind-friendly palette
  theme +
  
  #scale_x_discrete(breaks = c(1, seq(0, 96, by = 8))) +
  scale_y_continuous(breaks = seq(-10, 10, by = 2), limits = c(-4, 6)) +
    theme(axis.text.x = element_blank(),  
          axis.ticks.x = element_blank())



# rle with qc -------------------------------------------------------------
rle_report <- group_report %>%
  rename(sex = Gender, sample_no = no) %>% 
  mutate(log2 = log(PG.Quantity, 2)) %>%
  rename(age = Age.at.time.of.sample.collection ,
         group = Age.at.time.of.sample.collection.1) %>%
  select(sample_no, log2, PG.Genes, sex) %>% 
  filter()


rle_qc <- qc_report %>%
  filter(item %in% c('hela', 'tech')) %>%
  mutate(sample_no = case_when(
    qc_name == 'P4954_hela_1' ~ 0.1,
    qc_name == 'P4954_hela_2' ~ 32.5, 
    qc_name == 'P4954_hela_3' ~ 64.5,
    qc_name == 'P4954_hela_4' ~ 97,
    qc_name == 'P4954_tech_1' ~ 0.2,
    qc_name == 'P4954_tech_2' ~ 32.6, 
    qc_name == 'P4954_tech_3' ~ 64.6,
    qc_name == 'P4954_tech_4' ~ 97.6,
    TRUE ~ NA_real_  
  ),
  sex = item) %>%
  mutate(log2 = log(PG.Quantity,2)) %>%
  select(sample_no, log2, PG.Genes, sex)
 
rle <- rbind(rle_qc, rle_report)

rle %>%
  select(sample_no, log2, PG.Genes, sex) %>%
  mutate(injection = as.numeric(sample_no) + 1) %>%
  mutate(
         groupp = case_when(
         sex %in% c('Male', 'Female') ~ 'sample',
         injection == 2 ~ 'tech',
         sex == 'hela' ~ 'hela',
         sex == 'tech' ~ 'tech', 
         T ~ 'replicate'), 
         ) %>%
  filter() %>% 
  distinct() %>%
  group_by(PG.Genes) %>%
  mutate(RLE = log2 - median(log2, na.rm = TRUE), sample_no = as.numeric(sample_no) + 1) %>%  # Calculate RLE by subtracting the median
  ungroup() %>%
  # Reorder sample_no by group
  mutate(sample_no = factor(sample_no, levels = sort(unique(as.numeric(sample_no))))) %>%
  ggplot(aes(x = sample_no, y = RLE, fill = groupp)) + 
  #ggplot(aes(x = as.character(sample_no), y = RLE, fill =group)) + 
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.8), width = 0.7) + 
  labs(x = "Sample Number", y = "Relative Log Expression", fill = "Age Group") +  # Set legend title
  scale_fill_manual(values = c('#a2f2c9', '#f7e07b', '#87ceeb', '#EAB2EA')) +  # Use a colorblind-friendly palette
  theme +
  #scale_x_discrete(breaks = c(1,seq(0,103, by=8), 103))+
  scale_y_continuous(breaks=seq(-10,10, by=2), limits = c(-4,6)) +
  theme(axis.text.x = element_text( size=12, vjust=.5))

unique(rle$sex)

# check qcs  --------------------------------------------------------------

#beta gal and irts, as well as rle
#lacZ and iRT-Kit_WR_fusion
irt <- report %>%
  mutate(sample_no = factor(sample_no, c(seq(0,96,by=1)))) %>%
  filter(PG.ProteinAccessions %in% c("iRT-Kit_WR_fusion")) %>%
  select(PG.ProteinAccessions, sample_no, log2) %>%
  distinct()
cv_for_irt <- sd(irt$log2)/mean(irt$log2)
bgal <- report %>%
  mutate(sample_no = factor(sample_no, c(seq(0,96,by=1)))) %>%
  filter(PG.ProteinAccessions %in% c("P00722")) %>%
  select(PG.ProteinAccessions, sample_no, log2) %>%
  distinct() 

scale = 1

ggplot() + 
  geom_point(data = bgal, aes(x = sample_no, y = log2), color = 'black') +
  geom_point(data = irt, aes(x = sample_no, y = log2 * scale), color = '#ff7f0e') +  # Example scaling for secondary axis
  scale_y_continuous(name = expression(paste(beta, '-gal Intensity (log'[2]*')')),
  sec.axis = sec_axis(~./scale, name = expression('iRT Kit Intensity (log'[2]*')'))) +  # Adjust scaling for secondary axis
  labs(x = 'Sample') +
  theme + 
  scale_x_discrete(breaks = c(1,seq(0,96, by=16)))+
  coord_cartesian(ylim=c(0,20), xlim =c(-1,97),expand=F) + 
  theme(axis.text.y.right = element_text(colour='#ff7f0e'),
        axis.title.y.right = element_text(colour='#ff7f0e'),
        axis.line.y.right = element_line(colour='#ff7f0e'),
        axis.ticks.y.right = element_line(colour='#ff7f0e')) + 
  geom_vline(xintercept = c(32.5,64.5), linetype='dashed')
 

#hela's, and technical replicates (the frst sample)
#read in mapp  
qc_map <- read.delim('qc_mapping.csv', sep=',')
qc_report <- generate_report('qcs_Report.tsv')

# Original names
old_names <- c("P4954_hela_1", "P4954_hela_2", "P4954_hela_3", "P4954_hela_4", 
               "P4954_sample_1", "P4954_tech_1", "P4954_tech_2", "P4954_tech_3", 
               "P4954_sample_8", "P4954_sample_16", "P4954_sample_24", 
               "P4954_sample_32", "P4954_sample_40", "P4954_sample_48")

new_names <- c("P4954_hela_1", "P4954_hela_2", "P4954_hela_3", "P4954_hela_4", 
               "P4954_tech_1", "P4954_tech_2", "P4954_tech_3", "P4954_tech_4", 
               "P4954_pool_0", "P4954_pool_3", "P4954_pool_6", 
               "P4954_pool_9", "P4954_pool_12", "P4954_pool_15")
name_map <- data.frame(R.FileName = old_names, qc_name = new_names)

qc_report <- merge(qc_report, name_map, by='R.FileName', all.x=T)
qc_report <- qc_report %>% 
  separate(qc_name, c('projj', 'qc', 'qc_no'), remove=F)



# For distinct PG.Genes (proteins)
qc_genes <- qc_report %>%
  filter(qc %in% c('tech', 'hela')) %>%
  select(qc_no, qc, qc_name, PG.Genes) %>%
  distinct() %>%
  group_by(qc_no, qc_name, qc) %>%
  summarise(distinct_count = n_distinct(PG.Genes), .groups = 'drop') %>%
  mutate(type = "Proteins")

# For distinct EG.PrecursorId (precursors)
qc_precursors <- qc_report %>%
  filter(qc %in% c('tech', 'hela')) %>%
  select(qc_no, qc, qc_name, EG.PrecursorId) %>%
  distinct() %>%
  group_by(qc_no, qc_name, qc) %>%
  summarise(distinct_count = n_distinct(EG.PrecursorId), .groups = 'drop') %>%
  mutate(type = "Precursors")

# Combine both data frames
qc_combined <- bind_rows(qc_genes, qc_precursors)

# Plotting
qc_combined %>%
  filter(type == 'Precursors')



qc_precursors_filtered <- qc_report %>%
  filter(qc %in% c('tech', 'hela'), PEP.IsProteotypic=='True', PG.IsSingleHit == 'False') %>%
  select(qc_no, qc, qc_name, PG.Genes) %>%
  distinct() %>%
  group_by(qc_no, qc_name, qc) %>%
  summarise(distinct_count = n_distinct(PG.Genes), .groups = 'drop') %>%
  mutate(type = "Precursors")

# Split into plasma and Hela
plasma_precursors <- qc_precursors_filtered %>% filter(qc == 'tech')
hela_precursors <- qc_precursors_filtered %>% filter(qc == 'hela')
test <- rbind(plasma_precursors, hela_precursors)
test %>%
  ggplot(aes(x=qc_no, y=distinct_count)) +
  geom_point(size=3, aes(colour = qc)) + 
  scale_colour_manual(values=c("#DAA520", '#4169E1'), labels = c('Plasma', 'HeLa'), name ='QC Type') +
  theme
  

plasma_precursors %>%
  ggplot(aes(x=qc_no, y=distinct_count)) + 
  geom_point() + 
  geom_point(data=hela_precursors, aes(x=qc_no, y=distinct_count)) +
  theme

coeff =3
plasma_precursors %>%
  ggplot(aes(x=qc_no, y=distinct_count)) + 
  geom_point(color = "#DAA520", size = 3) + 
  geom_line(aes(group = 1), color = "#DAA520") +
  
  # Add Hela precursors on the secondary y-axis, scaling Hela counts
  geom_point(data = hela_precursors, aes(x = qc_no, y = distinct_count / coeff), color = "#4169E1", size = 3, show.legend=T) +
  geom_line(data = hela_precursors, aes(x = qc_no, y = distinct_count / coeff, group = 1), color = "#4169E1") +
  
  # Primary y-axis for plasma precursors
  scale_y_continuous(name = "Proteins Identified (Plasma)",
                     breaks= seq(0,20000, by=500), 
                       limits = c(0,3500),
                     sec.axis = sec_axis(~ . * coeff,
                                         name = "Proteins Identified (HeLa)")) +
  
  # X-axis labels
  labs(x = "QC Replicate") +
  
  coord_cartesian(expand=F,xlim = c(.7, 4.3)) +
  theme + 
  theme(axis.text.y.right = element_text(colour='#4169E1'),
        axis.title.y.right = element_text(colour='#4169E1'),
        axis.line.y.right = element_line(colour='#4169E1'),
        axis.ticks.y.right = element_line(colour='#4169E1'))



# Combine the two into a single data frame for plotting
combined_precursors <- merge(plasma_precursors, hela_precursors, by = "qc_no", suffixes = c("_plasma", "_hela"))

# Plotting with plasma on the primary y-axis and Hela on the secondary y-axis
combined_precursors %>%
  ggplot(aes(x = qc_no)) +
  geom_line(aes(y = distinct_count_plasma, group = 1, color = "Plasma")) +
  geom_point(aes(y = distinct_count_plasma, color = "Plasma"), size = 3) +
  
  # Add the Hela axis and points
  geom_line(aes(y = distinct_count_hela / max(distinct_count_hela) * max(distinct_count_plasma), group = 1, color = "Hela"), linetype = "dashed") +
  geom_point(aes(y = distinct_count_hela / max(distinct_count_hela) * max(distinct_count_plasma), color = "Hela"), size = 3) +
  
  # Set color scale
  scale_color_manual(values = c("Plasma" = "#1f77b4", "Hela" = "#ff7f0e")) +
  
  # Primary y-axis for plasma
  scale_y_continuous(name = "Distinct Plasma Precursors", 
                     sec.axis = sec_axis(~ . * max(distinct_count_hela) / max(distinct_count_plasma), 
                                         name = "Distinct Hela Precursors")) +
  
  labs(x = "QC No", color = "Sample Type") +
  theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))














ggplot(qc_combined, aes(x = qc_no, y = distinct_count, colour = qc, alpha = type, shape = type)) +
  geom_point(size = 3) +
  geom_line(aes(group = interaction(qc, type), linetype = type), size = .7, show.legend = F) +  # Different linetypes for proteins and precursors
  scale_alpha_manual(values = c(1, .6)) +  # Adjust alpha levels for distinction
  scale_colour_manual(values =c('#4169E1','#DAA520'), labels = c('HeLa', 'Plasma')) + 
  theme_minimal() +
  labs(
       x = "Replicate",
       y = expression('Distinct Counts (log'[10]*')'),
       colour = "QC Type",
       alpha = "QC Type",
       shape = "Metric") +  # Labels for the legend
  theme +
  #scale_y_continuous(breaks = seq(3,6, by=.5) )+
  coord_cartesian( expand=F, xlim=c(0.8, 4.2))


filt = list(
  list(pep_threshold = 0.05, pg_qval = 0.05), 
  list(pep_threshold = 0.01, pg_qval = 0.01))


sn_summary <- generate_summary_report(sn_report_group, filt)

sn_summary %>%
  select(method, replicate, IDstot, IDsfilt_level_1, IDsfilt_level_2) %>%
  pivot_longer(cols = c(IDstot, IDsfilt_level_1, IDsfilt_level_2),
               names_to = 'Metric',
               values_to = 'Values') %>%
  mutate(Metric = factor(Metric, c('IDstot', 'IDsfilt_level_1', 'IDsfilt_level_2'))) %>%
  mutate(mm = paste0(method,Metric)) %>%
  group_by(mm) %>% 
  mutate(mean=mean(Values)) %>% 
  
  ggplot(aes(x=method, y=mean, fill=method)) + 
  geom_bar(aes(alpha=Metric), colour='black', size=.5, stat='identity', position = 'dodge') + 
  geom_point(aes(x=method, y=Values, fill=method, group=Metric), shape=21,
             size=1.5, position=position_jitterdodge(
               jitter.width=.06,dodge.width=.9), show.legend = T) +
  theme


# jumana filtering --------------------------------------------------------
#add completness column

# completeness column  ----------------------------------------------------

sn_report <- read.delim('sn_report_groupsNEW.tsv', sep='\t')

report_filt <- report %>%
  filter(sample_no != 1)

# Filter the data for male and female samples
sn_report_filtered <- report_filt %>%
  filter(sex %in% c('Male', 'Female'))


groups <- unique(sn_report_filtered$group)
completeness_report <- data.frame()

for (i in groups) {
  total_samples <- sn_report_filtered %>%
    filter(group == i) %>%
    select(sample_no) %>%
    distinct() %>%
    nrow()
  
  group_completeness <- sn_report_filtered %>%
    filter(group == i) %>%
    group_by(PG.Genes) %>%
    summarise(
      samples_with_gene = n_distinct(sample_no),
      completeness = samples_with_gene / total_samples
    )

  group_completeness$group <- i
  completeness_report <- bind_rows(completeness_report, group_completeness)
}


sn_report_with_completeness <- report_filt %>%
  left_join(completeness_report, by = c("group", "PG.Genes"))


jf <- sn_report_with_completeness %>%
  filter(PG.IsSingleHit == 'False',
         PEP.IsProteotypic == 'True',
         completeness > 0.667)

#protein groups over time
jf_summary <- jf %>%
  group_by(sample_no) %>%
  mutate(hits = length(unique(PG.Genes))) %>%
  ungroup() 

jf_summary %>%
  arrange(group) %>%
  mutate(sample_no = factor(sample_no, levels = unique(sample_no))) %>%  # Reorder sample_no based on group
  select(sample_no, hits, group, sex) %>%
  distinct() %>%
  mutate(sex_cat = case_when(
    sex %in% c('Male', 'Female') ~ 'one',
    sex %in% c('pool ') ~ 'two', 
    TRUE ~ 'three')) %>%
  ggplot(aes(x = sample_no, y = hits)) + 
  geom_bar(aes(fill = as.factor(group), alpha = sex_cat), stat = 'identity', colour = 'black') + 
  theme + 
  ylim(0, 1000) +
  scale_alpha_manual(values = c(1, .1, .1))+
  #geom_vline(xintercept = c(32.5, 64.5)) +
  labs(x = "Sample No", y = "Hits", fill = "Group")








# pca ---------------------------------------------------------------------


pca_data <- report %>%
  select(sample_no, PG.Genes, log2, age, group, sex) %>%
  filter(sex %in% c('Male', 'Female'), sample_no != 1,
         group == 0) %>%# Select relevant columns
  distinct() %>%  # Remove duplicate rows
  spread(key = PG.Genes, value = log2)  # Convert to wide format

# Ensure rows are samples and columns are features, remove non-numeric columns
pca_matrix <- pca_data %>%
  select(-sample_no, -age, -group, -sex) %>%  # Exclude non-numeric columns
  select_if(~var(.) > 0) %>%  # Keep only columns with non-zero variance
  as.matrix()  # Convert to matrix

pca_result <- prcomp(pca_matrix, scale. = TRUE)
pca_df <- as.data.frame(pca_result$x)  # Get PCA results as a data frame
pca_df$sample_no <- pca_data$sample_no  # Add sample_no for identification
pca_df$age <- pca_data$age  # Add age for coloring
pca_df$group <- pca_data$group  # Add group for coloring
pca_df$sex <- pca_data$sex 

explained_variance <- pca_result$sdev^2 / sum(pca_result$sdev^2) * 100
pc1_variance <- round(explained_variance[1], 1)  # Variance for PC1
pc2_variance <- round(explained_variance[2], 1)  # Variance for PC2

# Update the PCA plot with variance in the axis labels
pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, fill = age*365)) +
  geom_point(size = 3, alpha = 1.4, shape=21, colour='black') +
  labs(
    
    x = paste0("PC1 (", pc1_variance, "%)"),  # PC1 label with variance
    y = paste0("PC2 (", pc2_variance, "%)"),  # PC2 label with variance
    fill = "Age\n(Days)",
    shape = "Sex"
  ) + 
  theme + 
  #scale_fill_manual(values = c('#a2f2c9', '#f7e07b', '#87ceeb', '#EAB2EA', '#d8e2dc', '#EA8484'))
  #scale_fill_manual(values = c('#8700F9', '#00C4AA'))
  scale_fill_continuous(type='viridis')
print(pca_plot)
ggsave('pcaplot_neonatesageindays.svg', plot=last_plot())

# levennes test  ----------------------------------------------------------


library(car)

# Define group levels
group_levels <- c(0, 3, 6, 9, 12, 15)

# Generate unique pairs of groups
group_pairs <- combn(group_levels, 2, simplify = FALSE)

# Initialize a data frame to store the results
test_lookup <- data.frame(
  group_1 = numeric(),
  group_2 = numeric(),
  test_type = character(),
  stringsAsFactors = FALSE
)

# Loop through each pair of groups
for (pair in group_pairs) {
  
  # Filter the data for the specific pair of groups
  data_for_test <- sn_report_group %>%
    filter(sex %in% c('Male', 'Female'),
           group %in% pair) %>%
    select(log2, group)
  
  # Perform Levene's test
  test_result <- leveneTest(log2 ~ as.factor(group), data = data_for_test)
  
  # Determine which test to use based on the p-value of the Levene's test
  p_value <- test_result$`Pr(>F)`[1]  # Extract p-value
  
  # If p-value > 0.05, use t-test; else, use Welch's test
  test_type <- ifelse(p_value > 0.05, "t-test", "Welch's test")
  
  # Add the result to the lookup table
  test_lookup <- rbind(
    test_lookup,
    data.frame(group_1 = pair[1], group_2 = pair[2], test_type = test_type)
  )
}

# View the lookup table
print(test_lookup)

# limit of detection ------------------------------------------------------

methods = seq(2,96, by=1)
limits_of_detection <- c()

# Loop through each method
for (method in methods) {
  print(method)
  # Filter the data for the current method
  method_data <- report %>%
    filter(sample_no == method) %>%
    select(PG.Genes, log2) %>%
    distinct() %>%
    arrange(log2)  # Order by PG.quant in ascending order
  
  # Extract the lowest 3 proteins
  lowest_3 <- method_data %>%
    slice(1:5) %>%  # Select the first 3 rows (lowest intensities)
    pull(log2)  # Get the PG.quant values
  
  # Calculate the average of the lowest 3 intensities
  avg_lowest_5 <- mean(lowest_5, na.rm = TRUE)
  print(avg_lowest_5)
  # Store the result in the limits_of_detection vector
  limits_of_detection <- c(limits_of_detection, avg_lowest_3)
}
lod_mu=mean(limits_of_detection)
lod_sd = sd(limits_of_detection)


# differnetial expression  ------------------------------------------------


library(plotly)




# Define the unique groups
groups <- c(0, 3, 6, 9, 12, 15)

# Function to perform DE analysis and generate interactive volcano plot for each group combination
generate_volcano_plots <- function(jf, groups, imp_comp, completeness, test_lookup) {
  # Get all unique combinations of groups, without doubling up
  group_combinations <- combn(groups, 2, simplify = TRUE)
  
  # Loop over each group combination
  for (i in 1:ncol(group_combinations)) {
    group1 <- group_combinations[1, i]
    group2 <- group_combinations[2, i]
    
    # Perform differential expression analysis for the current group combination
    test <- differential_expression_analysis(jf, 
                                             group1, group2, 
                                             imp_comp = imp_comp,
                                             completeness = completeness, 
                                             test_lookup)
    
    # Create the ggplot volcano plot
    volcano_plot <- test %>%
      mutate(imp = factor(imp, c(T,F))) %>%
      mutate(signif = case_when(
        -log(adjp, 10) > 1.3 ~ T,
        T ~ F)
      ) %>%
      mutate(signif= factor(signif, c(T,F))) %>%
      ggplot(aes(x = fc, y = -log10(pval), text = paste0(PG.Genes, '\n', round(fc,2)))) + 
      geom_point(aes(color = signif, shape=imp), size = 2) +  
      scale_shape_manual(values = c(8,16), labels = c('False', "True"), name = 'Unique\nIdentification')+
      labs(
        title = paste0(group1, ' vs. ', group2),
        x = "Fold Change (log2)",
        y = "-log10(Adjusted p-value)"
      ) +
      theme
    
    # Convert to interactive plot using plotly
    interactive_plot <- ggplotly(volcano_plot, tooltip = "text")
    
    # Save the plot as an HTML file
    html_file_name <- paste0("volcano_plot_", group1, "_vs_", group2, ".html")
    htmlwidgets::saveWidget(interactive_plot, html_file_name)
    
    # Print a message to indicate the plot has been saved
    print(paste("Saved interactive volcano plot for Group", group1, "vs Group", group2))
  }
}

# Example usage of the function
generate_volcano_plots(jf, groups, imp_comp = 0.667, completeness = 0.667, test_lookup = test_lookup)


groups <- c(0, 3, 6, 9, 12, 15)
test <- differential_expression_analysis(jf, 
                                         0, 3, 
                                         imp_comp = 0.667,
                                         completeness = 0.667, 
                                         test_lookup)
pval <- test$pval

adjp <- p.adjust(pval, method='fdr')
test$new <- adjp
test %>%
  mutate(imp = factor(imp, c(T,F))) %>%
  mutate(signif = case_when(
    new < 0.0081 ~ T,
    T ~ F)
  ) %>%
  ggplot(aes(x=fc, y=-log(pval, 10))) + 
  geom_point(aes(colour = signif))+ 
  ylim(2.5,3.5)
write.csv(test %>%
        filter(adjp<0.05) %>%
        pull(PG.Genes), '0v3.csv')

# Create the ggplot volcano plot
test %>%
  mutate(imp = factor(imp, c(T,F))) %>%
  mutate(signif = case_when(
    adjp < 0.05 ~ T,
    T ~ F)
  ) %>%
  mutate(signif= factor(signif, c(T,F))) %>%
  ggplot(aes(x = fc, y = -log10(adjp), text = paste0(PG.Genes, '\n', round(fc,2)))) + 
  geom_point(aes(color = signif, shape=imp), size = 2) +
  scale_colour_manual(values = c('red','grey'), name='Significant', labels = c('True', "False")) +
  scale_shape_manual(values = c(2,16), labels = c('True', "False"), name = 'Unique\nIdentification')+
  scale_x_continuous(breaks=seq(-21,20, by=3)) + 
  scale_y_continuous(breaks=seq(0,20, by=1)) +
  labs(
    #title = paste0( ' vs. '),
    x = expression("Fold Change (log"[2]*")"),
    y = expression("-log"[10]*"(p-value)")
  ) +
  theme

sample_group1 = 0
sample_group2 =3
differential_expression_analysis <- function(data, sample_group1, sample_group2, imp_comp=.7, completeness = .5, test_table = test_lookup) {
  # Filter data based on FDR thresholds for PEP and QVAL
  data <- jf
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
  
  
  
  print(df %>%
          filter(adjp<0.05) %>%
          nrow())
  
  
  # Create the ggplot volcano plot
  print(df %>%
          mutate(imp = factor(imp, c(T,F))) %>%
          mutate(signif = case_when(
            adjp < 0.05 ~ T,
            T ~ F)
          ) %>%
          mutate(signif= factor(signif, c(T,F))) %>%
          ggplot(aes(x = fc, y = -log10(pval), text = paste0(PG.Genes, '\n', round(fc,2)))) + 
          geom_point(aes(color = signif, shape=imp), size = 1.3) +
          scale_colour_manual(values = c('red','grey'), name='Significant', labels = c('True', "False")) +
          scale_shape_manual(values = c(16,16), labels = c('True', "False"), name = 'Unique\nIdentification')+
          scale_x_continuous(breaks=seq(-21,20, by=1)) + 
          scale_y_continuous(breaks=seq(0,20, by=1)) +
          labs(
            title = paste0(group1, ' vs. ', group2),
            x = expression("Fold Change (log"[2]*")"),
            y = expression("-log"[10]*"(p-value)")
          ) +
          theme)
  return(df)
  return(df)
}


####





# Define the unique groups
groups <- c(0, 3, 6, 9, 12, 15)

# Function to perform DE analysis and save FC and adjp for each group combination
save_de_results_to_csv <- function(jf, groups, imp_comp, completeness, test_lookup) {
  group_combinations <- combn(groups, 2, simplify = TRUE)
  
  results_list <- list()
  
  for (i in 1:ncol(group_combinations)) {
    group1 <- group_combinations[1, i]
    group2 <- group_combinations[2, i]
    
    test <- differential_expression_analysis(jf, 
                                             group1, group2, 
                                             imp_comp = imp_comp,
                                             completeness = completeness, 
                                             test_lookup)
    
    result_df <- test %>%
      select(PG.Genes, fc, pval, adjp) %>%
      mutate(Comparison = paste0(group1, "_vs_", group2)) %>%
      rename_with(~ paste0(., "_", group1, "_vs_", group2), c(fc, adjp, pval))
    
   results_list[[i]] <- result_df
  }
  
  merged_results <- Reduce(function(x, y) merge(x, y, by = "PG.Genes", all = TRUE), results_list)
  
  write.csv(merged_results, file = "differential_expression_results_no_imp_female.csv", row.names = FALSE)
  
  print("Saved differential expression results to 'differential_expression_results.csv'.")
}

# Example usage of the function
save_de_results_to_csv(jf, groups, imp_comp = 0.667, completeness = 0.667, test_lookup = test_lookup)





# de between sexes only ---------------------------------------------------

# Function to perform DE analysis and save FC and adjp for each age group (M vs. F)
save_de_results_to_csv <- function(jf, age_groups, imp_comp, completeness, test_lookup) {
  
  results_list <- list()
  
  # Loop over each age group
  for (gp in age_groups) {
    
   test <- differential_expression_analysis(jf, age_group = gp,
                                             imp_comp = imp_comp,
                                             completeness = completeness, 
                                             test_lookup)
    
     result_df <- test %>%
      select(PG.Genes, fc, pval, adjp) %>%
      mutate(Comparison = paste0(age_group, "_Male_vs_Female")) %>%
      rename_with(~ paste0(., "_", age_group, "_Male_vs_Female"), c(fc, adjp, pval))
    
    results_list[[age_group]] <- result_df
  }
  
  merged_results <- Reduce(function(x, y) merge(x, y, by = "PG.Genes", all = TRUE), results_list)
  
  write.csv(merged_results, file = "differential_expression_results_by_sex.csv", row.names = FALSE)
  
  print("Saved differential expression results to 'differential_expression_results_by_sex.csv'.")
}

age_groups = seq(0,15, by=3)
save_de_results_to_csv <- function(jf, age_groups, imp_comp, completeness, test_lookup) {
  
  results_list <- list()
  
  # Loop over each age group
  for (gp in age_groups) {
    
    test <- differential_expression_analysis(jf, age_group = gp,
                                             imp_comp = imp_comp,
                                             completeness = completeness, 
                                             test_lookup)
    
    if (nrow(test) > 0) {  # Ensure the test result has rows
      result_df <- test %>%
        select(PG.Genes, fc, pval, adjp) %>%
        mutate(Comparison = paste0(gp, "_Male_vs_Female")) %>%
        rename_with(~ paste0(., "_", gp, "_Male_vs_Female"), c(fc, adjp, pval))
      
      results_list[[as.character(gp)]] <- result_df
    } else {
      print(paste("No significant results for age group", gp))
    }
  }
  
  if (length(results_list) > 0) {
    # Merge all the data frames by PG.Genes
    merged_results <- Reduce(function(x, y) merge(x, y, by = "PG.Genes", all = TRUE), results_list)
    
    # Save the merged results to a CSV file
    write.csv(merged_results, file = "differential_expression_results_by_sex_imp.csv", row.names = FALSE)
    print("Saved differential expression results to 'differential_expression_results_by_sex.csv'.")
  } else {
    print("No results to save.")
  }
}


save_de_results_to_csv(jf, age_groups, imp_comp = 0.667, completeness = 0.667, test_lookup = test_lookup)



age_group = 15
#t-test for two sample groups for a MS-method adjusted p-vals
differential_expression_analysis <- function(data, age_group, imp_comp=.667, completeness = .667, test_table = test_lookup) {
  # Filter data based on FDR thresholds for PEP and QVAL
  
  
  data <- jf %>%
    filter(sex %in% c('Male','Female'), #since pool and tech are labelled in sex
           group == age_group,
           PG.IsSingleHit == 'False',
           PEP.IsProteotypic=='True') %>% 
    select(R.FileName, log2, PG.Genes, group, sample_no, sex) %>% 
    distinct()
  
  
  # Ensure sample completeness based on a threshold percentage
  group1_sample_size <- length(unique(data %>% filter(sex == 'Male') %>% pull(sample_no)))
  group2_sample_size <- length(unique(data %>% filter(sex == 'Female') %>% pull(sample_no)))
  
  group1_required <- ceiling(group1_sample_size * .667)  # Calculate required number of samples for group 1
  group2_required <- ceiling(group2_sample_size * .667)  # Calculate required number of samples for group 2
  
  genes <- unique(data$PG.Genes)
  
  FCvec <- c()
  pvec <- c()
  genevec <- c()
  impvec <- c()
  
  # Find which test to apply for the given sample groups
  test_type <-  "Welch's test" #everything is welch
  
  
  for (gene in genes) {
    #gene='AMH'
    gene_data <- data %>%
      
      filter(PG.Genes == gene)
    
    # Get data for group 1
    gene_data_group1 <- gene_data %>%
      filter(sex == 'Male') %>%
      select(R.FileName, log2, PG.Genes, group) %>% 
      distinct() %>%
      pull(log2)
    
    
    # Get data for group 2
    gene_data_group2 <- gene_data %>%
      filter(sex == 'Female') %>%
      select(R.FileName, log2, PG.Genes, group) %>% 
      distinct() %>%
      pull(log2)
    
    # Perform imputation if data in one group is less than required and completely missing in the other
    
    is_imp = F
    if (length(gene_data_group2) == 0 &length(gene_data_group1) >=.667 *group1_sample_size) {
      imputed_values <- impute_normal(gene_data_group1, group1_sample_size, group2_sample_size)  # Impute missing values for group 2
      gene_data_group1 <- imputed_values$g1_vec
      gene_data_group2 <- imputed_values$g2_vec
      is_imp = T
    }
    
    if (length(gene_data_group1) == 0 &length(gene_data_group2) >= .667 *group2_sample_size  ) {
      imputed_values <- impute_normal(gene_data_group2, group2_sample_size, group1_sample_size)  # Impute missing values for group 1
      gene_data_group1 <- imputed_values$g2_vec
      gene_data_group2 <- imputed_values$g1_vec
      is_imp = T
    }
    
    #If neither group has sufficient samples  skip
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
  df <- data.frame('PG.Genes' = genevec, 'fc' = FCvec, 'adjp' = padj, 'pval' = pvec, 'imp' = impvec)
  
  
 
  
  print(df %>%
          filter(adjp<0.05) %>%
          nrow())
  
  # Create the ggplot volcano plot
  print(df %>%
    mutate(imp = factor(imp, c(T,F))) %>%
    mutate(signif = case_when(
      adjp < 0.05 ~ T,
      T ~ F)
    ) %>%
    mutate(signif= factor(signif, c(T,F))) %>%
    ggplot(aes(x = fc, y = -log10(pval), text = paste0(PG.Genes, '\n', round(fc,2)))) + 
    geom_point(aes(color = signif, shape=imp), size = 1.3) +
    scale_colour_manual(values = c('red','grey'), name='Significant', labels = c('True', "False")) +
    scale_shape_manual(values = c(16,16), labels = c('True', "False"), name = 'Unique\nIdentification')+
    scale_x_continuous(breaks=seq(-21,20, by=1)) + 
    scale_y_continuous(breaks=seq(0,20, by=1)) +
    labs(
      title = paste0(group1, ' vs. ', group2),
      x = expression("Fold Change (log"[2]*")"),
      y = expression("-log"[10]*"(p-value)")
    ) +
    theme)
  return(df)
  
}

group1 = 15
test <- differential_expression_analysis(15, 
                                         imp_comp = 0.667,
                                         completeness = 0.667, 
                                         test_lookup)

print(test %>%
        filter(adjp<0.05) %>%
        nrow())

# Create the ggplot volcano plot
test %>%
  mutate(imp = factor(imp, c(T,F))) %>%
  mutate(signif = case_when(
    adjp < 0.05 ~ T,
    T ~ F)
  ) %>%
  mutate(signif= factor(signif, c(T,F))) %>%
  ggplot(aes(x = fc, y = -log10(adjp), text = paste0(PG.Genes, '\n', round(fc,2)))) + 
  geom_point(aes(color = signif, shape=imp), size = 1.3) +
  scale_colour_manual(values = c('red','grey'), name='Significant', labels = c('True', "False")) +
  scale_shape_manual(values = c(8,16), labels = c('True', "False"), name = 'Unique\nIdentification')+
  scale_x_continuous(breaks=seq(-21,20, by=3)) + 
  scale_y_continuous(breaks=seq(0,20, by=2)) +
  labs(
    title = paste0(group1, ' vs. ', group2),
    x = expression("Fold Change (log"[2]*")"),
    y = expression("-log"[10]*"(p-value)")
  ) +
  theme





unique_genes <- report %>%
  filter(PEP.IsProteotypic=='True') %>%
  select(PG.Genes) %>%
  distinct() %>% 
  pull(PG.Genes)


samples <- unique(report%>% 
                    filter(sex %in% c('Male', "Female")) %>% 
                    pull(sample_no))



sd_LOD = lod_sd  # log2
mu_LOD = lod_mu  # log2
mu = mu_LOD
sigma = sd_LOD * 0.3 

# Initialize an empty data frame to store the final results
final_df <- data.frame(PG.Genes = character(), 
                       sample_no = integer(), 
                       log2 = numeric(), 
                       IsImputed = character(), 
                       stringsAsFactors = FALSE)


# Loop over each sample number
for (sample in samples)  {
  
  tmp_no <- report %>%
    select(sample_no, log2, PG.Genes) %>%
    filter(sample_no == sample)
  
  # Loop over each unique PG.Gene
  for (gene in unique_genes) {
    tmp_gene <- tmp_no %>%
      filter(PG.Genes == gene) %>%
      select(log2) %>%
      distinct() %>% 
      pull(log2)
    
    is_imp = "No"  # Default to "No" for imputed flag
    
    # Check if log2 value is missing, and impute if necessary
    if (length(tmp_gene) == 0 || !is.numeric(tmp_gene) || (length(tmp_gene) == 1 && is.na(tmp_gene))) {
      tmp_gene <- rnorm(1, mu, sigma)  # Impute the missing value
      is_imp = "Yes"  # Mark as imputed
    }
    
    # Create a temporary data frame to store the results for the current gene and sample
    tmp_df <- data.frame(
      PG.Genes = gene, 
      sample_no = sample, 
      log2 = tmp_gene, 
      IsImputed = is_imp, 
      stringsAsFactors = FALSE
    )
    
    # Append the results to the final data frame
    final_df <- rbind(final_df, tmp_df)
  }
}

# Print or return the final data frame
t <- map%>% 
  rename(group = Age.at.time.of.sample.collection.1, sample_no = no, age = Age.at.time.of.sample.collection) %>%
  select(group, sample_no, Gender, age)
final_df_groups <- merge(final_df, t, all.x= T, by ='sample_no')
write.csv(final_df_groups, 'LS_table.csv')

groups = c(0,3,6,9,12,15)
for (gp in groups) {
  write.table(jf%>%
              filter(PEP.IsProteotypic=="True", completeness==1,
                     group==gp, 
                     PG.Genes != '',
                     PG.Genes !='lacZ')%>% 
              select(PG.Genes) %>% 
              distinct() %>% 
              pull(PG.Genes), paste0(gp, '_prots.csv'), col.names=F, row.names=F, sep=',')
}


# ids per group and total -------------------------------------------------



report %>% #all ids tot 
  group_by(sample_no, group) %>%  # Group by both sample_no and group
  summarise(uniq = n_distinct(PG.Genes), .groups = 'drop') %>%  # Calculate unique PG.Genes
  mutate(sample_no = as.numeric(sample_no)) %>%  # Ensure sample_no is numeric
  ggplot(aes(x = sample_no, y = uniq, fill = as.factor(group))) + 
  geom_bar(stat = 'identity', colour = 'black') + 
  scale_fill_manual(values = c('#a2f2c9', '#f7e07b', '#87ceeb', '#EAB2EA', '#d8e2dc', '#EA8484')) +
  labs(x = "Injection Number", y = "Number Of Proteins") +
  theme+
  scale_x_continuous(breaks = c(1, seq(16,96, by=16))) +
  scale_y_continuous(breaks = seq(0,2000, by=200)) +
  coord_cartesian(expand=F, ylim=c(0,1600))


sn_report_with_completeness %>%
  filter(PEP.IsProteotypic == 'True', PG.IsSingleHit == 'False', 
         completeness > 0.667) %>%
  group_by(sample_no, group) %>%  # Group by both sample_no and group
  summarise(uniq = n_distinct(PG.Genes), .groups = 'drop') %>%  # Calculate unique PG.Genes
  mutate(sample_no = as.numeric(sample_no)) %>%  # Ensure sample_no is numeric
  ggplot(aes(x = sample_no, y = uniq, fill = as.factor(group))) + 
  geom_bar(stat = 'identity', colour = 'black') + 
  scale_fill_manual(values = c('#a2f2c9', '#f7e07b', '#87ceeb', '#EAB2EA', '#d8e2dc', '#EA8484')) +
  labs(x = "Injection Number", y = "Number Of Proteins") +
  theme+
  scale_x_continuous(breaks = c(2, seq(16,96, by=16))) +
  scale_y_continuous(breaks = seq(0,2000, by=200)) +
  coord_cartesian(expand=F, ylim=c(0, 900))




filts = list(list())
report_filt <- sn_report_with_completeness %>% 
  filter(PEP.IsProteotypic == 'True', PG.IsSingleHit == 'False', 
         completeness > 0.667)
summary_report <- generate_summary_report(report, filts)



# Example usage:
plot_protein_precursor_errorbars(summary_report, mss = "Astral", coeff = 10, 
                                 fill_colors = c('#a2f2c9', '#f7e07b', '#87ceeb', '#EAB2EA', '#d8e2dc', '#EA8484'), 
                                 fill_labels = c(), 
                                 x_labels = c(),
                                 show_protein_error = T, 
                                 show_precursor_error = T)
coeff = 12
summary_report %>%
  group_by(method) %>%
  mutate(sem = sd(IDstot) / sqrt(14)) %>%  # Calculate standard error for each method
  mutate(sem_p = sd(avgprec)/sqrt(14),
         avg_p = mean(avgprec)) %>%
  ungroup() %>%
  mutate(method, factor(method, c(0, 3, 6, 9, 12, 15))) %>%
  
  ggplot(aes(x = as.factor(method), y = avgIDstot)) +
  geom_bar(aes(fill = as.factor(method)), alpha = 1, stat = 'identity', colour = 'black', position = 'dodge', show.legend=F) +
  geom_point(aes(y = IDstot, fill =as.factor(method)),
             colour='black', alpha = 1,shape = 21,  size = 2, position = position_jitter(width = 0.2), show.legend = F) + 
  geom_errorbar(aes(ymin = avgIDstot - sem, ymax = avgIDstot + sem), width = 0.2, size = 0.5, colour = 'black') +
  
  scale_fill_manual(values = c('#a2f2c9', '#f7e07b', '#87ceeb', '#EAB2EA', '#d8e2dc', '#EA8484')) +
  labs(x = "Age Group", y = "Number of Proteins") +
  theme + 
  scale_y_continuous(breaks = seq(0,1700, by=200)) +
  coord_cartesian(expand=F, ylim=c(0, 1000), xlim=c(0.3, 6.7))+
  
  
  geom_point(aes(x = as.factor(method), y = avg_p / coeff), size = 2.5, colour = 'red', show.legend = FALSE) +
  
  scale_y_continuous(name = "Number of Proteins",
                     sec.axis = sec_axis(~.*coeff, name = 'Number of Precursors', breaks = seq(0, 12000, by = 1500)),
                     breaks = seq(0, 1200, by = 100)) + 
  theme(axis.text.y.right = element_text(colour = 'red'),
        axis.title.y.right = element_text(colour = 'red'),
        axis.line.y.right = element_line(colour = 'red'),
        axis.ticks.y.right = element_line(colour = 'red')) + 
  

  geom_errorbar( 
                aes(x = as.factor(method), 
                    y = avg_p / coeff, 
                    ymin = (avg_p - sem_p) / coeff, 
                    ymax = (avg_p + sem_p) / coeff), 
                width = 0.08, size = 0.5, show.legend = FALSE)
  
  
  



# sample cohort characteristics  ------------------------------------------

t <- report %>%
  filter(sex %in% c('Male', 'Female')) %>%  # Filter for Male and Female
  select(sample_no, age, group,sex) %>% 
 # distinct() # Select relevant columns
  group_by(group, sex) %>%  # Group by both group and sex
  summarise(mean_age = mean(age, na.rm = TRUE),  # Calculate mean age
            sd_age = sd(age, na.rm = TRUE),      # Calculate standard deviation of age
            .groups = 'drop')



# gene ontology  ----------------------------------------------------------

#read gene scripts in http://127.0.0.1:8837/graphics/plot_zoom_png?width=1070&height=861

write.table( jf%>%
             filter(completeness==1, group ==15) %>%
             select(PG.Genes) %>%
               distinct(), 'prot_15.csv', row.names=F, col.names=F, sep=',')


go_process <- read.delim('go_comp_cluster4.txt')


go_process %>%
  arrange(Client.Text.Box.Input..FDR.) %>% 
  slice(2:13) %>% 
  ggplot(aes(x= Client.Text.Box.Input..fold.Enrichment., y = GO.cellular.component.complete, size=Client.Text.Box.Input..149.)) + 
  geom_point(aes(fill = -log(Client.Text.Box.Input..FDR.,10)), shape = 21, size =5) +
  theme


string <- read.delim('string_functional_annotations.tsv', sep = '\t') 


#file format go_[category]_[age].txt
cats = c('comp', 'process', 'function') 
clusters = c(1,2,3,4,5) 

cats= c( 'process') 
clusters =c(9)
for (cluster in clusters){
  for (cat in cats) {
    file_name = paste0('go_', cat,'_cluster',cluster,'.txt')
  
    go_process <- read.delim(file_name, header=F)
    go <- go_process %>%
      mutate(term = V1, fdr = V8, fold_enrich = V3/V4, count = V3) %>% 
      mutate(
        term = str_replace(term, " \\(GO:\\d+\\)", "")
      ) %>%
      filter(!str_detect(term, "Unclassified"))
    
    gplot <- go %>%
      arrange(fdr) %>%
      slice(1:9) %>%
      arrange(fold_enrich) %>%
      ggplot(aes(
        x = fold_enrich, 
        y = reorder(term, fold_enrich),  # Order GO terms by fold enrichment
        size = count,                     # Dot size represents count
        fill = -log10(fdr)                # Dot color represents -log10(FDR)
      )) + 
      geom_segment(aes(x = 0, 
                       xend = fold_enrich, 
                       y = term, 
                       yend = term, 
                       color = -log10(fdr)), 
                   size = 1.5, show.legend = F) +
      geom_point(shape = 21, color = "black") +
      scale_fill_gradient(low = "blue", high = "red", name = expression('-log'[10]*'(FDR)')) +
      scale_colour_gradient(low = "blue", high = "red") +
      scale_size_continuous(range = c(3, 10), name = "Genes") +
      coord_cartesian() + 
      labs(
        title = file_name,
        y = NULL,
        x = "Fold Enrichment"
        
      ) +
      
      theme
    ggsave(paste0(cat,cluster,".svg"), gplot)
  }
}

mean(sn_report_with_completeness %>%
  filter(sample_no != 1, PEP.IsProteotypic=='True') %>% #, PG.IsSingleHit=='False', completeness >0.667) %>%
  select(EG.PrecursorId, sample_no) %>% 
  distinct()%>%
   
  group_by(sample_no) %>%
  summarise(count = n_distinct(EG.PrecursorId)) %>% 
  pull(count))

# hierachal clustering --------------------------------------------------------------

LS_table <- final_df_groups
wide_data <- LS_table %>%
  filter(PG.Genes != c('', 'lacZ') , sample_no != 1) %>%
  mutate(id = paste0(group,'_',sample_no)) %>%
  pivot_wider(
    id_cols = PG.Genes ,          
    names_from = id,        
    values_from = log2            
  )

wide_data_ag <- wide_data %>%
  rowwise() %>%
  mutate('0' = mean(c_across(starts_with('0_'))),
         '3' = mean(c_across(starts_with('3_'))),
         '6' = mean(c_across(starts_with('6_'))),
         '9' = mean(c_across(starts_with('9_'))),
         '12' = mean(c_across(starts_with('12_'))),
         '15' = mean(c_across(starts_with('15_')))
         ) %>%
  ungroup() %>%
  select(PG.Genes, '0', '3', '6', '9', '12','15')
         

pg_genes <- wide_data_ag$PG.Genes 

data_for_clustering_scaled = as.matrix(scale(wide_data_ag[,-1]))

data_for_clustering <- data_for_clustering_scaled

dist_matrix <- dist(data_for_clustering, method = "euclidean")

hc <- hclust(dist_matrix, method = "average")

plot(hc, labels = pg_genes, main = "Hierarchical Clustering Dendrogram", xlab = "Proteins", sub = "", ylab = "Height")

clusters <- cutree(hc, k = 5)

data_with_clusters <- data.frame(PG.Genes = pg_genes, Cluster = clusters)

# Analyze and view proteins in each cluster
data_with_clusters %>%
  group_by(Cluster) %>%
  summarise(Proteins = paste(PG.Genes, collapse = ", "))

data_for_clustering_scaled$PG.Genes <- pg_genes
wide_data_ag_with_clusters <- wide_data_ag %>%
  left_join(data_with_clusters, by = "PG.Genes")

long_data_ag_with_clusters <- wide_data_ag_with_clusters %>%
  pivot_longer(cols = c('0', '3', '6', '9', '12', '15'),
               names_to = "Time",
               values_to = "log2_intensity") %>%
  mutate(Time = as.numeric(Time))

ggplot(long_data_ag_with_clusters, aes(x = Time, y = log2_intensity, color = as.factor(Cluster.x), group = Cluster.x)) +
  stat_summary(fun = mean, geom = "line", size = 1.2) +  # Plot the average line for each cluster
  labs(title = "Average Expression Profiles by Cluster Over Time", x = "Time", y = "Log2 Intensity", color = "Cluster") +
  theme_minimal()



# Step 1: Scale the data (excluding PG.Genes)
test_scaled <- scale(wide_data_ag[,-1])

# Step 2: Convert the scaled matrix to a data frame
test_scaled_df <- as.data.frame(test_scaled)

# Step 3: Add the PG.Genes column back to the data frame
test_scaled_df$PG.Genes <- pg_genes

# Step 4: Inspect the resulting data frame
head(scaled_data_ag_with_clusters)

scaled_data_ag_with_clusters <- test_scaled_df %>%
  left_join(data_with_clusters, by = "PG.Genes")

long_scaled_data_ag_with_clusters <- scaled_data_ag_with_clusters %>%
  pivot_longer(cols = c('0', '3', '6', '9', '12', '15'),
               names_to = "Time",
               values_to = "scaled_intensity") %>%
  mutate(Time = as.numeric(Time))

ggplot(long_scaled_data_ag_with_clusters, aes(x = Time, y = scaled_intensity, color = as.factor(Cluster), group = Cluster)) +
  stat_summary(fun = mean, geom = "line", size = 1.2) +  # Plot the average line for each cluster
  labs(title = "Average Scaled Expression Profiles by Cluster Over Time", x = "Time", y = "Scaled Intensity", color = "Cluster") +
  theme






long_data_ag <- data_with_clusters %>%
  gather(key = "Time", value = "log2_intensity", -PG.Genes, -Cluster) %>%
  mutate(Time = as.numeric(gsub("Time", "", Time)))

ggplot(long_data_ag, aes(x = Time, y = log2_intensity, color = as.factor(Cluster), group = Cluster)) +
  stat_summary(fun = mean, geom = "line", size = 1.2) +  # Plot the average line for each cluster
  labs(title = "Average Expression Profiles by Cluster", x = "Time", y = "Log2 Intensity", color = "Cluster") +
  theme



#plots the WSSS
# Function to calculate total within-cluster variance (sum of squares)
total_within_cluster_variance <- function(hc, k, data) {
  clusters <- cutree(hc, k = k)
  total_ss <- 0
  
  for (i in 1:k) {
    cluster_data <- data[clusters == i, , drop = FALSE]
    cluster_center <- colMeans(cluster_data)
    total_ss <- total_ss + sum(rowSums((cluster_data - cluster_center) ^ 2))
  }
  
  return(total_ss)
}

# Calculate total within-cluster variance for different numbers of clusters
k_values <- 1:10  # Try different numbers of clusters
wss_values <- sapply(k_values, function(k) total_within_cluster_variance(hc, k, data_for_clustering))

# Plot the WSS values (to find the elbow point)
wss_df <- data.frame(k = k_values, wss = wss_values)

ggplot(wss_df, aes(x = k, y = wss)) +
  geom_line() +
  geom_point() +
  scale_x_continuous(breaks=seq(0,10,by=1)) + 
  labs(x = "Number of Clusters", y = "Total Within-Cluster Variance", title = "Elbow Method for Hierarchical Clustering") +
  theme

# clustering, but relateive to zero ---------------------------------------

# Calculate differences from time 0 for each time point

LS_table <- final_df_groups
wide_data <- LS_table %>%
  filter(PG.Genes != c('', 'lacZ') , sample_no != 1) %>%
  mutate(id = paste0(group,'_',sample_no)) %>%
  pivot_wider(
    id_cols = PG.Genes ,          
    names_from = id,        
    values_from = log2            
  )

wide_data_ag <- wide_data %>%
  rowwise() %>%
  mutate('0' = mean(c_across(starts_with('0_'))),
         '3' = mean(c_across(starts_with('3_'))),
         '6' = mean(c_across(starts_with('6_'))),
         '9' = mean(c_across(starts_with('9_'))),
         '12' = mean(c_across(starts_with('12_'))),
         '15' = mean(c_across(starts_with('15_')))
  ) %>%
  ungroup() %>%
  select(PG.Genes, '0', '3', '6', '9', '12','15')


pg_genes <- wide_data_ag$PG.Genes 

wide_data_ag_diff <- wide_data_ag %>%
  rowwise() %>%
  mutate(  
          
         `6_diff` = `6` - `3`,
         `9_diff` = `9` - `3`,
         `12_diff` = `12` - `3`,
         `15_diff` = `15` - `3`,
         `3` = 3) %>%
  ungroup() %>%
  select(PG.Genes,`3`, `6_diff`, `9_diff`, `12_diff`, `15_diff`)

# Scaling the difference data
data_for_clustering_scaled_diff <- as.matrix(scale(wide_data_ag_diff[,-1]))

# Perform clustering on the scaled difference data
dist_matrix_diff <- dist(data_for_clustering_scaled_diff, method = "euclidean")
hc_diff <- hclust(dist_matrix_diff, method = "ward.D")

# Plot dendrogram
plot(hc_diff, labels = pg_genes, main = "Hierarchical Clustering Dendrogram (Differences)", xlab = "Proteins", sub = "", ylab = "Height")

for (i in 3:8) {
  # Cut into clusters
  clusters_diff <- cutree(hc_diff, k = 5)
  
  # Compute WSS for this clustering
  total_withinss <- 0
  for (j in 1:i) {
    cluster_data <- data_for_clustering_scaled_diff[clusters_diff == j, , drop = FALSE]
    if (nrow(cluster_data) > 0) {
      # Remove rows with NAs
      cluster_data <- cluster_data[complete.cases(cluster_data), , drop = FALSE]
      if (nrow(cluster_data) > 0) {
        cluster_center <- colMeans(cluster_data)
        ss <- sum(rowSums((cluster_data - cluster_center)^2))
        total_withinss <- total_withinss + ss
      } else {
        warning(paste("Cluster", j, "has no data after removing NAs."))
      }
    } else {
      warning(paste("Cluster", j, "is empty."))
    }
  }
  
  
  # Assign clusters and add to the difference data
  data_with_clusters_diff <- data.frame(PG.Genes = pg_genes, Cluster = clusters_diff)
  
  
  # Long format for plotting
  long_data_ag_with_clusters_diff <- wide_data_ag_diff %>%
    left_join(data_with_clusters_diff, by = "PG.Genes") %>%
    pivot_longer(cols = c(`3`, `6_diff`, `9_diff`, `12_diff`, `15_diff`),
                 names_to = "Time",
                 values_to = "log2_intensity_diff") %>%
    mutate(Time = as.numeric(gsub("_diff", "", Time)))
  
  # Plot the average profiles by cluster for the differences
  t <- ggplot(long_data_ag_with_clusters_diff, aes(x = Time, y = log2_intensity_diff, color = as.factor(Cluster), group = Cluster)) +
    stat_summary(fun = mean, geom = "line", size = 1.2) +  # Plot the average line for each cluster
    labs(title = "Average Difference Profiles by Cluster Over Time", x = "Time", y = "Log2 Intensity Difference", color = "Cluster") +
    theme
  print(t)
}


long_data_ag_with_clusters_diff %>%
  group_by(Cluster) %>%
  summarise(count = n_distinct(PG.Genes))




wide_data_ag_diff <- wide_data_ag %>%
  rowwise() %>%
  mutate(
    `3_diff` = `3` - `0`,  # Differences from time 0
    `6_diff` = `6` - `0`,
    `9_diff` = `9` - `0`,
    `12_diff` = `12` - `0`,
    `15_diff` = `15` - `0`
  ) %>%
  ungroup() %>%
  select(PG.Genes, `3_diff`, `6_diff`, `9_diff`, `12_diff`, `15_diff`)

# Scaling the difference data (excluding the first column)
data_for_clustering_scaled_diff <- as.matrix(scale(wide_data_ag_diff[,-1]))


# Function to calculate total within-cluster variance for hierarchical clustering
total_within_cluster_variance <- function(hc, k, dataa) {
  # Cut the dendrogram at k clusters
  clusters <- cutree(hc, k)
  
  # Initialize the total within-cluster variance
  total_variance <- 0
  
  # Loop through each cluster and calculate the within-cluster variance
  for (i in 1:k) {
    cluster_data <- dataa[clusters == i, ]
    if (nrow(cluster_data) > 1) {
      cluster_center <- colMeans(cluster_data)
      total_variance <- total_variance + sum(rowSums((cluster_data - cluster_center)^2))
    }
  }
  
  return(total_variance)
}

# Calculate total within-cluster variance for different numbers of clusters
k_values <- 1:10  # Try different numbers of clusters
wss_values <- sapply(k_values, function(k) total_within_cluster_variance(hc_diff, k, data_for_clustering_scaled_diff))

# Plot the WSS values (to find the elbow point)
wss_df <- data.frame(k = k_values, wss = wss_values)

# Plot the elbow plot
ggplot(wss_df, aes(x = k, y = wss)) +
  geom_line() +
  geom_point() +
  scale_x_continuous(breaks = seq(1, 10, by = 1)) + 
  labs(x = "Number of Clusters", y = "Total Within-Cluster Variance", title = "Elbow Method for Hierarchical Clustering") +
  theme



# upsets ------------------------------------------------------------------


library(ggupset)

# Filter the data based on specified criteria
upsetss <- sn_report_with_completeness %>%
  filter(PEP.IsProteotypic == 'True', PG.IsSingleHit == 'False') %>%
  filter(completeness > .667)

# Extract unique PG.Genes per group
gene_list <- upsetss %>%
  mutate(group = factor(group, levels = c('0','3','6','9','12','15')))%>%
  select(group, PG.Genes) %>%
  distinct() 

# Prepare the data for ggupset by creating a list-column
upset_data <- gene_list %>%
  group_by(PG.Genes) %>%
  summarize(groups = list(as.character(group))) 

# Create the UpSet plot
upset_data %>%
  ggplot(aes(x=groups)) +
  geom_bar(colour='black', fill='#FFFDD0') +
  scale_x_upset(n_intersections = 10) +  # Creates the combination matrix for the UpSet plot
  geom_text(stat='count', aes(label=after_stat(count)), vjust=-1) +
  theme +
  labs(
   
    x = "Groups",
    y = "Number of Proteins"
  )  


write.csv(length(upset_data %>% 
  filter(groups == '0') %>% 
  distinct() %>%
  pull(PG.Genes)), 'test.csv')


# missingness boxplot -----------------------------------------------------------------
pg_filt = unique(jf$PG.Genes)

missing <- final_df_groups %>%
  #filter(PG.Genes %in% pg_filt) %>%
  select(group, PG.Genes, IsImputed, sample_no) %>%
  distinct() %>%
  group_by(group, PG.Genes) %>%
  mutate(missingness = sum(IsImputed == 'Yes') / n()) %>%
  ungroup()

missing %>% 
  ggplot(aes(x=as.factor(group), y=missingness*100)) + 
  geom_boxplot(aes(fill=as.factor(group)), outliers = F, show.legend=F) + 
  ylim(0,100)+
  scale_fill_manual(values = c('#a2f2c9', '#f7e07b', '#87ceeb', '#EAB2EA', '#d8e2dc', '#EA8484')) +
  theme + 
  labs(y='Missingness (%)', x='Age Group')
 


# CVs violin plot ---------------------------------------------------------

jf #filterd dataset and only clinical smapoles

cv_report <- jf %>%
  select(PG.Genes, group, sample_no, log2) %>% 
  distinct() %>%
  group_by(group, PG.Genes) %>%
  mutate(cv = sd(log2)/mean(log2)) %>%
  ungroup()

cv_report %>%
  ggplot(aes(x=as.factor(group), y=cv*100)) + 
  geom_violin(aes(fill=as.factor(group)), show.legend=F)+
  scale_fill_manual(values = c('#a2f2c9', '#f7e07b', '#87ceeb', '#EAB2EA', '#d8e2dc', '#EA8484')) +
  labs(x='Age Group', y='CV (%)') +
  theme


# box  plot for protein ---------------------------------------------------
final_df_groups
library(ggpubr)

prot = 'HAP1'
comps = list(c('0', '3'), c('0', '15'), c('0','9'))

tmp <- final_df_groups %>%
  filter(PG.Genes == prot) %>% 
  select(PG.Genes, log2, sample_no, group) %>%
  distinct() %>%
  ggplot(aes(x = as.factor(group), y = log2)) +
  geom_boxplot(aes(fill=as.factor(group)), show.legend = F) + 
  stat_compare_means(label = 'p.format', 
                     comparisons = comps, 
                     method = 't.test', 
                     symnum.args = list(cutpoints = c(0,0.001, 0.01, 0.05, Inf), symbols = c('***', '**', '*', 'ns'))) +
  theme + 
  scale_fill_manual(values =c('#a2f2c9', '#f7e07b', '#87ceeb', '#EAB2EA', '#d8e2dc', '#EA8484')) +
  labs(y=expression('Protein Intensity (log'[2]*')'), x='Age Group')

# Print the plot
print(tmp)


plot_protein_comparisons <- function(protein, comparisons, sexes = c('Male', 'Female'), yl = NULL) {
  data <- final_df_groups
  tmp <- data %>%
    filter(PG.Genes == protein, Gender %in% sexes) %>% 
    select(PG.Genes, log2, sample_no, group) %>%
    distinct()
  
  plot <- ggplot(tmp, aes(x = as.factor(group), y = log2)) +
    geom_boxplot(fill = 'white', show.legend = F)+ 
    stat_compare_means(label = 'p.format', 
                       comparisons = comparisons, 
                       method = 't.test', 
                       symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, Inf), symbols = c('***', '**', '*', 'ns')),
                       #vjust=.5, 
                       step.increase=0.1) +
    #scale_fill_manual(values = c('#a2f2c9', '#f7e07b', '#87ceeb', '#EAB2EA', '#d8e2dc', '#EA8484')) +
    labs(y = expression('Protein Intensity (log'[2]*')'), x = 'Age Group', title = protein) +
    scale_y_continuous(limits = yl) +
    theme
  
  return(plot)
}

plot_protein_comparisons_sex <- function(protein, comparisons, yl = NULL) {
  data <- final_df_groups
  tmp <- data %>%
    filter(PG.Genes == protein, Gender == 'Female') %>% 
    select(PG.Genes, log2, sample_no, group,Gender) %>%
    distinct()
  
  plot <- ggplot(tmp, aes(x = as.factor(group), y = log2)) +
    geom_boxplot(aes(fill = as.factor(Gender)), show.legend = F) + 
    stat_compare_means(label = 'p.format', 
                       comparisons = comparisons, 
                       method = 't.test', 
                       symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, Inf), symbols = c('***', '**', '*', 'ns')),
                       #vjust=.5, 
                       ) +
    scale_fill_manual(values = c('#8700F9')) +
    labs(y = expression('Protein Intensity (log'[2]*')'), x = 'Age Group', title = paste0(protein, " - Female")) +
    scale_y_continuous(limits = yl) +
    theme
  print(plot)
  
  tmp <- data %>%
    filter(PG.Genes == protein, Gender =='Male') %>% 
    select(PG.Genes, log2, sample_no, group,Gender) %>%
    distinct() 
  plot <- ggplot(tmp, aes(x = as.factor(group), y = log2)) +
    geom_boxplot(aes(fill = as.factor(Gender)), show.legend = F) + 
    stat_compare_means(label = 'p.format', 
                       comparisons = comparisons, 
                       method = 't.test', 
                       symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, Inf), symbols = c('***', '**', '*', 'ns')),
                       #vjust=.5, 
                       ) +
    scale_fill_manual(values = c( '#00C4AA')) +
    labs(y = expression('Protein Intensity (log'[2]*')'), x = 'Age Group', title = paste0(protein, " - Male")) +
    scale_y_continuous(limits = yl) +
    theme
  print(plot)
  
}


# plot sex proteins -------------------------------------------------------

plot_protein_comparisons_sex('SPP1',comparisons=list(c('12','15'),c('0','12'),c('0','15')), c(6,11.5))
plot_protein_comparisons_sex('DPP4',list(c('12','15'),c('0','12'),c('0','15')), c(7.5,10.5))
plot_protein_comparisons_sex('CHL1',list(c('6','15'),c('0','6'),c('0','15')), c(7.5,10.5))
plot_protein_comparisons_sex('AMH',list(c('3','15'),c('0','3'),c('0','15')), c(1,10))

#while hormbonoe bidnign adn trasnport is sex specific, porteins have similar levels bewteen groups
plot_protein_comparisons('SHBG',list(c('12','15'),c('3','12'),c('0','3')))
plot_protein_comparisons('SPP2',list(c('12','15'),c('3','12'),c('0','3')))

# anova clusters ----------------------------------------------------------

first <- c('HAP1', 'USF3')
plot_protein_comparisons('HAP1',list(c('12','15'), c('3','9'), c('12','3'), c('0','12'))) #linkeind to n
plot_protein_comparisons('USF3',list(c('0','3'), c('3','9'), c('9','15')))

second <- 'IGFBP1'
plot_protein_comparisons('IGFBP1',list(c('0','3'), c('3','9'), c('9','15')))

third <- 'PM20D1'
plot_protein_comparisons('PM20D1',list(c('0','3'), c('3','9'), c('9','15')))

library(ggpubr)
#pd paper prots
pd <- c('PPP3CB', 'DKK3', 'HSPA5', 'SERPINF2', 'SERPINA3', 'SERPING1', 'C3', 'MASP2', 'GRN')
pd <- c('GRN', 'DKK3', 'C3', 'SERPINA3', 'HPX', 'SERPINF2', 'CAPN2', 'SERPING1', 'SELE')
for (prot in pd) {
  print(plot_protein_comparisons(prot,list(c('0','3'), c('3','9'), c('9','15'))))
}

#some other intersing ones
unique_prots <- c('CDH1', 'CDH6', 'IGF1')
for (prot in unique_prots) {
  print(plot_protein_comparisons(prot,list(c('0','3'), c('3','9'), c('3','15'))))
}

signif_genes


# immune prots ------------------------------------------------------------


immune_prots <- c('IGHV1-18', 'IGHV1-24', 'IGHV1-46', 'IGHV1-58', 'IGHV1-69', 'IGHV3-38', 'IGHV4-34', 'IGHV5-51',
                  'IGLV1-40', 'IGLV3-1', 'IGLV3-19', 'IGLV4-69')
antibody_immune_prots <- c( 'IGHV3-38','IGLV3-1','IGLV4-69', 'IGHG1')
for (prot in antibody_immune_prots) {
  print(plot_protein_comparisons(prot,list(c('0','3'), c('3','6'),c('3','15'), c('6','15'), c('0','15'))))
}
''
immune_prots <- c('IGKC', 'CFD', 'IGHG1')
for (prot in immune_prots) {
  print(plot_protein_comparisons(prot,list(c('0','3'), c('3','9'), c('0','15'))))
}

comp_immpune_prots <- c('CFD', 'SERPING1', 'C3' )
for (prot in comp_immpune_prots) {
  print(plot_protein_comparisons(prot,list(c('0','3'), c('3','9'), c('0','15'))))
}
print(plot_protein_comparisons('PTGDS',list(c('0','3'), c('3','9'), c('0','15'),c('0','9'))))

inflamation_proteins <- c('SERPINA5', 'SERPINA6', 'PTGDS', 'APCS', 'A2M')
inflamation_proteins <- c( 'SERPINA6', 'PTGDS')
for (prot in inflamation_proteins) {
  print(plot_protein_comparisons(prot,list(c('0','3'), c('3','9'), c('0','15'))))
}

complement_system <- c('C3', 'C5', 'C9', 'CFD', 'CFH', 'CFHR5')
complement_system <- c('C5', 'C9', 'CFD',  'CFHR5')
for (prot in complement_system) {
  print(plot_protein_comparisons(prot,list(c('0','3'), c('3','6'),c('3','15'),  c('0','15'))))
}

other_immune_proteins <- c('SERPINA6', 'PTGDS', 'SELL', 'IGFBP3')
for (prot in other_immune_proteins) {
  print(plot_protein_comparisons(prot,list(c('0','3'), c('3','6'),c('3','15'), c('6','15'), c('0','15'))))
}
# plot heamoproteins ------------------------------------------------------



plot_protein_comparisons('HAP1',list(c('0','3'), c('3','12'), c('3','15'), c('0','15')))

chantal <- c('F5', 'F2', 'F12', 'F10', 'PROS1','PROC','PROCR', 'F9', 'F11', 'F13A1','F13B', 'F8')
for (prot in chantal) {
  print(plot_protein_comparisons(prot,list(c('0','3'), c('3','15'))))
  
}

coag_factors <- c('FGA', 'FGB', 'FGG')
for (prot in coag_factors) {
  print(plot_protein_comparisons(prot,list(c('0','3'), c('3','15'), c('0','15'))))
}
plot_protein_comparisons('VWF',list(c('0','3'), c('3','12'), c('3','15'), c('0','15')))

clot_formation <-c('PROC', 'PROS1', 'F5', 'F8')
for (prot in clot_formation) {
  print(plot_protein_comparisons(prot,list(c('0','3'), c('3','15'), c('0','15'))))
}
print(plot_protein_comparisons('F8',list(c('0','3'), c('3','15'), c('0','15'))))

clot_breakdown <- c('PLG', 'SERPINF2')
for (prot in clot_breakdown) {
  print(plot_protein_comparisons(prot,list(c('0','3'), c('3','15'), c('0','15'))))
}

thrombins <- c('THBS1','THBS2','THBS3','THBS4')
for (prot in thrombins) {
  print(plot_protein_comparisons(prot,list(c('0','3'), c('3','15'), c('0','15'), c('0','9'))))
}

#coag inhibitor
plot_protein_comparisons('SERPINC1',list(c('0','3'), c('3','9'), c('3','12'), c('0','9')))
# anova -------------------------------------------------------------------

#how much of the varation can we explain across age, for a given protein
#is there a dif in log2 intesnities between the age group w/r to log2 intesntiteis
# Load necessary libraries
library(ggplot2)
library(reshape2)

perform_anova <- function(data) {
  # Create an empty list to store p-values for each protein
  anova_results <- data.frame(PG.Genes = character(), p_value = numeric(), stringsAsFactors = FALSE)
  
  # Loop through each unique protein
  for (protein in unique(data$PG.Genes)) {
    
   anova_tmp <- data %>%
      filter(PG.Genes == protein) %>%
      select(sample_no, group, log2) %>% 
      distinct()
    
   anova_tmp$group <- as.factor(anova_tmp$group)
    
    anova_res <- aov(log2 ~ group, data = anova_tmp)
    anova_summary <- summary(anova_res)
    p_value <- anova_summary[[1]]["Pr(>F)"][1,1]
    
   anova_results <- rbind(anova_results, data.frame(PG.Genes = protein, p_value = p_value))
  }
 anova_results$adjusted_p_value <- p.adjust(anova_results$p_value, method = "fdr")
  
  return(anova_results)
}
test <- perform_anova(final_df_groups)


signif_genes <- test %>%
  filter(adjusted_p_value <0.05) %>%
  arrange(adjusted_p_value) %>%
  #slice(1:100) %>%
  pull(PG.Genes)

normalized_data <- final_df_groups %>%
  filter(PG.Genes %in% signif_genes) %>%
  group_by(PG.Genes) %>%
  mutate(log2_norm = log2 - median(log2, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(group = factor(group, seq(0,15,by=3)))



heatmap_data <- normalized_data %>%
  mutate(group_id = paste0(group,'_',sample_no))%>%
  select(PG.Genes, group_id, log2_norm) %>%
  distinct() %>%
  spread(key = group_id, value = log2_norm)

colnames(heatmap_data)[-1] <- gsub("_.*", "", colnames(heatmap_data)[-1])
col_order <- order(as.numeric(gsub("[^0-9]", "", colnames(heatmap_data)[-1])))
heatmap_data <- heatmap_data[, c(1, col_order + 1)]

library(pheatmap)

# Convert heatmap_data to matrix
heatmap_matrix <- as.matrix(heatmap_data[,-1])  # Remove the PG.Genes column
rownames(heatmap_matrix) <- heatmap_data$PG.Genes

# Plot the heatmap using pheatmap
pheatmap(heatmap_matrix, 
         color = colorRampPalette(c("blue", "white", "red"))(50),  # Color scale
         cluster_rows = TRUE,  # Cluster proteins
         cluster_cols = F,  # Cluster samples
         fontsize_row = 8, 
         fontsize_col = 8,
         border_color = NA)  # Font size for colu



test_no0 <- perform_anova(final_df_groups %>%
                            filter(group != 0))
testies <- setdiff(signif_genes_no, signif_genes)

signif_genes_no <- test_no0 %>%
  filter(adjusted_p_value <0.05) %>%
  arrange(adjusted_p_value) %>%
  arrange(PG.Genes) %>%
  pull(PG.Genes)

normalized_data <- final_df_groups %>%
  filter(PG.Genes %in% signif_genes, group !=0) %>%
  group_by(PG.Genes) %>%
  mutate(log2_norm = log2 - median(log2, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(group = factor(group, seq(0,15,by=3)))



heatmap_data <- normalized_data %>%
  mutate(group_id = paste0(group,'_',sample_no))%>%
  select(PG.Genes, group_id, log2_norm) %>%
  distinct() %>%
  spread(key = group_id, value = log2_norm)

colnames(heatmap_data)[-1] <- gsub("_.*", "", colnames(heatmap_data)[-1])
col_order <- order(as.numeric(gsub("[^0-9]", "", colnames(heatmap_data)[-1])))
heatmap_data <- heatmap_data[, c(1, col_order + 1)]

library(pheatmap)

# Convert heatmap_data to matrix
heatmap_matrix <- as.matrix(heatmap_data[,-1])  # Remove the PG.Genes column
rownames(heatmap_matrix) <- heatmap_data$PG.Genes

# Plot the heatmap using pheatmap
pheatmap(heatmap_matrix, 
         color = colorRampPalette(c("blue", "white", "red"))(50),  # Color scale
         cluster_rows = TRUE,  # Cluster proteins
         cluster_cols = F,  # Cluster samples
         fontsize_row = 8, 
         fontsize_col = 8,
         border_color = NA)  # Font size for colu



jf %>%
  filter(group ==0) %>%
  arrange(age) %>%
  select(age,sample_no) %>%
  distinct()

plot_protein_comparisons('AMH', )



for (prot in signif_genes) {
  print(plot_protein_comparisons(prot,list(c('0','3'), c('3','9'), c('9','15'))))
}
# log log de --------------------------------------------------------------

log_log_de <- function( midpoint) {
  
  de1 <- differential_expression_analysis(jf, 0,midpoint,
                                          imp_comp = 0.667,
                                          completeness = 0.667, 
                                          test_lookup) #
  de2 <- differential_expression_analysis(jf, midpoint,15,
                                          imp_comp = 0.667,
                                          completeness = 0.667, 
                                          test_lookup) #
  combined_df   <- merge(de1, de2, by='PG.Genes', all=T, suffix = c("meth1", "meth2"))
  
  
  
  combined_df <- combined_df %>%
    mutate(meth1 = sign(fcmeth1) * -log(adjpmeth1,10), 
           meth2 = sign(fcmeth2) * -log(adjpmeth2,10))
  
  #plotting
  plot <- combined_df %>%
    mutate(sig = case_when(
      meth1 > 2 & meth2 < -2 ~ '2', 
      meth2 > 2 & meth1 < -2 ~ '2',
      meth1 > 1.3 & meth2 < -1.3 ~ '1', 
      meth2 > 1.3 & meth1 < -1.3 ~ '1',
      T ~ '0'
     
    )) %>%
    ggplot(aes(x=meth1, y=meth2)) + 
    geom_point(aes(colour=sig), show.legend=F) + 
    scale_colour_manual(values =c('grey', 'red','red')) +
    scale_x_continuous(breaks=seq(-18, 17, by= 2)) + 
    scale_y_continuous(breaks = seq(-18,17, by=2)) + 
    geom_vline(xintercept = c(-1.3, +1.3), linetype = 2) +
    geom_vline(xintercept = 0)+ 
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = c(-1.3, +1.3), linetype = 2)+
    geom_abline(slope =1) +
    labs(x=paste0('0 vs. ', midpoint), y= paste0(midpoint,' vs. 15')) +
    geom_text_repel(aes(
      label = ifelse(sig=='2', PG.Genes, "")), 
      size = 4, color = "red", 
      max.overlaps = Inf, 
      min.segment.length = 0, 
      box.padding = 1,  # Increased padding around text box
      point.padding = 0.5  # Padding between points and text labels
    ) +
    
    theme 
  print(plot)
}
midpoints = seq(3,12, by=3)
for (mp in midpoints) {
  log_log_de(mp)
}
plot_protein_comparisons('AMH', list(), c('Male'))

jf, 
group1, group2, 
imp_comp = 0.667,
completeness = 0.667, 
test_looku


