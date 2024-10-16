#DIA windows 
library('ggplot2')
library('tidyverse')
library('stringr')
library(RColorBrewer)

theme <- theme_classic() + theme(
  text = element_text((family = "Arial"), size = 12, colour ="black"),  # Set all text to Arial, size 12
  plot.title = element_text((family = "Arial"), size = 14),  # Adjust title size (optional)
  #strip.text.x = element_text((family = "Arial"), size = 12),  # Adjust strip text size (optional)
  axis.text = element_text((family = "Arial"), size = 12, colour ="black"),
  legend.text = element_text((family = "Arial"), size = 12, colour ="black")
  # ... adjust other text elements as needed
)



# script ------------------------------------------------------------------


#read in file
file = 'initial screen/report.tsv'
agc_report <- generate_report(file, is_diann =T, map_concs =F)

filters = list(
  list(pep_threshold = 0.05, pg_qval = 0.05),
  list(pep_threshold = 0.01, pg_qval = 0.01)
)

agc_report <- agc_report %>%
  filter(PEP.IsProteotypic == 'True')

agc_summary <- generate_summary_report(agc_report, filters)


#get prec mz values
fr = 'initial screen/report-lib.tsv'
mz_report <- read.csv(fr, sep = '\t')
mz_report <- mz_report %>%
  select(PrecursorMz, transition_group_id) %>%
  distinct() %>%
  rename(Precursor.Id = transition_group_id)


  
# detection - inital screen ---------------------------------------------------------------
plot_protein_precursor_errorbars <- function(summary_d, mss = "Astral", coeff = 9.3, fill_colors, fill_labels, x_labels, 
                                             show_protein_error = TRUE, show_precursor_error = TRUE) {
  
  desired_order <- c("Astral_2Th_20_5", "Astral_2Th_20_10", 
                     "Astral_2Th_40_5", "Astral_2Th_40_10", 
                     "Astral_2Th_100_5", "Astral_2Th_100_10")
  # Generate summary stats for protein data
  summary_stats <- summary_d %>%
    mutate(method = factor(method, desired_order)) %>%
    pivot_longer(cols = c(IDstot, IDsfilt_level_1, IDsfilt_level_2, precs),
                 names_to = 'Metric',
                 values_to = 'Values') %>%
    group_by(method, Metric) %>%
    summarise(mean = mean(Values),
              sd = sd(Values),
              se = sd(Values) / sqrt(n()),  # Calculate standard error
              .groups = 'drop') %>%
    separate(method, c('ms', 'wind', 'agc', 'it'), remove = F) %>%
    filter(ms == "Astral")
    
  
  # Generate data for plotting precursors
  line <- summary_d %>%
    select(method, avgprec, precs) %>%
    mutate(method = factor(method, desired_order)) %>%
    group_by(method) %>%
    summarise(avgprec = mean(precs),
              se_prec = sd(precs) / sqrt(n()),  # Calculate standard error
              .groups = 'drop') %>%
    separate(method, c('ms', 'wind', 'agc', 'it'), sep = "_", remove = F) %>%
    filter(ms == "Astral")
  
  # Main plot
  t <- summary_d %>%
    select(method, IDstot, IDsfilt_level_1, IDsfilt_level_2, replicate) %>%
    pivot_longer(cols = c(IDstot, IDsfilt_level_1, IDsfilt_level_2),
                 names_to = 'Metric',
                 values_to = 'Values') %>%
    mutate(Metric = factor(Metric, c('IDstot', 'IDsfilt_level_1', 'IDsfilt_level_2', 'precs'))) %>%
    mutate(method = factor(method, desired_order)) %>%
    separate(method, c('ms', 'wind', 'agc', 'it'), sep = "_", remove = F) %>%
    mutate(mm = paste0(method, Metric)) %>%
    group_by(mm) %>% 
    mutate(mean = mean(Values)) %>% 
    mutate(agc = factor(agc, c(20,40,100))) %>%
    filter(ms == "Astral") %>%
    
    ggplot(aes(x = method, y = mean, fill = agc)) + 
    geom_bar(aes(alpha = Metric), colour = 'black', size = 0.5, stat = 'identity', position = 'dodge') + 
    
    geom_point(aes(x = method, y = Values, fill = agc, group = Metric), shape = 21,
               size = 1.5, position = position_jitterdodge(jitter.width = 0.06, dodge.width = 0.9), show.legend = F) +
    
    scale_fill_manual(values = fill_colors,
                      labels = fill_labels,
                      name = "AGC Level") + 
    
    scale_alpha_manual(values = c(0.1, 0.3, 1),
                       labels = c('Proteotypic', "FDR < 5%", "FDR < 1%"),
                       name = "Protein ID's") +
    
    scale_x_discrete(labels = x_labels, name = 'Maximun Injection Time (ms)') +
    
    coord_cartesian(expand = FALSE, xlim = c(0.3, 6.7), ylim = c(0, 1500)) +
    
    theme +
    
    geom_line(data = line, aes(x = method, y = avgprec / coeff, group = agc), linewidth = 0.5) + 
    geom_point(data = line, aes(x = method, y = avgprec / coeff), size = 2.5, colour = 'red', show.legend = FALSE) +
    
    scale_y_continuous(name = "Number of Proteins",
                       sec.axis = sec_axis(~.*coeff, name = 'Number of Precursors', breaks = seq(0, 12000, by = 1500)),
                       breaks = seq(0, 2000, by = 100)) + 
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
plot_protein_precursor_errorbars(agc_summary, mss = "Astral", coeff = 5.5, 
                                 fill_colors = c('#00BCD4', "#CCFF00", "#ff00ff"), 
                                 fill_labels = c("20%", "40%", "100%"), 
                                 x_labels = c('5', '10', '5', '10', '5', '10'),
                                 show_protein_error = F, 
                                 show_precursor_error = T)

prot = 'APOA1'
test1 <- agc_report %>%
  select(R.FileName, PG.Genes, log2) %>% 
  distinct() %>%
  filter(PG.Genes == prot)

test2 <- opt_report %>%
  select(R.FileName, PG.Genes, log2) %>% 
  distinct() %>%
  filter(PG.Genes == prot)

test3 <- final_report %>%
  select(R.FileName, PG.Genes, log2) %>% 
  distinct() %>%
  filter(PG.Genes == prot)
test<-rbind(test1,test2) 
test <- rbind(test, test3)
test %>% 
  mutate(num = row_number()) %>%
  ggplot(aes(x=num, y=log2)) + 
  geom_point() +ylim(30,35) +
  geom_vline(xintercept = c(12.5,22.5), linetype='dashed') +
  scale_x_continuous(breaks = seq(1, 100, by=50)) +
  labs(y=expression('Protein Group Intensity (log'[2]*')'), x='Sample Type')+
  theme
t <- final_summary %>%
  separate(method, c('instrument', 'window', 'agc', 'it')) %>%
  select(window, agc, it, replicate, precs, IDstot, IDsfilt_level_1, IDsfilt_level_2)
# rank baundance - inital screen ---------------------------------------------------------


# Filter IDs for the 20_10 method
ids_20_10 <- agc_report %>%
  filter(R.FileName %in% c('P4899_2Th_20_10_R1','P4899_2Th_20_10_R2'),
         PEP.IsProteotypic == "True"
         #PG.PEP..Run.Wise. < 0.01,
         #PG.QValue..Run.Wise.<0.01,
  ) %>%
  select(PG.Genes) %>%
  distinct()

# Filter IDs for the 100_10 method
ids_100_10 <- agc_report %>%
  filter(R.FileName %in% c('P4899_2Th_100_5_R1','P4899_2Th_100_5_R1'),
         PEP.IsProteotypic == "True"
         #PG.PEP..Run.Wise. < 0.01,
         #PG.QValue..Run.Wise.<0.01
  ) %>%
  select(PG.Genes) %>%
  distinct()

# Get the different proteins between 20_10 and 100_10
diff_proteins <- setdiff(ids_20_10, ids_100_10)
diff_proteins <- diff_proteins$PG.Genes

# Generate the ranked abundance plot
agc_report_rank <- agc_report %>%
  filter(R.FileName %in% c('P4899_2Th_20_10_R1','P4899_2Th_20_10_R2'),
         PEP.IsProteotypic=='True'
         #PG.QValue..Run.Wise.<0.01
  ) %>%
  select(PG.Genes, PG.Quantity, log2) %>%
  distinct() %>%
  group_by(PG.Genes) %>%
  summarize(average_log2 = mean(log2)) %>%
  arrange(desc(average_log2)) %>%
  mutate(rank = row_number())

t <- agc_report_rank %>%
  mutate(
    test = case_when(PG.Genes %in% diff_proteins ~ T))

# Highlight QC proteins
t <- agc_report_rank %>%
  filter(average_log2 !=-Inf) %>%
  highlight_qc_proteins(diff_proteins) %>% #set what porteins you want highlifhted
  arrange(is_qc_protein) %>%
  mutate(is_qc_protein = factor(is_qc_protein, c(T,F))) %>%
  ggplot(aes(x = rank, y = average_log2)) +
  geom_point(aes(colour = is_qc_protein), size = 1.5, show.legend = F) +  
  labs(x = "Protein Rank", y = expression("Protein Intensity (log"[2]*")"), 
       color = "GPE") +
  scale_colour_manual(values = c('#00BCD4','black'), labels=c('Gas Phase Enriched', '')) +
  #scale_x_continuous(breaks = seq(0,1000, by=200), limits = c(0,1000))+
  scale_y_continuous(breaks = seq(0,50, by =5), limits=c(5,35)) +
  #coord_cartesian(expand = F, xlim=c(-20, 1630), ylim=c(7,35.2) ) + 
  theme

print(t)



# spectral counts - inital screen -----------------------------------------
counts <- read.delim('initial screen/counts.csv', sep = ',')
# Step 1: Summarize the data to calculate mean and standard error
summary_counts <- counts %>%
  group_by(agc, it) %>%
  summarise(
    mean_counts = mean(scan.counts),
    se_counts = sd(scan.counts) / sqrt(n())  # Standard Error
  )

# Step 2: Plot the averages with error bars
summary_counts %>%
  ggplot(aes(x = agc, y = mean_counts/1e5, colour = as.factor(it))) + 
  geom_point(size = 3) +  # Points for mean values
  geom_line(aes(group = it), size = 1) +  # Line connecting the mean values
  geom_errorbar(aes(ymin = (mean_counts - se_counts)/1e5, ymax = (mean_counts + se_counts)/1e5), 
                width = 1.5, colour='black') +  # Error bars
  theme + 
  scale_colour_manual(values = c('#B2DF6E','#FF4500'),
                      labels = c('5ms','10ms'))+
  scale_y_continuous(breaks = seq(0,3, by=.25))+
  coord_cartesian(expand=F, xlim=c(17,103),ylim=c(1,2.5)) +
  labs(x = "AGC Level (%)", y = expression("Total MS2 Scans (10"^5*")"), colour = "Max IT") 
   


test <- agc_summary %>% 
  select(method, precs, IDstot, IDsfilt_level_1, IDsfilt_level_2)



# ion traces  -------------------------------------------------------------

files = c(
  'initial screen/P4899_2Th_20_5_R1.csv',
  'initial screen/P4899_2Th_20_10_R1.csv',
  'initial screen/P4899_2Th_40_10_R2.csv',
  'initial screen/P4899_2Th_100_10_R1.csv'
  ) 

for (file in files) {
  test <- read.delim(file, sep = ',', header =F)
  cycle <- test %>% 
    mutate(cycle = case_when(V3==980 ~ T)) %>%
    select(V1, cycle) %>%
    filter(cycle ==T)
  mz_breaks = cycle$V1
    
    
  plot <- test %>%
    ggplot(aes(x=V1,y=V2))  +
    geom_line(size=.6) + 
    xlim(20, 20.05) +
    ylim(0,10) +
    geom_vline(xintercept = mz_breaks, colour = 'red',  size =1) +
    theme+ 
    labs(title=file)
    
  
  print(plot)
  
  
  # Create the lollipop plot
   plot <- test %>%
      ggplot(aes(x=V1,y=V2))  +
      geom_line(size=.6) + 
     
      scale_y_continuous(breaks=seq(0,15, by=5), limits=c(0,11))+
      geom_vline(xintercept = mz_breaks, colour = '#FF4500',  size =1) +
      theme+ 
      scale_x_continuous(breaks=seq(0,25,by=0.01), limits=c(x_start, x_end+.002))+
      labs( x= NULL, y= 'IT (ms)') + 
      coord_cartesian(expand=F)
  # Print the plot
  print(plot)
  
}



set_base_plot_theme <- function() {
  
  par(
    family = "Arial",        # Set font to Arial
    ## Remove box.
    bty = "n",
    ## Remove default x and y axis.
    xaxt = "n", yaxt = "n",
    bty = "l",
    ## Increase width of box lines.
    lwd = .1
    
  )
}

# Example usage
set_base_plot_theme()  # Apply the theme

# Now create your plots as usual
# Function to plot the cycle data based on the theme
plot_cycle_data_base <- function(files, start_rt) {
  set_base_plot_theme()  # Apply the theme before plotting
  
  for (file in files) {
    # Read the data file
    test <- read.delim(file, sep = ',', header = FALSE)
    
    # Identify the cycles
    cycle <- test %>%
      mutate(cycle = case_when(V3 == 980 ~ TRUE)) %>%
      select(V1, cycle) %>%
      filter(cycle == TRUE)
    mz_breaks <- cycle$V1
    
    # Find the closest mz break to the specified start RT
    closest_start <- mz_breaks[which.min(abs(mz_breaks - start_rt))]
    
    # Ensure there are at least two cycles to plot
    if (length(mz_breaks) < 2) {
      next
    }
    
    # Find the position of the closest start and the next mz break
    start_index <- which(mz_breaks == closest_start)
    if (start_index >= length(mz_breaks)) {
      next
    }
    x_start <- closest_start - 0.005
    x_end <- mz_breaks[start_index + 2] + 0.005

    
    plot <- test %>%
      ggplot(aes(x=V1,y=V2))  +
      geom_line(size=.6) + 
      #geom_point(size = 1, colour='black') +
      scale_y_continuous(breaks=seq(0,15, by=5), limits=c(0,11))+
      geom_vline(xintercept = mz_breaks, colour = '#FF4500',  size =1) +
      theme+ 
      scale_x_continuous(breaks=seq(0,25,by=0.01), limits=c(x_start, x_end+.002))+
      labs( x= NULL, y= 'IT (ms)') + 
      coord_cartesian(expand=F)
    
    
    print(plot)
  }
}

# Call the function to plot the data
files = c(
  'initial screen/P4899_2Th_20_5_R1.csv',
  'initial screen/P4899_2Th_20_10_R1.csv',
  'initial screen/P4899_2Th_40_10_R2.csv',
  'initial screen/P4899_2Th_100_10_R1.csv'
)
files = c('opt/P4899_1ex_20_10_R1.csv',
          'opt/P4899_2ex_20_10_R2.csv',
          'opt/P4899_3Th_20_10_R2.csv')

plot_cycle_data_base( files, start_rt=8)

#16.5,
# mz vs rt  ---------------------------------------------------------------


#mz_report and agc_report

it_report <- agc_report %>%
  filter(R.FileName == 'P4899_2Th_20_10_R1') %>%
  select(Precursor.Id, RT, log2) %>%
  distinct()

it_report <- merge(it_report, mz_report, by='Precursor.Id')

it_report %>% 
  filter(log2>22.80567) %>% #filtering for log2
  select(PrecursorMz) %>% 
  ggplot(aes(x=PrecursorMz)) + 
  geom_histogram( fill='#8A8D91', binwidth = 10, colour='black', size=.5, show.legend = T, boundary = 380) + 
  #scale_x_continuous(breaks = seq(380, 980, by=100), limits=c(380,980)) +
  scale_y_continuous(breaks= seq(0,300, by=50), limits = c(0,150)) +
  theme + 
  labs(x=expression(italic('m/z')), y='Count') +
  coord_cartesian(expand=F, ylim = c(0,150), xlim=c(380,988)) 

agc_report %>%
  filter(R.FileName == 'P4899_2Th_20_10_R1') %>%
  select(EG.PrecursorId, log2) %>%
  distinct() %>%
  ggplot(aes(x=log2)) + 
  xlim(5,35) + 
  geom_histogram(fill='#8A8D91', binwidth =.5) +
  geom_vline(aes(xintercept = 22.8), 
             color = "red", linetype = "dashed", size = 1) +
  labs(x=expression('Protein Group Intensity (log'[2]*')'), y = 'Count') +
  coord_cartesian(expand = F) + 
  theme

mean_value <- median(agc_report %>%
                     filter(R.FileName == 'P4899_2Th_20_10_R1') %>%
                     select(EG.PrecursorId, log2) %>%
                     distinct() %>%
                     pull(log2) %>%
                     na.omit() %>%  # Remove NA values
                     .[is.finite(.)],  # Keep only finite values (remove -Inf, Inf)
                   na.rm = TRUE)  # Ensure NA handling

# Print the result
print(mean_value)


p <- it_report %>%
  ggplot(aes(y=RT, x=PrecursorMz)) +
  geom_density_2d_filled(contour_var = "count") +
  scale_fill_viridis_d() +  
  theme + 
  scale_y_continuous(breaks = seq(0, 25, by=5), limits=c(2, 25)) +
  scale_x_continuous(breaks = seq(380, 980, by=100), limits=c(380, 980)) +
  coord_cartesian(expand = F) + 
  labs(x=expression(italic('m/z')), y='Retention Time (mins)', fill = 'Precursor\n Density')
print(p)
ggsave(filename = "plot.svg", plot = p, device = "svg", dpi = 1400)

# clustering of prec mz  -----------------------------------------------------

library(Ckmeans.1d.dp)

data <- it_report 
set.seed(123)  # For reproducibility
k.max <- 8
wss <- sapply(1:k.max, function(k) {
  kmeans(data[, c("RT", "PrecursorMz")], centers = k, nstart = 10)$tot.withinss
})

elbow_data <- data.frame(Clusters = 1:k.max, WSS = wss)

elbow_data %>%
  ggplot(aes(x = Clusters, y = WSS/1e8)) +
  geom_point(size = 3) +
  geom_line() +
  theme + 
  labs(x = "Number of Clusters",
       y = expression("Total WSS (10"^8*")")) +
  scale_x_continuous(breaks=seq(0,10,by=1)) +
  scale_y_continuous(breaks = seq(0,2.5, by=.5)) +
  coord_cartesian(expand =F, ylim = c(0,2.5), xlim = c(.7,8.3))


calculate_rt_midpoints <- function(data, optimal_clusters) {
  # Perform k-means clustering
  kmeans_result <- kmeans(data[, c("RT", "PrecursorMz")], centers = optimal_clusters, nstart = 10)
  data$cluster <- factor(kmeans_result$cluster)
  
  # Calculate midpoints for each cluster based on RT
  rt_midpoints <- data %>%
    group_by(cluster) %>%
    summarise(rt_midpoint = median(RT))
  
  return(rt_midpoints$rt_midpoint)
}

test <- calculate_rt_midpoints(it_report, 1)

# Generate the plot
rt_midpoint <- calculate_rt_midpoints(it_report, 1) #13.8261

it_report %>%
  mutate(V1 = case_when(
    RT> 13.8261 ~ T, 
    T ~ F #so everything else is flase
  )) %>% 
  ggplot(aes(y=RT, x=PrecursorMz)) +
  geom_point(aes(colour = V1), shape=19, alpha=1, size=1, show.legend=T) +
  scale_colour_manual(values = c('#8A8D91', '#00BCD4'))+
  theme + 
  scale_y_continuous(breaks = seq(0, 25, by=5), limits=c(2, 25)) +
  scale_x_continuous(breaks = seq(380, 980, by=100), limits=c(375, 980)) +
  coord_cartesian(expand = F, xlim=c(), ylim=c())+
  labs(x=expression(italic('m/z')), y='Retention Time (mins)') + 
  geom_hline(yintercept = 13.8261, colour ='black')


it_report %>%
  mutate(V1 = case_when(
    RT< 13.8261 ~ T, 
    T ~ F #so everything else is flase
  )) %>%
  filter(RT >13.8261) %>%
  ggplot(aes(x=PrecursorMz)) + 
  geom_histogram( fill='#00BCD4', binwidth = 10, colour='black', size=.5, show.legend = T, boundary = 380) + 
  scale_x_continuous(breaks = seq(380, 980, by=100), limits=c(380,980)) +
  scale_y_continuous(breaks= seq(0,300, by=50)) +
  theme + 
  labs(x=expression(italic('m/z')), y='Count') +
  coord_cartesian(expand=F, ylim = c(0,150), xlim=c(380,988))


# dyanmic window gneeration  ----------------------------------------------


test <- calculate_rt_midpoints(it_report, 1) # 13.821
test <- calculate_rt_midpoints(it_report, 2)  #12.4167, 15.96
# Generate dynamic windows based on precursor density
window_borders = list(c(), #empy, do the whole thing 
                      c(13.821), # two clustesr
                      c(12.41, 15.96)) # three clusters
#grab optimal window broders from asttral
astral1Th <- astral1Th %>%
  separate(col = 1, into = c("part1", "part2"), sep = "-")
vec <- astral1Th$part1
vec <- c(vec, 980.695655)
optimised_window_borders <- as.numeric(vec)


#### k=2, 2 groups
tmp_2a <- it_report %>%
  filter(RT <= 13.821) 
tmp_2b <- it_report  %>%
  filter(RT > 13.821) 

#### k=3, 3 groups
tmp_3a <- it_report  %>%
  filter(RT <= 12.41) 
tmp_3b <- it_report  %>%
  filter(12.41 < RT &RT <= 15.96) 
tmp_3c <- it_report  %>%
  filter(RT > 15.96) 

# Define the datasets and their names
ds <- list(it_report, tmp_2a, tmp_2b, tmp_3a, tmp_3b, tmp_3c)
ds_names <- c('it_report', "tmp_2a", "tmp_2b", "tmp_3a", "tmp_3b", "tmp_3c")
# Number of windows
var_methods <- 201 #-1 since add 2 windows (top and bottom)

# Iterate through each dataset
for (i in seq_along(ds)) {
  vec <- ds[[i]]$PrecursorMz
  mz_report_vec <- sort(vec)
  
  vec_astral <- c(380.422805)  # Starting point
  
  bins <- round(length(mz_report_vec) / var_methods)  # Calculate bins
  
  for (count in seq_along(mz_report_vec)) {
    if (count %% bins == 0) {
      num <- mz_report_vec[count]
      
      # Find the nearest value in 'optimised_window_borders'
      num_opt3 <- which.min(abs(optimised_window_borders - num))
      vec_astral <- c(vec_astral, optimised_window_borders[num_opt3])
    }
  }
  
  vec_astral <- c(vec_astral, 980.695655)  # Ending point
  vec_astral <- unique(vec_astral)  # Remove duplicates
  
  print(length(vec_astral))
  
  # Create m/z ranges in the astral format
  mz_list <- sapply(1:(length(vec_astral) - 1), function(j) {
    paste0(vec_astral[j], '-', vec_astral[j + 1])
  })
  
  # Write the results to a CSV file
  file_name <- paste0("vec_astral_", ds_names[i], ".csv")
  mz_df <- data.frame('m/z range' = mz_list)
  print(plot_mz_windows(mz_df))
  write.csv(mz_df, file = file_name, row.names = FALSE)
}



# agc opt script ----------------------------------------------------------

file= 'opt/report_new.tsv'
opt_report <- generate_report(file, is_diann = T, map_concs=F)
opt_report <- opt_report %>%
  filter(PEP.IsProteotypic == 'True')
opt_summary <- generate_summary_report(opt_report, filters)



# agc opt chart -----------------------------------------------------------

plot_protein_precursor_errorbars <- function(summary_d, mss = "Astral", coeff = 9.3, fill_colors, fill_labels, x_labels, 
                                             show_protein_error = TRUE, show_precursor_error = TRUE) {
  
  desired_order <- c("Astral_2Th_20_10", "Astral_3Th_20_10","Astral_1ex_20_10",
                     "Astral_2ex_20_10", "Astral_3ex_20_10")
  summary_d <- opt_summary
  # Generate summary stats for protein data
  summary_stats <- summary_d %>%
    mutate(method = factor(method, desired_order)) %>%
    pivot_longer(cols = c(IDstot, IDsfilt_level_1, IDsfilt_level_2, precs),
                 names_to = 'Metric',
                 values_to = 'Values') %>%
    group_by(method, Metric) %>%
    summarise(mean = mean(Values),
              sd = sd(Values),
              se = sd(Values) / sqrt(n()),  # Calculate standard error
              .groups = 'drop') %>%
    separate(method, c('ms', 'wind', 'agc', 'it'), remove = F, sep='_') %>%
    filter(ms == "Astral") %>%
    mutate(Metric = factor(Metric, c('IDstot', 'IDsfilt_level_1', 'IDsfilt_level_2', 'precs')),
           exp_no = case_when(
             wind == '3ex' ~ 'three',
             wind == '2ex' ~ 'two', 
             wind %in% c('2Th', '1ex', '3Th') ~ 'one'
           ))
  
  # Generate data for plotting precursors
  line <- summary_d %>%
    select(method, avgprec, precs) %>%
    mutate(method = factor(method, desired_order)) %>%
    group_by(method) %>%
    summarise(avgprec = mean(precs),
              se_prec = sd(precs) / sqrt(n()),  # Calculate standard error
              .groups = 'drop') %>%
    separate(method, c('ms', 'wind', 'agc', 'it'), sep = "_", remove = F) %>%
    mutate(exp_no = case_when(
             wind == '3ex' ~ 'three',
             wind == '2ex' ~ 'two', 
             wind %in% c('2Th', '1ex', '3Th') ~ 'one'
           )) %>%
    filter(ms == "Astral")
  
  # Main plot
  t <- summary_d %>%
    select(method, IDstot, IDsfilt_level_1, IDsfilt_level_2, replicate) %>%
    pivot_longer(cols = c(IDstot, IDsfilt_level_1, IDsfilt_level_2),
                 names_to = 'Metric',
                 values_to = 'Values') %>%
    mutate(Metric = factor(Metric, c('IDstot', 'IDsfilt_level_1', 'IDsfilt_level_2', 'precs'))) %>%
    mutate(method = factor(method, desired_order)) %>%
    separate(method, c('ms', 'wind', 'agc', 'it'), sep = "_", remove = F) %>%
    mutate(mm = paste0(method, Metric)) %>%
    mutate(exp_no = case_when(
      wind == '3ex' ~ 'three',
      wind == '2ex' ~ 'two', 
      wind %in% c('2Th', '1ex', '3Th') ~ 'one'
    )) %>%
    group_by(mm) %>% 
    mutate(mean = mean(Values)) %>% 
    mutate(exp_no = factor(exp_no, c('one', 'two', 'three'))) %>%
    filter(ms == "Astral") %>%
    
    ggplot(aes(x = method, y = mean, fill = exp_no)) + 
    geom_bar(aes(alpha = Metric), colour = 'black', size = 0.5, stat = 'identity', position = 'dodge') + 
    
    geom_point(aes(x = method, y = Values, fill = exp_no, group = Metric), shape = 21,
               size = 1.5, position = position_jitterdodge(jitter.width = 0.06, dodge.width = 0.9), show.legend = F) +
    
    scale_fill_manual(values = fill_colors,
                      labels = fill_labels,
                      name = "DIA Experiments") + 
    
    scale_alpha_manual(values = c(0.1, 0.3, 0.9),
                       labels = c('Proteotypic', "5%", "1%"),
                       name = "Filtered") +
    
    scale_x_discrete(labels = x_labels, name = NULL) +
    
    coord_cartesian(expand = FALSE, xlim = c(0.3, 5.7), ylim = c(0, 1400)) +
    
    theme +
    
    geom_line(data = line, aes(x = method, y = avgprec / coeff, group = exp_no), linewidth = 0.5) + 
    geom_point(data = line, aes(x = method, y = avgprec / coeff), size = 2.5, colour = 'red', show.legend = FALSE) +
    
    scale_y_continuous(name = "Number of Proteins",
                       sec.axis = sec_axis(~.*coeff, name = 'Number of Precursors', breaks = seq(0, 12000, by = 1500)),
                       breaks = seq(0, 2000, by = 100)) + 
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
plot_protein_precursor_errorbars(opt_summary, mss = "Astral", coeff = 5.6, 
                                 fill_colors = c('#00BCD4', "#007FFF", "#FF8C00"), 
                                 fill_labels = c("1", "2", "3"), 
                                 x_labels = c('2Th\n(300)', '3Th\n(200)', 
                                              'Variable\n(200)', 'Variable\n(200)', 'Variable\n(200)'),
                                 show_protein_error = F, 
                                 show_precursor_error = T)



# hidden scasn  -----------------------------------------------------------
file= 'op2/report_hidden.tsv'
hidden_report <- generate_report(file, is_diann = T, map_concs=F)
hidden_report <- hidden_report %>%
  filter(PEP.IsProteotypic == 'True')
hidden_summary <- generate_summary_report(hidden_report, filters)



# hidden scasn detec ------------------------------------------------------

plot_protein_precursor_errorbars <- function(summary_d, mss = "Astral", coeff = 9.3, fill_colors, fill_labels, x_labels, 
                                             show_protein_error = TRUE, show_precursor_error = TRUE) {
  
  desired_order <- c("Astral_2Th-10_20_10", "Astral_2Th-6_20_10","Astral_2Th-4_20_10")
  summary_d <- hidden_summary
  # Generate summary stats for protein data
  summary_stats <- summary_d %>%
    mutate(method = factor(method, desired_order)) %>%
    pivot_longer(cols = c(IDstot, IDsfilt_level_1, IDsfilt_level_2, precs),
                 names_to = 'Metric',
                 values_to = 'Values') %>%
    group_by(method, Metric) %>%
    summarise(mean = mean(Values),
              sd = sd(Values),
              se = sd(Values) / sqrt(n()),  # Calculate standard error
              .groups = 'drop') %>%
    separate(method, c('ms', 'wind', 'agc', 'it'), remove = F, sep='_') %>%
    filter(ms == "Astral") %>%
    mutate(Metric = factor(Metric, c('IDstot', 'IDsfilt_level_1', 'IDsfilt_level_2', 'precs')),
           exp_no = case_when(
             wind == '2Th-4' ~ '4',
             wind == '2Th-6' ~ '6', 
             wind =='2Th-10' ~ '10'
           ))
  
  # Generate data for plotting precursors
  line <- summary_d %>%
    select(method, avgprec, precs) %>%
    mutate(method = factor(method, desired_order)) %>%
    group_by(method) %>%
    summarise(avgprec = mean(precs),
              se_prec = sd(precs) / sqrt(n()),  # Calculate standard error
              .groups = 'drop') %>%
    separate(method, c('ms', 'wind', 'agc', 'it'), sep = "_", remove = F) %>%
    mutate(exp_no = case_when(
             wind == '2Th-4' ~ '4',
             wind == '2Th-6' ~ '6', 
             wind =='2Th-10' ~ '10'
           )) %>%
    filter(ms == "Astral")
  
  # Main plot
  t <- summary_d %>%
    select(method, IDstot, IDsfilt_level_1, IDsfilt_level_2, replicate) %>%
    pivot_longer(cols = c(IDstot, IDsfilt_level_1, IDsfilt_level_2),
                 names_to = 'Metric',
                 values_to = 'Values') %>%
    mutate(Metric = factor(Metric, c('IDstot', 'IDsfilt_level_1', 'IDsfilt_level_2', 'precs'))) %>%
    mutate(method = factor(method, desired_order)) %>%
    separate(method, c('ms', 'wind', 'agc', 'it'), sep = "_", remove = F) %>%
    mutate(mm = paste0(method, Metric)) %>%
    mutate(exp_no = case_when(
      wind == '2Th-4' ~ '4',
      wind == '2Th-6' ~ '6', 
      wind =='2Th-10' ~ '10'
    )) %>%
    group_by(mm) %>% 
    mutate(mean = mean(Values)) %>% 
    mutate(exp_no = factor(exp_no, c('one', 'two', 'three'))) %>%
    filter(ms == "Astral") %>%
    
    ggplot(aes(x = method, y = mean)) + 
    geom_bar(aes(alpha = Metric),fill='#00BCD4', colour = 'black', size = 0.5, stat = 'identity', position = 'dodge') + 
    
    geom_point(aes(x = method, y = Values, fill=Metric, group = Metric), shape = 21,
               size = 1.5, position = position_jitterdodge(jitter.width = 0.06, dodge.width = 0.9), show.legend = F) +
    
    scale_fill_manual(values= c('#00BCD4','#00BCD4','#00BCD4')) + 
    
    scale_alpha_manual(values = c(0.1, 0.3, 0.9),
                       labels = c('Proteotypic', "5%", "1%"),
                       name = "Filtered") +
    
    scale_x_discrete(labels = x_labels, name = 'AGC Survery Scan Rate') +
    
    coord_cartesian(expand = FALSE, xlim = c(0.3, 3.7), ylim = c(0, 1000)) +
    
    theme +
    
    geom_line(data = line, aes(x = method, y = avgprec / coeff, group = exp_no), linewidth = 0.5) + 
    geom_point(data = line, aes(x = method, y = avgprec / coeff), size = 2.5, colour = 'red', show.legend = FALSE) +
    
    scale_y_continuous(name = "Number of Proteins",
                       sec.axis = sec_axis(~.*coeff, name = 'Number of Precursors', breaks = seq(0, 12000, by = 1500)),
                       breaks = seq(0, 2000, by = 100)) + 
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
plot_protein_precursor_errorbars(hidden_summary, mss = "Astral", coeff = 5.6, 
                                 fill_colors = c('#00BCD4', "#1E90FF", "#4169E1"), 
                                 fill_labels = c("1", "2", "3"), 
                                 x_labels = c('2 Hz', '3.3 Hz', '5 Hz'),
                                 show_protein_error = F, 
                                 show_precursor_error = T)


# agc sweep 1ex -----------------------------------------------------------

file= 'sweep/report.tsv'
sweep_report <- generate_report(file, is_diann = T, map_concs=F)
sweep_report <- sweep_report %>%
  filter(PEP.IsProteotypic == 'True')
sweep_summary <- generate_summary_report(sweep_report, filters)





plot_protein_precursor_errorbars <- function(summary_d, mss = "Astral", coeff = 9.3, fill_colors, fill_labels, x_labels, 
                                             show_protein_error = TRUE, show_precursor_error = TRUE) {
  
  desired_order <- c("Astral_1ex_5_10", "Astral_1ex_10_10","Astral_1ex_15_10",
                     "Astral_1ex_20_10", "Astral_1ex_40_10")
  summary_d <- sweep_summary
  # Generate summary stats for protein data
  summary_stats <- summary_d %>%
    mutate(method = factor(method, desired_order)) %>%
    pivot_longer(cols = c(IDstot, IDsfilt_level_1, IDsfilt_level_2, precs),
                 names_to = 'Metric',
                 values_to = 'Values') %>%
    group_by(method, Metric) %>%
    summarise(mean = mean(Values),
              sd = sd(Values),
              se = sd(Values) / sqrt(n()),  # Calculate standard error
              .groups = 'drop') %>%
    separate(method, c('ms', 'wind', 'agc', 'it'), remove = F, sep='_') %>%
    filter(ms == "Astral") %>%
    mutate(Metric = factor(Metric, c('IDstot', 'IDsfilt_level_1', 'IDsfilt_level_2', 'precs')))
  
  # Generate data for plotting precursors
  line <- summary_d %>%
    select(method, avgprec, precs) %>%
    mutate(method = factor(method, desired_order)) %>%
    group_by(method) %>%
    summarise(avgprec = mean(precs),
              se_prec = sd(precs) / sqrt(n()),  # Calculate standard error
              .groups = 'drop') %>%
    separate(method, c('ms', 'wind', 'agc', 'it'), sep = "_", remove = F) %>%
    filter(ms == "Astral")
  
  # Main plot
  t <- summary_d %>%
    select(method, IDstot, IDsfilt_level_1, IDsfilt_level_2, replicate) %>%
    pivot_longer(cols = c(IDstot, IDsfilt_level_1, IDsfilt_level_2),
                 names_to = 'Metric',
                 values_to = 'Values') %>%
    mutate(Metric = factor(Metric, c('IDstot', 'IDsfilt_level_1', 'IDsfilt_level_2', 'precs'))) %>%
    mutate(method = factor(method, desired_order)) %>%
    separate(method, c('ms', 'wind', 'agc', 'it'), sep = "_", remove = F) %>%
    mutate(mm = paste0(method, Metric)) %>%
    group_by(mm) %>% 
    mutate(mean = mean(Values)) %>% 
    mutate(agc = factor(agc, c(5,10,15,20,40))) %>%
    filter(ms == "Astral") %>%
    
    ggplot(aes(x = method, y = mean)) + 
    geom_bar(aes(alpha = Metric,fill=agc), colour = 'black', size = 0.5, stat = 'identity', position = 'dodge') + 
    
    geom_point(aes(x = method, y = Values, fill=agc, group = Metric), shape = 21,
               size = 1.5, position = position_jitterdodge(jitter.width = 0.06, dodge.width = 0.9), show.legend = F) +
    
    scale_fill_manual(values= fill_colors) + 
    
    scale_alpha_manual(values = c(0.1, 0.3, 0.9),
                       labels = c('Proteotypic', "5%", "1%"),
                       name = "Filtered") +
    
    scale_x_discrete(labels = x_labels, name = 'AGC Target') +
    
    coord_cartesian(expand = FALSE, xlim = c(0.3, 5.7), ylim = c(0, 922)) +
    
    theme +
    
    geom_line(data = line, aes(x = method, y = avgprec / coeff, group = 1), linewidth = 0.5) + 
    geom_point(data = line, aes(x = method, y = avgprec / coeff), size = 2.5, colour = 'red', show.legend = FALSE) +
    
    scale_y_continuous(name = "Number of Proteins",
                       sec.axis = sec_axis(~.*coeff, name = 'Number of Precursors', breaks = seq(0, 12000, by = 1500)),
                       breaks = seq(0, 2000, by = 100)) + 
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
plot_protein_precursor_errorbars(hidden_summary, mss = "Astral", coeff = 5.6, 
                                 fill_colors = c('#7BC8F6',  '#00BFFF','#30D5C8', '#00BCD4', '#CCFF00'), 
                                 fill_labels = c(), 
                                 x_labels = c('5%', '10%', '15%', '20%', '40%'),
                                 show_protein_error = F, 
                                 show_precursor_error = T)



# finale ------------------------------------------------------------------


file= 'final/report.tsv'
file= 'final/report_mag.tsv'
file= 'final/sn_report.tsv'
file= 'final/report_neat.tsv'
final_report <- generate_report(file, is_diann =T, map_concs=F)
final_report <- final_report %>%
  filter(PEP.IsProteotypic == 'True')


filts = list(
  list(pep_threshold=0.05, pg_qval =0.05),
  list(pep_threshold=0.01, pg_qval =0.01))
final_summary <- generate_summary_report(final_report, filts)
final_summary <- final_summary %>% 
  rename(method_tmp = method) %>%
  mutate(method = substr(method_tmp, 1,9))





plot_protein_precursor_errorbars <- function(summary_d, mss = "Astral", coeff = 9.3, fill_colors, fill_labels, x_labels, 
                                             show_protein_error = TRUE, show_precursor_error = TRUE) {
  
  desired_order <- c("Astral_je", "Astral_op")
  summary_d <- final_summary
  # Generate summary stats for protein data
  summary_stats <- summary_d %>%
    mutate(method = factor(method, desired_order)) %>%
    pivot_longer(cols = c(IDstot, IDsfilt_level_1, IDsfilt_level_2, precs),
                 names_to = 'Metric',
                 values_to = 'Values') %>%
    group_by(method, Metric) %>%
    summarise(mean = mean(Values),
              sd = sd(Values),
              se = sd(Values) / sqrt(n()),  # Calculate standard error
              .groups = 'drop') %>%
    separate(method, c('ms', 'meth'), remove = F, sep='_') %>%
    filter(ms == "Astral") %>%
    mutate(Metric = factor(Metric, c('IDstot', 'IDsfilt_level_1', 'IDsfilt_level_2', 'precs')))
  
  # Generate data for plotting precursors
  line <- summary_d %>%
    select(method, avgprec, precs) %>%
    mutate(method = factor(method, desired_order)) %>%
    group_by(method) %>%
    summarise(avgprec = mean(precs),
              se_prec = sd(precs) / sqrt(n()),  # Calculate standard error
              .groups = 'drop') %>%
    separate(method, c('ms', 'meth'), sep = "_", remove = F) %>%
    filter(ms == "Astral")
  
  # Main plot
  t <- summary_d %>%
    select(method, IDstot, IDsfilt_level_1, IDsfilt_level_2, replicate) %>%
    pivot_longer(cols = c(IDstot, IDsfilt_level_1, IDsfilt_level_2),
                 names_to = 'Metric',
                 values_to = 'Values') %>%
    mutate(Metric = factor(Metric, c('IDstot', 'IDsfilt_level_1', 'IDsfilt_level_2', 'precs'))) %>%
    mutate(method = factor(method, desired_order)) %>%
    separate(method, c('ms', 'meth'), sep = "_", remove = F) %>%
    mutate(mm = paste0(method, Metric)) %>%
    group_by(mm) %>% 
    mutate(mean = mean(Values)) %>% 
    mutate(meth = factor(meth, c('je', 'op'))) %>%
    filter(ms == "Astral") %>%
    
    ggplot(aes(x = method, y = mean)) + 
    geom_bar(aes(alpha = Metric,fill=meth), colour = 'black', size = 0.5, stat = 'identity', position = 'dodge') + 
    
    geom_point(aes(x = method, y = Values, fill=meth, group = Metric), shape = 21,
               size = 1.5, position = position_jitterdodge(jitter.width = 0.06, dodge.width = 0.9), show.legend = F) +
    
    scale_fill_manual(values= fill_colors, labels = fill_labels, name = "Method") + 
    
    scale_alpha_manual(values = c(0.1, 0.3, 0.9),
                       labels = c('Proteotypic', "5%", "1%"),
                       name = "Filtered") +
    
    scale_x_discrete(labels = x_labels, name = 'Method') +
    
    coord_cartesian(expand = FALSE, xlim = c(0.3, 2.7), ylim = c(0, 1100)) +
    
    theme +
    
    geom_line(data = line, aes(x = method, y = avgprec / coeff, group = 1), linewidth = 0.5) + 
    geom_point(data = line, aes(x = method, y = avgprec / coeff), size = 2.5, colour = 'red', show.legend = FALSE) +
    
    scale_y_continuous(name = "Number of Proteins",
                       sec.axis = sec_axis(~.*coeff, name = 'Number of Precursors', breaks = seq(0, 12000, by = 1500)),
                       breaks = seq(0, 2000, by = 100)) + 
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
plot_protein_precursor_errorbars(hidden_summary, mss = "Astral", coeff = 5.8, 
                                 fill_colors = c('#8F8F1F', "#007FFF"), 
                                 fill_labels = c('A1', 'Optimised'), 
                                 x_labels = c('A1', 'Optimised'),
                                 show_protein_error = F, 
                                 show_precursor_error = T)


