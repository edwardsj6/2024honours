
# init --------------------------------------------------------------------

install.packages("ggpmisc")
library(tidyverse)





# SCRIPT read files------------------------------------------------------

report_files <- list(
  "P4854/15mins/Report.tsv",
  "P4854/22mins/Report.tsv",
  "P4854/30mins/Report.tsv"
)


###global report for workspace
Global_Report <- data.frame()

#stich files together
for (file in report_files) {
  data <- read.delim(file, sep = "\t")
  Global_Report <- rbind(Global_Report, data)
}


useful_columns <- c('R.FileName', "PG.Genes", "PG.ProteinAccessions", 
                    "PG.Quantity", 'PG.IsSingleHit', 'PEP.IsProteotypic',
                    'EG.PrecursorId','PG.QValue..Run.Wise.', 'PG.PEP..Run.Wise.')
Report <- Global_Report[,useful_columns]

Report$log2 <- log(Report$PG.Quantity)
Report$sample <- substring(Report$R.FileName, 7,9)
Report$grad <- substring(Report$R.FileName, 11,12)
Report$dia <- substring(Report$R.FileName, 14,14)
Report$method <- substring(Report$R.FileName, 11,14)
Report$ID <- paste0(Report$MS, '_', substring(Report$R.FileName, 7,14))




# CV calculations ---------------------------------------------------------
##init 
Report$ID_gene <- paste0(Report$ID, "_", Report$PG.Genes)

#remove dupes from fragements
cv_report <- Report[,c('R.FileName', "PG.Genes", "log2", "method", "ID_gene")]
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
  
  data  <- cv_report[cv_report$ID_gene == ID_gene,]
  data <- na.omit(data)
  completeness <- nrow(data)/3
  prot_cv <- sd(data$log2)/mean(data$log2)
  
  #update vecs 
  t1 <- c(t1, prot_cv)
  t2 <- c(t2, ID_gene)
  t3 <- c(t3, completeness)
}

#create cv mapping and merge datasets
cv_map <- data.frame('cv' = t1, 'ID_gene' = t2, 'completeness' = t3)
Report <- merge(Report, cv_map, by = 'ID_gene', all.x = TRUE)





#CV reporting
print(length(unique(Report$ID)))
g <- ggplot(Report) + 
  geom_violin(aes(x=ID, y=cv, fill = sample)) + 
  stat_summary(aes(x = ID, y = cv),fun = "mean", width = 0.5,colour = "red", geom = "crossbar") +
  ylim(0,1)

print(g)

#================================ flag 

#map instruments
Report$Proj <- substring(Report$R.FileName, 1,5)
MS_map <- data.frame("Proj" = c('P4854', 'P4777', 'P4778'), "MS" = c('Astral', 'TimsTOF', 'Eclipse'))
Report <- merge(Report, MS_map, by = 'Proj', all.x = TRUE)                    

Report_filt <- Report %>% 
  filter(PEP.IsProteotypic == 'True')

# beta-gal ----------------------------------------------------------------

#map prp concs (and bgal concs by extension) to df

conc_map <- data.frame(sample = c('101', '102','103', '104'), conc_PRP = c(1, .5, .1,0))
Report_filt <- merge(Report_filt, conc_map, by = 'sample', all.x = TRUE)


bgal <- Report_filt[grepl("P00722", Report_filt$PG.ProteinAccessions), ]
bgal$MS_grad_sample <- paste0(bgal$MS, '_', bgal$grad, '_',bgal$sample) 
bgal$MS_grad <- paste0(bgal$MS, '_', bgal$grad, '-', bgal$dia)
MS_grads <- unique(bgal$MS_grad) 

for (MS in MS_grads) {
  tmp <- bgal[bgal$MS_grad == MS,]
  tmp <- aggregate(PG.Quantity ~ MS_grad_sample, data = tmp, FUN = mean)
  max_mean <- max(tmp$PG.Quantity)
  tmp$relative_means <- round(tmp$PG.Quantity/max_mean,3)
  
  bgal_concs <- list(tmp$relative_means)
  
  print(paste(bgal_concs, max_mean, MS))
}



bgal_qc <- Report_filt[grepl("P00722", Report_filt$PG.ProteinAccessions), ]
bgal_qc$unique_id <- paste0(bgal_qc$R.FileName, bgal_qc$PG.Quantity)
bgal_qc <- bgal_qc[!duplicated(bgal_qc$unique_id),]
bgal_qc$MS_grad_dia <- paste0(bgal_qc$MS, '_', bgal_qc$grad,'_', bgal_qc$dia)

sets <- unique(bgal_qc$MS_grad_dia)
print(sets)
for (set in sets) {
  bgal_data <- bgal_qc[bgal_qc$MS_grad_dia==set,]
  tmp_title <- paste0(set, ' - ', nrow(bgal_data))
  g <- ggplot(bgal_data, aes(x=conc_PRP, y=log2)) + geom_point(size =2) + 
    labs(title = tmp_title, y='bgal (log2)')
  print(g)
  l <- ggplot(bgal_data, aes(x=conc_PRP , y=PG.Quantity)) + geom_point(size =2)+ 
    labs(title = tmp_title, y='bgal (raw)')
  print(l) 
}



# correlation caclulations ------------------------------------------------



#correlation calculations
quant_qc <- Report_filt

#create report without dupes - select columns  ERROR PRONEERROR PRONE
qc_prot_report <- quant_qc[, c('conc_PRP', 'R.FileName','method', 'PG.Quantity', "PG.Genes", 'cv', 'MS')]
qc_prot_report$id <- paste(qc_prot_report$PG.Quantity, qc_prot_report$PG.Genes, qc_prot_report$R.FileName)
qc_prot_report <- qc_prot_report[!duplicated(qc_prot_report$id),]

#define methods - grad and windwon ("grad_window")
qc_prot_report$MS_method <- paste0(qc_prot_report$MS, '_', qc_prot_report$method)
methods <- unique(qc_prot_report$MS_method)

#define qc proteins
qc_prot <- read.delim("P4854/qc_prot.csv", sep = ',')
qc_proteins <- qc_prot$Gene.names
print(qc_proteins)

#remove GSN from QC prots
qc_proteins <- qc_proteins[qc_proteins != 'GSN']
print(qc_proteins)


### rmeove TT data since only two repsCHANGE
qc_prot_report_CHANGE <- qc_prot_report[qc_prot_report$MS!= 'TimsTOF',]
methods <- unique(qc_prot_report_CHANGE$MS_method)


#filter qc proteins detected in all 12 samples
prot_qc1 <- c() #universla protein detection across meths
prot_qc2 <- c() #unique protein detection across meths
qc_prot_hit <- data.frame(gene =qc_proteins)

for (method in methods) {
  method_data <- qc_prot_report_CHANGE[qc_prot_report_CHANGE$MS_method == method,] 
  c_vec <- c()
  for (prot in qc_proteins) {
    prot_data <- method_data[method_data$PG.Genes == prot,] 
    
    d <- NA
    #do all four samples have three replicates? If not, NA
    if (nrow(prot_data) == 12) {
       d <- 'Yes'
    }
    c_vec <- c(c_vec, d)
  }
  qc_prot_hit <- cbind(qc_prot_hit, method= c_vec)
  
}

#remove NAs to leave prots ID'd in ALL runs
tmp_data <- na.omit(qc_prot_hit)
qc_prot1 <- unique(tmp_data$gene)
qc_prot2 <- qc_proteins




#### prot 1 calculations
#calc avg R2 for methods, count no of prots wiht R2 > 0.7 and full data completeness
method_vec <- c()
R2_1_vec <- c()
prot1_count_vec <- c()

methods <- unique(qc_prot_report$MS_method)
for (method in methods) {
  method_data <- qc_prot_report[qc_prot_report$MS_method == method,]
  
  prot1_R2_avg <- 0
  prot1_R2_sum <- 0
  prot1_count <- 0
  
  for (prot in qc_prot1) {
    filtered_data <- method_data[method_data$PG.Genes == prot,] 
    filtered_data$log_2 <- log(filtered_data$PG.Quantity, 2)
    filtered_data$replicate <- substring(filtered_data$R.FileName, 17,17)
    
    #for R2, use raw intensity
    #defeine linear model and plot === ERROR PRONE === ERROR PRONE ===
    model <- lm(PG.Quantity ~ conc_PRP, data = filtered_data)
    g <- ggplot(filtered_data, aes(x=conc_PRP, y=PG.Quantity))
    g <- g + geom_point(aes(colour=replicate), size =5) + 
      stat_smooth(method = "lm", se = FALSE)
    
    
    #feedback variables - corr
    r1 <- (summary(model)$r.squared) #R2 
    r2 <- summary(model)$df[2] +2 #no. of samplpes
    r3 <- summary(model)$coefficients[2,4] #p-val

    #ONLY USE FOR DEBUGGING
    gg <- paste0(prot,"_", method,  r1)
    g <- g + labs(title = gg)
    #print(g)
    
    if (r1 > 0.7) {
      prot1_count <- prot1_count + 1 
    }
    
    prot1_R2_sum <- prot1_R2_sum + r1
  }
  
  count <- length(qc_prot1)
  prot1_R2_avg <- prot1_R2_sum/count
  #update vecs
  
  R2_1_vec <- c(R2_1_vec, prot1_R2_avg)
  method_vec <- c(method_vec, method)
  prot1_count_vec <- c(prot1_count_vec, prot1_count)
}

prot1 <- data.frame()
prot1 <- data.frame(
  method =method_vec,
  R2_1  = R2_1_vec,
  count = prot1_count_vec
)


### calc for prot 2 ###
#count no of prots with 100% completeness, adn r2 >0.7 acrpss all qc proteins
method_vec <- c()
prot2_count_vec <- c()

methods <- unique(qc_prot_report$MS_method)
r2_df <- data.frame("Genes" = qc_prot2)
for (method in methods) {
  method_data <- qc_prot_report[qc_prot_report$MS_method == method,]
  
  prot2_count <- 0 
  
  prot_vec <- c()
  r1_vec <- c()
  
  for (prot in qc_prot2) {
    filtered_data <- method_data[method_data$PG.Genes == prot,] 
    
    #check plasma protein is found/100% completeness
    if (nrow(filtered_data) ==12) {
      filtered_data$log_2 <- log(filtered_data$PG.Quantity, 2)
      filtered_data$replicate <- substring(filtered_data$R.FileName, 17,17)
      
      #defeine linear model and plot === ERROR PRONE === ERROR PRONE ===
      model <- lm(log_2 ~ PG.Quantity, data = filtered_data)
      g <- ggplot(filtered_data, aes(x=conc_PRP, y=PG.Quantity))
      g <- g + geom_point(aes(colour=replicate), size =5) + 
        stat_smooth(method = "lm", se = FALSE)
      
      #feedback variables - corr
      r1 <- (summary(model)$r.squared) #R2 
      r2 <- summary(model)$df[2] +2 #no. of samplpes
      r3 <- summary(model)$coefficients[2,4] #p-val
      
      r1_vec <- c(r1_vec, r1) 
      prot_vec <- c(prot_vec, prot)
      
      #ONLY USE FOR DEBUGGING
      gg <- paste0(prot,"_", method, '   AIC:   ' ,round(AIC(model),3), '  BIC:   ', 
                   round(BIC(model),2), '    R2:  ', r1)
      g <- g + labs(title = gg)
      #print(g)
      
      if (r1 > 0.7 & nrow(filtered_data==12)) {
        prot2_count <- prot2_count + 1
      }
    }
  }
  
  method_vec <- c(method_vec, method)
  prot2_count_vec <- c(prot2_count_vec, prot2_count) 
  
  r1 <- data.frame('Genes' = prot_vec, 'v1' = r1_vec)
  names(r1)[2] <- method
  
  r2_df <- merge(r2_df, r1, by="Genes", all=T)
}

prot2 <- data.frame()
prot2 <- data.frame(
  method =method_vec,
  prot2_count = prot2_count_vec
)


# detection ---------------------------------------------------------------

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
    select(n_unique_genes, avg_unique_genes, sd_genes)
  
  return(list(
    IDs = result$n_unique_genes,
    avgIDs = result[[1,2]],
    sdIDs = result[[1,3]]
  ))
}

#count total ids for a method
Report_filt$MS_method <- paste0(Report_filt$MS, '_', Report_filt$method)
methods <- unique(Report_filt$MS_method)
print(methods)

methods_vec <- c()

IDstot_vec <- c()
IDsfilt_vec <- c()
IDs20_vec <- c()

avgIDstot_vec <- c()
avgIDsfilt_vec <- c()
avgIDs20_vec <- c()

sdIDstot_vec <- c()
sdIDsfilt_vec <- c()
sdIDs20_vec <- c()

avgprec_vec <- c()

for (meth in methods) {
  #we have filtered data for proteypticity only, we need filtered (single hits, PEP) and CV <20%
  data <- Report_filt[Report_filt$MS_method == meth & Report_filt$sample == 104,] 
  
  dIDs <- calculate_unique_genes(data)
  IDstot_vec <- c(IDstot_vec, dIDs$IDs)
  avgIDstot_vec <- c(avgIDstot_vec, rep(dIDs$avgIDs, length(dIDs$IDs)))
  sdIDstot_vec <- c(sdIDstot_vec, rep(dIDs$sdIDs, length(dIDs$IDs)))
  
  dIDsfilt <- calculate_unique_genes(data, pep_threshold = 0.05, pg_qval = 0.05)
  IDsfilt_vec <- c(IDsfilt_vec, dIDsfilt$IDs)
  avgIDsfilt_vec <- c(avgIDsfilt_vec, rep(dIDsfilt$avgIDs, length(dIDsfilt$IDs)))
  sdIDsfilt_vec <- c(sdIDsfilt_vec, rep(dIDsfilt$sdIDs, length(dIDsfilt$IDs)))
  
  dIDs20 <- calculate_unique_genes(data, pep_threshold = 0.01, pg_qval = 0.01)
  IDs20_vec <- c(IDs20_vec, dIDs20$IDs)
  avgIDs20_vec <- c(avgIDs20_vec, rep(dIDs20$avgIDs, length(dIDs20$IDs)))
  sdIDs20_vec <- c(sdIDs20_vec, rep(dIDs20$sdIDs, length(dIDs20$IDs)))
  methods_vec <- c(methods_vec, rep(meth, 3))
  print(meth)
  
  davgprec <-data %>%
    group_by(R.FileName) %>%
    summarise(n_unique_genes = n_distinct(EG.PrecursorId)) %>%
    mutate(avg_unique_genes = mean(n_unique_genes), 
           sd_genes = sd(n_unique_genes)) %>%
    select(n_unique_genes, avg_unique_genes, sd_genes)
  avgprec_vec <- c(avgprec_vec, davgprec$n_unique_genes)
  
}


plasma_summary <- data.frame(
  method = methods_vec,
  IDstot = IDstot_vec ,
  IDsfilt = IDsfilt_vec ,
  IDs20 = IDs20_vec ,
  
  avgIDstot = avgIDstot_vec ,
  avgIDsfilt =avgIDsfilt_vec ,
  avgIDs20 =avgIDs20_vec ,
  
  sdIDstot = sdIDstot_vec ,
  sdIDsfilt =sdIDsfilt_vec ,
  sdIDs20 =sdIDs20_vec, 
  
  avgprec = avgprec_vec
)

summary <- plasma_summary

# SCRIPT - SUMMARY --------------------------------------------------------

summary <- data.frame()
summary <- merge(plasma_summary, prot1, by = 'method')
summary <- merge(summary, prot2, by = 'method')


############# notes for jack ###############
'''
expaloin why gsn is removed. panel figure of gsn conc across all methods - showing
for ms PPP and PRP asnalusis it is not a good protein.

think about how you can show these methods quant perofrance, of all 7 proteins acorss each run, 
do you show avg corr of all 7, or no of IDs where R2 > 0.7. 

then for qc prot 2, how could  represent that, maube as simple as showing no of 
additional ids wiht r2 > 0.7

'''


# SCRIPT - other factors  -------------------------------------------------

#explore elution times for peptides FIX FOR NEW REPORT FORMAT
g <- ggplot(Report, aes(x=EG.ApexRT)) + geom_histogram()
print(g)

  

             

