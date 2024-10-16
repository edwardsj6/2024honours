
# formal analysis  --------------------------------------------------------
install.packages("wesanderson")
install.packages("RColorBrewer")
install.packages("tidyverse")
install.packages("ggrepel")
library("ggplot2")
library("RColorBrewer")
library(ggrepel)
library("ggridges")

library('tidyverse')

display.brewer.all(colorblindFriendly = T)
#take dark 2, 15=2 22=3 30=4 45=5

theme <- theme_classic() + theme(
  text = element_text((family = "Arial"), size = 12, colour ="black"),  # Set all text to Arial, size 12
  plot.title = element_text((family = "Arial"), size = 14),  # Adjust title size (optional)
  #strip.text.x = element_text((family = "Arial"), size = 12),  # Adjust strip text size (optional)
  axis.text = element_text((family = "Arial"), size = 12, colour ="black"),
  legend.text = element_text((family = "Arial"), size = 12, colour ="black")
  # ... adjust other text elements as needed
)

#another test
g 
# correlation analysis ----------------------------------------------------



library(ggplot2)
library(reshape)


#### global correlation 

cor_summary <- na.omit(summary[, -c(1)])
cormat <- round(cor(cor_summary, method='spearman'),2)
# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

upper_tri <- get_upper_tri(cormat)



melted_cormat <- melt(upper_tri)

ggheatmap <- ggplot(melted_cormat, aes(X2, X1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Spearman\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()


ggheatmap <- ggheatmap + 
  geom_text(aes(X2, X1, label = value), color = "black", size = 4) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position.inside  = c(0.6, 0.7),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5)) +
  labs(title = 'GLobal')
ggheatmap



#### correlation by isntrument
MSs <- c('Eclipse', 'Astral', 'TimsTOF')

for (MS in MSs) {
  cor_summary <- summary[grepl(MS, summary$method),]
  cor_summary <- na.omit(cor_summary[,-c(1)])
  cormat <- round(cor(cor_summary, method='pearson'),2)
  
  # Get lower triangle of the correlation matrix
  get_lower_tri<-function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
  }
  # Get upper triangle of the correlation matrix
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }
  
  upper_tri <- get_upper_tri(cormat)
  
  
  
  melted_cormat <- melt(upper_tri)
  
  ggheatmap <- ggplot(melted_cormat, aes(X2, X1, fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(-1,1), space = "Lab", 
                         name="Spearman\nCorrelation") +
    theme_minimal()+ # minimal theme
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 12, hjust = 1))+
    coord_fixed()
  
  
  ggheatmap <- ggheatmap + 
    geom_text(aes(X2, X1, label = value), color = "black", size = 4) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank(),
      legend.justification = c(1, 0),
      legend.position.inside  = c(0.6, 0.7),
      legend.direction = "horizontal")+
    guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                                 title.position = "top", title.hjust = 0.5)) +
    labs(title = MS)
  ggheatmap
  
}

# CVs ---------------------------------------------------------------------

#make cv pltos for all 

g <- ggplot(Report) + 
  geom_violin(aes(x=ID, y=cv.y, fill = sample)) + 
  stat_summary(aes(x = ID, y = cv.y),fun = "mean", width = 0.5,colour = "red", geom = "crossbar") +
  ylim(0,1)

print(g)


#make for each gradient
Report$MS_grad <- paste0(Report$MS, '_', Report$grad)
MS_grads <- unique(Report$MS_grad)
for (grad in MS_grads) {
  data <- Report[Report$MS_grad == grad,]
  
  g <- ggplot(data) + 
    geom_violin(aes(x=ID, y=cv, fill = sample)) + 
    stat_summary(aes(x = ID, y = cv),fun = "mean", width = 0.5,colour = "red", geom = "crossbar") +
    ylim(0,1) + labs(title = grad)+
    theme(axis.text.x = element_text(angle=45))
  
  print(g)
  
}


#make for each Instrument

MSs <- unique(Report$MS)
for (MS in MSs) {
  data <- Report[Report$MS == MS,]
  
  g <- ggplot(data) + 
    geom_violin(aes(x=ID, y=cv, fill = sample)) + 
    stat_summary(aes(x = ID, y = cv),fun = "mean", width = 0.5,colour = "red", geom = "crossbar") +
    ylim(0,1) + labs(title = MS)+
    theme(axis.text.x = element_text(angle=45))
  
  print(g)
  
}


#make for PLASMA ONLY

data <- Report[Report$sample == 104,]
data$lab <- substring(data$method)
  
g <- ggplot(data) +
  geom_violin(aes(x=ID, y=cv, fill = MS)) + 
  stat_summary(aes(x = ID, y = cv),fun = "mean", width = 0.5,colour = "red", geom = "crossbar") +
  ylim(0,1) + labs(title = 'Plasma')  +
  theme(axis.text.x = element_text(angle=45))
  
print(g)
  
#mutate(label = substring(ID, nchar(ID) -4, nchar(ID)))


# heatmap for summary -----------------------------------------------------

data <- summary

data$MS <-substring(data$method, 1,nchar(data$method)-5)
data$grad <- substring(data$method, nchar(data$method)-3,nchar(data$method)-2)                     
data$dia <- substring(data$method, nchar(data$method),nchar(data$method))

data$MS_grad<- substring(data$method, 1,nchar(data$method)-2)



#### all IDs 
data_d <- data
data_d$X3 <- data_d$IDs
data_d <- data_d[,c("MS_grad", "dia", "X3")]

melted_cormat <- data_d
ggheatmap <- ggplot(melted_cormat, aes(MS_grad, dia, fill = X3))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",midpoint = 500,
                       limits = c(0, 1200),
                       space = "Lab", 
                       name="") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()


ggheatmap <- ggheatmap + 
  geom_text(aes(MS_grad, dia, label = X3), color = "black", size = 4, fontface='bold') +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position.inside  = c(0.6, 0.7),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5)) +
  labs(title = 'IDs10')
ggheatmap




#### IDs10 

data_d <- data
data_d$X3 <- data_d$IDs10
data_d <- data_d[,c("MS_grad", "dia", "X3")]

melted_cormat <- data_d
ggheatmap <- ggplot(melted_cormat, aes(MS_grad, dia, fill = X3))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",midpoint = 500,
                       limits = c(0, 1200),
                       space = "Lab", 
                       name="") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()


ggheatmap <- ggheatmap + 
  geom_text(aes(MS_grad, dia, label = X3), color = "black", size = 4, fontface='bold') +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position.inside  = c(0.6, 0.7),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5)) +
  labs(title = 'IDs10')
ggheatmap


#### IDs5 
data_d <- data
data_d$X3 <- data_d$IDs5
data_d <- data_d[,c("MS_grad", "dia", "X3")]

melted_cormat <- data_d
ggheatmap <- ggplot(melted_cormat, aes(MS_grad, dia, fill = X3))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",midpoint = 500,
                       limits = c(0, 1200),
                       space = "Lab", 
                       name="") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()


ggheatmap <- ggheatmap + 
  geom_text(aes(MS_grad, dia, label = X3), color = "black", size = 4, fontface='bold') +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position.inside  = c(0.6, 0.7),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5)) +
  labs(title = 'IDs5')
ggheatmap


#### correlation 
data_d <- data
data_d$corr <- round(data_d$corr, 2)
data_d$X3 <- data_d$corr
data_d <- data_d[,c("MS_grad", "dia", "X3")]

melted_cormat <- data_d
ggheatmap <- ggplot(melted_cormat, aes(MS_grad, dia, fill = X3))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",midpoint = .9,
                       limits = c(.8, 1),
                       space = "Lab", 
                       name="") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()


ggheatmap <- ggheatmap + 
  geom_text(aes(MS_grad, dia, label = X3), color = "black", size = 4, fontface='bold') +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position.inside  = c(0.6, 0.7),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5)) +
  labs(title = 'corr')
ggheatmap

#### ids of qc proteins
data_d <- data

data_d$X3 <- data_d$no_id
data_d <- data_d[,c("MS_grad", "dia", "X3")]

melted_cormat <- data_d
ggheatmap <- ggplot(melted_cormat, aes(MS_grad, dia, fill = X3))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",midpoint = 15,
                       limits = c(10, 20),
                       space = "Lab", 
                       name="") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()


ggheatmap <- ggheatmap + 
  geom_text(aes(MS_grad, dia, label = X3), color = "black", size = 4, fontface='bold') +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position.inside  = c(0.6, 0.7),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5)) +
  labs(title = 'corr')
ggheatmap




# PCA ---------------------------------------------------------------------

install.packages('ggfortify')
library(ggfortify)


pca <- na.omit(summary)
pca$ms <- substring(pca$method, 1,nchar(pca$method)-5)
pca$grad <- substring(pca$method,nchar(pca$method)-3,nchar(pca$method)-2)
pca$dia <- substring(pca$method, nchar(pca$method),nchar(pca$method))
pca_res <- pca[2:6]
pca_res  <- prcomp(pca_res, scale. = T)
autoplot(pca_res, data = pca, colour = 'ms')






# CV graphs ---------------------------------------------------------------

#from orginal dataset, by instrument and rank for each cv
#remove dupes gene_cv_method for only plasma

rCV <- Report %>%
  filter(MS=='Astral', sample==104) %>%
  select(cv,method,PG.Genes,MS) %>%
  distinct() %>%
  na.omit() 

rCV %>% 
  mutate(cv1 = 100*cv) %>%
  ggplot(aes(x=cv1, y=method, fill= method)) +
  scale_fill_manual(
    values = c(
      brewer.pal(n=8, name = 'OrRd')[c(4,6,7)], #grad 15, 1,2,3
      brewer.pal(n=8, name = 'RdPu')[c(4,6,7)], #grad 22, 1,2,3
      brewer.pal(n=8, name = 'Greens')[c(4,6,7)] #grad 30, 1,2,3
    )
  )+
  geom_density_ridges(scale=3) +
  xlim(0,15) + 
  coord_cartesian(expand=F)+
  theme


#cycle through methods to calc rank for cvs
methods <- unique(rCV$method)
rCV_final <- data.frame()
for (meth in methods) {
  tmp <- rCV %>%
    filter(method ==meth) %>%
    arrange(cv) %>%
    mutate(rank=row_number())
  rCV_final <- rbind(rCV_final, tmp)
  
}


rCV_final %>%
  mutate(dia = substring(method,4,4), grad = substring(method,1,2)) %>%
  ggplot(aes(x=cv,y=rank, colour=interaction(dia,grad))) + 
  geom_point()+
  xlim(0,.15) +
  labs(x= 'cv 0.09') +
  theme + 
  scale_colour_manual(
    values = c(
      brewer.pal(n=8, name = 'OrRd')[c(4,6,7)], #grad 15, 1,2,3
      brewer.pal(n=8, name = 'RdPu')[c(4,6,7)], #grad 22, 1,2,3
      brewer.pal(n=8, name = 'Greens')[c(4,6,7)] #grad 30, 1,2,3
    )
    #labels = c("15", "22" , "30")
  ) 


#### %%%%%%%%%%%
#### Eclipse
#### %%%%%%%%%%%

rCV <- Report %>%
  filter(MS=='Eclipse', sample==104) %>%
  select(cv,method,PG.Genes,MS) %>%
  distinct() %>%
  na.omit() 

rCV %>% 
  mutate(cv1 = 100*cv) %>%
  ggplot(aes(x=cv1, y=method, fill= method)) +
  scale_fill_manual(
    values = c(
      brewer.pal(n=8, name = 'OrRd')[c(4,6,7)], #grad 15, 1,2,3
      brewer.pal(n=8, name = 'Greens')[c(4,6,7)], #grad 22, 1,2,3
      brewer.pal(n=8, name = 'YlOrBr')[c(2,3,4)] #grad 30, 1,2,3
    )
  )+
  geom_density_ridges(scale=3) +
  xlim(0,15) + 
  coord_cartesian(expand=F)+
  theme

#cycle through methods to calc rank for cvs
methods <- unique(rCV$method)
rCV_final <- data.frame()
for (meth in methods) {
  tmp <- rCV %>%
    filter(method ==meth) %>%
    arrange(cv) %>%
    mutate(rank=row_number())
  rCV_final <- rbind(rCV_final, tmp)
  
}


rCV_final %>%
  mutate(dia = substring(method,4,4), grad = substring(method,1,2)) %>%
  ggplot(aes(x=cv,y=rank, colour=interaction(dia,grad))) + 
  geom_point()+
  xlim(0,.15) +
  labs(x= 'cv 0.09') +
  theme + 
  scale_colour_manual(
    values = c(
      brewer.pal(n=8, name = 'OrRd')[c(4,6,7)], #grad 15, 1,2,3
      brewer.pal(n=8, name = 'RdPu')[c(4,6,7)], #grad 22, 1,2,3
      brewer.pal(n=8, name = 'Greens')[c(4,6,7)] #grad 30, 1,2,3
    )
    #labels = c("15", "22" , "30")
  ) 








#### %%%%%%%%%%%
#### TimsTOF
#### %%%%%%%%%%%

rCV <- Report %>%
  filter(MS=='TimsTOF', sample==104) %>%
  select(cv,method,PG.Genes,MS) %>%
  distinct() %>%
  na.omit() 

rCV %>% 
  mutate(cv1 = 100*cv) %>%
  ggplot(aes(x=cv1, y=method, fill= method)) +
  scale_fill_manual(
    values = c(
      brewer.pal(n=8, name = 'OrRd')[c(4,6,7)], #grad 15, 1,2,3
      brewer.pal(n=8, name = 'Greens')[c(4,6,7)], #grad 22, 1,2,3
      brewer.pal(n=8, name = 'YlOrBr')[c(2,3,4)] #grad 30, 1,2,3
    )
  )+
  geom_density_ridges(scale=3) +
  xlim(0,15) + 
  coord_cartesian(expand=F)+
  theme

#cycle through methods to calc rank for cvs
methods <- unique(rCV$method)
rCV_final <- data.frame()
for (meth in methods) {
  tmp <- rCV %>%
    filter(method ==meth) %>%
    arrange(cv) %>%
    mutate(rank=row_number())
  rCV_final <- rbind(rCV_final, tmp)
  
}


rCV_final %>%
  mutate(dia = substring(method,4,4), grad = substring(method,1,2)) %>%
  ggplot(aes(x=cv,y=rank, colour=interaction(dia,grad))) + 
  geom_point()+
  xlim(0,.15) +
  labs(x= 'cv 0.09') +
  theme + 
  scale_colour_manual(
    values = c(
      brewer.pal(n=8, name = 'OrRd')[c(4,6,7)], #grad 15, 1,2,3
      brewer.pal(n=8, name = 'RdPu')[c(4,6,7)], #grad 22, 1,2,3
      brewer.pal(n=8, name = 'Greens')[c(4,6,7)] #grad 30, 1,2,3
    )
    #labels = c("15", "22" , "30")
  ) 



# detection plasma  --------------------------------------------------------
brewer.pal(n = 8, name = "Dark2")
color_value <- brewer.pal(n = 8, name = "Set2")[1]
print(color_value)



ttt <- read.csv("P4777/Book1.csv", header=F)
ttt %>% 
  filter(V1=='timstof') %>%
  mutate(method = paste(V1, V2, V3),
         grad=as.factor(V2)) %>%
  ggplot(aes(x=method, y=V4))+
  geom_line(aes(group=grad, colour = grad), show.legend = F) +
  geom_point(aes(colour = grad),show.legend = F, size = 3) +
  scale_colour_manual(values = c('#C9AA81', "#8A360F",'#4B0082'), 
                      labels = c('15 mins', '22 mins', '30 mins'),
                      name = 'Gradient') +
  scale_y_continuous(limits=c(8,14), breaks = seq(0,20, by =1)) +
  theme +
  labs(y='Precursors/Protein')+
  scale_x_discrete(labels = c('T1', 'T2', 'T3', 'T1', 'T2', 'T3', 'T1', 'T2', 'T3'), name = NULL) 


fill_colors = c('#C9AA81', "#8A360F",'#4B0082')
values = c('#C9AA81', "#808000", "#8A360F")
prop_report %>%
  separate(method, c('ms', 'grad', 'dia'), remove = FALSE) %>%
  mutate(dia, factor(dia, c(1,2,3))) %>%
  filter(ms == 'Astral') %>%
  
  # Plot with points and lines for level2
  ggplot(aes(x = method, y = level3)) + 
  geom_line(aes(group = grad, colour = grad)) +
  geom_point(aes(colour = grad), size = 3) +
  
  # Add color scale for gradient
  scale_colour_manual(values = c('#C9AA81', "#808000", "#8A360F"), 
                      labels = c('15 mins', '22 mins', '30 mins'),
                      name = 'Gradient') +
  scale_x_discrete(labels = c('A1', 'A2', 'A3', 'A1', 'A2', 'A3', 'A1', 'A2', 'A3'), name = NULL) 
  
  




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

test_report <- generate_report(files, is_diann=F, map_concs=T)
test_report_backup <- test_report
unique(test_report$MS)
filters <- list(
  list(pep_threshold = .05, pg_qval = 0.05), 
  list(pep_threshold = 0.01, pg_qval = 0.01)
)
summary <- generate_summary_report(quant_report, filters)

result <- test_report %>%
  filter(conc_PRP %in% c(0, 1)) %>%  # Filter for conc_PRP 0 and 1
  select(PG.Genes, R.FileName, conc_PRP) %>%
  group_by(R.FileName, conc_PRP) %>%
  summarise(count = n_distinct(PG.Genes)) %>%
  ungroup() %>%
  separate(R.FileName, c('v1', 'v2', 'v3', 'v4', 'v5')) %>%
  mutate(method = paste(v1, v3)) %>%  # Create method identifier
  group_by(method, conc_PRP) %>%
  summarise(count = mean(count)) %>%  # Sum counts for each method and conc_PRP
  spread(conc_PRP, count) %>%  # Reshape to have conc_PRP values as columns
  mutate(ratio = `1` / `0`)  # Calculate the ratio (conc_PRP 1 / conc_PRP 0)

# Print the result
print(mean(result$ratio))
  
library(dplyr)
library(eulerr)

# Filter the data for conc_PRP == 0 and three methods (Astral, timsTOF, Eclipse)
venn_data <- test_report %>%
  filter(conc_PRP == 0, MS %in% c('Astral', 'TimsTOF', 'Eclipse'),
         PEP.IsProteotypic=="True") %>%
  select(PG.Genes, MS) %>%
  distinct()

ttt <- test_report %>%
  filter(conc_PRP == 0, PG.Genes %in% astral_proteins, PEP.IsProteotypic=='True') %>%
  arrange(log2) %>%
  select(PG.Genes, log2) %>%
  distinct()

l <- test_report %>% 
  filter(PG.Genes %in% only_astral) %>%
  select(PG.Genes, log2) %>%
  distinct() %>%
  group_by(PG.Genes) %>%
  summarise(mean = mean(log2))
  
# Create the sets for each method
astral_proteins <- venn_data %>% filter(MS == 'Astral') %>% pull(PG.Genes)
timstof_proteins <- venn_data %>% filter(MS == 'TimsTOF') %>% pull(PG.Genes)
eclipse_proteins <- venn_data %>% filter(MS == 'Eclipse') %>% pull(PG.Genes)

# Identify the overlaps between the methods
only_astral <- setdiff(astral_proteins, union(timstof_proteins, eclipse_proteins))
only_timstof <- setdiff(timstof_proteins, union(astral_proteins, eclipse_proteins))
only_eclipse <- setdiff(eclipse_proteins, union(astral_proteins, timstof_proteins))

astral_timstof <- intersect(astral_proteins, timstof_proteins) %>% setdiff(eclipse_proteins)
astral_eclipse <- intersect(astral_proteins, eclipse_proteins) %>% setdiff(timstof_proteins)
timstof_eclipse <- intersect(timstof_proteins, eclipse_proteins) %>% setdiff(astral_proteins)

all_three <- intersect(intersect(astral_proteins, timstof_proteins), eclipse_proteins)

# Create a named numeric vector for eulerr
protein_counts <- c(
  "Astral" = length(only_astral),
  "timsTOF" = length(only_timstof),
  "Eclipse" = length(only_eclipse),
  "Astral&timsTOF" = length(astral_timstof),
  "Astral&Eclipse" = length(astral_eclipse),
  "timsTOF&Eclipse" = length(timstof_eclipse),
  "Astral&timsTOF&Eclipse" = length(all_three)
)

# Generate the Euler diagram
fit <- euler(protein_counts)

# Plot the Euler diagram
plot(fit, quantities = TRUE, fill = c("#FFDDC1", "#C1EFFF", "#FFC1DD"),
     labels = NULL, main = NULL)

unique_genes = only_timstof
test_report %>% 
  filter(MS=="TimsTOF", PEP.IsProteotypic=="True", 
         PG.IsSingleHit == 'False', conc_PRP==0) %>%
  select(PG.Genes, log2) %>%
  group_by(PG.Genes) %>% 
  summarise(avg = mean(log2, na.rm = T)) %>%
  arrange(desc(avg)) %>%
  mutate(rank = row_number(), unique = case_when(
    PG.Genes %in% unique_genes ~ T, 
    T ~F
  )) %>%
  arrange(unique) %>%
  ggplot(aes(x=rank, y=avg)) + 
  geom_point(aes(colour=unique)) + 
  scale_colour_manual(values = c('black','red'), labels = c('False', 'True'))+
  geom_text_repel(aes(
    label = ifelse(PG.Genes %in% unique_genes, PG.Genes, "")), 
    size = 4, color = "red", 
    max.overlaps = Inf, 
    min.segment.length = 0, 
    box.padding = .5,  # Increased padding around text box
    point.padding = 1,
    force = 1# Padding between points and text labels
  )+
  theme  + 
  scale_y_continuous(breaks=seq(0,30,by=5)) + 
  scale_x_continuous(breaks=seq(0,2000, by=250))+
  labs(x = "Rank", 
       y = expression("Protein Intensity (log"[2]*")"), color = "Unique\nIdentification")





test_report %>%
  distinct() %>% 
  separate(R.FileName, c('v1', 'v2', 'v3', 'v4', 'v5')) %>%
  mutate(method_test = paste0(v3,'_', v4,'_', v5)) %>%
  group_by(method_test, conc_PRP) %>%
  summarise(avg_log2 = mean(PG.Quantity)) %>%
  ungroup() %>%
  
  # Calculate relative expression by setting conc_PRP = 0 to 1
  group_by(method_test) %>%
  mutate(relative_expression = ifelse(conc_PRP == 0, 1, avg_log2 / avg_log2[conc_PRP == 0])) %>%
  ungroup() %>%
  
  # Plot relative expression
  ggplot(aes(x = conc_PRP, y = relative_expression colour = as.factor(v3)) +
  geom_point() + 
  geom_abline(slope = -1, intercept = 1, linetype = "dashed", color = "red") +
  labs(y = 'Relative Protein Expression', x = 'Proportion PRP', colour = 'Method') + 
  theme
  

summary_d <- test_summary
##### %%%%%%%%%%%%%%
##### ASTRAL Plasma
##### %%%%%%%%%%%%%%

summary_stats <- summary_d %>%
  pivot_longer(cols = c(IDstot, IDsfilt_level_1, IDsfilt_level_2, precs),
               names_to = 'Metric',
               values_to = 'Values') %>%
  group_by(method, Metric) %>%
  summarise(mean = mean(Values),
            sd = sd(Values),
            se = sd(Values) / sqrt(n()),  # Calculate standard error
            .groups = 'drop') %>%
  separate(method, c('ms', 'grad', 'dia'), remove =F)  %>%
  filter(ms=='Astral')  %>%
  mutate(Metric = factor(Metric, c('IDstot', 'IDsfilt_level_1', 'IDsfilt_level_2', 'precs')))



# Generate data for plotting
line <- summary_d %>%
  select(method, avgprec, precs) %>%
  group_by(method) %>%
  summarise(avgprec = mean(precs),
            se_prec = sd(precs) / sqrt(n()),  # Calculate standard error
            .groups = 'drop') %>%
  separate(method, c('ms', 'grad', 'dia'), sep = "_", remove = F) %>%
  filter(ms == "Astral")

coeff = 9.3

# Ensure Metric is consistently available and correctly grouped
t <- summary_d %>%
  select(method, IDstot, IDsfilt_level_1, IDsfilt_level_2, replicate) %>%
  pivot_longer(cols = c(IDstot, IDsfilt_level_1, IDsfilt_level_2),
               names_to = 'Metric',
               values_to = 'Values') %>%
  mutate(Metric = factor(Metric, c('IDstot', 'IDsfilt_level_1', 'IDsfilt_level_2'))) %>%
  separate(method, c('ms', 'grad', 'dia'), sep="_",remove=F) %>%
  mutate(mm = paste0(method,Metric)) %>%
  group_by(mm) %>% 
  mutate(mean=mean(Values)) %>% 
  filter(ms == "Astral") %>%
  
  
  ggplot(aes(x=method, y=mean, fill=grad)) + 
  geom_bar(aes(alpha=Metric), colour='black', size=.5, stat='identity', position = 'dodge') + 
  
  geom_point(aes(x=method, y=Values, fill=grad, group=Metric), shape=21,
             size=1.5, position=position_jitterdodge(
               jitter.width=.06,dodge.width=.9), show.legend = T) +
  geom_errorbar(data = summary_stats, 
                aes(x = as.numeric(as.factor(method)) + (as.numeric(Metric) - 2) * 0.3, 
                   group = Metric, ymin = mean - se, ymax = mean + se),
               width = 0.13, size=.5) +
  scale_fill_manual(
    values = c('#D95F02', "#E7298A", "#66A61E" ),
    labels = c("22 mins", "30 mins", "30 mins"),
    name = "Gradient"
  )+ 
  
  
  
  scale_alpha_manual(
    values  = c(0.1, 0.25, 0.5),
    labels = c('Total', "Filtered, FDR < 5%", "Filtered, FDR < 1%"),
    name = "Protein ID's"
  )+
  scale_x_discrete(
    labels = c('A1', 'A2', 'A3',
               'A1', 'A2', 'A3',
               'A1', 'A2', 'A3'),
    name = NULL
  ) +
  labs( title = "Astral")+
  
  coord_cartesian(expand = F,xlim=c(0.3, 9.7), ylim=c(0,1200)) +
  theme + 
  #add precursros sec axis
  geom_line(data=line, aes(x=method, y=avgprec/coeff, group =grad), linewidth =.5)+ 
  geom_point(data=line, aes(x=method, y=avgprec/coeff), size =2.5, colour ='red', show.legend = F) +
  # Add the precursor error bars to the secondary y-axis
  geom_errorbar(data = line, 
               aes(x = method, 
                    y = avgprec / coeff, 
                    ymin = (avgprec - se_prec) / coeff, 
                  ymax = (avgprec + se_prec) / coeff), 
               width = 0.08, size=.5, show.legend = F) +
  
  scale_y_continuous(
    name = "Number of Proteins",
    sec.axis = sec_axis(~.*coeff, name = 'Number of Precursors', breaks = seq(0,12000, by=1500)),
    breaks = seq(0, 1200, by = 100)
  ) + 
  labs()+
  theme(axis.text.y.right = element_text(colour='red'),
        axis.title.y.right = element_text(colour='red'),
        axis.line.y.right = element_line(colour='red'),
        axis.ticks.y.right = element_line(colour='red'))

print(t)


# teast gpt ---------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(tidyr)

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
               size = 1.5, position = position_jitterdodge(jitter.width = 0.06, dodge.width = 0.9), show.legend = F) +
    
    scale_fill_manual(values = fill_colors,
                      labels = fill_labels,
                      name = "Gradient") + 
    
    scale_alpha_manual(values = c(0.1, 0.25, 0.5),
                       labels = c('Proteotypic', "5%", "1%"),
                       name = "Filtered") +
    
    scale_x_discrete(labels = x_labels, name = NULL) +
    
    labs(title = ms) +
    
    coord_cartesian(expand = FALSE, xlim = c(0.3, 9.7), ylim = c(0, 1200)) +
    
    theme +
    
    geom_line(data = line, aes(x = method, y = avgprec / coeff, group = grad), linewidth = 0.5) + 
    geom_point(data = line, aes(x = method, y = avgprec / coeff), size = 2.5, colour = 'red', show.legend = FALSE) +
    
    scale_y_continuous(name = "Number of Proteins",
                       sec.axis = sec_axis(~.*coeff, name = 'Number of Precursors', breaks = seq(0, 12000, by = 1500)),
                       breaks = seq(0, 1200, by = 100)) + 
    labs(title=NULL) +
    
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
                                 fill_colors = c('#C9AA81', "#808000", "#8A360F"), 
                                 fill_labels = c("15 mins", "22 mins", "30 mins"), 
                                 x_labels = c('A1', 'A2', 'A3', 'A1', 'A2', 'A3', 'A1', 'A2', 'A3'),
                                 show_protein_error = F, 
                                 show_precursor_error = TRUE)

plot_protein_precursor_errorbars(summary, ms = "Eclipse", coeff = 9.3, 
                                 fill_colors = c('#C9AA81', "#8A360F",'#4B0082'), 
                                 fill_labels = c("15 mins", "30 mins", "45 mins"), 
                                 x_labels = c('E1', 'E2', 'E3', 'E1', 'E2', 'E3', 'E1', 'E2', 'E3'),
                                 show_protein_error = F, 
                                 show_precursor_error = TRUE)


test <- create_plot(
  summary_d = summary_d,
  ms_filter = "Astral",
  coeff = 9.3,
  fill_colors = c('#D95F02', "#E7298A", "#66A61E"),
  fill_labels = c("22 mins", "30 mins", "30 mins"),
  x_labels = c('A1', 'A2', 'A3', 'A1', 'A2', 'A3', 'A1', 'A2', 'A3'),
  add_precursor_error = TRUE,
  add_protein_error = TRUE
)



# cont --------------------------------------------------------------------



##### %%%%%%%%%%%%%%
###### ECLIPSE Plasma
##### %%%%%%%%%%%%%%
line <- summary_d %>%
  select(method, avgprec) %>%
  distinct() %>% #since avg values
  separate(method, c('ms', 'grad', 'dia'), sep="_",remove=F) %>%
  filter(ms == "Eclipse")
coeff = 10

summary_d %>%
  select(method, IDstot, IDsfilt, IDs20) %>%
  pivot_longer(cols = c(IDstot, IDsfilt, IDs20),
               names_to = 'Metric',
               values_to = 'Values') %>%
  mutate(Metric = factor(Metric, c('IDstot', 'IDsfilt', 'IDs20'))) %>%
  separate(method, c('ms', 'grad', 'dia'), sep="_",remove=F) %>%
  mutate(mm = paste0(method,Metric)) %>%
  group_by(mm) %>% 
  mutate(mean=mean(Values)) %>% 
  filter(ms == "Eclipse") %>%
  
  
  ggplot(aes(x=method, y=mean, fill=grad)) + 
  geom_bar(aes(alpha=Metric), colour='black', size=.5, stat='identity', position = 'dodge') + 
  geom_point(aes(x=method, y=Values, fill=grad, group=Metric), 
             shape=21, size=1.5, position=position_jitterdodge(
               jitter.width=.06,dodge.width=.9), show.legend = F) +
  scale_fill_manual(
    values = c('#D95F02', "#66A61E", "#E6AB02"),
    labels = c("15 mins", "30 mins", "45 mins"),
    name = "Gradient"
  )+ 
  scale_alpha_manual(
    values  = c(0.1, 0.25, 0.5),
    labels = c('Total', "Filtered", "CV < 20%"),
    name = "Protein ID's"
  )+
  scale_x_discrete(
    labels = c('E1', 'E2', 'E3',
               'E1', 'E2', 'E3',
               'E1', 'E2', 'E3'),
    name = NULL
  ) +
  labs( title = "Eclipse")+
  
  coord_cartesian(expand = F,xlim=c(0.3, 9.7), ylim=c(0,1200)) +
  theme + 
  #add precursros sec axis
  geom_line(data=line, aes(x=method, y=avgprec/coeff, group =grad), linewidth =.5)+ 
  geom_point(data=line, aes(x=method, y=avgprec/coeff), size =2.5, colour ='red', show.legend = F) +
  
  scale_y_continuous(
    name = "Number of Proteins",
    sec.axis = sec_axis(~.*coeff, name = 'Number of Precursors', breaks = seq(0,12000, by=1500)),
    breaks = seq(0, 1200, by = 100)
  ) + 
  labs()+
  theme(axis.text.y.right = element_text(colour='red'),
        axis.title.y.right = element_text(colour='red'),
        axis.line.y.right = element_line(colour='red'),
        axis.ticks.y.right = element_line(colour='red'))





##### %%%%%%%%%%%%%%
###### TimsTOF Plasma 
##### %%%%%%%%%%%%%%

line <- summary_d %>%
  select(method, avgprec) %>%
  distinct() %>% #since avg values
  separate(method, c('ms', 'grad', 'dia'), sep="_",remove=F) %>%
  filter(ms == "TimsTOF")
coeff = 10

summary_d %>%
  select(method, IDstot, IDsfilt, IDs20) %>%
  pivot_longer(cols = c(IDstot, IDsfilt, IDs20),
               names_to = 'Metric',
               values_to = 'Values') %>%
  mutate(Metric = factor(Metric, c('IDstot', 'IDsfilt', 'IDs20'))) %>%
  separate(method, c('ms', 'grad', 'dia'), sep="_",remove=F) %>%
  mutate(mm = paste0(method,Metric)) %>%
  group_by(mm) %>% 
  mutate(mean=mean(Values)) %>% 
  filter(ms == "TimsTOF") %>%
  
  
  ggplot(aes(x=method, y=mean, fill=grad)) + 
  geom_bar(aes(alpha=Metric), colour='black', linewidth=.5, stat='identity', position = 'dodge') + 
  geom_point(aes(x=method, y=Values, fill=grad, group=Metric), 
             shape=21, size=1.5, position=position_jitterdodge(
               jitter.width=.06,dodge.width=.9), show.legend = F) +
  scale_fill_manual(
    values = c('#D95F02', "#66A61E", "#E6AB02"),
    labels = c("15 mins", "30 mins", "45 mins"),
    name = "Gradient"
  )+ 
  scale_alpha_manual(
    values  = c(0.1, 0.25, 0.5),
    labels = c('Total', "Filtered", "CV < 20%"),
    name = "Protein ID's"
  )+
  scale_x_discrete(
    labels = c('T1', 'T2', 'T3',
               'T1', 'T2', 'T3',
               'T1', 'T2', 'T3'),
    name = NULL
  ) +
  labs( title = "TimsTOF")+
  
  coord_cartesian(expand = F,xlim=c(0.3, 9.7), ylim=c(0,1200)) +
  theme + 
  #add precursros sec axis
  geom_line(data=line, aes(x=method, y=avgprec/coeff, group =grad), linewidth =.5)+ 
  geom_point(data=line, aes(x=method, y=avgprec/coeff), size =2.5, colour ='red', show.legend = F) +
  
  scale_y_continuous(
    name = "Number of Proteins",
    sec.axis = sec_axis(~.*coeff, name = 'Number of Precursors', breaks = seq(0,12000, by=1500)),
    breaks = seq(0, 1200, by = 100)
  ) + 
  labs()+
  theme(axis.text.y.right = element_text(colour='red'),
        axis.title.y.right = element_text(colour='red'),
        axis.line.y.right = element_line(colour='red'),
        axis.ticks.y.right = element_line(colour='red'))


# rank abundance plots ----------------------------------------------------

#highlight qc_proteins
highlight_qc_proteins <- function(data, qc_proteins_input = NULL) {
  # Define a default list of QC proteins
  default_qc_proteins <- c('FLNA', 'TLN1', 'MYH9', 'ACTB', 'VCL', 'ACTN1', 'TPM4', 'THBS1', 'TUBB1', 'YWHAZ', 'GSN', 'TUBA1B', 'ITGA2B', 'F13A1', 'PFN1', 'ITGB3', 'TAGLN2', 'FERMT3', 'RAP1B', 'PLEK', 'PPBP', 'GAPDH', 'MMRN1', 'MYL6', 'CFL1', 'PARVB', 'SDPR', 'TUBB4B', 'TMSB4X', 'PKM')
  
  # Use the default QC proteins if none are provided
  if (is.null(qc_proteins_input)) {
    qc_proteins_input <- default_qc_proteins
  }
  
  # Flag QC proteins in the data
  data %>%
    mutate(is_qc_protein = PG.Genes %in% qc_proteins_input)
}


methods <- unique(report_q$MS_method)

# Initialize an empty dataframe to store results
correlation_results <- data.frame(Method = character(), Correlation = numeric(), PValue = numeric(), stringsAsFactors = FALSE)

# Loop through each method and calculate correlation
for (meth in methods) {
  result <- correlation_analysis(Report, meth, proteins = c('lacZ'), samples = c(102, 103, 104))
  
  # Store results in the dataframe
  correlation_results <- rbind(correlation_results, data.frame(
    Method = meth,
    Correlation = result$correlation,
    PValue = result$p_value,
    adjp = result$p_adjusted
  ))
}

# View the results
print(correlation_results)
correlation_results %>%
  separate(Method, c('ms', 'grad', 'dia'), remove=F) %>%
  mutate(logpval = -log(adjp,10)) %>%
  filter(ms=='TimsTOF') %>% 
  ggplot(aes(x=Method, y=logpval)) + geom_point(aes(colour = grad)) + 
  geom_line(aes(group = grad, colour = grad))+ theme + 
  scale_x_discrete(labels=c(rep(c('1','2','3'), 3)))

Report %>%
  filter(MS_method == 'Astral_20_1', PG.Genes  == 'lacZ') %>% 
  select(R.FileName, conc_PRP, PG.Quantity) %>%
  ggplot(aes(x=conc_PRP, y=PG.Quantity)) + geom_point()  + 
  theme  +ylim(0,35000)

#### %%%%%%%%%%%%%%
#### Astral 
#### %%%%%%%%%%%%%%
qc_proteins <- c('FLNA', 'TLN1', 'MYH9', 'ACTB', 'VCL', 'ACTN1', 'TPM4', 'THBS1', 'TUBB1', 'YWHAZ', 'GSN', 'TUBA1B', 'ITGA2B', 'F13A1', 'PFN1', 'ITGB3', 'TAGLN2', 'FERMT3', 'RAP1B', 'PLEK', 'PPBP', 'GAPDH', 'MMRN1', 'MYL6', 'CFL1', 'PARVB', 'SDPR', 'TUBB4B', 'TMSB4X', 'PKM')



Report_rank_Astral <- test_report %>%
  filter(MS == 'Astral',method =='20_3', sample==103 ) %>%
  select(PG.Genes, PG.Quantity, log2) %>%
  distinct() %>%
  group_by(PG.Genes) %>%
  summarize(average_log2 = mean(log2),  cv = sd(log2) / mean(log2), count = n()) %>%
  arrange(desc(average_log2)) %>%
  mutate(rank=row_number())

#highlight 1c proteins
Report_rank_Astral %>%
  highlight_qc_proteins(qc_proteins) %>%#set what porteins you want highlifhted
  arrange(is_qc_protein) %>% #so red dots shown above black dots 
  ggplot(aes(x=rank, y=average_log2)) +
  geom_point(aes(colour=is_qc_protein), size=2, show.legend = F) +  
  labs(x = "Rank", y = expression("Protein Intensity (log"[2]*")"), color = "Platelet\nProtein") +
  scale_colour_manual(values = c('black', 'red'), labels = c('False', 'True'))+
  geom_text_repel(aes(
    label = ifelse(PG.Genes %in% c('F13A1', 'PFN1'), PG.Genes, "")), 
    size = 4, color = "red", 
    max.overlaps = Inf, 
    min.segment.length = 0, 
    box.padding = 1,  # Increased padding around text box
    point.padding = 0.5  # Padding between points and text labels
  ) + theme +
  coord_cartesian(expand=F, 
                  xlim=c(-80, max(Report_rank_Astral$rank)+80) , 
                  ylim = c(0,22))+ 
  scale_y_continuous(breaks= seq(-0, 30, by=2)) + 
  scale_x_continuous(breaks=seq(0,1500, by=300)) 

#colour by cv
Report_rank_Astral %>%
  #na.omit() %>%          # filter out NAs here
  #filter(count>=2) %>%
  ggplot(aes(x=rank, y=average_log2)) +
  geom_point(aes(colour=cv), size=2) +  
  labs(x = "Rank", y = "Intensity (log2)", color = "CV") +
  scale_colour_distiller(palette='Spectral', limits = c(0,.2), na.value='red')+
  theme

#colour by correlation
#prep df
corr_df <- r2_df %>%
  select(Genes, R2 = 'Astral_20_2')
Report_rank_Astral <- Report_rank_Astral %>%
  rename(Genes = PG.Genes)
Report_rank_Astral <- merge(Report_rank_Astral, corr_df, by="Genes", all=T)

#create plot
Report_rank_Astral %>%
  arrange(desc(R2)) %>%
  ggplot(aes(x=rank, y=average_log2)) +
  geom_point(aes(colour=R2), size=1.8, ) +
  geom_point(data = subset(Report_rank_Astral, R2 >= 0), aes(colour=R2), size=4)+
  labs(x = "Rank", y = "Intensity (log2)", color = "QC Protein") +
  scale_colour_distiller(palette='Spectral',direction="horizontal", limits = c(.5,1), na.value='grey') +
  theme +
  coord_cartesian(expand=F, 
                  xlim=c(-80, max(Report_rank_Astral$rank)+120) , 
                  ylim = c(0,15))+ 
  scale_y_continuous(breaks= seq(-1, 15, by=1)) +
  labs(x='Rank', y='Average log2')




#### %%%%%%%%%%%%%%
#### Eclipse 
#### %%%%%%%%%%%%%%
Report_rank_Eclipse <- Report %>%
  filter(MS_method == 'Eclipse_30_3', sample==104) %>%
  select(PG.Genes, PG.Quantity, log2) %>%
  distinct() %>%
  group_by(PG.Genes) %>%
  summarize(average_log2 = mean(log2),  cv = sd(log2) / mean(log2), count = n()) %>%
  arrange(desc(average_log2)) %>%
  mutate(rank=row_number())

#highlight 1c proteins
Report_rank_Eclipse %>%
  highlight_qc_proteins(qc_proteins) %>% #set what porteins you want highlifhted
  arrange(is_qc_protein) %>% #so red dots shown above black dots 
  ggplot(aes(x=rank, y=average_log2)) +
  geom_point(aes(colour=is_qc_protein), size=2) +  
  labs(x = "Rank", y = "Intensity (log2)", color = "QC Protein") +
  scale_colour_manual(values = c('black', 'red'))+
  geom_text_repel(aes(
    label = ifelse(is_qc_protein, PG.Genes, "") #ifelse(condition, True, False)
  ), size =4, color = "red", max.overlaps = Inf, min.segment.length = 0, box.padding = .5) +  
  theme

#colour by cv
Report_rank_Eclipse %>%
  #na.omit() %>%                 # filter out NAs here
  #filter(count>=2) %>%
  ggplot(aes(x=rank, y=average_log2)) +
  geom_point(aes(colour=cv), size=2) +  
  labs(x = "Rank", y = "Intensity (log2)", color = "CV") +
  scale_colour_distiller(palette='Spectral', limits = c(0,.2), na.value='red')+
  theme

#colour by corr
corr_df <- r2_df %>%
  select(Genes, R2 = 'Eclipse_30_3')
Report_rank_Eclipse <- Report_rank_Eclipse %>%
  rename(Genes = PG.Genes)
Report_rank_Eclipse <- merge(Report_rank_Eclipse, corr_df, by="Genes", all=T)

#create plot
Report_rank_Eclipse %>%
  arrange(desc(R2)) %>%
  ggplot(aes(x=rank, y=average_log2)) +
  geom_point(aes(colour=R2), size=1.8, ) +
  geom_point(data = subset(Report_rank_Eclipse, R2 >= 0), aes(colour=R2), size=4)+
  labs(x = "Rank", y = "Intensity (log2)", color = "QC Protein") +
  scale_colour_distiller(palette='Spectral',direction="horizontal", limits = c(.5,1), na.value='grey') +
  theme +
  coord_cartesian(expand=F, 
                  xlim=c(-80, max(Report_rank_Eclipse$rank)) +40, 
                  ylim = c(0,20))+ 
  scale_y_continuous(breaks= seq(-1, 20, by=1)) +
  labs(x='Rank', y='Average log2')





#### %%%%%%%%%%%%%%
#### TimsTof 
#### %%%%%%%%%%%%%%
Report_rank_TimsTof <- Report %>%
  filter(MS_method == 'TimsTOF_30_2', sample==104) %>%
  select(PG.Genes, PG.Quantity, log2) %>%
  distinct() %>%
  group_by(PG.Genes) %>%
  summarize(average_log2 = mean(log2),  cv = sd(log2) / mean(log2), count = n()) %>%
  arrange(desc(average_log2)) %>%
  mutate(rank=row_number())

#highlight 1c proteins
Report_rank_TimsTof %>%
  highlight_qc_proteins(qc_proteins) %>% #set what porteins you want highlifhted
  arrange(is_qc_protein) %>% #so red dots shown above black dots 
  ggplot(aes(x=rank, y=average_log2)) +
  geom_point(aes(colour=is_qc_protein), size=2) +  
  labs(x = "Rank", y = "Intensity (log2)", color = "QC Protein") +
  scale_colour_manual(values = c('black', 'red'))+
  geom_text_repel(aes(
    label = ifelse(is_qc_protein, PG.Genes, "") #ifelse(condition, True, False)
  ), size =4, color = "red", max.overlaps = Inf, min.segment.length = 0, box.padding = .5) +  
  theme

#colour by cv
Report_rank_TimsTof %>%
  #na.omit() %>%                 # filter out NAs here
  #filter(count>=2) %>%
  ggplot(aes(x=rank, y=average_log2)) +
  geom_point(aes(colour=cv), size=2) +  
  labs(x = "Rank", y = "Intensity (log2)", color = "CV") +
  scale_colour_distiller(palette='Spectral', limits = c(0,.2), na.value='red')+
  theme

#colour by corr
corr_df <- r2_df %>%
  select(Genes, R2 = 'TimsTOF_30_2')
Report_rank_TimsTof <- Report_rank_TimsTof %>%
  rename(Genes = PG.Genes)
Report_rank_TimsTof <- merge(Report_rank_TimsTof, corr_df, by="Genes", all=T)

#create plot
Report_rank_TimsTof %>%
  arrange(desc(R2)) %>%
  ggplot(aes(x=rank, y=average_log2)) +
  geom_point(aes(colour=R2), size=1.8, ) +
  geom_point(data = subset(Report_rank_TimsTof, R2 >= 0), aes(colour=R2), size=4)+
  labs(x = "Rank", y = "Intensity (log2)", color = "QC Protein") +
  scale_colour_distiller(palette='Spectral',direction="horizontal", limits = c(.5,1), na.value='grey') +
  theme +
  coord_cartesian(expand=F, 
                  xlim=c(-10, max(Report_rank_TimsTof$rank)) , 
                  ylim = c(0,12))+ 
  scale_y_continuous(breaks= seq(-1, 15, by=1)) +
  labs(x='Rank', y='Average log2')




# quant charts ------------------------------------------------------------

data <- summary %>%
  select(method, prot2_count, R2_1) %>%
  distinct() %>%
  separate(method, c('ms', 'grad', 'dia'), remove=F)

Report %>%
  filter(MS == 'TimsTOF',method =='15_3', sample %in% c(103)) %>%
  select(PG.Genes, PG.Quantity, log2) %>%
  distinct() %>%
  highlight_qc_proteins(qc_proteins) %>%
  filter(is_qc_protein==T) %>%
  nrow()/3

#highlight 1c proteins
Report_rank_Astral %>%
  highlight_qc_proteins(qc_proteins) %>%#set what porteins you want highlifhted
  arrange(is_qc_protein) %>% #so red dots shown above black dots 
  ggplot(aes(x=rank, y=average_log2)) +
  geom_point(aes(colour=is_qc_protein), size=2, show.legend = F) +  
  labs(x = "Rank", y = expression("Protein Intensity (log"[2]*")"), color = "Platelet\nProtein") +
  scale_colour_manual(values = c('black', 'red'), labels = c('False', 'True'))+
  geom_text_repel(aes(
    label = ifelse(PG.Genes %in% c('F13A1', 'PFN1'), PG.Genes, "")), 
    size = 4, color = "red", 
    max.overlaps = Inf, 
    min.segment.length = 0, 
    box.padding = 1,  # Increased padding around text box
    point.padding = 0.5  # Padding between points and text labels
  ) + theme +
  coord_cartesian(expand=F, 
                  xlim=c(-80, max(Report_rank_Astral$rank)+80) , 
                  ylim = c(0,22))+ 
  scale_y_continuous(breaks= seq(-0, 30, by=2)) + 
  scale_x_continuous(breaks=seq(0,1500, by=300)) 

  
#### %%%%%%%%%%%%%%
#### Astral 
#### %%%%%%%%%%%%%%
data %>%
  filter(ms == 'Astral') %>%
  ggplot(aes(x=method)) + 
  geom_bar(aes(fill=grad, y=prot2_count), stat = 'identity', colour = 'black', alpha = 0.5) + 
    scale_fill_manual(
    values = c('#D95F02', "#E7298A", "#66A61E" ),
    labels = c("15 mins", "22 mins", "30 mins"),
    name = "Gradient (mins)"
  )+ 
  scale_alpha_manual(
    values  = c(0.1, 0.25, 0.5),
    labels = c('Total', "CV < 20%", "CV < 5%"),
    name = "Protein ID's"
  )+
  scale_x_discrete(
    labels = c('2 Da, 7 ms', '4 Da, 10 ms', 'Facility',
               '2 Da, 7 ms', '4 Da, 10 ms', 'Facility',
               '2 Da, 7 ms', '4 Da, 10 ms', 'Facility'),
    name = NULL
  ) +
  theme +
  ylim(0, 19) +
  #add R2 sec axis
  geom_line(aes(y=R2_1*16, group=1), size=1)+
  geom_point(aes(y=R2_1*16), colour = 'red', size =3) +
  scale_y_continuous(
    name = 'Proteins (Count)',
    sec.axis = sec_axis(~./16, name = 'R^2 (avg)', )
  ) + 
  theme(axis.text.y.right = element_text(colour='red'),
        axis.title.y.right = element_text(colour='red'),
        axis.line.y.right = element_line(colour='red')) +
  coord_cartesian(expand=F, xlim = c(0.3, 9.7))
  



#### %%%%%%%%%%%%%%
#### Eclipse 
#### %%%%%%%%%%%%%%

data %>%
  filter(ms == 'Eclipse') %>%
  ggplot(aes(x=method)) + 
  geom_bar(aes(fill=grad, y=prot2_count), stat = 'identity', colour = 'black', alpha = 0.5) + 
  scale_fill_manual(
    values = c('#D95F02', "#66A61E", "#E6AB02"),
    labels = c("15 mins", "30 mins", "45 mins"),
    name = "Gradient (mins)"
  )+
  scale_x_discrete(
    labels = c('10 Da, 31 ms', '20 Da, 31 ms', 'Facility',
               '10 Da, 31 ms', '20 Da, 31 ms', 'Facility',
               '10 Da, 31 ms', '20 Da, 31 ms', 'Facility'),
    name = NULL
  ) +
  theme +
  ylim(0, 19) +
  #add R2 sec axis
  geom_line(aes(y=R2_1*14, group=1), size=1)+
  geom_point(aes(y=R2_1*14), colour = 'red', size =3) +
  scale_y_continuous(
    name = 'Proteins (Count)',
    sec.axis = sec_axis(~./14, name = 'R^2 (avg)', )
  ) + 
  theme(axis.text.y.right = element_text(colour='red'),
        axis.title.y.right = element_text(colour='red'),
        axis.line.y.right = element_line(colour='red')) +
  coord_cartesian(expand=F, xlim = c(0.3, 9.7))














#### %%%%%%%%%%%%%%
#### TimsTOF 
#### %%%%%%%%%%%%%%

data %>%
  filter(ms == 'TimsTOF') %>%
  ggplot(aes(x=method)) + 
  geom_bar(aes(fill=grad, y=prot2_count), stat = 'identity', colour = 'black', alpha = 0.5) + 
  scale_fill_manual(
    values = c('#D95F02', "#66A61E", "#E6AB02"),
    labels = c("15 mins", "30 mins", "45 mins"),
    name = "Gradient (mins)"
  )+
  scale_x_discrete(
    labels = c('thinPASEF', 'py-diAID', 'Facility',
               'thinPASEF', 'py-diAID', 'Facility',
               'thinPASEF', 'py-diAID', 'Facility'),
    name = NULL
  ) +
  theme +
  ylim(0, 19) +
  #add R2 sec axis
  geom_line(aes(y=R2_1*14, group=1), size=1)+
  geom_point(aes(y=R2_1*14), colour = 'red', size =3) +
  scale_y_continuous(
    name = 'Proteins (Count)',
    sec.axis = sec_axis(~./14, name = 'R^2 (avg)', )
  ) + 
  theme(axis.text.y.right = element_text(colour='red'),
        axis.title.y.right = element_text(colour='red'),
        axis.line.y.right = element_line(colour='red')) +
  coord_cartesian(expand=F, xlim = c(0.3, 9.7))








# window images  ----------------------------------------------------------


plot_pasef_cycles <- function(file_path) {
  # Read the data
  data <- read.csv(file_path, header = TRUE)
  
  # Convert columns to numeric and factor
  data$Start.IM..1.K0. <- as.numeric(data$Start.IM..1.K0.)
  data$End.IM..1.K0. <- as.numeric(data$End.IM..1.K0.)
  data$Start.Mass..m.z. <- as.numeric(data$Start.Mass..m.z.)
  data$End.Mass..m.z. <- as.numeric(data$End.Mass..m.z.)
  data$Cycle.Id <- as.factor(data$Cycle.Id)
  data <- data[-c(1), -c(7)]
  #data <- data[c(4),]
  
  # Determine number of cycles for the rainbow palette
  n_cycles <- length(unique(data$Cycle.Id))
  
  # Plot
  data %>% 
    ggplot() + 
    geom_rect(aes(
      xmin = Start.Mass..m.z., xmax = End.Mass..m.z.,
      ymin = Start.IM..1.K0., ymax = End.IM..1.K0.,
      fill = Cycle.Id
    ), color = "black", alpha = 0.5, show.legend = F) +
    labs(
      
      x = "m/z",
      y = expression("Ion Mobility (1/K"[0]*")"),
      fill = "PASEF Cycle"
    ) + 
    scale_y_continuous(breaks = seq(0,2, by=.1)) +
    scale_x_continuous(breaks = seq(0,1400, by =200)) + 
    scale_fill_manual(values = rainbow(n_cycles)) +
    coord_cartesian(ylim = c(0.6, 1.45), xlim = c(400,1210)) +
    theme+ 
    guides(fill = guide_legend(ncol = 3))
  
}

files <- c(
  'windows/current.txt',
  'windows/neat_diaPASEF_method.txt',
  'windows/thinpasef_final.txt'
)
for (file in files) {
  print(plot_pasef_cycles(file))
}

plot_mz_windows <- function(mz_list) {
  
  # Determine the number of original MS2 windows
  num_windows <- length(mz_list)
  
  # Insert MS1 windows at the start and at the midpoint
  midpoint <- ceiling(num_windows / 2)
  mz_list <- append(mz_list, "380-980", after = midpoint)  # Add at midpoint
  mz_list <- append(mz_list, "380-980", after = 0)  # Add at the start
  
  # Convert mz_list to a data frame
  mz_df <- data.frame(m.z.range = mz_list)
  
  # Expand the m/z range into start and end points
  mz_df <- mz_df %>%
    separate(m.z.range, into = c("mz_start", "mz_end"), sep = "-", convert = TRUE)
  
  # Recalculate the number of windows and set the time increment
  num_windows <- nrow(mz_df)
  time_increment <- 1 / num_windows
  
  # Add time relative to the full cycle for all rectangles
  mz_df <- mz_df %>%
    mutate(
      time = seq(from = 0, to = 1 - time_increment, by = time_increment),
      type = ifelse(mz_start == 380 & mz_end == 980, "MS1", "MS2")  # Label MS1 and MS2
    )
  
  # Combine the data and calculate the ymin and ymax for each rectangle
  all_rects <- mz_df %>%
    mutate(
      ymin = time,
      ymax = time + time_increment
    )
  
  # Plotting the m/z ranges as rectangles with the y-axis flipped
  p <- ggplot(all_rects, aes(xmin = mz_start, xmax = mz_end, ymin = ymin, ymax = ymax, fill = type)) +
    geom_rect(alpha = 1, colour = 'black', size = 0.1) +
    scale_x_continuous(breaks = seq(300,1500, by = 300), limits = c(348, 1651)) +
    scale_y_reverse(name = "Time Relative to Full Cycle", breaks = seq(0, 1.01, by = 0.1)) +  # Reverse the y-axis
    scale_fill_manual(values = c("MS1" = "#00BCD4", "MS2" = "#8A8D91"), name = "Scan Level") +
    theme_minimal() +
    labs(x = expression(italic('m/z'))) +
    coord_cartesian(expand = FALSE) +
    theme
  
  return(p)
}




files = c('windows/Mass List Table-Dyanmic - Copy.csv')
for (file in files) {
  tt <- read.delim('windows/eclipse final.csv', sep=',')
  vec <- tt$m.z.range #Calculated.m.z.Window for elcipse,. m.z.range for astral
  plot = plot_mz_windows(vec)
  print(plot) 
}


# upset -------------------------------------------------------------------

library(ggupset)
library(ggplot2)
library(tidyr)
library(dplyr)

# Extract distinct protein genes for each method
distinct_genes_per_method <- report_q %>%
  filter(PEP.IsProteotypic=='True') %>%
  select(MS_method, PG.Genes) %>%
  distinct() %>%
  group_by(MS_method) %>%
  summarize(Genes = list(PG.Genes))  # Create a list of genes per method

# Take the first 10 methods to plot
distinct_genes_top10 <- distinct_genes_per_method %>%
  slice(1:10)  # Select the first 10 combinations

# Create a data frame with a list-column containing methods for each gene
genes_combined <- distinct_genes_per_method %>%
  unnest(Genes) %>%
  group_by(Genes) %>%
  summarize(Methods = list(MS_method))  # Group genes and combine methods

# Generate the upset plot
ggplot(genes_combined, aes(x = Methods)) +
  geom_bar() +
  scale_x_upset(n_intersections = 10) +  # Display the first 10 combinations
  labs(title = "Upset Plot of Distinct Protein Genes by Method",
       x = "Method Combinations",
       y = "Number of Genes") +
  theme




