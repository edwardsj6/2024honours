
# formal analysis  --------------------------------------------------------
install.packages("wesanderson")
install.packages("RColorBrewer")
install.packages("tidyverse")
install.packages("ggrepel")
library("ggplot2")
library("RColorBrewer")

library("ggridges")

library('tidyverse')
library(ggrepel)
display.brewer.all(colorblindFriendly = T)
#take dark 2, 15=2 22=3 30=4 45=5

theme <- theme_classic() + theme(
  text = element_text((family = "Arial"), size = 10, colour ="black"),  # Set all text to Arial, size 12
  plot.title = element_text((family = "Arial"), size = 14),  # Adjust title size (optional)
  #strip.text.x = element_text((family = "Arial"), size = 12),  # Adjust strip text size (optional)
  axis.text = element_text((family = "Arial"), size = 10, colour ="black"),
  legend.text = element_text((family = "Arial"), size = 10, colour ="black")
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



#### all IDs --------
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




#### IDs10 --------

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


#### IDs5 --------

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


#### correlation ------------
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

#### ids of qc proteins -----
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

summary_d <- summary
##### %%%%%%%%%%%%%%
##### ASTRAL Plasma
##### %%%%%%%%%%%%%%
#data for precs
line <- summary_d %>%
  select(method, avgprec) %>%
  distinct() %>% #since avg values
  separate(method, c('ms', 'grad', 'dia'), sep="_",remove=F) %>%
  filter(ms == "NA")
coeff = 10.5

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
  filter(ms == "NA") %>%
  
  
  ggplot(aes(x=method, y=mean, fill=grad)) + 
  geom_bar(aes(alpha=Metric), colour='black', size=.5, stat='identity', position = 'dodge') + 
  geom_point(aes(x=method, y=Values, fill=grad, group=Metric), 
             shape=21, size=1.5, position=position_jitterdodge(
               jitter.width=.06,dodge.width=.9), show.legend = F) +
  scale_fill_manual(
    values = c('#D95F02', "#E7298A", "#66A61E" ),
    labels = c("15 mins", "22 mins", "30 mins"),
    name = "Gradient (mins)"
  )+ 
  
  
  
  scale_alpha_manual(
    values  = c(0.1, 0.25, 0.5),
    labels = c('Total', 'Filtered', 'CV < 20%'),
    name = "Protein ID's"
  )+
  scale_x_discrete(
    labels = c('2 Th', '4 Th', 'Facility',
               '2 Th', '4 Th', 'Facility',
               '2 Th', '4 Th', 'Facility'),
    name = NULL
  ) +
  labs( title = "Astral")+
  
  coord_cartesian(expand = F,xlim=c(0.3, 9.7), ylim=c(0,1800)) +
  theme + 
  #add precursros sec axis
  geom_line(data=line, aes(x=method, y=avgprec/coeff, group =1), linewidth =.5)+ 
  geom_point(data=line, aes(x=method, y=avgprec/coeff), size =2.5, colour ='red', show.legend = F) +
  
  scale_y_continuous(
    name = "Number of Proteins",
    sec.axis = sec_axis(~.*coeff, name = 'Number of Precursors', breaks = seq(0,20000, by=1500)),
    breaks = seq(0, 2000, by = 100)
  ) + 
  labs()+
  theme(axis.text.y.right = element_text(colour='red'),
        axis.title.y.right = element_text(colour='red'),
        axis.line.y.right = element_line(colour='red'),
        axis.ticks.y.right = element_line(colour='red'))


# Venn diagram ------------------------------------------------------------
install.packages('VennDiagram')
library(VennDiagram)

# Base data
data_base <- Report_filt[Report_filt$MS_method == 'Astral_30_1' & Report_filt$sample == 104, ]

# No filtering
genes_base <- data_base$PG.Genes

# Filter for PEP <= 0.01
data_pep <- data_base %>% filter(PG.PEP..Run.Wise. > 0.01) %>% distinct()
genes_pep <- data_pep$PG.Genes

# Filter for PG QValue <= 0.01
data_pg_qval <- data_base %>% filter(PG.QValue..Run.Wise.> 0.01)%>% distinct()
genes_pg_qval <- data_pg_qval$PG.Genes

# Filter for single hits == FALSE
data_single_hits <- data_base %>% filter(PG.IsSingleHit == 'True')%>% distinct()
genes_single_hits <- data_single_hits$PG.Genes


gene_list <- list(

  "PEP > 0.01" = genes_pep,
  "PG QValue > 0.01" = genes_pg_qval,
  " Is Single Hit" = genes_single_hits
)

venn_plot <- ggVennDiagram(gene_list,force_upset = TRUE) 
                           label_alpha = 0.8) + 
  theme_void() +  # Use a simple theme
  theme(legend.position = "bottom") 
print(venn_plot)

# rank abundance plots ----------------------------------------------------

#highlight qc_proteins
highlight_qc_proteins <- function(data, qc_proteins_input) {
  data %>%
    mutate(is_qc_protein = PG.Genes %in% qc_proteins_input)  # Flag QC proteins
}



#### %%%%%%%%%%%%%%
#### Astral 
#### %%%%%%%%%%%%%%
Report_rank_Astral <- Report %>%
  filter(R.FileName == 'P4776_104_20_2_R2') %>%
  select(PG.Genes, PG.Quantity, log2) %>%
  distinct() %>%
  group_by(PG.Genes) %>%
  summarize(average_log2 = mean(log2),  cv = sd(log2) / mean(log2), count = n()) %>%
  arrange(desc(average_log2)) %>%
  mutate(rank=row_number())

#highlight 1c proteins
Report_rank_Astral %>%
  highlight_qc_proteins(qc_proteins) %>% #set what porteins you want highlifhted
  arrange(is_qc_protein) %>% #so red dots shown above black dots 
  ggplot(aes(x=rank, y=average_log2)) +
  geom_point(aes(colour=is_qc_protein), size=2) +  
  labs(x = "Rank", y = "Intensity (log2)", color = "QC Protein") +
  scale_colour_manual(values = c('black', 'red'))+
  geom_text_repel(aes(
    label = ifelse(is_qc_protein, PG.Genes, "") #ifelse(condition, True, False)
    ), size =4, color = "red", max.overlaps = Inf, min.segment.length = 0, box.padding = .5) +  
  theme +
  coord_cartesian(expand=F, 
                  xlim=c(-80, max(Report_rank_Astral$rank)) +40, 
                  ylim = c(0,15))+ 
  scale_y_continuous(breaks= seq(-1, 15, by=1)) 


#colour by cv
Report_rank_Astral %>%
  #na.omit() %>%          # filter out NAs here
  #filter(count>=2) %>%
  ggplot(aes(x=rank, y=average_log2)) +
  geom_point(aes(colour=cv), size=2) +  
  labs(x = "Rank", y = "Intensity (log2)", color = "CV") +
  scale_colour_distiller(palette='Spectral', limits = c(0,.2), na.value='red')+
  theme





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





# magnet  -----------------------------------------------------------------
report_files <- list(
  "P4854/15mins/Report.tsv",
  "P4854/22mins/Report.tsv",
  "P4854/30mins/Report.tsv"
)


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
  
  
  
  
  #conc map 
  if (map_concs == T) {
    print('mapping concs')
    conc_map <- data.frame(sample = c('101', '102','103', '104'), conc_PRP = c(1, .5, .1,0))
    Report <- merge(Report, conc_map, by = 'sample', all.x = TRUE)
  }
  
  
  
  return(Report) 
}

generate_summary_report <- function(data_in, filter_levels) {
  # Filter proteotypic data
  Report_filt <- data_in %>%
    filter(PEP.IsProteotypic == 'True')
  
  # Count total IDs for a method, filtered for proteotypic
  Report_filt$MS_method <- paste0(Report_filt$MS, '_', Report_filt$method)
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
    methods_vec <- c(methods_vec, rep(meth, length(dIDs$Replicate)))
    
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
    prec_vec <- c(prec_vec, davgprec$n_unique_genes)
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


filts <- list(
  list(pep_threshold = .05, pg_qval = 0.05), 
  list(pep_threshold = 0.01, pg_qval = 0.01)
)

rep = generate_report(report_files)

summary <- generate_summary_report(rep, filts)

plot_protein_precursor_errorbars <- function(summary_d, ms = "Astral", coeff = 9.3, fill_colors, fill_labels, x_labels, 
                                             show_protein_error = TRUE, show_precursor_error = TRUE) {
  summary_d <- plasma_summary
  # Generate summary stats for protein data
  summary_stats <- summary_d %>%
    pivot_longer(cols = c(IDstot, IDsfilt, IDs20, avgprec),
                 names_to = 'Metric',
                 values_to = 'Values') %>%
    group_by(method, Metric) %>%
    summarise(mean = mean(Values),
              sd = sd(Values),
              se = sd(Values) / sqrt(n()),  # Calculate standard error
              .groups = 'drop') %>%
    separate(method, c('ms', 'grad', 'dia'), remove = F) %>%
    filter(ms == !!ms) %>%
    mutate(Metric = factor(Metric, c('IDstot', 'IDsfilt', 'IDs20', 'avgprec')))
  
  # Generate data for plotting precursors
  line <- summary_d %>%
    select(method, avgprec) %>%
    group_by(method) %>%
    summarise(se_prec = sd(avgprec) / sqrt(3), 
              avgprec = mean(avgprec),
               # Calculate standard error
              .groups = 'drop') %>%
    separate(method, c('ms', 'grad', 'dia'), sep = "_", remove = F) 
  
  # Main plot
  t <- summary_d %>%
    select(method, IDstot, IDsfilt, IDs20) %>%
    pivot_longer(cols = c(IDstot, IDsfilt, IDs20),
                 names_to = 'Metric',
                 values_to = 'Values') %>%
    mutate(Metric = factor(Metric, c('IDstot', 'IDsfilt', 'IDs20'))) %>%
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
    
    coord_cartesian(expand = FALSE, xlim = c(0.3, 9.7), ylim = c(0, 1900)) +
    
    theme +
    
    geom_line(data = line, aes(x = method, y = avgprec / coeff, group = grad), linewidth = 0.5) + 
    geom_point(data = line, aes(x = method, y = avgprec / coeff), size = 2.5, colour = 'red', show.legend = FALSE) +
    
    scale_y_continuous(name = "Number of Proteins",
                       sec.axis = sec_axis(~.*coeff, name = 'Number of Precursors', breaks = seq(0, 29000, by = 1500)),
                       breaks = seq(0, 2000, by = 100)) + 
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
plot_protein_precursor_errorbars(plasma_summary, ms = "NA", coeff = 9.6, 
                                 fill_colors = c('#C9AA81', "#808000", "#8A360F"), 
                                 fill_labels = c("15 mins", "22 mins", "30 mins"), 
                                 x_labels = c('A1', 'A2', 'A3', 'A1', 'A2', 'A3', 'A1', 'A2', 'A3'),
                                 show_protein_error = F, 
                                 show_precursor_error = TRUE)








# density plots -----------------------------------------------------------

neat <- generate_report('neat/dia1-20.tsv')
magnet_rep <- Report%>%
  filter(ID == '_104_20_1') %>%
  filter(PEP.IsProteotypic == 'True', PG.IsSingleHit == 'False') %>%
  select(PG.Genes, EG.PrecursorId, PG.Quantity, R.FileName) %>%
  distinct() %>%
  group_by(PG.Genes) %>% 
  summarise(avg_log10 = mean(log(PG.Quantity, 2))) %>%
  ungroup()

neat_rep <- neat %>%
  filter(PEP.IsProteotypic == 'True', PG.IsSingleHit == 'False') %>%
  select(PG.Genes, EG.PrecursorId, PG.Quantity, R.FileName) %>%
  distinct() %>%
  group_by(PG.Genes) %>% 
  summarise(avg_log10 = mean(log(PG.Quantity, 2))) %>%
  ungroup()
 

# Combine both datasets into one for easy plotting
combined_rep <- bind_rows(
  mutate(magnet_rep, group = "Magnet"),
  mutate(neat_rep, group = "Neat")
)

# Plot the density
p <- combined_rep %>%
  mutate(group = factor(group, c('Neat', 'Magnet'))) %>%
  ggplot(aes(x = avg_log10, fill = group)) +
  geom_density(alpha = 0.5) +  # Add transparency for overlapping densities
  scale_fill_manual(values = c( '#008080','#C71585')) +  # Color palette
  theme_minimal() +  # Clean theme
  labs(x = expression('Protein Intensity (log'[2]*')'),
       y = "Density",
       fill = "Sample Preparation") +
  
  theme + 
  coord_cartesian(expand=F)

# Print the plot
print(p)


library(eulerr)

magnet_prots <- Report%>%
  filter(ID == '_104_20_1') %>%
  filter(PEP.IsProteotypic == 'True', PG.IsSingleHit == 'False') %>%
  select(PG.Genes) %>% 
  distinct() %>%
  pull(PG.Genes) 
neat_rep <- neat %>%
  filter(PEP.IsProteotypic == 'True', PG.IsSingleHit == 'False') %>%
  select(PG.Genes) %>% 
  distinct() %>%
  pull(PG.Genes) 

# Filter the data for conc_PRP == 0 and three methods (Astral, timsTOF, Eclipse)


# Create the sets for each method
neat_proteins <- neat_rep
magnet_proteins <- magnet_prots

# Identify unique and overlapping proteins
only_neat <- setdiff(neat_proteins, magnet_proteins)
only_magnet <- setdiff(magnet_proteins, neat_proteins)
neat_and_magnet <- intersect(neat_proteins, magnet_proteins)

# Create a named numeric vector for eulerr
protein_counts <- c(
  "Neat" = length(only_neat),
  "Magnet" = length(only_magnet),
  "Neat&Magnet" = length(neat_and_magnet)
)

# Generate the Euler diagram
fit <- euler(protein_counts)

# Plot the Euler diagram
plot(fit, quantities = TRUE, fill = c( '#90D1D1','#DE9DC6'),
     labels = list(font = 2), main = NULL)

write.csv(only_neat, 'neatonly.csv')

unique_genes = only_neat
neat %>%
  filter(PEP.IsProteotypic == 'True', PG.IsSingleHit == 'False') %>%
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
  scale_colour_manual(values = c('black','#008080'), labels = c('False', 'True'))+
    theme  + 
  scale_y_continuous(breaks=seq(0,30,by=3)) + 
  scale_x_continuous(breaks=seq(0,2000, by=250))+
  labs(x = "Rank", 
       y = expression("Protein Intensity (log"[2]*")"), color = "Unique\nIdentification")




