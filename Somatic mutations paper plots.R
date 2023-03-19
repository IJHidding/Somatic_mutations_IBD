###### BOXplot or dots or whatever to show effect of filtering on the number of somatic mutations in each samples
library(ggplot2)
library(dplyr)
library(hrbrthemes)
library(viridis)

# create a dataset
filter_1_i <- read.table('/Users/iwan/Research/Somatic_mutations/plot_files/inflamed_biopsies_mutect_counts.txt')
filter_1_n <- read.table('/Users/iwan/Research/Somatic_mutations/plot_files/non_inflamed_biopsies_mutect_counts.txt')

filter_2_i <- read.table('/Users/iwan/Research/Somatic_mutations/plot_files/inflamed_biopsies_post_rare_filter_counts.txt')
filter_2_n <- read.table('/Users/iwan/Research/Somatic_mutations/plot_files/non_inflamed_biopsies_post_rare_filter_counts.txt')

filter_3_i <- read.table('/Users/iwan/Research/Somatic_mutations/plot_files/inflamed_biopsies_post_func_filter_counts.txt')
filter_3_n <- read.table('/Users/iwan/Research/Somatic_mutations/plot_files/non_inflamed_biopsies_post_func_filter_counts.txt')

#for (missing_rows in 1:(nrow(filter_1) - nrow(filter_2))){
#  filter_2[nrow(filter_2)+1,] <- NA
#}

#for (missing_rows in 1:(nrow(filter_1) - nrow(filter_3))){
#  filter_3[nrow(filter_3)+1,] <- NA
#}

data <- data.frame(
  name=c( rep("mutect calls non-inflamed",nrow(filter_1_n)), rep("mutect calls inflamed",nrow(filter_1_i)),
          rep("rarefaction filter non-inflamed",nrow(filter_2_n)), rep("rarefaction filter inflamed",nrow(filter_2_i)),
          rep("functional_filter non-inflamed",nrow(filter_3_n)), rep("functional_filter inflamed",nrow(filter_3_i))),
  value=c( filter_1_n$V1, filter_1_i$V1, filter_2_n$V1, filter_2_i$V1, filter_3_i$V1, filter_3_n$V1)
)

# sample size
sample_size = data %>% group_by(name) %>% dplyr::summarize(num=n())

# Plot
data %>%
  ggplot( aes(x=name, y=value, fill=name)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("A boxplot with jitter") +
  xlab("")








# beatufiul circular plot to show the base changes and numbers




# library
library(tidyverse)

# Create dataset#, rep('Insertion', 4), rep('Deletion', 4)) #, inser:(), del:())
# groups :T, G, A, C, insert, deletion
# individuals: T, G, A, C, 
data <- data.frame(
  individual=c('T', 'A', 'C', 'G'),
  group=c(rep('T, n=30882', 4), rep('A, n=33551', 4),rep('C, n=26164', 4), rep('G, n=27740', 4), rep('Del, n=18015', 4), rep('Ins, n=87799', 4)), 
  value=c(0, 0.15, 0.74, 0.11, 0.14, 0, 0.10, 0.76, 0.62, 0.20, 0, 0.18, 0.18, 0.63, 0.19, 0, 0.32, 0.31,0.19, 0.18, 0.40, 0.39, 0.12, 0.09)
)
data <- subset(data, value != 0)
rownames(data) <- NULL 

data <- data %>%
  mutate(group = factor(group, levels=c("T, n=30882", "A, n=33551", "C, n=26164", "G, n=27740", "Del, n=18015", "Ins, n=87799")))
  
#mutate(name = factor(name, levels=c("north", "north-east", "east", "south-east", "south", "south-west", "west", "north-west"))) %>%
  

#data <- data.frame(
#  individual=paste( "Mister ", seq(1,60), sep=""),
#  group=c( rep('A', 10), rep('B', 30), rep('C', 14), rep('D', 6)) ,
#  value=sample( seq(10,100), 60, replace=T)
#)

# Set a number of 'empty bar' to add at the end of each groupnlevels(data$group) levels(data$group
empty_bar <- 1
to_add <- data.frame( matrix(NA, empty_bar*6, ncol(data)) )
colnames(to_add) <- colnames(data)
to_add$group <- rep(c("A, n=33551", "T, n=30882" ,"C, n=26164", "G, n=27740", "Del, n=18015", "Ins, n=87799"), each=empty_bar)
data <- rbind(data, to_add)
data <- data %>% arrange(group)
data$id <- seq(1, nrow(data))

# Get the name and the y position of each label
label_data <- data
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

# prepare a data frame for base lines
base_data <- data %>% 
  group_by(group) %>% 
  dplyr::summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

# prepare a data frame for grid (scales)
grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]

# Make the plot
p <- ggplot(data, aes(x=as.factor(group), y=value, fill=individual)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  
  geom_bar(aes(x=as.factor(group), y=value, fill=individual),  position = 'stack', stat="identity", alpha=0.5) +
  

  # Add text showing the value of each 100/75/50/25 lines
  annotate("text", x = rep(max(data$id),4)+0.5, y = c(0.10, 0.25, 0.50, 0.75), label = c("0.10", "0.25", "0.50", "0.75") , color="grey", size=3 , angle=0, fontface="bold", hjust=1) +

  
  geom_bar(aes(x=as.factor(group), y=value, fill=individual), position = 'stack', stat="identity", alpha=0.5) +
  ylim(-0.5,1) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm")
  ) +
  #coord_polar() +
  #geom_text(data=label_data, aes(x=id, y=value+0.1, label=individual, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=2.5, angle= label_data$angle, inherit.aes = FALSE ) +
  #ylim(-10000, 40000)
  # Add base line information
  geom_segment(data=base_data, aes(x = start, y = -0.04, xend = end, yend = -0.04), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE ) +
  geom_text(data=base_data, aes(x = title, y = 0, label=group), hjust=c(0.2,-0.1,0.25,1, 1.4, 1.25), vjust=c(-18.9, -13, -4.8,-2.6,-10 ,-18 ),colour = "black", alpha=0.8, size=3, fontface="bold", inherit.aes = FALSE) +
  # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
  geom_segment(data=grid_data, aes(x = 0.4, y = 0.75, xend = 25, yend = 0.75), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = 0.4, y = 0.50, xend = 25, yend = 0.50), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = 0.4, y = 0.25, xend = 25, yend = 0.25), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = 0.4, y = 0.10, xend = 25, yend = 0.10), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  
  ggtitle('Somatic mutations single base change overview')

p

specie <- c(rep("sorgho" , 3) , rep("poacee" , 3) , rep("banana" , 3) , rep("triticum" , 3) )
condition <- rep(c("normal" , "stress" , "Nitrogen") , 4)
value <- abs(rnorm(12 , 0 , 15))
test_data <- data.frame(specie,condition,value)
ggplot(test_data, aes(fill=condition, y=value, x=specie)) + 
  geom_bar(position="stack", stat="identity")


p <- ggplot(data, aes(x=group, y=value, fill=individual)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  geom_bar(position = 'stack', stat="identity", alpha=0.5)
p
install.packages('ggstatsplot')
library(ggstatsplot)


library(tidyverse)
library(reshape2)

#data("penguins", package = "palmerpenguins")

#penguins <- drop_na(penguins)


df <- as.data.frame(read_csv('/Users/iwan/Research/Somatic_mutations/plot_dataframes/top_18_genes.csv'))
#rownames(df) <- df$...1
#df$...1 <- NULL

gene_length <- c(1494.0, 1098.0, 1056.0,  3966.0,  3785.0,  2583.0, 1527.0, 1047.0, 645.0, 327.0, 1902.0, 1677.0, 276.0, 1782.0, 1032.0, 735.0, 651.0)

df_corrected <- df[, -1] / gene_length 
df_corrected$...1 <- df$...1

#df$gene_length <- gene_length
#df[1:17,:]
melted_df <- melt(df, id='...1')
melted_df <- melted_df %>%
  mutate(...1 = factor(...1, levels=c("XIAP", "ZDHHC20", "TMEM236", "SPAG9", "AC073869.1", "MPHOSPH8"
                                      , "SLC19A3", "HFE", "MREG", "FKBP1C", "CLPX"
                                      , "GLUD2", "UQCRHL", "ARSD", "APOL6", "METTL7A"
                                      , "RAB2B")))


library(pals)
cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
plt <- ggbetweenstats(
  data = melted_df,
  x = ...1,
  y = value,
  pairwise.comparisons = FALSE,
  package = "pals",
  palette = "stepped",
  point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.3), alpha =
                      0.7, size = 2, stroke = 0),
  centrality.plotting = FALSE,
  results.subtitle = FALSE
  #centrality.point.args = list(size = 2, color = "darkred"),
  #centrality.label.args = list(size = 0.0001, nudge_x = 0.4, segment.linetype = 4,
                            #   min.segment.length = 0)
)

plt <- plt + 
  # Add labels and title
  labs(
    x = "Genes",
    y = "Number of somatic mutations per biopsy",
    title = "Number of somatic mutations per gene in genes with an average > 1"
  )


 show(plt)

 ## correct dataframe mutations for gene length 
 
 melted_df_corrected <- melt(df_corrected, id='...1')
 melted_df_corrected <- melted_df_corrected %>%
   mutate(...1 = factor(...1, levels=c("XIAP", "ZDHHC20", "TMEM236", "SPAG9", "AC073869.1", "MPHOSPH8"
                                       , "SLC19A3", "HFE", "MREG", "FKBP1C", "CLPX"
                                       , "GLUD2", "UQCRHL", "ARSD", "APOL6", "METTL7A"
                                       , "RAB2B")))
 
 
 library(pals)
 cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
           "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
 plt <- ggbetweenstats(
   data = melted_df_corrected,
   x = ...1,
   y = value,
   pairwise.comparisons = FALSE,
   package = "pals",
   palette = "stepped",
   point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.3), alpha =
                       0.7, size = 2, stroke = 0),
   centrality.plotting = FALSE,
   results.subtitle = FALSE
   #centrality.point.args = list(size = 2, color = "darkred"),
   #centrality.label.args = list(size = 0.0001, nudge_x = 0.4, segment.linetype = 4,
   #   min.segment.length = 0)
 )
 
 plt <- plt + 
   # Add labels and title
   labs(
     x = "Genes",
     y = "Number of somatic mutations per biopsy per base in the gene",
     title = "Number of somatic mutations per gene in genes with an average > 1 corrected for gene length"
   )
 
 
 show(plt) 
 
 
 
# correct mutation number by gene length: 
 
 





# create a dataset
filter_1_i <- read.table('/Users/iwan/Research/Somatic_mutations/plot_files/inflamed_biopsies_mutect_counts.txt')
filter_1_n <- read.table('/Users/iwan/Research/Somatic_mutations/plot_files/non_inflamed_biopsies_mutect_counts.txt')

filter_2_i <- read.table('/Users/iwan/Research/Somatic_mutations/plot_files/inflamed_biopsies_post_rare_filter_counts.txt')
filter_2_n <- read.table('/Users/iwan/Research/Somatic_mutations/plot_files/non_inflamed_biopsies_post_rare_filter_counts.txt')

filter_3_i <- read.table('/Users/iwan/Research/Somatic_mutations/plot_files/inflamed_biopsies_post_func_filter_counts.txt')
filter_3_n <- read.table('/Users/iwan/Research/Somatic_mutations/plot_files/non_inflamed_biopsies_post_func_filter_counts.txt')

#for (missing_rows in 1:(nrow(filter_1) - nrow(filter_2))){
#  filter_2[nrow(filter_2)+1,] <- NA
#}

#for (missing_rows in 1:(nrow(filter_1) - nrow(filter_3))){
#  filter_3[nrow(filter_3)+1,] <- NA
#}

# mutect files still have the header so we remove the number of header rows
filter_1_i$V1 <- filter_1_i$V1 - 3413
filter_1_n$V1 <- filter_1_n$V1 - 3413

data <- data.frame(
  name=c( rep("mutect calls non-inflamed",nrow(filter_1_n)), rep("mutect calls inflamed",nrow(filter_1_i)),
          rep("rarefaction filter non-inflamed",nrow(filter_2_n)), rep("rarefaction filter inflamed",nrow(filter_2_i)),
          rep("functional_filter non-inflamed",nrow(filter_3_n)), rep("functional_filter inflamed",nrow(filter_3_i))),
  value=c( filter_1_n$V1, filter_1_i$V1, filter_2_n$V1, filter_2_i$V1, filter_3_i$V1, filter_3_n$V1)
)

data <- data %>%
  mutate(name = factor(name, levels=c("mutect calls inflamed", "mutect calls non-inflamed", "rarefaction filter inflamed", "rarefaction filter non-inflamed", "functional_filter inflamed", "functional_filter non-inflamed")))

plt <- ggbetweenstats(
  data = data,
  x = name,
  y = value,
  pairwise.comparisons = TRUE,
  palette = "Dark2",
  point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.3), alpha =
                      0.7, size = 2, stroke = 0),
  centrality.plotting = TRUE,
  results.subtitle = FALSE,
  pairwise.display = 'ns'
  #centrality.point.args = list(size = 2, color = "darkred"),
  #centrality.label.args = list(size = 0.0001, nudge_x = 0.4, segment.linetype = 4,
  #   min.segment.length = 0)
)

plt <- plt + 
  # Add labels and title
  labs(
    x = "Filter levels",
    y = "Predicted number of somatic mutations",
    title = "Predicted number of somatic mutations per biopsy over different filters"
  )


show(plt)



# Assuming you have a vector containing the DNA bases in the 3' UTR
utr_bases <- c("A", "T", "C", "G")

# Assuming you have a vector containing the counts of each base
utr_counts <- c(10, 12, 8, 15)

# Creating the barplot
barplot(utr_counts, names.arg = utr_bases, xlab = "Base", ylab = "Count", main = "Base Distribution in 3' UTR")




#A to T : 5433
#A to C : 4182
#A to G : 29145
#T to A: 5584
#T to C: 25823
#T to G: 3911
#C to A: 5957
#C to G: 5517
#C to T: 18568
#G to A: 20143
#G to C: 5916
#G to T: 5780
# g inserted : 2672
# t inserted : 34293
# a inserted : 34806
# c inserted : 4234
# G deleted :16507
# T deleted :15098
# A deleted : 13805
# C deleted : 16455



# non-circular stacked barplot for somatic mutation spread 
data <- data.frame(
  individual=c('T', 'A', 'C', 'G'),
  group=c(rep('T, n=30882', 4), rep('A, n=33551', 4),rep('C, n=26164', 4), rep('G, n=31839', 4), rep('Del, n=15203', 4), rep('Ins, n=76005', 4)), 
  value=c(0, 0.16, 0.73, 0.11, 0.14, 0, 0.11, 0.75, 0.62, 0.20, 0, 0.18, 0.18, 0.63, 0.19, 0, 0.24, 0.22,0.27, 0.27, 0.44, 0.48, 0.05, 0.03)
)

data <- subset(data, value != 0)
rownames(data) <- NULL 
data <- data %>%
  mutate(group = factor(group, levels=c("T, n=30882", "A, n=33551", "C, n=26164", "G, n=31839", "Ins, n=76005", "Del, n=15203")))


#ggplot(data, aes(x=group, y=value, fill=individual)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
#  geom_bar(position = 'stack', stat="identity", alpha=0.5)


library(ggpattern)
stacked_bar <- ggplot(data, aes(x=group, y=value, fill=individual, pattern=group)) +
  geom_bar_pattern(position = 'stack', stat="identity", alpha=1,
  pattern_color = 'white',
  pattern_fill = "white",
  pattern_angle = 45,
  pattern_density = 0.05,
  pattern_spacing = 0.025,
  pattern_key_scale_factor = 0.6) + 
  scale_fill_manual(values=c("#7E549E", "#C2549D", "#FC8370", "#FECB3E")) +
  scale_pattern_manual(values = c("none", "none", "none", "none", "stripe", "stripe"))+
  theme_classic()+
  theme(
    axis.title.y = element_text(size=14, face="bold", colour = "black"),    
    axis.title.x = element_blank(),
    axis.text.x = element_text(size=14, face="bold", colour = "black"),
    axis.text.y = element_text(size=12, face="bold", colour = "black", hjust = 0),
    legend.text = element_text(face="bold"),
    legend.title = element_text(face = "bold"),
    plot.title=element_text(face = "bold"))+
  #theme(legend.position = "none")+
  ylab("Fraction of all mutations")+
  xlab("Type of mutation")
  #scale_color_manual(values = c())
stacked_bar





subs_dels <- c()
for (mut in data$group)
{
  if(mut ==  "T, n=30882" || mut ==  "A, n=33551" || mut ==  "C, n=26164" || mut ==  "G, n=31839"){
    subs_dels <- c(subs_dels, "Substitution")
  } else {
    subs_dels <- c(subs_dels, "Indels")
  }
}
data$subs_dels <- subs_dels

stacked_bar <- ggplot(data, aes(x=group, y=value, fill=individual, pattern=group)) +
  geom_bar(position = 'stack', stat="identity",  alpha=0.8) +
  scale_fill_manual(values=c("#7E549E", "#C2549D", "#FC8370", "#FECB3E")) +
  theme_bw()+
  theme(
    axis.title.y = element_text(size=14, face="bold", colour = "black"),    
    axis.title.x = element_blank(),
    axis.text.x = element_text(size=14, face="bold", colour = "black"),
    axis.text.y = element_text(size=12, face="bold", colour = "black", hjust = 0),
    legend.text = element_text(face="bold"),
    legend.title = element_text(face = "bold"),
    plot.title=element_text(face = "bold"))+
    
  #theme(legend.position = "none")+
  facet_grid( ~ subs_dels, scales = "free", space='free') +
  ylab("Fraction of all mutations")+
  xlab("Type of mutation")
#scale_color_manual(values = c())
stacked_bar







# mutation overview 


# horizontal bar to cover exons # see if ic an get broken axis or not 
# no broken axis so multiple subplots 


# Add stems with different colors and shapes for different mutations ? not bad



# add density for miRNA binding 







plot_df <- read.csv("/Users/iwan/Research/Somatic_mutations/output_data/somatic_mutations_post_filter.csv")
for (gene in relevant_p_filter$A){
  p <- ggplot(plot_df, aes_string(x="Diagnosis", y=gene,  fill="Diagnosis")) + 
    geom_violin(aes(color = Diagnosis), trim = FALSE, position = position_dodge(0.9), alpha=.10, size=0.8) + stat_n_text() +
    geom_boxplot(aes(color=Diagnosis), alpha=0, width=0.3, size=0.8)+
    theme_bw()+
    labs(y = "Number of somatic mutations per sample") +
    geom_jitter(aes(color=Diagnosis), alpha=0.4, width=0.1, height = 0.1) +
    ggtitle(gene)  
  #$#geom_signif(comparisons = list(c("Yes", "No")), map_signif_level=TRUE)
  print(p)
}


plot_df$Inflammation

plot_df <- plot_df %>%
  mutate(Inflammation = factor(Inflammation, levels=c("Yes", "No")))



library(ggpubr)
library(ggsignif)



edaradd_plot_df <- plot_df[! is.na(plot_df$EDARADD_gene_expression),]
p <- ggplot(edaradd_plot_df, aes_string(x="Inflammation", y="EDARADD",  fill="Inflammation")) + 
  geom_violin(aes(color = Inflammation), trim=FALSE, alpha=0.7, size=0.1, color='white', width=0.8, scale='width')  +
  scale_fill_manual(values=c("#DD2314", "#2A893B")) +
  geom_boxplot(alpha=1, width=0.09, size=0.4, fill='white', color='white', outlier.alpha=0)+
  theme_classic()+
  labs(y = "Number of somatic mutations per sample") +
  
  #geom_jitter(aes(color=Inflammation), alpha=0.8, width=0.1, height = 0.1, size=3) +
  ylim(0,5) +
  ggtitle("EDARADD") + guides(fill="none") + 
  geom_signif(comparisons = list(c("Yes", "No")), map_signif_level=TRUE)
p


#

inflammation_plot <- function(gene, max){
  
  p <- ggplot(plot_df, aes_string(x="Inflammation", y=gene,  fill="Inflammation")) + 
    geom_violin(aes(color = Inflammation), trim=FALSE, alpha=0.7, size=0.1, color='white', width=0.8, scale='width')  +
    scale_fill_manual(values=c("#DD2314", "#2A893B")) +
    geom_boxplot(alpha=1, width=0.09, size=0.4, fill='white', color='white', outlier.alpha=0)+
    theme_classic()+
    labs(y = "Number of somatic mutations per sample") +
    
    #geom_jitter(aes(color=Inflammation), alpha=0.8, width=0.1, height = 0.1, size=3) +
    ylim(0,max) +
    ggtitle(gene) + guides(fill="none") +
    geom_signif(comparisons = list(c("Yes", "No")), map_signif_level=TRUE) +
    scale_x_discrete(breaks=c("Yes","No"),
                     labels=c("Yes, n=192", "No, n=338"))
  return(p)
}




#creating all inflammation plots
inflammation_plot("CXCL1", 9)
inflammation_plot("SLC19A3", 12)
inflammation_plot("CXCL3", 4)
# not working
inflammation_plot("IGKV1.8", 45)

inflammation_plot("IL1B", 5)
inflammation_plot("IGLV3.27", 45)

# technically ns
inflammation_plot("DRAM1", 10)

inflammation_plot("FPR1", 3)
inflammation_plot("DUOX2", 12)
# does not exist
inflammation_plot("GBP1", 10)

inflammation_plot("VNN3", 3)

# whyyy
inflammation_plot("EDARADD", 30)

inflammation_plot("IGKV5.2", 45)
inflammation_plot("REEP5", 18)
inflammation_plot("ELN", 7)


location <- c()
for (loc in plot_df$Location_rough)
{
  if(loc ==  "ileum"){
    location <- c(location, "Ileum")
  } else {
    location <- c(location, "Colon")
  }
}
plot_df$Location <- location

plot_df <- plot_df %>%
  mutate(Location = factor(Location, levels=c("Ileum", "Colon")))




location_plot <- function(gene, max){
  #gene_plot_df <- plot_df[! is.na(plot_df[,glue::glue("{gene}_gene_expression")])]
  
  p <- ggplot(plot_df, aes_string(x="Location", y=gene,  fill="Location")) + 
    scale_fill_manual(values=c("#FFCF20", "#FF7020")) +
    geom_violin(aes(color = Inflammation), trim=FALSE, alpha=0.7, size=0.1, color='white', width=0.8, scale='width') +
    geom_boxplot(alpha=1, width=0.09, size=0.4, fill='white', color='white', outlier.alpha=0)+
    theme_classic()+
    labs(y = "Number of somatic mutations per sample") +
    ylim(0,max) +
    ggtitle(gene) + guides(fill="none") + 
    geom_signif(comparisons = list(c("Ileum", "Colon")), map_signif_level=TRUE) +
    scale_x_discrete(breaks=c("Ileum","Colon"),
                     labels=c("Ileum, n=181", "Colon, n=349"))
  return(p)
}


# all working location plots
location_plot("BMX", 10)
location_plot("REG3A", 6)
location_plot("APOA4", 8)
location_plot("APOB", 12)
location_plot("BPNT1", 18)
location_plot("MS4A12", 11)
location_plot("ZDHHC20", 70)
location_plot("SLC7A9", 20)
location_plot("DBNL", 13)
location_plot("CASP10", 7)
location_plot("TMEM236", 20)
location_plot("ITGA7", 10)
location_plot("TNK2", 15)
location_plot("SCP2", 13)
location_plot("MRPL30", 14)
location_plot("PYGB", 7)
location_plot("PREPL", 13)
location_plot("ENPP7", 3)

location_plot("CEACAM1", 15)
location_plot("CEACAM7", 7)
location_plot("NDUFS1", 14)
location_plot("SLC29A2", 15)
location_plot("SLC19A3", 11)
location_plot("MRPS16", 11)
location_plot("PEPD", 7)
location_plot("PEBP1P2", 11)
location_plot("SLC7A7", 7)
location_plot("CEACAM5", 17)
location_plot("CUBN", 8)
location_plot("ALDOB", 9)
location_plot("GGA2", 7)
location_plot("NOP14", 15)
location_plot("SLC28A2", 11)
location_plot("FGFBP1", 3)

# not signif on plot
location_plot("TLX1", 10)



## all not working location plots
location_plot("DEFA5", 10)
##

gene_plot_df <- plot_df[! is.na(plot_df[,("DEFA5_gene_expression")])]
gene_plot_df <- plot_df[! is.na(plot_df$DEFA5_gene_expression),]
gene_plot_df <- plot_df[! is.na(plot_df[,glue::glue("DEFA5_gene_expression")])]
gene_plot_df <- plot_df[! is.na(plot_df[,glue::glue("SLC28A2_gene_expression")])]

  
p <- ggplot(gene_plot_df, aes_string(x="Location", y="DEFA5",  fill="Location")) + 
  geom_violin(aes(color = Inflammation), trim=FALSE, alpha=0.7, size=0.1, color='white', width=0.8, scale='width')  +
  scale_fill_manual(values=c("#FFCF20", "#FF7020")) +
  geom_boxplot(alpha=1, width=0.09, size=0.4, fill='white', color='white', outlier.alpha=0)+
  theme_classic()+
  labs(y = "Number of somatic mutations per sample") +
  ylim(0,10) +
  ggtitle("gene") + guides(fill="none") + 
  geom_signif(comparisons = list(c("Ileum", "Colon")), map_signif_level=TRUE)
p


gene_plot_df <- plot_df[! is.na(plot_df$ACE_gene_expression),]
p <- ggplot(gene_plot_df, aes_string(x="Location", y="DEFA5",  fill="Location")) + 
  geom_violin(aes(color = Inflammation), trim=FALSE, alpha=0.7, size=0.1, color='white', width=0.8, scale='width')  +
  scale_fill_manual(values=c("#FFCF20", "#FF7020")) +
  geom_boxplot(alpha=1, width=0.09, size=0.4, fill='white', color='white', outlier.alpha=0)+
  theme_classic()+
  labs(y = "Number of somatic mutations per sample") +
  ylim(0,10) +
  ggtitle("gene") + guides(fill="none") + 
  geom_signif(comparisons = list(c("Ileum", "Colon")), map_signif_level=TRUE)
p
p <- ggplot(gene_plot_df, aes_string(x="Location", y="ACE",  fill="Location")) + 
  geom_violin(aes(color = Inflammation), trim=FALSE, alpha=0.7, size=0.1, color='white', width=0.8, scale='width')  +
  scale_fill_manual(values=c("#FFCF20", "#FF7020")) +
  geom_boxplot(alpha=1, width=0.09, size=0.4, fill='white', color='white', outlier.alpha=0)+
  theme_classic()+
  labs(y = "Number of somatic mutations per sample") +
  ylim(0,400) +
  ggtitle("gene") + guides(fill="none") + 
  geom_signif(comparisons = list(c("Ileum", "Colon")), map_signif_level=TRUE)
p


location_plot("ACE", 10)
##
location_plot("C11orf54", 10)
##
location_plot("NT5C2", 10)
##
location_plot("FSIP2", 10)
##
location_plot("BDH2P1", 10)
##
location_plot("TSGA10", 10)
##
location_plot("TM4SF5", 10)
##
location_plot("LY75.CD302", 10)
##
location_plot("FGFBP1", 10)
##
location_plot("COX15", 10)





library(ggbreak) 
library(patchwork)


pathogenic_overview <- read.csv('/Users/iwan/Research/Somatic_mutations/plot_dataframes/patho_plot_df', sep =',', header = TRUE)
#colnames(pathogenic_overview) <- c('Counts', 'Genes')

pathogenic_overview$X <-NULL
# stacked bar with colors of type of mutation # easy 

#barplot(sort(pathogenic_overview$Counts, decreasing = TRUE))

indel_mut <- pathogenic_overview[pathogenic_overview$Mutation_type == 'indel',]
point_mut <- pathogenic_overview[pathogenic_overview$Mutation_type == 'point',]
pathogenic_overview <- pathogenic_overview %>%
  mutate(Mutation_type = factor(Mutation_type, levels=c("indel", "point")))

point_mut <- point_mut %>%
  mutate(Genes = factor(Genes, levels=c('ARV1', 'ZMPSTE24', 'CNGA1', 'OFD1', 'CD3G', 'ANO10', 'GNPTAB',
                                        'RAD50', 'ALG6', 'IL7R', 'NBN', 'AHDC1', 'CD2AP', 'P3H1', 'MUT',
                                        'PALB2', 'FBXL4','PROM1', 'PTEN', 'TRIO', 'ANKRD11', 'MTM1', 
                                        'SGSH', 'TOP3A', 'ALMS1', 'EXT1', 'SLC26A2', 'AHCY', 'ARSA', 'BTD', 'CEP290',
                                        'CRYAB', 'DNM2', 'EIF2B5', 'EXT2', 'FBN1', 'FECH', 'HMGCL', 'HPS1', 'HPS4',
                                        'IDUA', 'IL10RA', 'IQSEC2', 'KIAA1109', 'KMT2B', 'KMT2E', 'MAN1B1', 'MTHFR',
                                        'MYRF', 'NSD2', 'PLA2G6', 'POLG', 'PTPN11', 'SGCE', 'SHANK3', 'SLC25A20', 'SOS1',
                                        'STIM1', 'TUBGCP6', 'UNC13D', 'ZSWIM6')))

indel_mut  <- indel_mut %>%
  mutate(Genes = factor(Genes, levels=c('ARV1', 'ZMPSTE24', 'CNGA1', 'OFD1', 'CD3G', 'ANO10', 'GNPTAB',
                                        'RAD50', 'ALG6', 'IL7R', 'NBN', 'AHDC1', 'CD2AP', 'P3H1', 'MUT',
                                        'PALB2', 'FBXL4','PROM1', 'PTEN', 'TRIO', 'ANKRD11', 'MTM1', 
                                        'SGSH', 'TOP3A', 'ALMS1', 'EXT1', 'SLC26A2', 'AHCY', 'ARSA', 'BTD', 'CEP290',
                                        'CRYAB', 'DNM2', 'EIF2B5', 'EXT2', 'FBN1', 'FECH', 'HMGCL', 'HPS1', 'HPS4',
                                        'IDUA', 'IL10RA', 'IQSEC2', 'KIAA1109', 'KMT2B', 'KMT2E', 'MAN1B1', 'MTHFR',
                                        'MYRF', 'NSD2', 'PLA2G6', 'POLG', 'PTPN11', 'SGCE', 'SHANK3', 'SLC25A20', 'SOS1',
                                        'STIM1', 'TUBGCP6', 'UNC13D', 'ZSWIM6')))


check<-merge(point_mut, indel_mut, by.x = "Genes", by.y = "Genes", all = T)
check$Value.x  <- check$Value.x + 24
check$Value.x  <- check$Value.x - 3

p1 <- ggplot(check, aes( x=Genes)) + 
  geom_bar(aes(fill=Mutations.y, y=Value.y), position="stack", stat="identity") +  
  scale_fill_manual(values=c("#7E549E", "#C2549D", "#FC8370", "7E549E")) +
  labs(y = "Number of observations", x='') +
  ggtitle('Observations of pathogenic mutations in IBD biopsies') +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 50, vjust = 1, hjust=1)) +
  geom_point(data=check, aes(x=Genes, y=Value.x, color=Mutations.x), show.legend=FALSE) +
  scale_color_manual(values=c("#7E549E", "#C2549D", "#FC8370", "7E549E")) 
  
  

p1
pg <- p1 + coord_cartesian(ylim(c(0,150))) #+ coord_cartesian(ylim(c(0,10)))
#pg
pg +
  facet_grid(Mutation_type~.,
             scales="free_y",
             space ="free_y") +



p2 <- ggplot(point_mut, aes(fill=Mutations, y=Value, x=reorder(Genes, -Value))) + 
  geom_bar(position="stack", stat="identity") +  
  scale_fill_manual(values=c("#7E549E", "#C2549D", "#FC8370", "7E549E")) +
  labs(y = "<= point mutations, indels =>", x='') +
  ggtitle('Observations of pathogenic mutations') +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 50, vjust = 1, hjust=1)) 

p1 + p2

#"#C2549D", "7E549E", "#FC8370", "#FECB3E", "#7E549E"








# table one and two code
table_2_df <- read.csv("/Users/iwan/Research/Somatic_mutations/output_data/somatic_mutations_post_filter.csv")
library(table1)
table_2_df$Sex <- 
  factor(table_2_df$Sex, levels=c(1,2),
         labels=c("Male", 
                  "Female"))

table_2_df$Diagnosis <- 
  factor(table_2_df$Diagnosis, levels=c("CD","UC"),
         labels=c("CD", 
                  "UC"))

table_2_df$Location_rough <- 
  factor(table_2_df$Location_rough, levels=c("colon","ileum"),
         labels=c("Colon", 
                  "Ileum"))

table_2_df$tnf_before <- 
  factor(table_2_df$tnf_before, levels=c(0,1),
         labels=c("No", 
                  "Yes"))

table_2_df$tnf_during <- 
  factor(table_2_df$tnf_during, levels=c(0,1),
         labels=c("No", 
                  "Yes"))

table_2_df$Montreal_B <- 
  factor(table_2_df$Montreal_B, levels=c(1,2,3),
         labels=c("Non stricturing, nonpenetrating", 
                  "Stricturing",
                  "Penetrating"))

table_2_df$Montreal_L <- 
  factor(table_2_df$Montreal_L, levels=c(1,2,3),
         labels=c("Terminal Ileum", 
                  "Colon",
                  "Ileocolon"))

table_2_df$Inflammation <- 
  factor(table_2_df$Inflammation, levels=c("Yes", "No"),
         labels=c("Inflamed", 
                  "Non-inflamed"))

label(table_2_df$tnf_before) <- "Anti-tnf alpha use before biopsy"
label(table_2_df$tnf_during) <- "Anti-tnf alpha use at time of biopsy"
label(table_2_df$Diagnosis) <- "Diagnosis"
label(table_2_df$Sex) <- "Sex"
label(table_2_df$Montreal_B) <- "Montreal Behavior"
label(table_2_df$Montreal_L) <- "Montreal Location"
label(table_2_df$Location_rough) <- "Location of biopsy"
label(table_2_df$Age_at_biopsy) <- "Age"
caption <- "Table 2, Sample description"
units(table_2_df$Age_at_biopsy) <- "Years"


table1(~  Sex + Age_at_biopsy + BMI + Diagnosis + 
         Location_rough + tnf_during + tnf_before + 
         Montreal_L + Montreal_B | Inflammation, data=table_2_df, caption=caption)









library(table1)


table_1_df <- read.csv("/Users/iwan/Research/Somatic_mutations/output_data/single_patients_somatic_mutations.csv")

table_1_df$Sex <- 
  factor(table_1_df$Sex, levels=c(1,2),
         labels=c("Male", 
                  "Female"))

table_1_df$Diagnosis <- 
  factor(table_1_df$Diagnosis, levels=c("CD","UC"),
         labels=c("CD", 
                  "UC"))

table_1_df$tnf_before <- 
  factor(table_1_df$tnf_before, levels=c(0,1),
         labels=c("No", 
                  "Yes"))

table_1_df$tnf_during <- 
  factor(table_1_df$tnf_during, levels=c(0,1),
         labels=c("No", 
                  "Yes"))

table_1_df$Montreal_B <- 
  factor(table_1_df$Montreal_B, levels=c(1,2,3),
         labels=c("Non stricturing, nonpenetrating", 
                  "Stricturing",
                  "Penetrating"))

table_1_df$Montreal_L <- 
  factor(table_1_df$Montreal_L, levels=c(1,2,3),
         labels=c("Terminal Ileum", 
                  "Colon",
                  "Ileocolon"))


label(table_1_df$tnf_before) <- "Anti-tnf alpha use before biopsy"
label(table_1_df$tnf_during) <- "Anti-tnf alpha use at time of biopsy"
label(table_1_df$Diagnosis) <- "Diagnosis"
label(table_1_df$Sex) <- "Sex"
label(table_1_df$Montreal_B) <- "Montreal Behavior"
label(table_1_df$Montreal_L) <- "Montreal Location"
label(table_1_df$Location_rough) <- "Location of biopsy"
label(table_1_df$Age_at_biopsy) <- "Age"
caption_1 <- "Table 1, Patient description"
units(table_1_df$Age_at_biopsy) <- "Years"


table1(~  Sex + Age_at_biopsy + BMI + tnf_during + tnf_before + 
         Montreal_L + Montreal_B | Diagnosis, data=table_1_df, caption=caption)



pvalue <- function(x, ...) {
  # Construct vectors of data y, and groups (strata) g
  y <- unlist(x)
  g <- factor(rep(1:length(x), times=sapply(x, length)))
  if (is.numeric(y)) {
    # For numeric variables, perform a standard 2-sample t-test
    p <- t.test(y ~ g)$p.value
  } else {
    # For categorical variables, perform a chi-squared test of independence
    p <- chisq.test(table(y, g))$p.value
  }
  # Format the p-value, using an HTML entity for the less-than sign.
  # The initial empty string places the output on the line below the variable label.
  c("", sub("<", "&lt;", format.pval(p, digits=3, eps=0.001)))
}

table1(~  Sex + Age_at_biopsy + BMI + tnf_during + tnf_before + 
         Montreal_L + Montreal_B | Diagnosis, data=table_1_df, caption=caption_1,overall=F, extra.col=list(`P-value`=pvalue))


table1(~  Sex + Age_at_biopsy + BMI + Diagnosis + 
         Location_rough + tnf_during + tnf_before + 
         Montreal_L + Montreal_B | Inflammation, data=table_2_df, caption=caption,overall=F, extra.col=list(`P-value`=pvalue))





# stacked bar for mutations in XIAP, one column for point, one column for ins and
# dels coloured by bases split ins and dels also show a distribution of bases in UTR for comparison




xiap_muts <- read.csv('/Users/iwan/Research/Somatic_mutations/plot_dataframes/XIAP_bases_plot_df.csv')


#subs_dels <- c()
#for (mut in data$group)
#{
#  if(mut ==  "T, n=30882" || mut ==  "A, n=33551" || mut ==  "C, n=26164" || mut ==  "G, n=31839"){
#    subs_dels <- c(subs_dels, "Substitution")
#  } else {
#    subs_dels <- c(subs_dels, "Indels")
#  }
#$}
#data$subs_dels <- subs_dels
xiap_muts

xiap_muts <- xiap_muts %>%
  mutate(Mutation_type = factor(Mutation_type, levels=c("XIAP bases", "Insertion", "A", "T", "C", "G")))

xiap_muts <- xiap_muts %>%
  mutate(Base = factor(Base, levels=c("A", "C", "G", "T")))

muts <- c()
for (mut in xiap_muts$Mutation_type)
{
  print(mut)
  if(mut ==  "Insertion" || mut ==  "A" || mut ==  "C" || mut ==  "G" || mut ==  "T"){
    muts <- c(muts, "Mutation Distribution")
    }
  else {
      muts <- c(muts, "XIAP base distribution")
  }
}
xiap_muts$muts <- muts


muts <- c()
for (mut in xiap_muts$Mutation_type)
{
  print(mut)
  if(mut ==  "Insertion"  ){
    muts <- c(muts,"Insertion, n=501")
  } else if (mut == "A"){
    muts <- c(muts,"A, n=7020")
    
  } else if (mut == "C"){
    muts <- c(muts, "C, n=17")
  }else if (mut =="G" ){
    muts <- c(muts,"G, n=16")
  }else if (mut =="T" ){
    muts <- c(muts,"T, n=34")
  }
  else {
    muts <- c(muts, "XIAP bases, n=8413")
  }
}
xiap_muts$Mutation_type_label <- muts

xiap_muts <- xiap_muts %>%
  mutate(Mutation_type_label = factor(Mutation_type_label, levels=c("XIAP bases, n=8413", "Insertion, n=501", "A, n=7020", "T, n=34", "C, n=17", "G, n=16")))


stacked_bar <- ggplot(xiap_muts, aes(x=Mutation_type_label, y=count, fill=Base)) +
  geom_bar(position = 'stack', stat="identity",  alpha=0.8) +
  scale_fill_manual(values=c("#7E549E", "#C2549D", "#FC8370", "#FECB3E")) +
  theme_bw()+
  theme(
    axis.title.y = element_text(size=14, face="bold", colour = "black"),    
    axis.title.x = element_blank(),
    axis.text.x = element_text(size=14, face="bold", colour = "black"),
    axis.text.y = element_text(size=12, face="bold", colour = "black", hjust = 0),
    legend.text = element_text(face="bold"),
    legend.title = element_text(face = "bold"),
    plot.title=element_text(face = "bold"))+
  
  #theme(legend.position = "none")+
  facet_grid( ~ muts, scales = "free", space='free') +
  ylab("Observed number of mutations")+
  xlab("Type of mutation")
#scale_color_manual(values = c())
stacked_bar


stacked_bar <- ggplot(xiap_muts, aes(x=Mutation_type_label, y=fraction, fill=Base)) +
  geom_bar(position = 'stack', stat="identity",  alpha=0.8) +
  scale_fill_manual(values=c("#7E549E", "#C2549D", "#FC8370", "#FECB3E")) +
  theme_bw()+
  theme(
    axis.title.y = element_text(size=14, face="bold", colour = "black"),    
    axis.title.x = element_blank(),
    axis.text.x = element_text(size=14, face="bold", colour = "black"),
    axis.text.y = element_text(size=12, face="bold", colour = "black", hjust = 0),
    legend.text = element_text(face="bold"),
    legend.title = element_text(face = "bold"),
    plot.title=element_text(face = "bold"))+
  
  #theme(legend.position = "none")+
  facet_grid( ~ muts, scales = "free", space='free') +
  ylab("Fraction of all mutations")+
  xlab("Type of mutation")
#scale_color_manual(values = c())
stacked_bar

