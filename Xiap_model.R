#R script for xiap model analysis
library(lme4) # for multilevel models


df <- read.csv('/Users/iwan/Research/Somatic_mutations/output_data/XIAP_mutations_for_R_post_filter.csv')


encode_ordinal <- function(x, order = unique(x)) {
  x <- as.numeric(factor(x, levels = order, exclude = NULL))
  x
}
df[["Inflammation"]] <- encode_ordinal(df[["Inflammation"]]) 
df[["Diagnosis"]] <- encode_ordinal(df[["Diagnosis"]]) 
df[["Location_rough"]] <- encode_ordinal(df[["Location_rough"]])
df[["Sequencing_Year"]] <- encode_ordinal(df[["Sequencing_Year"]])
df$Diagnosis <- df$Diagnosis -1
df$Inflammation <- df$Inflammation -1
df$Location_rough <- df$Location_rough -1
df$Sequencing_Year <- df$Sequencing_Year -1
df$Sex <- df$Sex -1

sumsum <- summary(lmer('XIAP ~ Inflammation + XIAP_gene_expression + tnf_combined + Sequencing_Year + Montreal_L + tnf_before + Montreal_B + Location_rough + Sex + Diagnosis+ Age_at_biopsy + BMI + (1|umcg_id)', data=df, na.action = na.omit))
coefs <- data.frame(coef(summary(sumsum)))
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs


newmodel <- summary(lmer('XIAP_gene_expression ~ Inflammation*XIAP+ tnf_combined + Sequencing_Year + Montreal_L + tnf_before + Montreal_B + Location_rough + Sex + Diagnosis+ Age_at_biopsy + BMI + (1|umcg_id)', data=df, na.action = na.omit))
newcoefs <- data.frame(coef(summary(newmodel)))
newcoefs$p.z <- 2 * (1 - pnorm(abs(newcoefs$t.value)))
newcoefs

newmodel <- summary(lmer('Inflammation ~ XIAP + XIAP_gene_expression + tnf_combined + Sequencing_Year + Montreal_L + tnf_before + Montreal_B + Location_rough + Sex + Diagnosis+ Age_at_biopsy + BMI + (1|umcg_id)', data=df, na.action = na.omit))
newcoefs <- data.frame(coef(summary(newmodel)))
newcoefs$p.z <- 2 * (1 - pnorm(abs(newcoefs$t.value)))
newcoefs


my_model <- lmer('XIAP ~ Inflammation + XIAP_gene_expression + XIAP : XIAP_gene_expression + tnf_combined + Sequencing_Year + Montreal_L + tnf_before + Montreal_B + Location_rough + Sex + Diagnosis+ Age_at_biopsy + BMI + (1|umcg_id)', data=df, na.action = na.omit)
anova(my_model) 
drop1(my_model,test="Chisq")

inflammation_plot <- function(gene, max){
  gene_plot_df <- plot_df[! is.na(plot_df[,glue::glue("{gene}_gene_expression")])]
  p <- ggplot(gene_plot_df, aes_string(x="Inflammation", y=gene,  fill="Inflammation")) + 
    geom_violin(aes(color = Inflammation), trim=FALSE, alpha=0.7, size=0.1, color='white', width=0.8, scale='width')  +
    scale_fill_manual(values=c("#DD2314", "#2A893B")) +
    geom_boxplot(alpha=1, width=0.09, size=0.4, fill='white', color='white', outlier.alpha=0)+
    theme_classic()+
    labs(y = "Number of somatic mutations per sample") +
    
    #geom_jitter(aes(color=Inflammation), alpha=0.8, width=0.1, height = 0.1, size=3) +
    ylim(0,max) +
    ggtitle(gene) + guides(fill="none") + 
    geom_signif(comparisons = list(c("Yes", "No")), map_signif_level=TRUE)
  return(p)
}

p <- ggplot(df, aes_string(x="Location_rough", y="XIAP",  fill="Location_rough")) + 
  geom_violin(aes(color = Inflammation), trim=FALSE, alpha=0.7, size=0.1, color='white', width=0.8, scale='width')  +
  scale_fill_manual(values=c("#DD2314", "#2A893B")) +
  geom_boxplot(alpha=1, width=0.09, size=0.4, fill='white', color='white', outlier.alpha=0)+
  theme_classic()+
  labs(y = "Number of somatic mutations per sample") +
  
  #geom_jitter(aes(color=Inflammation), alpha=0.8, width=0.1, height = 0.1, size=3) +
  ylim(0,70) +
  ggtitle("XIAP") + guides(fill="none") 
p
