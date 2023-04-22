## Loading packages

library(devtools)
library(car)
library(rlang)
library(dplyr)
library(phyloseq)
library(tidyverse)
library(plyr)
library(vegan)
library(stats)
library(ggplot2); theme_set(theme_bw()+
                              theme(panel.grid = element_blank(),
                                    strip.background = 
                                      element_rect(fill="white"),
                                    axis.text.y = 
                                      element_text(size = 12, colour = "black"),
                                    axis.title = 
                                      element_text(size=15, face="bold"),
                                    strip.text = 
                                      element_text(color="black", size=10),
                                    legend.text=element_text(size=10),
                                    axis.line = element_line(colour = "black"),
                                    axis.text.x = element_blank(),))




## Setting working directory

setwd("C:/Users/arach/OneDrive/Desktop/R")
filepath = "C:/Users/arach/OneDrive/Desktop/R"





## read data

PAE_filtered = readRDS("Rarefied_PAE_Dataset.RDS")

PAE_filtered@sam_data$number_of_reads = sample_sums(PAE_filtered)



## Setting dephyloseq function 

dephyloseq = function(phylo_obj){
  meta = as.data.frame(as.matrix(phylo_obj@sam_data))
  metacols = ncol(meta)+1
  otu = as.data.frame(t(as.matrix(phylo_obj@otu_table)))
  mo = merge(meta, otu, by=0)
  tax = as.data.frame(phylo_obj@tax_table)
  tax = tax %>% rownames_to_column(var="asv_id")
  mo = mo %>% pivot_longer(cols = -c(1:metacols), 
                           names_to = "asv_id", 
                           values_to="asv_abundance")
  mot = full_join(mo, tax)
  output = mot
}

## Further filtering

PAE_filtered_ratio = subset_taxa(PAE_filtered, 
                                 Rank2 %in% c("Bacteroidota", "Firmicutes"))


## Get filtered PAE data out of phyloseq

mot_filtered = dephyloseq(PAE_filtered_ratio)

mot_filtered$relative_abundance = 
  as.numeric(mot_filtered$asv_abundance)/
  as.numeric(mot_filtered$number_of_reads)

mot_filtered$plot_names = paste0(mot_filtered$Rank2, "; ", mot_filtered$Rank3)

## set days to numeric values

mot_filtered$Age[mot_filtered$Age == "P8"] <- "8"
mot_filtered$Age[mot_filtered$Age == "P22"] <- "22"
mot_filtered$Age[mot_filtered$Age == "P38"] <- "38"
mot_filtered$Age = as.numeric(as.character(mot_filtered$Age))

view(mot_filtered)


###############################################################################

## 2-factor ANOVA - 8 day

PAE_filtered_ratio_day8 = subset(mot_filtered, Age=="8")

pp8=ggplot(data=PAE_filtered_ratio_day8, aes(x=Group, y=asv_abundance))+
  geom_boxplot()
pp8

ps8=ggplot(data=PAE_filtered_ratio_day8, aes(x=Sex, y=asv_abundance))+
  geom_boxplot()
ps8

a8 <- aov(asv_abundance ~ Sex + Group, data = PAE_filtered_ratio_day8)
summary(a8)

#               Df Sum Sq Mean Sq F value  Pr(>F)   
# Sex            1    204  203.50   7.727 0.00547 **
# Group          1    148  148.18   5.626 0.01775 * 
# Residuals   3285  86516   26.34    

plot(a8)

leveneTest(asv_abundance ~ Sex*Group, data = PAE_filtered_ratio_day8)

# Levene's Test for Homogeneity of Variance (center = median)
#         Df F value    Pr(>F)    
# group    3  6.6042 0.0001907 ***
#       3284              

###############################################################################

## 2-factor ANOVA - 22 day

PAE_filtered_ratio_day22 = subset(mot_filtered, Age=="22")

pp22=ggplot(data=PAE_filtered_ratio_day22, aes(x=Group, y=asv_abundance))+
  geom_boxplot()
pp22

ps22=ggplot(data=PAE_filtered_ratio_day22, aes(x=Sex, y=asv_abundance))+
  geom_boxplot()
ps22

a22 <- aov(asv_abundance ~ Sex + Group, data = PAE_filtered_ratio_day22)
summary(a22)

#               Df   Sum Sq Mean Sq F value Pr(>F)
# Sex            1      498     498   0.128  0.721
# Group          1      688     688   0.177  0.674
# Residuals   4381 17078208    3898 

plot(a22)

leveneTest(asv_abundance ~ Sex*Group, data = PAE_filtered_ratio_day22)

# Levene's Test for Homogeneity of Variance (center = median)
#         Df F value Pr(>F)
# group    3  0.1533 0.9276
#       4380   

###############################################################################

## 2-factor ANOVA - 38 day 

PAE_filtered_ratio_day38 = subset(mot_filtered, Age=="38")

pp38=ggplot(data=PAE_filtered_ratio_day38, aes(x=Group, y=asv_abundance))+
  geom_boxplot()
pp38

ps38=ggplot(data=PAE_filtered_ratio_day38, aes(x=Sex, y=asv_abundance))+
  geom_boxplot()
ps38

a38 <- aov(asv_abundance ~ Sex + Group, data = PAE_filtered_ratio_day38)
summary(a38)

#               Df   Sum Sq Mean Sq F value Pr(>F)
# Sex            1        6       6   0.001  0.973
# Group          1       32      32   0.006  0.938
# Residuals   3970 21255572    5354     

plot(a38)

leveneTest(asv_abundance ~ Sex*Group, data = PAE_filtered_ratio_day38)

# Levene's Test for Homogeneity of Variance (center = median)
#         Df F value Pr(>F)
# group    3   0.082 0.9699
#       3969      

