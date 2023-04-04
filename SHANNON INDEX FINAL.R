##### WORKSPACE SETUP #####
library(phyloseq)
library(vegan)
library(stats)
library(ggplot2)
library(tidyverse)

### Setting working directory ###
setwd("/Users/yoovn/Desktop/BIOL_403/BIOL403_LAB/PAE_Project")

###VARIABLES###
readRDS("Rarefied_PAE_Dataset.RDS")
PAE_data = readRDS("Rarefied_PAE_Dataset.RDS")

#DEPHYLOSEQ
metadata = as.data.frame(as.matrix(PAE_data@sam_data))
otu = as.data.frame(t(as.matrix(PAE_data@otu_table)))

otu_table1 <- otu_table(otu, taxa_are_rows = FALSE)
View(otu_table1)

# merge the metadata and otu table together
metaasv = merge(metadata, otu, by=0)

View(metaasv)
View(metadata)


###
Day8_only_asv <- subset(metaasv, Age=="P8")
Day22_only_asv <- subset(metaasv, Age=="P22")
Day38_only_asv <- subset(metaasv, Age=="P38")

#Day 8 Only

PAE_Day8 <- Day8_only_asv[,-ncol(Day8_only_asv)]
PAE_meta_col <-c(1,2,3,4,5,6,7)
PAE_otu_col <-8:172 
view(PAE_otu_col)

#Day 22

PAE_Day22 <- Day22_only_asv[,-ncol(Day22_only_asv)]
PAE_meta_col <-c(1,2,3,4,5,6,7)
PAE_otu_col <-8:172 
view(PAE_otu_col)

#Day 38

PAE_Day38 <- Day38_only_asv[,-ncol(Day38_only_asv)]
PAE_meta_col <-c(1,2,3,4,5,6,7)
PAE_otu_col <-8:172 
view(PAE_otu_col)

######### RICHNESS ######## 
metadata$obs_richness = estimate_richness(otu_table1, split = TRUE, measures = c("Observed"))
view(metadata)
metadata$alpha <- diversity(otu, index = "shannon")

day8_only = subset(metadata, Age=="P8")
day22_only = subset(metadata, Age=="P22")
day38_only = subset(metadata, Age=="P38")

############ SHANNON'S INDEX #############

###All Days###
all_days_graph <- ggplot(metadata, aes(x = Group, y = alpha, fill=Sex)) + 
  geom_boxplot() + xlab("Group") +ylab("Shannon's Index") + 
  labs(title = " Figure 1: Alpha Diversity (α) of All Days Combined") +
  theme(plot.title = element_text(hjust = 0.5)) #+ geom_signif(
all_days_graph

###Day 8###
day8_only_graph <-ggplot(day8_only, aes(x = Group, y = alpha, fill= Sex)) + xlab("Group") +ylab("Shannon's Index") + 
  labs(title = " Figure 2: Alpha Diversity (α) at Day 8 of Treatment") +
  theme(plot.title = element_text(hjust = 0.5)) +geom_boxplot()
day8_only_graph

###Day 22###
day22_only_graph <-ggplot(day22_only, aes(x = Group, y = alpha, fill=Sex)) + 
  geom_boxplot() + xlab("Group") +ylab("Shannon's Index") + 
  labs(title = " Figure 3: Alpha Diversity (α) at Day 22 of Treatment") +
  theme(plot.title = element_text(hjust = 0.5)) 
day22_only_graph

###Day 38###
day38_only_graph <-ggplot(day38_only, aes(x = Group, y = alpha, fill=Sex)) + 
  geom_boxplot() + xlab("Group") +ylab("Shannon's Index") + 
  labs(title = " Figure 4: Alpha Diversity (α) at Day 38 of Treatment") +
  theme(plot.title = element_text(hjust = 0.5)) 
day38_only_graph

#SAVE FILE 1 PAGE#
pdf("PAE Alpha Plots no Title.pdf", width = 10, height = 8) # Set the name and size of the output PDF
colormodel = "cmyk"
gridExtra::grid.arrange(all_days_graph, day8_only_graph, day22_only_graph, day38_only_graph, ncol=2) # Use grid.arrange() function to arrange the plots on a grid
dev.off()


####SAVE FILE MULTIPLE PAGES####
pdf("PAE Alpha Diversity (α) Plots.pdf",         
    width = 8, height = 7, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk",    # Color model (cmyk is required for most publications)
    paper = "letter")          # Paper size
all_days_graph
day8_only_graph
day22_only_graph
day38_only_graph

### CLOSE ###
dev.off() 

#SAVE FILE 1 FIGURE#
# Combine all plots into a single figure
combined_plots <- all_days_graph + day8_only_graph + day22_only_graph + day38_only_graph + plot_layout(ncol = 2)

# Save the combined figure to a PDF
ggsave("combined_Alpha_plots.pdf", combined_plots)


##PERMANOVA 2-FACTOR##

##Run the PERMANOVA
p1 = adonis2(metaasv[,-c(1:7)]~Group + Sex, data=metaasv)
print(p1)

p2 = adonis2(Day8_only_asv[,-c(1:7)]~Group + Sex, data=Day8_only_asv)
print(p2)

p3 = adonis2(Day22_only_asv[,-c(1:7)]~Group + Sex, data=Day22_only_asv)
print(p3)

p4 = adonis2(Day38_only_asv[,-c(1:7)]~Group + Sex, data=Day38_only_asv)
print(p4)

