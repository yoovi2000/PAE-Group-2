
### Loading packages ###

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




### Setting working directory ###

setwd("C:/Users/arach/OneDrive/Desktop/R")
filepath = "C:/Users/arach/OneDrive/Desktop/R"



## Further filtering for the box plots

# read data

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


PAE_filtered_ratio = subset_taxa(PAE_filtered, 
                                 Rank2 %in% c("Bacteroidota", "Firmicutes"))


# Get filtered PAE data out of phyloseq

mot_filtered = dephyloseq(PAE_filtered_ratio)

mot_filtered$relative_abundance = 
  as.numeric(mot_filtered$asv_abundance)/
  as.numeric(mot_filtered$number_of_reads)

mot_filtered$plot_names = paste0(mot_filtered$Rank2, "; ", mot_filtered$Rank3)

# set days to numeric values

mot_filtered$Age[mot_filtered$Age == "P8"] <- "8"
mot_filtered$Age[mot_filtered$Age == "P22"] <- "22"
mot_filtered$Age[mot_filtered$Age == "P38"] <- "38"
mot_filtered$Age = as.numeric(as.character(mot_filtered$Age))

view(mot_filtered)
 
###############################################################################

## Processing data for PAE box plot - Day 8

# filter the day

PAE_filtered_ratio_day8 = subset(mot_filtered, Age=="8")


# making datafile to add 'other firmicutes' group


mot.group8 = ddply(PAE_filtered_ratio_day8, c("Group"),
                  summarise,
                  N=length(relative_abundance),
                  total = sum(relative_abundance))


mot.sum8 = ddply(PAE_filtered_ratio_day8, c("Group", "plot_names"),
                summarise,
                N=length(relative_abundance),
                sum = sum(relative_abundance))

samplegroups8 = unique(mot.sum8$Group)

sorted8 = mot.sum8[order(-mot.sum8$sum),]

firm.df8 = NULL

# start loop

for(i in samplegroups8) {
  for(j in i) {
    ## subset dataframe by samples
    sample8 = subset(sorted8, sorted8$Group %in% c(j))
    
    ## get other firmicutes
    firm8 = sample8[c(1:2),]
    
    ## save list of firmicutes abundance taxa
    t.tmp8 <- firm8
    firm.df8 <- rbind.fill(firm.df8, t.tmp8)
    
    ## close loop 
  }
}

firm.df8$place = "keep"

mot.firm8 = full_join(PAE_filtered_ratio_day8, firm.df8)

mot.firm8$place = replace(mot.firm8$place, is.na(mot.firm8$place), "firm")

mot.firm8[mot.firm8$place == "firm",]$plot_names <- "Other Firmicutes"


# making datafile to create ratios

mot.group.pup8 = ddply(mot.firm8, c("Group","Pup_number"),
                       summarise,
                       N=length(relative_abundance),
                       C_Sum = sum(ifelse(plot_names == 
                                            "Firmicutes; Clostridia", 
                                          relative_abundance, 0)),
                       B_Sum = sum(ifelse(plot_names == 
                                            "Bacteroidota; Bacteroidia", 
                                          relative_abundance, 0)),
                       F_Sum = sum(ifelse(plot_names == 
                                            "Other Firmicutes", 
                                          relative_abundance, 0)),
                       CB_Ratio = C_Sum/B_Sum,
                       BF_Ratio = B_Sum/(F_Sum+C_Sum),
                       CF_Ratio = C_Sum/F_Sum)


## Making PAE boxplots - Day 8

# Bacteroidota:Clostridia (CB Ratio)

day8_PAE_CB_plot <-ggplot(mot.group.pup8, aes(x=as.character(Group), 
                                            y=as.numeric(CB_Ratio), 
                                            fill=as.factor(Group))) +
                        geom_boxplot() +
                        labs(x=" ", y="Ratio", fill="Group") +
                        labs(title = 
                               "Bacteroidota:Clostridia Ratio for 8 Day Pups") +
                        theme(plot.title = element_text(hjust = 0.5))

day8_PAE_CB_plot

# Bacteroidota:All Firmicutes (BF Ratio)

day8_PAE_BF_plot <-ggplot(mot.group.pup8, aes(x=as.character(Group), 
                                            y=as.numeric(BF_Ratio), 
                                            fill=as.factor(Group))) +
  geom_boxplot() +
  labs(x=" ", y="Ratio", fill="Group") +
  labs(title = "Bacteroidota:All Firmicutes Ratio for 8 Day Pups") +
  theme(plot.title = element_text(hjust = 0.5))

day8_PAE_BF_plot

# Clostridia:Other Firmicutes  (CF Ratio)

day8_PAE_CF_plot <-ggplot(mot.group.pup8, aes(x=as.character(Group), 
                                              y=as.numeric(CF_Ratio), 
                                              fill=as.factor(Group))) +
  geom_boxplot() +
  labs(x=" ", y="Ratio", fill="Group") +
  labs(title = "Clostridia:Other Firmicutes Ratio for 8 Day Pups") +
  theme(plot.title = element_text(hjust = 0.5))

day8_PAE_CF_plot



## Processing data for Sex box plot - Day 8

# making datafile to add 'other firmicutes' group


mot.group.sex8 = ddply(PAE_filtered_ratio_day8, c("Sex"),
                   summarise,
                   N=length(relative_abundance),
                   total = sum(relative_abundance))


mot.sum.sex8 = ddply(PAE_filtered_ratio_day8, c("Sex", "plot_names"),
                 summarise,
                 N=length(relative_abundance),
                 sum = sum(relative_abundance))

samplegroups.sex8 = unique(mot.sum.sex8$Sex)

sorted.sex8 = mot.sum.sex8[order(-mot.sum.sex8$sum),]

firm.df.sex8 = NULL

# start loop

for(i in samplegroups.sex8) {
  for(j in i) {
    ## subset dataframe by samples
    sample.sex8 = subset(sorted.sex8, sorted.sex8$Sex %in% c(j))
    
    ## get other firmicutes
    firm.sex8 = sample.sex8[c(1:2),]
    
    ## save list of firmicutes abundance taxa
    t.tmp.sex8 <- firm.sex8
    firm.df.sex8 <- rbind.fill(firm.df.sex8, t.tmp.sex8)
    
    ## close loop 
  }
}

firm.df.sex8$place = "keep"

mot.firm.sex8 = full_join(PAE_filtered_ratio_day8, firm.df.sex8)

mot.firm.sex8$place = replace(mot.firm.sex8$place, 
                              is.na(mot.firm.sex8$place), "firm")

mot.firm.sex8[mot.firm.sex8$place == "firm",]$plot_names <- "Other Firmicutes"


# making datafile to create ratios

mot.group.pup.sex8 = ddply(mot.firm.sex8, c("Sex","Pup_number"),
                       summarise,
                       N=length(relative_abundance),
                       C_Sum = sum(ifelse(plot_names == 
                                            "Firmicutes; Clostridia", 
                                          relative_abundance, 0)),
                       B_Sum = sum(ifelse(plot_names == 
                                            "Bacteroidota; Bacteroidia", 
                                          relative_abundance, 0)),
                       F_Sum = sum(ifelse(plot_names == 
                                            "Other Firmicutes", 
                                          relative_abundance, 0)),
                       CB_Ratio = C_Sum/B_Sum,
                       BF_Ratio = B_Sum/(F_Sum+C_Sum),
                       CF_Ratio = C_Sum/F_Sum)


## Making Sex Boxplots - Day 8

# Bacteroidota:Clostridia (CB Ratio)

day8_Sex_CB_plot <-ggplot(mot.group.pup.sex8, aes(x=as.character(Sex), 
                                              y=as.numeric(CB_Ratio), 
                                              fill=as.factor(Sex))) +
  geom_boxplot() +
  labs(x=" ", y="Ratio", fill="Sex") +
  labs(title = "Bacteroidota:Clostridia Ratio for 8 Day Pups") +
  theme(plot.title = element_text(hjust = 0.5))

day8_Sex_CB_plot

# Bacteroidota:All Firmicutes (BF Ratio)

day8_Sex_BF_plot <-ggplot(mot.group.pup.sex8, aes(x=as.character(Sex), 
                                              y=as.numeric(BF_Ratio), 
                                              fill=as.factor(Sex))) +
  geom_boxplot() +
  labs(x=" ", y="Ratio", fill="Sex") +
  labs(title = "Bacteroidota:All Firmicutes Ratio for 8 Day Pups") +
  theme(plot.title = element_text(hjust = 0.5))

day8_Sex_BF_plot

# Clostridia:Other Firmicutes  (CF Ratio)

day8_Sex_CF_plot <-ggplot(mot.group.pup.sex8, aes(x=as.character(Sex), 
                                              y=as.numeric(CF_Ratio), 
                                              fill=as.factor(Sex))) +
  geom_boxplot() +
  labs(x=" ", y="Ratio", fill="Sex") +
  labs(title = "Clostridia:Other Firmicutes Ratio for 8 Day Pups") +
  theme(plot.title = element_text(hjust = 0.5))

day8_Sex_CF_plot



# 2-factor ANOVA - 8 day

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


## Day 22

## Processing data for PAE box plot - Day 8

# filter the day

PAE_filtered_ratio_day22 = subset(mot_filtered, Age=="22")


# making datafile to add 'other firmicutes' group


mot.group22 = ddply(PAE_filtered_ratio_day22, c("Group"),
                   summarise,
                   N=length(relative_abundance),
                   total = sum(relative_abundance))


mot.sum22 = ddply(PAE_filtered_ratio_day22, c("Group", "plot_names"),
                 summarise,
                 N=length(relative_abundance),
                 sum = sum(relative_abundance))

samplegroups22 = unique(mot.sum22$Group)

sorted22 = mot.sum22[order(-mot.sum22$sum),]

firm.df22 = NULL

# start loop

for(i in samplegroups22) {
  for(j in i) {
    ## subset dataframe by samples
    sample22 = subset(sorted22, sorted22$Group %in% c(j))
    
    ## get other firmicutes
    firm22 = sample22[c(1:2),]
    
    ## save list of firmicutes abundance taxa
    t.tmp22 <- firm22
    firm.df22 <- rbind.fill(firm.df22, t.tmp22)
    
    ## close loop 
  }
}

firm.df22$place = "keep"

mot.firm22 = full_join(PAE_filtered_ratio_day22, firm.df22)

mot.firm22$place = replace(mot.firm22$place, is.na(mot.firm22$place), "firm")

mot.firm22[mot.firm22$place == "firm",]$plot_names <- "Other Firmicutes"


# making datafile to create ratios

mot.group.pup22 = ddply(mot.firm22, c("Group","Pup_number"),
                       summarise,
                       N=length(relative_abundance),
                       C_Sum = sum(ifelse(plot_names == 
                                            "Firmicutes; Clostridia", 
                                          relative_abundance, 0)),
                       B_Sum = sum(ifelse(plot_names == 
                                            "Bacteroidota; Bacteroidia", 
                                          relative_abundance, 0)),
                       F_Sum = sum(ifelse(plot_names == 
                                            "Other Firmicutes", 
                                          relative_abundance, 0)),
                       CB_Ratio = C_Sum/B_Sum,
                       BF_Ratio = B_Sum/(F_Sum+C_Sum),
                       CF_Ratio = C_Sum/F_Sum)


## Making PAE boxplots - Day 22

# Bacteroidota:Clostridia (CB Ratio)

day22_PAE_CB_plot <-ggplot(mot.group.pup22, aes(x=as.character(Group), 
                                              y=as.numeric(CB_Ratio), 
                                              fill=as.factor(Group))) +
  geom_boxplot() +
  labs(x=" ", y="Ratio", fill="Group") +
  labs(title = "Bacteroidota:Clostridia Ratio for 22 Day Pups") +
  theme(plot.title = element_text(hjust = 0.5))

day22_PAE_CB_plot

# Bacteroidota:All Firmicutes (BF Ratio)

day22_PAE_BF_plot <-ggplot(mot.group.pup22, aes(x=as.character(Group), 
                                              y=as.numeric(BF_Ratio), 
                                              fill=as.factor(Group))) +
  ylim(0, 18) +
  geom_boxplot() +
  labs(x=" ", y="Ratio", fill="Group") +
  labs(title = "Bacteroidota:All Firmicutes Ratio for 22 Day Pups") +
  theme(plot.title = element_text(hjust = 0.5))

day22_PAE_BF_plot

# Clostridia:Other Firmicutes  (CF Ratio)

day22_PAE_CF_plot <-ggplot(mot.group.pup22, aes(x=as.character(Group), 
                                              y=as.numeric(CF_Ratio), 
                                              fill=as.factor(Group))) +
  ylim(0, 22) +
  geom_boxplot() +
  labs(x=" ", y="Ratio", fill="Group") +
  labs(title = "Clostridia:Other Firmicutes Ratio for 22 Day Pups") +
  theme(plot.title = element_text(hjust = 0.5))

day22_PAE_CF_plot



## Processing data for Sex box plot - Day 22

# making datafile to add 'other firmicutes' group


mot.group.sex22 = ddply(PAE_filtered_ratio_day22, c("Sex"),
                       summarise,
                       N=length(relative_abundance),
                       total = sum(relative_abundance))


mot.sum.sex22 = ddply(PAE_filtered_ratio_day22, c("Sex", "plot_names"),
                     summarise,
                     N=length(relative_abundance),
                     sum = sum(relative_abundance))

samplegroups.sex22 = unique(mot.sum.sex22$Sex)

sorted.sex22 = mot.sum.sex22[order(-mot.sum.sex22$sum),]

firm.df.sex22 = NULL

# start loop

for(i in samplegroups.sex22) {
  for(j in i) {
    ## subset dataframe by samples
    sample.sex22 = subset(sorted.sex22, sorted.sex22$Sex %in% c(j))
    
    ## get other firmicutes
    firm.sex22 = sample.sex22[c(1:2),]
    
    ## save list of firmicutes abundance taxa
    t.tmp.sex22 <- firm.sex22
    firm.df.sex22 <- rbind.fill(firm.df.sex22, t.tmp.sex22)
    
    ## close loop 
  }
}

firm.df.sex22$place = "keep"

mot.firm.sex22 = full_join(PAE_filtered_ratio_day22, firm.df.sex22)

mot.firm.sex22$place = replace(mot.firm.sex22$place, 
                               is.na(mot.firm.sex22$place), "firm")

mot.firm.sex22[mot.firm.sex22$place == "firm",]$plot_names <- "Other Firmicutes"


# making datafile to create ratios

mot.group.pup.sex22 = ddply(mot.firm.sex22, c("Sex","Pup_number"),
                           summarise,
                           N=length(relative_abundance),
                           C_Sum = sum(ifelse(plot_names == 
                                                "Firmicutes; Clostridia", 
                                              relative_abundance, 0)),
                           B_Sum = sum(ifelse(plot_names == 
                                                "Bacteroidota; Bacteroidia", 
                                              relative_abundance, 0)),
                           F_Sum = sum(ifelse(plot_names == 
                                                "Other Firmicutes", 
                                              relative_abundance, 0)),
                           CB_Ratio = C_Sum/B_Sum,
                           BF_Ratio = B_Sum/(F_Sum+C_Sum),
                           CF_Ratio = C_Sum/F_Sum)


## Making Sex Boxplots - Day 22

# Bacteroidota:Clostridia (CB Ratio)

day22_Sex_CB_plot <-ggplot(mot.group.pup.sex22, aes(x=as.character(Sex), 
                                                  y=as.numeric(CB_Ratio), 
                                                  fill=as.factor(Sex))) +
  geom_boxplot() +
  labs(x=" ", y="Ratio", fill="Sex") +
  labs(title = "Bacteroidota:Clostridia Ratio for 22 Day Pups") +
  theme(plot.title = element_text(hjust = 0.5))

day22_Sex_CB_plot

# Bacteroidota:All Firmicutes (BF Ratio)

day22_Sex_BF_plot <-ggplot(mot.group.pup.sex22, aes(x=as.character(Sex), 
                                                  y=as.numeric(BF_Ratio), 
                                                  fill=as.factor(Sex))) +
  ylim(0, 18) +
  geom_boxplot() +
  labs(x=" ", y="Ratio", fill="Sex") +
  labs(title = "Bacteroidota:All Firmicutes Ratio for 22 Day Pups") +
  theme(plot.title = element_text(hjust = 0.5))

day22_Sex_BF_plot

# Clostridia:Other Firmicutes  (CF Ratio)

day22_Sex_CF_plot <-ggplot(mot.group.pup.sex22, aes(x=as.character(Sex), 
                                                  y=as.numeric(CF_Ratio), 
                                                  fill=as.factor(Sex))) +
  ylim(0, 35) +
  geom_boxplot() +
  labs(x=" ", y="Ratio", fill="Sex") +
  labs(title = "Clostridia:Other Firmicutes Ratio for 22 Day Pups") +
  theme(plot.title = element_text(hjust = 0.5))

day22_Sex_CF_plot



# 2-factor ANOVA - 22 day

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

## Day 38

## Processing data for PAE box plot - Day 8

# filter the day

PAE_filtered_ratio_day38 = subset(mot_filtered, Age=="38")


# making datafile to add 'other firmicutes' group


mot.group38 = ddply(PAE_filtered_ratio_day38, c("Group"),
                   summarise,
                   N=length(relative_abundance),
                   total = sum(relative_abundance))


mot.sum38 = ddply(PAE_filtered_ratio_day38, c("Group", "plot_names"),
                 summarise,
                 N=length(relative_abundance),
                 sum = sum(relative_abundance))

samplegroups38 = unique(mot.sum38$Group)

sorted38 = mot.sum38[order(-mot.sum38$sum),]

firm.df38 = NULL

# start loop

for(i in samplegroups38) {
  for(j in i) {
    ## subset dataframe by samples
    sample38 = subset(sorted38, sorted38$Group %in% c(j))
    
    ## get other firmicutes
    firm38 = sample38[c(1:2),]
    
    ## save list of firmicutes abundance taxa
    t.tmp38 <- firm38
    firm.df38 <- rbind.fill(firm.df38, t.tmp38)
    
    ## close loop 
  }
}

firm.df38$place = "keep"

mot.firm38 = full_join(PAE_filtered_ratio_day38, firm.df38)

mot.firm38$place = replace(mot.firm38$place, is.na(mot.firm38$place), "firm")

mot.firm38[mot.firm38$place == "firm",]$plot_names <- "Other Firmicutes"


# making datafile to create ratios

mot.group.pup38 = ddply(mot.firm38, c("Group","Pup_number"),
                       summarise,
                       N=length(relative_abundance),
                       C_Sum = sum(ifelse(plot_names == 
                                            "Firmicutes; Clostridia", 
                                          relative_abundance, 0)),
                       B_Sum = sum(ifelse(plot_names == 
                                            "Bacteroidota; Bacteroidia", 
                                          relative_abundance, 0)),
                       F_Sum = sum(ifelse(plot_names == 
                                            "Other Firmicutes", 
                                          relative_abundance, 0)),
                       CB_Ratio = C_Sum/B_Sum,
                       BF_Ratio = B_Sum/(F_Sum+C_Sum),
                       CF_Ratio = C_Sum/F_Sum)


## Making PAE boxplots - Day 38

# Bacteroidota:Clostridia (CB Ratio)

day38_PAE_CB_plot <-ggplot(mot.group.pup38, aes(x=as.character(Group), 
                                              y=as.numeric(CB_Ratio), 
                                              fill=as.factor(Group))) +
  geom_boxplot() +
  labs(x=" ", y="Ratio", fill="Group") +
  labs(title = "Bacteroidota:Clostridia Ratio for 38 Day Pups") +
  theme(plot.title = element_text(hjust = 0.5))

day38_PAE_CB_plot

# Bacteroidota:All Firmicutes (BF Ratio)

day38_PAE_BF_plot <-ggplot(mot.group.pup38, aes(x=as.character(Group), 
                                              y=as.numeric(BF_Ratio), 
                                              fill=as.factor(Group))) +
  ylim(0, 65) +
  geom_boxplot() +
  labs(x=" ", y="Ratio", fill="Group") +
  labs(title = "Bacteroidota:All Firmicutes Ratio for 38 Day Pups") +
  theme(plot.title = element_text(hjust = 0.5))

day38_PAE_BF_plot

# Clostridia:Other Firmicutes  (CF Ratio)

day38_PAE_CF_plot <-ggplot(mot.group.pup38, aes(x=as.character(Group), 
                                              y=as.numeric(CF_Ratio), 
                                              fill=as.factor(Group))) +
  geom_boxplot() +
  labs(x=" ", y="Ratio", fill="Group") +
  labs(title = "Clostridia:Other Firmicutes Ratio for 38 Day Pups") +
  theme(plot.title = element_text(hjust = 0.5))

day38_PAE_CF_plot



## Processing data for Sex box plot - Day 38

# making datafile to add 'other firmicutes' group


mot.group.sex38 = ddply(PAE_filtered_ratio_day38, c("Sex"),
                       summarise,
                       N=length(relative_abundance),
                       total = sum(relative_abundance))


mot.sum.sex38 = ddply(PAE_filtered_ratio_day38, c("Sex", "plot_names"),
                     summarise,
                     N=length(relative_abundance),
                     sum = sum(relative_abundance))

samplegroups.sex38 = unique(mot.sum.sex38$Sex)

sorted.sex38 = mot.sum.sex38[order(-mot.sum.sex38$sum),]

firm.df.sex38 = NULL

# start loop

for(i in samplegroups.sex38) {
  for(j in i) {
    ## subset dataframe by samples
    sample.sex38 = subset(sorted.sex38, sorted.sex38$Sex %in% c(j))
    
    ## get other firmicutes
    firm.sex38 = sample.sex38[c(1:2),]
    
    ## save list of firmicutes abundance taxa
    t.tmp.sex38 <- firm.sex38
    firm.df.sex38 <- rbind.fill(firm.df.sex38, t.tmp.sex38)
    
    ## close loop 
  }
}

firm.df.sex38$place = "keep"

mot.firm.sex38 = full_join(PAE_filtered_ratio_day38, firm.df.sex38)

mot.firm.sex38$place = replace(mot.firm.sex38$place, 
                               is.na(mot.firm.sex38$place), "firm")

mot.firm.sex38[mot.firm.sex38$place == "firm",]$plot_names <- "Other Firmicutes"


# making datafile to create ratios

mot.group.pup.sex38 = ddply(mot.firm.sex38, c("Sex","Pup_number"),
                           summarise,
                           N=length(relative_abundance),
                           C_Sum = sum(ifelse(plot_names == 
                                                "Firmicutes; Clostridia", 
                                              relative_abundance, 0)),
                           B_Sum = sum(ifelse(plot_names == 
                                                "Bacteroidota; Bacteroidia", 
                                              relative_abundance, 0)),
                           F_Sum = sum(ifelse(plot_names == 
                                                "Other Firmicutes", 
                                              relative_abundance, 0)),
                           CB_Ratio = C_Sum/B_Sum,
                           BF_Ratio = B_Sum/(F_Sum+C_Sum),
                           CF_Ratio = C_Sum/F_Sum)


## Making Sex Boxplots - Day 38

# Bacteroidota:Clostridia (CB Ratio)

day38_Sex_CB_plot <-ggplot(mot.group.pup.sex38, aes(x=as.character(Sex), 
                                                  y=as.numeric(CB_Ratio), 
                                                  fill=as.factor(Sex))) +
  geom_boxplot() +
  labs(x=" ", y="Ratio", fill="Sex") +
  labs(title = "Bacteroidota:Clostridia Ratio for 38 Day Pups") +
  theme(plot.title = element_text(hjust = 0.5))

day38_Sex_CB_plot

# Bacteroidota:All Firmicutes (BF Ratio)

day38_Sex_BF_plot <-ggplot(mot.group.pup.sex38, aes(x=as.character(Sex), 
                                                  y=as.numeric(BF_Ratio), 
                                                  fill=as.factor(Sex))) +
  ylim(0, 65) +
  geom_boxplot() +
  labs(x=" ", y="Ratio", fill="Sex") +
  labs(title = "Bacteroidota:All Firmicutes Ratio for 38 Day Pups") +
  theme(plot.title = element_text(hjust = 0.5))

day38_Sex_BF_plot

# Clostridia:Other Firmicutes  (CF Ratio)

day38_Sex_CF_plot <-ggplot(mot.group.pup.sex38, aes(x=as.character(Sex), 
                                                  y=as.numeric(CF_Ratio), 
                                                  fill=as.factor(Sex))) +
  geom_boxplot() +
  labs(x=" ", y="Ratio", fill="Sex") +
  labs(title = "Clostridia:Other Firmicutes Ratio for 38 Day Pups") +
  theme(plot.title = element_text(hjust = 0.5))

day38_Sex_CF_plot

# 2-factor ANOVA - 38 day 

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

