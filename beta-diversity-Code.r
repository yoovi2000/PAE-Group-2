
### Loading packages ###

library(phyloseq)
library(tidyverse)
library(plyr)
library(vegan)
library(stats)
library(devtools)
library(pairwiseAdonis)

### Setting working directory ###

setwd("C:/Users/admin/Downloads")

getwd()


### Viewing Ranks ###

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


## getting the PAE data out of phyloseq

PAE_unfiltered = readRDS("Feb_12_2023_Developmental_phyloseq_P8_P22_P38_unfiltered.RDS")

mot = dephyloseq(PAE_unfiltered)
##turning age column into a numeric ##
mot$Age_numeric <- as.numeric(sub("P", "", mot$Age))

mot$relative_abundance = as.numeric(mot$asv_abundance)/as.numeric(mot$Age_numeric)

mot$plot_names = paste0(mot$Rank4, "; ", mot$Rank6)

view(mot)



### View the ASVs ###

## View the last few rows of the data

metadata = as.data.frame(as.matrix(PAE_unfiltered@sam_data))
otu = as.data.frame(as.matrix(PAE_unfiltered@otu_table))

samdata = as.data.frame(as.matrix(PAE_unfiltered@sam_data))


## View the last few rows of the data

metadata2 = as.data.frame(as.matrix(PAE_unfiltered@sam_data))
otu2 = as.data.frame(as.matrix(PAE_unfiltered@otu_table))

samdata2 = as.data.frame(as.matrix(PAE_unfiltered@sam_data))



### Filtering Data ###

## FILTERING - REMOVE OFF-TAGET TAXA AND MAJOR CONTAMINATION FROM TAXONOMY FILE

df = subset_taxa(PAE_unfiltered, Rank1!="Eukaryota" &
                                 Rank1!="Unassigned")

                   
## FILTERING - REMOVE SAMPLES WITH LESS THAN N READS

sample_sums(df)
plot(sort(sample_sums(df)))                 

df@sam_data$read_depth_noofftargets = sample_sums(df)

which(df@sam_data$read_depth_noofftargets < 1000) 

df.pruned <- prune_samples(sample_sums(df) >= 1000, df)

df.below1000 <- prune_samples(sample_sums(df) < 1000, df)
df.below1000 = as.matrix(df.below1000@sam_data)
write_rds(df.below1000, "PAE_samples_less_than_1000.csv")


## REMOVE INDIVIDUAL ASVS WITH LESS THAN N READS

otu.pruned <- as.data.frame(as.matrix(otu_table(df.pruned)))

otu.pruned$rowsum = rowSums(otu.pruned)

otu.pruned = subset(otu.pruned, otu.pruned$rowsum>100)

otu = subset(otu.pruned, select=-c(rowsum))


## FILTERING - REMOVE ASVs FOUND IN N SAMPLES OR LESS ##

richness = function(x){return(sum(x>0))}

otu$richness = apply(otu,1,richness)
summary(otu$richness)

otu = subset(otu, otu$richness> 3)

summary(otu$richness)

otu = subset(otu, select=-c(richness))


## DENOISING - MAKE ALL THE CELLS IN THE OTU TABLE WITH VALUES N OR LESS 0 ##

otu <- mutate_all(otu, funs(ifelse(. < 8, 0, .)))


## CREATE AND READ BACK IN FILTERED BUT NOT RAREFIED PHYLOSEQ OBJECT

otu_mat = as.matrix(otu)

OTU = otu_table(otu_mat, taxa_are_rows = TRUE)

tax_mat = df.pruned@tax_table
TAX = tax_table(tax_mat)  

samples = as.data.frame(as.matrix(df.pruned@sam_data))

samples = sample_data(samples)

df.prerarefaction = phyloseq(OTU, TAX, samples)

df.prerarefaction
rank_names(df.prerarefaction)
sample_variables(df.prerarefaction)


## FINAL DENOISING AND SAVE FILTERED DATA

df.prerarefaction = prune_taxa(taxa_sums(df.prerarefaction)>0, df.prerarefaction)

df.prerarefaction@sam_data$read_depth_filtered = sample_sums(df.prerarefaction)

histogram(df.prerarefaction@sam_data$read_depth_filtered, breaks=100)
# df.prerarefaction = subset_samples(df.prerarefaction, read_depth_filtered<N)
# Not necessary, only goes up to 25 000
histogram(df.prerarefaction@sam_data$read_depth_filtered, breaks=100)

write_rds(df.prerarefaction, "PAE_filtered.RDS")


## SIMPLE RAREFACTION - NORMALIZING THE DATA

plot(sort(sample_sums(df.prerarefaction)))
summary(sample_sums(df.prerarefaction))

sort(sample_sums(df.prerarefaction))

View(df.prerarefaction@sam_data)

# run "View(df.prerarefaction@sam_data)", then choose the second value from the top
# and set the cutoff to less than than number (289-1=288) to include it.
# later, we standardize to this number.

rare.threshold = which(sample_sums(df.prerarefaction) <288) 
rare.threshold

samples.lost <- prune_samples(sample_sums(df.prerarefaction) <= 288, df.prerarefaction)
metadata.lost = as.data.frame(samples.lost@sam_data)
write.csv(metadata.lost, "PAE_samples_below_rare_threshold.csv")

require(stats)

set.seed(5)

sample(10)
sample(10)

df.rarefied <- rarefy_even_depth(df.prerarefaction, sample.size = 288) 


## CHECK THAT RAREFACTION WORKED

abundance = function(x){
  return(sum(x,na.rm=TRUE))}

otu = as.data.frame(df.rarefied@otu_table)

otu=as.matrix(otu)

otu = t(otu)
otu = as.data.frame(otu)

otu$abundance = apply(otu,1,abundance)

summary(otu$abundance) 

write_rds(df.rarefied, "Rarefied_PAE_Dataset.RDS")


### END OF FILTERING ###




## set your seed for repoducibility 
set.seed(2)

## read in the phyloseq object 
Pae_filt = readRDS("Rarefied_PAE_Dataset.RDS")
##### RUN A PERMANOVA #####

## betadispersion test - see if your within group variation is the same
#set up for betadispersion test
project_bray.rarefied <- phyloseq::distance(Pae_filt, method = "bray")
sample_df <- data.frame(sample_data(Pae_filt))

beta.FACTOR1 <- betadisper(project_jaccard.rarefied, sample_df$Group) 
b1.1 = permutest(beta.FACTOR1) 
b1.1

# get your data out of phyloseq 
metadata = as.data.frame(as.matrix(Pae_filt@sam_data))
otu = as.data.frame(as.matrix(Pae_filt@otu_table))

# merge the metadata and otu table together
## by can be a column but b=0 joins by rownames. It's very handy 
metaotu = merge(metadata, otu, by=0)

## count the number of metadata columns you have. Do +1 because the rownames become a columns with merge
metacols = ncol(metadata)+1

##Run the PERMANOVA
p1 = adonis2(metaotu[,-c(1:metacols)]~Group, data=metaotu)
print(p1)

## Conduct a post-hoc test
pair.id = pairwise.adonis(metaotu[,-c(1:metacols)], metaotu$Group)
print(pair.id)

## this printed output sucks. get something better
sig.pairs = subset(pair.id, pair.id$p.adjusted<=0.1)
sig.pairs


##### MAKE AND NMDS PLOT #####
## ordinate your data to make the NMDS plot 
Pae_filt.ord = ordinate(Pae_filt, "NMDS","jaccard")

## see what the fit looks like
stressplot(Pae_filt.ord)

## See what the stress is for this solution 
## stress below 2 is considered low 
Pae_filt.ord

## make the NMDS plot
plot_ordination(Pae_filt, Pae_filt.ord, color="Group")+
  geom_point(size=3)+
  annotate("text", x = 1, y = 1.1, label = "Bray", size=4, color="black")+
  annotate("text", x = 1, y = 1.2, label = "stress = 0.178211", size=4, color="black")+
  scale_color_manual(values=c("navy","dodgerblue4","dodgerblue2", "skyblue4","skyblue1", "cyan3","cyan",
                              "forestgreen","green2","purple3","mediumorchid","pink4","pink1"))+
  theme_bw()



sample_df$Group <- ifelse(sample_df$Group == "Control", 0, 1) 
sample_df$Age <- as.numeric(gsub("P", "", sample_df$Age))
# Subset 1: Age=P38, Group=0
subset1 <- subset(sample_df, Age == 38 & Group == 0)
subset1[, 3:ncol(subset1)] <- lapply(subset1[, 3:ncol(subset1)], as.numeric)

# Subset 2: Age=P8, Group=1
subset2 <- subset(sample_df, Age == 8 & Group == 1)
subset2[, 3:ncol(subset2)] <- lapply(subset2[, 3:ncol(subset2)], as.numeric)
dist_mat <- vegdist(rbind(subset1[, -c(1,2)], subset2[, -c(1,2)]), method = "bray")
dist_mat
# Perform PERMANOVA test
res <- adonis(dist_mat ~ subset1$Age, permutations = 999)
res