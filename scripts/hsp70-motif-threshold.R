## Sammy Keyport
## June 27th, 2018
## Hsp70 Motif Threshold
## Purpose: Continuously plot the presence of an Hsp70 binding motif with increasing calling score

library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)
library(readr)

#read in al files in the folder that are output for motif finder code
files <- list.files(path="C:/Users/brandon/Desktop/research/hsp70-binding/data", pattern="*r10.txt", 
                    full.names=T, recursive=FALSE)
#set above to your path into hsp70-binding/data folder which contains all of the motif files


# function that takes the motif file, gets the name and the top score, and outputs it as a list
top_score <- function(motifs) {
  filename <- strsplit(motifs, "/")
  proteinfile <- strsplit(filename[[1]][8], "-")
  protein <- proteinfile[[1]][1]
  data <- read_tsv(motifs, comment="#")
  max.score <- max(data$score, NA.rm = TRUE)
  return(c(protein, as.numeric(max.score)))
}

# run top score function for all files, and change it into a data frame, give it nice column names
top.scores.list <- lapply(files, top_score)
top.scores.df <- as.data.frame(do.call(rbind, top.scores.list))
top.scores.df <- mutate_all(top.scores.df, funs(tolower))
names(top.scores.df) <- c("Name", "Score")

# load in file containing info about aggregation rate of protein of interest (here, just GTPases)
gtpase <- as.data.frame(read_tsv("../../hsp70-binding/analysis/gtpase_list.txt", col_names = F))
gtpase <- mutate_all(gtpase, funs(tolower))
names(gtpase) <- "Name"

#for all aggregating proteins, load in file and make sure the formatting is correct
all.agg.rate <- as.data.frame(read_tsv("../../heat_aggregators.txt", comment = "#", col_names = T))
# for calling aggregation, using psup_dryad file from 2015 paper with slight modifications. 
# Using a delta p value of 0.3 (was 0.4) at a difference from Time 0 to Time 4 
# (instead of from Time 0 to Time 2)

# give it names, and change to agg (aggregator) or superagg (super aggregator)
names(all.agg.rate) <- (c("ORF", "Name", "Category"))
all.agg.rate$Category[all.agg.rate$Category == "FALSE"] <- "agg"
all.agg.rate$Category[all.agg.rate$Category == "TRUE"] <- "superagg"

#need to make it all lowercase so downstream data frame joining will work
all.agg.rate <- mutate_all(all.agg.rate, funs(tolower))
all.agg.rate <- subset(all.agg.rate, select = c("Name", "Category"))

# load in file containing the names of proteins associated with ribosome biogenesis
ribo.bio <- as.data.frame(read_tsv("../analysis/gene_list_ribo_bio.txt", col_names = FALSE))
ribo.bio <- mutate_all(ribo.bio, funs(tolower))
names(ribo.bio) <- "Name"

# load in file containing the names of all proteins in the proteome
proteome <- as.data.frame(read_tsv("../analysis/yeast_names.txt", col_names = FALSE))
proteome <- mutate_all(proteome, funs(tolower)) # make lowercase
names(proteome) <- "Name"

# now, make new dfs of proteins from each group that aggregate
ribo.bio.agg <- inner_join(ribo.bio, all.agg.rate, by = "Name")
gtpase.agg <- inner_join(gtpase, all.agg.rate, by = "Name")
proteome.agg <- inner_join(proteome, all.agg.rate, by = "Name")

# join top score data frame with data frame separated by protein of interest by the protein name
all.agg.list <- inner_join(top.scores.df, all.agg.rate, by = "Name")
ribo.bio.list <- inner_join(top.scores.df, ribo.bio, by = "Name")
ribo.bio.agg <- inner_join(top.scores.df, ribo.bio.agg, by = "Name")
gtpase.list <- inner_join(top.scores.df, gtpase, by = "Name")
agg.gtpase.list <- inner_join(top.scores.df, gtpase.agg, by = "Name")
proteome.list <- inner_join(top.scores.df, proteome, by = "Name")
agg.proteome <- inner_join(top.scores.df, proteome.agg, by = "Name")

# Now we want to use a varying threshold - need the score column to do this
# initialize data frame and set sequence for varying the threshold
max.score.df <- data.frame(NAME = c(), SCORE = c())
in_list <- seq(from = 0, to = 20, by = 0.1)

# Function that takes a threshold value (will be from previously defined "in_list" - see above) 
# and a sub.frame that will be defined by aggregation category if applicable.
# Find what proportion of proteins in a group have a max score above the given threshold

Vary_Threshold <- function(threshold, sub.frame) {
  filter_data <- filter(sub.frame, as.numeric(Score) >= threshold)
  return(nrow(filter_data))
}

# run Vary_Threshold for each aggregation category (including all categories)
all.agg <- lapply(in_list, Vary_Threshold, sub.frame = 
                filter(all.agg.list, Category %in% c("agg", "superagg")))
all.super.agg <- lapply(in_list, Vary_Threshold, sub.frame = 
                          filter(all.agg.list, Category == "superagg"))
ribo.bio <- lapply(in_list, Vary_Threshold, sub.frame = ribo.bio.list)
agg.ribo.bio <- lapply(in_list, Vary_Threshold, sub.frame = ribo.bio.agg)
gtpase.all <- lapply(in_list, Vary_Threshold, sub.frame = gtpase.list)
gtpase.agg <- lapply(in_list, Vary_Threshold, sub.frame = 
                       filter(agg.gtpase.list, Category %in% c("agg", "superagg")))
gtpase.super.agg <- lapply(in_list, Vary_Threshold, sub.frame = 
                       filter(agg.gtpase.list, Category == "superagg"))
all.proteome <- lapply(in_list, Vary_Threshold, sub.frame = proteome.list)
proteome.agg <- lapply(in_list, Vary_Threshold, sub.frame = 
                         filter(agg.proteome, Category %in% c("agg", "superagg")))
proteome.super.agg <- lapply(in_list, Vary_Threshold, sub.frame = 
                         filter(agg.proteome, Category == "superagg"))

# change each list into a data frame 
agg.df <- as.data.frame(do.call(rbind, agg))
super.agg.df <- as.data.frame(do.call(rbind, super.agg))
ribo.bio.df <- as.data.frame(do.call(rbind, ribo.bio))
agg.ribo.bio.df <- as.data.frame(do.call(rbind, agg.ribo.bio))
gtpase.all.df <- as.data.frame(do.call(rbind, gtpase.all))
gtpase.agg.df <- as.data.frame(do.call(rbind, gtpase.agg))
gtpase.super.agg.df <- as.data.frame(do.call(rbind, gtpase.super.agg))
proteome.df <- as.data.frame(do.call(rbind, all.proteome))


# make it into a nice table before plotting
comp_table <- cbind(in_list, agg.df, super.agg.df, ribo.bio.df, agg.ribo.bio.df, 
                    gtpase.all.df, gtpase.agg.df, gtpase.super.agg.df
                    #, proteome.df
                    )

# make table into proportions, then combine back with threshold
props <- comp_table[1,2:length(comp_table)]
shorter <- comp_table[,2:length(comp_table)]
new_table <- mapply(`/`, data.frame(shorter), props)
comp_table2 <- as.data.frame(cbind(comp_table$in_list, new_table))
colnames(comp_table2) <- c("Threshold", 
                           "All Aggregating Proteins", "Super Aggregating Proteins", 
                           "Ribosome Biogenesis Proteins", "Aggregating RiboBio Proteins",
                           "All GTPases", "Aggregating GTPases", "Super Aggregating GTPases",
                           "All Proteins")

# make data into a form that ggplot likes :)
dd = melt(comp_table2, id=c("Threshold"))
dd1 = melt(comp_table, id=c("Threshold"))

# graphing
p <- ggplot(dd) + geom_line(aes(x=Threshold, y=value, colour=variable), size = 2) +
  scale_colour_manual(values=c("red","green","blue", "pink", "orange", "purple", "yellow", "dark green")) +
  ggtitle("Proteins with Hsp70 Motif") + ylab("Proportion with Motif") + 
  xlab("Threshold Score for Calling Motif") + geom_vline(xintercept=12, size = 1) +
  labs(colour = "Protein Type") 
p

# Trying to use stat_ecdf function. Did I do this right?
# Must include everything except "All Proteins" otherwise it looks horrible
comp_table <- cbind(in_list, agg.df, super.agg.df, ribo.bio.df, agg.ribo.bio.df, 
                    gtpase.all.df, gtpase.agg.df, gtpase.super.agg.df, proteome.df)
colnames(comp_table) <- c("Threshold", 
                          "All Aggregating Proteins", "Super Aggregating Proteins", 
                          "Ribosome Biogenesis Proteins", "Aggregating RiboBio Proteins",
                          "All GTPases", "Aggregating GTPases", "Super Aggregating GTPases")
p2 <- ggplot(dd1, aes(x = value, group = variable, color = variable)) + stat_ecdf() +
  ggtitle("Proteins with Hsp70 Motif") + ylab("Proportion with Motif") + 
  xlab("Count") + labs(colour = "Protein Type") 
p2

# find significance with KS test - weird results!
ks.test(comp_table2$All.Aggregating.Proteins, comp_table2$Aggregating.RiboBio.Proteins)
ks.test(comp_table2$Ribosome.Biogenesis.Proteins, comp_table2$Aggregating.RiboBio.Proteins)
ks.test(comp_table2$All.Aggregating.Proteins, comp_table2$Super.Aggregating.GTPases)
ks.test(comp_table2$Aggregating.RiboBio.Proteins, comp_table2$Ribosome.Biogenesis.Proteins)
ks.test(comp_table2$Aggregating.GTPases, comp_table2$Super.Aggregating.GTPases)
ks.test(comp_table2$All.GTPases, comp_table2$Super.Aggregating.GTPases)
ks.test(comp_table2$All.GTPases, comp_table2$Aggregating.GTPases)
ks.test(comp_table2$Super.Aggregating.Proteins, comp_table2$Ribosome.Biogenesis.Proteins)


# Wilcox test is better, but still some weird results. See the p-val when comparing 
# all RiboBio Proteins to just aggregating RiboBio proteins... I expect it to be significant when 
# looking at plot p

wilcox.test(comp_table2$All.Aggregating.Proteins, comp_table2$Aggregating.RiboBio.Proteins)
wilcox.test(comp_table2$Ribosome.Biogenesis.Proteins, comp_table2$Aggregating.RiboBio.Proteins)
wilcox.test(comp_table2$Super.Aggregating.GTPases, comp_table2$All.Proteins)
wilcox.test(comp_table2$Super.Aggregating.GTPases, comp_table2$All.GTPases)
wilcox.test(comp_table2$Super.Aggregating.GTPases, comp_table2$Super.Aggregating.Proteins)
wilcox.test(comp_table2$Super.Aggregating.GTPases, comp_table2$Aggregating.GTPases)
wilcox.test(comp_table2$Super.Aggregating.GTPases, comp_table2$Aggregating.RiboBio.Proteins)
