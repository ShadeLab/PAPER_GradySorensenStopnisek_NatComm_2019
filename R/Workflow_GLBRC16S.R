# Total R Code for GLBRC 
### Start Common
### DIRECTIONS, SET WORKING DIRECTORY TO "PAPER_GradySorensenStopnisek_InPrep/R"
### Before working on Analysis, run the entire Common block of Code (all the way until "### End Common")
library(dplyr)
library(tidyr)
library(reshape2)
library(RSQLite)
library(stringr)
# Read in OTU table
otu <- read.table("InputFiles/table_combined_merged_trimmed_otus.txt",sep="\t", header=TRUE, stringsAsFactors = FALSE, row.names=1)

#Preping the environmental metadata
glbrc <- dbConnect(SQLite(), dbname="InputFiles/GLBRC_bioenergy_db.db" )

#Content of the DB
dbListTables(glbrc)

#Content of selected tables in the DB
dbListFields(glbrc, 'plant')
dbListFields(glbrc, 'sequencing')
dbListFields(glbrc, 'soil')
dbListFields(glbrc, 'nucleic_acids')

#get tables into R
glbrc_NA <- dbGetQuery(glbrc, "select * from nucleic_acids") 
glbrc_soil <- dbGetQuery(glbrc, 'select * from soil')
glbrc_plant <- dbGetQuery(glbrc, 'select * from plant')
glbrc_plot <- dbGetQuery(glbrc, 'select * from plot')
glbrc_sampling <- dbGetQuery(glbrc, 'select * from sampling')
glbrc_sequncing <- dbGetQuery(glbrc, 'select * from sequencing')
glbrc_sequncing <- glbrc_sequncing %>% 
  mutate(nucleic_acid_name = str_trim(glbrc_sequncing$nucleic_acid_name, side = "both"))

#joining tables to create complete map file
metadata <- full_join(glbrc_plot, glbrc_sampling, by='plotID')
metadata <- full_join(metadata, glbrc_soil, by='sampleID')
metadata <- full_join(metadata, glbrc_plant, by='sampleID')
metadata <- full_join(metadata, glbrc_NA, by='sampleID')
metadata <- full_join(metadata, glbrc_sequncing, by='nucleic_acid_name')

unique(metadata$nucleic_acid_name)
map <- metadata
dim(map)

map_16 <- subset(map, map$primers == 'EMP V4')
dim(map_16)
head(map_16)
map_16$help_name = as.character(lapply(strsplit(as.character(map_16$nucleic_acid_name), split="D"), "[", 1))
unique(map_16$help_name)

#finding duplicates
n_occur <- data.frame(table(map_16$help_name))
n_occur[n_occur$Freq > 1,]
duplicate_df <- map_16[map_16$help_name %in% n_occur$Var1[n_occur$Freq > 1],]
list_dupli <- duplicate_df$sequence_name #list of duplicate samples

D1 <- duplicate_df[grep('D1', duplicate_df$sequence_name),]
D1$removing <- 'remove'
dim(D1)

map_full <- full_join(map, D1) 

#creating numeric time column
map_full$sampling_date <- paste0(map_full$month,'-', map_full$day,'-',map_full$year) 
map_full$sampling_date <- as.POSIXct(map_full$sampling_date, format='%m-%d-%Y')
time <- as.POSIXct(map_full$sampling_date, format='%m-%d-%Y')
time_numeric <- as.numeric(time)
map_time <- cbind(map_full, time_numeric)

#adding weather data
weather <- read.csv('InputFiles/kbs_weather_09212017.csv', encoding = 'UTF-8', na.strings= c("NA", " ", ""))
dim(weather)
head(weather)
weather$sampling_date <- as.POSIXct(weather$date, format='%d.%m.%y')

sub_weather <- weather[weather$sampling_date %in% map_time$sampling_date,] #subsetting weather file for sample dates

#merging dataframes - map file and weather
map_complete <- full_join(map_time, weather)

#Subset to 16S Samples
map_16S <- subset(map_complete, map_complete$primers=="EMP V4")
# Get the name of the samples we have data for
# Remove any duplicates that we don't want to use
map_small <- map_16S[map_16S$exclude_from_analysis=="N",]

#map_small_v1 <- subset(map_16S_small, map_16S_small$removing == 'remove') 
#map_small <- map_16S_small[!(rownames(map_16S_small) %in% c(33,108, 185,259,37,418,467,487,571)),]

samples <- colnames(otu)
# put taxonomy into its own variable
taxonomy <- read.csv('InputFiles/taxonomy_combined_merged_trimmed_otus.csv', header = T, row.names = 1, na.strings= c("NA", " ", ""))
# subset the map to include only those samples we have sequence data
map_small <- map_small[map_small$sequence_name %in% samples,]


# Remove rows that have N/A for the sequence name
map_small<- map_small[complete.cases(map_small$sequence_name),]

# Subset the samples to those we want to analyse (IE remove the duplicates)
samples <- samples[samples %in% map_small$sequence_name]
#Subset the OTU table to only the samples we want to analyze(IE remove the duplicates)
otu_sub <- otu[,colnames(otu) %in% samples]
# Order the samples
otu_sub <- otu_sub[,order(colnames(otu_sub))]
# Order the samples of the map the same way
map_small <- map_small[order(map_small$sequence_name),]
# Check to make sure they all match with each other
colnames(otu_sub) == map_small$sequence_name

# Map file that has duplicates of samples removed
map_16S <- map_small
# OTU table that removes duplicates of single samples
otu <- otu_sub
otu_CM <- otu
#removing Eukaryota from the OTU table
tax_short <- taxonomy[!grepl("Mitochondria", taxonomy$Family),]
tax_short <- tax_short[!grepl("Chloroplast", tax_short$Class),]

otu <- otu[rownames(otu) %in% rownames(tax_short),]

taxonomy_full <- taxonomy
taxonomy <- taxonomy[rowSums(otu_sub)>0,]
otu <- otu[rowSums(otu)>0,]
otu_soil <- otu[,map_16S$source=="soil"]

library(vegan)
set.seed(13)
otu_rare <- t(rrarefy(t(otu), 1000))
otu_soil_rare <- t(rrarefy(t(otu_soil), min(colSums(otu_soil))))
otu_rare <- otu_rare[,colSums(otu_rare)>999]
map_16S <- map_16S[map_16S$sequence_name%in%colnames(otu_rare),]

rarecurve(t(otu_rare), step=20, sample = min(colSums(otu_rare)), label = FALSE)

otu_rare <- otu_rare[,colSums(otu_rare)>999]

# ### Multiplot code taken from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


### End Common

### Start Nejc Analysis
setwd('Documents/git/phyllosphere/Analyses')

library(ggplot2)
library(dplyr)
library(tidyr)
library(vegan)
source("http://bioconductor.org/biocLite.R")
biocLite("Heatplus")
biocLite('limma')
library(Heatplus)
library("limma")
library(ggrepel)
library(codyn)
library(gridExtra)
library(grid)
install.packages('egg')
library(egg)

##################################

otu_col <- rep('darkgreen',length(names(otu)))
otu_col[map_16S$plant == 'switchgrass' & map_16S$source == 'phyllosphere' & map_16S$year==2016] <- 'darkolivegreen3'
otu_col[map_16S$plant == 'switchgrass' & map_16S$source == 'soil' & map_16S$year==2016] <- 'burlywood'
otu_col[map_16S$plant == 'miscanthus' & map_16S$source == 'soil' & map_16S$year==2016] <- 'burlywood4'
otu_col[map_16S$plant == 'switchgrass' & map_16S$source == 'phyllosphere' & map_16S$year==2017] <- 'darkolivegreen1'
otu_col[map_16S$plant == 'switchgrass' & map_16S$source == 'soil' & map_16S$year==2017] <- 'brown'

rarecurve(t(otu), step=1000, label=F, col=otu_col)

#remove empty species
otu_rare=otu_rare[rowSums(otu_rare)>0,]

misc_otu <- otu_rare[,map_16S$plant=="miscanthus" & (map_16S$source=='phyllosphere')]
swit16_otu <- otu_rare[,map_16S$plant=="switchgrass" & (map_16S$source=='phyllosphere') & (map_16S$year=='2016')]
swit17_otu <- otu_rare[,map_16S$plant=="switchgrass" & (map_16S$source=='phyllosphere') & (map_16S$year=='2017')]

misc_otu_soil <- otu[,map_16S$plant=="miscanthus" & (map_16S$source=='soil')]
swit_otu_soil <- otu[,map_16S$plant=="switchgrass" & (map_16S$source=='soil')]

soil_otu <- otu_rare[,map_16S$source=="soil"] #both soils: switchgrass and miscanthus

phyllo <- otu_rare[,map_16S$source=="phyllosphere"]
phyllo_otu <- phyllo[rowSums(phyllo)>0,]
nrow(phyllo_otu)
rownames(phyllo_otu)

# Venn Analysis of Switchgrass and miscanthus
#############################################

#Phyllosphere only
# make presence absence list from soil and plant into 1 & 0
swit16_otu <- 1*(rowSums(swit16_otu)>0)
swit17_otu <- 1*(rowSums(swit17_otu)>0)
misc16_otu <- 1*(rowSums(misc_otu)>0)

venn_phyllo_data <- cbind(swit16_otu,swit17_otu, misc16_otu)
colnames(venn_phyllo_data) <- c("Switchgrass 2016", "Switchgrass 2017", "Miscanthus 2016")
venn_phyllo_data=venn_phyllo_data[rowSums(venn_phyllo_data)>0,]
v_phyllo=vennCounts(venn_phyllo_data)
v_phyllo_2=round(v_phyllo[,"Counts"]/sum(v_phyllo[,"Counts"]),2) #calculate percentage of each group
vennDiagram(v_phyllo, circle.col = c('darkolivegreen3', 'darkolivegreen1','darkgreen'), lwd=6, cex=1.2, scale=F)

#####################################
####Contextual data plotting
#####################################

#Weather report
ggplot(weather[weather$Year == 2016 & weather$sampling_date>'2016-04-01 EDT' & weather$sampling_date<'2016-11-10 EDT',], 
       aes(x=sampling_date, y=precipitation))+
  geom_vline(xintercept = c(unique(map_time$sampling_date)),linetype="dashed") +
  geom_point(aes(y=RH), colour='blue', size=.8) +
  geom_point(aes(y=Air_temp_mean), colour='red', size=.8) +
  geom_point() +
  labs(y='Precipitation (mm) - black;\nMean temperature (°C)- red;\nRelative humidity (%) - blue', x='Date')
ggsave('Figures/weather_report.eps', height=4, width=5)

#Contextual data analysis
treat_colors <- rep("green", length(map_16S))
treat_colors[map_16S$treatment=="standard fertilization"]<-"black"

#Switchgrass
par(mfrow=c(2,2))
plot(map_16S$sampling_date[map_16S$source=='phyllosphere' & map_16S$plant=='switchgrass'], map_16S$height_mean_cm[map_16S$source=='phyllosphere' & map_16S$plant=='switchgrass'], 
     ylab='Plant mean height (cm)',
     xlab=NA,
     main= 'Switchgrass - plant',
     col=treat_colors)
legend('bottomright', col=c('green', 'black'), legend=c('N-free', 'N'), bty='n', pch=21, title='Treatment')
plot(map_16S$sampling_date[map_16S$source=='phyllosphere' & map_16S$plant=='switchgrass'], map_16S$carbon_per_nitrogen[map_16S$source=='phyllosphere' & map_16S$plant=='switchgrass'], 
     ylab='Leaf CN ratio',
     xlab=NA,
     main= 'Switchgrass - plant',
     col=treat_colors)
plot(map_16S$sampling_date[map_16S$source=='soil' & map_16S$plant=='switchgrass'], map_16S$NO3N_ppm[map_16S$source=='soil' & map_16S$plant=='switchgrass'], 
     ylab='Soil NO3-N (ppm)',
     xlab='Sampling date (Month)',
     main= 'Switchgrass - soil',
     col=treat_colors)
plot(map_16S$sampling_date[map_16S$source=='soil' & map_16S$plant=='switchgrass'], map_16S$organic_matter[map_16S$source=='soil' & map_16S$plant=='switchgrass'], 
     ylab='Soil organic matter (ppm)',
     xlab='Sampling date (Month)',
     main= 'Switchgrass - soil',
     col=treat_colors)

#Miscanthus
par(mfrow=c(2,2))
plot(map_16S$sampling_date[map_16S$source=='phyllosphere' & map_16S$plant=='miscanthus'], map_16S$height_mean_cm[map_16S$source=='phyllosphere' & map_16S$plant=='miscanthus'], 
     ylab='Plant mean height (cm)',
     xlab=NA,
     main= 'Miscanthus - plant',
     col=treat_colors)
legend('bottomright', col=c('green', 'black'), legend=c('N-free', 'N'), bty='n', pch=21, title='Treatment')
plot(map_16S$sampling_date[map_16S$source=='phyllosphere' & map_16S$plant=='miscanthus'], map_16S$carbon_per_nitrogen[map_16S$source=='phyllosphere' & map_16S$plant=='miscanthus'], 
     ylab='Leaf CN ratio',
     xlab=NA,
     main= 'Miscanthus - plant',
     col=treat_colors)
plot(map_16S$sampling_date[map_16S$source=='soil' & map_16S$plant=='miscanthus'], map_16S$NO3N_ppm[map_16S$source=='soil' & map_16S$plant=='miscanthus'], 
     ylab='Soil NO3-N (ppm)',
     xlab='Sampling date (Month)',
     main= 'Miscanthus - soil',
     col=treat_colors)
plot(map_16S$sampling_date[map_16S$source=='soil' & map_16S$plant=='miscanthus'], map_16S$organic_matter[map_16S$source=='soil' & map_16S$plant=='miscanthus'], 
     ylab='Soil organic matter (ppm)',
     xlab='Sampling date (Month)',
     main= 'Miscanthus - soil',
     col=treat_colors)

par(mfrow=c(1,1))

# #Using Occupancy data to retrieve OTUs with most frequent/common observations (found in more than 80% of samples)
# #Example for switchgrass
# names(Occ_swit)
# 
# occ_swit_df<- data.frame(otu=names(Occ_swit), occ=Occ_swit) %>% filter(occ > .4)
# abund_otu_switch <- data.frame(otu=names(Mean_abund_swit), abund=log(Mean_abund_swit)) %>% filter(abund > -6)
# switch_occ_abund <- full_join(occ_swit_df, abund_otu_switch)
# 
# occ_misc_df<- data.frame(otu=names(Occ_misc), occ_misc=Occ_misc) %>% filter(occ_misc > .4)
# abund_otu_misc <- data.frame(otu=names(Mean_abund_misc), abund_misc=log(Mean_abund_misc)) %>% filter(abund_misc > -6)
# misc_occ_abund <- full_join(occ_misc_df, abund_otu_misc)
# 
# phyllo_occ_abund <- full_join(switch_occ_abund, misc_occ_abund)
# 
# occ_swit_df$otu <- as.character(occ_swit_df$otu)
# 
# rel_otu_rare <- decostand(otu_rare, method="total", MARGIN=2)
# 
# occ_rel_otu_rare <- rel_otu_rare[rownames(rel_otu_rare) %in% occ_swit_df$otu,]
# occ_rel_otu_rare[5,]
# 
# tax_occ <- taxonomy[as.numeric(occ_swit_df$otu)]

# ########################################################
# #Heatmap using relative abundance data (rarefied counts)
# ########################################################
# #Switchgrass
# swit_otu_rel_presence=swit_otu_rel[rowSums(swit_otu_rel)>0,]
# 
# # determine the maximum relative abundance for each column
# maxab <- apply(swit_otu_rel_presence, 1, max)
# head(maxab)
# 
# # remove the genera with less than .001% as their maximum relative abundance
# n1 <- names(which(maxab < 0.01))
# swit_otu_rel.001 <- swit_otu_rel_presence[-which(rownames(swit_otu_rel_presence) %in% n1),]
# dim(swit_otu_rel_presence)
# 
# scaleyellowred <- colorRampPalette(c("white", "black"), space = "rgb")(100)
# 
# # calculate the Bray-Curtis dissimilarity matrix on the full dataset:
# data.dist_col <- vegdist(t(swit_otu_rel), method = "bray")
# data.dist_row <- vegdist(swit_otu_rel.001, method = "bray")
# # Do average linkage hierarchical clustering 
# cl.clus <- hclust(data.dist_col, "aver")
# row.clus <- hclust(data.dist_row, "aver")
# 
# # make the heatmap
# #pdf("switchgrass_heatmap.pdf", height = 7, width = 7)
# heatmap(as.matrix(swit_otu_rel.001), 
#         Colv = as.dendrogram(cl.clus), 
#         Rowv = as.dendrogram(row.clus), 
#         col = scaleyellowred, 
#         margins = c(10, 3),
#         main='Switchgrass OTU relative anundance'
# )
# #dev.off()
# 
# #Miscanthus
# # determine the maximum relative abundance for each column
# maxab_misc <- apply(misc_otu_rel, 1, max)
# head(maxab_misc)
# 
# # remove the genera with less than .001% as their maximum relative abundance
# n1_misc <- names(which(maxab_misc < 0.01))
# misc_otu_rel.001 <- misc_otu_rel[-which(rownames(misc_otu_rel) %in% n1_misc),]
# dim(misc_otu_rel.001)
# 
# # calculate the Bray-Curtis dissimilarity matrix on the full dataset:
# data.dist_col_misc <- vegdist(t(misc_otu_rel), method = "bray")
# data.dist_row_misc <- vegdist(misc_otu_rel.001, method = "bray")
# 
# # Do average linkage hierarchical clustering 
# cl.clus_misc <- hclust(data.dist_col_misc, method = "complete")
# row.clus_misc <- hclust(data.dist_row_misc, "average")
# 
# # make the heatmap with Rowv = as.dendrogram(row.clus)
# #pdf('miscanthus_heatmap.pdf', width = 7, height = 8)
# heatmap(as.matrix(misc_otu_rel.001), 
#         Colv = as.dendrogram(cl.clus_misc), 
#         Rowv = as.dendrogram(row.clus_misc), 
#         col = scaleyellowred, 
#         margins = c(10, 3),
#         main='Miscanthus OTU relative anundance')
# #dev.off()

# #************************************************************************************************************************************
# #Replotting the Occ_Abund data - using occupancy and relative abundance only. Also, use coler coding for unique OTUs and shared OTUs. 
# #Mark the quadrant used for the analysis - occ>.4 and rel abund > log10>-2.5
# #************************************************************************************************************************************
# 
#####2016 and 2017 data separated
swit_otu <- otu_rare[,map_16S$plant=="switchgrass" & (map_16S$source=='phyllosphere') & (map_16S$year==2016)]
swit_otu_17 <- otu_rare[,map_16S$plant=="switchgrass" & (map_16S$source=='phyllosphere') & (map_16S$year==2017)]

#Comulative Occ_Abund
#2016 data
swit_otu_PA <- 1*((swit_otu>0)==1)
swit_otu_PA <- swit_otu_PA[rowSums(swit_otu_PA)>0,]
Occ_swit <- rowSums(swit16_otu_PA)/ncol(swit_otu_PA)
swit_otu <- swit_otu[rowSums(swit_otu)>0,]
swit_otu_rel <- decostand(swit_otu, method="total", MARGIN=2)
com_abund_swit <- rowSums(swit_otu_rel)

#2017 data
swit17_otu_PA <- 1*((swit_otu_17>0)==1)
swit17_otu_PA <- swit17_otu_PA[rowSums(swit17_otu_PA)>0,]
Occ_swit17 <- rowSums(swit17_otu_PA)/ncol(swit17_otu_PA)
swit17_otu <- swit_otu_17[rowSums(swit_otu_17)>0,]
swit17_otu_rel <- decostand(swit17_otu, method="total", MARGIN=2)
com_abund_swit17 <- rowSums(swit17_otu_rel)
#nrow(swit17_otu)

#Miscanthus Occ_Abund
misc_otu <- misc_otu[rowSums(misc_otu)>0,]
misc_otu_PA <- 1*((misc_otu>0)==1)
misc_otu_PA <- misc_otu_PA[rowSums(misc_otu)>0,]
Occ_misc <- rowSums(misc_otu_PA)/ncol(misc_otu_PA)
misc_otu_rel <- decostand(misc_otu, method="total", MARGIN=2)
com_abund_misc <- rowSums(misc_otu_rel)

color_misc_top <- com_abund_misc
color_misc_top[] <- 'black' 

#Mean Occ_Abund
#Switch 2016
swit_otu_PA <- swit_otu_PA[rowSums(swit_otu_PA)>0,]
Mean_swit_otu_PA <- apply(swit_otu_PA, 1, mean)
Mean_Occ_swit <- rowSums(swit_otu_PA)/ncol(swit_otu_PA)
swit_otu <- swit_otu[rowSums(swit_otu)>0,]
swit_otu_rel <- decostand(swit_otu, method="total", MARGIN=2)
Mean_abund_swit <- apply(swit_otu_rel, 1, mean)

#2017
swit17_otu <- swit_otu_17[rowSums(swit_otu_17)>0,]
swit17_otu_PA <- swit17_otu_PA[rowSums(swit17_otu_PA)>0,]
Mean_swit17_otu_PA <- apply(swit17_otu_PA, 1, mean)
Mean_Occ_swit17 <- rowSums(swit17_otu_PA)/ncol(swit17_otu_PA)
swit17_otu_rel <- decostand(swit17_otu, method="total", MARGIN=2)
Mean_abund_swit17 <- apply(swit17_otu_rel, 1, mean)

#Miscanthus Occ_Abund
misc_otu_PA <- 1*((misc_otu>0)==1)
misc_otu_PA <- misc_otu_PA[rowSums(misc_otu)>0,]
Occ_misc <- rowSums(misc_otu_PA)/ncol(misc_otu_PA)
misc_otu <- misc_otu[rowSums(misc_otu)>0,]
misc_otu_rel <- decostand(misc_otu, method="total", MARGIN=2)
Mean_abund_misc <- apply(misc_otu_rel, 1, mean)

#Creating df for plotting with ggplot and adding color code for the shared and unique OTUs
swit_df_occ <- data.frame(otu=names(Mean_Occ_swit), occ=Mean_Occ_swit) 
swit_df_abun <- data.frame(otu=names(Mean_abund_swit), abun=log10(Mean_abund_swit))
switc_col <- data.frame(otu=names(color_swit_top), col=color_swit_top)
swit_occ_abun <- left_join(swit_df_abun, swit_df_occ, by='otu')
swit_occ_abun <- left_join(swit_occ_abun, switc_col, by='otu')
swit_occ_abun$unique <- 'shared'

swit17_df_occ <- data.frame(otu=names(Mean_Occ_swit17), occ=Mean_Occ_swit17) 
swit17_df_abun <- data.frame(otu=names(Mean_abund_swit17), abun=log10(Mean_abund_swit17))
switc17_col <- data.frame(otu=names(color_swit17_top), col=color_swit17_top)
swit17_occ_abun <- left_join(swit17_df_abun, swit17_df_occ, by='otu')
swit17_occ_abun <- left_join(swit17_occ_abun, switc17_col, by='otu')
swit17_occ_abun$unique <- 'shared'

misc_df_occ <- data.frame(otu=names(Occ_misc), occ=Occ_misc) 
misc_df_abun <- data.frame(otu=names(Mean_abund_misc), abun=log10(Mean_abund_misc))
misc_col <- data.frame(otu=names(color_misc_top), col=color_misc_top)
misc_occ_abun <- left_join(misc_df_abun, misc_df_occ, by='otu')
misc_occ_abun <- left_join(misc_occ_abun, misc_col, by='otu')
misc_occ_abun$unique <- 'shared'

shared_v1 <- swit_occ_abun$otu[swit_occ_abun$otu %in% misc_occ_abun$otu]
shared <- shared_v1[shared_v1 %in% swit17_occ_abun$otu] 
shared_season_switch <- swit_occ_abun$otu[swit_occ_abun$otu %in% swit17_occ_abun$otu] #187 shared between the 2016-2017 (298 OTUs total in 2017)

shared_season_switch_only <- shared_season_switch[!(shared_season_switch %in% shared)]

#swit_occ_abun$unique[swit_occ_abun$otu %in% shared] <- 'shared (n=144)'
swit_occ_abun$unique[!(swit_occ_abun$otu %in% shared_season_switch)] <- 'Switchgrass 2016 (n=563)'
#swit_occ_abun$unique[swit_occ_abun$otu %in% shared_season_switch_only] <- '2016-2017 (n=43)'

misc_occ_abun$unique[!(misc_occ_abun$otu %in% shared)] <- 'Miscanthus 2016 (n=486)' 

#swit17_occ_abun$unique[swit17_occ_abun$otu %in% shared] <- 'shared (n=144)'
swit17_occ_abun$unique[!(swit17_occ_abun$otu %in% shared_season_switch)] <- 'Switchgrass 2017 (n=111)'
#swit17_occ_abun$unique[swit17_occ_abun$otu %in% shared_season_switch_only] <- '2016-2017 (n=43)'


#####################################
#Figure 3A - Occupancy abundance plot
#####################################

FigA<- ggplot(data=swit_occ_abun, aes(x=abun, y=occ, fill=unique)) +
  theme_bw()+
  geom_point(size=3, pch=21, alpha=.5) +
  scale_fill_manual(breaks=unique, values=c('white','darkolivegreen3')) +
  labs(x=paste('log(mean relative abundace per OTU)\n (n=',nrow(swit_df_abun),'OTUs)',sep=' '), y=paste('Mean occupancy (n=',ncol(swit_otu_PA),')', sep=' '), title='A', fill=NULL) +
  geom_text_repel(data=swit_occ_abun[swit_occ_abun$otu %in% shared_otus$otu,], aes(label=otu), box.padding = unit(0.45, "lines")) +
  geom_segment(aes(x = -2.5, y = .4, xend = -2.5, yend = Inf), linetype= 'dashed') +
  geom_segment(aes(x = -2.5, y = .4, yend = .4, xend = Inf), linetype= 'dashed') +
  theme(legend.position = 'right',
        legend.background = element_rect(fill=alpha(0.1))) +
  scale_y_continuous(breaks=seq(0,1,.2)) +
  xlim(-4,-0.5) +
  guides(fill = guide_legend(override.aes = list(alpha = 1)))
ggsave(filename = 'Figures/switchgrass_occ_abund_2016.pdf', device = 'pdf', width = 6.5, height = 4)

FigB <- ggplot(data=swit17_occ_abun, aes(x=abun, y=occ, fill=unique)) +
  theme_bw()+
  geom_point(size=3, pch=21, alpha=.5) +
  scale_fill_manual(breaks=unique, values=c('white', 'darkolivegreen1')) +
  labs(x=paste('log(mean relative abundace per OTU)\n (n=',nrow(swit17_occ_abun),' OTUs)',sep=''), y=paste('Mean occupancy (n=',ncol(swit17_otu_PA),')', sep=''), title='B', fill=NULL) +
  geom_text_repel(data=swit17_occ_abun[swit17_occ_abun$otu %in% shared_otus$otu,], aes(label=otu), box.padding = unit(0.45, "lines")) +
  geom_segment(aes(x = -2.5, y = .4, xend = -2.5, yend = Inf), linetype= 'dashed') +
  geom_segment(aes(x = -2.5, y = .4, yend = .4, xend = Inf), linetype= 'dashed') +
  theme(legend.position = 'right',
        legend.background = element_rect(fill=alpha(0.1))) +
  xlim(-4,-0.5) +
  scale_y_continuous(breaks=seq(0,1,.2)) +
  #scale_x_continuous(breaks=seq(-4,0,1)) +
  guides(fill = guide_legend(override.aes = list(alpha = 1)))
ggsave(filename = 'Figures/switchgrass_occ_abund_2017.pdf', device = 'pdf', width = 6.5, height = 4)

FigC <- ggplot(data=misc_occ_abun, aes(x=abun, y=occ, fill=unique)) +
  theme_bw()+
  geom_point(size=3, pch=21, alpha=.5) +
  scale_fill_manual(breaks=unique, values=c('darkgreen', 'white')) +
  labs(x=paste('log(mean relative abundace per OTU)\n (n=',nrow(misc_df_abun),' OTUs)',sep=''), y=paste('Mean occupancy (n=',ncol(misc_otu_PA),')', sep=''), title='C', fill=NULL) +
  geom_text_repel(data=misc_occ_abun[misc_occ_abun$otu %in% shared_otus$otu,], aes(label=otu), box.padding = unit(0.45, "lines")) +
  geom_segment(aes(x = -2.5, y = .4, xend = -2.5, yend = Inf), linetype= 'dashed') +
  geom_segment(aes(x = -2.5, y = .4, yend = .4, xend = Inf), linetype= 'dashed') +
  theme(legend.position = 'right',
        legend.background = element_rect(fill=alpha(0.1))) +
  scale_y_continuous(breaks=seq(0,1,.2)) +
  xlim(-4,-0.5) +
  guides(fill = guide_legend(override.aes = list(alpha = 1)))
ggsave(filename = 'Figures/miscanthus_occ_abund_2016.pdf', device = 'pdf', width = 6.5, height = 4)

###Venn diagram for comparing the datasets
install.packages('VennDiagram')
library(VennDiagram)
grid.newpage()
venn_plot <- draw.triple.venn(area1 = 750, area2 = 298,area3= 630, n12 = 187, n23 = 164, n13 = 301, n123 = 144, 
                              fill= c('darkolivegreen3', 'darkolivegreen1','darkgreen'),
                              category = c('A', 'B', 'C'))




#collecting all the presence/absence data for each rep plot
#Collecting the OTUs that fall in the marked quadrant in the Occ_Abund plots
occ_swit_df <- data.frame(otu=names(Mean_Occ_swit), occ_2016=Mean_Occ_swit) %>% filter(occ_2016 > .4)
abund_otu_switch <- data.frame(otu=names(Mean_abund_swit), abun_2016=log10(Mean_abund_swit)) %>% filter(abun_2016 > -2.5)
switch_occ_abund <- full_join(occ_swit_df, abund_otu_switch)
swit_uniq_otu <- swit_occ_abun[swit_occ_abun$otu %in% switch_occ_abund$otu,]
switch_occ_abund_v2 <- left_join(switch_occ_abund, swit_uniq_otu, by='otu')
switch_occ_abund_v2 <- switch_occ_abund_v2[,-c(4,5,6)]
names(switch_occ_abund_v2)[4] <- 'unique_swit_2016'

occ_swit17_df <- data.frame(otu=names(Mean_Occ_swit17), occ_2017=Mean_Occ_swit17) %>% filter(occ_2017 > .4)
abund_otu_switch17 <- data.frame(otu=names(Mean_abund_swit17), abun_2017=log10(Mean_abund_swit17)) %>% filter(abun_2017 > -2.5)
switch17_occ_abund <- full_join(occ_swit17_df, abund_otu_switch17)
swit17_uniq_otu <- swit17_occ_abun[swit17_occ_abun$otu %in% switch17_occ_abund$otu,]
switch17_occ_abund_v2 <- left_join(switch17_occ_abund, swit17_uniq_otu, by='otu')
switch17_occ_abund_v2 <- switch17_occ_abund_v2[,-c(4,5,6)]
names(switch17_occ_abund_v2)[4] <- 'unique_swit_2017'

occ_misc_df<- data.frame(otu=names(Occ_misc), occ_misc=Occ_misc) %>% filter(occ_misc > .4)
abund_otu_misc <- data.frame(otu=names(Mean_abund_misc), abund_misc=log10(Mean_abund_misc)) %>% filter(abund_misc > -2.5)
misc_occ_abund <- full_join(occ_misc_df, abund_otu_misc)
misc_uniq_otu <- misc_occ_abun[misc_occ_abun$otu %in% misc_occ_abund$otu,]
misc_occ_abund_v2 <- left_join(misc_occ_abund, misc_uniq_otu, by='otu')
misc_occ_abund_v2 <- misc_occ_abund_v2[,-c(4,5,6)]
names(misc_occ_abund_v2)[4] <- 'unique_misc_2016'

phyllo_occ_abund <- full_join(switch_occ_abund_v2, misc_occ_abund_v2)
phyllo_occ_abund <- full_join(phyllo_occ_abund,switch17_occ_abund_v2)

#OTUs shared between the plants:
phyllo_occ_abund[phyllo_occ_abund=='NA'] <- NA
shared_otus <- phyllo_occ_abund[complete.cases(phyllo_occ_abund),] #15 OTUs shared between the plants
dim(shared_otus)

core_switch16 <- swit_occ_abun[swit_occ_abun$abun>-2.5 & swit_occ_abun$occ>.4,]
core_switch17 <- swit17_occ_abun[swit17_occ_abun$abun>-2.5 & swit17_occ_abun$occ>.4,]
core_misc <- misc_occ_abun[misc_occ_abun$abun>-2.5 & misc_occ_abun$occ>.4,]

#write.table('~/Dropbox/GLBRC_16S/Supplemental/core_switch')
dim(core_misc)
length(core_misc$otu[!(core_misc$otu %in% shared_otus$otu)])

#Switchgrass 2016 all:
switch_selected <- swit_occ_abun[swit_occ_abun$abun>-2.5 & swit_occ_abun$occ>.4,] #total 22 OTUs in the quadrant
#Switchgrass specific:
switch_spec_selected <- core_switch16[!(core_switch16$otu %in% shared_otus$otu),] #7 OTUs unique to switch
tax_short[rownames(tax_short) %in% switch_selected$otu,]

#Switchgrass 2017 all:
switch17_selected <- swit17_occ_abun[swit17_occ_abun$abun>-2.5 & swit17_occ_abun$occ>.4,] #total 25 OTUs in the quadrant
#Switchgrass specific:
switch17_spec_selected <- core_switch17[!(core_switch17$otu %in% shared_otus$otu),] #10 OTUs unique to switch17
tax_short[rownames(tax_short) %in% switch17_selected$otu,]

#Miscanthus all:
misc_selected <- misc_occ_abun[misc_occ_abun$abun>-2.5 & misc_occ_abun$occ>.4,] # 17 OTUs
#Miscanthus specific:
misc_spec_selected <- core_misc[!(core_misc$otu %in% shared_otus$otu),] #6 OTUs
tax_short[rownames(tax_short) %in% misc_selected$otu,]

combined_selected_otu <- full_join(switch_selected,switch17_selected)
combined_selected_otu <- full_join(combined_selected_otu,misc_selected)

dim(combined_selected_otu)

rel_otu_rare <- decostand(otu_rare, method="total", MARGIN=2)
abund_selected_all <- rel_otu_rare[rownames(rel_otu_rare) %in% combined_selected_otu$otu,]

abund_selected_all[order(rowSums(abund_selected_all))]

tax_short$otu <- rownames(tax_short)

selected_otus <- data.frame(otu = as.factor(row.names(abund_selected_all)), abund_selected_all) %>% gather(sequence_name, abun, -otu) %>% 
  left_join(map_16S[, c('sequence_name','rep','time_numeric', 'treatment' ,'source', 'plant', 'month', 'sampling_date', 
                        'carbon_per_nitrogen', 'nitrogen_percent', 'carbon_percent', 'year')], by = 'sequence_name') %>%
  left_join(tax_short, by='otu')

lastValue <- function(x) tail(x[!is.na(x)], 1)
last_taxons<- apply(selected_otus, 1, lastValue)
selected_otus$last_taxon <- last_taxons
head(selected_otus)
selected_otus$final_names <- paste(selected_otus$last_taxon, selected_otus$otu, sep=' - ')
unique(selected_otus$final_names[selected_otus$plant=='switchgrass'])

selected_otus_switch <- selected_otus[selected_otus$otu %in% factor(switch_selected$otu),]
selected_otus_misc <- selected_otus[selected_otus$otu %in% factor(misc_selected$otu),]
selected_otus_switch17 <- selected_otus[selected_otus$otu %in% factor(switch17_selected$otu),]

length(unique(selected_otus_switch$otu))
length(unique(selected_otus_switch17$otu))
length(unique(selected_otus_misc$otu))


##############################################################################
###Figure 3B - PA data representation for selected OTUs per replicate and time
##############################################################################
#Switchgrass
#df with stats for the time points with all replicates and per OTU
selected_otus_switch %>%
  filter(source == 'phyllosphere' & plant == 'switchgrass' & year == 2016) %>%
  group_by(sampling_date, final_names, Class, Order, Family, Genus) %>%
  dplyr::summarise(n=sum(abun>0)/length(abun),
                   all=length(abun),
                   rep_ab=mean(abun),
                   sd_rep=sd(abun)
  ) %>%
  filter(n>0) -> temp

temp %>%
  group_by(sampling_date) %>%
  summarise(all_n=unique(all))

#df with stats for the whole dataset per OTU
selected_otus_switch %>%
  filter(source == 'phyllosphere' & plant == 'switchgrass' & year == 2016) %>%
  group_by(final_names, otu) %>%
  dplyr::summarise(
    all_ab=mean(abun),
    all_sd=sd(abun)
  ) -> temp2

swit_selected_otu_list <- unique(selected_otus_switch$otu)

#combining df and calculating the z-score
z_df <- left_join(temp, temp2)
z_df %>% arrange(otu) %>%
  mutate(
    z_score=(rep_ab-all_ab)/all_sd
  ) %>%
  arrange(Class, Family) ->tmp3

#Miscanthus
#df with stats for the time points with all replicates and per OTU
selected_otus_misc %>%
  filter(source == 'phyllosphere' & plant == 'miscanthus') %>%
  group_by(sampling_date, final_names, Class, Order, Family, Genus) %>%
  dplyr::summarise(n=sum(abun>0)/length(abun),
                   all=length(abun),
                   rep_ab=mean(abun),
                   sd_rep=sd(abun)
  ) %>%
  filter(n>0) -> temp4

temp4 %>%
  group_by(sampling_date) %>%
  summarise(all_n=unique(all))

#df with stats for the whole dataset per OTU
selected_otus_misc %>%
  filter(source == 'phyllosphere' & plant == 'miscanthus') %>%
  group_by(final_names, otu) %>%
  dplyr::summarise(
    all_ab=mean(abun),
    all_sd=sd(abun)
  ) -> temp5

misc_selected_otu_list <- unique(selected_otus_misc$otu)

#combining df and calculating the z-score
z_df_misc <- left_join(temp4, temp5)
z_df_misc %>% arrange(otu) %>%
  mutate(
    z_score=(rep_ab-all_ab)/all_sd
  ) %>%
  arrange(Class, Family) -> temp6

#Switchgrass 2017
#OTU mean abundance per sampling time 
selected_otus_switch17 %>%
  filter(source == 'phyllosphere' & plant == 'switchgrass' & year == 2017) %>%
  group_by(sampling_date, final_names, Class, Order, Family, Genus) %>%
  dplyr::summarise(n=sum(abun>0)/length(abun),
                   all=length(abun),
                   rep_ab=mean(abun),
                   sd_rep=sd(abun)
  ) %>%
  filter(n>0) -> temp7

#df with stats for the whole dataset per OTU
selected_otus_switch17 %>%
  filter(source == 'phyllosphere' & plant == 'switchgrass' & year == 2017) %>%
  group_by(final_names) %>%
  dplyr::summarise(
    all_ab=mean(abun),
    all_sd=sd(abun)
  ) -> temp8

swit17_selected_otu_list <- unique(selected_otus_switch17$otu)

#combining df and calculating the z-score
z_df <- left_join(temp7, temp8)
z_df %>% arrange(final_names) %>%
  mutate(
    z_score=(rep_ab-all_ab)/all_sd
  ) %>%
  arrange(Class, Family) -> tmp9

same_across_data <- tmp3[tmp3$final_names %in% temp6$final_names,]
same_across_data <- same_across_data[same_across_data$final_names %in% tmp9$final_names,]
unique(same_across_data$otu)

combined_otus_misc <- temp6[temp6$final_names %in% same_across_data$final_names,]
uniq_misc <- temp6[!(temp6$final_names %in% same_across_data$final_names),]
length(unique(uniq_misc$final_names))

combined_otus_switch <- tmp3[tmp3$final_names %in% same_across_data$final_names,] #for switchgrass16 dataset
uniq_switch <- tmp3[!(tmp3$final_names %in% same_across_data$final_names),]

combined_otus_switch17 <- tmp9[tmp9$final_names %in% same_across_data$final_names,] #for switchgrass17 dataset
uniq_switch17 <- tmp9[!(tmp9$final_names %in% same_across_data$final_names),]


plot3 <- ggplot(combined_otus_misc,aes(x=factor(sampling_date), y=factor(combined_otus_misc$final_names, levels=rev(unique(combined_otus_misc$final_names))), 
                                       size=n)) +
  theme_bw()+
  geom_point(pch=21, colour='black', aes(fill=z_score)) +
  scale_fill_gradient2(low = "#b2182b", high ="#2166ac") +
  labs(x=NULL, y= NULL, title='Miscanthus (2016)') +
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        axis.text.y = element_text(size=10),
        axis.text.x=element_blank())

plot1 <- ggplot(combined_otus_switch,aes(x=factor(sampling_date), y=factor(combined_otus_switch$final_names, levels=rev(unique(combined_otus_switch$final_names))), 
                                         size=n)) +
  theme_bw()+
  geom_point(pch=21, colour='black', aes(fill=z_score)) +
  scale_fill_gradient2(low = "#b2182b", high ="#2166ac") +
  labs(x=NULL, y= NULL, title='Switchgrass (2016)', fill='z-score', size='Occupancy') +
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        axis.text.y = element_text(size=10),
        axis.text.x=element_blank(),
        legend.position = 'bottom')

plot5 <- ggplot(combined_otus_switch17,aes(x=factor(sampling_date), y=factor(combined_otus_switch17$final_names, 
                                                                             levels=rev(unique(combined_otus_switch17$final_names))), 
                                           size=n)) +
  theme_bw()+
  geom_point(pch=21, colour='black', aes(fill=z_score)) +
  scale_fill_gradient2(low = "#b2182b", high ="#2166ac") +
  labs(x=NULL, y= NULL, title='Switchgrass (2017)', fill='z-score', size='Occupancy') +
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        axis.text.y = element_text(size=10),
        axis.text.x=element_blank())

plot4 <- ggplot(uniq_misc,aes(x=factor(sampling_date), y=factor(uniq_misc$final_names, levels=rev(unique(uniq_misc$final_names))), 
                              size=n)) +
  theme_bw()+
  geom_point(pch=21, colour='black', aes(fill=z_score)) +
  scale_fill_gradient2(low = "#b2182b", high ="#2166ac") +
  labs(x="Sampling times", y= NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust=1, size=12),
        axis.text.y = element_text(size=10))

plot2 <- ggplot(uniq_switch,aes(x=factor(sampling_date), y=factor(uniq_switch$final_names, levels=rev(unique(uniq_switch$final_names))), 
                                size=n)) +
  theme_bw()+
  geom_point(pch=21, colour='black', aes(fill=z_score)) +
  scale_fill_gradient2(low = "#b2182b", high ="#2166ac") +
  labs(x="Sampling times", y= NULL, fill='z-score', size='Occupancy') +
  theme(axis.text.y = element_text(size=10),
        axis.text.x = element_text(angle = 45, hjust=1, size=12))

plot6 <- ggplot(uniq_switch17,aes(x=factor(sampling_date), y=factor(uniq_switch17$final_names, levels=rev(unique(uniq_switch17$final_names))), 
                                  size=n)) +
  theme_bw()+
  geom_point(pch=21, colour='black', aes(fill=z_score)) +
  scale_fill_gradient2(low = "#b2182b", high ="#2166ac") +
  labs(x="Sampling times", y= NULL, fill='z-score', size='Occupancy') +
  theme(axis.text.y = element_text(size=10),
        axis.text.x = element_text(angle = 45, hjust=1, size=12))

#function for extracting the legend in the ggplot figure now in package 'lemon'
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend <- g_legend(plot1)

plot1 <- plot1 + theme(legend.position="none")
plot3 <- plot3 + theme(legend.position="none")
plot2 <- plot2 + theme(legend.position="none")
plot4 <- plot4 + theme(legend.position="none")
plot5 <- plot5 + theme(legend.position="none")
plot6 <- plot6 + theme(legend.position="none")

library(egg)
setEPS()
postscript('Figures/Figure3B-2_v2.eps', width = 16, height = 9)
grid.draw(ggarrange(plot3,
                    plot1,
                    plot5,
                    plot4,
                    plot2,
                    plot6, heights = c(1.5,.8)))
dev.off()
####Replotting the Fig 3
combined_otus_switch17$tag <- 'switch17'
combined_otus_switch$tag <- 'switch16'
combined_otus_misc$tag <- 'misc16'
uniq_misc$tag <- 'misc16'
uniq_switch$tag <- 'switch16'
uniq_switch17$tag <- 'switch17'

combined_all <- rbind(combined_otus_switch17,combined_otus_switch)
combined_all <- rbind(combined_all,combined_otus_misc)
combined_all <- rbind(combined_all,uniq_misc)
combined_all <- rbind(combined_all,uniq_switch)
combined_all <- rbind(combined_all,uniq_switch17)

combined_all$sampling_week <- 0
combined_all$sampling_week[combined_all$sampling_date == '2017-05-15 EDT'] <- 1 
combined_all$sampling_week[combined_all$sampling_date == '2016-05-09 EDT'] <- 1

combined_all$sampling_week[combined_all$sampling_date == '2016-05-31 EDT'] <- 2
combined_all$sampling_week[combined_all$sampling_date == '2017-06-05 EDT'] <- 2 

combined_all$sampling_week[combined_all$sampling_date == '2016-06-20 EDT'] <- 3
combined_all$sampling_week[combined_all$sampling_date == '2017-06-26 EDT'] <- 3

combined_all$sampling_week[combined_all$sampling_date == '2016-07-12 EDT'] <- 4 
combined_all$sampling_week[combined_all$sampling_date == '2017-07-17 EDT'] <- 4

combined_all$sampling_week[combined_all$sampling_date == '2016-08-01 EDT'] <- 5 
combined_all$sampling_week[combined_all$sampling_date == '2017-08-07 EDT'] <- 5 

combined_all$sampling_week[combined_all$sampling_date == '2016-08-22 EDT'] <- 6 
combined_all$sampling_week[combined_all$sampling_date == '2017-08-28 EDT'] <- 6 

combined_all$sampling_week[combined_all$sampling_date == '2016-09-12 EDT'] <- 7 
combined_all$sampling_week[combined_all$sampling_date == '2017-09-18 EDT'] <- 7 

combined_all$sampling_week[combined_all$sampling_date == '2016-10-03 EDT'] <- 8

combined_all$sampling_week[combined_all$sampling_date == '2016-11-07 EST'] <- 9

class(combined_all$tag) 

#Figure 3
ggplot(combined_all,aes(x=factor(sampling_week), y=as.factor(tag), 
                        size=n)) +
  theme_bw()+
  geom_point(pch=21, colour='black', aes(fill=z_score)) +
  scale_fill_gradient2(low = "#b2182b", high ="#2166ac") +
  labs(x='Sampling week', y= NULL, title='Core') +
  theme(
        axis.text.y = element_text(size=10),
        strip.text = element_text(size=7)) + 
  facet_wrap(~final_names, labeller = label_wrap_gen(), ncol=4) +
  scale_y_discrete(limits=rev(levels(factor(combined_all$tag))))

###############################################################
#Suppl Fig - Dynamics of the OTUs classified as Alpha and Gamma proteobacteria
###############################################################
rel_otu_rare <- decostand(otu_rare, method="total", MARGIN=2)

proteo_plot <- data.frame(otu = as.factor(row.names(otu_rare)), otu_rare) %>% gather(sequence_name, abun, -otu) %>% 
  left_join(map_16S[, c('sequence_name','rep','treatment' ,'source', 'plant', 'sampling_date', 'year'
  )], by = 'sequence_name') %>%
  left_join(tax_short, by='otu') %>% 
  filter(grepl('Alphaproteobacteria|Betaproteobacteria|Deltaproteobacteria|Gammaproteobacteria', Class)) %>%
  filter(source == 'phyllosphere') %>%
  mutate(plant = factor(plant, levels = c('switchgrass', 'miscanthus'))) %>%
  group_by(plant, Class, sampling_date) %>%
  summarise(n=sum(abun),
            n_reps=length(unique(sequence_name))) %>%
  group_by(plant, sampling_date) %>%
  mutate(total_reads=sum(n),
         rel_abun=n/total_reads) %>%
  ggplot(aes(x=as.factor(sampling_date), y=rel_abun, fill=Class)) +
  geom_bar(color='black', stat = 'identity') +
  theme_classic() +
  facet_grid(~plant, scales='free_x') +
  labs(x='Sampling time', y='Normalized relative abundance')+
  theme(axis.text.x = element_text(angle = 45, hjust=1, size=8),
        legend.position = 'bottom') +
  guides(fill = guide_legend(ncol = 3, title=NULL))

#data organisation for plotting the Pseudomonas, Sphingomonas and Methylobacteria dynamics
selected_otus
selected_otus$sampling_week <- 0
selected_otus$sampling_week[map_16S$sampling_date == '2017-05-15 EDT'] <- 1 
selected_otus$sampling_week[map_16S$sampling_date == '2016-05-09 EDT'] <- 1
map_16S$sampling_week[map_16S$sampling_date == '2016-05-31 EDT'] <- 2
map_16S$sampling_week[map_16S$sampling_date == '2017-06-05 EDT'] <- 2 
map_16S$sampling_week[map_16S$sampling_date == '2016-06-20 EDT'] <- 3
map_16S$sampling_week[map_16S$sampling_date == '2017-06-26 EDT'] <- 3
map_16S$sampling_week[map_16S$sampling_date == '2016-07-12 EDT'] <- 4 
map_16S$sampling_week[map_16S$sampling_date == '2017-07-17 EDT'] <- 4
map_16S$sampling_week[map_16S$sampling_date == '2016-08-01 EDT'] <- 5 
map_16S$sampling_week[map_16S$sampling_date == '2017-08-07 EDT'] <- 5 
map_16S$sampling_week[map_16S$sampling_date == '2016-08-22 EDT'] <- 6 
map_16S$sampling_week[map_16S$sampling_date == '2017-08-28 EDT'] <- 6 
map_16S$sampling_week[map_16S$sampling_date == '2016-09-12 EDT'] <- 7 
map_16S$sampling_week[map_16S$sampling_date == '2017-09-18 EDT'] <- 7 
map_16S$sampling_week[map_16S$sampling_date == '2016-10-03 EDT'] <- 8
map_16S$sampling_week[map_16S$sampling_date == '2016-11-07 EST'] <- 9

methylobac<- combined_all[combined_all$Genus =='g:Methylobacterium',]
pseudomonas<- combined_all[combined_all$Genus =='g:Pseudomonas',]
sphingomonas <- combined_all[combined_all$Genus =='g:Sphingomonas',]

filtered_methylo <- sphingomonas %>% 
  #mutate(plant = factor(plant, levels = c('switchgrass', 'miscanthus')))  %>%
  filter(!is.na(c(sampling_week))) 

# sphing_aov <- aov(abun ~ factor(sampling_date)+plant,  data = filtered_methylo[filtered_methylo$otu=='2',])
# summary(sphing_aov)
# TukeyHSD(sphing_aov)
library(broom)

sphingo_result = filtered_methylo %>%
  filter(!(sampling_week %in% c(8,9))) %>%
  group_by(sampling_week, final_names) %>%
  do(tidy(t.test(.$abun~.$tag))) %>%
  mutate(FDR=p.adjust(p.value, method='fdr')) %>%
  filter(p.value < .05)

methylo_result$sampling_date <- as.factor(methylo_result$sampling_date)
methylo_result
sphingo_result
pseudo_result

sphingo <- ggplot(filtered_methylo, aes(x = as.factor(sampling_date), y = abun, fill=plant)) + 
  geom_boxplot() +
  theme_bw()+
  #scale_fill_manual(values=c('darkgreen','darkolivegreen3')) +
  labs(x="Sampling times", y= "Relative abundance") +
  facet_grid(plant ~ final_names, scale = 'free_y') +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        plot.title = element_text(hjust = 0.5),
        legend.position='bottom', 
        legend.background = element_rect(fill="transparent"),
        strip.text.x = element_blank()) +
  guides(fill = guide_legend(ncol=3, title=NULL))

install.packages("ggpubr")
library(ggpubr)

setEPS()
postscript('Figures/sphingo_dynamics.eps', width = 6, height = 5)
sphingo
dev.off()

big.phyllo.taxa <- data.frame(otu = as.factor(row.names(otu_rare)), otu_rare) %>% gather(sequence_name, abun, -otu) %>% 
  left_join(map_16S[, c('sequence_name','treatment' ,'source', 'plant', 'sampling_date', 'month'
  )], by = 'sequence_name') %>%
  left_join(tax_short, by='otu') %>% 
  filter(source=='phyllosphere') %>% 
  group_by(plant,Kingdom,Phylum,Class,Family,Genus,otu, month) %>%
  summarise(n=sum(abun)) %>%
  arrange(desc(n)) 
  
plant_otu <- otu_rare[,map_16S$source=="phyllosphere"]
plant_otu <- plant_otu[rowSums(plant_otu)>0,]
soil_otu <- otu_rare[,map_16S$source=="soil"] 
soil_otu <- soil_otu[rowSums(soil_otu)>0,]

total_sum_soil_otu <- as.data.frame(rowSums(soil_otu))
sum_soil_otu <- sum(total_sum_soil_otu$`rowSums(soil_otu)`)
total_sum_soil_otu$rel_abun <- total_sum_soil_otu$`rowSums(soil_otu)`/sum_soil_otu
#total_sum_soil_otu <- arrange(total_sum_soil_otu,desc(rel_abun))
sum(total_sum_soil_otu$rel_abun)
phylloINsoil=total_sum_soil_otu[rownames(total_sum_soil_otu) %in% rownames(plant_otu),]
sum(phylloINsoil$rel_abun) #OTUs find in phyllospere make up to 55.6% of soil community size 

#Dynamics of the unique taxa
####################################
switch_unique_otus <- rel_otu_rare[rownames(rel_otu_rare) %in% swit_occ_abun$otu[swit_occ_abun$unique=='unique'],]
(switch_unique <- data.frame(otu = as.factor(row.names(switch_unique_otus)), switch_unique_otus) %>% gather(sequence_name, abun, -otu) %>% 
    left_join(map_16S[, c('sequence_name','rep','treatment' ,'source', 'plant', 'sampling_date', 'month'
    )], by = 'sequence_name') %>%
    left_join(tax_short, by='otu') %>% 
    filter(source == 'phyllosphere' & plant == 'switchgrass' & month >= 8) %>%
    group_by(Genus) %>%
    filter(abun>0) %>%
    ggplot(aes(x=as.factor(sampling_date), y=abun)) +
    geom_boxplot() +
    theme_classic() +
    #facet_wrap( ~ as.factor(Class), scale = 'free_y') +
    labs(x='Sampling time', y='Mean abundance')+
    theme(axis.text.x = element_text(angle = 45, hjust=1, size=6),
          legend.position = 'bottom'))

setEPS()
postscript('Figures/switchgrass_unique_dynamics.eps', height = 9, width = 12)
switch_unique
dev.off()

misc_unique_otus <- rel_otu_rare[rownames(rel_otu_rare) %in% misc_occ_abun$otu[misc_occ_abun$unique=='unique'],]
(misc_unique <- data.frame(otu = as.factor(row.names(misc_unique_otus)), misc_unique_otus) %>% gather(sequence_name, abun, -otu) %>% 
    left_join(map_16S[, c('sequence_name','rep','treatment' ,'source', 'plant', 'sampling_date'
    )], by = 'sequence_name') %>%
    left_join(tax_short, by='otu') %>% 
    filter(source == 'phyllosphere' & plant == 'miscanthus') %>%
    group_by(sampling_date) %>%
    mutate(n=mean(abun)) %>%
    ggplot(aes(x=as.factor(sampling_date), y=abun)) +
    geom_boxplot() +
    theme_classic() +
    facet_wrap( ~ as.factor(Class), scale = 'free_y') +
    labs(x='Sampling time', y='Mean abundance')+
    theme(axis.text.x = element_text(angle = 45, hjust=1, size=6),
          legend.position = 'bottom'))

setEPS()
postscript('~/Dropbox/GLBRC_16S/Supplemental/miscanthus_unique_dynamics.eps', height = 9, width = 12)
misc_unique
dev.off()

####################################
#eLSA
###############################################################
#Data prep for the eLSA analysis
length(unique(selected_otus$otu))

swit_selected_otu_list <- unique(selected_otus_switch$otu)
misc_selected_otu_list <- unique(selected_otus_misc$otu)
swit17_selected_otu_list <- unique(selected_otus_switch17$otu)

slected_switch_otu_abundance<- otu_rare[rownames(otu_rare) %in% swit_selected_otu_list,]
slected_misc_otu_abundance<- otu_rare[rownames(otu_rare) %in% misc_selected_otu_list,]
slected_switch17_otu_abundance<- otu_rare[rownames(otu_rare) %in% swit17_selected_otu_list,]

switch_selected_otus <- selected_otus_switch%>%
  filter(source == 'phyllosphere' & plant == 'switchgrass' & year == 2016)
misc_selected_otus <- selected_otus_misc%>%
  filter(source == 'phyllosphere' & plant == 'miscanthus' & year == 2016)
switch17_selected_otus <- selected_otus_switch17 %>%
  filter(source == 'phyllosphere' & plant == 'switchgrass' & year == 2017)

#misc_selected_otus %>% 
#  group_by(sampling_date) %>%
#  summarise(n=length(unique(sequence_name)))
#length(unique(misc_selected_otus$sampling_date))

switch_timed_samples <- switch_selected_otus[order(switch_selected_otus$sampling_date),]
switch_timed_samples <- unique(switch_timed_samples$sequence_name)

misc_timed_samples <- misc_selected_otus[order(misc_selected_otus$sampling_date),]
misc_timed_samples <- unique(misc_timed_samples$sequence_name)

switch17_timed_samples <- switch17_selected_otus[order(switch17_selected_otus$sampling_date),]
switch17_timed_samples <- unique(switch17_timed_samples$sequence_name)

#organise abundance df so that the samples are time organised
LSA_switch_data <- slected_switch_otu_abundance[,colnames(slected_switch_otu_abundance) %in% switch_timed_samples]
LSA_switch_data <- LSA_switch_data[,switch_timed_samples]

LSA_misc_data <- slected_misc_otu_abundance[,colnames(slected_misc_otu_abundance) %in% misc_timed_samples]
LSA_misc_data <- LSA_misc_data[ ,misc_timed_samples]

LSA_switch17_data <- slected_switch17_otu_abundance[,colnames(slected_switch17_otu_abundance) %in% switch17_timed_samples]
LSA_switch17_data <- LSA_switch17_data[,switch17_timed_samples]
# nrow(LSA_misc_data)
# nrow(LSA_switch_data)
# nrow(LSA_switch17_data)

#write tables for the eLSA analysis
write.table(LSA_switch_data, 'InputFiles/LSAswitch16_input.txt', sep='\t')
write.table(LSA_misc_data, 'InputFiles/LSAmisc16_input.txt', sep='\t')
write.table(LSA_switch17_data, 'InputFiles/LSAswitch17_input.txt', sep='\t')
#######
#Import eLSA output tables for plotting 
s16.ELSA <- read.table('InputFiles/LSAswitch16_output.txt', header=T, sep='\t')
m.ELSA <- read.table('InputFiles/LSAmisc16_output.txt', header=T, sep='\t')
s17.ELSA <- read.table('InputFiles/LSAswitch17_output.txt', header=T, sep='\t')
# 
# swit16_plot_connection <- s16.ELSA %>% 
#   filter(P<0.05 & !(LS>-0.35 & LS<.35)) %>%
#   group_by(X) %>%
#   summarise(edgesNo=length(Y))%>%
#   arrange(desc(edgesNo)) %>%
#   ggplot(aes(x=(reorder(X, -edgesNo)), y=edgesNo)) +
#   geom_point() +
#   theme_bw()+
#   labs(x=NULL, y= 'Number of connections', title='Switchgrass 2016') + 
#   theme(plot.title = element_text(hjust = 0.5, size = 12),
#         axis.text.x = element_text(angle = 45, hjust=1, size=8.5))
# 
# swit17_plot_connection <- s17.ELSA %>% 
#   filter(P<0.05 & !(LS>-0.35 & LS<.35)) %>%
#   group_by(X) %>%
#   summarise(edgesNo=length(Y))%>%
#   arrange(desc(edgesNo)) %>%
#   ggplot(aes(x=(reorder(X, -edgesNo)), y=edgesNo)) +
#   geom_point() +
#   theme_bw()+
#   labs(x=NULL, y= 'Number of connections', title='Switchgrass 2017') + 
#   theme(plot.title = element_text(hjust = 0.5, size = 12),
#         axis.text.x = element_text(angle = 45, hjust=1, size=8.5))
# 
# misc_plot_connection <- m.ELSA %>% 
#   filter(P<0.05 & !(LS>-0.35 & LS<.35)) %>%
#   group_by(X) %>%
#   summarise(edgesNo=length(Y))%>%
#   arrange(desc(edgesNo)) %>%
#   ggplot(aes(x=(reorder(X, -edgesNo)), y=edgesNo)) +
#   geom_point() +
#   theme_bw()+
#   labs(x=NULL, y= 'Number of connections', title='Miscanthus 2016') + 
#   theme(plot.title = element_text(hjust = 0.5, size = 12),
#         axis.text.x = element_text(angle = 45, hjust=1, size=8.5))

s16.ELSA$plant <- 'Switchgrass 2016'
s17.ELSA$plant <- 'Switchgrass 2017'
m.ELSA$plant <- 'Miscanthus 2016'
eLSA_combined <- rbind(s16.ELSA, s17.ELSA, m.ELSA)

all_eLSA_plot_X_Y <- eLSA_combined %>% 
  filter(P<0.05 & !(LS>-0.35 & LS<.35)) %>%
  group_by(plant, X) %>%
  summarise(edgesNo=length(Y))%>%
  arrange(desc(edgesNo)) %>%
  ggplot(aes(x=(reorder(X, -edgesNo)), y=edgesNo)) +
  geom_point() +
  theme_bw()+
  facet_wrap(~plant, scales = "free_x") +
  labs(x=NULL, y= '# of connections building') + 
  theme(axis.text.x = element_text(angle = 90, hjust=1, size=8.5))

all_eLSA_plot_Y_X <- eLSA_combined %>% 
  filter(P<0.05 & !(LS>-0.35 & LS<.35)) %>%
  group_by(plant, Y) %>%
  summarise(edgesNo=length(X))%>%
  arrange(desc(edgesNo)) %>%
  ggplot(aes(x=(reorder(Y, -edgesNo)), y=edgesNo)) +
  geom_point() +
  theme_bw()+
  facet_wrap(~plant, scales = "free_x") +
  labs(x=NULL, y= '# of connections receiving') + 
  theme(axis.text.x = element_text(angle = 90, hjust=1, size=8.5))

names(all_eLSA_plot_X_Y)[2] <- 'otu'
names(all_eLSA_plot_Y_X)[2] <- 'otu'

connection_df <- rbind(all_eLSA_plot_X_Y, all_eLSA_plot_Y_X)
tax_otus_short <- unique(tax_short[,c('otu','Class')])
connection_df <- left_join(connection_df,tax_otus_short)

conn_plot <- connection_df %>%
  group_by(plant, otu, Class) %>%
  summarise(n=sum(edgesNo)) %>%
  arrange(desc(n)) %>%
  ggplot(aes(x=(reorder(otu, -n)), y=n, fill=Class)) +
  geom_point(pch=21, size=3) +
  theme_bw()+
  facet_wrap(~plant, scales = "free_x") +
  labs(x=NULL, y= '# of connections') +
  theme(axis.text.x = element_text(angle = 90, hjust=1, size=12, vjust=.5),
        legend.position = 'bottom')

switch16_core <- tmp3[,c(1,11,14)]
switch16_core <- as.data.frame(switch16_core)
#switch16_core_wide <- reshape(switch16_core, idvar='otu', timevar='sampling_date', direction='wide')
switch16_core_wide_df <- spread(switch16_core, key='sampling_date', value='z_score')
switch16_core_wide[is.na(switch16_core_wide)] <- 0
rownames(switch16_core_wide) <- switch16_core_wide$otu
switch16_core_wide$otu <-  NULL

set.seed(20)
clusters_switch16 <- hclust(dist(switch16_core_wide),'complete')
memb <- cutree(clusters_switch16, k=3)
plot(clusters_switch16)

switch17_core <- tmp9[,c(1,2,13)]
#switch17_core <- as.data.frame(switch17_core)
#switch16_core_wide <- reshape(switch16_core, idvar='otu', timevar='sampling_date', direction='wide')
switch17_core_wide <- spread(switch17_core, key='sampling_date', value='z_score')
switch17_core_wide[is.na(switch17_core_wide)] <- 0
rownames(switch17_core_wide) <- switch17_core_wide$final_names
switch17_core_wide$final_names <-  NULL
set.seed(21)
clusters_switch17 <- hclust(dist(switch17_core_wide),'average')
memb <- cutree(clusters_switch17, k=3)
plot(clusters_switch17)

misc16_core <- temp6[,c(1,2,14)]
#switch17_core <- as.data.frame(switch17_core)
#switch16_core_wide <- reshape(switch16_core, idvar='otu', timevar='sampling_date', direction='wide')
misc16_core_wide <- spread(misc16_core, key='sampling_date', value='z_score')
misc16_core_wide[is.na(misc16_core_wide)] <- 0
rownames(misc16_core_wide) <- misc16_core_wide$final_names
misc16_core_wide$final_names <-  NULL
set.seed(21)
clusters_misc16 <- hclust(dist(misc16_core_wide),'complete')
memb <- cutree(clusters_misc16, k=3)
plot(clusters_misc16)

#Plotting the abundance dynamics of the selected OTUs for Switchgrass and Miscanthus together
ggplot(selected_otus[selected_otus$source=='phyllosphere',], aes(x = as.factor(sampling_date), y = abun, fill = plant)) + 
  geom_boxplot() +
  scale_fill_manual(values=c('darkgreen','darkolivegreen3')) +
  labs(x="Sampling times", y= "Relative abundance", title='Relative abundance dynamics of the selected phyllosphere OTUs (n=35)') +
  facet_wrap( ~ as.factor(final_names), scale = 'free_y') +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        plot.title = element_text(hjust = 0.5))

swit_only<- selected_otus[selected_otus$source=='phyllosphere' & selected_otus$plant =='switchgrass',]

swit_only %>% 
  select(final_names, otu, Class) %>%
  filter(otu==) %>%
  arrange(Class)

group2_LSA_swit <- c('357','101', '30')
group1_LSA_swit <- c('10','21', '641','920', '1751')
group3_LSA_swit <- c('8','20', '9','53')
group4_LSA_swit <- c('12','6', '4','9', '19','108', '25', '2', '3')
antigroup1_LSA_swit <- c('12', '1751', '4', '641')
antigroup2_LSA_swit <- c('80', '63', '117', '37')

unique(swit_only$Class[swit_only$otu==c('10','2','61','40','44','157')])
swit_only_subset <- swit_only[(swit_only$otu %in% group4_LSA_swit),]

swit_LSA_4 <- ggplot(swit_only_subset, 
                     aes(x = as.factor(sampling_date), y = abun, fill=final_names)) + 
  geom_boxplot() +
  #scale_fill_manual(values=c('darkgreen','darkolivegreen3')) +
  labs(x="Sampling times", y= "Relative abundance") +
  #facet_wrap( ~ as.factor(final_names), scale = 'free_y') +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        plot.title = element_text(hjust = 0.5),
        legend.position='bottom', 
        legend.background = element_rect(fill="transparent")) +
  guides(fill = guide_legend(ncol = 2, title=NULL))

antigroup2 <- ggplot(swit_only_subset, 
                     aes(x = as.factor(sampling_date), y = abun, fill=final_names)) + 
  geom_boxplot() +
  #scale_fill_manual(values=c('darkgreen','darkolivegreen3')) +
  labs(x="Sampling times", y= "Relative abundance") +
  #facet_wrap( ~ as.factor(final_names), scale = 'free_y') +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        plot.title = element_text(hjust = 0.5),
        legend.position='bottom', 
        legend.background = element_rect(fill="transparent")) +
  guides(fill = guide_legend(ncol = 2, title=NULL))

setEPS()
postscript('Figures/eLSA_OTU_dynamics_group4.eps', width = 8, height = 6)
swit_LSA_4
dev.off()

misc_only<- selected_otus[selected_otus$source=='phyllosphere' & selected_otus$plant =='miscanthus',]
group1_LSA_misc <- c('10','30', '63','101', '641', '920')
group2_LSA_misc <- c('8','27', '16','32', '149', '83', '53', '40', '22','4')
group3_LSA_misc <- c('8','27', '16','32')
antigroup1_misc <- c('8','27', '16','32', '30', '10')
antigroup2_misc <- c('2','3', '10')
alphaVSgamma <- c('8','27', '16','32','22', '4', '9', '10', '641', '920', '101', '63')

misc_only_subset <- misc_only[(misc_only$otu %in% antigroup1_misc),]

antigroup <- ggplot(misc_only_subset, 
                    aes(x = as.factor(sampling_date), y = abun, fill=final_names)) + 
  geom_boxplot() +
  labs(x="Sampling times", y= "Relative abundance") +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        plot.title = element_text(hjust = 0.5),
        legend.position='bottom', 
        legend.background = element_rect(fill="transparent")) +
  guides(fill = guide_legend(ncol = 2, title=NULL))

setEPS()
postscript('Figures/eLSA_OTU_dynamics_antigroup1_misc.eps', width = 5, height = 5)
antigroup_1_misc
dev.off()

#Plotting the abundance dynamics of the selected OTUs for Switchgrass and Miscanthus together
ggplot(selected_otus[selected_otus$source=='phyllosphere',], aes(x = as.factor(sampling_date), y = abun, fill = plant)) + 
  geom_boxplot() +
  scale_fill_manual(values=c('darkgreen','darkolivegreen3')) +
  labs(x="Sampling times", y= "Relative abundance", title='Relative abundance dynamics of the selected phyllosphere OTUs (n=35)') +
  facet_wrap( ~ as.factor(final_names), scale = 'free_y') +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        plot.title = element_text(hjust = 0.5))

#Calculating relative abundance per groups based on taxonomy
ggplot(selected_otus[selected_otus$source=='phyllosphere',], aes(x = as.factor(sampling_date), y = abun, fill = plant)) + 
  geom_boxplot() +
  scale_fill_manual(values=c('darkgreen','darkolivegreen3')) +
  labs(x="Sampling times", y= "Relative abundance", title='Relative abundance dynamics of the selected phyllosphere OTUs taxonomically clustered to family (n=19)') +
  facet_wrap( ~ as.factor(Family), scale = 'free_y') +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        plot.title = element_text(hjust = 0.5))

# phyllo_outs <- otu_rare[,colnames(otu_rare) %in% selected_otus$sequence_name[selected_otus$source=='phyllosphere']]
# plant_otus <- data.frame(otu = as.factor(row.names(phyllo_outs)), phyllo_outs) %>% gather(sequence_name, abun, -otu) %>% 
#   left_join(map_16S[, c('sequence_name', 'plant', 'sampling_date', 'month')], by = 'sequence_name')
# head(plant_otus)

#Time dependent OTUs
#Who is there after the communities structures get less variable?
switch_unique_otus <- rel_otu_rare[rownames(rel_otu_rare) %in% swit_occ_abun$otu[swit_occ_abun$unique=='Switchgrass 2016'],]
switch_unique_late_stage <- data.frame(otu = as.factor(row.names(switch_unique_otus)), switch_unique_otus) %>% 
  gather(sequence_name, abun, -otu) %>% 
  left_join(map_16S[, c('sequence_name','rep','treatment' ,'source', 'plant', 'sampling_date', 'month', 'year'
  )], by = 'sequence_name') %>%
  left_join(tax_short, by='otu') %>% 
  filter(source == 'phyllosphere' & plant == 'switchgrass' & year == 2016) %>%
  group_by(sampling_date,plant,Kingdom,Phylum, Class,Order, Genus) %>%
  filter(abun>0) %>%
  summarise(n_otu=length(unique(otu))) %>%
  arrange(desc(n_otu))
sum(switch_unique_late_stage$n_otu) 

ggplot(switch_unique_late_stage, aes(x=as.factor(sampling_date), y=n_otu)) +
  theme_bw()+
  geom_bar(aes(fill=Phylum), stat='identity') +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        plot.title = element_text(hjust = 0.5))+
  labs(x='Date', y='Number of OTUs')

switch17_unique_otus <- rel_otu_rare[rownames(rel_otu_rare) %in% swit17_occ_abun$otu[swit17_occ_abun$unique=='Switchgrass 2017'],]
switch17_unique_late_stage <- data.frame(otu = as.factor(row.names(switch17_unique_otus)), switch17_unique_otus) %>% 
  gather(sequence_name, abun, -otu) %>% 
  left_join(map_16S[, c('sequence_name','rep','treatment' ,'source', 'plant', 'sampling_date', 'month', 'year'
  )], by = 'sequence_name') %>%
  left_join(tax_short, by='otu') %>% 
  filter(source == 'phyllosphere' & plant == 'switchgrass' & year == 2017) %>%
  group_by(sampling_date,plant,Kingdom,Phylum, Class,Order, Genus) %>%
  filter(abun>0) %>%
  summarise(n_otu=length(unique(otu))) %>%
  arrange(desc(n_otu))
sum(switch_unique_late_stage$n_otu) 

ggplot(switch17_unique_late_stage, aes(x=as.factor(sampling_date), y=n_otu)) +
  geom_bar(aes(fill=Phylum), stat='identity') +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        plot.title = element_text(hjust = 0.5))+
  labs(x='Date', y='Number of OTUs')

misc_unique_otus <- rel_otu_rare[rownames(rel_otu_rare) %in% misc_occ_abun$otu[misc_occ_abun$unique=='Miscanthus 2016'],]
dim(misc_unique_otus)
misc_unique_late_stage <- data.frame(otu = as.factor(row.names(misc_unique_otus)), misc_unique_otus) %>% 
  gather(sequence_name, abun, -otu) %>% 
  left_join(map_16S[, c('sequence_name','rep','treatment' ,'source', 'plant', 'sampling_date', 'month', 'year')], 
            by = 'sequence_name') %>%
  left_join(tax_short, by='otu') %>% 
  filter(source == 'phyllosphere' & plant == 'miscanthus' & year == 2016) %>%
  group_by(sampling_date,plant,Kingdom, Phylum, Class,Order, Genus) %>%
  filter(abun>0) %>%
  summarise(n_otu=length(unique(otu))) %>%
  arrange(desc(n_otu))
sum(misc_unique_late_stage$n_otu) 

ggplot(misc_unique_late_stage, aes(x=as.factor(sampling_date), y=n_otu)) +
  geom_bar(aes(fill=Phylum), stat='identity') +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        plot.title = element_text(hjust = 0.5))+
  labs(x='Date', y='Number of OTUs')

combined_phyllo_uniqes <- rbind(misc_unique_late_stage,switch_unique_late_stage)
combined_phyllo_uniqes <- rbind(combined_phyllo_uniqes,switch17_unique_late_stage)

ggplot(combined_phyllo_uniqes, aes(x=as.factor(sampling_date), y=n_otu)) +
  theme_bw()+
  geom_bar(aes(fill=Phylum), stat='identity') +
  facet_wrap(~plant) +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        plot.title = element_text(hjust = 0.5))+
  labs(x='Date', y='Number of OTUs')
ggsave('Figures/FigureS5.eps', height=4.5, width=7)

#overview of the unique taxa - taxonomy
unique_tax <- combined_phyllo_uniqes %>%
  group_by(Phylum, Order, Genus) %>%
  summarise(n=sum(n_otu)) %>%
  arrange(desc(n))

write.table(unique_tax, 'Figures/TableS1.txt', sep = '\t')


### End Nejc Analysis


##################################
##### Start Jackson Analysis #####
##################################
library(ggplot2)
library(scales)

#########################
#### Alpha Diversity ####
#########################
otu_rare.PA <- 1*(otu_rare>0)

s <- specnumber(otu_rare,MARGIN=2)
h <- diversity(t(otu_rare), index = "shannon")
pielou=h/log(s)
map_16S$sampling_date <- gsub(pattern = " EDT", replacement = "", map_16S$sampling_date)
map_16S$sampling_Rdate <- as.Date(map_16S$sampling_date)

otu_rare.phyllo.PA <- otu_rare.PA[,map_16S$source=="phyllosphere"]

sum(rowSums(otu_rare.phyllo.PA)>0)


### Settin up Contextual Data Maps
map.2017 <- map_16S[map_16S$year=="2017",]
map.2016 <- map_16S[map_16S$year=="2016",]

map.soil <- map_16S[map_16S$source=="soil",]
map.soil.2016 <- map_16S[map_16S$source=="soil"&map_16S$year==2016,]
map.soil.2017 <- map_16S[map_16S$source=="soil"&map_16S$year==2017,]

map.plant <- map_16S[map_16S$source=="phyllosphere",]
map.plant.2016 <- map.plant[map.plant$year==2016,]
map.plant.2017 <- map.plant[map.plant$year==2017,]

map.div <- map_16S
map.div$Richness <- s
map.div$Shannon <- h
map.div$Pielou <- pielou 
map.div.plant <- map.div[map.div$source=="phyllosphere",]
map.div.mis <- map.div.plant[map.div.plant$plant=="miscanthus",]
map.div.swg <- map.div.plant[map.div.plant$plant=="switchgrass",]

map.alpha <- melt(map.div, id.vars=c("sequence_name","treatment", "source", "plant", "time_numeric", "sampling_Rdate","year"), measure.vars=c("Richness", "Shannon", "Pielou"))

map.alpha.plant <- map.alpha[map.alpha$source=="phyllosphere",]

map.alpha.plant.2016 <- map.alpha[map.alpha$source=="phyllosphere"&map.alpha$year==2016,]

map.alpha.plant.2017 <- map.alpha[map.alpha$source=="phyllosphere"&map.alpha$year==2017,]

map.div.soil.2016 <- map.soil.2016
map.div.soil.2016$Richness <- ss.2016

map.alpha.soil <- melt(map.div.soil.2016, id.vars=c("sequence_name","treatment", "source", "plant", "time_numeric", "sampling_Rdate"), measure.vars=c("Richness"))


misc.map.2016 <- map.plant.2016[map.plant.2016$plant=="miscanthus",]
switch.map <- map.plant[map.plant$plant=="switchgrass",]
switch.map.2016 <- map.plant[map.plant$plant=="switchgrass"&map.plant$year==2016,]
switch.map.2017 <- map.plant[map.plant$plant=="switchgrass"&map.plant$year==2017,]

map.soil.misc <- map.soil[map.soil$plant=="miscanthus",]
map.soil.switch.2016 <- map.soil[map.soil$plant=="switchgrass"&map.soil$year==2016,]
map.soil.switch.2017 <- map.soil[map.soil$plant=="switchgrass"&map.soil$year==2017,]

### Soil rarefied datasets (down to 19,967 reads per sample)
otu_soil_rare.2016 <- otu_soil_rare[,map.soil$year==2016]
otu_soil_rare.2017 <- otu_soil_rare[,map.soil$year==2017]


ss <- specnumber(otu_soil_rare, MARGIN=2)
ss.2016 <- specnumber(otu_soil_rare.2016, MARGIN=2) 
ss.2017 <- specnumber(otu_soil_rare.2017, MARGIN=2)


cor.test(map.div.swg$time_numeric, map.div.swg$Richness)
cor.test(map.div.swg$time_numeric, map.div.swg$Shannon)
cor.test(map.div.swg$time_numeric, map.div.swg$Pielou)

cor.test(map.div.mis$time_numeric, map.div.mis$Richness)
cor.test(map.div.mis$time_numeric, map.div.mis$Shannon)
cor.test(map.div.mis$time_numeric, map.div.mis$Pielou)


plant_rare_otu <- otu_rare[,map_16S$source=="phyllosphere"]
plant_rare_otu.2016 <- otu_rare[,map_16S$source=="phyllosphere"&map_16S$year==2016]
plant_rare_otu.2017 <- otu_rare[,map_16S$source=="phyllosphere"&map_16S$year==2017]
soil_rare_otu.2016 <- otu_rare[,map_16S$source=="soil"&map_16S$year==2016]


### Fore variance partitioning of the core taxa
plant_rare_otu.2016.coretaxa <- plant_rare_otu.2016[core_taxa.df$x,]
plant_rare_otu.2016.coretaxa.2016 <- plant_rare_otu.2016[core_taxa_2016.df,]
plant_rare_otu.2016.coretaxa.rel <- decostand(plant_rare_otu.2016.coretaxa, MARGIN = 2, method = "total")
plant_rare_otu.2016.coretaxa.2016.rel <- decostand(plant_rare_otu.2016.coretaxa.2016, MARGIN=2, method="total")


plant_rare_otu.2016.coretaxa.2016.dist <- vegdist(t(plant_rare_otu.2016.coretaxa.2016), method="bray")

plant_rare_otu.2016.coretaxa.rel.dist <- vegdist(t(plant_rare_otu.2016.coretaxa.rel), method="bray")

plant_rare_otu.2016.coretaxa.2016.rel.dist <- vegdist(t(plant_rare_otu.2016.coretaxa.2016.rel), method="bray")


otu_rare.2017 <- otu_rare[,map_16S$year=="2017"]
colnames(otu_rare.2017) == map.2017$sequence_name

misc_rare_otu.2016 <- plant_rare_otu.2016[,map.plant.2016$plant=="miscanthus"]
switch_rare_otu <- plant_rare_otu[,map.plant$plant=="switchgrass"]
switch.rare.otu.2016 <- plant_rare_otu[,map.plant$plant=="switchgrass"&map.plant$year==2016]
switch.rare.otu.2017 <- plant_rare_otu[,map.plant$plant=="switchgrass"&map.plant$year==2017]

misc_unique_times <- unique(misc.map.2016$time_numeric)[order(unique(misc.map.2016$time_numeric))]

### Phyllosphere switchgrass unique times
switch_unique_times <- unique(switch.map$time_numeric)[order(unique(switch.map$time_numeric))]
switch_unique_times.2016 <- unique(switch.map.2016$time_numeric)[order(unique(switch.map.2016$time_numeric))]
switch_unique_times.2017 <-  unique(switch.map.2017$time_numeric)[order(unique(switch.map.2017$time_numeric))]

soil_misc_rare_otu.2016 <- soil_rare_otu.2016[,map.soil.2016$plant=="miscanthus"]
soil_switch_rare_otu.2016 <- soil_rare_otu.2016[,map.soil.2016$plant=="switchgrass"]

soil_misc_unique_date <- unique(map.soil.misc$sampling_Rdate)[order(unique(map.soil.misc$sampling_Rdate))]

soil_switch_unique_date.2016 <- unique(map.soil.switch.2016$sampling_Rdate)[order(unique(map.soil.switch.2016$sampling_Rdate))]
soil_switch_unique_date.2017 <- unique(map.soil.switch.2017$sampling_Rdate)[order(unique(map.soil.switch.2017$sampling_Rdate))]


# Look at species accumulation
misc_accumulation <- rep(1, length(misc_unique_times))
z <- NULL
for( i in 1:length(unique(misc.map.2016$time_numeric))){
  x <- misc_rare_otu.2016[,misc.map.2016$time_numeric==misc_unique_times[i]]
  z <- c(z, row.names(x[rowSums(x)!=0,]))
  misc_accumulation[i] <- length(unique(z))
}

soil_misc_accumulation <- rep(1, length(soil_misc_unique_date))
z <- NULL
for( i in 1:length(soil_misc_accumulation)){
  x <- soil_misc_rare_otu.2016[,map.soil.misc$sampling_Rdate==soil_misc_unique_date[i]]
  z <- c(z, row.names(x[rowSums(x)!=0,]))
  soil_misc_accumulation[i] <- length(unique(z))
}

switch_accumulation.2016 <- rep(1, length(switch_unique_times.2016))
z <- NULL
for( i in 1:length(unique(switch.map.2016$time_numeric))){
  x <- switch.rare.otu.2016[,switch.map.2016$time_numeric==switch_unique_times.2016[i]]
  z <- c(z, row.names(x[rowSums(x)!=0,]))
  switch_accumulation.2016[i] <- length(unique(z))
}

switch_accumulation.2017 <- rep(1, length(switch_unique_times.2017))
z <- NULL
for( i in 1:length(unique(switch.map.2017$time_numeric))){
  x <- switch.rare.otu.2017[,switch.map.2017$time_numeric==switch_unique_times.2017[i]]
  z <- c(z, row.names(x[rowSums(x)!=0,]))
  switch_accumulation.2017[i] <- length(unique(z))
}

switch_accumulation <- rep(1, length(switch_unique_times))
z <- NULL
for( i in 1:length(unique(switch.map$time_numeric))){
  x <- switch_rare_otu[,switch.map$time_numeric==switch_unique_times[i]]
  z <- c(z, row.names(x[rowSums(x)!=0,]))
  switch_accumulation[i] <- length(unique(z))
}

soil_switch_accumulation.2016 <- rep(1, length(soil_switch_unique_date.2016))
z <- NULL
for( i in 1:length(soil_switch_accumulation.2016)){
  x <- soil_switch_rare_otu.2016[,map.soil.switch.2016$sampling_Rdate==soil_switch_unique_date.2016[i]]
  z <- c(z, row.names(x[rowSums(x)!=0,]))
  soil_switch_accumulation.2016[i] <- length(unique(z))
}


Spec_accum.2016 <- NULL
Spec_accum.2016$Species_Accumulation <- c(misc_accumulation, switch_accumulation.2016)
Spec_accum.2016$Date <- c(unique(misc.map.2016$sampling_Rdate)[order(unique(misc.map.2016$sampling_Rdate))], unique(switch.map.2016$sampling_Rdate)[order(unique(switch.map.2016$sampling_Rdate))])
Spec_accum.2016 <- as.data.frame(Spec_accum.2016)
Spec_accum.2016$Plant <- c(rep("Miscanthus", 9), rep("Switchgrass", 8))
Spec_accum.2016

Spec_accum.2017 <- NULL
Spec_accum.2017$Species_Accumulation <- switch_accumulation.2017
Spec_accum.2017$Date <- unique(switch.map.2017$sampling_Rdate)[order(unique(switch.map.2017$sampling_Rdate))]
Spec_accum.2017 <- as.data.frame(Spec_accum.2017)
Spec_accum.2017$Plant <- rep("Switchgrass", nrow(Spec_accum.2017)) 

Spec_accum <- NULL
Spec_accum$Species_Accumulation <- switch_accumulation
Spec_accum$Date <- unique(switch.map$sampling_Rdate)[order(unique(switch.map$sampling_Rdate))]
Spec_accum <- as.data.frame(Spec_accum)
Spec_accum$Plant <- rep("Switchgrass", nrow(Spec_accum))

NewDates_Switch2017 <- c(as.Date("2016-05-15"), as.Date("2016-06-05"), as.Date("2016-06-26"), as.Date("2016-07-17"),as.Date("2016-08-07"), as.Date("2016-08-28"), as.Date("2016-09-18"))

Spec_accum_total <- NULL
Spec_accum_total$Species_Accumulation <- c(misc_accumulation, switch_accumulation.2016, switch_accumulation.2017)
Spec_accum_total$Date <- c(unique(misc.map.2016$sampling_Rdate)[order(unique(misc.map.2016$sampling_Rdate))], unique(switch.map.2016$sampling_Rdate)[order(unique(switch.map.2016$sampling_Rdate))], NewDates_Switch2017)
 Spec_accum_total <- as.data.frame(Spec_accum_total)
 Spec_accum_total$Plant <- c(rep("Miscanthus_2016", 9), rep("Switchgrass_2016", 8), rep("Switchgrass_2017", 7))

Spec_soil_accum <- NULL
Spec_soil_accum$Species_Accumulation <- c(soil_misc_accumulation, soil_switch_accumulation.2016)
Spec_soil_accum$Date <- c(soil_misc_unique_date, soil_switch_unique_date.2016)
Spec_soil_accum <- as.data.frame(Spec_soil_accum)
Spec_soil_accum$Plant <- c(rep("Miscanthus",10), rep("Switchgrass", 9))
Spec_soil_accum


#Plot Soil Miscanthus and Switchgrass Richness Through Season (rarefied to 76,488)
soil.misc.Richness <- map.alpha.soil[map.alpha.soil$plant=="miscanthus"&map.alpha.soil$variable=="Richness",]

A <- ggplot(soil.misc.Richness, aes(x=sampling_Rdate, y=value)) +
  geom_boxplot(mapping=aes(group=sampling_Rdate), fill="darkgreen") +
  labs(y="Richness", x="Date", title="Miscanthus") +
  scale_x_date(date_breaks = "1 month", labels=date_format("%b"), limits= as.Date(c('2016-04-01', '2016-12-01')))

soil.switch.Richness <- map.alpha.soil[map.alpha.soil$plant=="switchgrass"&map.alpha.soil$variable=="Richness",]

B <- ggplot(soil.switch.Richness, aes(x=sampling_Rdate, y=value)) +
  geom_boxplot(mapping=aes(group=sampling_Rdate), fill="darkolivegreen3") +
  labs(y="Richness", x="Date", title="Switchgrass")+
  scale_x_date(date_breaks = "1 month", labels=date_format("%b"), limits= as.Date(c('2016-04-01', '2016-12-01')))

C <- ggplot(Spec_soil_accum, aes(x=Date, y=Species_Accumulation, color=Plant)) + 
  geom_point() +
  geom_line()+
  labs(y="Total Observed Species", title ="Species Accumulation") + 
  scale_color_manual(values=c("darkgreen", "darkolivegreen3"))  + 
  guides(color=FALSE)+
  scale_x_date(date_breaks = "1 month", labels=date_format("%b"), limits= as.Date(c('2016-04-01', '2016-12-01')))

multiplot(A,B,C, cols=3)

# Code of Figure 2, Phyllosphere richness and accumulation through the season

### Figure 2A
switch.Richness.2016 <- map.alpha.plant.2016[map.alpha.plant.2016$plant=="switchgrass"&map.alpha.plant.2016$variable=="Richness",]
Fig2A <- ggplot(switch.Richness.2016, aes(x=sampling_Rdate, y=value)) +  
  geom_boxplot(mapping=aes(group=sampling_Rdate), fill="darkolivegreen3")+  
  #geom_point(color="darkolivegreen3") +
  theme_bw()+
  annotate("text", x = as.Date("2016-11-07"), y = 50, label = "NA")+
  labs(y="Richness", x="Date") + 
  scale_y_continuous(limits = c(0, 100))+
  scale_x_date(date_breaks = "1 month", labels=date_format("%b"), limits= as.Date(c('2016-05-01', '2016-11-15')))
#Are the samples different between each other?
fit_model <- aov(value~factor(sampling_Rdate),switch.Richness.2016)
summary(fit_model)
TukeyHSD(fit_model)

misc.Richness <- map.alpha.plant[map.alpha.plant$plant=="miscanthus"&map.alpha.plant$variable=="Richness",]
fit_model <- aov(value~factor(sampling_Rdate),misc.Richness)
summary(fit_model)
TukeyHSD(fit_model)

### Figure 2B
Fig2B <- ggplot(misc.Richness, aes(x=sampling_Rdate, y=value)) + 
  geom_boxplot(mapping=aes(group=sampling_Rdate), fill="darkgreen") + 
  #geom_point(color="darkgreen") + 
  theme_bw()+
  labs(y="Richness", x="Date") + 
  scale_y_continuous(limits = c(0, 100))+
  scale_x_date(date_breaks = "1 month", labels=date_format("%b"), limits= as.Date(c('2016-05-01', '2016-11-15'))) 

### Figure 2C
switch.Richness.2017 <- map.alpha.plant.2017[map.alpha.plant.2017$plant=="switchgrass"&map.alpha.plant.2017$variable=="Richness",]
Fig2C <- ggplot(switch.Richness.2017, aes(x=sampling_Rdate, y=value)) +  
  geom_boxplot(mapping=aes(group=sampling_Rdate), fill="darkolivegreen3")+  
  theme_bw()+
  #geom_point(color="darkolivegreen3") +
  annotate("text", x = as.Date("2017-11-07"), y = 50, label = "NA")+
  labs(y="Richness", x="Date") + 
  scale_y_continuous(limits = c(0, 100))+
  scale_x_date(date_breaks = "1 month", labels=date_format("%b"), limits= as.Date(c('2017-05-01', '2017-11-15')))


### Figure 2D 
Spec_accum_shapes.total <- c(rep(15,9), rep(16,15))
Spec_accum_lines.total <- c(rep(1,9), rep(2,8), rep(3,7))
Fig2D <- ggplot(Spec_accum_total, aes(x=Date, y=Species_Accumulation)) + 
  geom_point(aes(color=Plant, shape=Spec_accum_shapes.total)) +
  scale_shape_identity()+
  theme_bw()+
  scale_linetype_identity()+
  geom_line(aes(color=Plant, linetype=Spec_accum_lines.total))+
  labs(y="Total Observed Taxa") + 
  scale_color_manual(values=c("darkgreen","darkolivegreen3", "darkolivegreen3"))  + 
  guides(color=FALSE)+
  scale_y_continuous(limits = c(0,750))+
  scale_x_date(date_breaks = "1 month", labels=date_format("%b"), limits= as.Date(c('2016-05-01', '2016-11-15')))






setEPS()
postscript("Figure2.eps", width=6, height=4, paper="special")
par(ps = 8, cex = 1, cex.main = 1)
multiplot(Fig2A,Fig2C,Fig2B,Fig2D, cols=2)
dev.off()

# Plot Miscanthus Evenness Through Season
switch.Pielou.2016 <- map.alpha.plant.2016[map.alpha.plant.2016$plant=="switchgrass"&map.alpha.plant.2016$variable=="Pielou",]
A <- ggplot(switch.Pielou.2016, aes(x=sampling_Rdate, y=value)) +  
  geom_boxplot(mapping=aes(group=sampling_Rdate), fill="darkolivegreen3")+  
  #geom_point(color="darkolivegreen3") +
  annotate("text", x = as.Date("2016-11-07"), y = .5, label = "NA")+
  labs(y="Pielou's Evenness", x="Date", title ="A") + 
  scale_y_continuous(limits = c(0, 1))+
  scale_x_date(date_breaks = "1 month", labels=date_format("%b"), limits= as.Date(c('2016-05-01', '2016-11-15')))
A

switch.Pielou.2017 <- map.alpha.plant.2017[map.alpha.plant.2017$plant=="switchgrass"&map.alpha.plant.2017$variable=="Pielou",]
A_2017 <- ggplot(switch.Pielou.2017, aes(x=sampling_Rdate, y=value)) +  
  geom_boxplot(mapping=aes(group=sampling_Rdate), fill="darkolivegreen3")+  
  #geom_point(color="darkolivegreen3") +
  annotate("text", x = as.Date("2017-11-07"), y = .5, label = "NA")+
  labs(y="Pielou's Evenness", x="Date", title ="A") + 
  scale_y_continuous(limits = c(0, 1))+
  scale_x_date(date_breaks = "1 month", labels=date_format("%b"), limits= as.Date(c('2017-05-01', '2017-11-15')))

A_2017

misc.Pielou.2016 <- map.alpha.plant.2016[map.alpha.plant.2016$plant=="miscanthus"&map.alpha.plant.2016$variable=="Pielou",]

B <- ggplot(misc.Pielou.2016, aes(x=sampling_Rdate, y=value)) + 
  geom_boxplot(mapping=aes(group=sampling_Rdate), fill="darkgreen") + 
  #geom_point(color="darkgreen") + 
  labs(y="Pielou's Evenness", x="Date", title ="B") + 
  scale_y_continuous(limits = c(0, 1))+
  scale_x_date(date_breaks = "1 month", labels=date_format("%b"), limits= as.Date(c('2016-05-01', '2016-11-15'))) 

multiplot(A,B, cols=2)

setEPS()
postscript("FigureS3_Pielou.eps", width=7.3, height=3, paper="special")
par(ps = 8, cex = 1, cex.main = 1)
multiplot(A,B,cols=2)
dev.off()

t.test(map.alpha.plant.2016[map.alpha.plant.2016$variable=="Richness"&map.alpha.plant.2016$treatment=="standard fertilization",]$value, map.alpha.plant.2016[map.alpha.plant.2016$variable=="Richness"&map.alpha.plant.2016$treatment=="nitrogen free",]$value)


################################
#### Beta Diversity Analyses####
################################

NoMay_otu.2016 <- otu[,map_16S$sampling_date!="2016-05-31"&map_16S$sampling_date!="2016-05-09"&map_16S$year==2016]

NoLowSamples <- otu[,colSums(otu)>499]
NoLowSamples.146 <- rrarefy(t(NoLowSamples), 146)
NoLowSamples.500 <- rrarefy(t(NoLowSamples), 500)

NoLowSamples.146.dist <- vegdist(NoLowSamples.146, method="bray")
NoLowSamples.500.dist <- vegdist(NoLowSamples.500, method="bray")

mantel(NoLowSamples.146.dist, NoLowSamples.500.dist)

NoMay_rare_otu_500.2016 <- rrarefy(t(NoMay_otu.2016), 500)

NoMay_rare_otu_500.2016 <- t(NoMay_rare_otu_500.2016)

NoMay_map.2016 <- map_16S[map_16S$sampling_date!="2016-05-31"&map_16S$sampling_date!="2016-05-09"&map_16S$year==2016,]

NoMay_switch_map.2016 <- NoMay_map.2016[NoMay_map.2016$plant=="switchgrass",]



NoMay_rare_otu_146.2016 <- otu_rare[,map_16S$sampling_date!="2016-05-31"&map_16S$sampling_date!="2016-05-09"&map_16S$year==2016]

NoMay_switch_rare_otu_146.2016 <- NoMay_rare_otu_146.2016[,NoMay_map.2016$plant=="switchgrass"&NoMay_map.2016$source=="phyllosphere"]

NoMay_switch_rare_otu_500.2016 <- NoMay_rare_otu_500.2016[,NoMay_map.2016$plant=="switchgrass"&NoMay_map.2016$source=="phyllosphere"]


NoMay_switch_rare_146.2016.dist <- vegdist(t(NoMay_switch_rare_otu_146.2016), method = "bray")
NoMay_switch_rare_500.2016.dist <- vegdist(t(NoMay_switch_rare_otu_500.2016), method="bray")

mantel(NoMay_switch_rare_146.2016.dist, NoMay_switch_rare_500.2016.dist)




NoMay_rare_otu_500.2016.dist <- vegdist(t(NoMay_rare_otu_500.2016), method = "bray")

NoMay_rare_otu_146.2016.dist <- vegdist(t(NoMay_rare_otu_146.2016), method="bray")


NoMay_rare_otu_500.2016.pcoa <- cmdscale(NoMay_rare_otu_500.2016.dist, eig=TRUE)

NoMay_rare_otu_146.2016.pcoa <- cmdscale(NoMay_rare_otu_146.2016.dist, eig=TRUE)

plot(NoMay_rare_otu_146.2016.pcoa$points[,1], NoMay_rare_otu_146.2016.pcoa$points[,2])

plot(NoMay_rare_otu_500.2016.pcoa$points[,1], NoMay_rare_otu_500.2016.pcoa$points[,2])

NoMay_Colors <- rep("darkolivegreen3", nrow(NoMay_map.2016))
NoMay_Colors[NoMay_map.2016$source=="phyllosphere"&NoMay_map.2016$plant=="miscanthus"] <- "darkgreen"
NoMay_Colors[NoMay_map.2016$source=="soil"&NoMay_map.2016$plant=="miscanthus"] <- "burlywood4"
NoMay_Colors[NoMay_map.2016$source=="soil"&NoMay_map.2016$plant=="switchgrass"] <- "burlywood"


whole_plot_shapes <- rep(15, nrow(map.whole.2016))
whole_plot_shapes[map.whole.2016$plant=="switchgrass"] <- 16
whole_plot_shapes[map.whole.2016$plant=="switchgrass"&map.whole.2016$source=="soil"&map.whole.2016$treatment=="nitrogen free"] <- 1 

whole_plot_shapes[map.whole.2016$plant=="miscanthus"&map.whole.2016$source=="soil"&map.whole.2016$treatment=="nitrogen free"] <- 0 

NoMay_Symbols <- rep(15, nrow(NoMay_map.2016))
NoMay_Symbols[NoMay_map.2016$plant=="switchgrass"] <- 16
NoMay_Symbols[NoMay_map.2016$plant=="switchgrass"&NoMay_map.2016$source=="soil"&NoMay_map.2016$treatment=="nitrogen free"] <- 1
NoMay_Symbols[NoMay_map.2016$plant=="miscanthus"&NoMay_map.2016$source=="soil"&NoMay_map.2016$treatment=="nitrogen free"] <- 0

ax1.146 <- NoMay_rare_otu_146.2016.pcoa$eig[1]/sum(NoMay_rare_otu_146.2016.pcoa$eig)
ax2.146 <- NoMay_rare_otu_146.2016.pcoa$eig[2]/sum(NoMay_rare_otu_146.2016.pcoa$eig)

ax1.500 <- NoMay_rare_otu_500.2016.pcoa$eig[1]/sum(NoMay_rare_otu_500.2016.pcoa$eig)
ax2.500 <- NoMay_rare_otu_500.2016.pcoa$eig[2]/sum(NoMay_rare_otu_500.2016.pcoa$eig)


setEPS()
postscript(file="SupplementaryPCoA_RarefactionComparison.ps",  width=8, height=4, paper="special" )
par(mfrow=c(1,2))
plot(NoMay_rare_otu_146.2016.pcoa$points[,1], NoMay_rare_otu_146.2016.pcoa$points[,2], pch=NoMay_Symbols, col=NoMay_Colors, main="Rarefied to 146 Sequences",xlab=paste("PCoA1: ",100*round(ax1.146,3),"% var. explained",sep=""), ylab=paste("PCoA2: ",100*round(ax2.146,3),"% var. explained",sep=""))

plot(NoMay_rare_otu_500.2016.pcoa$points[,1], NoMay_rare_otu_500.2016.pcoa$points[,2], pch=NoMay_Symbols, col=NoMay_Colors, main="Rarefied to 500 Sequences",xlab=paste("PCoA1: ",100*round(ax1.500,3),"% var. explained",sep=""), ylab=paste("PCoA2: ",100*round(ax2.500,3),"% var. explained",sep=""))
dev.off()



whole.dist <- vegdist(t(otu_rare), method="bray")
whole.dist.2016 <- vegdist(t(otu_rare[,map_16S$year=="2016"]), method="bray")

plant.dist <- vegdist(t(otu_rare[,map_16S$source=="phyllosphere"]), method="bray")

plant.dist.2016 <- vegdist(t(otu_rare[,map_16S$source=="phyllosphere"&map_16S$year=="2016"]), method="bray")
miscanthus.dist.2016 <- vegdist(t(otu_rare[,map_16S$source=="phyllosphere"&map_16S$plant=="miscanthus"&map_16S$year=="2016"]), method="bray")


switch.dist <- vegdist(t(otu_rare[,map_16S$source=="phyllosphere"&map_16S$plant=="switchgrass"]), method="bray")

switch.dist.2016 <- vegdist(t(otu_rare[,map_16S$source=="phyllosphere"&map_16S$plant=="switchgrass"&map_16S$year=="2016"]), method="bray")

switch.dist.2017 <- vegdist(t(otu_rare[,map_16S$source=="phyllosphere"&map_16S$plant=="switchgrass"&map_16S$year=="2017"]), method="bray")

map.whole.2016 <- map_16S[map_16S$year==2016,]




### Make Soil maps, otu table and dissimilarity matrices for each year
map.soil.2016 <- map.soil[map.soil$year==2016,]
map.soil.2017 <- map.soil[map.soil$year==2017,]

otu_soil_rare.2016 <- otu_soil_rare[,map.soil$year==2016]
otu_soil_rare.2017 <- otu_soil_rare[,map.soil$year==2017]


soil.dist <- vegdist(t(otu_soil_rare), method="bray")
soil.dist.2016 <- vegdist(t(otu_soil_rare.2016), method="bray")
soil.dist.2017 <- vegdist(t(otu_soil_rare.2017), method="bray")


###################
### For PROTEST ###
###################


switch.otu.2016 <- otu_rare[,map_16S$source=="phyllosphere"&map_16S$plant=="switchgrass"&map_16S$year==2016&map_16S$month!=10]

switch.otu.2017 <- otu_rare[,map_16S$source=="phyllosphere"&map_16S$plant=="switchgrass"&map_16S$year==2017]

switch.map.2016.p <- map_16S[map_16S$source=="phyllosphere"&map_16S$plant=="switchgrass"&map_16S$year==2016&map_16S$month!=10,]
switch.map.2017.p <- map_16S[map_16S$source=="phyllosphere"&map_16S$plant=="switchgrass"&map_16S$year==2017,]

switch.map.2016.p$timepoint <- rep(1,nrow(switch.map.2016.p))
switch.map.2016.p$timepoint[switch.map.2016.p$sampling_date=="2016-05-31"]<- 2
switch.map.2016.p$timepoint[switch.map.2016.p$sampling_date=="2016-06-20"]<- 3
switch.map.2016.p$timepoint[switch.map.2016.p$sampling_date=="2016-07-12"]<- 4
switch.map.2016.p$timepoint[switch.map.2016.p$sampling_date=="2016-08-01"]<- 5
switch.map.2016.p$timepoint[switch.map.2016.p$sampling_date=="2016-08-22"]<- 6
switch.map.2016.p$timepoint[switch.map.2016.p$sampling_date=="2016-09-12"]<- 7


switch.map.2017.p$timepoint <- rep(1,nrow(switch.map.2017))
switch.map.2017.p$timepoint[switch.map.2017.p$sampling_date=="2017-06-05"] <- 2
switch.map.2017.p$timepoint[switch.map.2017.p$sampling_date=="2017-06-26"] <- 3
switch.map.2017.p$timepoint[switch.map.2017.p$sampling_date=="2017-07-17"] <- 4
switch.map.2017.p$timepoint[switch.map.2017.p$sampling_date=="2017-08-07"] <- 5
switch.map.2017.p$timepoint[switch.map.2017.p$sampling_date=="2017-08-28"] <- 6
switch.map.2017.p$timepoint[switch.map.2017.p$sampling_date=="2017-09-18"] <- 7

switch.collapsed.otu.2016 <- data.frame(Time1=rep(NA, nrow(switch.otu.2016)), Time2=rep(NA, nrow(switch.otu.2016)), Time3=rep(NA, nrow(switch.otu.2016)), Time4=rep(NA, nrow(switch.otu.2016)), Time5=rep(NA, nrow(switch.otu.2016)), Time6=rep(NA, nrow(switch.otu.2016)), Time7=rep(NA, nrow(switch.otu.2016)))

switch.collapsed.otu.2017 <- data.frame(Time1=rep(NA, nrow(switch.otu.2017)), Time2=rep(NA, nrow(switch.otu.2017)), Time3=rep(NA, nrow(switch.otu.2017)), Time4=rep(NA, nrow(switch.otu.2017)), Time5=rep(NA, nrow(switch.otu.2017)), Time6=rep(NA, nrow(switch.otu.2017)), Time7=rep(NA, nrow(switch.otu.2017)))


for (i in 1:7){
  temp_2016 <- switch.otu.2016[,switch.map.2016.p$timepoint==i]
  temp_2017 <- switch.otu.2017[, switch.map.2017.p$timepoint==i]
  switch.collapsed.otu.2016[,i] <- rowSums(temp_2016)/ncol(temp_2016)
  switch.collapsed.otu.2017[,i] <- rowSums(temp_2017)/ncol(temp_2017)
}

temp.2016 <- switch.collapsed.otu.2016
temp.2017 <- switch.collapsed.otu.2017

temp.2016$OTU <- row.names(switch.collapsed.otu.2016)
temp.2017$OTU <- row.names(switch.collapsed.otu.2017)

switch.collapsed.otu.matched <- inner_join(temp.2016, temp.2017, by="OTU")

switch.collapsed.otu.matched$OTU <- NULL

switch.collapsed.dist.2016 <- vegdist(t(switch.collapsed.otu.2016), method="bray")
switch.collapsed.dist.2017 <- vegdist(t(switch.collapsed.otu.2017), method="bray")
switch.collapsed.dist.both <- vegdist(t(switch.collapsed.otu.matched), method="bray")


switch.collapsed.pcoa.2016 <- cmdscale(switch.collapsed.dist.2016, eig=TRUE)
switch.collapsed.pcoa.2017 <- cmdscale(switch.collapsed.dist.2017, eig=TRUE)

protest(switch.collapsed.pcoa.2016, switch.collapsed.pcoa.2017)
protest(switch.collapsed.dist.2016, switch.collapsed.dist.2017)

adonis(switch.collapsed.dist.both~c(rep(2016, 7), rep(2017,7)))


switch.phyllo.map <- map_16S[map_16S$plant=="switchgrass"&map_16S$source=="phyllosphere",]

switch.dist.2017 <- vegdist(t(switch.otu.2017), method = "bray")
switch.pcoa.2017 <- cmdscale(switch.dist.2017, eig = TRUE)


plot(switch.pcoa.2017$points[,1], switch.pcoa.2017$points[,2], cex=0)
text(switch.pcoa.2017$points[,1], switch.pcoa.2017$points[,2], labels = switch.map.2017$month)

# Making matching pcoas for miscanthus and switchgrass, for use in PROTEST
otu_miscanthus_sub <- otu_rare[,map_16S$source=="phyllosphere"&map_16S$plant=="miscanthus"&map_16S$year=="2016"] 
otu_switchgrass_sub <- otu_rare[,map_16S$source=="phyllosphere"&map_16S$plant=="switchgrass"&map_16S$year=="2016"]

map_miscanthus_sub <- map_16S[map_16S$source=="phyllosphere"&map_16S$plant=="miscanthus"&map_16S$year=="2016",]
map_switchgrass_sub <- map_16S[map_16S$source=="phyllosphere"&map_16S$plant=="switchgrass"&map_16S$year=="2016",]

otu_miscanthus_matched_dates <- otu_miscanthus_sub[,map_miscanthus_sub$sampling_date%in%map_switchgrass_sub$sampling_date]
otu_switchgrass_matched_dates <- otu_switchgrass_sub[,map_switchgrass_sub$sampling_date%in%map_switchgrass_sub$sampling_date]

map_miscanthus_matched_dates <- map_miscanthus_sub[map_miscanthus_sub$sampling_date%in%map_switchgrass_sub$sampling_date,]
map_switchgrass_matched_dates <- map_switchgrass_sub[map_switchgrass_sub$sampling_date%in%map_switchgrass_sub$sampling_date,]


matched_dates <- unique(map_miscanthus_matched_dates$sampling_date)
day_misc <- NULL
day_switch <- NULL
for(i in 1:length(matched_dates)){
  temp_misc <- otu_miscanthus_matched_dates[,map_miscanthus_matched_dates$sampling_date==matched_dates[i]]
  temp_switch <- otu_switchgrass_matched_dates[,map_switchgrass_matched_dates$sampling_date==matched_dates[i]]
  day_misc <- rbind(day_misc,(rowSums(temp_misc)/ncol(temp_misc)))
  day_switch <- rbind(day_switch,(rowSums(temp_switch)/ncol(temp_switch)))
}
day_misc <- t(day_misc)
day_switch<-(t(day_switch))
row.names(day_misc) <- row.names(otu_miscanthus_matched_dates)
row.names(day_switch) <- row.names(otu_switchgrass_matched_dates)
colnames(day_misc) <- matched_dates
colnames(day_switch) <- matched_dates

misc_matched_dates.dist <- vegdist(t(day_misc), method="bray")
switch_matched_dates.dist <- vegdist(t(day_switch), method="bray")

misc_matched_dates.pcoa <- cmdscale(misc_matched_dates.dist, eig=TRUE)
switch_matched_dates.pcoa <- cmdscale(switch_matched_dates.dist, eig=TRUE)

protest(misc_matched_dates.pcoa, switch_matched_dates.pcoa)

########################################
#### Beta Dispersion of Phyllosphere####
########################################
misc_disper <-betadisper(miscanthus.dist.2016, group=misc.map.2016$sampling_date)

switch.disper.2016 <- betadisper(switch.dist.2016, group=switch.map.2016$sampling_date)

switch.disper.2017 <- betadisper(switch.dist.2017, group=switch.map.2017$sampling_date)

switch.disper.dates.2017 <- data.frame(distance=switch.disper.2017$distances, dates=as.Date(switch.map.2017$sampling_date))

misc_disperion_dates <- data.frame(distance=misc_disper$distances, dates=as.Date(misc.map.2016$sampling_date), sample=misc.map.2016$sequence_name)

switch_dispersion_dates.2016 <- data.frame(distance=switch.disper.2016$distances, dates=as.Date(switch.map.2016$sampling_date), sample=switch.map.2016$sampling_date)


A<- ggplot(switch_dispersion_dates.2016, aes(x=dates, y=distance))+
  geom_boxplot(aes(group=dates),fill="darkolivegreen3")+
  annotate("text", x = as.Date("2016-11-07"), y = 0.25, label = "NA")+
  labs(y="Distance to Median", x="Date", title="A") + 
  scale_y_continuous(limits = c(0, .5))+
  scale_x_date(date_breaks = "1 month", labels=date_format("%b"), limits= as.Date(c('2016-05-01', '2016-11-15')))
B<- ggplot(misc_disperion_dates, aes(x=dates, y=distance)) +
  geom_boxplot(aes(group=dates),fill="darkgreen")+
  labs(y="Distance to Median", x="Date", title="B") + 
  scale_y_continuous(limits = c(0, .5))+
  scale_x_date(date_breaks = "1 month", labels=date_format("%b"), limits= as.Date(c('2016-05-01', '2016-11-15'))) 
  


A_2017 <- ggplot(switch.disper.dates.2017, aes(x=dates, y=distance))+
  geom_boxplot(aes(group=dates),fill="darkolivegreen3")+
  annotate("text", x = as.Date("2017-11-07"), y = 0.25, label = "NA")+
  labs(y="Distance to Median", x="Date", title="A_2017") + 
  scale_y_continuous(limits = c(0, .5))+
  scale_x_date(date_breaks = "1 month", labels=date_format("%b"), limits= as.Date(c('2017-05-01', '2017-11-15')))
 
FigureS3 <- multiplot(A,A_2017,B, cols=2)
setEPS()
postscript("FigureS3_BetaDispersion.eps", width=6, height=4, paper="special")
par(ps = 8, cex = 1, cex.main = 1)
multiplot(A,A_2017,B, cols=2)
dev.off()

#Correlations between average distance to median and time for switchgrass and miscanthus
switch_dispersion_avg.2016 <- switch_dispersion_dates.2016 %>%
  group_by(dates) %>%
  summarise(mean = mean(distance), n = n())

unique(switch.map.2016$time_numeric)[order(unique(switch.map.2016$time_numeric))]
cor.test(switch_dispersion_avg.2016$mean, unique(switch.map.2016$time_numeric)[order(unique(switch.map.2016$time_numeric))]  )


misc_dispersion_avg <- misc_disperion_dates %>%
  group_by(dates) %>%
  summarise(mean = mean(distance), n = n())

unique(misc.map.2016$time_numeric)[order(unique(misc.map.2016$time_numeric))]
cor.test(misc_dispersion_avg$mean, unique(misc.map.2016$time_numeric)[order(unique(misc.map.2016$time_numeric))])



#########################################################
#### Plotting PCoA's for Soil, both and phyllosphere ####
#########################################################



############################
### Figure 1B Whole PCoA ###
############################
map_16S$timepoint <- rep(1,nrow(map_16S))
map_16S$timepoint[map_16S$sampling_date=="2016-05-31"]<- 2
map_16S$timepoint[map_16S$sampling_date=="2017-06-05"]<- 2
map_16S$timepoint[map_16S$sampling_date=="2016-06-20"]<- 3
map_16S$timepoint[map_16S$sampling_date=="2017-06-26"]<- 3
map_16S$timepoint[map_16S$sampling_date=="2016-07-12"] <- 4
map_16S$timepoint[map_16S$sampling_date=="2017-07-17"]<- 4
map_16S$timepoint[map_16S$sampling_date=="2016-08-01"] <- 5
map_16S$timepoint[map_16S$sampling_date=="2017-08-07"] <- 5
map_16S$timepoint[map_16S$sampling_date=="2016-08-22"] <- 6
map_16S$timepoint[map_16S$sampling_date=="2017-08-28"]<- 6
map_16S$timepoint[map_16S$sampling_date=="2016-09-12"] <- 7
map_16S$timepoint[map_16S$sampling_date=="2017-09-18"]<- 7
map_16S$timepoint[map_16S$sampling_date=="2016-10-03"]<- 8
map_16S$timepoint[map_16S$sampling_date=="2016-11-07"]<- 9


whole.pcoa.2016 <- cmdscale(whole.dist.2016, eig=TRUE)
ax1.whole.2016 <- whole.pcoa.2016$eig[1]/sum(whole.pcoa.2016$eig)
ax2.whole.2016 <- whole.pcoa.2016$eig[2]/sum(whole.pcoa.2016$eig)

whole_plot_colors <- rep("darkolivegreen3", nrow(map.whole.2016))
whole_plot_colors[map.whole.2016$plant=="miscanthus"] <- "darkgreen"
whole_plot_colors[map.whole.2016$source=="soil"] <- "burlywood4"
whole_plot_colors[map.whole.2016$source=="soil"&map.whole.2016$plant=="switchgrass"] <- "burlywood"
whole_plot_shapes <- rep(15, nrow(map.whole.2016))
whole_plot_shapes[map.whole.2016$plant=="switchgrass"] <- 16
whole_plot_shapes[map.whole.2016$plant=="switchgrass"&map.whole.2016$source=="soil"&map.whole.2016$treatment=="nitrogen free"] <- 1 

whole_plot_shapes[map.whole.2016$plant=="miscanthus"&map.whole.2016$source=="soil"&map.whole.2016$treatment=="nitrogen free"] <- 0 

Unique_Dates <- unique(map.whole.2016$sampling_date)[order(unique(map.whole.2016$sampling_date))]
Unique_Dates.2017 <- unique(map.plant.2017$sampling_date)[order(unique(map.plant.2017$sampling_date))]

Date_Size <- rep(1, nrow(map.whole.2016))
point_sizes_whole <- seq(from=1, to=3, length.out = length(Unique_Dates))
point.sizes_whole.2017 <- seq(from=1, to=3, length.out=length(unique(map.plant.2017$sampling_Rdate)))
for(i in 1:length(Unique_Dates)){
  Date_Size[map_16S$sampling_date==Unique_Dates[i]] <- point_sizes_whole[i]
}
unique(Date_Size)
U_Date_Size <- unique(Date_Size)[order(unique(Date_Size))]

setEPS() 
postscript("Figure1B_Whole_Size.eps", width=5, height=5,pointsize=10, paper="special",)
par(mar = c(5.1, 5.1, 2.1, 2.1))
plot(whole.pcoa.2016$points[,1], whole.pcoa.2016$points[,2],cex=Date_Size, cex.axis=1.6,cex.lab=1.6,  pch= whole_plot_shapes, col=whole_plot_colors, xlab=paste("PCoA1: ",100*round(ax1.whole.2016,3),"% var. explained",sep=""), ylab=paste("PCoA2: ",100*round(ax2.whole.2016,3),"% var. explained",sep=""), main="Soil and Phyllosphere PCoA")#+
#legend("bottomleft",legend = unique(map_16S$sampling_Rdate)[order(unique(map_16S$sampling_Rdate))], pt.cex=U_Date_Size, pch=c(rep(1,length(U_Date_Size)), rep(0,length(U_Date_Size))))
#legend("bottomright", legend=unique(map_16S$sampling_Rdate)[order(unique(map_16S$sampling_Rdate))], pch=rep(0,length(U_Date_Size)), pt.cex=U_Date_Size)
dev.off()

setEPS() 
postscript("Figure1B_CircleLegend.eps", width=5, height=5,pointsize=10, paper="special",)
par(mar = c(5.1, 5.1, 2.1, 2.1))
plot(whole.pcoa.2016$points[,1], whole.pcoa.2016$points[,2],cex=Date_Size, cex.axis=1.6,cex.lab=1.6,  pch= whole_plot_shapes, col=whole_plot_colors, xlab=paste("PCoA1: ",100*round(ax1.whole.2016,3),"% var. explained",sep=""), ylab=paste("PCoA2: ",100*round(ax2.whole.2016,3),"% var. explained",sep=""), main="Soil and Phyllosphere PCoA")#+
legend("bottomleft",legend = unique(map.whole.2016$sampling_Rdate)[order(unique(map.whole.2016$sampling_Rdate))], pt.cex=U_Date_Size, pch=c(rep(1,length(U_Date_Size)), rep(0,length(U_Date_Size))), cex=1.3)
#legend("top", legend=unique(map_16S$sampling_Rdate)[order(unique(map_16S$sampling_Rdate))], pch=rep(0,length(U_Date_Size)), pt.cex=U_Date_Size, cex=1.3)
dev.off()

setEPS() 
postscript("Figure1B_SquareLegend.eps", width=5, height=5,pointsize=10, paper="special",)
par(mar = c(5.1, 5.1, 2.1, 2.1))
plot(whole.pcoa.2016$points[,1], whole.pcoa.2016$points[,2],cex=Date_Size, cex.axis=1.6,cex.lab=1.6,  pch= whole_plot_shapes, col=whole_plot_colors, xlab=paste("PCoA1: ",100*round(ax1.whole.2016,3),"% var. explained",sep=""), ylab=paste("PCoA2: ",100*round(ax2.whole.2016,3),"% var. explained",sep=""), main="Soil and Phyllosphere PCoA")#+
#legend("bottomleft",legend = unique(map_16S$sampling_Rdate)[order(unique(map_16S$sampling_Rdate))], pt.cex=U_Date_Size, pch=c(rep(1,length(U_Date_Size)), rep(0,length(U_Date_Size))), cex=1.3)
legend("bottomleft", legend=unique(map.whole.2016$sampling_Rdate)[order(unique(map.whole.2016$sampling_Rdate))], pch=rep(0,length(U_Date_Size)), pt.cex=U_Date_Size, cex=1.3)
dev.off()




############################
### Figure 1A Plant PCoA ###
############################


plant.pcoa <- cmdscale(plant.dist, eig=TRUE)

switch.pcoa <- cmdscale(switch.dist, eig=TRUE)

plant.pcoa.2016 <- cmdscale(plant.dist.2016, eig=TRUE)
ax1.plant.2016 <- plant.pcoa.2016$eig[1]/sum(plant.pcoa.2016$eig)
ax2.plant.2016 <- plant.pcoa.2016$eig[2]/sum(plant.pcoa.2016$eig)
ax1.plant <- plant.pcoa$eig[1]/sum(plant.pcoa$eig)
ax2.plant <- plant.pcoa$eig[2]/sum(plant.pcoa$eig)

ax1.switch <- switch.pcoa$eig[1]/sum(switch.pcoa$eig)
ax2.switch <- switch.pcoa$eig[2]/sum(switch.pcoa$eig)

### Make collapsed PCoA Coordinates

Sample_Dates <- unique(map.plant$sampling_Rdate)
Sample_Dates <- Sample_Dates[order(Sample_Dates)]
Sample_Dates.2016 <- unique(map.plant.2016$sampling_Rdate)
plant_coordinates.2016 <- plant.pcoa.2016$points
plant_coordinates <- plant.pcoa$points

switch_coordinates <- switch.pcoa$points
Sample_Dates.switch <- unique(switch.map$sampling_Rdate)[order(unique(switch.map$sampling_Rdate))]  
  
misc_collapse.2016 <- data.frame(ax1.average=rep(NA,9),ax2.average=rep(NA,9), ax1.sd=rep(NA,9), ax2.sd=rep(NA,9) )
switch_collapse.2016 <- data.frame(ax1.average=rep(NA,9),ax2.average=rep(NA,9), ax1.sd=rep(NA,9), ax2.sd=rep(NA,9) )
switch_collapse.2017 <- data.frame(ax1.average=rep(NA,9), ax2.average=rep(NA,9), ax1.sd=rep(NA,9), ax2.sd=rep(NA,9))

switch_collapse <- data.frame(ax1.average=rep(NA,15),ax2.average=rep(NA,15), ax1.sd=rep(NA,15), ax2.sd=rep(NA,15) )

for (i in 1:length(unique(map.plant.2016$sampling_Rdate))){
  b <- plant_coordinates.2016[map.plant.2016$sampling_Rdate==Sample_Dates.2016[i],]
  dm <- map.plant[map.plant.2016$sampling_Rdate==Sample_Dates.2016[i],]
  m <- b[dm$plant=="miscanthus",]
  s <- b[dm$plant=="switchgrass",]
  misc_collapse.2016[i,1] <-mean(m[,1])
  misc_collapse.2016[i,2] <- mean(m[,2])
  misc_collapse.2016[i,3] <-sd(m[,1])
  misc_collapse.2016[i,4] <-sd(m[,1])
  switch_collapse.2016[i,1] <-mean(s[,1])
  switch_collapse.2016[i,2] <- mean(s[,2])
  switch_collapse.2016[i,3] <-sd(s[,1])
  switch_collapse.2016[i,4] <-sd(s[,1])
}

switch.2016.point_size <- seq(from=1, to=3, length.out = 8)
switch.2017.point_size <- seq(from=1, to=3, length.out = 7)

for (i in 1:length(unique(switch.map$sampling_Rdate))){
  b <- switch_coordinates[switch.map$sampling_Rdate==Sample_Dates.switch[i],]
  switch_collapse[i,1] <-mean(b[,1])
  switch_collapse[i,2] <- mean(b[,2])
  switch_collapse[i,3] <-sd(b[,1])
  switch_collapse[i,4] <-sd(b[,1])
}
point_sizes_switch.collapse <- c(switch.2016.point_size, switch.2017.point_size) 
point_colors_switch.collapse <- c(rep("darkolivegreen3", 8), rep("darkolivegreen1",7))

switch.weather <- switch.map[,c("precipitation", "Air_temp_mean", "air_temp_max", "Air_Temp_Min", "Air_Pressure", "RH", "AH", "Wind_Speed_Mean", "Solar_Radiation", "PAR", "soil_temp_5_cm_bare_avg", "sampling_Rdate", "time_numeric")]
switch.weather <- unique(switch.weather)

switch.leaf.chemistry <- switch.map[,c("LDMC_mg_per_g", "nitrogen_percent", "carbon_percent", "carbon_per_nitrogen", "height_mean_cm")]
switch.soil.chemsitry <- switch.map[,c("pH", "P_ppm", "K_ppm", "Ca_ppm", "Mg_ppm", "organic_matter", "NO3N_ppm", "NH4_ppm", "soil_moisture_percent", "soil_temp_10cm")]
for (i in 1:ncol(switch.soil.chemsitry)){
  switch.soil.chemsitry[,i] <- as.numeric(switch.soil.chemsitry[,i])
}
switch.LC.envfit <- envfit(switch.pcoa, switch.leaf.chemistry)
switch.LC.envfit.R40 <- envfit(switch.pcoa, switch.leaf.chemistry[,c("LDMC_mg_per_g", "nitrogen_percent", "carbon_per_nitrogen", "height_mean_cm")])

switch.SC.envfit <- envfit(switch.pcoa, switch.soil.chemsitry)
switch.SC.envfit.R40 <- envfit(switch.pcoa, switch.soil.chemsitry[,c("soil_moisture_percent")])
switch.weather.envfit <- envfit(switch_collapse[,1:2], switch.weather)

switch.LC.envfit$vectors


setEPS()
postscript("Switch_Phyllo_PCoa.eps", width=5, height=5,pointsize=10, paper="special",)
par(mar = c(5.1, 5.1, 2.1, 2.1))
plot(switch_collapse$ax1.average, switch_collapse$ax2.average, cex=point_sizes_switch.collapse, col=point_colors_switch.collapse, pch=16, xlim=c(-.5,.4), ylim=c(-.25,.4), xlab=paste("PCoA1: ",100*round(ax1.switch,3),"% var. explained",sep=""), ylab=paste("PCoA2: ",100*round(ax2.switch,3),"% var. explained",sep=""))+
arrows(switch_collapse$ax1.average, switch_collapse$ax2.average-switch_collapse$ax1.sd, switch_collapse$ax1.average, switch_collapse$ax2.average+switch_collapse$ax2.sd, angle=90,length= 0.05, code=3,lwd=0.75, col=point_colors_switch.collapse)+
arrows(switch_collapse$ax1.average- switch_collapse$ax1.sd, switch_collapse$ax2.average, switch_collapse$ax1.average + switch_collapse$ax1.sd, switch_collapse$ax2.average, angle=90,length= 0.05, code=3,lwd=0.75, col= point_colors_switch.collapse)+
  plot(switch.LC.envfit.R40, p=0.05, col="black")+
  plot(switch.SC.envfit.R40, p=0.05, col="black", labels="soil_moisture")+
  plot(switch.weather.envfit, p=0.05, col="black")
dev.off()

SwitchData_Pcoa <- data.frame(Ax1=switch_collapse$ax1.average, Ax2=switch_collapse$ax2.average, Year= factor(c(rep(2016,8), rep(2017,7)), levels=c(2016,2017)), PointSize =c(switch.2016.point_size, switch.2017.point_size))


switch.LC.envfit.R40
switch.LC.envfit.R40.df<-as.data.frame(switch.LC.envfit.R40$vectors$arrows*sqrt(switch.LC.envfit.R40$vectors$r))
switch.LC.envfit.R40.df$species<-rownames(switch.LC.envfit.R40.df)

spp.scrs <- as.data.frame(scores(switch.LC.envfit.R40, display = "vectors"))
spp.scrs <- cbind(spp.scrs, Species = rownames(spp.scrs))


vegan:::ordiArrowMul(switch.LC.envfit.R40)

##### Trying to figure out ggplot EnvFit
ggplot(SwitchData_Pcoa, aes(x=Ax1, y=Ax2)) +
  geom_point(size=SwitchData_Pcoa$PointSize*2, aes(color=Year)) +
  scale_color_manual(values =c("darkolivegreen3", "darkolivegreen1")) +
  coord_fixed()+
  geom_segment(data=switch.LC.envfit.R40.df,aes(x=0,xend=Dim1,y=0,yend=Dim2))
  
ggplot(SwitchData_Pcoa, aes(x=Ax1, y=Ax2)) +
  geom_point(size=SwitchData_Pcoa$PointSize*2, aes(color=Year)) +
  scale_color_manual(values =c("darkolivegreen3", "darkolivegreen1")) +
  coord_fixed()+
  geom_segment(data=spp.scrs, aes(x=0,xend=Dim1*0.4075751,y=0,yend=Dim2*0.4075751))



switch.whole.map <- map_16S[map_16S$plant=="switchgrass",]
switch.whole.otu <- otu_rare[,map_16S$plant=="switchgrass"]
switch.whole.dist <- vegdist(t(switch.whole.otu), method = "bray")
switch.whole.pcoa <- cmdscale(switch.whole.dist, eig=TRUE)

ax1.switch.whole <- switch.whole.pcoa$eig[1]/sum(switch.whole.pcoa$eig)
ax2.switch.whole <- switch.whole.pcoa$eig[2]/sum(switch.whole.pcoa$eig)

switch.whole.point_color <- rep("darkolivegreen3", nrow(switch.whole.map))
switch.whole.point_color[switch.whole.map$Year==2017] <- "darkolivegreen1"
switch.whole.point_color[switch.whole.map$Year==2016&switch.whole.map$source=="soil"] <- "burlywood"
switch.whole.point_color[switch.whole.map$Year==2017&switch.whole.map$source=="soil"] <- "brown"

switch.whole.point.size.2016 <- seq(from=1, to=3, length.out=length(unique(switch.whole.map$sampling_Rdate[switch.whole.map$Year==2016])))
switch.whole.point.size.2017 <- seq(from=1, to=3, length.out=length(unique(switch.whole.map$sampling_Rdate[switch.whole.map$Year==2017])))

switch.dates.2016 <- unique(switch.whole.map$sampling_Rdate[switch.whole.map$Year==2016])[order(unique(switch.whole.map$sampling_Rdate[switch.whole.map$Year==2016]))]
switch.dates.2017 <- unique(switch.whole.map$sampling_Rdate[switch.whole.map$Year==2017])[order(unique(switch.whole.map$sampling_Rdate[switch.whole.map$Year==2017]))]

switch.whole.point.sizes <- rep(1, nrow(switch.whole.map))

for(i in 1:length(switch.dates.2016)){
  switch.whole.point.sizes[switch.whole.map$sampling_Rdate==switch.dates.2016[i]] <- switch.whole.point.size.2016[i]  
}
for(i in 1:length(switch.dates.2017)){
  switch.whole.point.sizes[switch.whole.map$sampling_Rdate==switch.dates.2017[i]] <- switch.whole.point.size.2017[i]  
}


setEPS()
postscript("Switch_Whole_PCoa.eps", width=5, height=5,pointsize=10, paper="special",)
par(mar = c(5.1, 5.1, 2.1, 2.1))
plot(switch.whole.pcoa$points[,1],switch.whole.pcoa$points[,2], pch=16, col=switch.whole.point_color, cex=switch.whole.point.sizes, xlab=paste("PCoA1: ",100*round(ax1.switch.whole,3),"% var. explained",sep=""), ylab=paste("PCoA2: ",100*round(ax2.switch.whole,3),"% var. explained",sep=""))
dev.off()


switch.phyllo.otu <- switch.whole.otu[,switch.whole.map$source=="phyllosphere"]

switch.phyllo.ptaxa.otu <- switch.phyllo.otu[rowSums(switch.phyllo.otu)>0,]

switch.soil.otu <- switch.whole.otu[,switch.whole.map$source=="soil"]

switch.soil.ptaxa.otu <- switch.soil.otu[rowSums(switch.phyllo.otu)>0,]

switch.phyllo.avg.abund <- rowSums(switch.phyllo.ptaxa.otu)/sum(colSums(switch.phyllo.ptaxa.otu))
switch.soil.avg.abund <- rowSums(switch.soil.ptaxa.otu)/sum(colSums(switch.soil.ptaxa.otu))



Avg.abund.data <- data.frame(Average_Relative_Abundance_Phyllosphere=switch.phyllo.avg.abund, Average_Relative_Abundance_Soil=switch.soil.avg.abund)
AvgAbund_colors <- rep("black", nrow(Avg.abund.data))
AvgAbund_colors[row.names(Avg.abund.data)%in%Switch_core] <- "green"
AvgAbund_colors[row.names(Avg.abund.data)%in%Switch_core.2016] <- "darkolivegreen3"
AvgAbund_colors[row.names(Avg.abund.data)%in%Switch_core.2017] <- "darkolivegreen1"

AvgAbund_size <- rep(1, nrow(Avg.abund.data) )
AvgAbund_size[row.names(Avg.abund.data)%in%Switch_core] <- 3
AvgAbund_size[row.names(Avg.abund.data)%in%Switch_core.2016] <- 3
AvgAbund_size[row.names(Avg.abund.data)%in%Switch_core.2017] <- 3

FigX <-ggplot(Avg.abund.data, aes(y=log10(Average_Relative_Abundance_Phyllosphere), x=log10(Average_Relative_Abundance_Soil))) + geom_point(color=AvgAbund_colors, size=AvgAbund_size)
ggsave(filename ="RelativeABundancePlot.eps",FigX, device="eps", width=5, height = 5, units = "in" )





for (i in 1:length(unique(map.plant$sampling_Rdate))){
  b <- plant_coordinates[map.plant$sampling_Rdate==Sample_Dates[i],]
  dm <- map.plant[map.plant$sampling_Rdate==Sample_Dates[i],]
  m <- b[dm$plant=="miscanthus",]
  s.2016 <- b[dm$plant=="switchgrass"&dm$year=="2016",]
  s.2017 <- b[dm$plant=="switchgrass"&dm$year=="2017",]
  misc_collapse.2016[i,1] <-mean(m[,1])
  misc_collapse.2016[i,2] <- mean(m[,2])
  misc_collapse.2016[i,3] <-sd(m[,1])
  misc_collapse.2016[i,4] <-sd(m[,1])
  switch_collapse.2016[i,1] <-mean(s.2016[,1])
  switch_collapse.2016[i,2] <- mean(s.2016[,2])
  switch_collapse.2016[i,3] <-sd(s.2016[,1])
  switch_collapse.2016[i,4] <-sd(s.2016[,1])
  switch_collapse.2017[i,1] <-mean(s.2017[,1])
  switch_collapse.2017[i,2] <- mean(s.2017[,2])
  switch_collapse.2017[i,3] <-sd(s.2017[,1])
  switch_collapse.2017[i,4] <-sd(s.2017[,1])
}

misc_collapse.2016$SampleDates <- Sample_Dates
switch_collapse.2016$SampleDates <- Sample_Dates
switch_collapse.2017$SampleDates <- Sample_Dates

#weather_map2 <- map_16S[,c("precipitation", "Air_temp_mean", "air_temp_max", "Air_Temp_Min", "Air_Pressure", "RH", "AH", "Wind_Speed_Mean", "Solar_Radiation", "PAR", "soil_temp_5_cm_bare_avg", "soil_temp_5_cm_sod_avg", "sampling_Rdate", "time_numeric")]

weather_map2 <- map_16S[,c("precipitation", "Air_temp_mean", "air_temp_max", "Air_Temp_Min", "Air_Pressure", "RH", "AH", "Wind_Speed_Mean", "Solar_Radiation", "PAR", "soil_temp_5_cm_bare_avg", "sampling_Rdate", "time_numeric")]

# Leaf chemistry envfit on Plant PCoA
Leaf_Chemistry.2016 <- map.plant.2016[,c("LDMC_mg_per_g", "nitrogen_percent", "carbon_percent", "carbon_per_nitrogen", "height_mean_cm")]
LC.env.2016 <- envfit(plant.pcoa.2016, Leaf_Chemistry.2016)

R40_LC.env.2016 <- envfit(plant.pcoa.2016, Leaf_Chemistry.2016[,"height_mean_cm", drop=FALSE])

Leaf_Chemistry <- map.plant[,c("LDMC_mg_per_g", "nitrogen_percent", "carbon_percent", "carbon_per_nitrogen")]
LC.env <- envfit(plant.pcoa, Leaf_Chemistry)

R40_LC.env <- envfit(plant.pcoa, Leaf_Chemistry[,"nitrogen_percent"])


# Soil chemistry envfit on plant pcoa
# Don't include lime index because some values are missing. 
Soil_Chemistry_phyllo.2016<- map.plant.2016[,c("pH", "P_ppm", "K_ppm", "Ca_ppm", "Mg_ppm", "organic_matter", "NO3N_ppm", "NH4_ppm", "soil_moisture_percent", "soil_temp_10cm")]
for (i in 1:ncol(Soil_Chemistry_phyllo.2016)){
  Soil_Chemistry_phyllo.2016[,i] <- as.numeric(Soil_Chemistry_phyllo.2016[,i])
}

Soil_Chemistry_phyllo <- map.plant[,c("pH", "P_ppm", "K_ppm", "Ca_ppm", "Mg_ppm", "organic_matter", "NO3N_ppm", "NH4_ppm")]
for (i in 1:ncol(Soil_Chemistry_phyllo)){
  Soil_Chemistry_phyllo[,i] <- as.numeric(Soil_Chemistry_phyllo[,i])
}

SC_v_Phyllo.env.2016 <- envfit(plant.pcoa.2016, Soil_Chemistry_phyllo.2016)
SC_v_Phyllo.env <- envfit(plant.pcoa, Soil_Chemistry_phyllo)



plot_colors.2016 <- c(rep("darkgreen", 9), rep("darkolivegreen3", 8))
plot_colors <- c(rep("darkgreen", 9), rep("darkolivegreen3", 8), rep("darkolivegreen1",7))
collapsed_data.pcoa.2016 <- rbind(misc_collapse.2016, switch_collapse.2016[complete.cases(switch_collapse.2016),])
collapsed_data.pcoa <-  rbind(misc_collapse.2016[complete.cases(misc_collapse.2016),], switch_collapse.2016[complete.cases(switch_collapse.2016),], switch_collapse.2017[complete.cases(switch_collapse.2017),])


# Weather EnvFit on plant PCoA
weather_map2 <- unique(weather_map2)
nrow(weather_map2)
weather_map2 <- subset(weather_map2, sampling_Rdate%in%map.plant$sampling_Rdate)

weather_map_bothyears <- weather_map2[order(weather_map2$sampling_Rdate),]

weather_map3 <- rbind(weather_map2, weather_map2[1:5,])

weather_map_forEnvFit_bothyears <- rbind(weather_map_bothyears[1:9,], weather_map_bothyears[c(1:8, 10:16),])

weather_map_forEnvFit_bothyears$sampling_Rdate==collapsed_data.pcoa$SampleDates

# Remove Sampling Date Column
weather_map4 <- weather_map3[,-13]

collapsed.ef.2016 <- envfit(collapsed_data.pcoa.2016[,1:2], weather_map4)
collapsed.weather.ef <- envfit(collapsed_data.pcoa[,1:2], weather_map_forEnvFit_bothyears)
R40_collapsed.ef.2016 <- envfit(collapsed_data.pcoa.2016[,1:2], weather_map4[,c("RH", "time_numeric", "Wind_Speed_Mean")])

R40_collapsed.weather.ef <- envfit(collapsed_data.pcoa[,1:2], weather_map_forEnvFit_bothyears[,c("RH", "Wind_Speed_Mean")])

plant_shapes.2016 <- c(rep(15, 9), rep(16,8))
plant_shapes <- c(rep(15,9), rep(16,8), rep(1,7))


Phyllosphere_Sampling_Dates.2016 <- unique(collapsed_data.pcoa.2016$SampleDates)[order(unique(collapsed_data.pcoa.2016$SampleDates))]
phyllosphere_point_size.2016 <- rep(1, nrow(collapsed_data.pcoa.2016))
phyllosphere_point_size <- rep(1, nrow(collapsed_data.pcoa))


for(i in 1:nrow(collapsed_data.pcoa.2016)){
  phyllosphere_point_size.2016[collapsed_data.pcoa.2016$SampleDates==Unique_Dates[i]] <- point_sizes_whole[i]
}

for(i in 1:nrow(collapsed_data.pcoa)){
  phyllosphere_point_size[collapsed_data.pcoa$SampleDates==Unique_Dates[i]] <- point_sizes_whole[i]
}

for(i in 1:nrow(collapsed_data.pcoa)){
  phyllosphere_point_size[collapsed_data.pcoa$SampleDates==Unique_Dates.2017[i]] <- point_sizes_whole[2:8][i]
}


# Plot centroid and standard deviation (collapsing all biological reps and nitrogen into a single point per date)
# Make plot without labels on EnvFit vectors so they don't overlapp
setEPS()
postscript("Figure1A_Plant_NoLabels.eps", width=5, height=5,pointsize=10, paper="special",)
par(mar = c(5.1, 5.1, 2.1, 2.1))
plot(collapsed_data.pcoa.2016$ax1.average, collapsed_data.pcoa.2016$ax2.average, pch=plant_shapes.2016, cex=phyllosphere_point_size.2016, cex.axis=1.6, cex.lab=1.6, ylim = c(-0.3,0.4), xlim=c(-0.5,0.5), col= plot_colors.2016, xlab=paste("PCoA1: ",100*round(ax1.plant.2016,3),"% var. explained",sep=""), ylab=paste("PCoA2: ",100*round(ax2.plant.2016,3),"% var. explained",sep="")) +
arrows(collapsed_data.pcoa.2016$ax1.average, collapsed_data.pcoa.2016$ax2.average-collapsed_data.pcoa.2016$ax1.sd, collapsed_data.pcoa.2016$ax1.average, collapsed_data.pcoa.2016$ax2.average+collapsed_data.pcoa.2016$ax2.sd, angle=90,length= 0.05, code=3,lwd=0.75, col= plot_colors) +
arrows(collapsed_data.pcoa.2016$ax1.average- collapsed_data.pcoa.2016$ax1.sd, collapsed_data.pcoa.2016$ax2.average, collapsed_data.pcoa.2016$ax1.average + collapsed_data.pcoa.2016$ax1.sd, collapsed_data.pcoa.2016$ax2.average, angle=90,length= 0.05, code=3,lwd=0.75, col= plot_colors) +
plot(R40_collapsed.ef.2016, p=0.05, col="black", labels=NA)+
plot(R40_LC.env.2016, p=0.05, col="black", labels= NA)
#plot(SC_v_Phyllo.env, p=0.05, col="grey", labels=NA)+
#text(collapsed_data.pcoa$ax1.average, collapsed_data.pcoa$ax2.average,c(8,10,5,7,9,6,8,5,11), font=2, cex=1, pos=4)
dev.off()


### Combined 2016 and 2017 Phyllosphere PCoAs
setEPS()
postscript("Phyllosphere_PCoA_BothYears.eps", width=5, height=5,pointsize=10, paper="special",)
par(mar = c(5.1, 5.1, 2.1, 2.1))
plot(collapsed_data.pcoa$ax1.average, collapsed_data.pcoa$ax2.average, pch=plant_shapes, cex=phyllosphere_point_size, cex.axis=1.6, cex.lab=1.6, ylim = c(-0.4,0.4), xlim=c(-0.5,0.5), col= plot_colors, xlab=paste("PCoA1: ",100*round(ax1.plant,3),"% var. explained",sep=""), ylab=paste("PCoA2: ",100*round(ax2.plant,3),"% var. explained",sep="")) +
  arrows(collapsed_data.pcoa$ax1.average, collapsed_data.pcoa$ax2.average-collapsed_data.pcoa$ax1.sd, collapsed_data.pcoa$ax1.average, collapsed_data.pcoa$ax2.average+collapsed_data.pcoa$ax2.sd, angle=90,length= 0.05, code=3,lwd=0.75, col= plot_colors) +
  arrows(collapsed_data.pcoa$ax1.average- collapsed_data.pcoa$ax1.sd, collapsed_data.pcoa$ax2.average, collapsed_data.pcoa$ax1.average + collapsed_data.pcoa$ax1.sd, collapsed_data.pcoa$ax2.average, angle=90,length= 0.05, code=3,lwd=0.75, col= plot_colors)+
  plot(R40_collapsed.weather.ef, p=0.05, col="black")+
  plot(R40_LC.env, p=0.05, col="black", labels="nitrogen_percent")
dev.off()


setEPS()
postscript("Figure1A_Plant_Labels.eps", width=5, height=5,pointsize=10, paper="special",)
par(mar = c(5.1, 5.1, 2.1, 2.1))
plot(collapsed_data.pcoa.2016$ax1.average, collapsed_data.pcoa.2016$ax2.average, pch=plant_shapes.2016, cex=phyllosphere_point_size.2016, cex.axis=1.6, cex.lab=1.6, ylim = c(-0.3,0.4), xlim=c(-0.5,0.5), col= plot_colors.2016, xlab=paste("PCoA1: ",100*round(ax1.plant.2016,3),"% var. explained",sep=""), ylab=paste("PCoA2: ",100*round(ax2.plant.2016,3),"% var. explained",sep="")) +
  arrows(collapsed_data.pcoa.2016$ax1.average, collapsed_data.pcoa.2016$ax2.average-collapsed_data.pcoa.2016$ax1.sd, collapsed_data.pcoa.2016$ax1.average, collapsed_data.pcoa.2016$ax2.average+collapsed_data.pcoa.2016$ax2.sd, angle=90,length= 0.05, code=3,lwd=0.75, col= plot_colors.2016) +
  arrows(collapsed_data.pcoa.2016$ax1.average- collapsed_data.pcoa.2016$ax1.sd, collapsed_data.pcoa.2016$ax2.average, collapsed_data.pcoa.2016$ax1.average + collapsed_data.pcoa.2016$ax1.sd, collapsed_data.pcoa.2016$ax2.average, angle=90,length= 0.05, code=3,lwd=0.75, col= plot_colors.2016) +
  plot(R40_collapsed.ef.2016, p=0.05, col="black")+
  plot(R40_LC.env.2016, p=0.05, col="black")
#plot(SC_v_Phyllo.env, p=0.05, col="grey", labels=NA)+
#text(collapsed_data.pcoa$ax1.average, collapsed_data.pcoa$ax2.average,c(8,10,5,7,9,6,8,5,11), font=2, cex=1, pos=4)
dev.off()

###################################
### Table 1 and S1 Plant Envfit ###
###################################
Weather_plant_EnvFit.2016 <- as.data.frame(collapsed.ef.2016$vectors$arrows)
colnames(Weather_plant_EnvFit.2016) <- c("Dim1", "Dim2")
Weather_plant_EnvFit.2016$Rsquared <- collapsed.ef.2016$vectors$r
Weather_plant_EnvFit.2016$pvalue <- collapsed.ef.2016$vectors$pvals
Weather_plant_EnvFit.2016$DataType <- rep("Weather", nrow(Weather_plant_EnvFit.2016))


LC_plant_EnvFit.2016 <- as.data.frame(LC.env.2016$vectors$arrows)
LC_plant_EnvFit.2016$Rsquared <- LC.env.2016$vectors$r
LC_plant_EnvFit.2016$pvalue <- LC.env.2016$vectors$pvals
LC_plant_EnvFit.2016$DataType <- rep("Leaf Chemistry", nrow(LC_plant_EnvFit.2016))

SC_plant_EnvFit.2016 <- as.data.frame(SC_v_Phyllo.env.2016$vectors$arrows)
SC_plant_EnvFit.2016$Rsquared <-SC_v_Phyllo.env.2016$vectors$r
SC_plant_EnvFit.2016$pvalue <- SC_v_Phyllo.env.2016$vectors$pvals
SC_plant_EnvFit.2016$DataType <- rep("Soil Chemistry", nrow(SC_plant_EnvFit.2016))

All_plant_EnvFit <- rbind(Weather_plant_EnvFit.2016, LC_plant_EnvFit.2016, SC_plant_EnvFit.2016)

Sig_plant_EnvFit <- All_plant_EnvFit[All_plant_EnvFit$pvalue<0.05,]
NonSig_plant_EnvFit <- All_plant_EnvFit[All_plant_EnvFit$pvalue>0.05,]

write.table(file = "Table1_PhyllosphereEnvFit.txt", x=Sig_plant_EnvFit, sep="\t", quote=FALSE, col.names = TRUE, row.names = TRUE)
write.table(file="TableS1_PhyllosphereEnvFit.txt", x=NonSig_plant_EnvFit, sep="\t", quote=FALSE, col.names = TRUE, row.names=TRUE)

###########################
### Figure 1C Soil PCoA ###
###########################

soil.pcoa <- cmdscale(soil.dist, eig=TRUE)
ax1.soil <- soil.pcoa$eig[1]/sum(soil.pcoa$eig)
ax2.soil <- soil.pcoa$eig[2]/sum(soil.pcoa$eig)

soil.pcoa.2016 <- cmdscale(soil.dist.2016, eig=TRUE)
ax1.soil.2016 <- soil.pcoa.2016$eig[1]/sum(soil.pcoa.2016$eig)
ax2.soil.2016 <- soil.pcoa.2016$eig[2]/sum(soil.pcoa.2016$eig)

Soil_Chemistry_soil <- map.soil[,c("pH", "P_ppm", "K_ppm", "Ca_ppm", "Mg_ppm", "organic_matter", "NO3N_ppm", "NH4_ppm", "soil_moisture_percent", "soil_temp_10cm")]
for (i in 1:ncol(Soil_Chemistry_soil)){
  Soil_Chemistry_soil[,i] <- as.numeric(Soil_Chemistry_soil[,i])
}

soil_coordinates <- soil.pcoa$points



Soil_Chemistry_soil.2016 <- map.soil.2016[,c("pH", "P_ppm", "K_ppm", "Ca_ppm", "Mg_ppm", "organic_matter", "NO3N_ppm", "NH4_ppm", "soil_moisture_percent", "soil_temp_10cm")]
for (i in 1:ncol(Soil_Chemistry_soil.2016)){
  Soil_Chemistry_soil.2016[,i] <- as.numeric(Soil_Chemistry_soil.2016[,i])
}

soil_coordinates.2016 <- soil.pcoa.2016$points

### Make collapsed PCoA for Weather EnvFit

dates <- unique(map.soil$sampling_Rdate)[order(unique(map.soil$sampling_Rdate))]

dates.2016 <- unique(map.soil.2016$sampling_Rdate)[order(unique(map.soil.2016$sampling_Rdate))]

misc.soil_NF.2016 <- data.frame(ax1.average=rep(NA,10), ax1.sd = rep(NA,10), ax2.average=rep(NA,10), ax2.sd=rep(NA,10))
misc.soil_MAIN.2016 <- data.frame(ax1.average=rep(NA,10), ax1.sd = rep(NA,10), ax2.average=rep(NA,10), ax2.sd=rep(NA,10))

switch.soil_NF.2016 <- data.frame(ax1.average=rep(NA,9), ax1.sd = rep(NA,9), ax2.average=rep(NA,9), ax2.sd=rep(NA,9))
switch.soil_MAIN.2016 <- data.frame(ax1.average=rep(NA,9), ax1.sd = rep(NA,9), ax2.average=rep(NA,9), ax2.sd=rep(NA,9))

for (i in 1:length(dates)){
 w <- soil_coordinates.2016[map.soil.2016$sampling_Rdate==dates[i]&map.soil.2016$plant=="miscanthus"&map.soil.2016$treatment=="standard fertilization",]
 x <- soil_coordinates.2016[map.soil.2016$sampling_Rdate==dates[i]&map.soil.2016$plant=="miscanthus"&map.soil.2016$treatment=="nitrogen free",]
y <- soil_coordinates.2016[map.soil.2016$sampling_Rdate==dates[i]&map.soil.2016$plant=="switchgrass"&map.soil.2016$treatment=="standard fertilization",]
z <- soil_coordinates.2016[map.soil.2016$sampling_Rdate==dates[i]&map.soil.2016$plant=="switchgrass"&map.soil.2016$treatment=="nitrogen free",]
misc.soil_MAIN.2016[i,] <- c(mean(w[,1]),sd(w[,1]), mean(w[,2]), sd(w[,2]))
misc.soil_NF.2016[i,] <- c(mean(x[,1]),sd(x[,1]), mean(x[,2]), sd(x[,2]))
switch.soil_MAIN.2016[i,] <- c(mean(y[,1]),sd(y[,1]), mean(y[,2]), sd(y[,2]))
switch.soil_NF.2016[i,] <- c(mean(z[,1]),sd(z[,1]), mean(z[,2]), sd(z[,2]))
}

Collapsed_Soil.2016 <- rbind(misc.soil_MAIN.2016, misc.soil_NF.2016, switch.soil_MAIN.2016[1:9,], switch.soil_NF.2016[1:9,])
Collapsed_Soil.2016$plant <- c(rep("miscanthus", 20), rep("switchgrass",18))
Collapsed_Soil.2016$Date <- c(dates,dates,dates[1:9], dates[1:9])


weather_soil.2016 <- map.soil.2016[,c("precipitation", "Air_temp_mean", "air_temp_max", "Air_Temp_Min", "Air_Pressure", "RH", "AH", "Wind_Speed_Mean", "Solar_Radiation", "PAR", "soil_temp_5_cm_bare_avg", "soil_temp_5_cm_sod_avg", "sampling_Rdate", "time_numeric")]
weather_soil2.2016 <- unique(weather_soil.2016)
weather_soil3.2016 <- weather_soil2.2016[order(weather_soil2.2016$sampling_Rdate),]
weather_soil4.2016 <- rbind(weather_soil3.2016, weather_soil3.2016, weather_soil3.2016[1:9,], weather_soil3.2016[1:9,])

weather_soil4.2016$sampling_Rdate==Collapsed_Soil.2016$Date

# Soil Chemistry EnvFit
R40_SoilChemistry_soil.2016 <- Soil_Chemistry_soil.2016[,"pH", drop=FALSE]
SC_Soil.env.2016 <- envfit(soil.pcoa.2016, Soil_Chemistry_soil.2016)
R40_SC_Soil.env.2016 <- envfit(soil.pcoa.2016, R40_SoilChemistry_soil.2016)

#Weather Soil EnvFit
R40_weather_soil.2016 <- weather_soil4.2016[,c("Wind_Speed_Mean","precipitation","time_numeric")]
weather_soil.env.2016 <- envfit(Collapsed_Soil.2016[,c(1,3)], weather_soil4.2016[,c(1:12,14)])
R40_weather_soil.env.2016 <- envfit(Collapsed_Soil.2016[,c(1,3)], R40_weather_soil.2016)


library(lubridate)

### Set point colors based on plant (Dark == Miscanthus, light == Switchgrass)
soil_plot_colors.2016 <- rep("burlywood4", nrow(map.soil.2016))
soil_plot_colors.2016[map.soil.2016$plant=="switchgrass"] <- "burlywood"

### Set point shapes based on plant (Squares == Miscanthus, Circles ==Switchgrass, open== Nitrogen Free)
soil_plot_shapes.2016 <- rep(15, nrow(map.soil.2016))
soil_plot_shapes.2016[map.soil.2016$plant=="switchgrass"] <- 16
soil_plot_shapes.2016[map.soil.2016$plant=="switchgrass"& map.soil.2016$treatment=="nitrogen free"] <- 1
soil_plot_shapes.2016[map.soil.2016$plant=="miscanthus"&map.soil.2016$treatment=="nitrogen free"] <- 0

### Set point sizes based on date (Smaller symbols earlier in season)
soil_plot_size.2016 <- rep(1, nrow(map.soil.2016))
for (i in 1:length(Unique_Dates)){
  soil_plot_size.2016[map.soil.2016$sampling_date==Unique_Dates[i]] <- point_sizes_whole[i]
}


setEPS()
postscript("Figure1C_Soil_Labels.eps", width=5, height=5,pointsize=10, paper="special",)
par(mar = c(5.1, 5.1, 2.1, 2.1))
Soil_PcoA_Plot.2016 <- plot(soil.pcoa.2016$points[,1], soil.pcoa.2016$points[,2], pch=soil_plot_shapes.2016, col=soil_plot_colors.2016, cex=soil_plot_size.2016, cex.axis=1.6, cex.lab=1.6, xlab=paste("PCoA1: ",100*round(ax1.soil.2016,3),"% var. explained",sep=""), ylab=paste("PCoA2: ",100*round(ax2.soil.2016,3),"% var. explained",sep=""))+
  #text(soil.pcoa$points[,1], soil.pcoa$points[,2], map.soil$month, pos=4,  font=2)+
plot(R40_SC_Soil.env.2016, p=0.05, col="black")+
plot(R40_weather_soil.env.2016, p=0.05, col="black")
dev.off()

setEPS()
postscript("Figure1C_Soil_NoLabels.eps", width=5, height=5,pointsize=10, paper="special",)
par(mar = c(5.1, 5.1, 2.1, 2.1))
Soil_PcoA_Plot.2016 <- plot(soil.pcoa.2016$points[,1], soil.pcoa.2016$points[,2], pch=soil_plot_shapes.2016, col=soil_plot_colors.2016, cex=soil_plot_size.2016, cex.axis=1.6, cex.lab=1.6, xlab=paste("PCoA1: ",100*round(ax1.soil.2016,3),"% var. explained",sep=""), ylab=paste("PCoA2: ",100*round(ax2.soil.2016,3),"% var. explained",sep=""))+
  #text(soil.pcoa$points[,1], soil.pcoa$points[,2], map.soil$month, pos=4,  font=2)+
  plot(R40_SC_Soil.env.2016, p=0.05, col="black", labels = FALSE)+
  plot(R40_weather_soil.env.2016, p=0.05, col="black", labels = FALSE)
dev.off()

############################
### Table S2 Soil Envfit ###
############################

SC_soil_EnvFit.2016 <- as.data.frame(SC_Soil.env.2016$vectors$arrows)
SC_soil_EnvFit.2016$Rsquared <- SC_Soil.env.2016$vectors$r
SC_soil_EnvFit.2016$pvalue <- SC_Soil.env.2016$vectors$pvals
SC_soil_EnvFit.2016$DataType <- rep("Soil Chemistry", nrow(SC_soil_EnvFit.2016))

weather_soil_EnvFit.2016 <- as.data.frame(weather_soil.env.2016$vectors$arrows)
weather_soil_EnvFit.2016$Rsquared <- weather_soil.env.2016$vectors$r
weather_soil_EnvFit.2016$pvalue <- weather_soil.env.2016$vectors$pvals
weather_soil_EnvFit.2016$DataType <- rep("Weather", nrow(weather_soil_EnvFit.2016))

colnames(weather_soil_EnvFit.2016) <- colnames(SC_soil_EnvFit.2016)

AllVar_soil_EnvFit.2016 <- rbind(SC_soil_EnvFit.2016, weather_soil_EnvFit.2016)


write.table(file = "TableS2_SoilEnvFit.txt", x=AllVar_soil_EnvFit.2016, sep="\t", quote=FALSE, col.names = TRUE, row.names = TRUE)


######################################
#### 2017 Whole Distance Matrices ####
######################################
whole.dist.2017 <- vegdist(t(otu_rare.2017), method = "bray")
whole.pcoa.2017 <- cmdscale(whole.dist.2017, eig = TRUE)
plot(whole.pcoa.2017$points[,1], whole.pcoa.2017$points[,2])


#############################################
#### Hypothesis Testing of Beta Diversity####
#############################################
adonis(whole.dist.2016~map.whole.2016$source)
  #Soil is distinct from plant

adonis(whole.dist.2017~map.2017$source)

adonis(plant.dist.2016~map.plant.2016$plant)
  # PLant type has significant effect on phyllosphere community

adonis(soil.dist.2016~map.soil.2016$plant)
  # Plant type has significant effect on soil community

adonis(soil.dist.2016~map.soil.2016$time_numeric)
  # Time has significant influence on soil community

adonis(plant.dist.2016~map.plant.2016$time_numeric)
  # Time has significant influence on phyllosphere community  

adonis(soil.dist.2016~map.soil.2016$treatment)
  # Fertilization has significant effect on Soil community

adonis(plant.dist.2016~map.plant.2016$treatment)
  # Fertilization has no significant effect on plant community
adonis(plant.dist.2016~map.plant.2016$time_numeric*map.plant.2016$treatment)


adonis(soil.dist.2017~map.soil.2017$time_numeric)
  # Soil exhibits seasonal patterns
adonis(soil.dist.2017~map.soil.2017$treatment)
  # Fertilization has significant effect on soil community
adonis(switch.dist.2017~switch.map.2017$time_numeric)
  # Phyllosphere exhibits seasonal patterns
adonis(switch.dist.2017~switch.map.2017$treatment)
  # Fertilization has no significant effect on phyllosphere community

adonis(switch.dist~switch.map$year)


permutest(betadisper(whole.dist.2016, map.2016$source))

##########################################################################
#### Occupancy and abundance of Share taxa in Soils and Phyllosphere ####
##########################################################################

plant_rare_taxa.2017 <- 2*(rowSums(plant_rare_otu.2017)>0)

plant_taxa.2016 <- (1*rowSums(plant_rare_otu.2016)>0)
soil_taxa.2016 <- 2*(rowSums(soil_rare_otu.2016)>0)

switch.rare.taxa.2016 <- 1*(rowSums(switch.rare.otu.2016)>0)
switch.rare.taxa.2017 <- 2*(rowSums(switch.rare.otu.2017)>0)
  
table(switch.rare.taxa.2016 + switch.rare.taxa.2017)  
  
# Abundance of shared taxa in soil samples
soil_both_taxa.2016 <- soil_rare_otu.2016[(plant_taxa.2016 + soil_taxa.2016)==3,]

#Abundance of shared taxa in phyllosphere samples
plant_both_taxa.2016 <- plant_rare_otu.2016[(plant_taxa.2016 + soil_taxa.2016)==3,]


colSums(soil_both_taxa.2016)/colSums(soil_rare_otu.2016)
range(colSums(soil_both_taxa.2016)/colSums(soil_rare_otu.2016))

# Relative abundance of all shared taxa in soil samples
RA_st_soil.2016 <- colSums(soil_both_taxa.2016)/colSums(soil_rare_otu.2016)

# Average relative abudnace of all shared taxa in soil samples
ARA_st_soil.2016 <- RA_st_soil.2016/nrow(soil_both_taxa.2016)

# Shared Taxa's average relative abundance in soil samples
ST_ARA_S.2016<-rowSums(soil_both_taxa.2016)/sum(colSums(soil_rare_otu.2016))

# Shared Taxa's average relative abundance in phyllosphere samples
ST_ARA_P.2016<-rowSums(plant_both_taxa.2016)/sum(colSums(plant_rare_otu.2016))

ST_ARA_plotdata <- data.frame(OTU=row.names(plant_both_taxa.2016), SoilAbundance=ST_ARA_S.2016, PhyllosphereAbundance=ST_ARA_P.2016)

plant_both_taxa_PA.2016 <- 1*(plant_both_taxa.2016>0)

soil_both_taxa_PA.2016 <- 1*(soil_both_taxa.2016>0)

ST_Occ_Plant.2016 <- rowSums(plant_both_taxa_PA.2016)/ncol(plant_both_taxa_PA.2016)

ST_Occ_Soil.2016 <- rowSums(soil_both_taxa_PA.2016)/ ncol(soil_both_taxa_PA.2016) 

ST_ARA_plotdata.2016 <- data.frame(OTU=row.names(plant_both_taxa.2016), SoilAbundance=ST_ARA_S.2016, PhyllosphereAbundance=ST_ARA_P.2016)

ST_Occ_plotdata.2016 <- data.frame(OTU=row.names(plant_both_taxa.2016), SoilOccupancy=ST_Occ_Soil.2016, PhyllosphereOccupancy=ST_Occ_Plant.2016)


core_taxa.df <- read.table("~/Downloads/core_list.csv", header = TRUE, row.names = 1, sep=",", stringsAsFactors = FALSE)
core_taxa_2016.df <- as.character(core_taxa.df$x[!core_taxa.df$x%in%c("OTU84", "OTU86", "OTU4223")])

Switch_core <- c("OTU5", "OTU47", "OTU14", "OTU2771", "OTU7", "OTU10","OTU18","OTU41","OTU4", "OTU23", "OTU2","OTU430", "OTU6", "OTU22", "OTU21", "OTU192", "OTU90", "OTU537")
Switch_core.2016 <- c("OTU995", "OTU15", "OTU32", "OTU1334", "OTU6096", "OTU842")
Switch_core.2017 <- c("OTU519", "OTU84", "OTU86", "OTU4223")


core_taxa <- c("53","117","4","16","27","32","8","9","20","19","2","25","3","21","37","30","6","44","63","920","10","101","641","80")

switchgrass_core_taxa <- c("104","12","1751","157","108","357","61")

miscanthus_core_taxa <- c("83","149","22","40")

write.table(file="Core_Taxa.txt",x=c(core_taxa,switchgrass_core_taxa, miscanthus_core_taxa), quote=FALSE, row.names = FALSE)


core_taxa.2016 <- c("OTU47","OTU10", "OTU18", "OTU21", "OTU23", "OTU2", "OTU430", "OTU6", "OTU41", "OTU7", "OTU14", "OTU2771", "OTU5", "OTU4", "OTU22", "OTU995", "OTU32")

core_taxa_misc.2016 <- c("OTU15", "OTU48", "OTU50", "OTU3994", "OTU995", "OTU519", "OTU1674", "OTU32")
core_taxa_switch.2016 <- c("OTU192", "OTU6096", "OTU537", "OTU90", "OTU1334", "OTU842")

FigS4_Colors <- rep("Black", length(ST_ARA_P.2016))
FigS4_Colors[row.names(plant_both_taxa_PA.2016)%in%core_taxa.2016] <- "Green"
FigS4_Colors[row.names(plant_both_taxa_PA.2016)%in%core_taxa_switch.2016] <- "darkolivegreen3"
FigS4_Colors[row.names(plant_both_taxa_PA.2016)%in%core_taxa_misc.2016] <- "darkgreen"

FigS4_Size <- rep(1, length(ST_ARA_P.2016))
FigS4_Size[row.names(plant_both_taxa_PA.2016)%in%c(core_taxa.2016, core_taxa_switch.2016, core_taxa_misc.2016)] <- 2



Fig3_A.2016 <- ggplot(data=ST_ARA_plotdata.2016, aes(x=log10(SoilAbundance), y=log10(PhyllosphereAbundance)))+
  geom_point(color=FigS4_Colors, size= FigS4_Size)+
  labs(x=expression('log'[10]*'(Relative Abundance Soil)'), y=expression('log'[10]*'(Relative Abundance Phyllosphere)'))
  

Fig3_B.2016 <- ggplot(data=ST_Occ_plotdata.2016, aes(x=SoilOccupancy, y=PhyllosphereOccupancy))+
  geom_point(color=FigS4_Colors, size=FigS4_Size)+
  labs(x="Soil Occupancy", y="Phyllosphere Occupancy")




ggsave("Figure3A.eps", plot=Fig3_A.2016, height = 4, width = 4, device="eps")

ggsave("Figre3B.eps", plot=Fig3_B.2016, height = 4, width=4, device = "eps")

# Average shared taxa abundance in phyllosphere samples




colSums(plant_both_taxa)/colSums(plant_rare_otu)

soil_full_otu.2016 <- otu[,map_16S$source=="soil"&map_16S$year==2016]
plant_full_otu.2016 <- otu[,map_16S$source=="phyllosphere"&map_16S$year==2016]

misc_soil_full_otu.2016 <- soil_full_otu.2016[,map.soil.2016$plant=="miscanthus"]
switch_soil_full_otu.2016 <- soil_full_otu.2016[,map.soil.2016$plant=="switchgrass"]

misc_soil_taxa.2016 <- 1*(rowSums(misc_soil_full_otu.2016)>0)
switch_soil_taxa.2016 <- 2*(rowSums(switch_soil_full_otu.2016)>0)

table(misc_soil_taxa.2016 + switch_soil_taxa.2016)
21795/(21795+2638+2304)  

misc_soil_rare_otu.2016 <- otu_soil_rare[,map.soil$plant=="miscanthus"&map.soil$year==2016]
switch_soil_rare_otu.2016 <- otu_soil_rare[,map.soil$plant=="switchgrass"&map.soil$year==2016]

misc_soil_rare_taxa.2016 <- 1*(rowSums(misc_soil_rare_otu.2016)>0)
switch_soil_rare_taxa.2016 <- 2*(rowSums(switch_soil_rare_otu.2016)>0)
table(misc_soil_rare_taxa.2016 + switch_soil_rare_taxa.2016)


soil_full_taxa.2016 <- (2*(rowSums(soil_full_otu.2016)>0))
plant_full_taxa <- (1*(rowSums(plant_full_otu)>0))

plant_rare_taxa.2016 <- (1*(rowSums(otu_rare[,map_16S$source=="phyllosphere"&map_16S$year==2016])>0))
table(plant_rare_taxa.2016+soil_full_taxa.2016)
916/(135+916)


# Relative abundance of all shared taxa in soil samples
RA_ST_full_soil <- colSums(soil_full_both_taxa)/colSums(soil_full_otu)

mean(RA_ST_full_soil)
range(RA_ST_full_soil)

nrow(soil_full_otu[(plant_full_taxa+ soil_full_taxa)==3,])/sum(rowSums(plant_full_otu)>0)

############################################################################
#### Making Table of Core Taxa Blast Results against Isolate Collection ####
############################################################################
Core_taxa_BlastResults <- read.table("~/GLBRC_Isoaltes/Core_Taxa_Isolate_Blast_Results_6.txt", sep="\t", header=FALSE, row.names=NULL, stringsAsFactors = FALSE)

colnames(Core_taxa_BlastResults) <- c("OTU", "isolate_name", "Percent_Identity", "Alignment_Length", "Mismatches","Gap_opens", "Query_start", "Query_end", "Subject_start", "Subject_end", "E_value", "Bit_score")

All_Core_taxa <- unique(Core_taxa_BlastResults$OTU)

Best_Match <- NULL
Best_Matches <- NULL
for (i in 1:length(All_Core_taxa)){
  x <- Core_taxa_BlastResults[Core_taxa_BlastResults$OTU==All_Core_taxa[i],]
  z <- x[x$E_value==min(x$E_value),]
  Best_Match <- rbind(Best_Match, z[1,])
  Best_Matches <- rbind(Best_Matches, z)
}


Match_97 <- Best_Matches[Best_Matches$Percent_Identity>=97,]

Isolate_Info <- read.table("GLBRC_db_isolatestable.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE, row.names = NULL)


Joined_Blast_Hits <- inner_join(Best_Match, Isolate_Info, by="isolate_name")

Isolate_Taxonomy <- read.table("fixrank_Trimmed_GLBRC_Isolate.fasta_classified.txt", sep=";", header=FALSE, stringsAsFactors = FALSE, row.names = NULL)

colnames(Isolate_Taxonomy) <- c("isolate_name", "strand", "Domain", "Domain_confidence", "Phylum", "Phylum_confidence", "Class", "Class_confidence", "Order", "Order_confidence", "Family", "Family_confidence", "Genus", "Genus_confidence")

Joined_Blast_RDPTaxonomy_Isolates <- inner_join(Joined_Blast_Hits, Isolate_Taxonomy)

simple_core <- data.frame(OTU=as.numeric(c(core_taxa, switchgrass_core_taxa,miscanthus_core_taxa)), BlastHit=rep("Yes", length(c(core_taxa, switchgrass_core_taxa,miscanthus_core_taxa))))

Joined_Blast_RDPTaxonomy_Isolates <- left_join(simple_core, Joined_Blast_RDPTaxonomy_Isolates, by="OTU")
Joined_Blast_RDPTaxonomy_Isolates <- Joined_Blast_RDPTaxonomy_Isolates[,-2]

Core_Taxa_Isolate_Summary <- Joined_Blast_RDPTaxonomy_Isolates[,c(1,2,4,3,11,39:42)]

Core_Taxa_Isolate_Summary$Taxonomy <- paste("g:",Core_Taxa_Isolate_Summary$Genus)
Core_Taxa_Isolate_Summary[1,10] <- paste("f:", Core_Taxa_Isolate_Summary$Family[1])
Core_Taxa_Isolate_Summary$Isolate_Taxonomy <- gsub("g: NA", "No Match", x=Core_Taxa_Isolate_Summary$Taxonomy)

Core_Taxa_Isolate_Summary <- Core_Taxa_Isolate_Summary[,c(1:5,10)]

write.table(file="CoreTaxa_BlastMatches.txt", x=Core_Taxa_Isolate_Summary, sep="\t", quote=FALSE, row.names=FALSE)

# Phyllosphere taxa
sum(rowSums(otu_rare.PA[,map_16S$source=="phyllosphere"])>0)
# Unique Switchgrass Taxa (537 are unique to switchgrass)
sum(rowSums(otu_rare.PA[,map_16S$source=="phyllosphere"&map_16S$plant=="switchgrass"])>0&rowSums(otu_rare.PA[,map_16S$source=="phyllosphere"&map_16S$plant=="miscanthus"])==0)
#unique to Miscanthus (295 are unique to miscanthus)
sum(rowSums(otu_rare.PA[,map_16S$source=="phyllosphere"&map_16S$plant=="miscanthus"])>0&rowSums(otu_rare.PA[,map_16S$source=="phyllosphere"&map_16S$plant=="switchgrass"])==0)
# Shared between the two (323 are shared) 
sum(rowSums(otu_rare.PA[,map_16S$source=="phyllosphere"&map_16S$plant=="miscanthus"])>0&rowSums(otu_rare.PA[,map_16S$source=="phyllosphere"&map_16S$plant=="switchgrass"])>0)
537 + 295 + 323


################################
### Switchgrass 2016 vs 2017 ###
################################

switch.core.2016 <- c("OTU7", "OTU47", "OTU192", "OTU23", "OTU2", "OTU22", "OTU4", "OTU2771", "OTU14", "OTU537", "OTU90", "OTU41", "OTU5", "OTU430", "OTU6", "OTU15", "OTU995", "OTU6096", "OTU1334", "OTU18", "OTU842", "OTU21", "OTU10", "OTU32")

switch.core.2017 <- c("OTU7", "OTU47", "OTU192", "OTU23", "OTU2", "OTU22", "OTU4", "OTU2771", "OTU86", "OTU14", "OTU84", "OTU537", "OTU90", "OTU4223", "OTU41", "OTU519", "OTU5", "OTU430", "OTU6", "OTU18", "OTU21", "OTU10")

switch.occ.2016 <- rowSums(1*(switch.rare.otu.2016>0))/ncol(switch.rare.otu.2016)

switch.occ.2017 <- rowSums(1*(switch.rare.otu.2017>0))/ncol(switch.rare.otu.2017)

switch.abund.2016 <- rowSums(switch.rare.otu.2016)/sum(colSums(switch.rare.otu.2016))

switch.abund.2017 <- rowSums(switch.rare.otu.2017)/sum(colSums(switch.rare.otu.2017))

plot(log10(switch.abund.2016), log10(switch.abund.2017))

plot(switch.occ.2016, switch.occ.2017)

Symbols.switch.core <- rep(1, nrow(switch.rare.otu.2016))
Symbols.switch.core[row.names(switch.rare.otu.2016)%in%switch.core.2016] <- 0
Symbols.switch.core[row.names(switch.rare.otu.2017)%in%switch.core.2017] <- 2
Symbols.switch.core[row.names(switch.rare.otu.2017)%in%switch.core.2017&row.names(switch.rare.otu.2016)%in%switch.core.2016] <- 5


plot(switch.occ.2016, switch.occ.2017, pch=Symbols.switch.core)

plot(log10(switch.abund.2016), log10(switch.abund.2017), pch=Symbols.switch.core)

setEPS()
postscript("FigureS#_SwitchgrassSeasonalOcc.eps", width=8, height=4, paper="special")
par(ps = 8, cex = 1, cex.main = 1, mfrow=c(1,2))
plot(switch.occ.2016, switch.occ.2017, pch=Symbols.switch.core)
plot(log10(switch.abund.2016), log10(switch.abund.2017), pch=Symbols.switch.core)
dev.off()

###############################
#### Variance Partitioning ####
###############################
nrow(plant_rare_otu)
colSums(plant_rare_otu)

# Make Big huge table for variance partitioning
var_part_map <- map_16S[map_16S$source=="phyllosphere"&map_16S$year==2016,c("precipitation", "Air_temp_mean", "air_temp_max", "Air_Temp_Min", "Air_Pressure", "RH", "AH", "Wind_Speed_Mean", "Solar_Radiation", "PAR", "soil_temp_5_cm_bare_avg", "soil_temp_5_cm_sod_avg", "LDMC_mg_per_g", "nitrogen_percent", "carbon_percent", "carbon_per_nitrogen", "height_mean_cm","pH", "P_ppm", "K_ppm", "Ca_ppm", "Mg_ppm", "organic_matter", "NO3N_ppm", "NH4_ppm", "soil_moisture_percent", "soil_temp_10cm")]


colinearity_time <- data.frame(Variable=rep(NA, ncol(var_part_map)), T=rep(NA, ncol(var_part_map)), pvalue=rep(NA, ncol(var_part_map)), estimate=rep(NA, ncol(var_part_map)))
for(i in 1:ncol(var_part_map)){
  x <- cor.test(map_16S$time_numeric[map_16S$source=="phyllosphere"&map_16S$year==2016], var_part_map[,i])
  colinearity_time$Variable[i] <- colnames(var_part_map)[i]
  colinearity_time$T[i] <- x$statistic 
  colinearity_time$pvalue[i] <- x$p.value
  colinearity_time$estimate[i] <- x$estimate
  }

sig_colinear <- colinearity_time[colinearity_time$pvalue<0.05, 1]

#Use below dataset for ABIOTIC
reduced_var_part_ABIOTIC <- var_part_map[,!colnames(var_part_map)%in%sig_colinear]
ABIOTIC <- reduced_var_part_ABIOTIC


reduced_sample_names <- unlist(strsplit(colnames(plant_rare_otu.2016), split="_"))

reduced_sample_names <- rep(NA, ncol(plant_rare_otu.2016))
for(i in 1:ncol(plant_rare_otu.2016)) {
  reduced_sample_names[i] <-  unlist(strsplit(colnames(plant_rare_otu.2016)[i], split="_"))[1] 
}

spatial <- read.table("LongSpatialDistanceMatrix.txt", sep="\t", header = TRUE, row.names = NULL, stringsAsFactors = FALSE)

Big_Spatial_Distance_Matrix <- NULL
plots <- unique(reduced_sample_names)
for(i in 1:length(plots)){
  z <- matrix(nrow = sum(reduced_sample_names==plots[i]), ncol=length(reduced_sample_names))
  row.names(z) <- colnames(plant_rare_otu.2016)[reduced_sample_names==plots[i]]
  spatial_x <- spatial[spatial$Site1==plots[i],]
  for (w in 1:length(plots)){
    z[,reduced_sample_names==plots[w]] <- spatial_x[w,3]
  }
  Big_Spatial_Distance_Matrix <- rbind(Big_Spatial_Distance_Matrix,z)
}
colnames(Big_Spatial_Distance_Matrix) <- row.names(Big_Spatial_Distance_Matrix)

# Use below dataset for SPACE
SPACE <- Big_Spatial_Distance_Matrix

### Is Space significant for community similarity? Answer== NO
mantel(Big_Spatial_Distance_Matrix, plant.dist.2016)

SPACE.pcnm <- pcnm(Big_Spatial_Distance_Matrix,threshold = .5)
SPACE.pcnm$vectors

# Use below for HOST
HOST <- model.matrix(~plant, map_16S[map_16S$source=="phyllosphere"&map_16S$year==2016])[,-1]

HOST <- 1*(map_16S[map_16S$source=="phyllosphere"&map_16S$year==2016, "plant"]=="switchgrass")

# Use below for TIME 
TIME <- map_16S[map_16S$source=="phyllosphere"&map_16S$year==2016, "time_numeric"]


parting_pcnm<-varpart(Y=t(plant_rare_otu.2016), HOST, ABIOTIC, TIME, transfo = "total")
parting_pcnm
plot(parting_pcnm, main="PCNM", Xnames=c("HOST", "ABIOTIC", "TIME"), digits=2, cutoff=-1)


CorePartition <- varpart(Y=plant_rare_otu.2016.coretaxa.rel.dist, HOST, ABIOTIC, TIME)
Core2016Partition <- varpart(Y=plant_rare_otu.2016.coretaxa.2016.rel.dist, HOST, ABIOTIC, TIME)

Core2016Partition.norel <- varpart(Y=plant_rare_otu.2016.coretaxa.2016.dist, HOST, ABIOTIC, TIME)


plot(CorePartition, Xnames=c("HOST", "ABIOTIC", "TIME"), cutoff=-1)
plot(Core2016Partition, Xnames=c("HOST", "ABIOTIC", "TIME"), cutoff=-1)
plot(Core2016Partition.norel, Xnames=c("HOST", "ABIOTIC", "TIME"), cutoff=-1)


parting.dist.pcnm <- varpart(Y=plant.dist.2016, HOST, ABIOTIC, TIME)

setEPS()
postscript("FigureS5_VariancePartitioning.eps", width=5, height=5, paper="special")
par(ps = 8, cex = 1, cex.main = 1)
plot(parting.dist.pcnm,Xnames=c("HOST", "ABIOTIC", "TIME"), digits=2, cutoff=-1)
dev.off()

### Dead Code
Spec_accum_shapes.2016 <- c(rep(15,9), rep(16,8))
Spec_accum_lines.2016 <- c(rep(2,9), rep(1,8))
C <-  ggplot(Spec_accum.2016, aes(x=Date, y=Species_Accumulation)) + 
  geom_point(aes(color=Plant, shape=Spec_accum_shapes.2016)) +
  scale_shape_identity()+
  scale_linetype_identity()+
  geom_line(aes(color=Plant, linetype=Spec_accum_lines.2016))+
  labs(y="Total Observed Taxa", title ="C") + 
  scale_color_manual(values=c("darkgreen", "darkolivegreen3"))  + 
  guides(color=FALSE)+
  scale_x_date(date_breaks = "1 month", labels=date_format("%b"), limits= as.Date(c('2016-05-01', '2016-11-15')))
C

Spec_accum_shapes.2017 <- rep(16, nrow(Spec_accum.2017))
Spec_accum_lines.2017 <- rep(1, nrow(Spec_accum.2017))
C_2017 <- ggplot(Spec_accum.2017, aes(x=Date, y=Species_Accumulation)) + 
  geom_point(aes(color=Plant, shape=Spec_accum_shapes.2017)) +
  scale_shape_identity()+
  scale_linetype_identity()+
  geom_line(aes(color=Plant, linetype=Spec_accum_lines.2017))+
  labs(y="Total Observed Taxa", title ="C") + 
  scale_color_manual(values=c("darkolivegreen3", "darkgreen"))  + 
  guides(color=FALSE)+
  scale_x_date(date_breaks = "1 month", labels=date_format("%b"), limits= as.Date(c('2017-05-01', '2017-11-15')))


otu_soil_full <- otu[,map_16S$source=="soil"]
otu_soil_full.2016 <- otu[,map_16S$source=="soil"&map_16S$year=="2016"]
otu_soil_full.2017 <- otu[,map_16S$source=="soil"&map_16S$year=="2017"]

plant_rare_otu.PA.2017 <- 1*(rowSums(plant_rare_otu.2017)>0)

otu_soil_full.PA.2017 <- 2*(rowSums(otu_soil_full.2017)>0) 

table(plant_rare_otu.PA.2017 + otu_soil_full.PA.2017)

##### Nejc's Updates to code

#####
#Table S1 - EnvFit
#####


whole.dist <- vegdist(t(otu_rare), method="bray")
switch.dist <- vegdist(t(otu_rare[,map_16S$source=="phyllosphere"&map_16S$plant=="switchgrass"]), method="bray")
plant.dist <- vegdist(t(otu_rare[,map_16S$source=="phyllosphere"]), method="bray")
soil.dist <- vegdist(t(otu_rare[,map_16S$source=="soil"]), method="bray")

map_16S$timepoint <- rep(1,nrow(map_16S))
map_16S$timepoint[map_16S$sampling_date=="2016-05-31"]<- 2
map_16S$timepoint[map_16S$sampling_date=="2017-06-05"]<- 2
map_16S$timepoint[map_16S$sampling_date=="2016-06-20"]<- 3
map_16S$timepoint[map_16S$sampling_date=="2017-06-26"]<- 3
map_16S$timepoint[map_16S$sampling_date=="2016-07-12"] <- 4
map_16S$timepoint[map_16S$sampling_date=="2017-07-17"]<- 4
map_16S$timepoint[map_16S$sampling_date=="2016-08-01"] <- 5
map_16S$timepoint[map_16S$sampling_date=="2017-08-07"] <- 5
map_16S$timepoint[map_16S$sampling_date=="2016-08-22"] <- 6
map_16S$timepoint[map_16S$sampling_date=="2017-08-28"]<- 6
map_16S$timepoint[map_16S$sampling_date=="2016-09-12"] <- 7
map_16S$timepoint[map_16S$sampling_date=="2017-09-18"]<- 7
map_16S$timepoint[map_16S$sampling_date=="2016-10-03"]<- 8
map_16S$timepoint[map_16S$sampling_date=="2016-11-07"]<- 9

map.whole <- map_16S
map.switch <- map_16S[map_16S$source=="phyllosphere" & map_16S$plant=="switchgrass",]
map.plant <- map_16S[map_16S$source=="phyllosphere",]
map.soil <- map_16S[map_16S$source=="soil",]

rownames(t(otu_rare[,map_16S$source=="phyllosphere"])) == map.plant$sequence_name
rownames(t(otu_rare[,map_16S$source=="soil"])) == map.soil$sequence_name

Soil_Chemistry_whole<- map.whole[,c('timepoint',"pH", "P_ppm", "K_ppm", "Ca_ppm", "Mg_ppm", "organic_matter", "NO3N_ppm", "NH4_ppm")]
for (i in 1:ncol(Soil_Chemistry_whole)){
  Soil_Chemistry_whole[,i] <- as.numeric(Soil_Chemistry_whole[,i])
}
Leaf_Chemistry.whole <- map.whole[,c("precipitation", "Air_temp_mean", "air_temp_max", "Air_Temp_Min", "Air_Pressure", "RH", "AH", "Wind_Speed_Mean", "Solar_Radiation", "PAR", "soil_temp_5_cm_bare_avg", "sampling_Rdate")]
numeric_map.whole <- cbind(Soil_Chemistry_whole,Leaf_Chemistry.whole)

Soil_Chemistry_soil<- map.soil[,c('timepoint',"pH", "P_ppm", "K_ppm", "Ca_ppm", "Mg_ppm", "organic_matter", "NO3N_ppm", "NH4_ppm")]
for (i in 1:ncol(Soil_Chemistry_soil)){
  Soil_Chemistry_soil[,i] <- as.numeric(Soil_Chemistry_soil[,i])
}
Leaf_Chemistry.soil <- map.soil[,c("precipitation", "Air_temp_mean", "air_temp_max", "Air_Temp_Min", "Air_Pressure", "RH", "AH", "Wind_Speed_Mean", "Solar_Radiation", "PAR", "soil_temp_5_cm_bare_avg", "sampling_Rdate")]
numeric_map.soil <- cbind(Soil_Chemistry_soil,Leaf_Chemistry.soil)

Soil_Chemistry_plant<- map.plant[,c('timepoint',"pH", "P_ppm", "K_ppm", "Ca_ppm", "Mg_ppm", "organic_matter", "NO3N_ppm", "NH4_ppm")]
for (i in 1:ncol(Soil_Chemistry_plant)){
  Soil_Chemistry_plant[,i] <- as.numeric(Soil_Chemistry_plant[,i])
}
Leaf_Chemistry.plant <- map.plant[,c("precipitation", "Air_temp_mean", "air_temp_max", "Air_Temp_Min", "Air_Pressure", "RH", "AH", "Wind_Speed_Mean", "Solar_Radiation", "PAR", "soil_temp_5_cm_bare_avg", "sampling_Rdate",
                                     'carbon_percent', 'nitrogen_percent', 'carbon_per_nitrogen')]
numeric_map.plant <- cbind(Soil_Chemistry_plant,Leaf_Chemistry.plant)

Soil_Chemistry_switch<- map.switch[,c('timepoint',"pH", "P_ppm", "K_ppm", "Ca_ppm", "Mg_ppm", "organic_matter", "NO3N_ppm", "NH4_ppm")]
for (i in 1:ncol(Soil_Chemistry_switch)){
  Soil_Chemistry_switch[,i] <- as.numeric(Soil_Chemistry_switch[,i])
}
Leaf_Chemistry.switch <- map.switch[,c("precipitation", "Air_temp_mean", "air_temp_max", "Air_Temp_Min", "Air_Pressure", "RH", "AH", "Wind_Speed_Mean", "Solar_Radiation", "PAR", "soil_temp_5_cm_bare_avg", "sampling_Rdate",
                                       'carbon_percent', 'nitrogen_percent', 'carbon_per_nitrogen')]
numeric_map.switch <- cbind(Soil_Chemistry_switch,Leaf_Chemistry.switch)

plant.pcoa <- cmdscale(plant.dist, eig=TRUE)
soil.pcoa <- cmdscale(soil.dist, eig=TRUE)
switch.pcoa <- cmdscale(switch.dist, eig=TRUE)
whole.pcoa <- cmdscale(whole.dist, eig=TRUE)

env.plant <- envfit(plant.pcoa, numeric_map.plant)
env.soil <- envfit(soil.pcoa, numeric_map.soil)
env.switch <- envfit(switch.pcoa, numeric_map.switch)
env.whole <- envfit(whole.pcoa, numeric_map.whole)

EnvFit.plant <- as.data.frame(env.plant$vectors$arrows)
EnvFit.plant$Rsquared <- env.plant$vectors$r
EnvFit.plant$pvalue <- env.plant$vectors$pvals
EnvFit.plant$desc <- 'phyllosphere'
EnvFit.plant$variable <- rownames(EnvFit.plant)
rownames(EnvFit.plant) <- NULL
Sig_plant_EnvFit <- EnvFit.plant[EnvFit.plant$pvalue<0.05,]

EnvFit.soil <- as.data.frame(env.soil$vectors$arrows)
EnvFit.soil$Rsquared <- env.soil$vectors$r
EnvFit.soil$pvalue <- env.soil$vectors$pvals
EnvFit.soil$desc <- 'soil'
EnvFit.soil$variable <- rownames(EnvFit.soil)
rownames(EnvFit.soil) <- NULL
Sig_soil_EnvFit <- EnvFit.soil[EnvFit.soil$pvalue<0.05,]

EnvFit.switch <- as.data.frame(env.switch$vectors$arrows)
EnvFit.switch$Rsquared <- env.switch$vectors$r
EnvFit.switch$pvalue <- env.switch$vectors$pvals
EnvFit.switch$desc <- 'switchgrass'
EnvFit.switch$variable <- rownames(EnvFit.switch)
rownames(EnvFit.switch) <- NULL
Sig_switch_EnvFit <- EnvFit.switch[EnvFit.switch$pvalue<0.05,]

EnvFit.whole <- as.data.frame(env.whole$vectors$arrows)
EnvFit.whole$Rsquared <- env.whole$vectors$r
EnvFit.whole$pvalue <- env.whole$vectors$pvals
EnvFit.whole$desc <- 'all'
EnvFit.whole$variable <- rownames(EnvFit.whole)
rownames(EnvFit.whole) <- NULL
Sig_whole_EnvFit <- EnvFit.whole[EnvFit.whole$pvalue<0.05,]

All_EnvFit <- rbind(Sig_whole_EnvFit, Sig_soil_EnvFit, Sig_plant_EnvFit, Sig_switch_EnvFit)
library(reshape2)
envfit_data <- acast(All_EnvFit, factor(desc)~variable, value.var="Rsquared")
envfit_data

write.table(file = "~/Dropbox/GLBRC_16S/Table1_EnvFit.txt", x=envfit_data, sep="\t", quote=FALSE, col.names = TRUE, row.names = TRUE)


#####
#Figure1ABC
#####
soil.pcoa <- cmdscale(soil.dist, eig=TRUE)
whole.pcoa <- cmdscale(whole.dist, eig=TRUE)
ax1.whole <- whole.pcoa$eig[1]/sum(whole.pcoa$eig)
ax2.whole <- whole.pcoa$eig[2]/sum(whole.pcoa$eig)

Unique_Dates_soil <- unique(map.soil$sampling_week)[order(unique(map.soil$sampling_week))]
Date_Size_soil <- rep(1, nrow(map.soil))
point_sizes_soil <- seq(from=1, to=3, length.out = length(Unique_Dates_soil))
for(i in 1:length(Unique_Dates_soil)){
  Date_Size_soil[map_16S$sampling_week==Unique_Dates_soil[i]] <- point_sizes_soil[i]
}

soil_plot_colors <- rep("darkolivegreen3", nrow(map.soil))
soil_plot_colors[map.soil$plant=="miscanthus" & map.soil$year==2016] <- "darkgreen"
soil_plot_colors[map.soil$plant=="switchgrass" & map.soil$year==2017] <- "darkolivegreen1"
soil_plot_colors[map.soil$source=="soil" & map.soil$year==2016] <- "burlywood4"
soil_plot_colors[map.soil$source=="soil"&map.soil$plant=="switchgrass"&map.soil$year==2016] <- "burlywood"
soil_plot_colors[map.soil$source=="soil"&map.soil$year==2017] <- "brown"
soil_plot_shapes <- rep(15, nrow(map.soil))
soil_plot_shapes[map.soil$plant=="switchgrass"] <- 16
soil_plot_shapes[map.soil$plant=="switchgrass"& map.soil$source=="soil" & map.soil$treatment=="nitrogen free"] <- 1 
soil_plot_shapes[map.soil$plant=="miscanthus"& map.soil$source=="soil"& map.soil$treatment=="nitrogen free"] <- 0 

Unique_Dates_plant <- unique(map.plant$sampling_week)[order(unique(map.plant$sampling_week))]
Date_Size_plant <- rep(1, nrow(map.plant))
point_sizes_plant <- seq(from=1, to=3, length.out = length(Unique_Dates_plant))
for(i in 1:length(Unique_Dates_plant)){
  Date_Size_plant[map_16S$sampling_week==Unique_Dates_plant[i]] <- point_sizes_plant[i]
}


plant_plot_colors <- rep("darkolivegreen3", nrow(map.plant))
plant_plot_colors[map.plant$plant=="miscanthus" & map.plant$year==2016] <- "darkgreen"
plant_plot_colors[map.plant$plant=="switchgrass" & map.plant$year==2017] <- "darkolivegreen1"
plant_plot_colors[map.plant$source=="soil" & map.plant$year==2016] <- "burlywood4"
plant_plot_colors[map.plant$source=="soil"&map.plant$plant=="switchgrass"&map.plant$year==2016] <- "burlywood"
plant_plot_colors[map.plant$source=="soil"&map.plant$year==2017] <- "brown"
plant_plot_shapes <- rep(15, nrow(map.plant))
plant_plot_shapes[map.plant$plant=="switchgrass"] <- 16
plant_plot_shapes[map.plant$plant=="switchgrass"& map.plant$source=="soil" & map.plant$treatment=="nitrogen free"] <- 1 
plant_plot_shapes[map.plant$plant=="miscanthus"& map.plant$source=="soil"& map.plant$treatment=="nitrogen free"] <- 0 

env.whole <- envfit(whole.pcoa, numeric_map.whole[,-which(names(numeric_map.whole)=='sampling_Rdate')])

Unique_Dates_whole <- unique(map.whole$sampling_week)[order(unique(map.whole$sampling_week))]
Date_Size_whole <- rep(1, nrow(map.whole))
point_sizes_whole <- seq(from=1, to=3, length.out = length(Unique_Dates_whole))
for(i in 1:length(Unique_Dates_whole)){
  Date_Size_whole[map_16S$sampling_week==Unique_Dates_whole[i]] <- point_sizes_whole[i]
}

whole_plot_colors <- rep("darkolivegreen3", nrow(map.whole))
whole_plot_colors[map.whole$plant=="miscanthus" & map.whole$year==2016] <- "darkgreen"
whole_plot_colors[map.whole$plant=="switchgrass" & map.whole$year==2017] <- "darkolivegreen1"
whole_plot_colors[map.whole$source=="soil" & map.whole$year==2016] <- "burlywood4"
whole_plot_colors[map.whole$source=="soil"&map.whole$plant=="switchgrass"&map.whole$year==2016] <- "burlywood"
whole_plot_colors[map.whole$source=="soil"&map.whole$year==2017] <- "brown"
whole_plot_shapes <- rep(15, nrow(map.whole))
whole_plot_shapes[map.whole$plant=="switchgrass"] <- 16
whole_plot_shapes[map.whole$plant=="switchgrass"& map.whole$source=="soil" & map.whole$treatment=="nitrogen free"] <- 1 
whole_plot_shapes[map.whole$plant=="miscanthus"& map.whole$source=="soil"& map.whole$treatment=="nitrogen free"] <- 0 


plant.pcoa <- cmdscale(plant.dist, eig=TRUE)
soil.pcoa <- cmdscale(soil.dist, eig=TRUE)
env.plant <- envfit(plant.pcoa, numeric_map.plant[,-which(names(numeric_map.plant)=='sampling_Rdate')])
env.soil <- envfit(soil.pcoa, numeric_map.soil[,-which(names(numeric_map.soil)=='sampling_Rdate')])
ax1.soil <- soil.pcoa$eig[1]/sum(soil.pcoa$eig)
ax2.soil <- soil.pcoa$eig[2]/sum(soil.pcoa$eig)
ax1.plant <- plant.pcoa$eig[1]/sum(plant.pcoa$eig)
ax2.plant <- plant.pcoa$eig[2]/sum(plant.pcoa$eig)

env.soil.sig <- envfit(soil.pcoa, map.soil[,c("K_ppm", "pH", "Air_Temp_Min", "Ca_ppm", "organic_matter")])

#Fig1B
plot(soil.pcoa$points[,1], soil.pcoa$points[,2], cex = Date_Size_soil,  pch= soil_plot_shapes, col=soil_plot_colors, 
     xlab=paste("PCoA1: ",100*round(ax1.soil,3),"% var. explained",sep=""), 
     ylab=paste("PCoA2: ",100*round(ax2.soil,3),"% var. explained",sep="")) +
  title('B', adj=0) 

#Fig1A
plot(ylim = c(-0.4,0.4), xlim=c(-0.5,0.4),collapsed_data.pcoa$ax1.average, collapsed_data.pcoa$ax2.average, pch=plant_shapes, cex=Date_Size_plant, col= plot_colors, xlab=paste("PCoA1: ",100*round(ax1.plant,3),"% var. explained",sep=""), ylab=paste("PCoA2: ",100*round(ax2.plant,3),"% var. explained",sep="")) +
  arrows(collapsed_data.pcoa$ax1.average, collapsed_data.pcoa$ax2.average-collapsed_data.pcoa$ax1.sd, collapsed_data.pcoa$ax1.average, collapsed_data.pcoa$ax2.average+collapsed_data.pcoa$ax2.sd, angle=90,length= 0.05, code=3,lwd=0.75, col= plot_colors) +
  arrows(collapsed_data.pcoa$ax1.average- collapsed_data.pcoa$ax1.sd, collapsed_data.pcoa$ax2.average, collapsed_data.pcoa$ax1.average + collapsed_data.pcoa$ax1.sd, collapsed_data.pcoa$ax2.average, angle=90,length= 0.05, code=3,lwd=0.75, col= plot_colors)+
  plot(R40_collapsed.weather.ef, p=0.05, col="black")+
  plot(R40_LC.env, p=0.05, col="black", labels="nitrogen_percent") +
  title('A', adj=0)

#Fig1C
plot(whole.pcoa$points[,1], whole.pcoa$points[,2], cex= Date_Size_whole, pch= whole_plot_shapes, col=whole_plot_colors, 
     xlab=paste("PCoA1: ",100*round(ax1.whole.2016,3),"% var. explained",sep=""), 
     ylab=paste("PCoA2: ",100*round(ax2.whole.2016,3),"% var. explained",sep="")) +
  title('C', adj=0)
legend(-.342,-.19,box.lty=0,
       legend = c('Switchgrass soil (NF) 2016','Switchgrass soil (F) 2016',
                  'Switchgrass soil (NF) 2017','Switchgrass soil (NF) 2017',
                  'Switchgrass phyllosphere 2016','Switchgrass phyllosphere 2017',
                  'Miscanthus soil (NF) 2016','Miscanthus soil (F) 2016',
                  'Miscanthus phyllosphere 2016'), 
       pch=c(1,16,1,16,16,16,0,15,15), 
       col=c("burlywood", "burlywood", 
             "brown", "brown",
             "darkolivegreen3", "darkolivegreen1",
             'burlywood4', 'burlywood4',
             "darkgreen"),
       cex=0.75,
       pt.cex=1,
       title='Sample (treatment)') 
legend(-.28, -.1, box.lty = 0,
       pch=16,
       title='Sampling week',
       c('0', '3', '9'),
       horiz=TRUE,
       pt.cex=c(1, 1.666667, 3),
       cex=.75)

############
##FigS2
############
C <- ggplot(switch.disper.dates.2017, aes(x=dates, y=distance))+
  geom_boxplot(aes(group=dates),fill="darkolivegreen3")+
  them_bw()+
  annotate("text", x = as.Date("2017-11-07"), y = 0.25, label = "NA")+
  labs(y="Distance to Median", x="Date", title="C") + 
  scale_y_continuous(limits = c(0, .5))+
  scale_x_date(date_breaks = "1 month", labels=date_format("%b"), limits= as.Date(c('2017-05-01', '2017-11-15')))

B<- ggplot(switch_dispersion_dates.2016, aes(x=dates, y=distance))+
  geom_boxplot(aes(group=dates),fill="darkolivegreen3")+
  them_bw()+
  annotate("text", x = as.Date("2016-11-07"), y = 0.25, label = "NA")+
  labs(y="Distance to Median", x="Date", title="B") + 
  scale_y_continuous(limits = c(0, .5))+
  scale_x_date(date_breaks = "1 month", labels=date_format("%b"), limits= as.Date(c('2016-05-01', '2016-11-15')))

A<- ggplot(misc_disperion_dates, aes(x=dates, y=distance)) +
  geom_boxplot(aes(group=dates),fill="darkgreen")+
  them_bw()+
  labs(y="Distance to Median", x="Date", title="A") + 
  scale_y_continuous(limits = c(0, .5))+
  scale_x_date(date_breaks = "1 month", labels=date_format("%b"), limits= as.Date(c('2016-05-01', '2016-11-15'))) 

setEPS()
postscript("~/Dropbox/GLBRC_16S/Figures/FigureS2_v2.eps", width=4, height=8, paper="special")
par(ps = 8, cex = 1, cex.main = 1)
multiplot(A,B,C, cols=1)
dev.off()

#####
#FigureS4
#####
B <- ggplot(switch.Pielou.2016, aes(x=sampling_Rdate, y=value)) +  
  geom_boxplot(mapping=aes(group=sampling_Rdate), fill="darkolivegreen3")+  
  #geom_point(color="darkolivegreen3") +
  annotate("text", x = as.Date("2016-11-07"), y = .5, label = "NA")+
  labs(y="Pielou's Evenness", x="Date", title ="B") + 
  scale_y_continuous(limits = c(0, 1))+
  scale_x_date(date_breaks = "1 month", labels=date_format("%b"), limits= as.Date(c('2016-05-01', '2016-11-15')))

C <- ggplot(switch.Pielou.2017, aes(x=sampling_Rdate, y=value)) +  
  geom_boxplot(mapping=aes(group=sampling_Rdate), fill="darkolivegreen3")+  
  #geom_point(color="darkolivegreen3") +
  annotate("text", x = as.Date("2017-11-07"), y = .5, label = "NA")+
  labs(y="Pielou's Evenness", x="Date", title ="C") + 
  scale_y_continuous(limits = c(0, 1))+
  scale_x_date(date_breaks = "1 month", labels=date_format("%b"), limits= as.Date(c('2017-05-01', '2017-11-15')))

A <- ggplot(misc.Pielou.2016, aes(x=sampling_Rdate, y=value)) + 
  geom_boxplot(mapping=aes(group=sampling_Rdate), fill="darkgreen") + 
  #geom_point(color="darkgreen") + 
  labs(y="Pielou's Evenness", x="Date", title ="A") + 
  scale_y_continuous(limits = c(0, 1))+
  scale_x_date(date_breaks = "1 month", labels=date_format("%b"), limits= as.Date(c('2016-05-01', '2016-11-15'))) 

setEPS()
postscript("~/Dropbox/GLBRC_16S/Figures/FigureS4_v2.eps", width=4, height=8, paper="special")
par(ps = 8, cex = 1, cex.main = 1)
multiplot(A,B,C, cols=1)
dev.off()


#####
#Ternary plot
####
misc_rare_otu.2016 <- otu_rare[,map_16S$source=="phyllosphere"&map_16S$year==2016&map_16S$plant=='miscanthus']
switch_rare_otu.2016 <- otu_rare[,map_16S$source=="phyllosphere"&map_16S$year==2016&map_16S$plant=='switchgrass']

misc_taxa.2016 <- 3*(rowSums(misc_rare_otu.2016)>0)
switch_taxa.2016 <- 1*(rowSums(switch_rare_otu.2016)>0)
soil_taxa.2016 <- 2*(rowSums(soil_rare_otu.2016)>0)

table(misc_taxa.2016 + switch_taxa.2016 + soil_taxa.2016)  

# Abundance of shared taxa in soil samples
soil_both_taxa.2016 <- soil_rare_otu.2016[(misc_taxa.2016 + switch_taxa.2016 + soil_taxa.2016)==6,]

#Abundance of shared taxa in misc phyllosphere samples
misc_both_taxa.2016 <- misc_rare_otu.2016[(misc_taxa.2016 + switch_taxa.2016 + soil_taxa.2016)==6,]

#Abundance of shared taxa in switch phyllosphere samples
switch_both_taxa.2016 <- switch_rare_otu.2016[(misc_taxa.2016 + switch_taxa.2016 + soil_taxa.2016)==6,]


colSums(soil_both_taxa.2016)/colSums(soil_rare_otu.2016)
range(colSums(soil_both_taxa.2016)/colSums(soil_rare_otu.2016))

# Relative abundance of all shared taxa in soil samples
RA_st_soil.2016 <- colSums(soil_both_taxa.2016)/colSums(soil_rare_otu.2016)

# Average relative abudnace of all shared taxa in soil samples
ARA_st_soil.2016 <- RA_st_soil.2016/nrow(soil_both_taxa.2016)

# Shared Taxa's average relative abundance in soil samples
ST_ARA_Soil.2016<-rowSums(soil_both_taxa.2016)/sum(colSums(soil_rare_otu.2016))

# Shared Taxa's average relative abundance in misc phyllosphere samples
ST_ARA_M.2016<-rowSums(misc_both_taxa.2016)/sum(colSums(misc_rare_otu.2016))

# Shared Taxa's average relative abundance in switch phyllosphere samples
ST_ARA_S.2016<-rowSums(switch_both_taxa.2016)/sum(colSums(switch_rare_otu.2016))

ST_ARA_plotdata <- data.frame(OTU=row.names(soil_both_taxa.2016), SoilAbundance=ST_ARA_Soil.2016, MiscanthusAbundance=ST_ARA_M.2016, SwitchgrassAbundance=ST_ARA_S.2016)

misc_both_taxa_PA.2016 <- 1*(misc_both_taxa.2016>0)
switch_both_taxa_PA.2016 <- 1*(switch_both_taxa.2016>0)
soil_both_taxa_PA.2016 <- 1*(soil_both_taxa.2016>0)

ST_Occ_Misc.2016 <- rowSums(misc_both_taxa_PA.2016)/ncol(misc_both_taxa_PA.2016)
ST_Occ_Switch.2016 <- rowSums(switch_both_taxa_PA.2016)/ncol(switch_both_taxa_PA.2016)
ST_Occ_Soil.2016 <- rowSums(soil_both_taxa_PA.2016)/ ncol(soil_both_taxa_PA.2016) 

ST_Occ_plotdata.2016 <- data.frame(OTU=row.names(soil_both_taxa.2016), SoilOccupancy=ST_Occ_Soil.2016, MiscanthusOccupancy=ST_Occ_Misc.2016, SwitchgrassOccupancy=ST_Occ_Switch.2016)

core_taxa <- c("53","117","4","16","27","32","8","9","20","19","2","25","3","21","37","30","6","44","63","920","10","101","641","80")
switchgrass_core_taxa <- c("104","12","1751","157","108","357","61")
miscanthus_core_taxa <- c("83","149","22","40")

core_taxa.2016 <- c("OTU47","OTU10", "OTU18", "OTU21", "OTU23", "OTU2", "OTU430", "OTU6", "OTU41", "OTU7", "OTU14", "OTU2771", "OTU5", "OTU4", "OTU22", "OTU995", "OTU32")
core_taxa_misc.2016 <- c("OTU15", "OTU48", "OTU50", "OTU3994", "OTU995", "OTU519", "OTU1674", "OTU32")
core_taxa_switch.2016 <- c("OTU192", "OTU6096", "OTU537", "OTU90", "OTU1334", "OTU842")

FigS4_Colors <- rep("Black", length(ST_ARA_S.2016))
FigS4_Colors[row.names(switch_both_taxa_PA.2016)%in%core_taxa.2016] <- "Green"
FigS4_Colors[row.names(switch_both_taxa_PA.2016)%in%core_taxa_switch.2016] <- "darkolivegreen3"
FigS4_Colors[row.names(switch_both_taxa_PA.2016)%in%core_taxa_misc.2016] <- "darkgreen"

FigS4_Size <- rep(1, length(ST_ARA_P.2016))
FigS4_Size[row.names(plant_both_taxa_PA.2016)%in%c(core_taxa.2016, core_taxa_switch.2016, core_taxa_misc.2016)] <- 2

plant_occ <- data.frame(ST_Occ_Plant.2016)
plant_occ$OTU <- rownames(plant_occ)
rownames(plant_occ) <- c()
ST_Occ_plotdata.2016$OTU <- as.character(ST_Occ_plotdata.2016$OTU) 
plant_occ$OTU
Occ_ternary_data <- left_join(data.frame(ST_Occ_plotdata.2016), plant_occ)

#install.packages('ggtern')
library(ggtern)
ternaryOcc<- ggtern(data=Occ_ternary_data,aes(x=MiscanthusOccupancy, y=SoilOccupancy, z=SwitchgrassOccupancy, size=ST_Occ_Plant.2016, fill=ST_Occ_Plant.2016)) + 
  theme_classic()+
  theme_showarrows() +
  theme_rotate(180)+
  geom_point(pch=21)  + 
  theme(legend.position = 'top') +
  guides(size = 'none') +
  labs(x="Miscanthus",y="Soil",z="Switchgrass", fill = 'Phyllosphere occupancy') 

ST_ARA_plotdata$OTU <- as.character(ST_ARA_plotdata$OTU) 
plant_abun <- data.frame(ST_ARA_P.2016)
plant_abun$OTU <- rownames(plant_abun)
rownames(plant_abun) <- c()
Abun_ternary_data <- left_join(data.frame(ST_ARA_plotdata), plant_abun)

ternaryAbun <- ggtern(data=Abun_ternary_data,aes(x=MiscanthusAbundance, y=SoilAbundance, z=SwitchgrassAbundance, size=ST_ARA_P.2016, fill=ST_ARA_P.2016)) + 
  theme_classic()+
  theme_showarrows() +
  theme_rotate(180)+
  geom_point(pch=21)  + 
  theme(legend.position = 'top') +
  guides(size = 'none') +
  labs(x="Miscanthus",y="Soil",z="Switchgrass", fill = 'Phyllosphere \nrel. abundance') 

setEPS()
postscript("~/Dropbox/GLBRC_16S/Figures/ternary.eps", width=8, height=6, paper="special")
par(ps = 8, cex = 1, cex.main = 1)
multiplot(ternaryOcc,ternaryAbun,cols=2)
dev.off()

ST_Occ_plotdata.2016$SoilOccupancy+ST_Occ_plotdata.2016$MiscanthusOccupancy+ST_Occ_plotdata.2016$SwitchgrassOccupancy

ST_ARA_plotdata.2016 <- data.frame(OTU=row.names(plant_both_taxa.2016), SoilAbundance=ST_ARA_S.2016, PhyllosphereAbundance=ST_ARA_P.2016)
ST_Occ_plotdata.2016 <- data.frame(OTU=row.names(plant_both_taxa.2016), SoilOccupancy=ST_Occ_Soil.2016, PhyllosphereOccupancy=ST_Occ_Plant.2016)
ST_Occ_plotdata.2016$test <- .5

ggtern(data=ST_Occ_plotdata.2016,aes(x=SoilOccupancy, y=PhyllosphereOccupancy, z=test)) + 
  theme_bw()+
  theme_showarrows() +
  theme_rotate(180)+
  geom_point()  + 
  labs(x="Soil",y="Phyllo",z="test",title="Occupancy - example OTU4; soil=0.157,\n misc=0.769, switch=0.854")+
  scale_T_continuous(breaks = seq(0,1,0.2), labels = seq(0,1,0.2)) +
  scale_L_continuous(breaks = seq(0,1,0.2), labels = seq(0,1,0.2)) +
  scale_R_continuous(breaks = seq(0,1,0.2), labels = seq(0,1,0.2))



c <- as.data.frame(table(map.plant$sampling_date))
otu_146 <- otu_rare

otu_500 <- t(rrarefy(t(otu), 500))

otu_1000 <- t(rrarefy(t(otu), 1000))             

otu_5k <- t(rrarefy(t(otu)), 5000)
             
otu_10K <- t(rrarefy(t(otu), 10000))

otu_146_10k_Match <- otu_146[,colSums(otu_10K)>9999]
phyllo_otu_146_10k_Match <- otu_146[,colSums(otu_10K)>9999&map_16S$source=="phyllosphere"]

otu_500_clean <- otu_500[,colSums(otu_500)>499]
map_500_clean <- map_16S[colSums(otu_500)>499,]
map_500_clean_phyllo <- map_500_clean[map_500_clean$source=="phyllosphere",]
otu_500_clean_phyllo <- otu_500_clean[,colnames(otu_500_clean)%in%map_500_clean_phyllo$sequence_name]

d <- as.data.frame(table(map_500_clean_phyllo$sampling_date))


otu_1k_Clean <- otu_1000[,colSums(otu_1000)>999]
map_1k_clean <- map_16S[map_16S$sequence_name%in%colnames(otu_1k_Clean),]
map_1k_clean_phyllo <- map_1k_clean[map_1k_clean$source=="phyllosphere",]
otu_1k_Clean_phyllo <- otu_1k_Clean[,colnames(otu_1k_Clean)%in%map_1k_clean_phyllo$sequence_name]

map_1k_clean_phyllo_2016 <- map_1k_clean_phyllo[map_1k_clean_phyllo$Year==2016,]
Week_2016.df <- data.frame(sampling_date=unique(map_1k_clean_phyllo_2016$sampling_date)[order(unique(map_1k_clean_phyllo_2016$sampling_date))], Week=c(1:9))
map_1k_clean_phyllo_2016 <- left_join(map_1k_clean_phyllo_2016, Week_2016.df, by="sampling_date")

map_1k_clean_phyllo_2016_noMay <- map_1k_clean_phyllo_2016[map_1k_clean_phyllo_2016$Week>2,]

otu_1k_Clean_phyllo_2016_noMay <- otu_1k_Clean_phyllo[,colnames(otu_1k_Clean_phyllo)%in%map_1k_clean_phyllo_2016_noMay$sequence_name]

otu_1k_Clean_phyllo_2016_noMay.bc <- vegdist(t(otu_1k_Clean_phyllo_2016_noMay), method="bray")

adonis(otu_1k_Clean_phyllo_2016_noMay.bc~map_1k_clean_phyllo_2016_noMay$time_numeric)
adonis(otu_1k_Clean_phyllo_2016_noMay.bc~map_1k_clean_phyllo_2016_noMay$Week)
adonis(otu_1k_Clean_phyllo_2016_noMay.bc~map_1k_clean_phyllo_2016_noMay$plant)



anosim(otu_1k_Clean_phyllo_2016_noMay.bc, map_1k_clean_phyllo_2016_noMay$Week)
anosim(otu_1k_Clean_phyllo_2016_noMay.bc, map_1k_clean_phyllo_2016_noMay$plant)



otu_1k_Clean_phyllo.bc <- vegdist(t(otu_1k_Clean_phyllo), method = "bray")
a <- table(map_1k_clean_phyllo$sampling_date)
a <- as.data.frame(a)

otu_10K_Clean <- otu_10K[,colSums(otu_10K)>9999]
phyllo_otu_10K_Clean <- otu_10K[,colSums(otu_10K)>9999&map_16S$source=="phyllosphere"]


colnames(otu_146_10k_Match) == colnames(otu_10K_Clean)

otu_10K_Clean.bc <- vegdist(t(otu_10K_Clean), method = "bray")
otu_146_10k_Match.bc <- vegdist(t(otu_146_10k_Match), method="bray")

phyllo_otu_10K_Clean.bc <- vegdist(t(phyllo_otu_10K_Clean), method = "bray")
phyllo_otu_146_10k_Match.bc <- vegdist(t(phyllo_otu_146_10k_Match), method="bray")

mantel(otu_10K_Clean.bc, otu_146_10k_Match.bc)

mantel(phyllo_otu_10K_Clean.bc, phyllo_otu_146_10k_Match.bc)

map_10K_Clean <- map_16S[map_16S$sequence_name%in%colnames(otu_10K_Clean),]

map_10K_Clean_phyllo <- map_10K_Clean[map_10K_Clean$source=="phyllosphere",]

b <- as.data.frame(table(map_10K_Clean_phyllo$sampling_date))


Summary_rarefaction <- left_join(c,d, by="Var1")

Summary_rarefaction <- left_join(Summary_rarefaction, a, by="Var1")
Summary_rarefaction <- left_join(Summary_rarefaction, b, by="Var1")


colnames(Summary_rarefaction) <- c("Sampling_Date", "146_Reads","500_Reads", "1000_Reads", "10000_Reads")
write.table(x=Summary_rarefaction, file = "Summary_rarefaction_phyllo.txt", sep="\t", quote=FALSE, row.names = FALSE )

table(map_10K_Clean_phyllo$sampling_Rdate)

adonis(otu_10K_Clean.bc~map_10K_Clean$source)
anosim(otu_10K_Clean.bc,map_10K_Clean$source)

map_10K_Clean_phyllo_2016 <- map_10K_Clean_phyllo[map_10K_Clean_phyllo$Year==2016,]
unique(map_10K_Clean_phyllo_2016$sampling_date)[order(unique(map_10K_Clean_phyllo_2016$sampling_date))]

### Histogram of % mitochondria and chloroplasts in 2016 phyllosphere samples
hist(((1-(colSums(otu)/colSums(otu_CM)))*100)[map_16S$source=="phyllosphere"&map_16S$year==2016])
range(((1-(colSums(otu)/colSums(otu_CM)))*100)[map_16S$source=="phyllosphere"&map_16S$year==2016])


### Histogram of % mitochondria and chloroplasts in 2017 phyllosphere samples
hist(((1-(colSums(otu)/colSums(otu_CM)))*100)[map_16S$source=="phyllosphere"&map_16S$year==2017])
range(((1-(colSums(otu)/colSums(otu_CM)))*100)[map_16S$source=="phyllosphere"&map_16S$year==2017])


### Histogram of % mitochondria and chloroplasts in 2016 soil samples
hist(((1-(colSums(otu)/colSums(otu_CM)))*100)[map_16S$source=="soil"&map_16S$year==2016])
range(((1-(colSums(otu)/colSums(otu_CM)))*100)[map_16S$source=="soil"&map_16S$year==2016])


### Histogram of % mitochondria and chloroplasts in 2017 soil samples
hist(((1-(colSums(otu)/colSums(otu_CM)))*100)[map_16S$source=="soil"&map_16S$year==2017])
range(((1-(colSums(otu)/colSums(otu_CM)))*100)[map_16S$source=="soil"&map_16S$year==2017])




##################################
#####  End Jackson Analysis  #####
##################################
