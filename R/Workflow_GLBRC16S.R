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
otu <- read.table("R/InputFiles/table_combined_merged_trimmed_otus.txt",sep="\t", header=TRUE, stringsAsFactors = FALSE, row.names=1)

#Preping the environmental metadata
glbrc <- dbConnect(SQLite(), dbname="R/InputFiles/GLBRC_bioenergy_db.db" )

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
weather <- read.csv('R/InputFiles/kbs_weather_09212017.csv', encoding = 'UTF-8', na.strings= c("NA", " ", ""))
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
taxonomy <- read.csv('R/InputFiles/taxonomy_combined_merged_trimmed_otus.csv', header = T, row.names = 1, na.strings= c("NA", " ", ""))
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

tax_filtered <- tax_short%>%
  mutate(otu = rownames(tax_short)) %>%
  filter(!is.na(Phylum))

silva_bact_only <- read.csv('R/InputFiles/silva_bacteria_only_glbrc.csv', header=T)

keep_otus <- silva_bact_only %>%
  filter(lca_tax_slv != 'Unclassified;') %>%
  separate(lca_tax_slv, into=c("Kingdom", "Phylum", "Class", 
                               "Order", "Family", "Genus", "Species"), sep=";", remove=F) %>%
  mutate(taxonomy = lca_tax_slv) %>%
  filter(Kingdom!= 'Eukaryota') %>%
  select(-lca_tax_slv) %>%
  bind_rows(tax_filtered) %>%
  select(-taxonomy)

otu_filtered <- otu[rownames(otu) %in% keep_otus$otu,]

library(vegan)
set.seed(13)
#otu_rare <- t(rrarefy(t(otu), min(colSums(otu))))
otu_rare<- t(rrarefy(t(otu_filtered), 1000))
otu_soil_rare <- t(rrarefy(t(otu_soil), min(colSums(otu_soil))))
otu_141 <- t(rrarefy(t(otu_filtered), min(colSums(otu_filtered))))

otu.141 <- otu_141[,colSums(otu_141)>140]
otu_rare <- otu_rare[,colSums(otu_rare)>999]
map_16S <- map_16S[map_16S$sequence_name%in%colnames(otu_rare),]

map_16S.141 <- map_16S[map_16S$sample_name%in%colnames(otu_141)]

curve_colors <- rep("darkgreen", ncol(otu))
curve_colors[map_small$source=="phyllosphere"&map_small$plant=="switchgrass"] <- "darkolivegreen3"
curve_colors[map_small$source=="soil"&map_small$plant=="switchgrass"] <- "burlywood"
curve_colors[map_small$source=="soil"&map_small$plant=="miscanthus"] <- "burlywood4"

curve_colors_rare <- rep("darkgreen", ncol(otu))
curve_colors_rare[map_16S$source=="phyllosphere"&map_16S$plant=="switchgrass"] <- "darkolivegreen3"
curve_colors_rare[map_16S$source=="soil"&map_16S$plant=="switchgrass"] <- "burlywood"
curve_colors_rare[map_16S$source=="soil"&map_16S$plant=="miscanthus"] <- "burlywood4"

rarecurve(t(otu), step=1000, sample=max(colSums(otu)), label=FALSE, col = curve_colors)
rarecurve(t(otu), step=50, sample=10000, label=FALSE, col = curve_colors, xlim=c(0,1000))

rarecurve(t(otu_rare), step=50, label=FALSE, col= curve_colors_rare)

#rarecurve(t(otu_rare), step=20, sample = min(colSums(otu_rare)), label = FALSE)

#otu_rare <- otu_rare[,colSums(otu_rare)>999]

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
setwd('~/Documents/git/PAPER_GradySorensenStopnisek_InPrep/')
library(ggplot2)
library(dplyr)
library(tidyr)
library(vegan)
library(Heatplus)
library("limma")
library(ggrepel)
library(codyn)
library(gridExtra)
library(grid)
library(egg)

#Adding sampling week property
map_16S$sampling_week <- 0
map_16S$sampling_week[map_16S$sampling_date == '2017-05-15 EDT'] <- 1 
map_16S$sampling_week[map_16S$sampling_date == '2016-05-09 EDT'] <- 1

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

rownames(keep_otus) <- keep_otus$otu
keep_otus[] = lapply(keep_otus, blank2na, na.strings=c('','NA','na','N/A','n/a','NaN','nan'))
lastValue <- function(x) tail(x[!is.na(x)], 1)
last_taxons<- apply(keep_otus, 1, lastValue)
keep_otus$last_taxon <- last_taxons
head(keep_otus)
keep_otus$final_names <- paste(keep_otus$last_taxon, keep_otus$otu, sep=' - ')
tax_filtered <- keep_otus
###########################
#' Venn diagram - Figure S7
###########################
#' Subsetting the data to the crop and year
otu_rare=otu_rare[rowSums(otu_rare)>0,]

misc_otu <- otu_rare[,map_16S$plant=="miscanthus" & (map_16S$source=='phyllosphere')]
swit16_otu <- otu_rare[,map_16S$plant=="switchgrass" & (map_16S$source=='phyllosphere') & (map_16S$year==2016)]
swit17_otu <- otu_rare[,map_16S$plant=="switchgrass" & (map_16S$source=='phyllosphere') & (map_16S$year==2017)]

# make presence absence list from soil and plant into 1 & 0
swit16_otu_venn <- 1*(rowSums(swit16_otu)>0)
swit17_otu_venn<- 1*(rowSums(swit17_otu)>0)
misc16_otu_venn <- 1*(rowSums(misc_otu)>0)

#' plot Venn
venn_phyllo_data <- cbind(swit16_otu_venn,swit17_otu_venn, misc16_otu_venn)
colnames(venn_phyllo_data) <- c("Switchgrass 2016", "Switchgrass 2017", "Miscanthus 2016")
venn_phyllo_data=venn_phyllo_data[rowSums(venn_phyllo_data)>0,]
v_phyllo=vennCounts(venn_phyllo_data)
v_phyllo_2=round(v_phyllo[,"Counts"]/sum(v_phyllo[,"Counts"]),2) #calculate percentage of each group
vennDiagram(v_phyllo, circle.col = c('darkolivegreen3', 'darkolivegreen1','darkgreen'), lwd=6, cex=1.2, scale=F)

#############################
#Occupancy abundance analysis
#############################

#Comulative Occ_Abund
#Switch 2016 
swit16_otu <- swit16_otu[rowSums(swit16_otu)>0,]
swit_otu_PA <- 1*((swit16_otu>0)==1)
swit_otu_PA <- swit_otu_PA[rowSums(swit_otu_PA)>0,]
Occ_swit <- rowSums(swit_otu_PA)/ncol(swit_otu_PA)
swit_otu_rel <- decostand(swit16_otu, method="total", MARGIN=2)
Mean_abund_swit <- apply(swit_otu_rel, 1, mean)

#Creating df for plotting with ggplot and adding color code for the shared and unique OTUs
com_abund_swit <- rowSums(swit_otu_rel)
color_swit_top <- com_abund_swit
color_swit_top[] <- 'black'

swit_df_occ <- data.frame(otu=names(Occ_swit), occ=Occ_swit) 
swit_df_abun <- data.frame(otu=names(Mean_abund_swit), abun=log10(Mean_abund_swit))
switc_col <- data.frame(otu=names(color_swit_top), col=color_swit_top)
swit_occ_abun <- left_join(swit_df_abun, swit_df_occ, by='otu')
swit_occ_abun <- left_join(swit_occ_abun, switc_col, by='otu')

#Switch 2017 
swit17_otu <- swit17_otu[rowSums(swit17_otu)>0,]
swit17_otu_PA <- 1*((swit17_otu>0)==1)
swit17_otu_PA <- swit17_otu_PA[rowSums(swit17_otu_PA)>0,]
Mean_swit17_otu_PA <- apply(swit17_otu_PA, 1, mean)
Occ_swit17 <- rowSums(swit17_otu_PA)/ncol(swit17_otu_PA)
swit17_otu_rel <- decostand(swit17_otu, method="total", MARGIN=2)
Mean_abund_swit17 <- apply(swit17_otu_rel, 1, mean)

#Creating df for plotting with ggplot and adding color code for the shared and unique OTUs
com_abund_swit17 <- rowSums(swit17_otu_rel)
color_swit17_top <- com_abund_swit17
color_swit17_top[] <- 'black'

swit17_df_occ <- data.frame(otu=names(Occ_swit17), occ=Occ_swit17) 
swit17_df_abun <- data.frame(otu=names(Mean_abund_swit17), abun=log10(Mean_abund_swit17))
switc17_col <- data.frame(otu=names(color_swit17_top), col=color_swit17_top)
swit17_occ_abun <- left_join(swit17_df_abun, swit17_df_occ, by='otu')
swit17_occ_abun <- left_join(swit17_occ_abun, switc17_col, by='otu')

#Miscanthus 2016
misc_otu <- misc_otu[rowSums(misc_otu)>0,]
misc_otu_PA <- 1*((misc_otu>0)==1)
misc_otu_PA <- misc_otu_PA[rowSums(misc_otu)>0,]
Occ_misc <- rowSums(misc_otu_PA)/ncol(misc_otu_PA)
misc_otu <- misc_otu[rowSums(misc_otu)>0,]
misc_otu_rel <- decostand(misc_otu, method="total", MARGIN=2)
Mean_abund_misc <- apply(misc_otu_rel, 1, mean)

#Creating df for plotting with ggplot and adding color code for the shared and unique OTUs
com_abund_misc <- rowSums(misc_otu_rel)
color_misc_top <- com_abund_misc
color_misc_top[] <- 'black' 

misc_df_occ <- data.frame(otu=names(Occ_misc), occ=Occ_misc) 
misc_df_abun <- data.frame(otu=names(Mean_abund_misc), abun=log10(Mean_abund_misc))
misc_col <- data.frame(otu=names(color_misc_top), col=color_misc_top)
misc_occ_abun <- left_join(misc_df_abun, misc_df_occ, by='otu')
misc_occ_abun <- left_join(misc_occ_abun, misc_col, by='otu')

#' Creating a lsts of shared taxa
shared_v1 <- swit_occ_abun$otu[swit_occ_abun$otu %in% misc_occ_abun$otu]
shared <- shared_v1[shared_v1 %in% swit17_occ_abun$otu] 
shared_season_switch <- swit_occ_abun$otu[swit_occ_abun$otu %in% swit17_occ_abun$otu] 
shared_season_switch_only <- shared_season_switch[!(shared_season_switch %in% shared)]

##################################################################################
#Determining the core by the taxa contribution to the community structure varience
##################################################################################
sw <- swit16_otu
sw17<- swit17_otu
mis <- misc_otu

#Calculating OTU presence within the time points and presence in consecutive time points for each dataset
PresenceSum <- data.frame(otu = as.factor(row.names(sw)), sw) %>% 
  gather(sequence_name, abun, -otu) %>%
  left_join(map_16S[, c('sequence_name','source', 'plant', 'sampling_date', 'year')], by = 'sequence_name') %>%
  mutate(time_points=length(unique(sampling_date))) %>%
  group_by(sampling_date, otu, time_points) %>%
  summarise(date_presence=sum(abun>0)/length(abun),
            all_reps=length(abun),
            detect=ifelse(date_presence > 0, 1, 0)) %>%
  group_by(otu) %>%
  summarise(no_times=sum(detect),
            meanRep=mean(date_presence),
            Index=no_times*meanRep)

PresenceSum_sw17 <- data.frame(otu = as.factor(row.names(sw17)), sw17) %>% 
  gather(sequence_name, abun, -otu) %>%
  left_join(map_16S[, c('sequence_name','source', 'plant', 'sampling_date', 'year')], by = 'sequence_name') %>%
  mutate(time_points=length(unique(sampling_date))) %>%
  group_by(sampling_date, otu, time_points) %>%
  summarise(date_presence=sum(abun>0)/length(abun),
            all_reps=length(abun),
            detect=ifelse(date_presence > 0, 1, 0)) %>%
  group_by(otu) %>%
  summarise(no_times=sum(detect),
            meanRep=mean(date_presence),
            Index=no_times*meanRep)

PresenceSum_mi <- data.frame(otu = as.factor(row.names(mis)), mis) %>% 
  gather(sequence_name, abun, -otu) %>%
  left_join(map_16S[, c('sequence_name','source', 'plant', 'sampling_date', 'year')], by = 'sequence_name') %>%
  mutate(time_points=length(unique(sampling_date))) %>%
  group_by(sampling_date, otu, time_points) %>%
  summarise(date_presence=sum(abun>0)/length(abun),
            all_reps=length(abun),
            detect=ifelse(date_presence > 0, 1, 0)) %>%
  group_by(otu) %>%
  summarise(no_times=sum(detect),
            meanRep=mean(date_presence),
            Index=no_times*meanRep)

#' Creating OTU rankes by their abundacne and  weighted occupancy 
sw_ranked <- swit_occ_abun %>%
  left_join(PresenceSum, by='otu') %>%
  transmute(rank=(abun) + (occ*Index),
            otu=otu) %>%
  arrange(desc(rank))
sw_ranked$otu <- as.character(sw_ranked$otu) 

#' Calculating the B-C 
BCswAddition <- c()

for(i in sw_ranked$otu){
  otu_start=sw_ranked$otu[1]
  start_matrix <- as.matrix(sw[otu_start,])
  start_matrix <- t(start_matrix)
  x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]]- start_matrix[,x[2]]))/2000)
  x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse=' - '))
  df_s <- data.frame(x_names,x)
  names(df_s)[2] <- 1 
  BCswAddition <- rbind(BCswAddition,df_s)
  
  for(i in 2:length(sw_ranked$otu)){
    otu_add=sw_ranked$otu[i]
    add_matrix <- as.matrix(sw[otu_add,])
    add_matrix <- t(add_matrix)
    start_matrix <- rbind(start_matrix, add_matrix)
    x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]]-start_matrix[,x[2]]))/2000)
    x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse=' - '))
    df <- data.frame(x_names,x)
    names(df)[2] <- i 
    BCswAddition <- left_join(BCswAddition, df, by=c('x_names'))
  }
}

rownames(BCswAddition) <- BCswAddition$x_names
temp_BC <- BCswAddition
temp_BC$x_names <- NULL
temp_BC_matrix <- as.matrix(temp_BC)

plot_BC_ranked <- data.frame(rank = as.factor(row.names(t(temp_BC_matrix))),t(temp_BC_matrix)) %>% 
  gather(comparison, BC, -rank) %>%
  group_by(rank) %>%
  summarise(MeanBC=mean(BC)) %>%
  arrange(-desc(MeanBC)) %>%
  mutate(proportionBC=MeanBC/max(MeanBC))

Increase=plot_BC_ranked$MeanBC[-1]/plot_BC_ranked$MeanBC[-length(plot_BC_ranked$MeanBC)]
increaseDF <- data.frame(IncreaseBC=c(0,(Increase)), rank=factor(c(1:(length(Increase)+1))))
plot_BC_ranked <- left_join(plot_BC_ranked, increaseDF)

BC_plot <- ggplot(plot_BC_ranked[1:75,], aes(x=factor(plot_BC_ranked$rank[1:75], levels=plot_BC_ranked$rank[1:75]))) +
  geom_histogram(aes(y=IncreaseBC/5), stat='identity', col='darkolivegreen3', fill='darkolivegreen3', cex=1, alpha=.4) +
  geom_point(aes(y=MeanBC)) +
  labs(x=NULL, y='Mean Bray-Curtis distance') +
  theme_classic() + theme(strip.background = element_blank(),axis.text.x = element_text(size=7, angle=45)) +
  ylim(0,.57) +
  geom_vline(xintercept=c(37.5), lty=3, col='red', cex=.5) +
  scale_y_continuous(sec.axis = sec_axis(~.*5, name = "Increase in Bray-Curtis distance between points")) +
  
  theme_classic() + theme(strip.background = element_blank(), 
                          axis.text.x = element_text(size=7, angle=45))

#plot_BC_ranked[plot_BC_ranked$IncreaseBC>=1.02,]
sw_core_OTUs <- sw_ranked$otu[1:37]

#Switchgrass 2017
#' Creating OTU rankes by their abundacne and  weighted occupancy
sw17_ranked <- swit17_occ_abun %>%
  left_join(PresenceSum_sw17, by='otu') %>%
  transmute(rank=(abun) + occ*Index,
            otu=otu) %>%
  arrange(desc(rank))
sw17_ranked$otu <- as.character(sw17_ranked$otu) 

BCswAddition <- c()
for(i in sw17_ranked$otu){
  otu_start=sw17_ranked$otu[1]
  start_matrix <- as.matrix(sw17[otu_start,])
  start_matrix <- t(start_matrix)
  x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]]- start_matrix[,x[2]]))/2000)
  x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse=' - '))
  df_s <- data.frame(x_names,x)
  names(df_s)[2] <- 1 
  BCswAddition <- rbind(BCswAddition,df_s)
  
  for(i in 2:length(sw17_ranked$otu)){
    otu_add=sw17_ranked$otu[i]
    add_matrix <- as.matrix(sw17[otu_add,])
    add_matrix <- t(add_matrix)
    start_matrix <- rbind(start_matrix, add_matrix)
    x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]]-start_matrix[,x[2]]))/2000)
    x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse=' - '))
    df <- data.frame(x_names,x)
    names(df)[2] <- i 
    BCswAddition <- left_join(BCswAddition, df, by=c('x_names'))
  }
}

rownames(BCswAddition) <- BCswAddition$x_names
temp_BC <- BCswAddition
temp_BC$x_names <- NULL
temp_BC_matrix <- as.matrix(temp_BC)

plot_sw17_BC_ranked <- data.frame(rank = as.factor(row.names(t(temp_BC_matrix))),t(temp_BC_matrix)) %>% 
  gather(comparison, BC, -rank) %>%
  group_by(rank) %>%
  summarise(MeanBC=mean(BC)) %>%
  arrange(-desc(MeanBC)) %>%
  mutate(proportionBC=MeanBC/max(MeanBC))

Increase=plot_sw17_BC_ranked$MeanBC[-1]/plot_sw17_BC_ranked$MeanBC[-length(plot_sw17_BC_ranked$MeanBC)]
increaseDF <- data.frame(IncreaseBC=c(0,(Increase)), rank=factor(c(1:(length(Increase)+1))))
plot_sw17_BC_ranked <- left_join(plot_sw17_BC_ranked, increaseDF)

sw17_BC_plot <- ggplot(plot_sw17_BC_ranked[1:75,],aes(x=factor(plot_sw17_BC_ranked$rank[1:75], levels=plot_sw17_BC_ranked$rank[1:75]))) +
  geom_histogram(aes(y=IncreaseBC/5), stat='identity', col='darkolivegreen3', fill='white', cex=1) +
  geom_point(aes(y=MeanBC)) +  
  labs(x='Ranked OTUs', y=NULL) +
  theme_classic() + theme(strip.background = element_blank(), 
                          axis.text.x = element_text(size=7, angle=45))+
  ylim(0,.57)+
  geom_vline(xintercept=c(20.5), lty=3, col='red', cex=.5) +
  scale_y_continuous(sec.axis = sec_axis(~.*5))

plot_sw17_BC_ranked[plot_sw17_BC_ranked$IncreaseBC>=1.02,]
sw17_core_OTUs <- sw17_ranked$otu[1:20]

#Miscanthus 2016
#' Ranking the OTUs by their Occ and Abun' Creating OTU rankes by their abundacne and  weighted occupancy
mi_ranked <- misc_occ_abun %>%
  left_join(PresenceSum_mi, by='otu') %>%
  transmute(rank=(abun) + occ*Index,
            otu=otu) %>%
  arrange(desc(rank))
mi_ranked$otu <- as.character(mi_ranked$otu) 

BCswAddition <- c()
for(i in mi_ranked$otu){
  otu_start=mi_ranked$otu[1]
  start_matrix <- as.matrix(mis[otu_start,])
  start_matrix <- t(start_matrix)
  x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]]- start_matrix[,x[2]]))/2000)
  x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse=' - '))
  df_s <- data.frame(x_names,x)
  names(df_s)[2] <- 1 
  BCswAddition <- rbind(BCswAddition,df_s)
  
  for(i in 2:length(mi_ranked$otu)){
    otu_add=mi_ranked$otu[i]
    add_matrix <- as.matrix(mis[otu_add,])
    add_matrix <- t(add_matrix)
    start_matrix <- rbind(start_matrix, add_matrix)
    x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]]-start_matrix[,x[2]]))/2000)
    x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse=' - '))
    df <- data.frame(x_names,x)
    names(df)[2] <- i 
    BCswAddition <- left_join(BCswAddition, df, by=c('x_names'))
  }
}

rownames(BCswAddition) <- BCswAddition$x_names
temp_BC <- BCswAddition
temp_BC$x_names <- NULL
temp_BC_matrix <- as.matrix(temp_BC)

plot_mi_BC_ranked <- data.frame(rank = as.factor(row.names(t(temp_BC_matrix))),t(temp_BC_matrix)) %>% 
  gather(comparison, BC, -rank) %>%
  group_by(rank) %>%
  summarise(MeanBC=mean(BC)) %>%
  arrange(-desc(MeanBC)) %>%
  mutate(proportionBC=MeanBC/max(MeanBC))

Increase=plot_mi_BC_ranked$MeanBC[-1]/plot_mi_BC_ranked$MeanBC[-length(plot_mi_BC_ranked$MeanBC)]
increaseDF <- data.frame(IncreaseBC=c(0,(Increase)), rank=factor(c(1:(length(Increase)+1))))
plot_mi_BC_ranked <- left_join(plot_mi_BC_ranked, increaseDF)

mi_BC_plot <- ggplot(plot_mi_BC_ranked[1:75,],aes(x=factor(plot_mi_BC_ranked$rank[1:75], levels=plot_mi_BC_ranked$rank[1:75]))) +
  geom_histogram(aes(y=IncreaseBC/5), stat='identity', col='darkgreen', fill='darkgreen', alpha=.4) +
  geom_point(aes(y=MeanBC)) +    labs(x=NULL, y=NULL) +
  ylim(0,.57)+
  geom_vline(xintercept=c(24.5), lty=3, col='red', cex=.5)+
  theme_classic() + theme(strip.background = element_blank(), 
                          axis.text.x = element_text(size=7, angle=45),
                          axis.title.x = element_blank(), axis.title.y = element_blank()) +
  scale_y_continuous(sec.axis = sec_axis(~.*5))

plot_mi_BC_ranked[plot_mi_BC_ranked$IncreaseBC>=1.02,]
misc_core_OTUs <- mi_ranked$otu[1:24]

#*********************
# Figure S6
#*********************
setEPS()
postscript("~/Desktop/figure/figureS6.eps", width=6, height=4.5,pointsize=10, paper="special")
grid.draw(ggarrange(mi_BC_plot, BC_plot, sw17_BC_plot,nrow = 3))
dev.off()

###################
#' The core members
###################

core_list <- unique(c(misc_core_OTUs,sw_core_OTUs,sw17_core_OTUs))
list_of_core_taxa <- data.frame(colnames=c('plant', 'OTU'), OTU=misc_core_OTUs,sw_core_OTUs,sw17_core_OTUs)
list_of_core_taxa <- tax_filtered[tax_filtered$otu %in% core_list,]
list_of_core_taxa$color <- 'black'
list_of_core_taxa$color[list_of_core_taxa$Class == "c:Alphaproteobacteria"] <- 'green'
list_of_core_taxa$color[list_of_core_taxa$Class == "c:Betaproteobacteria"] <- 'blue'
list_of_core_taxa$color[list_of_core_taxa$Class == "c:Gammaproteobacteria"] <- 'purple'
list_of_core_taxa$color[list_of_core_taxa$Class == "c:Cytophagia"] <- 'orange'
list_of_core_taxa$color[list_of_core_taxa$Class == "c:Sphingobacteriia"] <- 'grey'
list_of_core_taxa$color[list_of_core_taxa$Class == "c:Actinobacteria"] <- 'red'
list_of_core_taxa$color[list_of_core_taxa$Class == "c:Flavobacteriia"] <- 'yellow'
list_of_core_taxa[c(4,10,11)]
write.csv(list_of_core_taxa[c(4,10,11)],'~/Desktop/core_colors.csv')


# core_list <- read.csv('../../../../../Desktop/core.csv')
# misc_core_OTUs <- as.character(core_list$otu[core_list$plant=='misc'])
# sw_core_OTUs <- as.character(core_list$otu[core_list$plant=='s16'])
# sw17_core_OTUs <- as.character(core_list$otu[core_list$plant=='s17'])

##################################################################
#Clustering to determine OTUs belonging to the plant growth stages
##################################################################
rel_otu_rare <- decostand(otu_rare, method="total", MARGIN=2)
abund_selected_all <- rel_otu_rare[rownames(rel_otu_rare) %in% core_list,]
selected_otus <- data.frame(otu = as.factor(row.names(abund_selected_all)), abund_selected_all) %>% gather(sequence_name, abun, -otu) %>% 
  left_join(map_16S[, c('sequence_name','rep','time_numeric', 'treatment' ,'source', 'plant', 'month', 'sampling_date', 
                        'carbon_per_nitrogen', 'nitrogen_percent', 'carbon_percent', 'year')], by = 'sequence_name') %>%
  left_join(tax_filtered, by='otu')

selected_otus_switch <- selected_otus[selected_otus$otu %in% sw_core_OTUs,]
selected_otus_misc <- selected_otus[selected_otus$otu %in% misc_core_OTUs,]
selected_otus_switch17 <- selected_otus[selected_otus$otu %in% sw17_core_OTUs,]

#Calculating the z score 
#***********************
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



#df with stats for the whole dataset per OTU
selected_otus_switch %>%
  filter(source == 'phyllosphere' & plant == 'switchgrass' & year == 2016) %>%
  group_by(final_names, otu) %>%
  dplyr::summarise(
    all_ab=mean(abun),
    all_sd=sd(abun)
  ) -> temp2

#combining df and calculating the z-score
z_df <- left_join(temp, temp2)
z_df %>% arrange(otu) %>%
  mutate(
    z_score=(rep_ab-all_ab)/all_sd
  ) %>%
  arrange(Class, Family) -> tmp3

#***************
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

#df with stats for the whole dataset per OTU
selected_otus_misc %>%
  filter(source == 'phyllosphere' & plant == 'miscanthus') %>%
  group_by(final_names, otu) %>%
  dplyr::summarise(
    all_ab=mean(abun),
    all_sd=sd(abun)
  ) -> temp5

#combining df and calculating the z-score
z_df_misc <- left_join(temp4, temp5)
z_df_misc %>% arrange(otu) %>%
  mutate(
    z_score=(rep_ab-all_ab)/all_sd
  ) %>%
  arrange(Class, Family) -> temp6

#*******************
#Switchgrass 2017
#OTU mean abundance per sampling time 
selected_otus_switch17 %>%
  filter(source == 'phyllosphere' & plant == 'switchgrass' & year == 2017) %>%
  group_by(sampling_date, otu, Class, Order, Family, Genus, final_names) %>%
  dplyr::summarise(n=sum(abun>0)/length(abun),
                   all=length(abun),
                   rep_ab=mean(abun),
                   sd_rep=sd(abun)
  ) %>%
  filter(n>0) -> temp7

#df with stats for the whole dataset per OTU
selected_otus_switch17 %>%
  filter(source == 'phyllosphere' & plant == 'switchgrass' & year == 2017) %>%
  group_by(otu, final_names) %>%
  dplyr::summarise(
    all_ab=mean(abun),
    all_sd=sd(abun)
  ) -> temp8

#combining df and calculating the z-score
z_df <- left_join(temp7, temp8)
z_df %>% arrange(otu) %>%
  mutate(
    z_score=(rep_ab-all_ab)/all_sd
  ) %>%
  arrange(Class, Family) -> tmp9

#Hierarchical clustering of the core taxa
switch16_core <- tmp3[,c(1,14,11)]
switch16_core <- as.data.frame(switch16_core)
switch16_core_wide <- spread(switch16_core, key='sampling_date', value='z_score')
switch16_core_wide[is.na(switch16_core_wide)] <- 0
rownames(switch16_core_wide) <- switch16_core_wide$otu
switch16_core_wide$otu <-  NULL
set.seed(20)
clusters_switch16 <- hclust(dist(switch16_core_wide),'complete')
memb_sw16 <- cutree(clusters_switch16, k=3)
sw16_dend <- plot(clusters_switch16, main=NULL)
rect.hclust(clusters_switch16, k=3)

switch17_core <- tmp9[,c(1,2,14)]
switch17_core <- as.data.frame(switch17_core)
switch17_core_wide <- spread(switch17_core, key='sampling_date', value='z_score')
switch17_core_wide[is.na(switch17_core_wide)] <- 0
rownames(switch17_core_wide) <- switch17_core_wide$otu
switch17_core_wide$otu <-  NULL
set.seed(21)
clusters_switch17 <- hclust(dist(switch17_core_wide),'complete')
memb_sw17 <- cutree(clusters_switch17, k=3)
sw17_dend <- plot(clusters_switch17, main=NULL)
rect.hclust(clusters_switch17, k=3)

misc16_core <- temp6[,c(1,11,14)]
misc16_core <- as.data.frame(misc16_core)
misc16_core_wide <- spread(misc16_core, key='sampling_date', value='z_score')
misc16_core_wide[is.na(misc16_core_wide)] <- 0
rownames(misc16_core_wide) <- misc16_core_wide$otu
misc16_core_wide$otu <-  NULL
set.seed(21)
clusters_misc16 <- hclust(dist(misc16_core_wide),'complete')
memb_misc <- cutree(clusters_misc16, k=4)
misc_dend <- plot(clusters_misc16, main=NULL)
rect.hclust(clusters_misc16, k=4)

#------------------------------------
# Figure S8

setEPS()
postscript('~/Desktop/figure/figureS8.eps', width=14, height=4, paper = 'special')
par(mfrow=c(1,3))
plot(clusters_misc16, main=NULL)
rect.hclust(clusters_misc16, k=4)
plot(clusters_switch16, main=NULL)
rect.hclust(clusters_switch16, k=3)
plot(clusters_switch17, main=NULL)
rect.hclust(clusters_switch17, k=3)
dev.off()
par(mfrow=c(1,1))

#**************************
# Occupancy abundance plots
#**************************

sw16_clusters <- data.frame(memb_sw16)
sw16_clusters$otu <-  rownames(sw16_clusters)
sw16_clusters$stage <- 'mid'
sw16_clusters$stage[sw16_clusters$memb_sw16==3] <- 'early'
sw16_clusters$stage[sw16_clusters$memb_sw16==2] <- 'late'

sw17_clusters <- data.frame(memb_sw17)
sw17_clusters$otu <- rownames(sw17_clusters)
sw17_clusters$stage[sw17_clusters$memb_sw17==2] <- 'late'
sw17_clusters$stage[sw17_clusters$memb_sw17==3] <- 'early'
sw17_clusters$stage[sw17_clusters$memb_sw17==1] <- 'mid'

misc_clusters <- data.frame(memb_misc)
misc_clusters$otu <-rownames(misc_clusters)
misc_clusters$stage[misc_clusters$memb_misc==2] <- 'late'
misc_clusters$stage[misc_clusters$memb_misc==3] <- 'early'
misc_clusters$stage[misc_clusters$memb_misc==1] <- 'mid'
misc_clusters$stage[misc_clusters$memb_misc==4] <- 'early'

swit_occ_abun$unique <- 'shared'
misc_occ_abun$unique <- 'shared'
swit17_occ_abun$unique <- 'shared'

swit_occ_abun$unique[!(swit_occ_abun$otu %in% shared_season_switch)] <- 'Switchgrass 2016 (n=418)'
swit_occ_abun$unique[(swit_occ_abun$otu %in% sw_core_OTUs)] <- 'Switchgrass 2016 (n=37)'
swit_occ_abun$unique[swit_occ_abun$otu %in% sw16_clusters$otu[sw16_clusters$stage=='mid']] <- 'mid'
swit_occ_abun$unique[swit_occ_abun$otu %in% sw16_clusters$otu[sw16_clusters$stage=='late']] <- 'late'
swit_occ_abun$unique[swit_occ_abun$otu %in% sw16_clusters$otu[sw16_clusters$stage=='early']] <- 'early'

misc_occ_abun$unique[!(misc_occ_abun$otu %in% shared)] <- 'Miscanthus 2016 (n=614)' 
misc_occ_abun$unique[(misc_occ_abun$otu %in% misc_core_OTUs)] <- 'Core (n=24)'
misc_occ_abun$unique[misc_occ_abun$otu %in% misc_clusters$otu[misc_clusters$stage=='mid']] <- 'mid'
misc_occ_abun$unique[misc_occ_abun$otu %in% misc_clusters$otu[misc_clusters$stage=='late']] <- 'late'
misc_occ_abun$unique[misc_occ_abun$otu %in% misc_clusters$otu[misc_clusters$stage=='early']] <- 'early'

swit17_occ_abun$unique[!(swit17_occ_abun$otu %in% shared_season_switch)] <- 'Switchgrass 2017 (n=357)'
swit17_occ_abun$unique[(swit17_occ_abun$otu %in% sw17_core_OTUs)] <- 'Core (n=20)'
swit17_occ_abun$unique[swit17_occ_abun$otu %in% sw17_clusters$otu[sw17_clusters$stage=='mid']] <- 'mid'
swit17_occ_abun$unique[swit17_occ_abun$otu %in% sw17_clusters$otu[sw17_clusters$stage=='late']] <- 'late'
swit17_occ_abun$unique[swit17_occ_abun$otu %in% sw17_clusters$otu[sw17_clusters$stage=='early']] <- 'early'

FigA <- ggplot(data=swit_occ_abun, aes(x=abun, y=occ, fill=unique)) +
  theme_bw()+
  geom_point(size=3, pch=21, alpha=.8) +
  scale_fill_manual(breaks=unique, values=c('red', 'blue', 'grey','white','darkolivegreen3')) +
  labs(x=paste('log(mean relative abundace per OTU)\n (n=',nrow(swit_df_abun),'OTUs)',sep=' '), y=paste('Occupancy (n=',ncol(swit_otu_PA),')', sep=' '), fill=NULL) +
  theme(legend.position = 'none',
        legend.background = element_rect(fill=alpha(0.1)),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_y_continuous(breaks=seq(0,1,.2)) +
  guides(fill = guide_legend(override.aes = list(alpha = 1)))

FigB <- ggplot(data=swit17_occ_abun, aes(x=abun, y=occ, fill=unique)) +
  theme_bw()+
  geom_point(size=3, pch=21, alpha=.8) +
  scale_fill_manual(breaks=unique, values=c('red', 'blue', 'grey','white', 'darkolivegreen1')) +
  labs(x=paste('log(mean relative abundace per OTU)\n (n=',nrow(swit17_occ_abun),' OTUs)',sep=''), y=paste('Occupancy (n=',ncol(swit17_otu_PA),')', sep=''), fill=NULL) +
  theme(legend.position = 'none',
        legend.background = element_rect(fill=alpha(0.1)),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_y_continuous(breaks=seq(0,1,.2)) +
  guides(fill = guide_legend(override.aes = list(alpha = 1)))

FigC <- ggplot(data=misc_occ_abun, aes(x=abun, y=occ, fill=unique)) +
  theme_bw()+
  geom_point(size=3, pch=21, alpha=.8) +
  scale_fill_manual(breaks=unique, values=c('red', 'blue', 'grey','darkgreen', 'white')) +
  labs(x=paste('log(mean relative abundace per OTU)\n (n=',nrow(misc_df_abun),' OTUs)',sep=''), y=paste('Occupancy (n=',ncol(misc_otu_PA),')', sep=''), fill=NULL) +
  theme(legend.position = 'none',
        legend.background = element_rect(fill=alpha(0.1)),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_y_continuous(breaks=seq(0,1,.2)) +
  guides(fill = guide_legend(override.aes = list(alpha = 1)))

grid.arrange(FigA, FigB, FigC, nrow=1)

#---------------------------------------------------
#Contribution of the core to community dissimilarity

#switch 2016
map_16S$sequence_name == colnames(otu_rare)
switch16full <- otu_rare[,map_16S$plant == 'switchgrass' & map_16S$year == 2016 & map_16S$source == 'phyllosphere']
switch16full <- switch16full[rowSums(switch16full)>0,]
switch16core <- switch16full[rownames(switch16full) %in% selected_otus_switch$otu,]

for(i in switch16full){
  week = c(3:8) 
  switch16_BC_full <- NULL
  for(i in week){
    samples = map_16S$sequence_name[map_16S$sampling_week== i]    #selecting the unique dates for pairwise calcualtions
    time_sub=switch16full[,colnames(switch16full) %in% samples]   #subsetting matrix to the selected dates
    x <- apply(combn(ncol(time_sub), 2), 2, function(x) sum(abs(time_sub[,x[1]]- time_sub[,x[2]]))/2000) 
    x_names <- apply(combn(ncol(time_sub), 2), 2, function(x) paste(colnames(time_sub)[x], collapse=' - '))
    df <- data.frame(x, x_names, i)
    switch16_BC_full <- rbind(switch16_BC_full,df)
  }
}

for(i in switch16core){
  week = c(3:8) 
  switch16_BC_core <- NULL
  for(i in week){
    samples = map_16S$sequence_name[map_16S$sampling_week== i]    #selecting the unique dates for pairwise calcualtions
    time_sub=switch16core[,colnames(switch16core) %in% samples]   #subsetting matrix to the selected dates
    x <- apply(combn(ncol(time_sub), 2), 2, function(x) sum(abs(time_sub[,x[1]]- time_sub[,x[2]]))/2000) 
    x_names <- apply(combn(ncol(time_sub), 2), 2, function(x) paste(colnames(time_sub)[x], collapse=' - '))
    df <- data.frame(x, x_names, i)
    switch16_BC_core <- rbind(switch16_BC_core,df)
  }
}

diffBC_switch16 <- left_join(switch16_BC_full,switch16_BC_core, by=c('x_names', 'i'))
diffBC_switch16$diff_BC <- diffBC_switch16$x.y/diffBC_switch16$x.x

switch16_BC_plot <- ggplot(diffBC_switch16, aes(x=factor(i), y=diff_BC)) +
  geom_violin(trim=FALSE,cex=1, color='darkolivegreen3',fill='darkolivegreen3', alpha=.2)+
  geom_jitter(width=0.1, color='darkolivegreen3')+
  ylim(0,1)+
  theme_bw()+
  labs(x=NULL, y=NULL)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text.x = element_blank(), panel.border = element_blank())
#switch 2017
switch17full <- otu_rare[,map_16S$plant == 'switchgrass' & map_16S$year == 2017 & map_16S$source == 'phyllosphere']
switch17full <- switch17full[rowSums(switch17full)>0,]
switch17core <- switch17full[rownames(switch17full) %in% selected_otus_switch17$otu,]
#unique(map_16S$sampling_week[map_16S$plant == 'switchgrass' & map_16S$year == 2017 & map_16S$source == 'phyllosphere'])

for(i in switch17full){
  week = c(1:7) 
  switch17_BC_full <- NULL
  for(i in week){
    samples = map_16S$sequence_name[map_16S$sampling_week== i]    #selecting the unique dates for pairwise calcualtions
    time_sub=switch17full[,colnames(switch17full) %in% samples]   #subsetting matrix to the selected dates
    x <- apply(combn(ncol(time_sub), 2), 2, function(x) sum(abs(time_sub[,x[1]]- time_sub[,x[2]]))/2000) 
    x_names <- apply(combn(ncol(time_sub), 2), 2, function(x) paste(colnames(time_sub)[x], collapse=' - '))
    df <- data.frame(x, x_names, i)
    switch17_BC_full <- rbind(switch17_BC_full,df)
  }
}

for(i in switch17core){
  week = c(1:7) 
  switch17_BC_core <- NULL
  for(i in week){
    samples = map_16S$sequence_name[map_16S$sampling_week== i]    #selecting the unique dates for pairwise calcualtions
    time_sub=switch17core[,colnames(switch17core) %in% samples]   #subsetting matrix to the selected dates
    x <- apply(combn(ncol(time_sub), 2), 2, function(x) sum(abs(time_sub[,x[1]]- time_sub[,x[2]]))/2000) 
    x_names <- apply(combn(ncol(time_sub), 2), 2, function(x) paste(colnames(time_sub)[x], collapse=' - '))
    df <- data.frame(x, x_names, i)
    switch17_BC_core <- rbind(switch17_BC_core,df)
  }
}

diffBC_switch17 <- left_join(switch17_BC_full,switch17_BC_core, by=c('x_names', 'i'))
diffBC_switch17$diff_BC <- diffBC_switch17$x.y/diffBC_switch17$x.x
mean(diffBC_switch17$diff_BC)
switch17_BC_plot <- ggplot(diffBC_switch17, aes(x=factor(i), y=diff_BC)) +
  geom_violin(trim=FALSE,cex=1, color='darkolivegreen3')+
  geom_jitter(width=0.1, color='darkolivegreen3', pch=21)+
  ylim(0,1)+
  theme_bw()+
  labs(x=NULL, y=NULL) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text.x = element_blank(), panel.border = element_blank()) 

#misc 2016
misc16full <- otu_rare[,map_16S$plant == 'miscanthus' & map_16S$year == 2016 & map_16S$source == 'phyllosphere']
misc16full <- misc16full[rowSums(misc16full)>0,]
misc16core <- misc16full[rownames(misc16full) %in% selected_otus_misc$otu,]
#unique(map_16S$sampling_week[map_16S$plant == 'miscanthus' & map_16S$year == 2016 & map_16S$source == 'phyllosphere'])

for(i in misc16full){
  week = c(3:9) 
  misc16_BC_full <- NULL
  for(i in week){
    samples = map_16S$sequence_name[map_16S$sampling_week== i]    #selecting the unique dates for pairwise calcualtions
    time_sub=misc16full[,colnames(misc16full) %in% samples]   #subsetting matrix to the selected dates
    x <- apply(combn(ncol(time_sub), 2), 2, function(x) sum(abs(time_sub[,x[1]]- time_sub[,x[2]]))/2000) 
    x_names <- apply(combn(ncol(time_sub), 2), 2, function(x) paste(colnames(time_sub)[x], collapse=' - '))
    df <- data.frame(x, x_names, i)
    misc16_BC_full <- rbind(misc16_BC_full,df)
  }
}

for(i in misc16core){
  week = c(3:9) 
  misc16_BC_core <- NULL
  for(i in week){
    samples = map_16S$sequence_name[map_16S$sampling_week== i]    #selecting the unique dates for pairwise calcualtions
    time_sub=misc16core[,colnames(misc16core) %in% samples]   #subsetting matrix to the selected dates
    x <- apply(combn(ncol(time_sub), 2), 2, function(x) sum(abs(time_sub[,x[1]]- time_sub[,x[2]]))/2000) 
    x_names <- apply(combn(ncol(time_sub), 2), 2, function(x) paste(colnames(time_sub)[x], collapse=' - '))
    df <- data.frame(x, x_names, i)
    misc16_BC_core <- rbind(misc16_BC_core,df)
  }
}

diffBC_misc16 <- left_join(misc16_BC_full,misc16_BC_core, by=c('x_names', 'i'))
diffBC_misc16$diff_BC <- diffBC_misc16$x.y/diffBC_misc16$x.x
mean(diffBC_misc16$diff_BC)
misc16_BC_plot <- ggplot(diffBC_misc16, aes(x=factor(i), y=diff_BC)) +
  geom_violin(trim=FALSE, cex=1, color='darkgreen', fill='darkgreen', 
              alpha=.2) +
  geom_jitter(width=0.1, color='darkgreen') +
  ylim(0,1) +
  expand_limits(x =factor(seq(1,9, by=1))) +
  theme_bw() +
  labs(x=NULL, y='Fraction of beta diversity\nwithin the time points') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text.x = element_blank(), panel.border = element_blank())

#----------------------------------------------------------
#Dynamics of the clusters assigned by the hclust function

abund_plot_sw16 <- data.frame(otu = as.factor(row.names(rel_otu_rare)), rel_otu_rare) %>% 
  gather(sequence_name, abun, -otu) %>%  
  left_join(map_16S[, c('sequence_name','source', 'plant', 'sampling_date', 'year', 'sampling_week','rep')], by = 'sequence_name') %>%
  filter(year == 2016, plant == 'switchgrass', source != 'soil') %>%
  filter(otu %in% sw_core_OTUs) %>%
  left_join(sw16_clusters, by='otu') %>%
  group_by(sampling_date) %>%
  mutate(sample_size = length(unique(sequence_name))) %>%
  group_by(sampling_week, sampling_date, memb_sw16, stage) %>%
  summarize(n_relabun=sum(abun)/unique(sample_size)) %>%
  filter(!is.na(memb_sw16)) %>%
  ggplot(aes(x=factor(sampling_date), y=n_relabun, color=stage, group=stage))+
  geom_line(size=2)+
  scale_color_manual(values = c('red','blue','grey'))+
  theme_classic() + theme(strip.background = element_blank(),
                          axis.text.x = element_text(angle=45, hjust = 1),
                          legend.position = 'none') +
  labs(x='Sampling date', y=NULL) +
  ylim(0,1)

abund_plot_sw17 <- data.frame(otu = as.factor(row.names(rel_otu_rare)), rel_otu_rare) %>% 
  gather(sequence_name, abun, -otu) %>%  
  left_join(map_16S[, c('sequence_name','source', 'plant', 'sampling_date', 'year', 'sampling_week','rep')], by = 'sequence_name') %>%
  left_join(sw17_clusters, by=c('otu')) %>%
  filter(otu %in% sw17_core_OTUs , year == 2017, plant == 'switchgrass', source != 'soil') %>%
  group_by(sampling_date)%>%
  mutate(sample_size = length(unique(sequence_name))) %>%
  group_by(sampling_week, sampling_date, memb_sw17, stage) %>%
  summarize(n_relabun=sum(abun)/unique(sample_size)) %>%
  filter(!is.na(memb_sw17)) %>%
  ggplot(aes(x=factor(sampling_date), y=n_relabun, color=stage, group=factor(memb_sw17)))+
  geom_line(size=2)+
  scale_color_manual(values = c('red','blue','grey'))+
  theme_classic() + theme(strip.background = element_blank(),
                          axis.text.x = element_text(angle=45, hjust = 1),
                          legend.position = 'none') +
  labs(x=NULL, y=NULL)+
  ylim(0,1)

abund_plot_mi16 <- data.frame(otu = as.factor(row.names(rel_otu_rare)), rel_otu_rare) %>% 
  gather(sequence_name, abun, -otu) %>%  
  left_join(map_16S[, c('sequence_name','source', 'plant', 'sampling_date', 'year', 'sampling_week','rep')], by = 'sequence_name') %>%
  filter(otu %in% unique(selected_otus_misc$otu) & year == 2016 & plant == 'miscanthus' & source != 'soil') %>%
  left_join(misc_clusters, by='otu') %>%
  group_by(sampling_date)%>%
  mutate(sample_size = length(unique(sequence_name))) %>%
  group_by(sampling_week, sampling_date, memb_misc, stage) %>%
  summarize(n_relabun=sum(abun)/unique(sample_size)) %>%
  ggplot(aes(x=factor(sampling_date), y=n_relabun, color=factor(stage), group=memb_misc))+
  geom_line(size=2)+
  scale_color_manual(values = c('red','blue','grey'))+
  theme_classic() + theme(strip.background = element_blank(),
                          axis.text.x = element_text(angle=45, hjust = 1),
                          legend.position = 'none') +
  labs(x=NULL, y='Relative\nabundance')+
  ylim(0,1)

grid.draw(ggarrange(FigC,FigA,FigB,
                    misc16_BC_plot,switch16_BC_plot,switch17_BC_plot, 
                    abund_plot_mi16,abund_plot_sw16,abund_plot_sw17,
                    nrow = 3))

#-------------------------------------------------------
#Heatmap - Figure S4
#-------------------------------------------------------
#' Preparing the data for plotting 
plant_otu <- otu_rare[,map_16S$source=="phyllosphere"]
plant_otu <- plant_otu[rowSums(plant_otu)>0,]
soil_otu <- otu_rare[,map_16S$source=="soil"] 
soil_otu <- soil_otu[rowSums(soil_otu)>0,]

soil_relabun <- decostand(soil_otu,method="total", MARGIN=2)
plant_relabun <- decostand(plant_otu,method="total", MARGIN=2)

soilOTU_df <- data.frame(otu=rownames(soil_relabun), soil_relabun) %>% 
  gather(sequence_name, relabun, -otu) %>%
  mutate(source='soil')
soilOTUsub_df <- soilOTU_df[soilOTU_df$otu %in% rownames(plant_relabun),]

leafOTU_df <- data.frame(otu=rownames(plant_relabun), plant_relabun) %>% 
  gather(sequence_name, relabun, -otu) %>%
  mutate(source='phyllosphere')

OTU_relabun_notime <- rbind(soilOTUsub_df, leafOTU_df) %>%
  left_join(map_16S[,c('sequence_name', 'sampling_date', 'year', 'plant')], by='sequence_name') %>%
  group_by(source, otu, year, plant) %>%
  summarise(rel_abun=sum(relabun)/length(unique(sequence_name)))

OTU_relabun <- rbind(soilOTUsub_df, leafOTU_df) %>%
  left_join(map_16S[,c('sequence_name', 'sampling_date', 'year', 'plant')], by='sequence_name') %>%
  group_by(source, otu, sampling_date,year, plant) %>%
  summarise(rel_abun=sum(relabun)/length(unique(sequence_name)))

OTU_relabun_notime$memb <- 'other'
OTU_relabun_notime$memb[OTU_relabun_notime$otu %in% unique(core_list)] <- 'core'

#' Calculating soil OTU ranks
misc_soil_otu <- otu_rare[,map_16S$plant=="miscanthus" & (map_16S$source=='soil')]
swit16_soil_otu <- otu_rare[,map_16S$plant=="switchgrass" & (map_16S$source=='soil') & (map_16S$year==2016)]
swit17_soil_otu <- otu_rare[,map_16S$plant=="switchgrass" & (map_16S$source=='soil') & (map_16S$year==2017)]
misc_soil_rank=sort(rowSums(misc_soil_otu)/sum(rowSums(misc_soil_otu)), decreasing=TRUE)
s16_soil_rank=sort(rowSums(swit16_soil_otu)/sum(rowSums(swit16_soil_otu)), decreasing=TRUE)
s17_soil_rank=sort(rowSums(swit17_soil_otu)/sum(rowSums(swit17_soil_otu)), decreasing=TRUE)

misc_leafinSoil <- misc_soil_rank[names(misc_soil_rank) %in% rownames(mis)]
s16_leafinSoil <- s16_soil_rank[names(s16_soil_rank) %in% rownames(sw)]
s17_leafinSoil <- s17_soil_rank[names(s17_soil_rank) %in% rownames(sw17)]

#Miscanthus
misc_soil_phyllo_df <- data.frame(rel_abun=misc_soil_rank, otu = names(misc_soil_rank))
misc_soil_phyllo_df$source <- 'soil'
misc_soil_phyllo_df$source[misc_soil_phyllo_df$otu %in% rownames(mis)] <- 'phyllo'
misc_soil_phyllo_df <- left_join(misc_soil_phyllo_df, OTU_relabun_notime[OTU_relabun_notime$source=='phyllosphere',], by='otu')
misc_soil_phyllo_df[is.na(misc_soil_phyllo_df)] <- 0
misc_soil_phyllo_df$otu <- factor(misc_soil_phyllo_df$otu, levels=misc_soil_phyllo_df$otu)

#' Using only the top 50 soil taxa
misc_soil_phyllo_df <- misc_soil_phyllo_df[misc_soil_phyllo_df$otu %in% rownames(mis),]
top50_miscTaxa <- as.character(unique(misc_soil_phyllo_df$otu)[1:50])

#' Creating abundance matrix for the heatmap
OTU_relabun_leaf_misc<- OTU_relabun[OTU_relabun$source == 'phyllosphere' & OTU_relabun$plant == 'miscanthus',]
rel_abun_matrix_misc<- acast(OTU_relabun_leaf_misc[OTU_relabun_leaf_misc$otu %in% top50_miscTaxa,], otu~sampling_date, value.var="rel_abun")
remaining_OTUs <- top50_miscTaxa[top50_miscTaxa %in% rownames(rel_abun_matrix_misc)]
rel_abun_matrix_misc <- rel_abun_matrix_misc[match(remaining_OTUs,rownames(rel_abun_matrix_misc)), ]
remaining_OTUs==rownames(rel_abun_matrix_misc)
scaleyellowred <- colorRampPalette(c("white", "black"), space = "rgb")(10)
sampleColor <- rep('grey', length(rownames(rel_abun_matrix_misc)))
sampleColor[rownames(rel_abun_matrix_misc) %in% misc_core_OTUs] <- 'red'
heatmap(rel_abun_matrix_misc,Colv = NA, Rowv = NA,scale="row",  
        ylab='50 top ranked soil OTUs', col = scaleyellowred,
        RowSideColors = sampleColor)

#Switchgrass 2016
s16_soil_phyllo_df <- data.frame(rel_abun=s16_soil_rank, otu = names(s16_soil_rank))
s16_soil_phyllo_df$source <- 'soil'
s16_soil_phyllo_df$source[s16_soil_phyllo_df$otu %in% rownames(sw)] <- 'phyllo'
s16_soil_phyllo_df <- left_join(s16_soil_phyllo_df, OTU_relabun_notime[OTU_relabun_notime$source=='phyllosphere',], by='otu')
s16_soil_phyllo_df[is.na(s16_soil_phyllo_df)] <- 0

#' Using only the top 50 soil taxa
s16_soil_phyllo_df <- s16_soil_phyllo_df[s16_soil_phyllo_df$otu %in% rownames(sw),]
top50_swTaxa <- as.character(unique(s16_soil_phyllo_df$otu)[1:50])

#' Creating abundance matrix for the heatmap
OTU_relabun_leaf_s16<- OTU_relabun[OTU_relabun$source == 'phyllosphere' & OTU_relabun$plant == 'switchgrass' & OTU_relabun$year == '2016',]
rel_abun_matrix_sw<- acast(OTU_relabun_leaf_s16[OTU_relabun_leaf_s16$otu %in% top50_swTaxa,], otu~sampling_date, value.var="rel_abun")
remaining_OTUs <- top50_swTaxa[top50_swTaxa %in% rownames(rel_abun_matrix_sw)]
rel_abun_matrix_sw <- rel_abun_matrix_sw[match(remaining_OTUs,rownames(rel_abun_matrix_sw)), ]
remaining_OTUs==rownames(rel_abun_matrix_sw)
scaleyellowred <- colorRampPalette(c("white", "black"), space = "rgb")(10)
sampleColor <- rep('grey', length(rownames(rel_abun_matrix_sw)))
sampleColor[rownames(rel_abun_matrix_sw) %in% sw_core_OTUs] <- 'red'
heatmap(rel_abun_matrix_sw,Colv = NA, Rowv = NA,scale="row",  
        ylab='50 top ranked soil OTUs', col = scaleyellowred,
        RowSideColors = sampleColor)

#Switchgrass 2017
s17_soil_phyllo_df <- data.frame(rel_abun=s17_soil_rank, otu = names(s17_soil_rank))
s17_soil_phyllo_df$source <- 'soil'
s17_soil_phyllo_df$source[s17_soil_phyllo_df$otu %in% rownames(sw17)] <- 'phyllo'
s17_soil_phyllo_df <- left_join(s17_soil_phyllo_df, OTU_relabun_notime[OTU_relabun_notime$source=='phyllosphere' & OTU_relabun_notime$year == 2017,], by='otu') %>% arrange(desc(rel_abun.x))
s17_soil_phyllo_df[is.na(s17_soil_phyllo_df)] <- 0

#' Using only the top 50 soil taxa found also in the phyllopshre
s17_soil_phyllo_df <- s17_soil_phyllo_df[s17_soil_phyllo_df$otu %in% rownames(sw17),]
top50_sw17_Taxa <- as.character(unique(s17_soil_phyllo_df$otu)[1:50])

#' Creating abundance matrix for the heatmap
OTU_relabun_leaf_s17<- OTU_relabun[OTU_relabun$source == 'phyllosphere' & OTU_relabun$plant == 'switchgrass' & OTU_relabun$year == '2017',]
rel_abun_matrix_sw17<- acast(OTU_relabun_leaf_s17[OTU_relabun_leaf_s17$otu %in% top50_sw17_Taxa,], otu~sampling_date, value.var="rel_abun")
remaining_OTUs <- top50_sw17_Taxa[top50_sw17_Taxa %in% rownames(rel_abun_matrix_sw17)]
rel_abun_matrix_sw17 <- rel_abun_matrix_sw17[match(remaining_OTUs,rownames(rel_abun_matrix_sw17)), ]
remaining_OTUs==rownames(rel_abun_matrix_sw17)
scaleyellowred <- colorRampPalette(c("white", "black"), space = "rgb")(10)
sampleColor <- rep('grey', length(rownames(rel_abun_matrix_sw17)))
sampleColor[rownames(rel_abun_matrix_sw17) %in% sw17_core_OTUs] <- 'red'
heatmap(rel_abun_matrix_sw17,Colv = NA, Rowv = NA,scale="row",  
        ylab='50 top ranked soil OTUs', col = scaleyellowred,
        RowSideColors = sampleColor)


###############################################################
#Suppl Fig - Dynamics of the OTUs classified as Alpha and Gamma proteobacteria
###############################################################
(proteo_plot <- data.frame(otu = as.factor(row.names(otu_rare)), otu_rare) %>% gather(sequence_name, abun, -otu) %>%
  left_join(map_16S[, c('sequence_name','rep','treatment' ,'source', 'plant', 'sampling_date', 'year','sampling_week'
  )], by = 'sequence_name') %>%
  left_join(tax_filtered, by='otu') %>%
  filter(grepl('Alphaproteobacteria|Betaproteobacteria|Deltaproteobacteria|Gammaproteobacteria', Class)) %>%
  filter(source == 'phyllosphere') %>%
  mutate(plant = factor(plant, levels = c('switchgrass', 'miscanthus'))) %>%
  mutate(class=if_else(Class=='Alphaproteobacteria','c:Alphaproteobacteria', Class)) %>%
  mutate(class=if_else(class=='Betaproteobacteria','c:Betaproteobacteria', class)) %>%
  mutate(class=if_else(class=='Gammaproteobacteria','c:Gammaproteobacteria', class)) %>%
  mutate(class=if_else(class=='Deltaproteobacteria','c:Deltaproteobacteria', class)) %>%
  group_by(plant, class, year,sampling_week) %>%
  summarise(n=sum(abun),
            n_reps=length(unique(sequence_name))) %>%
  group_by(plant, year,sampling_week) %>%
  mutate(total_reads=sum(n),
         rel_abun=n/total_reads) %>%
  ggplot(aes(x=as.factor(sampling_week), y=rel_abun, fill=class)) +
  geom_bar(color='black', stat = 'identity') +
  theme_classic() +
  facet_grid(~factor(plant, levels=c("miscanthus","switchgrass","switchgrass"), labels=c("miscanthus","switchgrass","switchgrass"))+year) +
  labs(x='Sampling week', y='Normalized relative\n abundance')+
  theme(axis.text.x = element_text(size=8),
        legend.position = 'bottom',
        legend.text=element_text(size=8),
        legend.key.size = unit(.3, "cm")) +
  guides(fill = guide_legend(ncol = 4, title=NULL)))

data.frame(otu = as.factor(row.names(otu_rare)), otu_rare) %>% gather(sequence_name, abun, -otu) %>%
  left_join(map_16S[, c('sequence_name','rep','treatment' ,'source', 'plant', 'sampling_date', 'year','sampling_week'
  )], by = 'sequence_name') %>%
  filter(source=='phyllosphere') %>%
  left_join(tax_filtered, by='otu') %>%
  filter(grepl('Alphaproteobacteria|Betaproteobacteria|Deltaproteobacteria|Gammaproteobacteria', Class)) %>%
  mutate(plant = factor(plant, levels = c('switchgrass', 'miscanthus'))) %>%
  mutate(class=if_else(Class=='Alphaproteobacteria','c:Alphaproteobacteria', Class)) %>%
  mutate(class=if_else(class=='Betaproteobacteria','c:Betaproteobacteria', class)) %>%
  mutate(class=if_else(class=='Gammaproteobacteria','c:Gammaproteobacteria', class)) %>%
  mutate(class=if_else(class=='Deltaproteobacteria','c:Deltaproteobacteria', class)) %>%
  group_by(Phylum, otu) %>%
  summarise(n=sum(abun)) %>% 
  filter(n>0) %>%
  summarise(n_sum=sum(n))


#Plotting the abundance dynamics of the selected OTUs for Switchgrass and Miscanthus together
selected_otus %>% 
  filter(source=='phyllosphere',
         otu %in% core_list
         #grepl('Alphaproteobacteria|Betaproteobacteria|Deltaproteobacteria|Gammaproteobacteria', Class)
         ) %>%
  ggplot(aes(x = as.factor(sampling_date), y = abun, color = plant, group=plant)) + 
  geom_point() +
  stat_smooth(method = "loess") +
  scale_color_manual(values=c('darkgreen','darkolivegreen3')) +
  labs(x="Sampling times", y= "Relative abundance", title='Relative abundance dynamics of the selected phyllosphere OTUs') +
  facet_wrap( ~ otu +year, scale = 'free', ncol=10) +
  theme_classic() + theme(strip.background = element_blank(),
                          axis.text.x = element_blank(),
                          legend.position = 'none')

(snippet_otu <- selected_otus %>% 
  filter(otu %in% c('OTU4','OTU2', 'OTU4223'),
         source=='phyllosphere') %>%
  mutate(members=if_else(year==2016 & plant =='miscanthus', 'mi', 'sw16'),
         members=if_else(year==2017 & plant =='switchgrass', 'sw17', members)) %>%
  left_join(map_16S[,c('sequence_name','sampling_week')], by = 'sequence_name') %>%
  ggplot(aes(x = as.factor(sampling_week), y = abun, color = plant, fill = plant, group=members)) + 
  geom_point() +
  stat_smooth(method = "loess") +
  scale_color_manual(values=c('darkgreen','darkolivegreen3')) +
  scale_fill_manual(values=c('darkgreen','darkolivegreen3')) +
  labs(x="Sampling times", y= "Relative abundance") +
  facet_wrap(year ~ as.factor(Class)+final_names, scale = 'free_y', ncol=3) +
  theme_classic() + theme(strip.background = element_blank(),
                          axis.text.x = element_text(size=8),
                          legend.position = 'none',
                          strip.text.x = element_blank()))

setEPS()
postscript('~/Desktop/figure/Figure5.eps', width = 6,height = 7, paper = 'special')
grid.draw(ggarrange(proteo_plot,
                    snippet_otu,
                    heights = 1:2.5))
dev.off()

selected_otus %>% 
  filter(source=='phyllosphere',
         Genus == 'g:Hymenobacter') %>%
  mutate(members=if_else(year==2016 & plant =='miscanthus', 'mi', 'sw16'),
         members=if_else(year==2017 & plant =='switchgrass', 'sw17', members)) %>%
  left_join(map_16S[,c('sequence_name','sampling_week')], by = 'sequence_name') %>%
  ggplot(aes(x = as.factor(sampling_week), y = abun, color = plant, fill = plant, group=members)) + 
  geom_point() +
  stat_smooth(method = "loess") +
  scale_color_manual(values=c('darkgreen','darkolivegreen3')) +
  scale_fill_manual(values=c('darkgreen','darkolivegreen3')) +
  labs(x="Sampling times", y= "Relative abundance") +
  facet_wrap(year ~ otu, scale = 'free_y', ncol=6) +
  theme_classic() + theme(strip.background = element_blank(),
                          axis.text.x = element_text(size=8),
                          legend.position = 'none',
                          strip.text.x = element_blank())


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

map.plant <- map_16S[map_16S$source=="phyllosphere",]
map.plant.2016 <- map.plant[map.plant$year==2016,]
map.plant.2017 <- map.plant[map.plant$year==2017,]
map.plant.switch <- map.plant.2016[map.plant.2016$plant=="switchgrass",] 
map.plant.misc <- map.plant.2016[map.plant.2016$plant=="miscanthus",]


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

misc.map.2016 <- map.plant.2016[map.plant.2016$plant=="miscanthus",]
switch.map <- map.plant[map.plant$plant=="switchgrass",]
switch.map.2016 <- map.plant[map.plant$plant=="switchgrass"&map.plant$year==2016,]
switch.map.2017 <- map.plant[map.plant$plant=="switchgrass"&map.plant$year==2017,]

cor.test(map.div.swg$time_numeric, map.div.swg$Richness)
cor.test(map.div.swg$time_numeric, map.div.swg$Shannon)
cor.test(map.div.swg$time_numeric, map.div.swg$Pielou)

cor.test(map.div.mis$time_numeric, map.div.mis$Richness)
cor.test(map.div.mis$time_numeric, map.div.mis$Shannon)
cor.test(map.div.mis$time_numeric, map.div.mis$Pielou)


plant_rare_otu <- otu_rare[,map_16S$source=="phyllosphere"]
plant_rare_otu.2016 <- otu_rare[,map_16S$source=="phyllosphere"&map_16S$year==2016]
plant_rare_otu.2017 <- otu_rare[,map_16S$source=="phyllosphere"&map_16S$year==2017]

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


# Look at species accumulation
misc_accumulation <- rep(1, length(misc_unique_times))
z <- NULL
for( i in 1:length(unique(misc.map.2016$time_numeric))){
  x <- misc_rare_otu.2016[,misc.map.2016$time_numeric==misc_unique_times[i]]
  y <- matrix(x, ncol=sum(misc.map.2016$time_numeric==misc_unique_times[i]))
  row.names(y) <- row.names(misc_rare_otu.2016)
  z <- c(z, row.names(misc_rare_otu.2016[rowSums(y)!=0,]))
  misc_accumulation[i] <- length(unique(z))
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


NewDates_Switch2017 <- c(as.Date("2016-05-15"), as.Date("2016-06-05"), as.Date("2016-06-26"), as.Date("2016-07-17"),as.Date("2016-08-07"), as.Date("2016-08-28"), as.Date("2016-09-18"))

Spec_accum_total <- NULL
Spec_accum_total$Species_Accumulation <- c(misc_accumulation, switch_accumulation.2016, switch_accumulation.2017)
Spec_accum_total$Date <- c(unique(misc.map.2016$sampling_Rdate)[order(unique(misc.map.2016$sampling_Rdate))], unique(switch.map.2016$sampling_Rdate)[order(unique(switch.map.2016$sampling_Rdate))], NewDates_Switch2017)
 Spec_accum_total <- as.data.frame(Spec_accum_total)
 Spec_accum_total$Plant <- c(rep("Miscanthus", length(misc_accumulation)), rep("Switchgrass", length(switch_accumulation.2016)), rep("Switchgrass", length(switch_accumulation.2017)))



# Code of Figure 2, Phyllosphere richness and accumulation through the season
switch.Richness.2016 <- map.alpha.plant.2016[map.alpha.plant.2016$plant=="switchgrass"&map.alpha.plant.2016$variable=="Richness",]
fit_model <- aov(value~factor(sampling_Rdate),switch.Richness.2016)
summary(fit_model)
TukeyHSD(fit_model)

misc.Richness <- map.alpha.plant[map.alpha.plant$plant=="miscanthus"&map.alpha.plant$variable=="Richness",]
fit_model <- aov(value~factor(sampling_Rdate),misc.Richness)
summary(fit_model)
TukeyHSD(fit_model)

switch.Richness.2017 <- map.alpha.plant.2017[map.alpha.plant.2017$plant=="switchgrass"&map.alpha.plant.2017$variable=="Richness",]


Richness.df <- rbind(switch.Richness.2016, switch.Richness.2017, misc.Richness)
Richness.df$Factor <- factor(c(rep("Switchgrass 2016", nrow(switch.Richness.2016)), rep("Switchgrass 2017", nrow(switch.Richness.2017)), rep("Miscanthus 2016", nrow(misc.Richness))), levels = c("Miscanthus 2016", "Switchgrass 2016", "Switchgrass 2017"))
Richness.df$FakeDate <- as.Date(gsub(x=Richness.df$sampling_Rdate, pattern = "2017", replacement = "2016" ))

  

### Spec Accumulation Curves
Spec_accum_total$Year <- c(rep("2016",length(misc_accumulation)), rep("2016",length(switch_accumulation.2016)), rep("2017",length(switch_accumulation.2017)))
Spec_accum_total
Spec_accum_total$PlantYear <- paste(Spec_accum_total$Plant, Spec_accum_total$Year)

Fig1B <- ggplot(Spec_accum_total, aes(x=Date, y=Species_Accumulation)) + 
  geom_point(aes(color=Plant, shape=PlantYear), size=2) +
  geom_line(aes(color=Plant, linetype=Year))+
  scale_shape_manual(values=c(15,16,1))+
  scale_linetype_manual(values = c(1,2))+
  theme_bw()+
  labs(y="Total Observed Taxa") + 
  scale_color_manual(values=c("darkgreen","darkolivegreen3"))  + 
  scale_x_date(date_breaks = "1 month", labels=date_format("%b"))+
  guides(shape=FALSE)

ggsave("Figures/Figure1B_SpeciesAccumulation.eps", Fig1B, device = "eps", width=6, height=4, units = "in")

# Plot Miscanthus Evenness Through Season
switch.Pielou.2016 <- map.alpha.plant.2016[map.alpha.plant.2016$plant=="switchgrass"&map.alpha.plant.2016$variable=="Pielou",]

switch.Pielou.2017 <- map.alpha.plant.2017[map.alpha.plant.2017$plant=="switchgrass"&map.alpha.plant.2017$variable=="Pielou",]

misc.Pielou.2016 <- map.alpha.plant.2016[map.alpha.plant.2016$plant=="miscanthus"&map.alpha.plant.2016$variable=="Pielou",]

t.test(map.alpha.plant.2016[map.alpha.plant.2016$variable=="Richness"&map.alpha.plant.2016$treatment=="standard fertilization",]$value, map.alpha.plant.2016[map.alpha.plant.2016$variable=="Richness"&map.alpha.plant.2016$treatment=="nitrogen free",]$value)


################################
#### Beta Diversity Analyses####
################################
Date_to_Week <- data.frame(sampling_date=unique(map_16S$sampling_date)[order(unique(map_16S$sampling_date))], Week=c(1:10,1:8))

map_16S <- left_join(map_16S, Date_to_Week, by="sampling_date")
map.2016 <- map_16S[map_16S$Year==2016,]
map.2017 <- map_16S[map_16S$Year==2017,]

map.soil <- map_16S[map_16S$source=="soil",]
map.soil.2016 <- map.soil[map.soil$Year=="2016",]
map.soil.2017 <- map.soil[map.soil$Year=="2017",]

otu.2016 <- otu_rare[,map_16S$Year==2016]
otu.2017 <- otu_rare[,map_16S$Year==2017]


otu.plant <- otu_rare[,map_16S$source=="phyllosphere"]

otu.plant.2016 <- otu_rare[,map_16S$source=="phyllosphere"&map_16S$Year==2016]
otu.plant.2017 <- otu_rare[,map_16S$source=="phyllosphere"&map_16S$Year==2017]

otu.soil.2016 <- otu_soil_rare[,map.soil$Year=="2016"]
otu.soil.2017 <- otu_soil_rare[,map.soil$Year=="2017"]

dist.otu <- vegdist(t(otu_rare), method="bray")
dist.otu.2016 <- vegdist(t(otu.2016), method="bray")
dist.otu.2017 <- vegdist(t(otu.2017), method="bray")

dist.plant <- vegdist(t(otu.plant), method="bray")
dist.plant.2016 <- vegdist(t(otu.plant.2016), method="bray")
dist.plant.2017 <- vegdist(t(otu.plant.2017), method="bray")

dist.soil.2016 <- vegdist(t(otu.soil.2016), method="bray")
dist.soil.2017 <- vegdist(t(otu.soil.2017), method="bray")

pcoa.plant.2016 <- cmdscale(dist.plant.2016, eig=TRUE)
pcoa.plant.2017 <- cmdscale(dist.plant.2017, eig=TRUE)
pcoa.plant <- cmdscale(dist.plant, eig=TRUE)


adonis(dist.plant.2016~map.plant.2016$time_numeric)
adonis(dist.plant.2017~map.plant.2017$time_numeric)
adonis(dist.soil.2016~map.soil.2016$time_numeric)
adonis(dist.soil.2017~map.soil.2017$time_numeric)

adonis(dist.plant.2016~map.plant.2016$treatment)
adonis(dist.plant.2017~map.plant.2017$treatment)

adonis(dist.soil.2016~map.soil.2016$treatment)
adonis(dist.soil.2017~map.soil.2017$treatment)


adonis(dist.plant.2016~map.plant.2016$plant)
adonis(dist.soil.2016~map.soil.2016$plant)

adonis(dist.otu.2016~map.2016$source)
adonis(dist.otu.2017~map.2017$source)


#######################################
### Setting up collapsed plant pcoa ###
#######################################

### Start by getting the dates for plant smaples
plant.dates <- unique(map.plant$sampling_date)[order(unique(map.plant$sampling_date))]
sg.dates <- unique(map.plant.switch$sampling_date)[order(unique(map.plant.switch$sampling_date))]
sg.dates <- c(sg.dates,unique(map.plant.2017$sampling_date)[order(unique(map.plant.2017$sampling_date))] )
ms.dates <- unique(map.plant.misc$sampling_date)[order(unique(map.plant.misc$sampling_date))]

### Set up the points from the "entire" plant PCoA

plant.points <- pcoa.plant$points
sg.points <- plant.points[map.plant$plant=="switchgrass",]
ms.points <- plant.points[map.plant$plant=="miscanthus",]

### determine the average location of each crop on the PCoA at each timepoint
sg.points.collapsed <- NULL
sg.samples.per.date <- NULL
for (i in 1:length(sg.dates)){
  x <- plant.points[map.plant$sampling_date==sg.dates[i]&map.plant$plant=="switchgrass",]
  sg.samples.per.date <- c(sg.samples.per.date, nrow(x))
  sg.points.collapsed <- rbind(sg.points.collapsed, c(colSums(x), sd(x[,1]), sd(x[,2])))
}
sg.points.collapsed <- as.data.frame(sg.points.collapsed)
colnames(sg.points.collapsed) <- c("Axis1", "Axis2", "sd_axis1", "sd_axis2")
sg.points.collapsed$sampling_date <- sg.dates
sg.points.collapsed[,1:2] <- sg.points.collapsed[,1:2]/sg.samples.per.date

ms.points.collapsed <- NULL
ms.samples.per.date <- NULL
for(i in 1:length(ms.dates)){
  x <- ms.points[map.plant.misc$sampling_date==ms.dates[i],]
  ms.samples.per.date <- c(ms.samples.per.date, nrow(x))
  ms.points.collapsed <- rbind(ms.points.collapsed, c(colSums(matrix(x, ncol=2)), sd(matrix(x,ncol=2)[,1]), sd(matrix(x, ncol=2)[,2])))
}

ms.points.collapsed <- as.data.frame(ms.points.collapsed)
colnames(ms.points.collapsed) <- c("Axis1", "Axis2", "sd_axis1", "sd_axis2")
ms.points.collapsed$sampling_date <- ms.dates
ms.points.collapsed[,1:2] <- ms.points.collapsed[,1:2]/ms.samples.per.date

plant.points.collapsed <- rbind(sg.points.collapsed, ms.points.collapsed)
plant.points.collapsed <- left_join(plant.points.collapsed, Date_to_Week, by="sampling_date")
plant.points.collapsed$plant <- c(rep("Switchgrass", 13), rep("Miscanthus", 9))
plant.points.collapsed$Week <- factor(x = plant.points.collapsed$Week, levels=c(1,2,3,4,5,6,7,8,9,10))
plant.points.collapsed$WeekNumeric <- as.numeric(plant.points.collapsed$Week)
plant.points.collapsed$Year <- c(rep(2016,6), rep(2017,7), rep(2016,9))
plant.points.collapsed$Year <- factor(plant.points.collapsed$Year, levels = c(2016,2017))
plant.points.collapsed$plant <- factor(plant.points.collapsed$plant, levels = c("Switchgrass", "Miscanthus"))
plant.points.collapsed <- plant.points.collapsed[order(plant.points.collapsed$sampling_date),]

### Set up the weather maps for Weather Environmental fits
switch.weather <- map.plant.switch[,c("precipitation", "Air_temp_mean", "air_temp_max", "Air_Temp_Min", "Air_Pressure", "RH", "AH", "Wind_Speed_Mean", "Solar_Radiation", "PAR", "soil_temp_5_cm_bare_avg", "sampling_date", "time_numeric")]
switch.weather <- unique(switch.weather)

switch.weather.2017 <- map.plant.2017[,c("precipitation", "Air_temp_mean", "air_temp_max", "Air_Temp_Min", "Air_Pressure", "RH", "AH", "Wind_Speed_Mean", "Solar_Radiation", "PAR", "soil_temp_5_cm_bare_avg", "sampling_date", "time_numeric")]
switch.weather.2017 <- unique(switch.weather.2017)

misc.weather <- map.plant.misc[,c("precipitation", "Air_temp_mean", "air_temp_max", "Air_Temp_Min", "Air_Pressure", "RH", "AH", "Wind_Speed_Mean", "Solar_Radiation", "PAR", "soil_temp_5_cm_bare_avg", "sampling_date", "time_numeric")]
misc.weather <- unique(misc.weather)


scale_arrow <- function(arrows, data, at = c(0, 0), fill = 0.75) {
  u <- c(range(data[,1], range(data[,2])))
  u <- u - rep(at, each = 2)
  r <- c(range(arrows[, 1], na.rm = TRUE), range(arrows[, 2], na.rm = TRUE))
  rev <- sign(diff(u))[-2]
  
  if (rev[1] < 0) {
    u[1:2] <- u[2:1]
  }
  if (rev[2] < 0) {
    u[3:4] <- u[4:3]
  }
  u <- u/r
  u <- u[is.finite(u) & u > 0]
  invisible(fill * min(u))
}


arrow_scaling <- scale_arrow(data = scores(pcoa.plant, display = "sites"), arrows = scores(envfit.phyllo.sc, display = "vectors"))

arrow_scaling <- scale_arrow(data = scores(pcoa.plant, display = "sites"), arrows = rbind(scores(envfit.phyllo.sc, display="vectors"), scores(envfit.phyllo.lc, display="vectors")))



### Set up EnvFit for ggplot (remembering to use vegan:::ordiArrowMul to adjust arrow sizes)
phyllo.weather <- rbind(switch.weather, switch.weather.2017, misc.weather)
phyllo.weather <- phyllo.weather[order(phyllo.weather$sampling_date),]
phyllo.weather$sampling_date == plant.points.collapsed$sampling_date
phyllo.weather$Week <- plant.points.collapsed$WeekNumeric

envfit.phyllo.weather <- envfit(plant.points.collapsed[,1:2], phyllo.weather[,c(1:11,14)])
envfit.phyllo.weather.df<-as.data.frame(scores(envfit.phyllo.weather, display = "vectors"))
colnames(envfit.phyllo.weather.df) <- c("Dim1", "Dim2")




phyllo.leaf.chemistry <- map.plant[,c("LDMC_mg_per_g", "nitrogen_percent", "carbon_percent", "carbon_per_nitrogen", "height_mean_cm")]
envfit.phyllo.lc <- envfit(pcoa.plant$points, phyllo.leaf.chemistry)
envfit.phyllo.lc.df<-as.data.frame(scores(envfit.phyllo.lc, display = "vectors"))


phylloi.soil.chemsitry <- map.plant[,c("pH", "P_ppm", "K_ppm", "Ca_ppm", "Mg_ppm", "organic_matter", "NO3N_ppm", "NH4_ppm", "soil_moisture_percent", "soil_temp_10cm")]
envfit.phyllo.sc <- envfit(pcoa.plant$points, phylloi.soil.chemsitry)
envfit.phyllo.sc.df<-as.data.frame(scores(envfit.phyllo.sc, display = "vectors"))

envfit.phyllo.total <- rbind(envfit.phyllo.weather.df, envfit.phyllo.lc.df, envfit.phyllo.sc.df)

envfit.phyllo.total$r2 <- c(envfit.phyllo.weather$vectors$r, envfit.phyllo.lc$vectors$r, envfit.phyllo.sc$vectors$r)
envfit.phyllo.total$pval <- c(envfit.phyllo.weather$vectors$pvals, envfit.phyllo.lc$vectors$pvals, envfit.phyllo.sc$vectors$pvals)

arrow_scaling <- scale_arrow(arrows=envfit.phyllo.total[,1:2], data=plant.points.collapsed)

envfit.phyllo.total$Axis1 <- envfit.phyllo.total$Dim1 * arrow_scaling
envfit.phyllo.total$Axis2 <- envfit.phyllo.total$Dim2 * arrow_scaling
envfit.phyllo.total$Variable <- row.names(envfit.phyllo.total)

### Subset to p <0.05 & Rsquared >0.4
envfit.phyllo.sub <- subset(envfit.phyllo.total, envfit.phyllo.total$r2 >0.3999) 
envfit.phyllo.sub <- subset (envfit.phyllo.sub , envfit.phyllo.sub$pval<0.05)


write.table(file = "Figures/TableS2_EnvFit.txt", x=envfit.phyllo.total[,1:4], sep="\t", quote=FALSE)



### Set up some of the plotting specifics
Point_Sizes <- seq(from=2, to=6, length.out = 10)
Ax1.plant <- pcoa.plant$eig[1]/sum(pcoa.plant$eig)
Ax2.plant <- pcoa.plant$eig[2]/sum(pcoa.plant$eig)

### Plot Plant PCoA
Fig2B <- ggplot(plant.points.collapsed, aes(x=Axis1, y=Axis2))+
  geom_point(aes(size=Week, color=plant, shape=Year)) +
  scale_color_manual(values=c("darkolivegreen3","darkgreen"))+
  scale_shape_manual(values = c(19,1))+
  scale_size_manual(values=Point_Sizes)+
  coord_fixed()+
  geom_segment(data=plant.points.collapsed, aes(x=Axis1,xend=Axis1+sd_axis1,y=Axis2,yend=Axis2,color=plant ))+
  geom_segment(data=plant.points.collapsed, aes(x=Axis1,xend=Axis1-sd_axis1,y=Axis2,yend=Axis2,color=plant))+
  geom_segment(data=plant.points.collapsed, aes(x=Axis1,xend=Axis1,y=Axis2,yend=Axis2+sd_axis2,color=plant))+
  geom_segment(data=plant.points.collapsed, aes(x=Axis1,xend=Axis1,y=Axis2,yend=Axis2-sd_axis2,color=plant))+
  xlab(label = paste(round(Ax1.plant,digits = 3)*100, "% Var. Explained", sep = ""))+
  ylab(label= paste(round(Ax2.plant,digits = 3)*100, "% Var. Explained", sep = ""))+
  geom_segment(data = envfit.phyllo.sub,
   aes(x = 0, xend = Axis1, y = 0, yend = Axis2),
   arrow = arrow(length = unit(0.25, "cm")), colour = "grey")+
  geom_text(data=envfit.phyllo.sub, aes(label=Variable))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.position = "none")

ggsave(filename = "Figures/Figure2B_PhyllospherePCoA.eps", Figure2B, device = "eps", height = 6, width = 6, units = "in")

### Make Whole PCoA
pcoa.whole <- cmdscale(dist.otu, eig=TRUE)

### Set up dataframe for plotting
points.whole <- pcoa.whole$points
points.whole <- as.data.frame(points.whole)
colnames(points.whole) <- c("Axis1", "Axis2")
points.whole$plant <- factor(map_16S$plant, levels=c("switchgrass", "miscanthus"))
points.whole$source <- factor(map_16S$source, levels=c("soil", "phyllosphere"))
points.whole$Year <- factor(map_16S$Year, levels=c(2016,2017))
points.whole$sampling_date <- map_16S$sampling_date
points.whole <- left_join(points.whole, Date_to_Week, by="sampling_date")
points.whole$Week <- factor(points.whole$Week, levels=c(1:10))
points.whole$PlantYear <- paste(points.whole$plant,points.whole$Year, sep="")
points.whole$PlantYear <- factor(points.whole$PlantYear, levels=c("miscanthus2016", "switchgrass2016", "switchgrass2017"))

points.whole$SourceFert <- paste(map_16S$source, map_16S$treatment)
points.whole$SourceFert[grepl(pattern = "phyllosphere", points.whole$SourceFert)] <- "phyllosphere"
points.whole$SourceFert <- factor(points.whole$SourceFert, levels = c("soil standard fertilization", "soil nitrogen free", "phyllosphere"))
points.whole$SampleType <- points.whole$SourceFert


### Determine % variation explained on each axis
Ax1.whole <- pcoa.whole$eig[1]/sum(pcoa.whole$eig)
Ax2.whole <- pcoa.whole$eig[2]/sum(pcoa.whole$eig)

### ggplot Whole PCoA
Fig2A <- ggplot(points.whole, aes(x=Axis1, y=Axis2))+
  coord_fixed()+
  geom_point(aes(shape=SampleType, color=plant, size=Week, fill=PlantYear))+
  scale_shape_manual(values = c(22,23,21))+
  scale_size_manual(values=Point_Sizes)+
  scale_color_manual(values=c("darkolivegreen3","darkgreen"))+
  scale_fill_manual(values=c("darkgreen", "darkolivegreen3","NA"))+
  xlab(label = paste(round(Ax1.whole, 3)*100, "% Var. Explained", sep=""))+
  ylab(label= paste(round(Ax2.whole, 3)*100, "% Var. Explained", sep=""))+
  theme(legend.position = "none")+
  guides(fill=guide_legend(override.aes=list(colour=c(Miscanthus2016="darkgreen",Switchgrass2016="darkolivegreen3", Switchgrass2017="white"))), color=FALSE)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggsave("Figures/Figure2A_WholePCoA_NoLegend.eps", Fig2A, width = 6, height = 6, device = "eps", units="in")


Fig2A_leg <- ggplot(points.whole, aes(x=Axis1, y=Axis2))+
  coord_fixed()+
  geom_point(aes(shape=SampleType, color=plant, size=Week, fill=PlantYear))+
  scale_shape_manual(values = c(22,23,21))+
  scale_size_manual(values=Point_Sizes)+
  scale_color_manual(values=c("darkolivegreen3","darkgreen"))+
  scale_fill_manual(values=c("darkgreen", "darkolivegreen3","NA"))+
  xlab(label = paste(round(Ax1.whole, 3)*100, "% Var. Explained", sep=""))+
  ylab(label= paste(round(Ax2.whole, 3)*100, "% Var. Explained", sep=""))+
  theme(legend.position = "right", legend.box = "horizontal")+
  guides(fill=guide_legend(override.aes=list(shape=21, color=c(miscanthus2016="darkgreen",switchgrass2016="darkolivegreen3", switchgrass2017="darkolivegreen3"),fill=c(miscanthus2016="darkgreen",switchgrass2016="darkolivegreen3", switchgrass2017="white"), size=5)), color=FALSE, shape=guide_legend(override.aes = list(size=5)))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

Fig2_Legend <- get_legend(Fig2A_leg)

setEPS()
postscript("Figures/Figure2.eps", width=10, height=10,pointsize=10, paper="special")
plot_grid(plot_grid(Fig2B, Fig2A, align = "h"), Fig2_Legend, ncol=1)
dev.off()



library(cowplot)

GenSciLegend <- ggplot(points.whole, aes(x=Axis1, y=Axis2))+
  coord_fixed()+
  geom_point(aes(shape=SampleType, color=plant, size=Week, fill=PlantYear))+
  scale_shape_manual(values = c(22,23,21))+
  scale_size_manual(values=Point_Sizes)+
  scale_color_manual(values=c("darkolivegreen3","darkgreen"))+
  scale_fill_manual(values=c("darkgreen", "darkolivegreen3","NA"))+
  xlab(label = paste(round(Ax1.whole, 3)*100, "% Var. Explained", sep=""))+
  ylab(label= paste(round(Ax2.whole, 3)*100, "% Var. Explained", sep=""))+
  theme(legend.position = "right", legend.box = "horizontal")+
  guides(fill=guide_legend(override.aes=list(shape=21, colour=c(Miscanthus2016="darkgreen",Switchgrass2016="darkolivegreen3", Switchgrass2017="white"), size=5)), color=FALSE, shape=guide_legend(override.aes = list(size=5)))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

legend <- get_legend(GenSciLegend)


setEPS()
postscript("~/Desktop/GenScipPCoAs.eps", width=10, height=10,pointsize=10, paper="special",)
plot_grid(Figure2B, FigS2B, Fig2A, legend)
dev.off()

###################
### For PROTEST ###
###################


switch.otu.2016 <- otu_rare[,map_16S$source=="phyllosphere"&map_16S$plant=="switchgrass"&map_16S$year==2016&map_16S$month!=10]

switch.otu.2017 <- otu_rare[,map_16S$source=="phyllosphere"&map_16S$plant=="switchgrass"&map_16S$year==2017]

switch.map.2016.p <- map_16S[map_16S$source=="phyllosphere"&map_16S$plant=="switchgrass"&map_16S$year==2016&map_16S$month!=10,]
switch.map.2017.p <- map_16S[map_16S$source=="phyllosphere"&map_16S$plant=="switchgrass"&map_16S$year==2017,]



switch.collapsed.otu.2016 <- data.frame(Time3=rep(NA, nrow(switch.otu.2016)), Time4=rep(NA, nrow(switch.otu.2016)), Time5=rep(NA, nrow(switch.otu.2016)), Time6=rep(NA, nrow(switch.otu.2016)), Time7=rep(NA, nrow(switch.otu.2016)))

switch.collapsed.otu.2017 <- data.frame(Time3=rep(NA, nrow(switch.otu.2017)), Time4=rep(NA, nrow(switch.otu.2017)), Time5=rep(NA, nrow(switch.otu.2017)), Time6=rep(NA, nrow(switch.otu.2017)), Time7=rep(NA, nrow(switch.otu.2017)))

a <- c(4:8)

for (i in 1:5){
  temp_2016 <- switch.otu.2016[,switch.map.2016.p$Week==a[i]]
  temp_2017 <- switch.otu.2017[, switch.map.2017.p$Week==a[i]]
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



set.seed(13)
c <- as.data.frame(table(map.plant$sampling_date))
set.seed(13)
otu_min <- t(rrarefy(t(otu_filtered), min(colSums(otu_filtered))))
set.seed(13)
otu_500 <- t(rrarefy(t(otu_filtered), 500))
otu_500 <- otu_500[,colSums(otu_500)>499]
set.seed(13)
otu_1000 <- t(rrarefy(t(otu_filtered), 1000))  
otu_1000 <- otu_1000[,colSums(otu_1000)>999]
set.seed(13)
otu_5k <- t(rrarefy(t(otu_filtered), 5000))
otu_5k <- otu_5k[,colSums(otu_5k)>4999]
set.seed(13)            
otu_10K <- t(rrarefy(t(otu_filtered), 10000))
otu_10K <- otu_10K[,colSums(otu_10K)>9999]

otu_min_10K <- otu_min[,colnames(otu_min)%in%colnames(otu_10K)]
otu_min_10k.dist <- vegdist(t(otu_min_10K), method="bray")
otu_10K.dist <- vegdist(t(otu_10K), method="bray")
mantel(otu_min_10k.dist, otu_10K.dist)

otu_500_10K <- otu_500[,colnames(otu_500)%in%colnames(otu_10K)]
otu_500_10k.dist <- vegdist(t(otu_500_10K), method="bray")
mantel(otu_500_10k.dist, otu_10K.dist)

otu_1000_10K <- otu_1000[,colnames(otu_1000)%in%colnames(otu_10K)]
otu_1000_10k.dist <- vegdist(t(otu_1000_10K), method="bray")
mantel(otu_1000_10k.dist, otu_10K.dist)

otu_5k_10K <- otu_5k[,colnames(otu_5k)%in%colnames(otu_10K)]
otu_5k_10k.dist <- vegdist(t(otu_5k_10K), method="bray")
mantel(otu_5k_10k.dist, otu_10K.dist)

otu_5k.dist <- vegdist(t(otu_5k), method = "bray")

otu_min_5k <- otu_min[,colnames(otu_min)%in%colnames(otu_5k)]
otu_min_5k.dist <- vegdist(t(otu_min_5k), method = "bray")
mantel(otu_min_5k.dist, otu_5k.dist)

otu_500_5k <- otu_500[,colnames(otu_500)%in%colnames(otu_5k)]
otu_500_5k.dist <- vegdist(t(otu_500_5k), method = "bray")
mantel(otu_500_5k.dist, otu_5k.dist)

otu_1000_5k <- otu_1000[,colnames(otu_1000)%in%colnames(otu_5k)]
otu_1000_5k.dist <- vegdist(t(otu_1000_5k), method = "bray")
mantel(otu_1000_5k.dist, otu_5k.dist)


otu_1000.dist <- vegdist(t(otu_1000), method = "bray")

otu_min_1000 <- otu_min[,colnames(otu_min)%in%colnames(otu_1000)]
otu_min_1000.dist <- vegdist(t(otu_min_1000), method = "bray")
mantel(otu_min_1000.dist, otu_1000.dist)

otu_500_1000 <- otu_500[,colnames(otu_500)%in%colnames(otu_1000)]
otu_500_1000.dist <- vegdist(t(otu_500_1000), method = "bray")
mantel(otu_500_1000.dist, otu_1000.dist)


otu_500.dist <- vegdist(t(otu_500), method="bray")

otu_min_500 <- otu_min[,colnames(otu_min)%in%colnames(otu_500)]
otu_min_500.dist <- vegdist(t(otu_min_500), method = "bray")
mantel(otu_min_500.dist, otu_500.dist)




adonis(otu_1k_Clean_phyllo_2016_noMay.bc~map_1k_clean_phyllo_2016_noMay$time_numeric)
adonis(otu_1k_Clean_phyllo_2016_noMay.bc~map_1k_clean_phyllo_2016_noMay$Week)
adonis(otu_1k_Clean_phyllo_2016_noMay.bc~map_1k_clean_phyllo_2016_noMay$plant)



anosim(otu_1k_Clean_phyllo_2016_noMay.bc, map_1k_clean_phyllo_2016_noMay$Week)
anosim(otu_1k_Clean_phyllo_2016_noMay.bc, map_1k_clean_phyllo_2016_noMay$plant)



otu_1k_10k_Clean <- otu_1000[,colnames(otu_1000)%in%colnames(otu_10K)]

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



otu_1k_10k_Clean <- otu_1000[,colnames(otu_1000)%in%colnames(otu_10K_Clean)]
colnames(otu_1k_10k_Clean)==colnames(otu_10K_Clean)

bc.1k_10K_Clean <- vegdist(t(otu_1k_10k_Clean), method = "bray")
bc.10K_Clean <- vegdist(t(otu_10K_Clean), method="bray")


mantel(bc.1k_10K_Clean, bc.10K_Clean)

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

#################
### Soil PCoA ###
#################

soil.context <- map.soil[,c("precipitation", "Air_temp_mean", "air_temp_max", "Air_Temp_Min", "Air_Pressure", "RH", "AH", "Wind_Speed_Mean", "Solar_Radiation", "PAR", "soil_temp_5_cm_bare_avg", "sampling_date", "Week.x","pH", "P_ppm", "K_ppm", "Ca_ppm", "Mg_ppm", "organic_matter", "NO3N_ppm", "NH4_ppm", "soil_moisture_percent", "soil_temp_10cm")]

dist.soil <- vegdist(t(otu_soil_rare), method="bray")
pcoa.soil <- cmdscale(dist.soil, eig = TRUE)

points.soil <- as.data.frame(pcoa.soil$points)
points.soil$Week <- map.soil$Week.x
points.soil$Plant <- map.soil$plant
points.soil$SourceFert <- paste(map.soil$source, map.soil$treatment)
points.soil$SourceFert <- factor(points.soil$SourceFert, levels=c("soil standard fertilization", "soil nitrogen free", "phyllosphere"))
points.soil$PlantYear <- factor(paste(map.soil$plant, map.soil$Year, sep=""), levels=c("miscanthus2016", "switchgrass2016", "switchgrass2017"))

points.soil$Week <- as.factor(points.soil$Week)


soil.envfit <- envfit(pcoa.soil$points, soil.context)

soil.ef.plot <- as.data.frame(soil.envfit$vectors$arrows)
soil.ef.plot$r2 <- soil.envfit$vectors$r
soil.ef.plot$pval <- soil.envfit$vectors$pvals
soil.ef.plot
soil.ef.plot$Variable <- row.names(soil.ef.plot)

arrow_scaling_soil <- scale_arrow(arrows=soil.ef.plot[,1:2], data=pcoa.soil$points)

soil.ef.plot$Axis1 <- soil.ef.plot$Dim1* arrow_scaling_soil
soil.ef.plot$Axis2 <- soil.ef.plot$Dim2* arrow_scaling_soil

soil.ef.plot.sub <- soil.ef.plot[soil.ef.plot$pval<0.05,]
soil.ef.plot.sub <- soil.ef.plot.sub[soil.ef.plot.sub$r2>.25,]

soil.ef.plot.sub$Variable <- gsub(pattern = "Week.x", replacement = "Week", soil.ef.plot.sub$Variable)


Ax1.soil <- pcoa.soil$eig[1]/sum(pcoa.soil$eig)
Ax2.soil <- pcoa.soil$eig[2]/sum(pcoa.soil$eig)
colnames(points.soil)<- c("Axis1", "Axis2", "Week", "Plant", "SourceFert", "PlantYear", "SampleType")

levels(points.soil$SampleType)

### ggplot Whole PCoA
FigureS2C <- ggplot(points.soil, aes(x=Axis1, y=Axis2))+
  coord_fixed()+
  geom_point(aes(shape=SourceFert, color=Plant, size=Week, fill=PlantYear))+
  scale_shape_manual(values = c(22,23,21))+
  scale_size_manual(values=Point_Sizes)+
  scale_color_manual(values=c("darkgreen", "darkolivegreen3"))+
  scale_fill_manual(values=c("darkgreen", "darkolivegreen3","NA"))+
  xlab(label = paste(round(Ax1.soil, 3)*100, "% Var. Explained", sep=""))+
  ylab(label= paste(round(Ax2.soil, 3)*100, "% Var. Explained", sep=""))+
  geom_segment(data = soil.ef.plot.sub,
               aes(x = 0, xend = Axis1, y = 0, yend = Axis2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "grey")+
  geom_text(data=soil.ef.plot.sub, aes(label=Variable))+
  theme(legend.position = "none")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))





#######################
### Analysis at 141 ###
#######################

otu_141.PA <- 1*(otu_141>0)

s <- specnumber(otu_141,MARGIN=2)
h <- diversity(t(otu_141), index = "shannon")
pielou=h/log(s)

map.141 <- map_small[map_small$sequence_name%in%colnames(otu_141),]

map.141$sampling_date <- gsub(pattern = " EDT", replacement = "", map.141$sampling_date)
map.141$sampling_Rdate <- as.Date(map.141$sampling_date)

otu_141.phyllo.PA <- otu_141.PA[,map.141$source=="phyllosphere"]

### Settin up Contextual Data Maps
map.141.2017 <- map.141[map.141$year=="2017",]
map.141.2016 <- map.141[map.141$year=="2016",]

map.141.plant <- map.141[map.141$source=="phyllosphere",]
map.141.plant.2016 <- map.141.plant[map.141.plant$year==2016,]
map.141.plant.2017 <- map.141.plant[map.141.plant$year==2017,]
map.141.plant.switch <- map.141.plant.2016[map.141.plant.2016$plant=="switchgrass",] 
map.141.plant.misc <- map.141.plant.2016[map.141.plant.2016$plant=="miscanthus",]


map.141.div <- map.141
map.141.div$Richness <- s
map.141.div$Shannon <- h
map.141.div$Pielou <- pielou 
map.141.div.plant <- map.141.div[map.141.div$source=="phyllosphere",]
map.141.div.mis <- map.141.div.plant[map.141.div.plant$plant=="miscanthus",]
map.141.div.swg <- map.141.div.plant[map.141.div.plant$plant=="switchgrass",]

map.141.alpha <- melt(map.141.div, id.vars=c("sequence_name","treatment", "source", "plant", "time_numeric", "sampling_Rdate","year"), measure.vars=c("Richness", "Shannon", "Pielou"))

map.141.alpha.plant <- map.141.alpha[map.141.alpha$source=="phyllosphere",]

map.141.alpha.plant.2016 <- map.141.alpha[map.141.alpha$source=="phyllosphere"&map.141.alpha$year==2016,]

map.141.alpha.plant.2017 <- map.141.alpha[map.141.alpha$source=="phyllosphere"&map.141.alpha$year==2017,]

map.141.misc.2016 <- map.141.plant.2016[map.141.plant.2016$plant=="miscanthus",]
map.141.switch <- map.141.plant[map.141.plant$plant=="switchgrass",]
map.141.switch.2016 <- map.141.plant[map.141.plant$plant=="switchgrass"&map.141.plant$year==2016,]
map.141.switch.2017 <- map.141.plant[map.141.plant$plant=="switchgrass"&map.141.plant$year==2017,]

cor.test(map.141.div.swg$time_numeric, map.141.div.swg$Richness)
cor.test(map.141.div.swg$time_numeric, map.141.div.swg$Shannon)
cor.test(map.141.div.swg$time_numeric, map.141.div.swg$Pielou)

cor.test(map.141.div.mis$time_numeric, map.141.div.mis$Richness)
cor.test(map.141.div.mis$time_numeric, map.141.div.mis$Shannon)
cor.test(map.141.div.mis$time_numeric, map.141.div.mis$Pielou)


plant.141.otu <- otu.141[,map.141$source=="phyllosphere"]
plant.141.otu.2016 <- otu.141[,map.141$source=="phyllosphere"&map.141$year==2016]
plant.141.otu.2017 <- otu.141[,map.141$source=="phyllosphere"&map.141$year==2017]

otu.141.2017 <- otu.141[,map.141$year=="2017"]

misc.141.otu.2016 <- plant.141.otu.2016[,map.141.plant.2016$plant=="miscanthus"]
switch.141.otu <- plant.141.otu[,map.141.plant$plant=="switchgrass"]
switch.141.otu.2016 <- plant.141.otu[,map.141.plant$plant=="switchgrass"&map.141.plant$year==2016]
switch.141.otu.2017 <- plant.141.otu[,map.141.plant$plant=="switchgrass"&map.141.plant$year==2017]

misc.141.unique_times <- unique(map.141.misc.2016$time_numeric)[order(unique(map.141.misc.2016$time_numeric))]

### Phyllosphere switchgrass unique times
switch.141.unique_times <- unique(map.141.switch$time_numeric)[order(unique(map.141.switch$time_numeric))]
switch.141.unique_times.2016 <- unique(map.141.switch.2016$time_numeric)[order(unique(map.141.switch.2016$time_numeric))]
switch.141.unique_times.2017 <-  unique(map.141.switch.2017$time_numeric)[order(unique(map.141.switch.2017$time_numeric))]

# Look at species accumulation
misc.141.accumulation <- rep(1, length(misc.141.unique_times))
z <- NULL
for( i in 1:length(unique(map.141.misc.2016$time_numeric))){
  x <- misc.141.otu.2016[,map.141.misc.2016$time_numeric==misc.141.unique_times[i]]
  y <- matrix(x, ncol=sum(map.141.misc.2016$time_numeric==misc.141.unique_times[i]))
  row.names(y) <- row.names(misc.141.otu.2016)
  z <- c(z, row.names(misc.141.otu.2016[rowSums(y)!=0,]))
  misc.141.accumulation[i] <- length(unique(z))
}

switch.141.accumulation.2016 <- rep(1, length(switch.141.unique_times.2016))
z <- NULL
for( i in 1:length(unique(map.141.switch.2016$time_numeric))){
  x <- switch.141.otu.2016[,map.141.switch.2016$time_numeric==switch.141.unique_times.2016[i]]
  z <- c(z, row.names(x[rowSums(x)!=0,]))
  switch.141.accumulation.2016[i] <- length(unique(z))
}

switch.141.accumulation.2017 <- rep(1, length(switch.141.unique_times.2017))
z <- NULL
for( i in 1:length(unique(map.141.switch.2017$time_numeric))){
  x <- switch.141.otu.2017[,map.141.switch.2017$time_numeric==switch.141.unique_times.2017[i]]
  z <- c(z, row.names(x[rowSums(x)!=0,]))
  switch.141.accumulation.2017[i] <- length(unique(z))
}

NewDates_Switch2017 <- c(as.Date("2016-05-15"), as.Date("2016-06-05"), as.Date("2016-06-26"), as.Date("2016-07-17"),as.Date("2016-08-07"), as.Date("2016-08-28"), as.Date("2016-09-18"))

Spec_accum_total.141 <- NULL
Spec_accum_total.141$Species_Accumulation <- c(misc.141.accumulation, switch.141.accumulation.2016, switch.141.accumulation.2017)
Spec_accum_total.141$Date <- c(unique(map.141.plant.misc$sampling_Rdate)[order(unique(map.141.plant.misc$sampling_Rdate))], unique(map.141.switch.2016$sampling_Rdate)[order(unique(map.141.switch.2016$sampling_Rdate))], NewDates_Switch2017)
Spec_accum_total.141 <- as.data.frame(Spec_accum_total.141)
Spec_accum_total.141$Plant <- c(rep("Miscanthus", length(misc.141.accumulation)), rep("Switchgrass", length(switch.141.accumulation.2016)), rep("Switchgrass", length(switch.141.accumulation.2017)))

Spec_accum_total.141$Year <- c(rep("2016", length(misc.141.accumulation)), rep("2016", length(switch.141.accumulation.2016)), rep("2017", length(switch.141.accumulation.2017)))

Spec_accum_total.141$PlantYear <- paste(Spec_accum_total.141$Plant, Spec_accum_total.141$Year)

switch.141.Richness.2016 <- map.141.alpha.plant.2016[map.141.alpha.plant.2016$plant=="switchgrass"&map.141.alpha.plant.2016$variable=="Richness",]
fit_model <- aov(value~factor(sampling_Rdate),switch.141.Richness.2016)
summary(fit_model)
TukeyHSD(fit_model)

misc.141.Richness <- map.141.alpha.plant[map.141.alpha.plant$plant=="miscanthus"&map.141.alpha.plant$variable=="Richness",]
fit_model <- aov(value~factor(sampling_Rdate),misc.Richness)
summary(fit_model)
TukeyHSD(fit_model)

switch.141.Richness.2017 <- map.141.alpha.plant.2017[map.141.alpha.plant.2017$plant=="switchgrass"&map.141.alpha.plant.2017$variable=="Richness",]
summary(fit_model)
TukeyHSD(fit_model)


Richness.df <- rbind(switch.Richness.2016, switch.Richness.2017, misc.Richness)
Richness.df$Factor <- factor(c(rep("Switchgrass 2016", nrow(switch.Richness.2016)), rep("Switchgrass 2017", nrow(switch.Richness.2017)), rep("Miscanthus 2016", nrow(misc.Richness))), levels = c("Miscanthus 2016", "Switchgrass 2016", "Switchgrass 2017"))
Richness.df$FakeDate <- as.Date(gsub(x=Richness.df$sampling_Rdate, pattern = "2017", replacement = "2016" ))

Richness.141.df <- rbind(switch.141.Richness.2016, switch.141.Richness.2017, misc.141.Richness)
Richness.141.df$Factor <- factor(c(rep("Switchgrass 2016", nrow(switch.141.Richness.2016)), rep("Switchgrass 2017", nrow(switch.141.Richness.2017)), rep("Miscanthus 2016", nrow(misc.141.Richness))), levels = c("Miscanthus 2016", "Switchgrass 2016", "Switchgrass 2017"))

Richness.141.df$FakeDate <- as.Date(gsub(x=Richness.141.df$sampling_Rdate, pattern = "2017", replacement = "2016" ))


Supplemental <- rbind(Richness.df, Richness.141.df)
Supplemental$RarefactionDepth <- c(rep("1000 Reads", nrow(Richness.df)), rep("141 Reads", nrow(Richness.141.df)))
Supplemental$RarefactionDepth <- factor(Supplemental$RarefactionDepth, levels=c("1000 Reads", "141 Reads"))




SupplementaryFigure1 <- ggplot(data = Supplemental, aes(x=FakeDate, y=value))+
  geom_boxplot(mapping=aes(group=FakeDate, fill=Factor, color=as.character(year)), width=8, position="dodge")+
  facet_grid(rows=vars(RarefactionDepth), cols=vars(Factor), scales = "free")+
  scale_fill_manual(values=c("darkgreen", "darkolivegreen3", "White"))+
  scale_color_manual(values=c("black","darkolivegreen3"))+
  scale_x_date(date_breaks = "1 month", labels=date_format("%b") ) +
  ylab(label = "Richness")+
  xlab(label="Date")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 60,vjust = 0.75), legend.position = "none")

ggsave("Figures/FigureS1_AlphaDiversity.eps", SupplementaryFigure1, device = "eps", width=6, height=4, units = "in")






### Soil Alpha Diversity
colSums(otu_soil_rare)

s.soil <- specnumber(otu_soil_rare,MARGIN=2)
h.soil <- diversity(t(otu_soil_rare), index = "shannon")
pielou.soil =h.soil/log(s.soil)

map.soil.div <- map.soil
map.soil.div$Richness <- s.soil
map.soil.div$Shannon <- h.soil
map.soil.div$Pielou <- pielou.soil 


map.soil.alpha <- melt(map.soil.div, id.vars=c("sequence_name","treatment", "source", "plant", "time_numeric", "sampling_Rdate","year"), measure.vars=c("Richness", "Shannon", "Pielou"))


Richness.soil.df <- map.soil.alpha[map.soil.alpha$variable=="Richness",]
Combined.Richness <- rbind(Richness.soil.df, Richness.df[,1:9])

Combined.Richness$Factor <- factor(c(rep("Soil 19,967 Reads", nrow(Richness.soil.df)), rep("Phyllosphere 1,000 Reads", nrow(Richness.df))), levels=c("Phyllosphere 1,000 Reads", "Soil 19,967 Reads"))
Combined.Richness$FakeDate <- as.Date(gsub(pattern="2017", replacement = "2016", Combined.Richness$sampling_Rdate))
Combined.Richness$PlantYear <- factor(paste(Combined.Richness$plant, Combined.Richness$year), levels=c("miscanthus 2016", "switchgrass 2016", "switchgrass 2017"))

Richness.soil.df %>%
  group_by(sampling_Rdate) %>%
  dplyr::summarise(n=sum(abun>0)/length(abun),
                   all=length(abun),
                   rep_ab=mean(abun),
                   sd_rep=sd(abun)
  ) 

max(colSums(otu_soil_rare>0))
  
  


GenSciMeeting.alpha <- ggplot(Combined.Richness, aes(x=FakeDate, y=value))+
  facet_grid(rows=vars(Factor), cols=vars(PlantYear), scales="free_y")+
  geom_boxplot(mapping = aes(group=FakeDate, fill=PlantYear, color=as.character(year)), , width=8, position="dodge")+
  scale_fill_manual(values = c("darkgreen", "darkolivegreen3", "white"))+
  scale_color_manual(values=c("black", "darkolivegreen3"))+
  ylab(label = "Richness")+
  theme_bw()+
  scale_x_date(date_breaks = "1 month", labels=date_format("%b"))+
  theme(axis.text.x = element_text(angle=60, vjust = .25), axis.title.x = element_blank())+
  guides(fill=FALSE, color=FALSE)

ggsave(filename = "~/Desktop/AlphaDiversity_GenSci.eps", GenSciMeeting.alpha, device = "eps", width=6, height=4, units = "in" )

Spec_accum_total$Plant
### Spec Accumulation Curves Phyllosphere 141 Reads

SupplementalAccumCurve <- ggplot(Spec_accum_total.141, aes(x=Date, y=Species_Accumulation)) + 
  geom_point(aes(color=Plant, shape=PlantYear), size=2) +
  geom_line(aes(color=Plant, linetype=Year))+
  scale_shape_manual(values=c(15,16,1))+
  scale_linetype_manual(values = c(1,2))+
  theme_bw()+
  labs(y="Total Observed Taxa") + 
  scale_color_manual(values=c("darkgreen","darkolivegreen3"))  + 
  scale_x_date(date_breaks = "1 month", labels=date_format("%b"))+
  guides(shape=FALSE)



soil.2016.unique.times <- unique(map.soil.2016$sampling_Rdate)
soil.2016.unique.times <- soil.2016.unique.times[order(soil.2016.unique.times)]

soil.2017.unique.times <- unique(map.soil.2017$sampling_Rdate)
soil.2017.unique.times <- soil.2017.unique.times[order(soil.2017.unique.times)]

soil.switch.2016.accum <- rep(1, length(soil.2016.unique.times))
soil.misc.2016.accum <- rep(1, length(soil.2016.unique.times))
soil.switch.2017.accum <- rep(1, length(soil.2017.unique.times))

z.misc <- NULL
z.switch <- NULL
for( i in 1:length(soil.2016.unique.times)){
  x.misc <- otu.soil.2016[,map.soil.2016$sampling_Rdate==soil.2016.unique.times[i]&map.soil.2016$plant=="miscanthus"]
  x.switch <- otu.soil.2016[,map.soil.2016$sampling_Rdate==soil.2016.unique.times[i]&map.soil.2016$plant=="switchgrass"]
  z.misc <- c(z.misc, row.names(x.misc[rowSums(x.misc)!=0,]))
  z.switch <- c(z.switch, row.names(x.switch[rowSums(x.switch)!=0,]))
  soil.switch.2016.accum[i] <- length(unique(z.switch))
  soil.misc.2016.accum[i] <- length(unique(z.misc))
}

z.switch <- NULL
for( i in 1:length(soil.2017.unique.times)){
  x.switch <- otu.soil.2017[,map.soil.2017$sampling_Rdate==soil.2017.unique.times[i]&map.soil.2017$plant=="switchgrass"]
  z.switch <- c(z.switch, row.names(x.switch[rowSums(x.switch)!=0,]))
  soil.switch.2017.accum[i] <- length(unique(z.switch))
}

Soil.Spec.Accum <- data.frame(NumberSpecies = c(soil.switch.2016.accum, soil.misc.2016.accum, soil.switch.2017.accum), Date=c(soil.2016.unique.times, soil.2016.unique.times, soil.2017.unique.times))
Soil.Spec.Accum$Plant <- c(rep("Switchgrass", length(soil.switch.2016.accum)), rep("Miscanthus", length(soil.misc.2016.accum)), rep("Switchgrass", length(soil.switch.2017.accum)))
Soil.Spec.Accum$Year <- c(rep("2016", length(soil.switch.2016.accum)), rep("2016", length(soil.misc.2016.accum)), rep("2017", length(soil.switch.2017.accum)))

Soil.Spec.Accum$FakeDate <- as.Date(gsub(pattern = "2017", replacement = "2016", Soil.Spec.Accum$Date))
Soil.Spec.Accum$PlantYear <- paste(Soil.Spec.Accum$Plant, Soil.Spec.Accum$Year)

Soil_Accum_Curve <- ggplot(Soil.Spec.Accum, aes(x=FakeDate, y=NumberSpecies)) + 
  geom_point(aes(color=Plant,shape=PlantYear), size=2) +
  scale_shape_manual(values = c(15,16,1))+
  theme_bw()+
  geom_path(aes(color=Plant, linetype=Year))+
  scale_linetype_manual(values=c(1,2))+
  labs(y="Total Observed Taxa", x="Date") + 
  scale_color_manual(values=c("darkgreen","darkolivegreen3", "darkolivegreen3"))  + 
  guides(shape=FALSE)+
  guides(color=FALSE, linetype=FALSE)
  scale_x_date(date_breaks = "1 month", labels=date_format("%b"))

ggsave("~/Desktop/SoilAccumCurves_NoGuides.eps", Soil_Accum_Curve, device = "eps", height = 4, width = 4, units = "in" )

rbind(Spec_accum_total, Soil.Spec.Accum)

temp_spec_accum <- Spec_accum_total

colnames(temp_spec_accum) <- c("NumberSpecies", "Date", "Plant", "Line", "Year", "PlantYear")

temp_spec_accum$FakeDate <- temp_spec_accum$Date
temp_spec_accum <- temp_spec_accum[,-4]
temp_spec_accum <- temp_spec_accum[,c(1:4,6,5)]

Combined_Spec_Accum <- rbind(temp_spec_accum, Soil.Spec.Accum)
Combined_Spec_Accum$Source <- c(rep("Phyllosphere 1,000 Reads", nrow(temp_spec_accum)), rep("Soil 19,967 Reads", nrow(Soil.Spec.Accum)))

Both_accum_curve <- ggplot(Combined_Spec_Accum, aes(x=FakeDate, y=NumberSpecies)) + 
  geom_point(aes(color=Plant,shape=PlantYear), size=2) +
  facet_grid(rows = vars(Source), scales="free_y")+
  scale_shape_manual(values = c(15,16,1))+
  theme_bw()+
  geom_path(aes(color=Plant, linetype=Year))+
  scale_linetype_manual(values=c(1,2))+
  labs(y="Total Observed Taxa", x="Date") + 
  scale_color_manual(values=c("darkgreen","darkolivegreen3")) + 
  guides(shape=FALSE)+
  guides(color=FALSE, linetype=FALSE)+
  scale_x_date(date_breaks = "1 month", labels=date_format("%b"))

ggsave("~/Desktop/Combined_accumCurves.eps", Both_accum_curve, height=4, width = 4, device = "eps", units = "in" )

Both_accum_curve_faceted <- ggplot(Combined_Spec_Accum, aes(x=FakeDate, y=NumberSpecies)) + 
  geom_point(aes(color=Plant,shape=PlantYear), size=2) +
  facet_grid(rows = vars(Source), cols=vars(PlantYear), scales="free_y") +
  scale_shape_manual(values = c(15,16,1))+
  theme_bw()+
  geom_path(aes(color=Plant))+
  labs(y="Total Observed Taxa", x="Date") + 
  scale_color_manual(values=c("darkgreen","darkolivegreen3")) + 
  guides(shape=FALSE)+
  guides(color=FALSE, linetype=FALSE)+
  scale_x_date(date_breaks = "1 month", labels=date_format("%b"))+
  theme(axis.text.x=element_text(angle=60, vjust = .25), axis.title.x = element_blank())

ggsave("~/Desktop/AccumulationCurves_Facted.eps", Both_accum_curve_faceted, height=4, width = 6, device = "eps", units = "in" )


###################
### Ordinations ###
###################

### Make Whole PCoA

dist.otu.141 <- vegdist(t(otu.141), method="bray")
pcoa.otu.141 <- cmdscale(dist.otu.141, eig=TRUE)

### Set up dataframe for plotting
points.141 <- pcoa.otu.141$points
points.141 <- as.data.frame(points.141)
colnames(points.141) <- c("Axis1", "Axis2")
points.141$plant <- factor(map.141$plant, levels=c("switchgrass", "miscanthus"))
points.141$source <- factor(map.141$source, levels=c("soil", "phyllosphere"))
points.141$Year <- factor(map.141$Year, levels=c(2016,2017))
points.141$sampling_date <- map.141$sampling_date
points.141 <- left_join(points.141, Date_to_Week, by="sampling_date")
points.141$Week <- factor(points.141$Week, levels=c(1:10))
points.141$SourceYear <- paste(points.141$source,points.141$Year, sep="")
points.141$SourceYear <- factor(points.141$SourceYear, levels=unique(points.141$SourceYear))
points.141$Treatment <- map.141$treatment


points.141$SampleType <- map.141$source
points.141$SampleType[points.141$SampleType=="soil"] <- paste(map.141$source[points.141$source=="soil"], map.141$treatment[points.141$source=="soil"]) 

points.141$SampleType <- factor(points.141$SampleType, levels=c("soil standard fertilization", "soil nitrogen free", "phyllosphere"))

points.141$PlantYear <- factor(paste(map.141$plant, map.141$year, sep=""), levels=c("miscanthus2016", "switchgrass2016", "switchgrass2017"))

### Determine % variation explained on each axis
Ax1.whole.141 <- pcoa.otu.141$eig[1]/sum(pcoa.otu.141$eig)
Ax2.whole.141 <- pcoa.otu.141$eig[2]/sum(pcoa.otu.141$eig)



### ggplot Whole PCoA
FigureS2B <- ggplot(points.141, aes(x=Axis1, y=Axis2))+
  coord_fixed()+
  geom_point(aes(shape=SampleType, color=plant, size=Week, fill=PlantYear))+
  scale_shape_manual(values = c(22,23,21))+
  scale_size_manual(values=Point_Sizes)+
  scale_color_manual(values=c("darkolivegreen3","darkgreen"))+
  scale_fill_manual(values=c("darkgreen", "darkolivegreen3", "white"))+
  xlab(paste(round(Ax1.whole.141, 3)*100, "% Var. Explained", sep=""))+
  ylab(paste(round(Ax2.whole.141, 3)*100, "% Var. Explained", sep=""))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.position = "none")

guides(fill=guide_legend(override.aes=list(shape=21, colour=c(Miscanthus2016="darkgreen",Switchgrass2016="darkolivegreen3", Switchgrass2017="white"), size=5)), color=FALSE, shape=guide_legend(override.aes = list(size=5)))

FigureS2_Legend <- ggplot(points.141, aes(x=Axis1, y=Axis2))+
  coord_fixed()+
  geom_point(aes(shape=SampleType, color=plant, size=Week, fill=PlantYear))+
  scale_shape_manual(values = c(22,23,21))+
  scale_size_manual(values=Point_Sizes)+
  scale_color_manual(values=c("darkolivegreen3","darkgreen"))+
  scale_fill_manual(values=c("darkgreen", "darkolivegreen3", "white"))+
  xlab(paste(round(Ax1.whole.141, 3)*100, "% Var. Explained", sep=""))+
  ylab(paste(round(Ax2.whole.141, 3)*100, "% Var. Explained", sep=""))+
  theme(legend.position = "right", legend.box = "horizontal")+
  guides(fill=guide_legend(override.aes=list(shape=21, color=c(miscanthus2016="darkgreen",switchgrass2016="darkolivegreen3", switchgrass2017="darkolivegreen3"),fill=c(miscanthus2016="darkgreen",switchgrass2016="darkolivegreen3", switchgrass2017="white"), size=5)), color=FALSE, shape=guide_legend(override.aes = list(size=5)))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

legend_S2 <- get_legend(FigureS2_Legend)

#ggsave("Figures/FigureS2A_WholePCoALegend.eps", FigS2A, width = 6, height = 6, device = "eps", units="in")


### For Figure S2B
otu.141.plant <- otu.141[,map.141$source=="phyllosphere"]
dist.141.plant <- vegdist(t(otu.141.plant), method = "bray")
pcoa.plant.141 <- cmdscale(dist.141.plant, eig = TRUE)

plant.dates.141 <- unique(map.141.plant$sampling_date)[order(unique(map.141.plant$sampling_date))]
sg.dates.141 <- unique(map.141.plant.switch$sampling_date)[order(unique(map.141.plant.switch$sampling_date))]
sg.dates.141 <- c(sg.dates.141,unique(map.141.plant.2017$sampling_date)[order(unique(map.141.plant.2017$sampling_date))] )
ms.dates.141 <- unique(map.141.plant.misc$sampling_date)[order(unique(map.141.plant.misc$sampling_date))]

### Set up the points from the "entire" plant PCoA

plant.points.141 <- pcoa.plant.141$points
sg.points.141 <- plant.points.141[map.141.plant$plant=="switchgrass",]
ms.points.141 <- plant.points.141[map.141.plant$plant=="miscanthus",]

### determine the average location of each crop on the PCoA at each timepoint
sg.points.collapsed.141 <- NULL
sg.samples.per.date.141 <- NULL
for (i in 1:length(sg.dates.141)){
  x <- plant.points.141[map.141.plant$sampling_date==sg.dates.141[i]&map.141.plant$plant=="switchgrass",]
  sg.samples.per.date.141 <- c(sg.samples.per.date.141, nrow(x))
  sg.points.collapsed.141 <- rbind(sg.points.collapsed.141, c(colSums(x), sd(x[,1]), sd(x[,2])))
}
sg.points.collapsed.141 <- as.data.frame(sg.points.collapsed.141)
colnames(sg.points.collapsed.141) <- c("Axis1", "Axis2", "sd_axis1", "sd_axis2")
sg.points.collapsed.141$sampling_date <- sg.dates.141
sg.points.collapsed.141[,1:2] <- sg.points.collapsed.141[,1:2]/sg.samples.per.date.141

ms.points.collapsed.141 <- NULL
ms.samples.per.date.141 <- NULL
for(i in 1:length(ms.dates.141)){
  x <- ms.points.141[map.141.plant.misc$sampling_date==ms.dates.141[i],]
  ms.samples.per.date.141 <- c(ms.samples.per.date.141, nrow(x))
  ms.points.collapsed.141 <- rbind(ms.points.collapsed.141, c(colSums(matrix(x, ncol=2)), sd(matrix(x,ncol=2)[,1]), sd(matrix(x, ncol=2)[,2])))
}

ms.points.collapsed.141 <- as.data.frame(ms.points.collapsed.141)
colnames(ms.points.collapsed.141) <- c("Axis1", "Axis2", "sd_axis1", "sd_axis2")
ms.points.collapsed.141$sampling_date <- ms.dates.141
ms.points.collapsed.141[,1:2] <- ms.points.collapsed.141[,1:2]/ms.samples.per.date.141


### Combine the dataframe for us in ggplot
plant.points.collapsed.141 <- rbind(sg.points.collapsed.141, ms.points.collapsed.141)
plant.points.collapsed.141 <- left_join(plant.points.collapsed.141, Date_to_Week, by="sampling_date")
plant.points.collapsed.141$plant <- c(rep("Switchgrass", nrow(sg.points.collapsed.141)), rep("Miscanthus", nrow(ms.points.collapsed.141)))
plant.points.collapsed.141$Week <- factor(x = plant.points.collapsed.141$Week, levels=c(1,2,3,4,5,6,7,8,9,10))
plant.points.collapsed.141$WeekNumeric <- as.numeric(plant.points.collapsed.141$Week)
plant.points.collapsed.141$Year <- c(rep(2016,8), rep(2017,7), rep(2016,9))
plant.points.collapsed.141$Year <- factor(plant.points.collapsed.141$Year, levels = c(2016,2017))
plant.points.collapsed.141$plant <- factor(plant.points.collapsed.141$plant, levels = c("Switchgrass", "Miscanthus"))
plant.points.collapsed.141 <- plant.points.collapsed.141[order(plant.points.collapsed.141$sampling_date),]

### Set up the weather maps for Weather Environmental fits
switch.weather.141 <- map.141.plant.switch[,c("precipitation", "Air_temp_mean", "air_temp_max", "Air_Temp_Min", "Air_Pressure", "RH", "AH", "Wind_Speed_Mean", "Solar_Radiation", "PAR", "soil_temp_5_cm_bare_avg", "sampling_date", "time_numeric")]
switch.weather.141 <- unique(switch.weather.141)

switch.weather.2017.141 <- map.141.plant.2017[,c("precipitation", "Air_temp_mean", "air_temp_max", "Air_Temp_Min", "Air_Pressure", "RH", "AH", "Wind_Speed_Mean", "Solar_Radiation", "PAR", "soil_temp_5_cm_bare_avg", "sampling_date", "time_numeric")]
switch.weather.2017.141 <- unique(switch.weather.2017.141)

misc.weather.141 <- map.141.plant.misc[,c("precipitation", "Air_temp_mean", "air_temp_max", "Air_Temp_Min", "Air_Pressure", "RH", "AH", "Wind_Speed_Mean", "Solar_Radiation", "PAR", "soil_temp_5_cm_bare_avg", "sampling_date", "time_numeric")]
misc.weather.141 <- unique(misc.weather.141)

phyllo.weather.141 <- rbind(switch.weather.141, switch.weather.2017.141, misc.weather.141)
phyllo.weather.141 <- phyllo.weather.141[order(phyllo.weather.141$sampling_date),]
phyllo.weather.141$sampling_date == plant.points.collapsed.141$sampling_date
phyllo.weather.141$Week <- plant.points.collapsed.141$WeekNumeric


### Set up EnvFit for each category of contextual data (weather, leaf chemistry, soil chemistry)
envfit.phyllo.weather.141 <- envfit(plant.points.collapsed.141[,1:2], phyllo.weather.141[,c(1:11,14)])
envfit.phyllo.weather.df.141<-as.data.frame(scores(envfit.phyllo.weather.141, display = "vectors"))
colnames(envfit.phyllo.weather.df.141) <- c("Dim1", "Dim2")




phyllo.leaf.chemistry.141 <- map.141.plant[,c("LDMC_mg_per_g", "nitrogen_percent", "carbon_percent", "carbon_per_nitrogen", "height_mean_cm")]
envfit.phyllo.lc.141 <- envfit(pcoa.plant.141$points, phyllo.leaf.chemistry.141)
envfit.phyllo.lc.df.141<-as.data.frame(scores(envfit.phyllo.lc.141, display = "vectors"))


phyllo.soil.chemsitry.141 <- map.141.plant[,c("pH", "P_ppm", "K_ppm", "Ca_ppm", "Mg_ppm", "organic_matter", "NO3N_ppm", "NH4_ppm", "soil_moisture_percent", "soil_temp_10cm")]
envfit.phyllo.sc.141 <- envfit(pcoa.plant.141$points, phyllo.soil.chemsitry.141)
envfit.phyllo.sc.df.141<-as.data.frame(scores(envfit.phyllo.sc.141, display = "vectors"))

envfit.phyllo.total.141 <- rbind(envfit.phyllo.weather.df.141, envfit.phyllo.lc.df.141, envfit.phyllo.sc.df.141)

envfit.phyllo.total.141$r2 <- c(envfit.phyllo.weather.141$vectors$r, envfit.phyllo.lc.141$vectors$r, envfit.phyllo.sc.141$vectors$r)
envfit.phyllo.total.141$pval <- c(envfit.phyllo.weather.141$vectors$pvals, envfit.phyllo.lc.141$vectors$pvals, envfit.phyllo.sc.141$vectors$pvals)

arrow_scaling <- scale_arrow(arrows=envfit.phyllo.total.141[,1:2], data=plant.points.collapsed.141)

envfit.phyllo.total.141$Axis1 <- envfit.phyllo.total.141$Dim1 * arrow_scaling
envfit.phyllo.total.141$Axis2 <- envfit.phyllo.total.141$Dim2 * arrow_scaling
envfit.phyllo.total.141$Variable <- row.names(envfit.phyllo.total.141)

### Subset to p <0.05 & Rsquared >0.4
envfit.phyllo.sub.141 <- subset(envfit.phyllo.total.141, envfit.phyllo.total.141$r2 >0.3999) 
envfit.phyllo.sub.141 <- subset (envfit.phyllo.sub.141 , envfit.phyllo.sub.141$pval<0.05)


write.table(file = "Figures/TableS2_141EnvFit.txt", x=envfit.phyllo.total.141[,1:4], sep="\t", quote=FALSE)

Point_Sizes.141 <- seq(from=2, to=6, length.out = 10)
Ax1.plant.141 <- pcoa.plant.141$eig[1]/sum(pcoa.plant.141$eig)
Ax2.plant.141 <- pcoa.plant.141$eig[2]/sum(pcoa.plant.141$eig)

plant.points.collapsed.141$PlantYear <- factor(paste(plant.points.collapsed.141$plant,plant.points.collapsed.141$Year, sep=""), levels=c("Miscanthus2016", "Switchgrass2016", "Switchgrass2017"))

### Plot Plant PCoA
FigureS2A <- ggplot(plant.points.collapsed.141, aes(x=Axis1, y=Axis2))+
  geom_point(aes(size=Week, color=plant, fill=PlantYear),shape=21) +
  scale_color_manual(values=c("darkolivegreen3","darkgreen"))+
  scale_fill_manual(values=c("darkgreen","darkolivegreen3", "white"))+
  scale_size_manual(values=Point_Sizes)+
  coord_fixed()+
  geom_segment(data=plant.points.collapsed.141, aes(x=Axis1,xend=Axis1+sd_axis1,y=Axis2,yend=Axis2,color=plant ))+
  geom_segment(data=plant.points.collapsed.141, aes(x=Axis1,xend=Axis1-sd_axis1,y=Axis2,yend=Axis2,color=plant))+
  geom_segment(data=plant.points.collapsed.141, aes(x=Axis1,xend=Axis1,y=Axis2,yend=Axis2+sd_axis2,color=plant))+
  geom_segment(data=plant.points.collapsed.141, aes(x=Axis1,xend=Axis1,y=Axis2,yend=Axis2-sd_axis2,color=plant))+
  xlab(label = paste(round(Ax1.plant.141,digits = 3)*100, "% Var. Explained", sep = ""))+
  ylab(label= paste(round(Ax2.plant.141,digits = 3)*100, "% Var. Explained", sep = ""))+
  geom_segment(data = envfit.phyllo.sub.141,
               aes(x = 0, xend = Axis1, y = 0, yend = Axis2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "grey")+
  geom_text(data=envfit.phyllo.sub.141, aes(label=Variable))+
  theme(panel.grid.major = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.position = "none")

#ggsave(filename = "Figures/FigureS2B_PhyllospherePCoA141.eps", FigureS2B, device = "eps", height = 6, width = 6, units = "in")

soil.dist <- vegdist(t(otu_soil_rare), method="bray")
soil.dist.2016 <- vegdist(t(otu.soil.rare.2016), method="bray")
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

Collapsed_Soil.2016 <- rbind(misc.soil_MAIN.2016[1:10,], misc.soil_NF.2016[1:10,], switch.soil_MAIN.2016[1:9,], switch.soil_NF.2016[1:9,])
Collapsed_Soil.2016$plant <- c(rep("miscanthus", 20), rep("switchgrass",18))
Collapsed_Soil.2016$Date <- c(dates[1:10],dates[1:10],dates[1:9], dates[1:9])


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

setEPS()
postscript("Figures/FigureS2.eps", width=10, height=10,pointsize=10, paper="special",)
plot_grid(FigureS2A,FigureS2B, FigureS2C,  legend_S2, axis="bl")
dev.off()

align = "h", label_size = 10, labels = c("A", "B", "C") 

plot_grid(plot_grid(FigureS2A, FigureS2B, ncol=2, align = "h"), plot_grid(FigureS2C, legend_S2, ncol=2, align="h"), ncol=1, align = "v")

####################################################
### Look at betadispersion of crops through time ###
####################################################

otu.switch.2016 <- otu.plant.2016[,map.plant.2016$plant=="switchgrass"]
map.plant.switch.2016 <- map.plant.switch[map.plant.switch$Year==2016,]

otu.switch.2017 <- otu.plant.2017
map.plant.switch.2017 <- map.plant.2017

dist.switch.2016 <- vegdist(t(otu.switch.2016), method="bray")
dist.switch.2017 <- vegdist(t(otu.switch.2017), method="bray")

Dispersion.switch.2016 <- betadisper(dist.switch.2016, map.plant.switch.2016$sampling_date)
Dispersion.switch.2017 <- betadisper(dist.switch.2017, map.plant.switch.2017$sampling_date)


names(Dispersion.switch.2016$distances)==map.plant.switch.2016$sequence_name
Dispersion.switch.2016.df <- data.frame(Distance_to_Median=Dispersion.switch.2016$distances, Date=map.plant.switch.2016$sampling_date)

names(Dispersion.switch.2017$distances)==map.plant.switch.2017$sequence_name
Dispersion.switch.2017.df <- data.frame(Distance_to_Median=Dispersion.switch.2017$distances, Date=map.plant.switch.2017$sampling_date)

SG_Disp.2016 <- ggplot(Dispersion.switch.2016.df, aes(x=Date, y=Distance_to_Median))+
  geom_boxplot()+
  geom_point()+
  ggtitle("2016 Switchgrass Beta Dispersion")+
  theme(axis.text.x = element_text(angle = 60))+
  ylim(c(0,1))

ggsave("Figures/2016SG_BetaDispersion_1000.eps", SG_Disp.2016, device = "eps", width = 4, height = 4, units = "in")
TukeyHSD(aov(data = Dispersion.switch.df, Distance_to_Median~Date))

SG_Disp.2017 <- ggplot(Dispersion.switch.2017.df, aes(x=Date, y=Distance_to_Median))+
  geom_boxplot()+
  geom_point()+
  ggtitle("2017 Switchgrass Beta Dispersion")+
  theme(axis.text.x = element_text(angle = 60))+
  ylim(c(0,1))

ggsave("Figures/2017SG_BetaDispersion_1000.eps", SG_Disp.2017, device = "eps", width = 4, height = 4, units = "in")


otu.misc.2016 <- otu.plant.2016[,map.plant.2016$plant=="miscanthus"]
dist.misc.2016 <- vegdist(t(otu.misc.2016), method = "bray")
Dispersion.misc.2016 <- betadisper(dist.misc.2016, map.plant.misc$sampling_date)
names(Dispersion.misc.2016$distances)==misc.map.2016$sequence_name
Dispersion.misc.2016.df <- data.frame(Distance_to_Median=Dispersion.misc.2016$distances, Date=misc.map.2016$sampling_date)

MS_Disp <- ggplot(Dispersion.misc.2016.df, aes(x=Date, y=Distance_to_Median))+
  geom_boxplot()+
  geom_point()+
  ggtitle("2016 Miscanthus Beta Dispersion")+
  theme(axis.text.x = element_text(angle = 60))+
  ylim(c(0,1))

ggsave("Figures/MS_BetaDispersion_1000.eps", device = "eps", width = 4, height = 4, units = "in")



otu.141.misc <- otu.141.plant[,map.141.plant$plant=="miscanthus"]
dist.misc.141 <- vegdist(t(otu.141.misc), method = "bray")
Dispersion.misc.141 <- betadisper(dist.misc.141, map.141.plant.misc$sampling_date)
names(Dispersion.misc.141$distances)==map.141.plant.misc$sequence_name
Dispersion.misc.141.df <- data.frame(Distance_to_Median=Dispersion.misc.141$distances, Date=map.141.plant.misc$sampling_date)

otu.switch.141.2016 <- otu.141.plant[,map.141.plant$plant=="switchgrass"&map.141.plant$Year==2016]
map.141.plant.switch.2016 <- map.141.plant.switch[map.141.plant.switch$Year==2016,]

otu.switch.141.2017 <- otu.141.plant[,map.141.plant$Year==2017]
map.141.plant.switch.2017 <- map.141.plant.2017

dist.switch.2016.141 <- vegdist(t(otu.switch.141.2016), method="bray")
dist.switch.2017.141 <- vegdist(t(otu.switch.141.2017), method="bray")

Dispersion.switch.2016.141 <- betadisper(dist.switch.2016.141, map.141.plant.switch.2016$sampling_date)
Dispersion.switch.2017.141 <- betadisper(dist.switch.2017.141, map.141.plant.switch.2017$sampling_date)


names(Dispersion.switch.2016.141$distances)==map.141.plant.switch.2016$sequence_name
Dispersion.switch.2016.df.141 <- data.frame(Distance_to_Median=Dispersion.switch.2016.141$distances, Date=map.141.plant.switch.2016$sampling_date)

names(Dispersion.switch.2017.141$distances)==map.141.plant.switch.2017$sequence_name
Dispersion.switch.2017.df.141 <- data.frame(Distance_to_Median=Dispersion.switch.2017.141$distances, Date=map.141.plant.switch.2017$sampling_date)

AllDispersion <- rbind(Dispersion.switch.2016.df,Dispersion.switch.2017.df, Dispersion.misc.2016.df, Dispersion.switch.2016.df.141, Dispersion.switch.2017.df.141, Dispersion.misc.141.df)
AllDispersion$Date <- gsub(pattern="2017", replacement = "2016", AllDispersion$Date)
AllDispersion$NumberSequences <- factor(c(rep("1000 Reads", 140), rep("141 Reads", 171)), levels=c("1000 Reads", "141 Reads"))

AllDispersion$Factor <- factor(x = c(rep("Switchgrass 2016", nrow(Dispersion.switch.2016.df)), rep("Switchgrass 2017", nrow(Dispersion.switch.2017.df)), rep("Miscanthus 2016", nrow(Dispersion.misc.2016.df)),rep("Switchgrass 2016", nrow(Dispersion.switch.2016.df.141)), rep("Switchgrass 2017", nrow(Dispersion.switch.2017.df.141)), rep("Miscanthus 2016", nrow(Dispersion.misc.141.df)) ), levels=c("Miscanthus 2016", "Switchgrass 2016", "Switchgrass 2017"))

AllDispersion$Year <- factor(x = c(rep("2016", nrow(Dispersion.switch.2016.df)), rep("2017", nrow(Dispersion.switch.2017.df)), rep("2016", nrow(Dispersion.misc.2016.df)),rep("2016", nrow(Dispersion.switch.2016.df.141)), rep("2017", nrow(Dispersion.switch.2017.df.141)), rep("2016", nrow(Dispersion.misc.141.df)) ), levels=c("2016", "2017"))

AllDispersion$Date <- as.Date(AllDispersion$Date)

SupplementaryFigure3 <- ggplot(AllDispersion, aes(x=Date, y=Distance_to_Median))+
  facet_grid(rows=vars(NumberSequences), cols=vars(Factor))+
  geom_boxplot(width=8, mapping=aes(group=Date, fill=Factor, color=Year))+
  scale_fill_manual(values=c("darkgreen", "darkolivegreen3", "White"))+
  scale_color_manual(values=c("black","darkolivegreen3"))+
  scale_x_date(date_breaks = "1 month", labels=date_format("%b") ) +
  ylab(label = "Beta Dispersion")+
  xlab(label="Date")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 60,vjust = 0.75), legend.position = "none")

ggsave("Figures/FigureS3_BetaDispersion.eps", SupplementaryFigure3, device = "eps", width=6, height=4, units = "in")

##### Mean Pairwise distance between crops through time

unique(map.plant.switch$sampling_date)

dates.2016 <- unique(map.plant.misc$sampling_date)[unique(map.plant.misc$sampling_date)%in%unique(map.plant.switch$sampling_date)]


pairwise.distance <- NULL
for (i in 1:length(dates.2016)){
  x <- otu_rare[,map_16S$sampling_date==dates.2016[i]&map_16S$source=="phyllosphere"]
  y <- map_16S[map_16S$sampling_date==dates.2016[i]&map_16S$source=="phyllosphere",]
  x.bc <- vegdist(t(x), method="bray")
  x.bc.matrix <- as.matrix(x.bc)
  x.bc.matrix.sub <- x.bc.matrix[y$plant=="switchgrass", y$plant=="miscanthus"]
  z <- melt(x.bc.matrix.sub)
  z$Date <- dates.2016[i]
  pairwise.distance <- rbind(pairwise.distance, z)
}

Phyllo_PWD <- ggplot(pairwise.distance, aes(x=Date, y=value))+
  geom_violin()+
  geom_point()+
  ylab(label = "Bray Curtis Dissimilarity between Miscanthus and Switchgrass")+
  theme(axis.text.x = element_text(angle=60, size=8), axis.title.x=element_text(size=10), axis.title.y = element_text(size=10))

ggsave("Figures/FigureS4_BCBtwnCrops.eps", Phyllo_PWD, device = "eps", width = 6, height = 6, units = "in")

library(limma)
venn_switch.2016 <- 1*rowSums(switch.rare.otu.2016)>0
venn_misc.2016 <-  1*rowSums(misc_rare_otu.2016)>0
venn_switch.2017 <- 1*rowSums(switch.rare.otu.2017)>0

venn_phyllo_data <- cbind(venn_switch.2016, venn_switch.2017, venn_misc.2016)
colnames(venn_phyllo_data) <- c("Switchgrass 2016", "Switchgrass 2017", "Miscanthus 2016")
venn_phyllo_data=venn_phyllo_data[rowSums(venn_phyllo_data)>0,]
v_phyllo=vennCounts(venn_phyllo_data)
v_phyllo_2=round(v_phyllo[,"Counts"]/sum(v_phyllo[,"Counts"]),2) #calculate percentage of each group

setEPS()
postscript("Figures/FigureS6_VennDiagram.eps", width=6, height=6,pointsize=10, paper="special",)
vennDiagram(v_phyllo, circle.col = c('darkolivegreen3', 'darkolivegreen3','darkgreen'), lwd=6, cex=1.2, scale=F)
dev.off()


dates.141.2016 <- unique(map.141.misc.2016$sampling_date)[unique(map.141.misc.2016$sampling_date)%in%unique(map.141.switch.2016$sampling_date)]


pairwise.141.distance <- NULL
for (i in 1:length(dates.141.2016)){
  x <- otu.141[,map.141$sampling_date==dates.141.2016[i]&map.141$source=="phyllosphere"]
  y <- map.141[map.141$sampling_date==dates.141.2016[i]&map.141$source=="phyllosphere",]
  x.bc <- vegdist(t(x), method="bray")
  x.bc.matrix <- as.matrix(x.bc)
  x.bc.matrix.sub <- x.bc.matrix[y$plant=="switchgrass", y$plant=="miscanthus"]
  z <- melt(x.bc.matrix.sub)
  z$Date <- dates.141.2016[i]
  pairwise.141.distance <- rbind(pairwise.141.distance, z)
}

Phyllo_PWD.141 <- ggplot(pairwise.141.distance, aes(x=Date, y=value))+
  geom_violin()+
  geom_point()+
  ylab(label = "Bray Curtis Dissimilarity between Miscanthus and Switchgrass")+
  theme(axis.text.x = element_text(angle=60, size=8), axis.title.x=element_text(size=10), axis.title.y = element_text(size=10))

ggsave("Figures/FigureS4_BCBtwnCrops.141.eps", Phyllo_PWD.141, device = "eps", width = 6, height = 6, units = "in")


pairwise.141.distance$Reads <- rep("141 Reads", nrow(pairwise.141.distance))
pairwise.distance$Reads <- rep("1,000 Reads", nrow(pairwise.distance))

pwd.combined <- rbind(pairwise.141.distance, pairwise.distance)

pwd.combined$Reads <- factor(pwd.combined$Reads, levels=c("141 Reads", "1,000 Reads"))

FigureS5 <- ggplot(pwd.combined, aes(x=Date, y=value))+
  geom_violin()+
  geom_point()+
  facet_grid(cols=vars(Reads))+
  ylab(label = "Bray Curtis Dissimilarity between Miscanthus and Switchgrass")+
  theme(axis.text.x = element_text(angle=60, size=8, vjust = 1, hjust = 1), axis.title.x=element_text(size=10), axis.title.y = element_text(size=10))

ggsave(filename = "Figures/FigureS5_BCDissimilarity.eps", width=8, height = 5, units = "in", device = "eps")

#################################
### Variance Partitioning 141 ###
#################################
Plant_Core_OTU.141 <- plant.141.otu.2016[row.names(plant.141.otu.2016)%in%core.taxa.16,]
Plant_Core_OTU.141.rel <- decostand(Plant_Core_OTU.141, MARGIN=2, method = "total")
Plant_Core.141.dist <- vegdist(t(Plant_Core_OTU.141.rel), method = "bray")

var_part_map.141 <- map.141[map.141$source=="phyllosphere"&map.141$year==2016,c("precipitation", "Air_temp_mean", "air_temp_max", "Air_Temp_Min", "Air_Pressure", "RH", "AH", "Wind_Speed_Mean", "Solar_Radiation", "PAR", "soil_temp_5_cm_bare_avg", "soil_temp_5_cm_sod_avg", "LDMC_mg_per_g", "nitrogen_percent", "carbon_percent", "carbon_per_nitrogen", "height_mean_cm","pH", "P_ppm", "K_ppm", "Ca_ppm", "Mg_ppm", "organic_matter", "NO3N_ppm", "NH4_ppm", "soil_moisture_percent", "soil_temp_10cm")]

colinearity_time.141 <- data.frame(Variable=rep(NA, ncol(var_part_map.141)), T=rep(NA, ncol(var_part_map.141)), pvalue=rep(NA, ncol(var_part_map.141)), estimate=rep(NA, ncol(var_part_map.141)))
for(i in 1:ncol(var_part_map.141)){
  x <- cor.test(map.141$time_numeric[map.141$source=="phyllosphere"&map.141$year==2016], var_part_map.141[,i])
  colinearity_time.141$Variable[i] <- colnames(var_part_map.141)[i]
  colinearity_time.141$T[i] <- x$statistic 
  colinearity_time.141$pvalue[i] <- x$p.value
  colinearity_time.141$estimate[i] <- x$estimate
}

sig_colinear.141 <- colinearity_time.141[colinearity_time.141$pvalue<0.05, 1]

reduced_var_part_ABIOTIC.141 <- var_part_map.141[,!colnames(var_part_map.141)%in%sig_colinear.141]
ABIOTIC.141 <- reduced_var_part_ABIOTIC.141

ABIOTIC_Bigger <- var_part_map[,!colnames(var_part_map)%in%sig_colinear.141]


reduced_sample_names.141 <- unlist(strsplit(colnames(plant.141.otu.2016), split="_"))

reduced_sample_names.141 <- rep(NA, ncol(plant.141.otu.2016))
for(i in 1:ncol(plant.141.otu.2016)) {
  reduced_sample_names.141[i] <-  unlist(strsplit(colnames(plant.141.otu.2016)[i], split="_"))[1] 
}

plots.141 <- unique(reduced_sample_names.141)
Big_Spatial_Distance_Matrix.141 <- NULL

for(i in 1:length(plots.141)){
  z <- matrix(nrow = sum(reduced_sample_names.141==plots.141[i]), ncol=length(reduced_sample_names.141))
  row.names(z) <- colnames(plant.141.otu.2016)[reduced_sample_names.141==plots.141[i]]
  spatial_x <- spatial[spatial$Site1==plots.141[i],]
  for (w in 1:length(plots.141)){
    z[,reduced_sample_names.141==plots.141[w]] <- spatial_x[w,3]
  }
  Big_Spatial_Distance_Matrix.141 <- rbind(Big_Spatial_Distance_Matrix.141,z)
}

colnames(Big_Spatial_Distance_Matrix.141) <- row.names(Big_Spatial_Distance_Matrix.141)

SPACE.141 <- Big_Spatial_Distance_Matrix.141
mantel(Plant_Core.141.dist, SPACE.141)

HOST.141 <- 1*(map.141[map.141$source=="phyllosphere"&map.141$year==2016, "plant"]=="switchgrass")

TIME.141 <- map.141[map.141$source=="phyllosphere"&map.141$year==2016, "time_numeric"]

cca_core.rel.141 <- varpart(Y=t(Plant_Core_OTU.141.rel), HOST.141, ABIOTIC.141, TIME.141)
plot(cca_core.rel.141, Xnames=c("HOST.rel", "ABIOTIC.rel", "TIME.rel"), cutoff=-1, main="Relativized 141 reads")

var.core.dist.141 <- varpart(Y=Plant_Core.141.dist, HOST.141, ABIOTIC.141, TIME.141)
plot(var.core.dist.141, Xnames=c("HOST", "ABIOTIC", "TIME"), cutoff=-1)


Plant_Core_OTU.141.rel.1kmatch <- Plant_Core_OTU.141.rel[,colnames(Plant_Core_OTU.141.rel)%in%colnames(Plant_Core_OTU)]

colnames(Plant_Core_OTU.141.rel.1kmatch) == map.plant.2016$sequence_name

ak <- varpart(Y=t(Plant_Core_OTU.141.rel.1kmatch), HOST, ABIOTIC, TIME)
plot(ak, Xnames=c("HOST", "ABIOTIC", "TIME"))



##################################
#####  End Jackson Analysis  #####
##################################