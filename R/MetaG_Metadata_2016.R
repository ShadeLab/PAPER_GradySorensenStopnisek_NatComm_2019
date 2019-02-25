setwd("~/Documents/Github/PAPER_GradySorensenStopnisek_InPrep")
library(dplyr)
library(tidyr)
library(reshape2)
library(RSQLite)
library(stringr)


# Preping the environmental metadata
glbrc <- dbConnect(RSQLite::SQLite(), "R/InputFiles/GLBRC_bioenergy_db.db" )

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
metadata <- full_join(glbrc_sampling, glbrc_plot, by='plotID')
metadata <- full_join(metadata, glbrc_soil, by='sampleID')
metadata <- full_join(metadata, glbrc_plant, by='sampleID')
metadata <- full_join(metadata, glbrc_NA, by='sampleID')
metadata <- full_join(metadata, glbrc_sequncing, by='nucleic_acid_name')

map_full <- metadata

# #creating numeric time column
map_full$sampling_date <- paste0(map_full$month,'-', map_full$day,'-',map_full$year)
map_full$sampling_date <- as.POSIXct(map_full$sampling_date, format='%m-%d-%Y')
time <- as.POSIXct(map_full$sampling_date, format='%m-%d-%Y')
time_numeric <- as.numeric(time)
map_time <- cbind(map_full, time_numeric)

#adding weather data
weather <- read.csv("R/InputFiles/kbs_weather_09212017.csv", encoding = 'UTF-8')
dim(weather)
head(weather)
weather$sampling_date <- as.POSIXct(weather$date, format='%d.%m.%y')

#subsetting weather file for sample dates
sub_weather <- weather[weather$sampling_date %in% map_time$sampling_date,] 

#merging dataframes - map file and weather
map_complete <- full_join(map_time, weather)


#subset metadata for metagenomes, unique, and 2016
metag_2016 <- subset(map_complete, exclude_from_analysis == "N" & sequencing_type == 'Illumina Hiseq Metagenome Sequencing' & year == "2016")
dim(metag_2016)

# write out file
write.table(metag_2016, file = "R/Outputs/2016_MetaG_map.txt", sep = "\t", row.names=FALSE, quote=FALSE)

