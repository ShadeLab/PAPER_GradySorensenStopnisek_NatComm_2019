#script by john guittar, guittarj@gmail.com, february 2019

#### setup ####

#setwd
wd <- 'C:\\Users\\John\\Documents\\msu\\phyllosphere_soil_metacommunity\\PAPER_GradySorensenStopnisek_InPrep\\R'
setwd(wd)
library(tidyverse)
library(data.table)
library(kableExtra)
library(vegan)
library(cowplot)

#load metadata
meta <- read.table('InputFiles\\Map_GLBRC_16S.txt', sep = '\t') %>% 
  transmute(
    plot = plotID, 
    sp = plant,
    source = source,
    date = as.Date(sampling_date),
    year = Year,
    sample = sequence_name) %>%
  mutate_if(is.factor, as.character)

otus_unrarefied <- read.table('InputFiles\\Filtered_OTU_Table_JacksonMarch11.txt', sep = '\t', header = TRUE)
otus_unrarefied <- data.frame(otu = row.names(otus_unrarefied), otus_unrarefied) %>%
  gather(sample, abun, -otu)

otus <- read.table('InputFiles\\Rarefied_OTU_JacksonMarch11.txt', sep = '\t', header = TRUE)
otus <- data.frame(otu = row.names(otus), otus) %>%
  gather(sample, abun, -otu)

#remove redundant samples; G5R3_NF_09MAY2016_LD2 is taken at the same time/place as G5R3_NF_09MAY2016_LD1; 
#G5R4_MAIN_12SEP2016_LD2
otus_unrarefied <- filter(otus_unrarefied, !sample %in% c('G5R3_NF_09MAY2016_LD2', 'G5R3_NF_09MAY2016_LD2'))
otus <- filter(otus, !sample %in% c('G5R3_NF_09MAY2016_LD2', 'G5R3_NF_09MAY2016_LD2'))
meta <- filter(meta, sample %in% otus$sample | sample %in% otus_unrarefied$sample)

#load cores from Nejc
cores <- read.csv('InputFiles\\core_list.csv', stringsAsFactors = FALSE)

#do you want to re-run simulations or use last saved run?
recalculate <- FALSE

#create soil pool -- all samples combined over both years
soilpool <- otus_unrarefied %>% 
  left_join(meta[, c('sample','source')], by = "sample") %>% 
  filter(source == 'soil') %>% 
  group_by(otu) %>% 
  summarise(abun = sum(abun)) %>%
  filter(abun > 0)

#### null model ####
###identify new arriving OTUS for each plot at each time point
#prepare dataframe columns used to compare pairs of sequential samples
imms <- otus %>% 
  left_join(meta, by = "sample") %>%
  filter(source == 'phyllosphere') %>%
  filter(abun > 0) %>%
  group_by(plot, year) %>%
  arrange(plot, date) %>%
  rename(date1 = date) %>%
  mutate(seqdate = as.numeric(as.factor(date1)),
         date0 = as.Date(NA),
         new = NA)

# for each otu at each plot...
for(i in 1:nrow(imms)) {
  
  #identify the prior sampling date
  #NA if it was the first sample that year
  seqdate0 <- imms$seqdate[i] - 1
  
  imms$date0[i] <- imms$date1[
    imms$plot == imms$plot[i] &
    imms$year == imms$year[i] &
    imms$seqdate == seqdate0][1]
  
  #Determine if each current OTU was present at the earlier sampling date 
  imms$new[i] <- ifelse(
    imms$otu[i] %in% imms$otu[imms$plot == imms$plot[i] &
    imms$year == imms$year[i] &
    imms$seqdate == seqdate0], 
    FALSE, TRUE)
}

#filter out first samples, and non-new spp
imms <- imms %>%
  filter(!is.na(date0)) %>%
  group_by(sample) %>%
  filter(new) %>%
  select(-new)

#count new spp per plot per date
newspp <- imms %>%
  group_by(plot, year, date0, date1, seqdate) %>%
  summarise(n = length(unique(otu))) %>%
  ungroup()

#calcute median abundance of immigrants after arrival
imm_abun <- median(imms$abun)

#calculate numbers of 0.1% incremental increases and decreases in abundance *of non-immigrants* for each plot between each pair of sequential samples
#importantly, all otu zeroes are included
#the null model will work on a per plot per sample basis
#takes 10 minutes because lazy coding...
if (recalculate) {
  abun_shifts <- otus %>%
    left_join(meta, by = "sample") %>%
    filter(source == 'phyllosphere') %>%
    group_by(plot, year, otu) %>%
    arrange(plot, otu, date) %>%
    do(data.frame(stringsAsFactors = FALSE,
      plot = .$plot[1],
      year = .$year[1],
      otu = .$otu[1],
      date0 = .$date[-length(.$date)],
      date1 = .$date[-1],
      diff = diff(.$abun))) %>%
    group_by(plot, year, date0, date1) %>%
    filter(!otu %in% imms$otu[imms$plot == plot[1] & imms$date1 == date1[1]]) %>%
    summarise(diff = sum(diff[diff > 0])) %>%
    ungroup()
  
  saveRDS(abun_shifts, file = 'Outputs\\abun_shifts.RDS')

}

#load abundance_shift data
abun_shifts <- readRDS('Outputs\\abun_shifts.RDS')

#generate list of randomly sampled "new arrivals" from the soil, which match the number of new tax observed at each time step
arrivals <- list()
for(i in 1:nrow(newspp)) {
    arrivals[[i]] <- data.frame(stringsAsFactors = FALSE,
      plot = newspp$plot[i],
      year = newspp$year[i],
      date0 = newspp$date0[i],
      date1 = newspp$date1[i], 
      otu = as.character(sample(soilpool$otu, newspp$n[i], replace = T, prob = soilpool$abun)),
      abun = imm_abun)
}
arrivals <- bind_rows(arrivals)

#reorganize data into a list of observed communities by date
comms <- otus %>%
  left_join(meta, by = "sample") %>%
  group_by(plot, date) %>%
  filter(source == 'phyllosphere') %>%
  filter(abun > 0) %>%
  ungroup() %>%
  mutate(by = paste(plot, date))

#save a version of initial observations for later comparison
comms0 <- comms

##split() breaks data into a list, which will be looped over and appended to in the simulation
#last data don't need to be here because... there is nothing to simulate.
comms <- comms %>%
  group_by(plot, year) %>%
  filter(date < max(date)) %>% 
  split(.$by)

# 1) loop over each plot-date combo (list item in comms), 
# 2) randomly identify x individuals of otus to increase/decrease by 0.1% increments of relative abundance, where x is equal to the number of increases/decreases observed between time points
# 3) randomly sample x otus from the soil to immigrate onto the leaf, where n is equal to the number of novel otus observed on that leaf at that sampling point. They arrive with the median observed abundance of immigrants. I ensure that the taxa are new, and not already present in the phyllosphere (a conservative measure)
# 4) append and loop over until all dates/plots are simulated

#for each community-plot
if (recalculate) {
  for(i in 1:length(comms)) {
  
    #print progress
    print(paste(i, 'of', length(comms))); flush.console()
    
    #rename for convenience and clean up
    j <- comms[[i]] %>% rename(date0 = date) %>% select(-by) %>% mutate_if(is.factor, as.character)
    
    #determine how many increases/decreases there will be
    diffs <- j %>%
      left_join(transmute(abun_shifts, plot, year, date0, diff), by = c("plot", "date0", "year")) %>%
      summarise(diff = length(diff)) %>%
      pull(unique(diff))
    
    #determine how many immigrants are slated to arrive
    n_imms <- j %>% 
      left_join(newspp, by = c("plot", "year", "date0")) %>%
      summarise(n = unique(n)) %>%
      pull(n)
    
    #loop through communities, increasing/decreasing one individual at a time until diff is reached
    #then summarise by OTU
    #this looping approach is necessary because it ensures no taxon drops below 0 abundance...
    for(ii in 1:diffs) {
      up <- sample(as.character(j$otu), 1, prob = j$abun)
      down <- sample(as.character(j$otu), 1, prob = j$abun)
      j <- bind_rows(j,
                     mutate(j[1, ], otu = up, abun = 1),
                     mutate(j[1, ], otu = down, abun = -1)) %>%
        group_by(plot, year, date0, otu) %>%
        summarise(abun = sum(abun))
    }
    
    #determine subsequent sampling date for the plot
    j$date1 <- j %>% 
      left_join(newspp, by = c("plot", "year", "date0")) %>% 
      summarise(d1 = unique(date1)) %>% 
      pull(d1)
      
    #determine how many immigrants are slated to arrive
    imms <- j %>% 
      left_join(newspp, by = c("plot", "year", "date0")) %>%
      summarise(n = unique(n)) %>%
      pull(n)
    
    #reduce community size by the number of immigrants arriving times their median initial observed abundance, 
    #this is necessary to maintain community sizes of 1000
    j <- j %>%
      do(data.frame(stringsAsFactors = FALSE,
                    plot = .$plot[1],
                    year = .$year[1],
                    date0 = .$date0[1],
                    date1 = .$date1[1],
                    otu = sample(rep(.$otu, times = .$abun), 1000 - imm_abun * n_imms),
                    abun = 1))
    
    #now add soil imms, resulting in community sizes of 1000
    j <- j %>%
      bind_rows(filter(arrivals, plot == j$plot[1] & date0 == j$date0[1])) %>%
      group_by(plot, year, date0, date1, otu) %>%
      summarise(abun = sum(abun))
    
    #save
    comms[[i]] <- j
  
  }  
  
  comms %>% 
    bind_rows() %>% 
    ungroup() %>%
    rename(date = date1) %>%
    select(-date0) %>%
    saveRDS("Outputs\\sim_comms.RDS")

}

comms <- readRDS("Outputs\\sim_comms.RDS")

#summarize neutral expectations and add metadata

nulls0 <- otus %>% 
  left_join(distinct(meta, sample, year, plot, date, source, sp), by = c('sample')) %>%
  filter(source == "phyllosphere") %>%
  group_by(plot, year, sp) %>%
  filter(date == min(date)) %>%
  select(-source)

nulls <- comms %>% 
  left_join(distinct(meta, sample, plot, date, sp), by = c('plot','date')) %>%
  bind_rows(nulls0) %>%
  mutate(data = 'neutral_prediction')


#### plotting ####

#plot illustrating that most phyllosphere taxa are present in the soil.
fig_soil_presence <- function(){}

#determine which phyllosphere taxa (rarefied) are in the soilpool (unrarefied)
j <- otus %>%
  filter(abun > 0) %>%
  left_join(meta[, c('sample','source','year','date')], by = "sample") %>%
  filter(source == 'phyllosphere') %>%
  group_by(year, otu) %>%
  mutate(origin = ifelse(otu %in% soilpool$otu, 'soil', 'other'))

#how many phyllosphere taxa were NEVER found in the soil?
j %>% 
  ungroup() %>% 
  distinct(otu, origin) %>% 
  summarise(
    phyllosphere_otus_never_in_soil_pool = length(otu[origin == 'other']),
    all_phyllosphere_otus = length(otu),
    percent_otus_never_in_soil_pool = phyllosphere_otus_never_in_soil_pool / all_phyllosphere_otus) %>%
  gather(var, val)

#calculate column for calendar year
meta$year_date <- meta$date
meta$year_date[meta$year == 2017] <- meta$date[meta$year == 2017] - 365

#calculate community and richness proportions present in soil
tmp <- j %>%
  left_join(meta[, c('sample','sp','plot','year','year_date','date')], 
            by = c("sample", "year", "date")) %>%
  group_by(sp, plot, year, date, year_date, sample) %>%
  summarise(Abundance = sum(abun[origin == 'soil'])/sum(abun),
            Richness = length(unique(otu[origin == 'soil'])) / length(unique(otu))) %>%
  gather(var, val, Abundance, Richness) %>%
  ungroup() %>%
  mutate(sp = ifelse(sp == 'miscanthus', 'Miscanthus', 'Switchgrass'),
         group = paste(sp, year)) 

#calculate temporary dataframe of means for plotting lines
tmp2 <- tmp %>% 
  group_by(sp, year, year_date, var) %>% 
  summarise(val = mean(val))

#create list of dates used for x-axis breaks
datez1 <- as.Date(c('2016-06-01','2016-08-01','2016-10-01'))
datez2 <- as.Date(c('2017-06-01','2017-08-01','2017-10-01'))

#general plot
p <- ggplot(tmp, aes(x = year_date, y = val, color = sp)) +
  stat_summary(fun.data = "mean_cl_boot", aes(shape = factor(year))) +
  scale_color_manual(values = c('#026400','#A3CD57'), name = '') +
  scale_shape_manual(values = c(16,1), name = '') +
  scale_linetype_discrete(name = '') +
  scale_y_continuous(labels = scales::percent) +
  scale_x_date(breaks = datez1, date_labels = "%b") +
  expand_limits(y = 0) +
  theme_classic() +
  theme(plot.margin = unit(c(0, 0, 0, 0), "pt"), legend.position = 'bottom') +
  labs(x = '')

#panel a
p1 <- p %+% filter(tmp, var == 'Richness') +
  geom_point(data = filter(tmp2, var == 'Richness' & sp == 'Miscanthus'), size = 2) +
  labs(x = 'Calendar day', y = '% OTUs in soil', tag = 'a') +
  geom_line(aes(group = interaction(year, sp)), size = 1, 
            data = filter(tmp2, var == 'Richness')) + 
  theme(legend.position = 'none')

#panel b
p2 <- p %+% filter(tmp, var == 'Abundance') +
  geom_point(data = filter(tmp2, var == 'Abundance' & sp == 'Miscanthus'), size = 2) +
  labs(x = 'Calendar day', y = '% community in soil', tag = 'b') +
  geom_line(aes(group = interaction(year, sp)), 
            data = filter(tmp2, var == 'Abundance')) + 
  guides(shape = guide_legend(order = 2, override.aes = list(color = '#A3CD57')),
         color = guide_legend(order = 1))

#save, then remove legend
myleg <- get_legend(p2)
p2 <- p2 + theme(legend.position = 'none')

#compare abundances/occupancy between cores
fig_abun_biplot <- function(){}

#use unrarefied OTUs? if so, necessary to calculate relative abundances within sample before pooling to account for variable sampling sizes
j <- otus_unrarefied %>%
  group_by(sample) %>%
  mutate(abun = abun / sum(abun)) %>%
  left_join(meta, by = 'sample') %>%
  mutate(group = ifelse(source == 'soil', source, sp)) %>%
  group_by(group, otu)

#calculate relative abundances across all samples for each plant species
#also calculate mean relative abundances for each OTU across both plant species
relabun <- j %>%
  summarise(relabun = sum(abun)) %>%
  group_by(group) %>%
  mutate(relabun = relabun / sum(relabun)) %>%
  ungroup() %>%
  mutate(group = paste(group, 'relabun', sep = '_')) %>%
  spread(group, relabun) %>%
  mutate(
    plant_relabun = (miscanthus_relabun + switchgrass_relabun) / 2,
    mcore = otu %in% cores$OTU[cores$plant_year == 'misc16'],
    score = otu %in% cores$OTU[cores$plant_year != 'misc16']) %>%
  select(-plant_relabun) %>%
  rename(miscanthus = miscanthus_relabun, switchgrass = switchgrass_relabun) %>%
  gather(sp, relabun, miscanthus, switchgrass)

p3 <- ggplot(relabun, aes(x = soil_relabun, y = relabun)) +
  geom_point(alpha = 0.2, size = 1) +
  geom_point(data = filter(relabun, score & sp == 'switchgrass'), shape = 21, fill = '#A3CD57', size = 2) +
  geom_point(data = filter(relabun, mcore & sp == 'miscanthus'), shape = 21, fill = '#026400', size = 2) +
  geom_point(data = filter(relabun, score & mcore), shape = 21, fill = '#A3CD57') +
  geom_point(aes(y = soil_relabun, x = relabun), alpha = 0) + 
  geom_abline(slope = 1, lty = 3) +
  stat_smooth(color = 'red', se = FALSE) +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = 'Rel. abun. in soil', y = 'Rel. abun. in phyllosphere', tag = 'c') +
  theme_classic() + 
  theme(
    legend.position = 'none',
    strip.background = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "pt"),
    plot.title = element_text(hjust = 0.5))

#violin plots comparing observed - predicted OTU richness (null/neutral model)
fig_violin_filtering <- function(){}

#combine observed and null prediction data
j <- otus %>%
  left_join(meta, by = 'sample') %>%
  filter(source == 'phyllosphere' & abun > 0) %>%
  transmute(plot, year, date, otu, abun, sp, data = 'observed') %>%
  bind_rows(nulls) %>%
  filter(abun > 0) %>% 
  group_by(plot, year) 

#create calendar year column called year_date
j$year_date <- j$date
j$year_date[j$year == 2017] <- j$year_date[j$year == 2017] - 365

#calculate richness for each observed and predicted sample
j <- j %>%
  mutate(sp = ifelse(sp == 'switchgrass', 'Switchgrass', 'Miscanthus')) %>%
  group_by(data, plot, year, year_date, date, sp) %>%
  summarise(rich = length(unique(otu))) %>%
  spread(data, rich) %>%
  mutate(group = paste(sp, year))

t0s <- j %>% 
  group_by(year, plot) %>%
  filter(date == min(date))

tmp <- j %>%
  group_by(year, plot) %>%
  filter(date != min(date))

#general plot
p <- ggplot(tmp, aes(y = observed - neutral_prediction, x = date, group = date)) + 
  geom_hline(yintercept = 0, lty = 3) +
  scale_y_continuous(breaks = c(-80,-60,-40,-20,0)) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    legend.position = 'none',
    plot.margin = unit(c(0, 0, 0, 0), "pt"),
    plot.title = element_text(hjust = 0.5))

#panel d
p4 <- p %+% filter(tmp, group == 'Miscanthus 2016') + 
  scale_x_date(breaks = datez1, date_labels = "%b '%y") +
  expand_limits(y = c(-92,16), x = range(j$year_date)) +
  geom_violin(width = 20, fill = '#026400', color = '#026400', alpha = 0.4) + 
  geom_jitter(color = '#026400', width = 2) +
  geom_point(color = '#026400', data = filter(t0s, group == 'Miscanthus 2016'), shape = 4) +
  labs(x = '', y = 'Observed - predicted OTUs', tag = 'd')

#panel b
p5 <- p %+% filter(tmp, group == 'Switchgrass 2016') + 
  scale_x_date(breaks = datez1, date_labels = "%b '%y") +
  expand_limits(y = c(-92,16), x = range(j$year_date)) +
  geom_violin(width = 20, fill = '#A3CD57', color = '#A3CD57', alpha = 0.4) + 
  geom_jitter(color = '#A3CD57', width = 2) +
  geom_point(color = '#A3CD57', data = filter(t0s, group == 'Switchgrass 2016'), shape = 4) +
  labs(x = '', y = '', tag = 'e')

#panel c
p6 <- p %+% filter(tmp, group == 'Switchgrass 2017') + 
  scale_x_date(breaks = datez2, date_labels = "%b '%y") +
  expand_limits(y = c(-92,16), x = range(j$year_date + 365)) +
  geom_violin(width = 20, fill = NA, color = '#A3CD57') + 
  geom_jitter(color = '#A3CD57', width = 2, shape = 1) +
  geom_point(color = '#A3CD57', data = filter(t0s, group == 'Switchgrass 2017'), shape = 4) +
  labs(x = '', y = '', tag = 'f')


#put it together and save PDF figure``
fig_combine <- function(){}

ps <- plot_grid(p1, p2, p3, p4, p5, p6, ncol = 3, align = 'hv')
pss <- plot_grid(ps, myleg, ncol = 1, rel_heights = c(1, .08))
pss
ggsave('Figures\\Fig3_SoilPhyllosphere.pdf', width = 8, height = 5)
