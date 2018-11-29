library(tidyr)
library(dplyr)
library(RColorBrewer)
library(knitr)
library(ggplot2)

#Setting up and reading out files for CheckM data, and gtdbtk taxonomy assignments

```{r read, warning=FALSE, message=FALSE}
arc_class <- read.table("~/MICB405_Project/gtdbtk.ar122.classification_pplacer.tsv", sep="\t")
bac_class <- read.table("~/MICB405_Project/gtdbtk.bac120.classification_pplacer.tsv", sep="\t")
gtdb_dat <- rbind(arc_class, bac_class) %>% 
  dplyr::rename(mag = V1) %>% 
  separate(V2, sep=';', into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")) %>%
  dplyr::select(mag, Kingdom, Phylum, Class, Order, Family)
checkm_dat <- read.table("~//MICB405_Project/MetaBAT2_SaanichInlet_10m_min1500_checkM_stdout.tsv",
                         header=TRUE,
                         sep="\t",
                         comment.char = '') %>% 
  dplyr::rename(mag = Bin.Id) %>% 
  dplyr::select(mag, Completeness, Contamination)
# Due to a bug in the renaming script we have to rename the bins. Its a bit hacky but works using tidyverse functions
metag_rpkm <- read.table("~/MICB405_Project/SaanichInlet_10m_binned.rpkm.csv", header=T, sep=',') %>% 
  mutate(Sequence = gsub('m_', 'm.', Sequence)) %>% 
  mutate(Sequence = gsub('Inlet_', 'Inlet.', Sequence)) %>% 
  separate(col=Sequence, into=c("mag", "contig"), sep='_', extra="merge") %>% 
  group_by(Sample, mag) %>% 
  summarise(g_rpkm = sum(RPKM)) %>% 
  mutate(mag = gsub('Inlet.', 'Inlet_', mag))
```
#Joining the metagenome RPKM, checkM, and GTDB-tk data frames 
#Then take the mean RPKM of each MAG across the different cruises --> Specific for the quality and abundance bubble plots

```{r }
rpkm_dat <- left_join(metag_rpkm, checkm_dat, by="mag") %>% 
  left_join(gtdb_dat, by="mag") %>% 
  group_by(mag, Kingdom, Phylum, Class, Completeness, Contamination) %>% 
  summarise(g_rpkm = mean(g_rpkm))

#Joining files and filtering names
rpkm_dat <- left_join(metag_rpkm, checkm_dat, by="mag") %>%
  left_join(gtdb_dat, by="mag") %>% 
  left_join(ko_F, by="mag")
rpkm_dat <- mutate(rpkm_dat, Family=substring(Family, 4)) 
names(rpkm_dat)[names(rpkm_dat)=="g_rpkm"] <- "Total_RPKM" 
names(rpkm_dat)[names(rpkm_dat)=="Sample"] <- "Cruise" 
rpkm_dat <- mutate(rpkm_dat, Cruise=substring(Cruise,0,5)) 

#Plotting in 'ggplot' --> Creating a plot of Contamination vs Completion
  ggplot(rpkm_dat, aes(x=Completeness, y=Contamination, col=Family)) +
  geom_jitter(aes(size=Total_RPKM)) +
  scale_size(range=c(1,10)) +
  xlim(c(50,100)) +
  ylim(c(0,100)) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_line(colour = "#bdbdbd", linetype = "dotted"),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(face="bold", angle = 90), 
        axis.title.x = element_text(face = "bold"), 
        axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = ("bottom")) + geom_hline(yintercept = 5, linetype = 2) + geom_hline(yintercept = 10, linetype = 2) + geom_hline(yintercept = 15, linetype = 2) +
    geom_vline(xintercept = 70, linetype = 2) + geom_vline(xintercept = 90, linetype = 2) +
  ylab("Contamination (%)") +
  xlab("Completion(%)")

```

#Plotting in 'ggplot' continued --> Determing the relative abundance of MAGs across the different cruises and time points for Orders involved in phothosynthetic pathways
names(rpkm_dat)[names(rpkm_dat)=="Cruise"] <- "Sample"
rpkm_dat <- filter(rpkm_dat, Completeness > 50 & Contamination < 10) # Filtering out low quality MAGs
rpkm_dat <- filter(rpkm_dat, Order %in% c("o__Synechococcales_A", "o__Rhodospirillales", "o__Rhizobiales")) #Filtering for Orders of Bacteria that are involved in photosynthetic pathways
rpkm_dat <- mutate(rpkm_dat, Order=substring(Order, 4))
names(rpkm_dat)[names(rpkm_dat)=="Sample"] <- "Cruise"
rpkm_dat <- mutate(rpkm_dat, mag=substring(mag, 18))

ggplot(rpkm_dat, aes(x=Cruise, y=mag, col=Order, shape=Order)) +
  geom_jitter(aes(size=Total_RPKM, color= Order)) +
  scale_size(range=c(1,10)) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_line(colour = "#bdbdbd", linetype = "dotted"),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(face="bold", angle = 90), 
        axis.title.x = element_text(face = "bold"), 
        axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = ("bottom")) +
  xlab("Cruise") +
  ylab("Metagenome Assembled Genome")
