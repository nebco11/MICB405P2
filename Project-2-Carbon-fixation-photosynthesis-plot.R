#Connor's initial script:

```{r read, warning=FALSE, message=FALSE}

arc_class <- read.table("~/MICB405_Project/gtdbtk.ar122.classification_pplacer.tsv", sep="\t")
bac_class <- read.table("~/MICB405_Project/gtdbtk.bac120.classification_pplacer.tsv", sep="\t")
gtdb_dat_bub <- rbind(arc_class, bac_class) %>% 
  dplyr::rename(mag = V1) %>% 
  separate(V2, sep=';', into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")) %>%
  dplyr::select(mag, Kingdom, Phylum, Class, Order, Family)
checkm_dat <- read.table("~/Desktop/Project 2/MetaBAT2_SaanichInlet_10m_min1500_checkM_stdout.tsv",
                         header=TRUE,
                         sep="\t",
                         comment.char = '') %>% 
  dplyr::rename(mag = Bin.Id) %>% 
  dplyr::select(mag, Completeness, Contamination)
# Due to a bug in the renaming script we have to rename the bins. Its a bit hacky but works using tidyverse functions
metag_rpkm <- read.table("~/Desktop/Project 2/SaanichInlet_10m_binned.rpkm.csv", header=T, sep=',') %>% 
  mutate(Sequence = gsub('m_', 'm.', Sequence)) %>% 
  mutate(Sequence = gsub('Inlet_', 'Inlet.', Sequence)) %>% 
  separate(col=Sequence, into=c("mag", "contig"), sep='_', extra="merge") %>% 
  group_by(Sample, mag) %>% 
  summarise(g_rpkm = sum(RPKM)) %>% 
  mutate(mag = gsub('Inlet.', 'Inlet_', mag))
```

#Assigning Gene Names for Bubble Plot (get gene names from first ggplot without gene names) (see line 84)
Photosynthesis_CF_Genes <- c("mdh", "mdh1", "maeB", "E1.1.1.82", "gapA", 
                            "gap2", "tktA, tktB", "pgk", "ppdk", "ppc", "ALDO", "FBA", "rpe", "TPI", "rpiA", "rpiB", "glpX", "fbp", "glpX-SEBP", "fbaB")


#Load KO data and merge with Prokka MAG Map
KOs <- read_csv("~/Downloads/KOs.csv") %>% dplyr::rename(zname = V1) %>% dplyr::rename(ko = V2)
KOs <- KOs[-c(1),] #remove first row of headers
ko_F <- merge(KOs, prokka_mag_map, by.x="zname") 

#Generate Carbon fixation in Photosynthetic organism Cruise Bubbles

#Step 1 - make rpkm_dat_all table

rpkm_dat_all <- left_join(metag_rpkm, checkm_dat, by="mag") %>% 
  left_join(gtdb_dat_bub, by="mag") %>% 
  left_join(ko_F, by="mag") %>%
  filter(Completeness > 50 & Contamination < 10) %>% 
  filter(ko %in% c("K00024", "K00025", "K00026", "K00028", "K00029", "K00051", "K00134", 
                   "K00150","K00615", "K00814", "K00855", "K00927", "K01006", "K01086", 
                   "K01100", "K01595", "K01601", "K01621", "K01623", "K01624", "K01783", 
                   "K01803", "K01807", "K01808", "K02446", "K03841", "K04041", "K05298", 
                   "K11214", "K11532", "K11645", "K14272", "K014454", "K14455")) #filter for all CF photosynthesis genes
 
#step 2: Run the following individually in order listed
  
  rpkm_dat_all <- rpkm_dat_all[-c(523:1044),] #Change the numbers in [-c(a:b)] depending on where SI073 samples start appearing (will need to go through the table)
  rpkm_dat_all <- mutate(rpkm_dat_all, mag=substring(mag, 18)) #Keep only MAG number
  names(rpkm_dat_all)[names(rpkm_dat_all)=="g_rpkm"] <- "Total_RPKM"#Rename g_rpkm to Total_RPKM
  names(rpkm_dat_all)[names(rpkm_dat_all)=="mag"] <- "MAG" #Rename mag to MAG
  names(rpkm_dat_all)[names(rpkm_dat_all)=="ko"] <- "KO" #Rename ko to KO
  names(rpkm_dat_all)[names(rpkm_dat_all)=="Sample"] <- "Cruise" #Rename Sample to Cruise
  rpkm_dat_all <- mutate(rpkm_dat_all, Cruise=substring(Cruise,0,5)) #Remove prefix of Cruise 
  rpkm_dat_all <- mutate(rpkm_dat_all, Family=substring(Family,4)) #remove prefix of Family

#step 3: Generate bubble plot
  
ggplot(data = rpkm_dat_all, aes(x=MAG, y=KO)) +
geom_jitter(aes(color= Family, size= Total_RPKM, shape = Cruise)) + #use Jitter to avoid overlapping data points on plot
theme(panel.background = element_blank(),
        panel.grid.major = element_line(colour = "#bdbdbd", linetype = "dotted"),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(face="bold", angle = 90),
        axis.title.x = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45,
                                   hjust = 1),
        legend.position = ("bottom")) +
        ylab("Carbon Fixation in Photosynthesis Organism Gene") +
        xlab("Metagenome Assembled Genome") +
        scale_y_discrete(labels=Photosynthesis_CF_Genes) #Hide line 84 in initial plot to identify KO values to find gene names in KEGG.(see Line 27)



