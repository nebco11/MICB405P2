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
CF_Genes <- c("mdh", "IDH1, IDH2, icd", "korB", "sdhA, frdA", "sdhB, frdB", "sdhC, frdC", "sdhD, frdD", "metF", "atoB", "ackA", "ppdK", "ppsA", "folD", "ppc", "fumA, fumB", "fumB", "fumC", "acnA", "acnB", "MUT", "mcmA1", "acs", "sucD", "sucC", "fhs", "pyc", "pycA", "pycB", "accC", "accA", "accD", "accB, bccP", "epi", "pta", "frdB")


#Load KO data and merge with Prokka MAG Map
KOs <- read.table("~/Downloads/KOs.csv") %>% dplyr::rename(zname = V1) %>% dplyr::rename(ko = V2) 
KOs <- KOs[-c(1),] #remove first row of headers
ko_F <- merge(KOs, prokka_mag_map, by.x="zname")


#Generate Carbon Fixation Cruise Bubbles

#Step 1 - make rpkm_dat_all table

rpkm_dat_all <- left_join(metag_rpkm, checkm_dat, by="mag") %>% 
  left_join(gtdb_dat_bub, by="mag") %>% 
  left_join(ko_F, by="mag") %>%
  filter(Completeness > 50 & Contamination < 10) %>% 
  filter(ko %in% c("K00024", "K00031",
         "K00169",
         "K00170",
         "K00171",
         "K00172",
         "K00174",
         "K00175",
         "K00176",
         "K00177",
         "K00194",
         "K00196",
         "K00197",
         "K00198",
         "K00239",
         "K00240",
         "K00241",
         "K00242",
         "K00244",
         "K00245",
         "K00246",
         "K00247",
         "K00297",
         "K00625",
         "K00626",
         "K00925",
         "K01006",
         "K01007",
         "K01491",
         "K01595",
         "K01648",
         "K01676",
         "K01677",
         "K01678",
         "K01679",
         "K01681",
         "K01682",
         "K01847",
         "K01848",
         "K01849",
         "K01895",
         "K01902",
         "K01903",
         "K01938",
         "K01958",
         "K01959",
         "K01960",
         "K01961",
         "K01962",
         "K01963",
         "K01964",
         "K02160",
         "K03737",
         "K05299",
         "K05606",
         "K08691",
         "K09709",
         "K13788",
         "K14138",
         "K14449",
         "K14465",
         "K14466",
         "K14467",
         "K14468",
         "K14469",
         "K14470",
         "K14471",
         "K14472",
         "K14534",
         "K15016",
         "K15017",
         "K15018",
         "K15019",
         "K15020",
         "K15022",
         "K15023",
         "K15024",
         "K15036",
         "K15037",
         "K15038",
         "K15039",
         "K15052",
         "K15230",
         "K15231",
         "K15232",
         "K15233",
         "K15234",
         "K18209",
         "K18210",
         "K18556",
         "K18557",
         "K18558",
         "K18559",
         "K18560",
         "K18593",
         "K18594",
         "K18602",
         "K18603",
         "K18604",
         "K18605",
         "K18859",
         "K18860", "K18861")) #filter for all CF genes

#step 2: Run the following individually in order listed
  
  rpkm_dat_all <- rpkm_dat_all[-c(1147:2292),] #Change the numbers in [-c(a:b)] depending on where SI073 samples start appearing (will need to go through the table)
  rpkm_dat_all <- mutate(rpkm_dat_all, mag=substring(mag, 18)) #Keep only MAG number
  names(rpkm_dat_all)[names(rpkm_dat_all)=="g_rpkm"] <- "Total_RPKM"#Rename g_rpkm to Total_RPKM
  names(rpkm_dat_all)[names(rpkm_dat_all)=="mag"] <- "MAG" #Rename mag to MAG
  names(rpkm_dat_all)[names(rpkm_dat_all)=="ko"] <- "KO" #Rename ko to KO
  names(rpkm_dat_all)[names(rpkm_dat_all)=="Sample"] <- "Cruise" #Rename Sample to Cruise
  rpkm_dat_all <- mutate(rpkm_dat_all, Cruise=substring(Cruise,0,5)) #Remove prefix of Cruise
  rpkm_dat_all <- mutate(rpkm_dat_all, Family=substring(Family,4)) #Remove prefix of Family

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
        ylab("Prokaryotic Carbon Fixation Gene") +
        xlab("Metagenome Assembled Genome") +
       scale_y_discrete(labels=CF_Genes) #Hide line 84 in initial plot to identify KO values to find gene names in KEGG.(see Line 27)
