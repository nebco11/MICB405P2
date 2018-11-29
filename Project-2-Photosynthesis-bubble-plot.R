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
Photosynthesis_Genes <- c("atpB",	"atpF",	"atpE",	"atpA",	"atpD",	"atpH",	"atpC",	"atpG",	"petC")


#Load KO data and merge with Prokka MAG Map
KOs <- read_csv("~/Downloads/KOs.csv")
ko_F <- merge(KOs, prokka_mag_map, by.x="zname")


#Generate Photosynthesis Cruise Bubbles

#Step 1 - make rpkm_dat_all table

rpkm_dat_all <- left_join(metag_rpkm, checkm_dat, by="mag") %>% 
  left_join(gtdb_dat_bub, by="mag") %>% 
  left_join(ko_F, by="mag") %>%
  filter(Completeness > 50 & Contamination < 10) %>% 
  filter(ko %in% c("K02108", "K02109", "K02110", "K02111", "K02112", "K02113", "K02114", 
                   "K02115","K02634", "K02635", "K02636", "K02637", "K02638", "K02639", 
                   "K02640", "K02641", "K02642", "K02643", "K02689", "K02690", "K02691", 
                   "K02692", "K02693", "K02694", "K02695", "K02696", "K02697", "K02698", 
                   "K02699", "K02700", "K02701", "K02702", "K02703", "K02704", "K02705", 
                   "K02706", "K02707", "K02708", "K02709", "K02710", "K02711", "K02712", 
                   "K02713", "K02714","K02716", "K02717", "K02718", "K027119", "K02720", 
                   "K02721", "K02722", "K02723", "K02724", "K03541", "K03542", "K03689", 
                   "K08901", "K08902", "K08903", "K08904","K08905", "K08906", "K14332" )) #filter for all photosynthesis genes
 
#step 2: Run the following individually in order listed
  
  rpkm_dat_all <- rpkm_dat_all[-c(340:678),] #Change the numbers in [-c(a:b)] depending on where SI073 samples start appearing (will need to go through the table)
  rpkm_dat_all <- mutate(rpkm_dat_all, mag=substring(mag, 18)) #Keep only MAG number
  names(rpkm_dat_all)[names(rpkm_dat_all)=="g_rpkm"] <- "Total_RPKM"#Rename g_rpkm to Total_RPKM
  names(rpkm_dat_all)[names(rpkm_dat_all)=="mag"] <- "MAG" #Rename mag to MAG
  names(rpkm_dat_all)[names(rpkm_dat_all)=="ko"] <- "KO" #Rename ko to KO
  names(rpkm_dat_all)[names(rpkm_dat_all)=="Sample"] <- "Cruise" #Rename Sample to Cruise
  rpkm_dat_all <- mutate(rpkm_dat_all, Cruise=substring(Cruise,0,5)) #Remove prefix of Cruise 

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
        scale_color_hue(labels = c("Unknown Family", "BACL11", "Crocinitomicaceae", 
                                   "Cryomorphaceae", "D2472", "Flavobacteriaceae", "Halieaceae", "HIMB59", 
                                   "HTCC2089", "Maricaulaceae", "Microbacteriaceae", "NAC60-12", "Pelagibacteraceae", 
                                   "Planctomycetaceae", "Porticoccaceae", "Pseudohongiellaceae", "Puniceispirillaceae", 
                                   "Rhodobacteraceae", "RS24", "SG8-40", "Thioglobaceae", "TMED25", "UA16", "UBA10009", 
                                   "UBA7434", "UBC8229", "UBA9320")) + #Hide lines 76-81 in initial plot; change name of labels depending on what species appear on initial plot
        ylab("Metabolic Genes") +
        xlab("Metagenome Assembled Genomes") +
        scale_y_discrete(labels=Photosynthesis_Genes) #Hide line 84 in initial plot to identify KO values to find gene names in KEGG.(see Line 27)



