library(tidyverse)
options(scipen = 999)

F1A = read.csv("~/PhD/ARG project/data/F1A_RGI_taxon.tsv", sep = "\t", row.names=1)
F1B = read.csv("~/PhD/ARG project/data/F1B_RGI_taxon.tsv", sep = "\t", row.names=1)
F1C = read.csv("~/PhD/ARG project/data/F1C_RGI_taxon.tsv", sep = "\t", row.names=1)
F1D = read.csv("~/PhD/ARG project/data/F1D_RGI_taxon.tsv", sep = "\t", row.names=1)
F2A = read.csv("~/PhD/ARG project/data/F2A_RGI_taxon.tsv", sep = "\t", row.names=1)
F2B = read.csv("~/PhD/ARG project/data/F2B_RGI_taxon.tsv", sep = "\t", row.names=1)
F2C = read.csv("~/PhD/ARG project/data/F2C_RGI_taxon.tsv", sep = "\t", row.names=1)
F2D = read.csv("~/PhD/ARG project/data/F2D_RGI_taxon.tsv", sep = "\t", row.names=1)
F3A = read.csv("~/PhD/ARG project/data/F3A_RGI_taxon.tsv", sep = "\t", row.names=1)
F3B = read.csv("~/PhD/ARG project/data/F3B_RGI_taxon.tsv", sep = "\t", row.names=1)
F3C = read.csv("~/PhD/ARG project/data/F3C_RGI_taxon.tsv", sep = "\t", row.names=1)
F3D = read.csv("~/PhD/ARG project/data/F3D_RGI_taxon.tsv", sep = "\t", row.names=1)
F4A = read.csv("~/PhD/ARG project/data/F4A_RGI_taxon.tsv", sep = "\t", row.names=1)
F4B = read.csv("~/PhD/ARG project/data/F4B_RGI_taxon.tsv", sep = "\t", row.names=1)
F4C = read.csv("~/PhD/ARG project/data/F4C_RGI_taxon.tsv", sep = "\t", row.names=1)
F4D = read.csv("~/PhD/ARG project/data/F4D_RGI_taxon.tsv", sep = "\t", row.names=1)
F6A = read.csv("~/PhD/ARG project/data/F6A_RGI_taxon.tsv", sep = "\t", row.names=1)
F6B = read.csv("~/PhD/ARG project/data/F6B_RGI_taxon.tsv", sep = "\t", row.names=1)
F6C = read.csv("~/PhD/ARG project/data/F6C_RGI_taxon.tsv", sep = "\t", row.names=1)
F6D = read.csv("~/PhD/ARG project/data/F6D_RGI_taxon.tsv", sep = "\t", row.names=1)
F8A = read.csv("~/PhD/ARG project/data/F8A_RGI_taxon.tsv", sep = "\t", row.names=1)
F8B = read.csv("~/PhD/ARG project/data/F8B_RGI_taxon.tsv", sep = "\t", row.names=1)
F8C = read.csv("~/PhD/ARG project/data/F8C_RGI_taxon.tsv", sep = "\t", row.names=1)
F8D = read.csv("~/PhD/ARG project/data/F8D_RGI_taxon.tsv", sep = "\t", row.names=1)

F1A$pond = "1A"
F1B$pond = "1B"
F1C$pond = "1C"
F1D$pond = "1D"
F2A$pond = "2A"
F2B$pond = "2B"
F2C$pond = "2C"
F2D$pond = "2D"
F3A$pond = "3A"
F3B$pond = "3B"
F3C$pond = "3C"
F3D$pond = "3D"
F4A$pond = "4A"
F4B$pond = "4B"
F4C$pond = "4C"
F4D$pond = "4D"
F6A$pond = "6A"
F6B$pond = "6B"
F6C$pond = "6C"
F6D$pond = "6D"
F8A$pond = "8A"
F8B$pond = "8B"
F8C$pond = "8C"
F8D$pond = "8D"

df = rbind(F1A,F1B,F1C,F1D,
           F2A,F2B,F2C,F2D,
           F3A,F3B,F3C,F3D,
           F4A,F4B,F4C,F4D,
           F6A,F6B,F6C,F6D,
           F8A,F8B,F8C,F8D)

#remove nudge genes :( loss from 30k to 1k
taxonomy = df[df$Nudged != "True", ]

# define the levels of taxa in "Taxonomy"
taxa = c("superkingdom","phylum","class","order","family","genus","species","strain")

# seperate the 7 levels of taxa + strain from taxonomy as new columns
ponds = taxonomy %>% separate(Taxonomy, taxa, "; ")

#remove strain column
ponds$strain = NULL
taxa = c("superkingdom","phylum","class","order","family","genus","species")

# remove any "NA" or blank spaces in taxa columns with "Unclassified"
ponds[taxa] = lapply(ponds[taxa], gsub, pattern = "NA|^$", replacement = "Unclassified")

# change <NA> to "Unclassified"
ponds[taxa][is.na(ponds[taxa])] = "Unclassified"

# if a gene conferes multiple resistances, split into new row copies with each resistance
ponds = ponds %>% mutate(Drug.Class = strsplit(as.character(Drug.Class), "; ")) %>% unnest(Drug.Class)

# remove the word "antibiotic" and any spaces from the Drug.Class column
ponds$Drug.Class = str_trim(gsub("antibiotic", "", ponds$Drug.Class), side = c("both"))

#pluralise Drug.Class and add specific used examples from social dataset to drug class
ponds$Drug.Class = ponds$Drug.Class %>% str_c("s") %>% 
  str_replace_all("penams", "penams (e.g. amoxicillin)") %>%
  str_replace_all("macrolides", "macrolides (e.g. erythromycin)") %>%
  str_replace_all("diaminopyrimidines", "diaminopyrimidines (e.g. trimethoprim)") %>%
  str_replace_all("tetracyclines", "tetracyclines (e.g. chlorotetracycline,\ndoxycycline, oxytetracycline)") %>%
  str_replace_all("fluroquinolones", "fluroquinolones (e.g. enrofloxacin, ciprofloxacin)") %>%
  str_replace_all("aminoglycosides", "aminoglycosides (e.g. neomycin)") %>%
  str_replace_all("sulfonamides", "sulfonamides (e.g. sulfadiazine,\nsulfamethoxazole)")

# rename Drug.Class colname as antibiotic class
ponds = ponds %>% rename("Antibiotic Class" = "Drug.Class")

# assign colours to drug classes
clrs = c('#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#42d4f4', '#f032e6', '#bfef45', '#fabed4', '#469990', '#dcbeff', '#800000', '#aaffc3', '#808000', '#000075', '#a9a9a9', '#000000')
names(clrs) = unique(ponds$`Antibiotic Class`)
colScale <- scale_fill_manual(name = "Antibiotic Class",values = clrs)

# trim dataset of useful information and replace genenames with shorten version
ARGs = ponds[c("Antibiotic Class","ARO","pond", "Best_Hit_ARO", "RPKM")]
ARGs$Best_Hit_ARO = gsub("Mycobacterium tuberculosis folC with mutation conferring resistance to para-aminosalicylic acid", "folc", ARGs$Best_Hit_ARO)
ARGs$Best_Hit_ARO = gsub("Acinetobacter baumannii AbaQ", "AbaQ", ARGs$Best_Hit_ARO)
ARGs$Best_Hit_ARO = gsub("Mycobacterium tuberculosis intrinsic murA conferring resistance to fosfomycin", "murA", ARGs$Best_Hit_ARO)
ARGs$Best_Hit_ARO = gsub("Mycobacterium tuberculosis rpsL mutations conferring resistance to Streptomycin", "rpsL", ARGs$Best_Hit_ARO)

# plot gene diversity
figure = ggplot(ARGs, aes(fill=`Antibiotic Class`, y=Best_Hit_ARO, x=pond)) +
  geom_tile() +
  theme_bw() +
  colScale +
  ylab("Name of gene conferring resistance") +
  xlab("Farm")

#save figure
ggsave(filename = "~/PhD/ARG project/figures/Figure_3.pdf", plot = figure)
