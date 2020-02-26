####Evaluation of NHS-Acetate and DEPC labelling for determination of solvent accessible amino acid residues in protein complexes####
#This R script was used to calculate a normalised intensity for each modified site.
#For this, MaxQuant output files containing information on the individual modifications 
#(e.g. acetylation of lysine or acetylation of serine, threonine, tyrosine) were combined 
#and filtered for the protein of interest. 
#Modified sites with a localisation probability < 0.75 and a peptide score < 80 were discarded. 
#The MaxQuant intensity of each modified site was then normalised 
#by applying a normalisation factor calculated from the sum of all corresponding peptides (i.e. all modified and unmodified peptides) containing the respective site.
#The following equation was used: ((Intensity of modified residue)*100)/(sum corresponding peptides (modified+unmodified)).
#Subsequently, mean and standard error were calculated.


# Install if not installed

#install.packages("openxlsx")
#install.packages("tidyverse")
#install.packages("ggforce")
#install.packages("reshape")
#install.packages('Rmisc', dependencies = TRUE)

# Required libraries

library("openxlsx")
library("tidyverse")
library("ggforce")
library("reshape")
library("Rmisc")

#### 1. Import data ####

# ModificationSpecificPeptides.txt (MaxQuant result file) contains information on peptides
peptides = read.table("modificationSpecificPeptides.txt", sep="\t", header = TRUE, stringsAsFactors = F)

# Peptides file (MaxQuant) is used for info of start and end of peptide sequence
peptides_start_end = read.table("peptides.txt", sep="\t", header = TRUE, stringsAsFactors = F)

# Modification files from Max Quant (identified modified sites)
mod_sites_1 = read.table("mod_name_1.txt", sep="\t", header = TRUE, stringsAsFactors = F)
mod_sites_2 = read.table("mod_name_2.txt", sep="\t", header = TRUE, stringsAsFactors = F)
mod_sites_3 = read.table("mod_name_3.txt", sep="\t", header = TRUE, stringsAsFactors = F)
mod_sites_4 = read.table("mod_name_4.txt", sep="\t", header = TRUE, stringsAsFactors = F)

# File Frac_ASA_protein_structure contains info of GETAREA output, "Ratio." is the relative solvent accessibility of the amino acid
Frac_asa <- read.xlsx("Frac_ASA_Protein_structure.xlsx")

#### 1.1 Add start and end position of every peptide to peptides table ####
peptides = select( my_data_info, Start.position, End.position, Proteins, Sequence)
peptides = merge(peptides, peptides_start_end, by=c("Sequence","Proteins"))

#### 1.2 Remove all contaminant entries and "no Name" entries ####
peptides = peptides[-grep("^CON", peptides$Proteins),]
peptides = peptides[peptides$Proteins != "",]

#### 1.3 Change structure of peptide table ####
#Name: Experiment_1 can be changed in experimental table within MaxQuant search

peptides = select(peptides, Start.position, End.position, Proteins, Intensity.Experiment_1_R1:Intensity.Experiment_5_R3)
peptides_new_structure <- melt(peptides, id=(c("Start.position", "End.position","Proteins")))

#### 2. Filter modified sites ####
#### 2.1 Add column with specific modification information ####
mod_sites_1$Modification = "mod_name_1"
mod_sites_2$Modification = "mod_name_2"
mod_sites_3$Modification = "mod_name_3"
mod_sites_4$Modification = "mod_name_4"

### 2.2 Merge modification files into one ####
mod_sites = merge(mod_sites_1,mod_sites_2, all=TRUE)
mod_sites = merge(mod_sites,mod_sites_3, all=TRUE)
mod_sites = merge(mod_sites,mod_sites_4, all=TRUE)

#### 2.3 Select rows without contaminats ####
mod_sites = mod_sites[-grep("^CON", mod_sites$Proteins),]
mod_sites = mod_sites2[mod_sites2$Proteins != "",]

#### 2.4 Select rows according to score and localisation probability ####
# Score selection can be adjusted
mod_sites = mod_sites[mod_sites$Score > 80,]

# Localization pobability filter can be adjusted
mod_sites = mod_sites[mod_sites$Localization.prob > 0.75,]

#### 2.5 Change structure of modification table ####
# Selection of protein, position, amino acid, modification and intensity of Experiment 1 to Experiment 5 of all replicates
mod_sites = select(mod_sites, Protein, Position, Amino.acid, Modification, Intensity.ADH_DEPC_0_0_mM_R1:Intensity.ADH_DEPC_5_0_mM_R3)
mod_sites_new_structure <- melt(mod_sites, id=(c("Protein","Position","Amino.acid","Modification")))

#### 3. Normalisation of Intensities ####
#### 3.1 Sum intensities of (unmodified+modifified) peptides of every modified site ####
#Pept: peptides in peptides_new_structure table that include the modification site position, the same protein and replicate number as in the mod-sites_new_structure table are selected

for(i in 1:nrow(mod_sites_new_structure)){
  mod_pos = mod_sites_new_structure$Position[i]
  mod_pos = as.integer(mod_pos)
  repl    = mod_sites_new_structure$variable[i]
  protein = mod_sites_new_structure$Protein[i]
  pept    = which(peptides_new_structure$variable == repl & as.integer(peptides_structure$Start.position) <= mod_pos & as.integer(peptides_new_structure$End.position) >= mod_pos & peptides_new_structure$Proteins == protein)
  mod_sites_new_structure$Sum_cor_pept[i] = sum(peptides_new_structure$value[pept], na.rm = TRUE)
}

#### 3.2 Build ratio modified / (unmodified+modified) ####
mod_sites_ratio = mod_sites_new_structure
mod_sites_ratio$Ratio <- (mod_sites_ratio$value *100)/(mod_sites_ratio$Sum_cor_pept)

#### 4. Add new columns with additional information ####
#### 4.1 Concentration (molar excess labelling reagent) ####
mod_sites_ratio$Concentration <- 0 #Experiment 1
mod_sites_ratio$Concentration[mod_sites_ratio$variable %in% c("Intensity.Experiment_2_R1",
                                                            "Intensity.Experiment_2_R2",
                                                            "Intensity.Experiment_2_R3")] <- 2 # value of molar excess
mod_sites_ratio$Concentration[mod_sites_ratio$variable %in% c("Intensity.Experiment_3_R1",
                                                            "Intensity.Experiment_3_R2",
                                                            "Intensity.Experiment_3_R3")] <- 10
mod_sites_ratio$Concentration[mod_sites_ratio$variable %in% c("Intensity.Experiment_4_R1",
                                                            "Intensity.Experiment_4_R2",
                                                            "Intensity.Experiment_4_R3")] <- 50
mod_sites_ratio$Concentration[mod_sites_ratio$variable %in% c("Intensity.Experiment_5_R1",
                                                            "Intensity.Experiment_5_R2",
                                                            "Intensity.Experiment_5_R3")] <- 100

#### 4.2 Replicate information ####
mod_sites_ratio$Replicate = mod_sites_melt$variable
mod_sites_ratio$Replicate = gsub("Intensity.Experiment_\d{1}_R", "",mod_sites_melt$Replicate)

#### 4.3 Add column with Frac_ASA info (obtained from GETAREA) and save file ####
Frac_asa$Position=Frac_asa$Residue
mod_sites_ratio_final=merge(mod_sites_ratio,Frac_asa, by = "Position", all.x= TRUE)

write.xlsx(mod_sites_ratio_final, "Name_of_new_file.xlsx")

#### 5.Calculate mean values and standard error for every concentration ####

mod_sites_ratio_final <- as_data_frame(mod_sites_ratio_final)
mod_sites_ratio_final <- summarySE(mod_sites_ratio_final, measurevar="Ratio", groupvars=c("Position","Amino.acid","Concentration","Modification","Sidechain","Ratio.","In.Out"))

#### 6. Plot results for every amino acid ####

# Create colour gradient 0-20 buried, above 50 solvent accessible
custom <- colorRampPalette(c("red4","red3","orangered3","orange1","orange","cyan3","deepskyblue2","dodgerblue1","dodgerblue3","blue","blue3"))

# Plot
plot_SE = ggplot(mod_sites_ratio_final, aes(x= as.numeric(Concentration), y=Ratio, order= Position, colour=Ratio., shape= Modification)) + 
  geom_errorbar(aes(ymin=Ratio-se, ymax=Ratio+se),size=1.5, width=.1) +
  geom_point(size= 4)+ 
  scale_y_continuous(name = "Percent labelled peptide")+
  scale_x_continuous(name = "Molar excess labelling reagent")+
  facet_wrap_paginate(~ Position + Modification, nrow = 9, ncol =7, scales ="free_x", page=1)+
  theme_grey(base_size = 30)+
  theme(panel.spacing = unit(2, "lines"))+
  scale_colour_gradientn(name = element_blank(),
                         breaks = c(20,50),
                         labels = c("Buried","Solvent accessible"),
                         colours = custom(100),
                         limits =c (0,100),
                         na.value = "grey50")+
  ylim (0,100)+ 
  theme(strip.background = element_blank(),
    strip.text.x = element_blank())+
  scale_shape_manual(values=c(15,19))+
  geom_text(
    size    = 9,
    data    = id1,
    mapping = aes(x = -Inf, y = Inf, label = Label),
    hjust   = -0.2,
    vjust   = 1.5,
    col     = "black")

x11()
plot_SE


