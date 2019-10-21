#Normalization of Labelling data (obtained with MaxQuant) using Peptide Intensitie sums
# for Normalisation of Labelling Intensities of individual labelled Residues

# clean up everything
rm(list=ls())
gc()

# Installing if required
install.packages("readr")
install.packages("tidyverse")
install.packages('Rmisc', dependencies = TRUE)
install.packages("dplyr")
install.packages("naniar")
install.packages("lme4")
install.packages("cowplot")
install.packages("openxlsx")



#####1. Loading required packages###
library("readr")
library("tidyverse")
library("ggplot2")
library("reshape")
library('Rmisc')
library("dplyr")
library("naniar")
library("lme4")
library("cowplot")
library("openxlsx")

####2. Import Peptide MaxQuant (change in excel to csv file)####
#other option 
#read txt from peptides
#my_data <- read_tsv(file.choose())

my_data <- read.csv("peptides.csv", header = T, stringsAsFactors = F)
#remove all contaminant entries and "no Name" entries
my_data = my_data[-grep("^CON", my_data$Proteins),]
my_data = my_data[my_data$Proteins != "",]

#calculate Intensity sum of each experiment
sum_0000_R1 = sum(my_data$`Intensity.ADH_NHS_Ac_0_R1`)
sum_0000_R2 = sum(my_data$`Intensity.ADH_NHS_Ac_0_R2`)
sum_0000_R3 = sum(my_data$`Intensity.ADH_NHS_Ac_0_R3`)

sum_0050_R1 = sum(my_data$`Intensity.ADH_NHS_Ac_0050_R1`)
sum_0050_R2 = sum(my_data$`Intensity.ADH_NHS_Ac_0050_R2`)
sum_0050_R3 = sum(my_data$`Intensity.ADH_NHS_Ac_0050_R3`)

sum_0100_R1 = sum(my_data$`Intensity.ADH_NHS_Ac_0100_R1`)
sum_0100_R2 = sum(my_data$`Intensity.ADH_NHS_Ac_0100_R2`)
sum_0100_R3 = sum(my_data$`Intensity.ADH_NHS_Ac_0100_R3`)

sum_0250_R1 = sum(my_data$`Intensity.ADH_NHS_Ac_0250_R1`)
sum_0250_R2 = sum(my_data$`Intensity.ADH_NHS_Ac_0250_R2`)
sum_0250_R3 = sum(my_data$`Intensity.ADH_NHS_Ac_0250_R3`)

sum_0500_R1 = sum(my_data$`Intensity.ADH_NHS_Ac_0500_R1`)
sum_0500_R2 = sum(my_data$`Intensity.ADH_NHS_Ac_0500_R2`)
sum_0500_R3 = sum(my_data$`Intensity.ADH_NHS_Ac_0500_R3`)

sum_1000_R1 = sum(my_data$`Intensity.ADH_NHS_Ac_1000_R1`)
sum_1000_R2 = sum(my_data$`Intensity.ADH_NHS_Ac_1000_R2`)
sum_1000_R3 = sum(my_data$`Intensity.ADH_NHS_Ac_1000_R3`)

sum_1500_R1 = sum(my_data$`Intensity.ADH_NHS_Ac_1500_R1`)
sum_1500_R2 = sum(my_data$`Intensity.ADH_NHS_Ac_1500_R2`)
sum_1500_R3 = sum(my_data$`Intensity.ADH_NHS_Ac_1500_R3`)



# open K and YST modifications csv file 

mod_sites <- read.csv("mod_sites.csv", header = T, stringsAsFactors = F)
dim(mod_sites)
head(mod_sites)

# select rows without contaminats
mod_sites2 = mod_sites[-grep("^CON", mod_sites$Proteins),]
mod_sites = mod_sites2

# select rows according to scroe and localisation probability
mod_sites2 = mod_sites[mod_sites$Score > 80,]
mod_sites = mod_sites2

mod_sites2 = mod_sites[mod_sites$Localization.prob > 0.75,]
mod_sites = mod_sites2

#Buit ratio

mod_sites$ratio_Intensity_ADH_NHS_Ac_0_R1 = mod_sites$Intensity.ADH_NHS_Ac_0_R1 / sum_0000_R1
mod_sites$ratio_Intensity_ADH_NHS_Ac_0_R2 = mod_sites$Intensity.ADH_NHS_Ac_0_R2 / sum_0000_R2
mod_sites$ratio_Intensity_ADH_NHS_Ac_0_R3 = mod_sites$Intensity.ADH_NHS_Ac_0_R3 / sum_0000_R3

mod_sites$ratio_Intensity_ADH_NHS_Ac_0050_R1 = mod_sites$Intensity.ADH_NHS_Ac_0050_R1 / sum_0050_R1
mod_sites$ratio_Intensity_ADH_NHS_Ac_0050_R2 = mod_sites$Intensity.ADH_NHS_Ac_0050_R2 / sum_0050_R2
mod_sites$ratio_Intensity_ADH_NHS_Ac_0050_R3 = mod_sites$Intensity.ADH_NHS_Ac_0050_R3 / sum_0050_R3

mod_sites$ratio_Intensity_ADH_NHS_Ac_0100_R1 = mod_sites$Intensity.ADH_NHS_Ac_0100_R1 / sum_0100_R1
mod_sites$ratio_Intensity_ADH_NHS_Ac_0100_R2 = mod_sites$Intensity.ADH_NHS_Ac_0100_R2 / sum_0100_R2
mod_sites$ratio_Intensity_ADH_NHS_Ac_0100_R3 = mod_sites$Intensity.ADH_NHS_Ac_0100_R3 / sum_0100_R3

mod_sites$ratio_Intensity_ADH_NHS_Ac_0250_R1 = mod_sites$Intensity.ADH_NHS_Ac_0250_R1 / sum_0250_R1
mod_sites$ratio_Intensity_ADH_NHS_Ac_0250_R2 = mod_sites$Intensity.ADH_NHS_Ac_0250_R2 / sum_0250_R2
mod_sites$ratio_Intensity_ADH_NHS_Ac_0250_R3 = mod_sites$Intensity.ADH_NHS_Ac_0250_R3 / sum_0250_R3

mod_sites$ratio_Intensity_ADH_NHS_Ac_0500_R1 = mod_sites$Intensity.ADH_NHS_Ac_0500_R1 / sum_0500_R1
mod_sites$ratio_Intensity_ADH_NHS_Ac_0500_R2 = mod_sites$Intensity.ADH_NHS_Ac_0500_R2 / sum_0500_R2
mod_sites$ratio_Intensity_ADH_NHS_Ac_0500_R3 = mod_sites$Intensity.ADH_NHS_Ac_0500_R3 / sum_0500_R3

mod_sites$ratio_Intensity_ADH_NHS_Ac_1000_R1 = mod_sites$Intensity.ADH_NHS_Ac_1000_R1 / sum_1000_R1
mod_sites$ratio_Intensity_ADH_NHS_Ac_1000_R2 = mod_sites$Intensity.ADH_NHS_Ac_1000_R2 / sum_1000_R2
mod_sites$ratio_Intensity_ADH_NHS_Ac_1000_R3 = mod_sites$Intensity.ADH_NHS_Ac_1000_R3 / sum_1000_R3

mod_sites$ratio_Intensity_ADH_NHS_Ac_1500_R1 = mod_sites$Intensity.ADH_NHS_Ac_1500_R1 / sum_1500_R1
mod_sites$ratio_Intensity_ADH_NHS_Ac_1500_R2 = mod_sites$Intensity.ADH_NHS_Ac_1500_R2 / sum_1500_R2
mod_sites$ratio_Intensity_ADH_NHS_Ac_1500_R3 = mod_sites$Intensity.ADH_NHS_Ac_1500_R3 / sum_1500_R3


####Melt data#####

#library(reshape)
#md <- melt(mydata, id=(c("id", "time")))
mod_sites$Position
mod_sites$Amino.acid

mod_sites_for_melt = mod_sites %>% select(Position, Amino.acid, ratio_Intensity_ADH_NHS_Ac_0_R1:ratio_Intensity_ADH_NHS_Ac_1500_R3)
mod_sites_melt <- melt(mod_sites_for_melt, id=(c("Position", "Amino.acid")))


#### add new column with concentration information#####
mod_sites_melt$Concentration <- 0
mod_sites_melt$Concentration[mod_sites_melt$variable %in% c("ratio_Intensity_ADH_NHS_Ac_0050_R1",
                                                            "ratio_Intensity_ADH_NHS_Ac_0050_R2",
                                                            "ratio_Intensity_ADH_NHS_Ac_0050_R3")] <- 50
mod_sites_melt$Concentration[mod_sites_melt$variable %in% c("ratio_Intensity_ADH_NHS_Ac_0100_R1",
                                                            "ratio_Intensity_ADH_NHS_Ac_0100_R2",
                                                            "ratio_Intensity_ADH_NHS_Ac_0100_R3")] <- 100
mod_sites_melt$Concentration[mod_sites_melt$variable %in% c("ratio_Intensity_ADH_NHS_Ac_0250_R1",
                                                            "ratio_Intensity_ADH_NHS_Ac_0250_R2",
                                                            "ratio_Intensity_ADH_NHS_Ac_0250_R3")] <- 250
mod_sites_melt$Concentration[mod_sites_melt$variable %in% c("ratio_Intensity_ADH_NHS_Ac_0500_R1",
                                                            "ratio_Intensity_ADH_NHS_Ac_0500_R2",
                                                            "ratio_Intensity_ADH_NHS_Ac_0500_R3")] <- 500
mod_sites_melt$Concentration[mod_sites_melt$variable %in% c("ratio_Intensity_ADH_NHS_Ac_1000_R1",
                                                            "ratio_Intensity_ADH_NHS_Ac_1000_R2",
                                                            "ratio_Intensity_ADH_NHS_Ac_1000_R3")] <- 1000
mod_sites_melt$Concentration[mod_sites_melt$variable %in% c("ratio_Intensity_ADH_NHS_Ac_1500_R1",
                                                            "ratio_Intensity_ADH_NHS_Ac_1500_R2",
                                                            "ratio_Intensity_ADH_NHS_Ac_1500_R3")] <- 1500



####sort positions
mod_sites_melt$Position <- as.integer(mod_sites_melt$Position)
mod_sites_melt$Position <- as.character(mod_sites_melt$Position)
mod_sites_melt$Position <- factor(mod_sites_melt$Position, levels = sort(unique(as.integer(mod_sites_melt$Position))))

###add Replicate information
mod_sites_melt$Replicate = mod_sites_melt$variable

mod_sites_melt$Replicate = gsub("ratio_Intensity_ADH_NHS_Ac_0_R", "", mod_sites_melt$Replicate)
mod_sites_melt$Replicate = gsub("ratio_Intensity_ADH_NHS_Ac_0050_R", "", mod_sites_melt$Replicate)
mod_sites_melt$Replicate = gsub("ratio_Intensity_ADH_NHS_Ac_0100_R", "", mod_sites_melt$Replicate)
mod_sites_melt$Replicate = gsub("ratio_Intensity_ADH_NHS_Ac_0250_R", "", mod_sites_melt$Replicate)
mod_sites_melt$Replicate = gsub("ratio_Intensity_ADH_NHS_Ac_0500_R", "", mod_sites_melt$Replicate)
mod_sites_melt$Replicate = gsub("ratio_Intensity_ADH_NHS_Ac_1000_R", "", mod_sites_melt$Replicate)
mod_sites_melt$Replicate = gsub("ratio_Intensity_ADH_NHS_Ac_1500_R", "", mod_sites_melt$Replicate)



###create data for linear models without 0 or NA
id1= mod_sites_melt[mod_sites_melt$value > 0,]
id1$Concentration2= id1$Concentration
id1$Concentration2= as.character(id1$Concentration2)
id1$Concentration2 =gsub("^0","1",id1$Concentration2)
id = id1



####linear model with looping#####

id_split <- split(id, f = id$Position)

l_mods <- list()
for(i in seq_along(id_split)){
  
  temp    <- id_split[[i]]
  temp    <- temp[!temp$value == 0, ]
  l_mod   <- lm(log2(value) ~ log2(as.numeric(Concentration2)), data = temp, singular.ok = TRUE)
  l_mods[[i]] <- l_mod
  
}
names(l_mods) <- names(id_split)

model_data <- lapply(l_mods, function(mod) {
  
  data.frame(intercept = summary(mod)$coefficients[1],
             slope     = summary(mod)$coefficients[2],
             intercept_se     = summary(mod)$coefficients[3],
             slope_se         = summary(mod)$coefficients[4],
             intercept_tvalue = summary(mod)$coefficients[5],
             slope_tvalue     = summary(mod)$coefficients[6],
             intercept_pvalue = summary(mod)$coefficients[7],
             slope_pvalue     = summary(mod)$coefficients[8],
             stringsAsFactors = FALSE
  )
  
})
model_data

#### combine model_data loop into a single data frame and save as xlxs ####
model_data <- do.call("rbind", model_data)
head(model_data)
nrow(model_data)
length(l_mods)

# name the positions
model_data$Position <- names(l_mods)
head(model_data)

#save output as excel file

write.xlsx(model_data, "20190515_ADH_NHS_Norm_Intensities_With_MQ_Peptide_sum")



####plot data#####

####plot id (data without 0 or NA)
x11()
g <- ggplot(id,
            aes(x = log2(as.numeric(Concentration2)), y = log2(value), color = Position, group = Position))
#g <- g + geom_point(aes(color = Experiment))
#g <- g + geom_line()
g <- g + stat_smooth(method = "lm", se= FALSE)
g <- g + theme(axis.text.x = element_text(angle = 30, hjust = 1))
g <- g + scale_y_continuous(name = "log2(Ratio)")
g <- g + scale_x_continuous(name = "log2(Concentration)")
g




gd1 <- as_data_frame(mod_sites_melt)
gd <- gd1 %>% 
  group_by(Position, Concentration) %>% 
  summarise(Ratio.mean = mean(value),
            Ratio.sd  = sd(value))

plot_mean = ggplot(data=gd)+
  geom_point(aes(x=log2(Concentration), y=log2(Ratio.mean), col=Position))
plot_mean

###statistics 
#id = replace_with_na(replace = list(x==0))
gd2 <- summarySE(id, measurevar="value", groupvars=c("Position","Concentration"))
gd2

# plot Standard error of the mean
x11()
plot_SE = ggplot(gd2, aes(x= log2(Concentration), y=log2(value), colour=Position)) + 
  geom_errorbar(aes(ymin=log2(value-se), ymax=log2(value+se)), width=.1) +
  
  geom_point()
plot_SE

#geom_line() +

# for single residues
single_res= gd2[gd2$Position ==  "28", ]

plot_single_res = ggplot(data=single_res) +
  geom_point(aes(x = log2(Concentration), y = log2(value)), size=2) +
  scale_y_continuous(name = "Ratio") +
  scale_x_continuous(name = "Concentration")+
  geom_smooth(aes(x = log2(Concentration), y = log2(value), method="lm", se=FALSE))

plot_single_res


single_res= id[id$Position ==  "28", ]

plot_single_res = ggplot(data=single_res) +
  geom_point(aes(x = log2(as.numeric(Concentration2)), y = log2(value)), size=2) +
  scale_y_continuous(name = "log2(Ratio)") +
  scale_x_continuous(name = "log2(Concentration)")+
  geom_smooth(aes(x = log2(Concentration2), y = log2(value), method="lm", se=FALSE))

plot_single_res




single_res2= gd1[gd1$Position == "315",]
head(single_res2)

plot_single_res2 = ggplot(data=single_res2) +
  geom_point(aes(x = log2(Concentration), y = log2(value)), col = single_res2$Replicate, group=single_res2$Replicate) +
  scale_y_continuous(name = "log2(Ratio)") +
  scale_x_continuous(name = "log2(Concentration)")+
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
  
x11()
plot_single_res2                 

plot_single_res2 = ggplot(data=single_res2) +
  geom_point(aes(x = Concentration, y = log2(value)), col = single_res2$Replicate, group=single_res2$Replicate, fill=single_res2$Replicate) +
  scale_y_continuous(name = "log2(Ratio)") +
  scale_x_continuous(name = "Concentration")+
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

x11()
plot_single_res2   



#plot for every residue
#plot.res <- ggplot(gd2, aes(x=log2(Concentration), y=log2(value))) + 
#geom_point() + facet_grid(. ~ Position) + stat_smooth(method = "lm") +
#background_grid(major = 'y', minor = "none") + # add thin horizontal lines 
#panel_border() # and a border around each panel
#x11(plot.res)


#linear models
#fits <- lmList(Concentration ~ value | Position, data=id, na.action=na.exclude)
#fits








####plots####
plot1 <- ggplot(data=mod_sites_melt) + 
  geom_point(aes(x = Concentration, y = log2(value), col = Position, shape = Amino.acid, group = Position), size=2) +
  geom_line(aes(x = Concentration, y = log2(value), col = Position), size=0.5)+
  scale_y_continuous(name = "Ratio") +
  scale_x_continuous(name = "Concentration") 

plot1

plot2 <- ggplot(data=mod_sites_melt) + 
  geom_point(aes(x = log2(Concentration), y = log2(value), col = Position, shape = Amino.acid), size=2) +
  geom_smooth(aes(x = log2(Concentration), y = log2(value), method="lm"))
  scale_y_continuous(name = "Ratio") +
  scale_x_continuous(name = "Concentration") 
plot2

plot3 <- ggplot(data=mod_sites_melt) + 
  geom_point(aes(x = log2(Concentration), y = log2(value), col = Position, shape = Amino.acid), size=2) +
  geom_line(aes(x = log2(Concentration), y = log2(value), col = Position), size=0.5)+
  scale_y_continuous(name = "Ratio") +
  scale_x_continuous(name = "Concentration")+
  geom_smooth(aes(x = log2(Concentration), y = log2(value), method="lm"))

plot3

data1= mod_sites_melt[mod_sites_melt$Position ==  "27", ]

plot4 = ggplot(data=data1) +
                 geom_point(aes(x = (Concentration), y = log2(value), col = Replicate, shape = Amino.acid, group=Replicate), size=2) +
                 geom_line(aes(x = (Concentration), y = log2(value), col = Replicate), size=0.5)+
                 scale_y_continuous(name = "Ratio") +
                 scale_x_continuous(name = "Concentration")+
                geom_smooth(aes(x = (Concentration), y = log2(value), method="lm"))

plot4






#library(tidyverse)

#mode_site_new <- mutate(mod_sites, ratio_Intensity_ADH_NHS_Ac_0_R1 = ifelse(Score.ADH_NHS_Ac_0050_R1 > 80 | NA |
#                                                                            Localization.prob.ADH_NHS_Ac_0050_R1>0.75| NA, 
#                                                                            ratio_Intensity_ADH_NHS_Ac_0_R1, NA))


#mode_site_new <- mutate(mod_sites, ratio_Intensity_ADH_NHS_Ac_0_R1 = ifelse(mod_sites$Localization.prob.ADH_NHS_Ac_0050_R1>0.75, ratio_Intensity_ADH_NHS_Ac_0_R1, NA))

#replace Ratio with NA if Identification typ !== by MS/MS








