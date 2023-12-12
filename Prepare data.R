library(this.path)
library(dplyr)
library(tidyverse)

setwd(this.dir()) # set the active working directory to the directory of this file
getwd()

data <- read.csv("20231027 Afkeurbeleid Donaties met leeftijd JK.csv", sep = "\t", head = T, fileEncoding="UTF-16LE")

test <- aggregate(data$Donatie.type.omschrijving, by = list(data$Donor.code, data$Donatie.type.omschrijving), FUN=length)

test2 <- data %>% group_by(Donor.code, Donatie.type.omschrijving) %>% summarise(n=n()) 
test2 <- test2 %>% pivot_wider(names_from = Donatie.type.omschrijving, values_from = n)
test3<- test2 %>% filter((Bloed > 5*Plasma | is.na(Plasma) ) & ( Bloed>5*Cyta | is.na(Cyta)))

data <- data %>% dplyr::select(Donor.code, Donor.geslacht.code, Datum.donatie, Donatie.hemoglobine.van.donor,  Donatie.HB.predonatie.screening, Donatie.type.omschrijving) %>% 
  filter(Donor.code %in% test3$Donor.code) %>% mutate(Datum.donatie=as.Date(Datum.donatie, format = "%Y-%m-%d"), Donatie.HB.predonatie.screening=as.numeric(Donatie.HB.predonatie.screening)) %>%
  filter(Donatie.type.omschrijving=="Bloed") %>% dplyr::select(-Donatie.type.omschrijving)

colnames(data) <- c("KeyID", "Sex","DonDate", "Hb", "pre_Hb")
data$Sex[data$Sex=="V"]<-"F"
data <- data %>% filter(!(is.na(data$Hb)&is.na(data$pre_Hb)))

#alter prescreening values that were likely a typo and remove unlikely values
data$pre_Hb[data$pre_Hb > 9 & data$pre_Hb < 20 & !is.na(data$Hb) & !is.na(data$pre_Hb)] <- data$pre_Hb[data$pre_Hb > 9 & data$pre_Hb < 20 & !is.na(data$Hb) & !is.na(data$pre_Hb)]*10
data$pre_Hb[data$pre_Hb <= 9 & !is.na(data$Hb)] <- NA



saveRDS(data, file="donations.rds")
