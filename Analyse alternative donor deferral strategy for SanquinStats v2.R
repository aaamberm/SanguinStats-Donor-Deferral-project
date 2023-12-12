#######################################################
# Analyse impact of alternative deferral strategy
#######################################################
# initialize R
rm(list=ls())
while (!is.null(dev.list())) dev.off()
cat("\014")  

# set the working path
library(this.path)
setwd(this.dir()) # set the active working directory to the directory of this file
getwd()

library("car")


# Set code date 
maincodedatestamp<-"20231115"

#############################
# Define user input
#############################

# set name of the datafile to use
FileToUse<-"donations.rds" # to be set by the USER 

# "FileToUSe" should contain the following data columns (data type in brackets):
# KeyID   : Unique identifier for each donor (integer)
# Sex     : indicator for Male (M) or Female (F) donor (factor)
# DonDate : date of donation (Date)
# Hb      : donor Hb at donation date (numeric)

# All relevant input and output information is stored in the "tosave" variable
# This information is stored in a file called "SavedDeferralData_DD-MM-YYYY.RDS"

# parameter to determine whether the plots go to a pdf file or to screen.
plot_to_pdf<-T # to be set by the USER 


# Set minimum acceptable Hb levels for males and females
dtm<-135  # to be set by the USER 
dtf<-125  # to be set by the USER 

# Are Hb levels given in g/L (T) or in mmol/L (F)
Hb_in_gpl<-T# to be set by the USER 

# European threshold (g/L)
135*0.06206 # 135g/L for males = 8.3781 mmol/L
round(135*0.06206,1) # 8.4 mmol/L
125*0.06206 # 125g/L for females = 7.7575 mmol/L
round(125*0.06206,1) # 7.8 mmol/L

# Is a change in anonymous donor IDs required?
changeIDs<-F # to be set by the USER 

# Set deferral percentile level
cutoffperc<-0.99 # to be set by the USER 

# set minimum size of the groups for aggregated data
mingroupsize<-20 # to be set by the USER 

# Note that the analyses may take a substantial amount of time. Therefore an analysis file 
# ("donations_analysis_data.RDS") will be created, which allows speeding up the analyses
# by reading in this file whenever the analyses are repeated (e.g. for making additional graphs)
# Do not forget to remove this file in case the original data file is updated!

#############################
# Initialize R
#############################
# load various R libraries required for the analyses
library("zoo")       # required for rollmean function
library("lubridate") # required for extracting year info from a date variable
library("Hmisc")     # to enable calculating splitpoints
library("fitdistrplus") # to fit and plot normal distribution fits to the data
library("stringr")   # to manipulate strings
library("dplyr")     # data wrangling

# restore old color palette
palette("R3")

# load support functions from the second code file
source("General_functions.R")

#############################
# read data
#############################
# but first write some of the input data to the tosave variable
# the first item is the name of the datafile used for the analyses
tosave<-list(FileToUse=FileToUse)

# save code datestamps
tosave<-append(tosave, list(maincodedatestamp=maincodedatestamp))
tosave<-append(tosave, list(generalfunctionscodedatestamp=generalfunctionscodedatestamp))

# save cutoff percentage applied
tosave<-append(tosave, list(cutoffperc=cutoffperc))

# now read in the data
data<-readRDS(FileToUse)
(classes<-sapply(data,class))
#     KeyID         Sex     DonDate          Hb    pre_Hb
# "integer" "character"      "Date"   "numeric"   "numeric"


# when was the first donation
daterange<-min(data$DonDate)
(daterange<-c(daterange,max(data$DonDate)))
nrrecs<-nrow(data)

# save info
tosave<-append(tosave, list(daterange=daterange))

# remove all missing records
sum(is.na(data$Hb))
data_complete <- data #Amber: hier is/was een probleem, namelijk donaties met een pre-donatiescreening onder de treshold hebben een NA Hb waarde (omdat er na de predonatie screening geen donatie was), maar we hebben in de verdere analyse wel beide nodig
data_prescreening <- data[!is.na(data$pre_Hb),] #Amber: make a separate dataset with all records where a pre-screening was conducted
data<-data[!is.na(data$Hb),] 
nrrecs<-c(nrrecs,nrow(data))

data<-data[!is.na(data$KeyID),]
nrrecs<-c(nrrecs,nrow(data))
data<-data[!is.na(data$DonDate),]
nrrecs<-c(nrrecs,nrow(data))
nrrecs<-c(nrrecs,length(unique(data$KeyID)))

# save changes in dataset
tosave<-append(tosave, list(nrrecs=nrrecs))

# Sort by date per donor
data<-data[order(data$KeyID,data$DonDate),] 
data_complete<-data_complete[order(data_complete$KeyID,data_complete$DonDate),] 
# Create index for nr of donations
data$numdons <- sequence(rle(data$KeyID)$lengths)
data_complete$numdons <- sequence(rle(data_complete$KeyID)$lengths)

# nr of donors
length(unique(data$KeyID)) 
# nr of donations
length((data$KeyID))       

# change KeyIDs if requested by the user
if (changeIDs) {
  donIDs<-unique(data$KeyID)
  donIDs<-cbind(donIDs, 1:length(donIDs))
  #donIDs<-as.data.frame(cbind(donIDs, 1:length(donIDs)))
  colnames(donIDs)<-c("KeyID", "NewID")
  data<-merge(data,donIDs, by="KeyID")
  data$KeyID<-NULL
  colnames(data)[which(colnames(data)=="NewID")]<-"KeyID"
}

# calculate distribution of nr of donations per donor
table(data$numdons, data$Sex) 
dontable <- table(data$numdons, data$Sex) 
tosave<-append(tosave, list(dontable=dontable))


# convert Hb levels if so required
if (!Hb_in_gpl){
  dtm<-dtm/0.06206 -1e-6 # subtract a small margin to compensate for rounding errors
  dtf<-dtf/0.06206 -1e-6
  data$Hb<-data$Hb/0.06206
  data$pre_Hb<-data$pre_Hb/0.06206
  data_complete$Hb <- data_complete$Hb/0.06206
  data_complete$pre_Hb <- data_complete$pre_Hb/0.06206
  data_prescreening$Hb <- data_prescreening$Hb/0.06206
  data_prescreening$pre_Hb <- data_prescreening$pre_Hb/0.06206
}

# calculate sex, mean Hb, Sd and nr of donations per donor
meanHb<-aggregate(data$Hb, by=list(data$KeyID), mean)
SdHb<-aggregate(data$Hb, by=list(data$KeyID), sd)
nHb<-aggregate(data$Hb, by=list(data$KeyID, data$Sex), length)
# aggregate these results
adata<-merge(meanHb,SdHb, by="Group.1")
adata<-merge(adata,nHb, by="Group.1")
adata$Group.1<-NULL
colnames(adata)<-c("Hb", "sd", "Sex", "Nrdon")

# calculate distributions for various subsets of nr of donations

if(plot_to_pdf) pdf(file="Hb_distribution_male.pdf")
# nrofquantiles are to be set by the USER 
malefits  <-fitHbdistributions(adata[adata$Sex=="M",],nrofquantiles=20)
if(plot_to_pdf) dev.off()

if(plot_to_pdf) pdf(file="Hb_distribution_female.pdf")
# nrofquantiles are to be set by the USER 
femalefits<-fitHbdistributions(adata[adata$Sex=="F",],nrofquantiles=20)
if(plot_to_pdf) dev.off()

# stop execution of groupsize is larger than required by the user
if(malefits$minsubset<mingroupsize | femalefits$minsubset<mingroupsize) {
  print("Code stopped because aggregated group size is smaller than specified by the user")
  print("Please decrease the nrofquantiles parameter in the fitHbdistributions functions (line 144/145)")
  print("or increase the mingroupsize (line 44)")
  stop("Change analysis code")
}

# set parameter for maximum follow-up
maxDons<-max(nHb$x)

# save distribution fits and maxDons
tosave<-append(tosave, list(malefits=malefits))
tosave<-append(tosave, list(femalefits=femalefits))
tosave<-append(tosave, list(maxDons=maxDons))

# Postdonation screenings 
# Set indicator for deferral
data$def<-ifelse(data$Hb<dtf,1,0)
data$def[data$Sex=="M"]<-ifelse(data$Hb[data$Sex=="M"]<dtm,1,0)

#create year variable
data$year<-year(data$DonDate)

# table deferrals per year 
with(data,table(year, Sex))
with(data,table(Sex,year,def))
deff<-with(data[data$Sex=="F",],table(year,def))
proportions(deff, margin=1)
defm<-with(data[data$Sex=="M",],table(year,def))
proportions(defm, margin=1)
tosave<-append(tosave, list(defm=defm))
tosave<-append(tosave, list(deff=deff))

# plot deferrals per year
maxdef<-max(c(proportions(defm, margin=1)[,2],proportions(deff, margin=1)[,2]))
if(plot_to_pdf) pdf(file="Deferrals_per_year.pdf")
plot(rownames(defm), proportions(defm, margin=1)[,2], ylim=c(0,maxdef*1.2), col="blue", type="l",
     ylab="Proportion deferred", xlab="Year")
lines(rownames(deff), proportions(deff, margin=1)[,2], col="red")
points(rownames(defm), proportions(defm, margin=1)[,2], pch=1, col="blue")
points(rownames(deff), proportions(deff, margin=1)[,2], pch=2, col="red")
legend("topright", c("Males", "Females"), pch=c(1,2), lty=c(1,1), col=c("blue","red"))
if(plot_to_pdf) dev.off()

if(plot_to_pdf) pdf(file="Distribution_of_donation_counts.pdf")
# plot distribution of number of donations per sex
with(data[data$Sex=="F",], plot(as.numeric(table(numdons)), type="l", col="red", xlim=c(1,maxDons), log='y',
                                xlab="Number of donations", ylab="number of donors"))
with(data[data$Sex=="M",], lines(as.numeric(table(numdons)), type="l", col="blue"))
with(data[data$Sex=="M",], points(as.numeric(table(numdons)), pch=1, col="blue"))
with(data[data$Sex=="F",], points(as.numeric(table(numdons)), pch=2, col="red"))
legend("topright", c("Males", "Females"), pch=c(1,2), lty=c(1,1), col=c("blue","red"))
if(plot_to_pdf) dev.off()

# Predonation screenings 
# Set indicator for deferral
data_prescreening$def<-ifelse(data_prescreening$pre_Hb<dtf,1,0)
data_prescreening$def[data_prescreening$Sex=="M"]<-ifelse(data_prescreening$pre_Hb[data_prescreening$Sex=="M"]<dtm,1,0)

#create year variable
data_prescreening$year<-year(data_prescreening$DonDate)

# table deferrals per year 
with(data_prescreening,table(year, Sex))
with(data_prescreening,table(Sex,year,def))
pre_deff<-with(data_prescreening[data_prescreening$Sex=="F",],table(year,def))
proportions(pre_deff, margin=1)
pre_defm<-with(data_prescreening[data_prescreening$Sex=="M",],table(year,def))
proportions(pre_defm, margin=1)
tosave<-append(tosave, list(pre_defm=pre_defm))
tosave<-append(tosave, list(pre_deff=pre_deff))

# plot deferrals per year
pre_maxdef<-max(c(proportions(pre_defm, margin=1)[,2],proportions(pre_deff, margin=1)[,2]))
if(plot_to_pdf) pdf(file="Predonationscreening_Deferrals_per_year.pdf")
plot(rownames(pre_defm), proportions(pre_defm, margin=1)[,2], ylim=c(0,pre_maxdef*1.2), col="blue", type="l",
     ylab="Proportion deferred", xlab="Year")
lines(rownames(pre_deff), proportions(pre_deff, margin=1)[,2], col="red")
points(rownames(pre_defm), proportions(pre_defm, margin=1)[,2], pch=1, col="blue")
points(rownames(pre_deff), proportions(pre_deff, margin=1)[,2], pch=2, col="red")
legend("topright", c("Males", "Females"), pch=c(1,2), lty=c(1,1), col=c("blue","red"))
if(plot_to_pdf) dev.off()

#all postdonation screenings and predonation screenings
# Set indicator for deferral
data_complete$def<-ifelse(data_complete$Hb<dtf,1,0)
data_complete$def[data_complete$Sex=="M"]<-ifelse(data_complete$Hb[data_complete$Sex=="M"]<dtm,1,0)
data_complete$predef<-ifelse(data_complete$pre_Hb<dtf,1,0)
data_complete$predef[data_complete$Sex=="M"]<-ifelse(data_complete$pre_Hb[data_complete$Sex=="M"]<dtm,1,0)
data_complete$def[is.na(data_complete$def)] <- 0
data_complete$predef[is.na(data_complete$predef)] <- 0 
data_complete$sumdef <- (data_complete$def | data_complete$predef)*1

#create year variable
data_complete$year<-year(data_complete$DonDate)

# table deferrals per year 
with(data_complete,table(year, Sex))
with(data_complete,table(Sex,year,sumdef))
all_deff<-with(data_complete[data_complete$Sex=="F",],table(year,sumdef))
proportions(all_deff, margin=1)
all_defm<-with(data_complete[data_complete$Sex=="M",],table(year,sumdef))
proportions(all_defm, margin=1)
tosave<-append(tosave, list(all_defm=all_defm))
tosave<-append(tosave, list(all_deff=all_deff))

# plot deferrals per year
all_maxdef<-max(c(proportions(all_defm, margin=1)[,2],proportions(all_deff, margin=1)[,2]))
if(plot_to_pdf) pdf(file="alldonationscreening_Deferrals_per_year.pdf")
plot(rownames(all_defm), proportions(all_defm, margin=1)[,2], ylim=c(0,all_maxdef*1.2), col="blue", type="l",
     ylab="Proportion deferred", xlab="Year")
lines(rownames(all_deff), proportions(all_deff, margin=1)[,2], col="red")
points(rownames(all_defm), proportions(all_defm, margin=1)[,2], pch=1, col="blue")
points(rownames(all_deff), proportions(all_deff, margin=1)[,2], pch=2, col="red")
legend("topright", c("Males", "Females"), pch=c(1,2), lty=c(1,1), col=c("blue","red"))
if(plot_to_pdf) dev.off()

# table deferrals per year for first time donors
with(data_complete[data_complete$numdons==1 & !is.na(data_complete$pre_Hb),],table(year, Sex))
with(data_complete[data_complete$numdons==1 & !is.na(data_complete$pre_Hb),],table(Sex,year,sumdef))
first_deff<-with(data_complete[data_complete$numdons==1 & !is.na(data_complete$pre_Hb) & data_complete$Sex=="F",],table(year,sumdef))
proportions(first_deff, margin=1)
first_defm<-with(data_complete[data_complete$numdons==1 & !is.na(data_complete$pre_Hb) & data_complete$Sex=="M",],table(year,sumdef))
proportions(first_defm, margin=1)
tosave<-append(tosave, list(first_defm=first_defm))
tosave<-append(tosave, list(first_deff=first_deff))

# plot deferrals per year
first_maxdef<-max(c(proportions(first_defm, margin=1)[,2],proportions(first_deff, margin=1)[,2]))
if(plot_to_pdf) pdf(file="First_donationscreening_Deferrals_per_year.pdf")
plot(rownames(first_defm), proportions(first_defm, margin=1)[,2], ylim=c(0,first_maxdef*1.2), col="blue", type="l",
     ylab="Proportion deferred", xlab="Year")
lines(rownames(first_deff), proportions(first_deff, margin=1)[,2], col="red")
points(rownames(first_defm), proportions(first_defm, margin=1)[,2], pch=1, col="blue")
points(rownames(first_deff), proportions(first_deff, margin=1)[,2], pch=2, col="red")
legend("topright", c("Males", "Females"), pch=c(1,2), lty=c(1,1), col=c("blue","red"))
if(plot_to_pdf) dev.off()


##################################
# Estimate measurement variation
##################################

# create an index for pointing to the previous donation
idx<-1:nrow(data)
precursor<-idx-1
precursor[1]<-nrow(data)

# calculate time since last donation
data$dt<-NA
data$dt<-as.numeric(data$DonDate-data$DonDate[precursor])
data$dt[data$KeyID != data$KeyID[precursor]]<-NA
# calculate log10 of the time difference
data$ldt<-log10(data$dt)

# calculate change in Hb from previous donation
data$dHb<-NA
data$dHb<-data$Hb-data$Hb[precursor]
data$dHb[data$KeyID != data$KeyID[precursor]]<-NA

# also make a column for the previous Hb value
data$prevHb<-NA
data$prevHb<-data$Hb[precursor]
data$prevHb[data$KeyID != data$KeyID[precursor]]<-NA


# set plot parameters # to be set by the USER 
intervaltoshow<-c(50, 730) # time interval to show
ylim<-c(-40, 40)           # interval in dHb to show
nrtoprint<-15000           # nr of observations to select for printing
rollmeanWidth<-1000        # width of the rolling window

# plot association between time and Hb change for females
set.seed(1)
self<-which(data$Sex=="F" & data$dt>=intervaltoshow[1] & !is.na(data$ldt))
sel<-sample(self, nrtoprint, replace=F)
sel<-sel[order(data$dt[sel])]
if(plot_to_pdf) pdf(file="Change_in_Hb_level_by_interval_female.pdf")
with(data[sel,], plot(dt, jitter(dHb, amount=.5), col="red", log="x", xlim=intervaltoshow, ylim=ylim,
                      xlab="Days between donations", ylab="Change in Hb level [g/L]"))
abline(h=0,col=8)
# add rolling mean
with(data[sel,], lines(dt, rollmean(dHb, rollmeanWidth, fill = list(NA, NULL, NA)),
                       col = 3, lwd = 3 ))
linfitf<-lm(dHb~ldt, data=data[self,])
summary(linfitf)
cx<-as.data.frame(log10(intervaltoshow))
colnames(cx)<-"ldt"
lines(intervaltoshow, predict(linfitf,cx), lwd=2)
if(plot_to_pdf) dev.off()
mean(data$Hb[self])   
sd(data$Hb[self])     
sd(data$dHb[self])    
sd(linfitf$residuals) 

# save info
tosave<-append(tosave, list(coeff=linfitf$coefficients))
tosave<-append(tosave, list(sdf=c(mean(data$Hb[self]), sd(data$Hb[self]), 
                                  sd(data$dHb[self]), sd(linfitf$residuals))))

# plot association between time and Hb change for males
intervaltoshow<-c(50, 730) # time interval to show
set.seed(1)
selm<-which(data$Sex=="M" & data$dt>=intervaltoshow[1] & !is.na(data$ldt))
sel<-sample(selm, nrtoprint, replace=F)
sel<-sel[order(data$dt[sel])]
if(plot_to_pdf) pdf(file="Change_in_Hb_level_by_interval_male.pdf")
with(data[sel,], plot(dt, jitter(dHb, amount=.5), col="blue", log="x", xlim=intervaltoshow, ylim=ylim,
                      xlab="Days between donations", ylab="Change in Hb level [g/L]"))
abline(h=0,col=8)
# add rolling mean
with(data[sel,], lines(dt, rollmean(dHb, rollmeanWidth, fill = list(NA, NULL, NA)),
                       col = 3, lwd = 3 ))
linfitm<-lm(dHb~ldt, data=data[selm,])
summary(linfitm)
lines(intervaltoshow, predict(linfitm,cx), lwd=2)

if(plot_to_pdf) dev.off()
mean(data$Hb[selm])   
sd(data$Hb[selm])     
sd(data$dHb[selm])    
sd(linfitm$residuals) 


# save info
tosave<-append(tosave, list(coefm=linfitm$coefficients))
tosave<-append(tosave, list(sdm=c(mean(data$Hb[selm]), sd(data$Hb[selm]), 
                                  sd(data$dHb[selm]), sd(linfitm$residuals))))


###########################################################
# Hb recovery extra analysis
###########################################################
#index for the previous donation
idx<-1:nrow(data_complete)
precursor<-idx-1
precursor[1]<-nrow(data_complete)

# check whether the previous observation was indeed a donation
idx_noprev_donation<-idx[data_complete$KeyID == data_complete$KeyID[precursor] & is.na(data_complete$Hb[precursor])]
precursor2<-precursor[idx_noprev_donation]-1
if(F){
  ### change some records to test 
  idx2<-which(data_complete$KeyID[idx_noprev_donation]==data_complete$KeyID[precursor2])
  data_complete[precursor2[idx2[1:10]],]<-data_complete[(precursor[idx_noprev_donation])[idx2[1:10]],]
  View(data_complete[data_complete$KeyID %in% data_complete$KeyID[idx_noprev_donation[idx2][1:10]],])
}
precisvalid<-data_complete$KeyID[idx_noprev_donation] == data_complete$KeyID[precursor2]
precisHb<-!is.na(data_complete$Hb[precursor2])
sum(precisvalid) # previous record was from the same donor
sum(precisHb) # previous record had a valid Hb
precisset<-precisvalid & precisHb # reference is set if both conditions are satisfied
sum(precisset)
i<-0
while ( sum(!precisvalid | precisset)!=length(precursor2)){ # loop until all previous records found are either invalid or contain an Hb value 
  i<-i+1
  precursor2[!precisset & precisvalid]<-precursor2[!precisset & precisvalid]-1
  precisvalid[!precisset & precisvalid]<-(data_complete$KeyID[idx_noprev_donation] == data_complete$KeyID[precursor2])[!precisset & precisvalid]
  precisHb[!precisset & precisvalid]<-!is.na(data_complete$Hb[precursor2])[!precisset & precisvalid]
  precisset[!precisset & precisvalid]<-(precisvalid & precisHb)[!precisset & precisvalid]
}
paste("Nr of loops to find a valid Hb:",i)
precursor[idx_noprev_donation]<-precursor2
precursor[idx_noprev_donation][!precisset]<-NA
sum(is.na(precursor))

# also make a column for the previous Hb value
data_complete$prevHb<-NA
data_complete$prevHb<-data_complete$Hb[precursor]
data_complete$prevHb[data_complete$KeyID != data_complete$KeyID[precursor]]<-NA

data_complete$prevHb_prescr<-NA
data_complete$prevHb_prescr<-data_complete$pre_Hb[precursor]
data_complete$prevHb_prescr[data_complete$KeyID != data_complete$KeyID[precursor]]<-NA

# calculate time since last donation
data_complete$dt<-NA
data_complete$dt<-as.numeric(data_complete$DonDate-data_complete$DonDate[precursor])
data_complete$dt[data_complete$KeyID != data_complete$KeyID[precursor]]<-NA

# calculate log10 of the time difference
data_complete$ldt<-log10(data_complete$dt)

# calculate change in Hb from previous donation
data_complete$dHb<-NA
data_complete$dHb<-data_complete$Hb-data_complete$prevHb

data_complete$dHb_prescr <- NA
data_complete$dHb_prescr <- data_complete$pre_Hb-data_complete$prevHb

#association for females:
#select all females with donations above the cut off value
self<-which(data_complete$Sex=="F" & data_complete$dt>=intervaltoshow[1] & data_complete$prevHb>dtf)
linfitf_notdeferred<-lm(dHb~ldt, data=data_complete[self,])
summary(linfitf_notdeferred)
tosave<-append(tosave, list(coeff_notdeferred=c(linfitf_notdeferred$coefficients, length(self), sd(linfitf_notdeferred$residuals))))

#select all females with donations below the cut off value and the capillary follow up measurement
self<-which(data_complete$Sex=="F" & data_complete$dt>=intervaltoshow[1] & data_complete$prevHb<dtf & !is.na(data_complete$pre_Hb))
linfitf_deferred_capillary<-lm(dHb_prescr~ldt, data=data_complete[self,])
summary(linfitf_deferred_capillary)
tosave<-append(tosave, list(coeff_deferred_capillary=c(linfitf_deferred_capillary$coefficients, length(self), sd(linfitf_deferred_capillary$residuals))))

#select all females with donations below the cut off value and the venous follow up measurement
self<-which(data_complete$Sex=="F" & data_complete$dt>=intervaltoshow[1] & data_complete$prevHb<dtf)
linfitf_deferred_venous<-lm(dHb~ldt, data=data_complete[self,])
summary(linfitf_deferred_venous)
tosave<-append(tosave, list(coeff_deferred_venous=c(linfitf_deferred_venous$coefficients, length(self), sd(linfitf_deferred_venous$residuals))))

#associations for males
#select all males with donations above the cut off value
selm<-which(data_complete$Sex=="M" & data_complete$dt>=intervaltoshow[1] & data_complete$prevHb>dtf)
linfitm_notdeferred<-lm(dHb~ldt, data=data_complete[selm,])
summary(linfitm_notdeferred)
tosave<-append(tosave, list(coefm_notdeferred=c(linfitm_notdeferred$coefficients, length(selm), sd(linfitm_notdeferred$residuals))))

#select all males with donations below the cut off value and the capillary follow up measurement
selm<-which(data_complete$Sex=="M" & data_complete$dt>=intervaltoshow[1] & data_complete$prevHb<dtf & !is.na(data_complete$pre_Hb))
linfitm_deferred_capillary<-lm(dHb_prescr~ldt, data=data_complete[selm,])
summary(linfitm_deferred_capillary)
tosave<-append(tosave, list(coefm_deferred_capillary=c(linfitm_deferred_capillary$coefficients, length(selm), sd(linfitm_deferred_capillary$residuals))))

#select all males with donations below the cut off value and the venous follow up measurement
selm<-which(data_complete$Sex=="M" & data_complete$dt>=intervaltoshow[1] & data_complete$prevHb<dtf)
linfitm_deferred_venous<-lm(dHb~ldt, data=data_complete[selm,])
summary(linfitm_deferred_venous)
tosave<-append(tosave, list(coefm_deferred_venous=c(linfitm_deferred_venous$coefficients, length(selm), sd(linfitm_deferred_venous$residuals))))

#save the dataset with the prescreening
data_prescreening <- data_complete[!is.na(data_complete$pre_Hb),]

data_complete$Hbfilled<-data_complete$Hb
data_complete$Hbfilled[is.na(data_complete$Hbfilled)]<-data_complete$pre_Hb[is.na(data_complete$Hbfilled)]
data_complete <- data_complete %>%
  group_by(KeyID) %>%
  mutate(cum_Hb= cumsum(Hbfilled)) %>% ungroup() %>% mutate(meanHb = cum_Hb/numdons)

idx<-1:nrow(data_complete)
precursor<-idx-1
precursor[1]<-nrow(data_complete)

data_complete$prevMeanHb<-NA
data_complete$prevMeanHb<-data_complete$meanHb[precursor]
data_complete$prevMeanHb[data_complete$KeyID != data_complete$KeyID[precursor]]<-NA

data_complete$prevMeanHb[is.na(data_complete$prevMeanHb)] <- data_complete$meanHb[is.na(data_complete$prevMeanHb)]

###########################################################
# Analyse capillary variation 
###########################################################

data_cap<-data_complete[!is.na(data_complete$pre_Hb),]
data_cap$cap_dHb<-data_cap$pre_Hb-data_cap$prevMeanHb
(sd_cap<-sd(data_cap$cap_dHb))
(sd_capf<-sd(data_cap$cap_dHb[data_cap$Sex=="F"]))
(sd_capm<-sd(data_cap$cap_dHb[data_cap$Sex=="M"]))

tosave<-append(tosave, list(sd_cap=c(sd_cap, sd_capf, sd_capm)))

qqPlot(data_cap$cap_dHb, main="Q-Q plot of residuals capillary measurement", ylab="residuals [g/L]" )
qqPlot(data_cap$cap_dHb[data_cap$Sex=="F"], main="Q-Q plot of residuals capillary measurement", ylab="residuals [g/L]" )
qqPlot(data_cap$cap_dHb[data_cap$Sex=="M"], main="Q-Q plot of residuals capillary measurement", ylab="residuals [g/L]" )

quantile(data_cap$cap_dHb, c(.01,.99))
plot(ecdf(data_cap$cap_dHb[data_cap$cap_dHb<75 & data_cap$cap_dHb>-75]))

qqPlot(data_cap$cap_dHb[data_cap$Sex=="F" & data_cap$cap_dHb<75 & data_cap$cap_dHb>-75], main="Q-Q plot of residuals capillary measurement", ylab="residuals [g/L]" )
qqPlot(data_cap$cap_dHb[data_cap$Sex=="M" & data_cap$cap_dHb<75 & data_cap$cap_dHb>-75], main="Q-Q plot of residuals capillary measurement", ylab="residuals [g/L]" )

data_cap <- data_cap[data_cap$cap_dHb<75 & data_cap$cap_dHb>-75,]
(sd_cap<-sd(data_cap$cap_dHb))
(sd_capf<-sd(data_cap$cap_dHb[data_cap$Sex=="F"]))
(sd_capm<-sd(data_cap$cap_dHb[data_cap$Sex=="M"]))

tosave<-append(tosave, list(sd_cap_outliersremoved=c(sd_cap, sd_capf, sd_capm)))

###########################################################
# Create analysis file 
###########################################################
# the analysis file contains per donor
# 1) KeyID
# 2) Sex
# 3) per donor per visit the 
#     Hb level (Hbi), 
#     Mean Hb level (MeanHbi) of all past donations, 
#     Nr of visits (nHbi), and 
#     Nr of measurements that were above the threshold value (HbOki)
# the new dataset is called datt

data <- data_complete[!is.na(data_complete$Hb),]

#make a column for the previous Mean Hb
idx<-1:nrow(data)
precursor<-idx-1
precursor[1]<-nrow(data)

# View(data[data$KeyID==30943,])
# View(data_complete[data_complete$KeyID==30943,])

############################
# Analyse deferrals
############################

# standard deviation of individual venous measurements
(malesd<-sd(linfitm$residuals)/sqrt(2))
(femalesd<-sd(linfitf$residuals)/sqrt(2))

data$d<-qnorm(cutoffperc)*femalesd
data$d[data$Sex=="M"]<-qnorm(cutoffperc)*malesd
# set thresholds per gender
data$th<-dtf
data$th[data$Sex=="M"]<-dtm
table(data$th-data$d, data$Sex)

data_complete$d <- qnorm(cutoffperc)*femalesd
data_complete$d[data_complete$Sex=="M"]<-qnorm(cutoffperc)*malesd
data_complete$th<-dtf
data_complete$th[data_complete$Sex=="M"]<-dtm
table(data_complete$th-data_complete$d, data_complete$Sex)

#for capillary measurements
data_complete$diff_prehb_mean <- data_complete$pre_Hb - data_complete$prevMeanHb

(malessd_pre <- sd(data_complete$diff_prehb_mean[data_complete$Sex=="M"],na.rm=T))
(malessd_pre <- sd_capm)
(femalessd_pre <- sd(data_complete$diff_prehb_mean[data_complete$Sex=="F"],na.rm=T))
(femalessd_pre <- sd_capf)

# of berekend op basis van verschil met vorig meting
sd(data_cap$dHb_prescr, na.rm = T)
sd_capm
(malessd_pre <- sqrt(var(data_cap$dHb_prescr[data_cap$Sex=="M" & !is.na(data_cap$dHb_prescr)])-malesd^2))
sd_capf
(femalessd_pre <- sqrt(var(data_cap$dHb_prescr[data_cap$Sex=="F" & !is.na(data_cap$dHb_prescr)])-femalesd^2))

data_complete$pre_d<-NA
data_complete$pre_d[data_complete$Sex=="M"] <- qnorm(cutoffperc)*malessd_pre
data_complete$pre_d[data_complete$Sex=="F"] <- qnorm(cutoffperc)*femalessd_pre

tosave<-append(tosave, list(malessd_pre=malessd_pre))
tosave<-append(tosave, list(femalessd_pre=femalessd_pre))

###################################################################
# now analyse what the new donor deferral policy would achieve 
# in terms of number of donors now allowed to donate
###################################################################
# define an output array with 8 columns, one row per subsequent donation
# the items stored in each row is explained below
outputsummarytable<-as.data.frame(matrix(0,maxDons,8))
colnames(outputsummarytable)<-c("Deferred", "Non-deferred", "ShouldNotDeferred", "ShouldNotDonate", 
                                "Should_not_have_donated","Missed_by_stopped_donor", "ShouldNotDeferred2", "RequiresReview")

stopped_KeyID <- NA
stopped2_KeyID <- NA
stopafter <- 3
maxDons <- max(data$numdons)

for(i in 1:maxDons){
  calculate <- data[data$numdons==i,]
  # 1 - Deferred donors
  # The number of donors deferred at step i are those with a Hb value that is  
  # below the deferral threshold
  outputsummarytable[i,1] <- sum(!is.na(calculate$Hb) & calculate$Hb<calculate$th)
  
  
  # 2 - non-Deferred donors
  # The number of donors not deferred at step i are those with a Hb value that   
  # is equal or larger than the deferral threshold
  outputsummarytable[i,2] <- sum(!is.na(calculate$Hb) & calculate$Hb >= calculate$th)
  
  # 3 - Donors that should not have been deferred as the Hb deviation relative to
  #     their mean Hb value does not provide sufficient evidence against donation
  # only count events from second donation onwards
  if(i>1) outputsummarytable[i,3] <- sum(!is.na(calculate$Hb) & calculate$Hb < calculate$th & calculate$Hb >= calculate$prevMeanHb-calculate$d, na.rm=T)
  
  # 4 - Donors that should not donate as their Hb is demonstrably below the eligibility threshold
  outputsummarytable[i,4]<-sum(!is.na(calculate$Hb) & calculate$meanHb < calculate$th - calculate$d/sqrt(i), na.rm=T)
  
  # 5 - Donors that should not have donated
  outputsummarytable[i,5] <- sum(!is.na(calculate$Hb) & calculate$prevMeanHb < calculate$th-calculate$d/sqrt(i-1) & calculate$Hb >= calculate$th, na.rm=T)
  
  if(i>stopafter) stopped_KeyID <- c(stopped_KeyID, calculate$KeyID[!is.na(calculate$Hb) & calculate$prevMeanHb < calculate$th - calculate$d/sqrt(i-1)])
  # index stopped2 indicates that the mean Hb level is below the eligibility threshold at previous donation
  if (i>stopafter) stopped2_KeyID <- c(stopped2_KeyID, calculate$KeyID[!is.na(calculate$Hb) & calculate$prevMeanHb < calculate$th])
  
  # 6 - donations missed as a result of new deferral rule
  if(i>stopafter) outputsummarytable[i,6] <- sum(!is.na(calculate$Hb) & calculate$Hb >= calculate$th & calculate$KeyID %in% stopped_KeyID, na.rm=T)
  
  # 7 - Donors that should not have been deferred as the Hb deviation relative to 
  #     the absolute Hb threshold is insufficient (see also evaluation 3 above)
  # this basically presumes that anyone with a Hb level over th-d may donate
  outputsummarytable[i,7]<- sum(!is.na(calculate$Hb) & calculate$Hb< calculate$th & calculate$Hb>=calculate$th-calculate$d, na.rm=T)
  
  # 8 - Identified as outlier, but not deferred
  if(i>1) outputsummarytable[i,8] <- sum(!is.na(calculate$Hb) & calculate$Hb < calculate$prevMeanHb - calculate$d & calculate$Hb >= calculate$th, na.rm=T)
}

length(unique(stopped_KeyID)) 
length(unique(stopped2_KeyID))


length(unique(stopped_KeyID))/length(stopped_KeyID) #proportion of stopped donors
length(unique(stopped2_KeyID))/length(stopped2_KeyID) #proportion of stopped2 donors

outputsummarytable$defprop<-outputsummarytable$Deferred/(outputsummarytable$Deferred+outputsummarytable$`Non-deferred`)
outputsummarytable$nondefprop<-outputsummarytable$ShouldNotDeferred/outputsummarytable$Deferred
outputsummarytable
tosave<-append(tosave, list(outputsummarytable=outputsummarytable))

########################################
# analysis of pre-donation screenings
########################################

outputtableprescreening <- as.data.frame(matrix(0,maxDons,4))
colnames(outputtableprescreening) <- c("Pre-donation screenings conducted", "Deferral during pre-donation screening", "Should not deferred", "Should not donate")

maxDons <- max(data_complete$numdons)

for(i in 1:maxDons){
  calculate <- data_complete[data_complete$numdons==i,]
  #1 - number of pre-screenings conducted
  outputtableprescreening[i,1] <- sum(!is.na(calculate$pre_Hb))
  
  #2 - number of pre-screenings that were a deferral
  outputtableprescreening[i,2] <- sum(calculate$predef==1)
  
  #3 - deferrals at pre-screening that were not necessary, because the pre-screening Hb was within expected variation of the mean
  outputtableprescreening[i,3] <- sum(!is.na(calculate$pre_Hb) & calculate$pre_Hb < calculate$th & calculate$pre_Hb >= calculate$prevMeanHb-calculate$pre_d, na.rm=T)
  
  #4 - donations done after the pre-screening measurment, that should not have been done because the pre-screening Hb was not sufficient
  outputtableprescreening[i,4]<- sum(!is.na(calculate$pre_Hb) & ((calculate$pre_Hb + calculate$prevMeanHb)/i) < calculate$th - calculate$pre_d/sqrt(i), na.rm=T)
}

tosave<-append(tosave, list(outputtableprescreening=outputtableprescreening))


###################################################################
# Write tosave data to datafile
###################################################################
saveRDS(tosave, paste0("SavedDeferralData_",Sys.Date(),".RDS"))

###################################################################
# now calculate various statistics of the updated policy 
###################################################################

# calculate sums per column
sms<-colSums(outputsummarytable)
(totn<-sms[1]+sms[2]) # number of measurements included

sms2<-colSums(outputsummarytable[2:nrow(outputsummarytable),])
(totn2<-sms2[1]+sms2[2]) # number of measurements included
sms2[1]/totn2   # proportion deferred 
sms2[3]/totn2   # proportion that should not have been deferred 
(sms2[1]-sms2[3])/totn2   # proportion with deviant Hb values 
sms2[3]/sms2[1] # proportion unnecessary deferrals within deviation of mean
sms2[7]/totn2   # proportion that Should not have been deferred 
(sms2[1]-sms2[7])/totn2   # proportion that should have been deferred 
sms2[7]/sms2[1] # proportion unnecessary deferrals
sms2[5]         # number of donations 
sms2[5]/totn2   # proportion that should not have donated 
sms2[6]         # number missed by deferred donors
sms2[6]/totn2   # proportion missed by deferred donors 
sms2[8]/totn2   # Reviewed for low relative Hb 

# however, if you only want to start this rule after two donations
sms3<-colSums(outputsummarytable[3:nrow(outputsummarytable),])
(totn3<-sms3[1]+sms3[2])  # number of measurements included
sms3[1]/totn3   # proportion deferred 
sms3[3]/totn3   # proportion that should not have been deferred 
(sms3[1]-sms3[3])/totn3   # proportion with deviant Hb values 
sms3[3]/sms3[1] # proportion unnecessary deferrals within deviation of mean
sms3[7]/totn3   # proportion that should not have been deferred 
(sms3[1]-sms3[7])/totn3   # proportion that should have been deferred 
sms3[7]/sms3[1] # proportion unnecessary deferrals
sms3[5]         # number of donations 
sms3[5]/totn3   # proportion that should not have donated
sms3[6]         # number missed by deferred donors
sms3[6]/totn3   # proportion missed by deferred donors
sms3[8]/totn3   # Reviewed for low relative Hb 



#######################################################
# plot some individual donor profiles
#######################################################
data_complete <- data_complete %>% group_by(KeyID) %>% mutate(totaldon = n())
data <- data %>% group_by(KeyID) %>% mutate(totaldon = n())
# for internal use only
# note that the number of 

maxplots<-3  # USER: Set the maximum number of graphs to plot in row/column of a matrix
# plotdonorprofile(KeyID[250], leg=T, ylim=c(0,190)) # nice illustration with a range of 10 unnecessary deferrals 
# plotdonorprofile(KeyID[250], ylim=c(75,210)) # nice illustration with a range of 10 unnecessary deferrals 

###########################
# Create a selection of donors that should have been deferred at donation 'def' but did donate 
# at least n times
def<-6 # to be set by the USER 
n<-9  # to be set by the USER 
# Nr of donors that fit the criterion
eval(parse(text=paste0("sum(data$numdons==",def, "& data$meanHb + data$d/sqrt(",def,")<data$th & data$Hb > data$th & !is.na(data$Hb) & data$totaldon >",n,", na.rm=T)")))
if (eval(parse(text=paste0("sum(data$numdons==",def, "& data$meanHb + data$d/sqrt(",def,")<data$th & data$Hb > data$th & !is.na(data$Hb)  & data$totaldon >",n,", na.rm=T)")))>0) {
  # set selection of donors
  eval(parse(text=paste0("selID<-data$KeyID[data$numdons==",def, "& data$meanHb + data$d/sqrt(",def,")<data$th & data$Hb > data$th & !is.na(data$Hb)  & data$totaldon >",n,"]")))
  print(selID[!is.na(selID)])
  # plot the donor profile in a matrix
  if(plot_to_pdf) pdf(file=gsub(" ","_",paste0("Deferral at donation ",def," minimum of ",n," donations.pdf")))
  plotmatrix(selID,maxplots=maxplots, ylim=c(70,160),)
  if(plot_to_pdf) dev.off()
}
###########################
# select donors with at least ndef deferrals at donation n and an average Hb level above the deferral threshold
data <- data %>% group_by(KeyID) %>% mutate(def_count = cumsum(def))
data_complete <- data %>% group_by(KeyID) %>% mutate(def_count = cumsum(def))

n<-29     # to be set by the USER 
ndef<-8  # to be set by the USER 
# Nr of donors selected
eval(parse(text=paste0("sum(data$numdons ==",n," & data$def_count >= ndef & data$meanHb > data$th & !is.na(data$Hb))")))
if(eval(parse(text=paste0("sum(data$numdons ==",n," & data$def_count >= ndef & data$meanHb > data$th & !is.na(data$Hb))")))>0){
  # set selection of donors
  eval(parse(text=paste0("selID<-data$KeyID[data$numdons ==",n," & data$def_count >= ndef & data$meanHb > data$th & !is.na(data$Hb)]")))
  print(selID)
  # plot the donor profile in a matrix
  if(plot_to_pdf) pdf(file=gsub(" ","_",paste0("Deferral at ",n,", but with an average Hb level above the threshold.pdf")))
  plotmatrix(selID,maxplots=maxplots)
  if(plot_to_pdf) dev.off()
}
###########################
# Create a selection of donors that were deferred at donation n but have 
# an average Hb level well above (a value "delta") the deferral threshold
n<-30    # to be set by the USER 
delta<-5 # to be set by the USER 
# Nr of donors selected
eval(parse(text=paste0("sum(data$numdons ==",n,"& data$meanHb>data$th+delta & data$Hb<data$th & !is.na(data$Hb))")))
if (eval(parse(text=paste0("sum(data$numdons ==",n,"& data$meanHb>data$th+delta & data$Hb<data$th & !is.na(data$Hb))")))>0){
  # set selection of donors
  eval(parse(text=paste0("selID<-data$KeyID[data$numdons ==",n,"& data$meanHb>data$th+delta & data$Hb<data$th & !is.na(data$Hb)]")))
  print(selID)
  # plot the donor profile in a matrix
  if(plot_to_pdf) pdf(file=gsub(" ","_",paste0("Deferral at ",n," but with an average of ",delta," above the threshold.pdf")))
  plotmatrix(selID, maxplots=maxplots)
  if(plot_to_pdf) dev.off()
  # plot another subset
  # plotmatrix(selID, 2, seedvalue=2)
}

