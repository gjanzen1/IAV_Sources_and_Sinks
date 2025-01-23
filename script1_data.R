
################
### This script reads raw data and updated clade classifications and generates cleaned, formatted datasets
### Datasets:
# data:               all data
# data_nona:          all data, no NAs (excludes samples with incomplete and missing constellations and clades)
# data_uniq:          all unique data from data_nona, with summary statistics
# data_usda:          per-state hog production values
################

setwd("C:/Users/garrett.janzen/OneDrive - USDA/Projects/IAV_Env_Eco")

##########################

library("ggplot2")
library("ggraph")
library("igraph")
library("gplots")
library("ggpmisc")
library("tidygraph")
library("networkD3")
library("viridis")
library("readxl")
library("readr")
library("reshape2")
library("maps")
library("tmap")
library("raster")
library("dismo")
library("sf")
library("spData")
library("magick")
library("grid")
library("data.table")
library("dplyr")
library("stringr")
library("lubridate")
library("network")
library("geodist")
library("pscl")
library("randomForest")
library("Metrics")
# library("reprtree")
library("gbm")
library("minerva")
library("emmeans")
library("zetadiv")
library("TTR")
library("zoo")
library("visNetwork")
library("markovchain")
library("lattice")
library("diagram")
library("ade4")
library("lme4")
library("vegan")

##########################
### Start with a clean memory (not required)
# rm(list = ls())

##########################
### Define some functions

`%!in%` <- Negate(`%in%`)
`%!like%` <- Negate(`%like%`)

eq <- function(x,y) {
  m <- lm(y ~ x)
  as.character(
    as.expression(
      substitute(~~italic(r)^2~"="~r2,
                 list(r2 = format(summary(m)$r.squared, digits = 3)))
    )
  )
}

##########################
### Read in the datasets to be merged.

data1 <- as.data.frame(read_excel("Data/swine-surveillance-data_8_27_24.xlsx"))
data1$Date <- as.Date(data1$Date)

data1$H1 <- str_remove(data1$H1, ".*\\|");table(data1$H1)
data1$H3 <- str_remove(data1$H3, ".*\\|");table(data1$H3)
data1$N1 <- str_remove(data1$N1, ".*\\|");table(data1$N1)
data1$N2 <- str_remove(data1$N2, ".*\\|");table(data1$N2)

##########################

data2 <- as.data.frame(read.table(file = 'Data/000-query-result-v3.tsv', sep = '\t', header = TRUE))
data2 <- data2[,c(1:16)] #drops columns for HA sequence and NA sequence
colnames(data2) <- c("Barcode","Strain","Host","Subtype","Year","Month","Day","Country","State","Source",
                     "H_Genbank","N_Genbank","Constellation","Ha_clade","Na_clade","GL_Clade");dim(data2);head(data2)
data2$Date <- as.Date(paste(data2$Year, data2$Month, data2$Day, sep="-"))

data2 <- data2[,c("Barcode","Strain","Host","Subtype","Year","Month","Day","Country","State","Source",
                  "H_Genbank","N_Genbank","Constellation","Na_clade")]

##########################
### Modify dataset structure before they are merged.

data2_reclassifications <- as.data.frame(read.table(file = 'Data/000-query-result-v3--octoFLU-classify.txt', sep = '\t', header = TRUE));head(data2_reclassifications);dim(data2_reclassifications)
colnames(data2_reclassifications) <- c("Strain", "Barcode","segment_name","space","segment_subtype","Clade_rec","Gl_clade_rec")
data2_reclassifications$segment_name <- ifelse(is.na(data2_reclassifications$segment_name), "NA", data2_reclassifications$segment_name);table(data2_reclassifications$segment_name)
data2_rec_HA <- data2_reclassifications[which(data2_reclassifications$segment_name == "HA"),];dim(data2_rec_HA)
data2_rec_NA <- data2_reclassifications[which(data2_reclassifications$segment_name == "NA"),];dim(data2_rec_NA)

data2_rec_HA <- data2_rec_HA[,c("Strain","Clade_rec","Gl_clade_rec")];dim(data2_rec_HA)
data2_rec_NA <- data2_rec_NA[,c("Strain","Clade_rec","Gl_clade_rec")];dim(data2_rec_NA)

data2_HA <- data2[which(data2$N_Genbank == ""),];dim(data2_HA)
data2_NA <- data2[which(data2$N_Genbank != ""),];dim(data2_NA)

data2m_HA <- merge(data2_HA, data2_rec_HA, by="Strain", all.x=TRUE, all.y=FALSE);head(data2m_HA);dim(data2m_HA)
data2m_NA <- merge(data2_NA, data2_rec_NA, by="Strain", all.x=TRUE, all.y=FALSE);head(data2m_NA);dim(data2m_NA)

colnames(data2m_HA)[15:16] <- c("HA_clade_rec", "HA_gl_clade_rec")

data2m_NA <- data2m_NA[,c("Strain","Clade_rec","Gl_clade_rec")]
colnames(data2m_NA) <- c("Strain","NA_clade_rec","NA_gl_clade_rec")

data2m <- merge(data2m_HA, data2m_NA, by="Strain", all.x=TRUE, all.y=FALSE);head(data2m);dim(data2m)
data2 <- data2m

data2$Na_clade <- NULL

colnames(data2)[14:17] <- c("Ha_clade_US","Ha_clade","Na_clade_US","Na_clade")

data2$N1 <- NA
data2$N2 <- NA
data2$H1 <- NA
data2$H3 <- NA
data2$N1 <- ifelse(substr(data2$Subtype, 3, 4) == "N1", data2$Na_clade, data2$N1)
data2$N2 <- ifelse(substr(data2$Subtype, 3, 4) == "N2", data2$Na_clade, data2$N2)
data2$H1 <- ifelse(substr(data2$Subtype, 1, 2) == "H1", data2$Ha_clade, data2$H1)
data2$H3 <- ifelse(substr(data2$Subtype, 1, 2) == "H3", data2$Ha_clade, data2$H3)
# data2$GL_Clade <- NULL
# data2$Ha_gl_clade <- NULL
# data2$US_Clade <- NULL
# data2$Ha_clade <- NULL
# data2$Na_clade <- NULL
data2$Haseq <- NULL
data2$Naseq <- NULL
data2$Host <- NULL
data2$Source <- NULL

# data2 does not have columns corresponding to the 8 genes, but we can pull that data from the Constellation column
# Genes: H, N, PB2, PB1, PA, NP, M, NS
data2$PB2 <- substr(data2$Constellation, 1, 1)
data2$PB1 <- substr(data2$Constellation, 2, 2)
data2$PA <- substr(data2$Constellation, 3, 3)
data2$NP <- substr(data2$Constellation, 4, 4)
data2$M <- substr(data2$Constellation, 5, 5)
data2$NS <- substr(data2$Constellation, 6, 6)
data2$PB2 <- ifelse(data2$PB2 == "T", "TRIG", data2$PB2)
data2$PB2 <- ifelse(data2$PB2 == "P", "pdm", data2$PB2)
data2$PB2 <- ifelse(data2$PB2 == "V", "LAIV", data2$PB2)
# data2$PB2 <- ifelse(data2$PB2 == "X", "humanSeasonal", data2$PB2)
data2$PB1 <- ifelse(data2$PB1 == "T", "TRIG", data2$PB1)
data2$PB1 <- ifelse(data2$PB1 == "P", "pdm", data2$PB1)
data2$PB1 <- ifelse(data2$PB1 == "V", "LAIV", data2$PB1)
# data2$PB1 <- ifelse(data2$PB1 == "X", "humanSeasonal", data2$PB1)
data2$PA <- ifelse(data2$PA == "T", "TRIG", data2$PA)
data2$PA <- ifelse(data2$PA == "P", "pdm", data2$PA)
data2$PA <- ifelse(data2$PA == "V", "LAIV", data2$PA)
# data2$PA <- ifelse(data2$PA == "X", "humanSeasonal", data2$PA)
data2$NP <- ifelse(data2$NP == "T", "TRIG", data2$NP)
data2$NP <- ifelse(data2$NP == "P", "pdm", data2$NP)
data2$NP <- ifelse(data2$NP == "V", "LAIV", data2$NP)
# data2$NP <- ifelse(data2$NP == "X", "humanSeasonal", data2$NP)
data2$M <- ifelse(data2$M == "T", "TRIG", data2$M)
data2$M <- ifelse(data2$M == "P", "pdm", data2$M)
data2$M <- ifelse(data2$M == "V", "LAIV", data2$M)
# data2$M <- ifelse(data2$M == "X", "humanSeasonal", data2$M)
data2$NS <- ifelse(data2$NS == "T", "TRIG", data2$NS)
data2$NS <- ifelse(data2$NS == "P", "pdm", data2$NS)
data2$NS <- ifelse(data2$NS == "V", "LAIV", data2$NS)
# data2$NS <- ifelse(data2$NS == "X", "humanSeasonal", data2$NS)

i <- NULL
for(i in 1:length(data2$State)){
  #i <- 1
  data2$State[i] <- state.abb[grep(data2$State[i], state.name)]
}

data2$Date <- as.Date(paste(data2$Year,data2$Month,data2$Day,sep="-"));head(data2$Date)

##########################
###data3
data3 <- as.data.frame(read.table(file = 'Data/2023-11-09.tsv', sep = '\t', header = FALSE));dim(data3)
data3 <- data3[,c(1:16)] #drops columns for HA sequence and NA sequence

colnames(data3) <- c("Barcode","Strain","Host","Subtype","Year","Month","Day","Country","State","Source",
                     "H_Genbank","N_Genbank","Constellation","Ha_clade","Na_clade","GL_Clade");dim(data3);head(data3)
data3$Date <- as.Date(paste(data3$Year, data3$Month, data3$Day, sep="-"))

data3 <- data3[,c("Barcode","Strain","Host","Subtype","Year","Month","Day","Country","State","Source",
                  "H_Genbank","N_Genbank","Constellation","Na_clade")]

data3_reclassifications <- as.data.frame(read.table(file = 'Data/2023-11-09--octoFLU-classify.txt', sep = '\t', header = TRUE));head(data3_reclassifications);dim(data3_reclassifications)
colnames(data3_reclassifications) <- c("Strain", "Barcode","segment_name","space","segment_subtype","Clade_rec","Gl_clade_rec")
data3_reclassifications$segment_name <- ifelse(is.na(data3_reclassifications$segment_name), "NA", data3_reclassifications$segment_name);table(data3_reclassifications$segment_name)
data3_rec_HA <- data3_reclassifications[which(data3_reclassifications$segment_name == "HA"),];dim(data3_rec_HA)
data3_rec_NA <- data3_reclassifications[which(data3_reclassifications$segment_name == "NA"),];dim(data3_rec_NA)

data3_rec_HA <- data3_rec_HA[,c("Strain","Clade_rec","Gl_clade_rec")];dim(data3_rec_HA)
data3_rec_NA <- data3_rec_NA[,c("Strain","Clade_rec","Gl_clade_rec")];dim(data3_rec_NA)

data3_HA <- data3[which(data3$N_Genbank == ""),];dim(data3_HA)
data3_NA <- data3[which(data3$N_Genbank != ""),];dim(data3_NA)

data3m_HA <- merge(data3_HA, data3_rec_HA, by="Strain", all.x=TRUE, all.y=FALSE);head(data3m_HA);dim(data3m_HA)
data3m_NA <- merge(data3_NA, data3_rec_NA, by="Strain", all.x=TRUE, all.y=FALSE);head(data3m_NA);dim(data3m_NA)

colnames(data3m_HA)[15:16] <- c("HA_clade_rec", "HA_gl_clade_rec")

data3m_NA <- data3m_NA[,c("Strain","Clade_rec","Gl_clade_rec")]
colnames(data3m_NA) <- c("Strain","NA_clade_rec","NA_gl_clade_rec")

data3m <- merge(data3m_HA, data3m_NA, by="Strain", all.x=TRUE, all.y=FALSE);head(data3m);dim(data3m)
data3 <- data3m

data3$Na_clade <- NULL

colnames(data3)[14:17] <- c("Ha_clade_US","Ha_clade","Na_clade_US","Na_clade")

data3$N1 <- NA
data3$N2 <- NA
data3$H1 <- NA
data3$H3 <- NA
data3$N1 <- ifelse(substr(data3$Subtype, 3, 4) == "N1", data3$Na_clade, data3$N1)
data3$N2 <- ifelse(substr(data3$Subtype, 3, 4) == "N2", data3$Na_clade, data3$N2)
data3$H1 <- ifelse(substr(data3$Subtype, 1, 2) == "H1", data3$Ha_clade, data3$H1)
data3$H3 <- ifelse(substr(data3$Subtype, 1, 2) == "H3", data3$Ha_clade, data3$H3)
# data3$GL_Clade <- NULL
# data3$Ha_gl_clade <- NULL
# data3$US_Clade <- NULL
# data3$Ha_clade <- NULL
# data3$Na_clade <- NULL
data3$Haseq <- NULL
data3$Naseq <- NULL
data3$Host <- NULL
data3$Source <- NULL

# data3 does not have columns corresponding to the 8 genes, but we can pull that data from the Constellation column
# Genes: H, N, PB2, PB1, PA, NP, M, NS
data3$PB2 <- substr(data3$Constellation, 1, 1)
data3$PB1 <- substr(data3$Constellation, 2, 2)
data3$PA <- substr(data3$Constellation, 3, 3)
data3$NP <- substr(data3$Constellation, 4, 4)
data3$M <- substr(data3$Constellation, 5, 5)
data3$NS <- substr(data3$Constellation, 6, 6)
data3$PB2 <- ifelse(data3$PB2 == "T", "TRIG", data3$PB2)
data3$PB2 <- ifelse(data3$PB2 == "P", "pdm", data3$PB2)
data3$PB2 <- ifelse(data3$PB2 == "V", "LAIV", data3$PB2)
# data3$PB2 <- ifelse(data3$PB2 == "X", "humanSeasonal", data3$PB2)
data3$PB1 <- ifelse(data3$PB1 == "T", "TRIG", data3$PB1)
data3$PB1 <- ifelse(data3$PB1 == "P", "pdm", data3$PB1)
data3$PB1 <- ifelse(data3$PB1 == "V", "LAIV", data3$PB1)
# data3$PB1 <- ifelse(data3$PB1 == "X", "humanSeasonal", data3$PB1)
data3$PA <- ifelse(data3$PA == "T", "TRIG", data3$PA)
data3$PA <- ifelse(data3$PA == "P", "pdm", data3$PA)
data3$PA <- ifelse(data3$PA == "V", "LAIV", data3$PA)
# data3$PA <- ifelse(data3$PA == "X", "humanSeasonal", data3$PA)
data3$NP <- ifelse(data3$NP == "T", "TRIG", data3$NP)
data3$NP <- ifelse(data3$NP == "P", "pdm", data3$NP)
data3$NP <- ifelse(data3$NP == "V", "LAIV", data3$NP)
# data3$NP <- ifelse(data3$NP == "X", "humanSeasonal", data3$NP)
data3$M <- ifelse(data3$M == "T", "TRIG", data3$M)
data3$M <- ifelse(data3$M == "P", "pdm", data3$M)
data3$M <- ifelse(data3$M == "V", "LAIV", data3$M)
# data3$M <- ifelse(data3$M == "X", "humanSeasonal", data3$M)
data3$NS <- ifelse(data3$NS == "T", "TRIG", data3$NS)
data3$NS <- ifelse(data3$NS == "P", "pdm", data3$NS)
data3$NS <- ifelse(data3$NS == "V", "LAIV", data3$NS)
# data3$NS <- ifelse(data3$NS == "X", "humanSeasonal", data3$NS)

i <- NULL
for(i in 1:length(data3$State)){
  #i <- 1
  data3$State[i] <- state.abb[grep(data3$State[i], state.name)]
}

# data3 <- data3[which(data3$N_Genbank != ""),]
# data3 <- data3[which(data3$Constellation != ""),]

sum(!is.na(data3$H1) & !is.na(data3$H3)) # number of rows that have both H1 and H3 values
sum(is.na(data3$H1) & is.na(data3$H3)) # number of rows that have neither H1 or H3 values
sum(!is.na(data3$N1) & !is.na(data3$N2)) # number of rows that have both N1 and N2 values
sum(is.na(data3$N1) & is.na(data3$N2)) # number of rows that have neither N1 or N2 values
sum(is.na(data3$N1) & is.na(data3$N2) & is.na(data3$H1) & is.na(data3$H3)) # number of rows that lack all four

### Rectify Ha clades (column name and formatting)
colnames(data1)[which(colnames(data1) == "GL_Clade")] <-    "Ha_clade"
colnames(data2)[which(colnames(data2) == "Ha_gl_clade")] <- "Ha_clade"
colnames(data3)[which(colnames(data3) == "Ha_clade")] <-    "Ha_clade"
# Some Ha clades have "3." appended to the string. This drops those characters.
data1$Ha_clade <- ifelse(substr(data1$Ha_clade, 1, 1)=="3", str_sub(data1$Ha_clade, 3), data1$Ha_clade);table(data1$Ha_clade)
data2$Ha_clade <- ifelse(substr(data2$Ha_clade, 1, 1)=="3", str_sub(data2$Ha_clade, 3), data2$Ha_clade);table(data2$Ha_clade)
data3$Ha_clade <- ifelse(substr(data3$Ha_clade, 1, 1)=="3", str_sub(data3$Ha_clade, 3), data3$Ha_clade);table(data3$Ha_clade)

data3$Date <- as.Date(paste(data3$Year,data3$Month,data3$Day,sep="-"));head(data3$Date)

##########################
# for troubleshooting:
data1$file <- "data1"
data2$file <- "data2"
data3$file <- "data3"

head(data1);head(data2);head(data3)
dim(data1);dim(data2);dim(data3)
data1$Date <- as.Date(data1$Date)
data2$Date <- as.Date(data2$Date)
data3$Date <- as.Date(data3$Date)

data2[which(data2$Strain %!in% data1$Strain),]
data3[which(data3$Strain %!in% data1$Strain),]
data3[which(data3$Strain %!in% data2$Strain),]

##############################
##############################
##############################
# Merge the datasets here.
data <- bind_rows(data1, data2, data3)
head(data);dim(data)
data <- data[order(data$Date),]


##############################
##############################
##############################
##############################
### Begin cleaning and formatting the merged dataset.

data$Ha_clade <- ifelse(data$Ha_clade == "Other-Human-2020", "Other-Human-2020", data$Ha_clade)
data$Ha_clade <- ifelse(data$Ha_clade == "1A.2-3-like", "1A.4", data$Ha_clade)
data$Ha_clade <- ifelse(data$Ha_clade == "1A.3.3.2-vaccine", "1A.3.3.2", data$Ha_clade)
data$Ha_clade <- ifelse(data$Ha_clade == "1990.4.1", "1990.4.a", data$Ha_clade)
data$Ha_clade <- ifelse(data$Ha_clade == "1990.4.2", "1990.4.b", data$Ha_clade)
data$Ha_clade <- ifelse(data$Ha_clade == "1990.4.3", "1990.4.c", data$Ha_clade)
data$Ha_clade <- ifelse(data$Ha_clade == "1990.4.4", "1990.4.d", data$Ha_clade)
data$Ha_clade <- ifelse(data$Ha_clade == "1990.4.5", "1990.4.e", data$Ha_clade)
data$Ha_clade <- ifelse(data$Ha_clade == "1990.4.6", "1990.4.f", data$Ha_clade)
data$Ha_clade <- ifelse(data$Ha_clade == "1990.4.7", "1990.4.g", data$Ha_clade)
data$Ha_clade <- ifelse(data$Ha_clade == "1990.4.8", "1990.4.h", data$Ha_clade)
data$Ha_clade <- ifelse(data$Ha_clade == "1990.4.9", "1990.4.i", data$Ha_clade)
data$Ha_clade <- ifelse(data$Ha_clade == "1990.4.10", "1990.4.j", data$Ha_clade)
data$Ha_clade <- ifelse(data$Ha_clade == "1990.4.11", "1990.4.k", data$Ha_clade)
data$Ha_clade <- ifelse(data$Ha_clade == "1990.4.12", "1990.4.l", data$Ha_clade)

data$Year <- ifelse(is.na(data$Year), substr(data$Date, 1, 4), data$Year)
data <- data[order(as.Date(data$Date, format="%Y/%m/%d")),]
data_dmin <- as.Date(min(data$Date))
data_dmax <- as.Date(max(data$Date))

data$Date_format <- format(data$Date, format="%m-%d")
data$Date_month <- NA
data$Date_month <- ifelse(substr(data$Date_format, start=1, stop=2) == "01", "January", data$Date_month)
data$Date_month <- ifelse(substr(data$Date_format, start=1, stop=2) == "02", "February", data$Date_month)
data$Date_month <- ifelse(substr(data$Date_format, start=1, stop=2) == "03", "March", data$Date_month)
data$Date_month <- ifelse(substr(data$Date_format, start=1, stop=2) == "04", "April", data$Date_month)
data$Date_month <- ifelse(substr(data$Date_format, start=1, stop=2) == "05", "May", data$Date_month)
data$Date_month <- ifelse(substr(data$Date_format, start=1, stop=2) == "06", "June", data$Date_month)
data$Date_month <- ifelse(substr(data$Date_format, start=1, stop=2) == "07", "July", data$Date_month)
data$Date_month <- ifelse(substr(data$Date_format, start=1, stop=2) == "08", "August", data$Date_month)
data$Date_month <- ifelse(substr(data$Date_format, start=1, stop=2) == "09", "September", data$Date_month)
data$Date_month <- ifelse(substr(data$Date_format, start=1, stop=2) == "10", "October", data$Date_month)
data$Date_month <- ifelse(substr(data$Date_format, start=1, stop=2) == "11", "November", data$Date_month)
data$Date_month <- ifelse(substr(data$Date_format, start=1, stop=2) == "12", "December", data$Date_month)
data$Year <- as.character(data$Year)
data$Month <- NULL; data$Day <- NULL
names(data)[names(data) == 'Year'] <- 'Date_year'

data$Date_season <- NA
data$Date_season <- ifelse(data$Date_month == "December" |
                             data$Date_month == "January" |
                             data$Date_month == "February", "winter", data$Date_season)
data$Date_season <- ifelse(data$Date_month == "March" |
                             data$Date_month == "April" |
                             data$Date_month == "May", "spring", data$Date_season)
data$Date_season <- ifelse(data$Date_month == "June" |
                             data$Date_month == "July" |
                             data$Date_month == "August", "summer", data$Date_season)
data$Date_season <- ifelse(data$Date_month == "September" |
                             data$Date_month == "October" |
                             data$Date_month == "November", "fall", data$Date_season)

##############################
##############################
##############################
### Getting rid of rows of duplicated observations.
#Strain names incorporate place, time, and subtype of the virus, and include a barcode sequence.
#Each strain name should appear once in the dataset. While constellation/clade combinations can recur,
#Each strain is a single isolate. Therefore, we can drop all rows with duplicated strain names.

repeated_strains <- names(table(data$Strain)[which(table(data$Strain) > 1)]);length(repeated_strains)
data_temp <- data[which(data$Strain %in% repeated_strains),]
data_keepers <- data[which(data$Strain %!in% repeated_strains),]

i <- NULL
for(i in 1:length(repeated_strains)){
  # i <- 1
  dat <- data_temp[which(data_temp$Strain == repeated_strains[i]),]
  if(length(unique(rowSums(!is.na(dat)))) < nrow(dat)){ #if there are instances where there are fewer row sums than rows, then there are ties.
    # print(paste0("i = ",i,", ties present, dropping missing Constellations and rechecking."))
    # print(paste0(rowSums(!is.na(dat))))
    dat <- dat[which(!is.na(dat$Constellation)),]
    if(length(unique(rowSums(!is.na(dat)))) < nrow(dat)){ #if there are instances where there are fewer row sums than rows, then there are ties.
      # print("Ties are still present. Dropping missing Ha and Na clades and rechecking.")
      dat <- dat[which(!is.na(dat$Ha_clade)),]
      dat <- dat[which(!is.na(dat$Na_clade)),]
    }
    # print("Exiting tie loop.")
  }
  dat <- dat[which.max(rowSums(!is.na(dat))),]
  data_keepers <- rbind(data_keepers, dat)
}

data <- data_keepers
data_keepers <- NULL
data <- data[order(as.Date(data$Date, format="%Y/%m/%d")),]

#########################################################
### Manually fixing some errors that can't be done procedurally.

#There are a few strange cases where a barcode is applied to two samples.
#See the conversation with Blake Inderski on Teams on 11/15/2023
# data <- data[which(data$Strain != "A/swine/Indiana/A02750668/2022"),]
table(table(data$Barcode)) 
tbl <- table(data$Barcode)[order(table(data$Barcode), decreasing=TRUE)];head(tbl)# not all are 1!
repeated_barcodes <- names(tbl)[2:4];repeated_barcodes # these are the ones we can fix
temp <- data[which(data$Barcode %in% repeated_barcodes),];temp
data <- data[which(data$Barcode %!in% repeated_barcodes),] #removal of these data, we add them back in after selecting desired rows

temp1 <- temp[which(temp$Barcode == "A02711864"),] #first one has typo, lowercase A in barcode, BUT also includes Constellation.
temp1$Constellation <- temp1$Constellation[1]
temp1 <- temp1[2,]
temp2 <- temp[which(temp$Barcode == "A02750668"),] #According to GenBank, the Indiana instance is the one to retain.
temp2 <- temp2[1,]
temp3 <- temp[which(temp$Barcode == "A02978463"),] #first one has typos, including doubled strain name, BUT also includes Constellation.
temp3$Constellation <- temp3$Constellation[1]
temp3 <- temp3[2,]

temp <- rbind(temp1,temp2,temp3)
data <- rbind(data, temp)
rm(temp, tbl, repeated_barcodes)

#These strains are erroneous, and are here fixed manually.
#See the conversation with Tavis Anderson on Teams on 11/22/2023
data[which(data$Strain == "A/swine/Tennessee/A01785435/2018"),]
data$Subtype <- ifelse(data$Strain == "A/swine/Tennessee/A01785435/2018",  "H1N1", data$Subtype)
data$Na_clade <- ifelse(data$Strain == "A/swine/Tennessee/A01785435/2018", "N1.C.3.2", data$Na_clade)
data$N1 <- ifelse(data$Strain == "A/swine/Tennessee/A01785435/2018",       "N1.C.3.2", data$N1)
data$N2 <- ifelse(data$Strain == "A/swine/Tennessee/A01785435/2018",       NA, data$N2)
data[which(data$Strain == "A/swine/Minnesota/A02245867/2021"),]
data$Subtype <- ifelse(data$Strain == "A/swine/Minnesota/A02245867/2021", "H3N2", data$Subtype)
data$H1 <- ifelse(data$Strain == "A/swine/Minnesota/A02245867/2021",      NA, data$H1)
data$H3 <- ifelse(data$Strain == "A/swine/Minnesota/A02245867/2021",      "1990.4.a", data$H3)
#These strains have a "mixed" subtype, which we can solve by looking at H1 and N1 columns.
data[which(data$Subtype == "mixed"),]
data$Subtype <- ifelse(data$Strain == "A/swine/Minnesota/A01244318/2012",  "H1N1", data$Subtype)
data$Subtype <- ifelse(data$Strain == "A/swine/Kansas/A01377621/2015",     "H1N1", data$Subtype)
#These strains have a "mixed" constellation or has "A", which we cannot solve, so set to ""
data[which(data$Constellation == "mixed"),]
i <- NULL
for (i in 1:nrow(data)){
  dat <- data[i,]
  if(!is.na(dat$Constellation)){
    if(dat$Constellation == "mixed" | dat$Constellation == "AAAAAA"){
      dat$PB2 <- dat$PB1 <- dat$PA <- dat$NP <- dat$M <- dat$Constellation <- ""
    }
  }
  data[i,] <- dat
}
# data[which(data$Constellation == "mixed"),]
# data[which(data$Constellation == "AAAAAA"),]

#####################
#These strains have two subtypes or two strains, we sort this out too.

temp <- data[grep(",", data$Strain), ]
data <- data[-grep(",", data$Strain), ] #taking out problem rows, to be added back in when fixed
setDT(temp)[, paste0("Strain",1:2) := tstrsplit(Strain, ",")]
temp$StrainTest <- "No"

for(i in 1:nrow(temp)){
  # i <- 1
  straini1 <- temp[i,]$Strain1
  straini2 <- temp[i,]$Strain2
  if(straini1 == straini2){
    temp[i,]$StrainTest == "Match"
  }
}

#As none match, we need to pick winners.
selection <- c(2, #spelling
               1,1,1,1,1,1, #state, rather than country
               2,2,2, #spelling of Illinois              #10
               1, #year AND spelling, will fixt after the for-loop
               1, #year
               1,1,1, #year #15
               1,2,2,1,1,2, #State spelling #21
               1 #year
)

for(i in 1:nrow(temp)){
  # i <- 1
  sel <- selection[i]
  temp$Strain[i] <- ifelse(sel == 1, temp$Strain1[i], temp$Strain2[i])
  temp$Selection[i] <- sel
}

temp$Strain <- ifelse(temp$Strain == "A/swine/Ilinois/A01240775/2011", "A/swine/Illinois/A01240775/2011", temp$Strain)
# temp

temp$Strain1 <- temp$Strain2 <- temp$StrainTest <- temp$Selection <- NULL
# ncol(data);ncol(temp)

data <- rbind(data, temp)
rm(temp, selection)

####################

#These constellations have two subtypes or two constellations, we sort this out too.
temp <- data[grep(",", data$Constellation), ]
data <- data[-grep(",", data$Constellation), ] #taking out problem rows, to be added back in when fixed
setDT(temp)[, paste0("Constellation",1:2) := tstrsplit(Constellation, ",")]
temp$ConstellationTest <- "No"

for(i in 1:nrow(temp)){
  # i <- 1
  constellationi1 <- temp[i,]$Constellation1
  constellationi2 <- temp[i,]$Constellation2
  if(constellationi1 == constellationi2){
    temp[i,]$ConstellationTest == "Match"
  }
}

#Rather than picking winners, we will blend the constellations for strains with two constellations:
temp$PB2 <- temp$PB1 <- temp$PA <- temp$NP <- temp$M <- temp$NS <- ""
temp$Constellation_merge <- "------"

for(i in 1:nrow(temp)){        #for each row of temp,
  # i <- 1
  Conm <- ""
  for(j in 1:6){                #for each gene in Constellation 1,
    # j <- 1
    Con1genej <- substr(temp$Constellation1[i], j, j)
    Con2genej <- substr(temp$Constellation2[i], j, j)
    Conmgenej <- "-"
    Conmgenej <- ifelse(Con1genej != "-", Con1genej, Conmgenej)
    Conmgenej <- ifelse(Con2genej != "-", Con2genej, Conmgenej)
    Conm <- paste0(Conm, Conmgenej)
  }
  temp$Constellation_merge[i] <- Conm
}

temp$Constellation <- temp$Constellation_merge

temp$Constellation1 <- temp$Constellation2 <- temp$ConstellationTest <- temp$Constellation_merge <- NULL

ncol(data);ncol(temp)
data <- rbind(data, temp)
rm(temp, Con1genej, Con2genej, Conm, Conmgenej, constellationi1, constellationi2, straini1, straini2)

#######
# Some constellations contain an X or A character, rather than the appropriate -, or are "NA"/NA. We fix that here.
data$Constellation <- gsub('X', '-', data$Constellation);table(data$Constellation)
data$Constellation <- gsub('A', '-', data$Constellation);table(data$Constellation)
data$Constellation <- ifelse(data$Constellation == "NA" | is.na(data$Constellation) | data$Constellation == "",
                             "------", data$Constellation)
data$PB2 <- ifelse(data$PB2 == "A" | data$PB2 == "", "-", data$PB2);table(data$PB2)
data$PB1 <- ifelse(data$PB1 == "A" | data$PB1 == "", "-", data$PB1);table(data$PB1)
data$PA <- ifelse(data$PA == "A" | data$PA == "", "-", data$PA);table(data$PA)
data$NP <- ifelse(data$NP == "A" | data$NP == "", "-", data$NP);table(data$NP)
data$M <- ifelse(data$M == "A" | data$M == "", "-", data$M);table(data$M)
data$NS <- ifelse(data$NS == "A" | data$NS == "", "-", data$NS);table(data$NS)
table(data$Constellation)
###################
###################
###################

data$H1 <- ifelse(data$H1 == "PDM", "pdm", data$H1)
data$H1 <- ifelse(data$H1 == "pdm-vaccine", "pdm", data$H1)

data$H3 <- ifelse(data$H3 == "2010-human like", "2010 human-like", data$H3)
data$H3 <- ifelse(data$H3 == "2016-human like", "2016 human-like", data$H3)
data$H3 <- ifelse(data$H3 == "Other-Human_2020", "2020 human-like", data$H3)
data$N1 <- ifelse(data$N1 == "PDM", "pdm", data$N1)

# #A H4N6 is included as an outgroup. We should remove it. #UPDATE: It's an outgroup, but still a real sample, so retain it.
# data <- data[data$GL_Clade != "Outgroup",]

### This has to be re-run, some of the internal gene columns are still deficient.
data$PB2 <- substr(data$Constellation, 1, 1)
data$PB1 <- substr(data$Constellation, 2, 2)
data$PA <- substr(data$Constellation, 3, 3)
data$NP <- substr(data$Constellation, 4, 4)
data$M <- substr(data$Constellation, 5, 5)
data$NS <- substr(data$Constellation, 6, 6)
data$PB2 <- ifelse(data$PB2 == "T", "TRIG", data$PB2)
data$PB2 <- ifelse(data$PB2 == "P", "pdm", data$PB2)
data$PB2 <- ifelse(data$PB2 == "V", "LAIV", data$PB2)
# data$PB2 <- ifelse(data$PB2 == "X", "humanSeasonal", data$PB2)
data$PB1 <- ifelse(data$PB1 == "T", "TRIG", data$PB1)
data$PB1 <- ifelse(data$PB1 == "P", "pdm", data$PB1)
data$PB1 <- ifelse(data$PB1 == "V", "LAIV", data$PB1)
# data$PB1 <- ifelse(data$PB1 == "X", "humanSeasonal", data$PB1)
data$PA <- ifelse(data$PA == "T", "TRIG", data$PA)
data$PA <- ifelse(data$PA == "P", "pdm", data$PA)
data$PA <- ifelse(data$PA == "V", "LAIV", data$PA)
# data$PA <- ifelse(data$PA == "X", "humanSeasonal", data$PA)
data$NP <- ifelse(data$NP == "T", "TRIG", data$NP)
data$NP <- ifelse(data$NP == "P", "pdm", data$NP)
data$NP <- ifelse(data$NP == "V", "LAIV", data$NP)
# data$NP <- ifelse(data$NP == "X", "humanSeasonal", data$NP)
data$M <- ifelse(data$M == "T", "TRIG", data$M)
data$M <- ifelse(data$M == "P", "pdm", data$M)
data$M <- ifelse(data$M == "V", "LAIV", data$M)
# data$M <- ifelse(data$M == "X", "humanSeasonal", data$M)
data$NS <- ifelse(data$NS == "T", "TRIG", data$NS)
data$NS <- ifelse(data$NS == "P", "pdm", data$NS)
data$NS <- ifelse(data$NS == "V", "LAIV", data$NS)
# data$NS <- ifelse(data$NS == "X", "humanSeasonal", data$NS)

# There may be cases where we can re-derive clades from dedicated columns
table(data$Na_clade);table(is.na(data$Na_clade))
data$N1 <- ifelse(is.na(data$N1), "", data$N1);head(data$N1)
data$N2 <- ifelse(is.na(data$N2), "", data$N2);head(data$N2)
i <- iclade <- dat <- NULL
for(i in 1:nrow(data)){
  dat <- data[i,]
  iclade <- dat$Na_clade
  iclade <- ifelse(is.na(iclade), paste0(dat$N1, dat$N2), iclade)
  dat$Na_clade <- iclade
  data[i,] <- dat
}

table(data$Ha_clade);table(is.na(data$Ha_clade))
data$H1 <- ifelse(is.na(data$H1), "", data$H1);head(data$H1)
data$H3 <- ifelse(is.na(data$H3), "", data$H3);head(data$H3)
i <- iclade <- dat <- NULL
for(i in 1:nrow(data)){
  dat <- data[i,]
  iclade <- dat$Ha_clade
  iclade <- ifelse(is.na(iclade), paste0(dat$H1, dat$H3), iclade)
  dat$Ha_clade <- iclade
  data[i,] <- dat
}

#################################

data$H_simple <- substr(data$Subtype, 1, 2)
data$H_complex <- NA
data$H_complex <- ifelse(data$H_simple == "H1", paste("H1-", data$Ha_clade, sep=""), data$H_complex)
data$H_complex <- ifelse(data$H_simple == "H3", paste("H3-", data$Ha_clade, sep=""), data$H_complex)
data$H_complex <- ifelse(data$H_simple == "H4", paste("H4-", data$Ha_clade, sep=""), data$H_complex)
data$H_complex <- ifelse(substring(data$H_complex,nchar(data$H_complex)-2+1,nchar(data$H_complex))=="NA", 
                         substr(data$H_complex,1,nchar(data$H_complex)-2),
                         data$H_complex)

data$Na_clade <- ifelse(substr(data$Na_clade, 1, 2) == "N1", 
                        substr(data$Na_clade, 2, nchar(data$Na_clade)), 
                        data$Na_clade)
data$N_simple <- substr(data$Subtype, 3, 4)
data$N_complex <- NA
data$N_complex <- ifelse(data$N_simple == "N1", paste("N1-", data$Na_clade, sep=""), data$N_complex)
data$N_complex <- ifelse(data$N_simple == "N2", paste("N2-", data$Na_clade, sep=""), data$N_complex)
data$N_complex <- ifelse(data$N_simple == "N6", paste("N6-", data$Na_clade, sep=""), data$N_complex)
data$N_complex <- ifelse(substring(data$N_complex,nchar(data$N_complex)-2+1,nchar(data$N_complex))=="NA", 
                         substr(data$N_complex,1,nchar(data$N_complex)-2),
                         data$N_complex)

#Values of $H_complex "H1" and "H3" are not suitable. They are essentially without a clade. So, let's redefine them so they don't count.
#The same is true of "N1" and "N2" clades in $N_complex, though I'll accept the "N6", as it's the only one for N6.
data$H_complex <- ifelse(data$H_complex == "H1-" | data$H_complex == "H3-" | data$H_complex == "H4-", "", data$H_complex)
data$N_complex <- ifelse(data$N_complex == "N1-" | data$H_complex == "N2-", "", data$N_complex)

any(is.na(data$Subtype))
any(is.na(data$Constellation))
data$Subtype <- ifelse(is.na(data$Subtype), "", data$Subtype)
data$Constellation <- ifelse(is.na(data$Constellation) | data$Constellation == "", "------", data$Constellation)

data$UID_simple <- paste0(data$Subtype, data$Constellation);table(data$UID_simple)
data$UID_complex <- paste0(data$H_complex, data$N_complex, data$Constellation)
data$H1 <- NULL; data$H3 <- NULL; data$H4 <- NULL; data$N1 <- NULL; data$N2 <- NULL; data$N6 <- NULL

##############################
##############################
###############################For this publication, we set an end date of the data at December 31, 2022, dropping all data from 2023 and beyond
dim(data)
data <- data[which(as.numeric(data$Date_year) < 2023),]
data <- data[order(data$Date),]
dim(data)
# data_dmax <- as.Date("2022-12-31") # this was written, treading the date cutoff above as if we have no knowledge after,
# but we in fact do, so we don't need to use this date for the window calculation.
##########################


##############################
##############################
##############################
### This section creates data_uniq.
data_uniq <- data[!duplicated(data$UID_complex),];dim(data);dim(data_uniq)

# This version of the for-loop does not count first/last dates if they are within a configurable window of days of the first or last day of the dataset.
# The buffer creates a window around the first and last dates of detection of a particular constellation [Currently experimental]
# The persistence_minimum allows me to later disregard samples that were detected for fewer than this number of days.
window <- 180
# For small persistence windows, outliers abound. Define here a lower limit to persistence.
persistence_minimum <- 0
df2 <- data.frame(matrix(ncol = 7, nrow = 0))
UIDs <- unique(data$UID_complex)
constellations <- unique(data$Constellation)
df <- data.frame(matrix(ncol = 3, nrow = 0))
list_sources <- list()
list_sinks <- list()
i <- NULL
# buffer <- 3.5
for(i in 1:length(UIDs)){
  # i <- 10
  if(i > 1){
    rm(date_min, date_max, state_first, con, prop)
  }
  con <- UIDs[i]
  dat <- data[which(data$UID_complex == con),]
  date_min <- min(dat$Date)
  date_max <- max(dat$Date)
  state_first <- dat[which(dat$Date == date_min),]$State[1]
  data_subset <- data[which(data$Date >= date_min & data$Date <= date_max),]
  # data_subset <- data[which(data$Date >= date_min-(buffer*24*60*60) & data$Date <= date_max+(buffer*24*60*60)),]
  prop <- nrow(dat)/nrow(data_subset)
  df2 <- rbind(df2, c(con, prop, paste(date_min), paste(date_max), paste(state_first), paste(nrow(dat)), paste(nrow(data_subset))))
}
colnames(df2) <- c("UID_complex","Proportion","Date_min","Date_max","State_first","Count","Count_Background")

df2$Count <- as.numeric(df2$Count)
df2$Count_per_day <- df2$Count/(as.numeric(as.Date(df2$Date_max) - as.Date(df2$Date_min))+1)
df2$Count_per_day_note <- ifelse(df2$Count_per_day < 1, "ok","too_few_days")
df2$Persistence <- as.numeric((as.Date(df2$Date_max) - as.Date(df2$Date_min))+1)
df2$Proportion <- ifelse(df2$Proportion == "Inf", NA, as.numeric(df2$Proportion))
df2$Proportion <- ifelse(df2$Count <= 1, NA, as.numeric(df2$Proportion))
df2$Proportion <- ifelse(as.Date(df2$Date_max) - as.Date(df2$Date_min) < persistence_minimum, NA, as.numeric(df2$Proportion))
df2$Date_min_note <- df2$Date_max_note <- df2$Persistence_note <- NA
df2$Date_min_note <- ifelse(df2$Date_min > data_dmin+window, "ok","too_early")
df2$Date_max_note <- ifelse(df2$Date_max < data_dmax-window, "ok","too_late")
df2$Persistence_note <- ifelse(df2$Persistence > 1, "ok","too_few") 

df3 <- merge(df2, data_uniq, by="UID_complex", all.y=FALSE);dim(df3) # df3 is a precursor to data_uniq

##############################
##############################
##############################
### Let's add in first state, last state, and maximal geographic distance between state centroids
# Download state centroids
# https://github.com/ajduberstein/us_centroids/blob/master/state.csv
centroids <- read.csv("Data/state.csv")

stats::dist(centroids[,1:2],diag=TRUE,upper=TRUE)
geodistm <- centroids[,2:1];colnames(geodistm) <- c("lon","lat");rownames(geodistm) <- centroids$postal_code
geodistm <- geodist(geodistm, measure="vincenty")
rownames(geodistm) <- colnames(geodistm) <- centroids$postal_code

dfstore <- data.frame(matrix(ncol = 4, nrow = 0))
i <- NULL

for(i in 1:length(UIDs)){
  # i <- 1
  con <- UIDs[i]
  dat <- data[which(data$UID_complex == con),]
  state_source <- dat$State[1]
  state_sink <- dat$State[length(dat$State)]
  states <- as.data.frame(unique(dat$State))
  colnames(states) <- c("postal_code")
  states <- merge(states, centroids, by="postal_code", all.x=TRUE, all.y=FALSE)
  colnames(states) <- c("state_code","y","x","state");states
  states <- geodist(states, measure="haversine") # distance in meters 
  Maxdist <- max(states)
  dfstore <- rbind(dfstore, as.data.frame(cbind(con, state_source, state_sink, Maxdist))
  )
}

df4 <- merge(df3, dfstore, by.x="UID_complex", by.y="con", all=TRUE)
col.num <- which(colnames(df4) %in% c("state_source","State_first"))
data_uniq <- df4

data_uniq$Persistence_scaled <- data_uniq$Persistence/mean(data_uniq$Persistence, na.rm=T)
data_uniq$Date_min_year <- as.factor(substr(data_uniq$Date_min, 1, 4))
data_uniq$Date_min_year <- factor(data_uniq$Date_min_year)

table(data_uniq$UID_complex)
##############################
##############################
##############################
########
### Read in data from USDA on hog production:
# There is not enough data resolution to get weekly or monthly state-level data on anything other than slaughter.
# https://quickstats.nass.usda.gov/
# Year: all since 2008 
# When the data are pulled up, click "spreadsheet" to download.
data_usda <- as.data.frame(read.csv("Data/usda_quick_stats_10_3_24.csv"))
data_usda$State_code <- state.abb[match(tolower(data_usda$State),tolower(state.name))]
data_usda <- data_usda[which(data_usda$Year != 2023),]

data_usda_hi <- data_usda[which(data_usda$Data.Item == "HOGS - INVENTORY"),]
data_usda_hi <- as.data.frame(data_usda_hi[,c("Year","Period","State_code","Value")])
data_usda_hi$Value = as.numeric(str_replace_all(data_usda_hi$Value,",",""))
data_usda_hi <- data_usda_hi[!is.na(data_usda_hi$State_code),]
data_usda_hi_agg <- aggregate(data_usda_hi$Value, by=list(data_usda_hi$State_code), FUN="mean", na.rm=TRUE)
colnames(data_usda_hi_agg) <- c("State_code","Mean_hog_inventory")

data_usda_hbi <- data_usda[which(data_usda$Data.Item == "HOGS, BREEDING - INVENTORY"),]
data_usda_hbi <- as.data.frame(data_usda_hbi[,c("Year","Period","State_code","Value")])
data_usda_hbi$Value = as.numeric(str_replace_all(data_usda_hbi$Value,",",""))
data_usda_hbi <- data_usda_hbi[!is.na(data_usda_hbi$State_code),]
data_usda_hbi_agg <- aggregate(data_usda_hbi$Value, by=list(data_usda_hbi$State_code), FUN="mean", na.rm=TRUE)
colnames(data_usda_hbi_agg) <- c("State_code","Mean_hog_inventory_breeding")

data_usda_agg <- merge(data_usda_hi_agg, data_usda_hbi_agg, by="State_code", all=TRUE)
plot(data_usda_agg$Mean_hog_inventory, data_usda_agg$Mean_hog_inventory_breeding)

df <- as.data.frame(table(data$State))
colnames(df) <- c("State_code","Cases")

data$Bin_quarter <- NA
data$Bin_quarter <- ifelse(month(as.POSIXlt(data$Date, format="%d/%m/%Y")) == 3 |
                             month(as.POSIXlt(data$Date, format="%d/%m/%Y")) == 4 |
                             month(as.POSIXlt(data$Date, format="%d/%m/%Y")) == 5,
                           "q1",data$Bin_quarter)
data$Bin_quarter <- ifelse(month(as.POSIXlt(data$Date, format="%d/%m/%Y")) == 6 |
                             month(as.POSIXlt(data$Date, format="%d/%m/%Y")) == 7 |
                             month(as.POSIXlt(data$Date, format="%d/%m/%Y")) == 8,
                           "q2",data$Bin_quarter)
data$Bin_quarter <- ifelse(month(as.POSIXlt(data$Date, format="%d/%m/%Y")) == 9 |
                             month(as.POSIXlt(data$Date, format="%d/%m/%Y")) == 10 |
                             month(as.POSIXlt(data$Date, format="%d/%m/%Y")) == 11,
                           "q3",data$Bin_quarter)
data$Bin_quarter <- ifelse(month(as.POSIXlt(data$Date, format="%d/%m/%Y")) == 12 |
                             month(as.POSIXlt(data$Date, format="%d/%m/%Y")) == 1 |
                             month(as.POSIXlt(data$Date, format="%d/%m/%Y")) == 2,
                           "q4",data$Bin_quarter)
table(data$Bin_quarter)

### Making data_nona

constellations_exclude <- c(unique(data$Constellation[grep("-", data$Constellation)]), "", "NA", NA)
UIDs_exclude <- c(unique(data$UID_simple[grep("-", data$UID_simple)]), "", "NA", NA)

hclades_exclude <- c("H1", "H3", "", " ", "NA", NA)
nclades_exclude <- c("N1","N2","", " ", "NA", NA)

data_nona <- data;dim(data_nona)
data_nona <- data_nona[data_nona$Constellation %!in% constellations_exclude,];dim(data_nona)
data_nona <- data_nona[data_nona$H_complex %!in% hclades_exclude,];dim(data_nona)
data_nona <- data_nona[data_nona$N_complex %!in% nclades_exclude,];dim(data_nona)

UIDs_include <- intersect(c(unique(data$UID_complex[grep("H", data$UID_complex)])),
                          c(unique(data$UID_complex[grep("V", data$UID_complex)])))

write.csv(data, "data.csv", quote=FALSE)
write.csv(data_nona, "data_nona.csv", quote=FALSE)
write.csv(data_uniq, "data_uniq.csv", quote=FALSE)
# data:               all data, useful for determining frequencies
# data_nona:          data, but excludes samples with incomplete and missing constellations and clades
# data_uniq:          all unique data from data_nona, with summary statistics

rm(data_usda_agg, data_usda_hbi, data_usda_hbi_agg, data_usda_hi, data_usda_hi_agg)
rm(dat, data_subset, data1, data2, df, df2, df3, df4, dfstore, list_sinks, list_sources, states, col.num)
rm(temp1, temp2, temp3, con, i, iclade, j, prop, sel, state_first, state_sink, state_source)
rm(data3_HA, data3_NA, data3_rec_HA, data3_rec_NA, data3_reclassifications, data3m, data3m_HA, data3m_NA, data3)
rm(data2_HA, data2_NA, data2_rec_HA, data2_rec_NA, data2_reclassifications, data2m, data2m_HA, data2m_NA)
rm(data_temp, data_keepers)

##############################
##############################
##############################
### We save the progress here.
### While the rest of the script produces some values cited in the paper and creates inputs for Microreact, nothing created needs to be saved in .RData.
save.image("IAV_Sources_and_Sinks.RData")
#######

##############################
##############################
##############################
#Pulling stats for Table 1 and first part of 3.2:

data_nona$Pairing <- substr(data_nona$UID_complex, 1, nchar(data_nona$UID_complex)-6)
data_nona_2022 <- data_nona[which(data_nona$Date_year == "2022"),]
data_nona_no2022 <- data_nona[which(data_nona$Date_year != "2022"),]

sort(table(data_nona_2022$Constellation), decreasing=TRUE)[1:9]
nrow(data_nona_2022)

data_nona_2022_TTTPPT <- data_nona_2022[which(data_nona_2022$Constellation == "TTTPPT"),]
nrow(data_nona_2022_TTTPPT);100*(nrow(data_nona_2022_TTTPPT)/nrow(data_nona_2022))
data_nona_2022_TTTTPT <- data_nona_2022[which(data_nona_2022$Constellation == "TTTTPT"),]
nrow(data_nona_2022_TTTTPT);100*(nrow(data_nona_2022_TTTTPT)/nrow(data_nona_2022))
data_nona_2022_TVVPPT <- data_nona_2022[which(data_nona_2022$Constellation == "TVVPPT"),]
nrow(data_nona_2022_TVVPPT);100*(nrow(data_nona_2022_TVVPPT)/nrow(data_nona_2022))
data_nona_2022_TVVTPT <- data_nona_2022[which(data_nona_2022$Constellation == "TVVTPT"),]
nrow(data_nona_2022_TVVTPT);100*(nrow(data_nona_2022_TVVTPT)/nrow(data_nona_2022))

data_nona_2022_TTPPPT <- data_nona_2022[which(data_nona_2022$Constellation == "TTPPPT"),]
nrow(data_nona_2022_TTPPPT);100*(nrow(data_nona_2022_TTPPPT)/nrow(data_nona_2022))
data_nona_2022_TTPTPT <- data_nona_2022[which(data_nona_2022$Constellation == "TTPTPT"),]
nrow(data_nona_2022_TTPTPT);100*(nrow(data_nona_2022_TTPTPT)/nrow(data_nona_2022))
data_nona_2022_TTVPPT <- data_nona_2022[which(data_nona_2022$Constellation == "TTVPPT"),]
nrow(data_nona_2022_TTVPPT);100*(nrow(data_nona_2022_TTVPPT)/nrow(data_nona_2022))
data_nona_2022_TTVTPT <- data_nona_2022[which(data_nona_2022$Constellation == "TTVTPT"),]
nrow(data_nona_2022_TTVTPT);100*(nrow(data_nona_2022_TTVTPT)/nrow(data_nona_2022))

#percentage of that constellation matching the pairing
head(sort(table(data_nona_2022_TTTPPT$Pairing), decreasing=TRUE), n=1);100*(sort(table(data_nona_2022_TTTPPT$Pairing), decreasing=TRUE)[[1]]/nrow(data_nona_2022_TTTPPT))
head(sort(table(data_nona_2022_TTTTPT$Pairing), decreasing=TRUE), n=1);100*(sort(table(data_nona_2022_TTTTPT$Pairing), decreasing=TRUE)[[1]]/nrow(data_nona_2022_TTTTPT))
head(sort(table(data_nona_2022_TVVPPT$Pairing), decreasing=TRUE), n=1);100*(sort(table(data_nona_2022_TVVPPT$Pairing), decreasing=TRUE)[[1]]/nrow(data_nona_2022_TVVPPT))
head(sort(table(data_nona_2022_TVVTPT$Pairing), decreasing=TRUE), n=1);100*(sort(table(data_nona_2022_TVVTPT$Pairing), decreasing=TRUE)[[1]]/nrow(data_nona_2022_TVVTPT))

head(sort(table(data_nona_2022_TTPPPT$Pairing), decreasing=TRUE), n=1);100*(sort(table(data_nona_2022_TTPPPT$Pairing), decreasing=TRUE)[[1]]/nrow(data_nona_2022_TTPPPT))
head(sort(table(data_nona_2022_TTPTPT$Pairing), decreasing=TRUE), n=1);100*(sort(table(data_nona_2022_TTPTPT$Pairing), decreasing=TRUE)[[1]]/nrow(data_nona_2022_TTPTPT))
head(sort(table(data_nona_2022_TTVPPT$Pairing), decreasing=TRUE), n=1);100*(sort(table(data_nona_2022_TTVPPT$Pairing), decreasing=TRUE)[[1]]/nrow(data_nona_2022_TTVPPT))
head(sort(table(data_nona_2022_TTVTPT$Pairing), decreasing=TRUE), n=1);100*(sort(table(data_nona_2022_TTVTPT$Pairing), decreasing=TRUE)[[1]]/nrow(data_nona_2022_TTVTPT))

#entire dataset
nrow(data_nona)
sort(table(data_nona$Constellation), decreasing=TRUE)[1:9]

data_nona_TTTPPT <- data_nona[which(data_nona$Constellation == "TTTPPT"),]
nrow(data_nona_TTTPPT);100*(nrow(data_nona_TTTPPT)/nrow(data_nona))
data_nona_TTTTPT <- data_nona[which(data_nona$Constellation == "TTTTPT"),]
nrow(data_nona_TTTTPT);100*(nrow(data_nona_TTTTPT)/nrow(data_nona))
data_nona_TTPPPT <- data_nona[which(data_nona$Constellation == "TTPPPT"),]
nrow(data_nona_TTPPPT);100*(nrow(data_nona_TTPPPT)/nrow(data_nona))
data_nona_TTTTTT <- data_nona[which(data_nona$Constellation == "TTTTTT"),]
nrow(data_nona_TTTTTT);100*(nrow(data_nona_TTTTTT)/nrow(data_nona))

data_nona_TTPTPT <- data_nona[which(data_nona$Constellation == "TTPTPT"),]
nrow(data_nona_TTPTPT);100*(nrow(data_nona_TTPTPT)/nrow(data_nona))
data_nona_PPPPPP <- data_nona[which(data_nona$Constellation == "PPPPPP"),]
nrow(data_nona_PPPPPP);100*(nrow(data_nona_PPPPPP)/nrow(data_nona))
data_nona_TVVPPT <- data_nona[which(data_nona$Constellation == "TVVPPT"),]
nrow(data_nona_TVVPPT);100*(nrow(data_nona_TVVPPT)/nrow(data_nona))
data_nona_TTPPPP <- data_nona[which(data_nona$Constellation == "TTPPPP"),]
nrow(data_nona_TTPPPP);100*(nrow(data_nona_TTPPPP)/nrow(data_nona))


#percentage of that constellation matching the pairing
head(sort(table(data_nona_TTTPPT$Pairing), decreasing=TRUE), n=1);100*(sort(table(data_nona_TTTPPT$Pairing), decreasing=TRUE)[[1]]/nrow(data_nona_TTTPPT))
head(sort(table(data_nona_TTTTPT$Pairing), decreasing=TRUE), n=1);100*(sort(table(data_nona_TTTTPT$Pairing), decreasing=TRUE)[[1]]/nrow(data_nona_TTTTPT))
head(sort(table(data_nona_TTPPPT$Pairing), decreasing=TRUE), n=1);100*(sort(table(data_nona_TTPPPT$Pairing), decreasing=TRUE)[[1]]/nrow(data_nona_TTPPPT))
head(sort(table(data_nona_TTTTTT$Pairing), decreasing=TRUE), n=1);100*(sort(table(data_nona_TTTTTT$Pairing), decreasing=TRUE)[[1]]/nrow(data_nona_TTTTTT))

head(sort(table(data_nona_TTPTPT$Pairing), decreasing=TRUE), n=1);100*(sort(table(data_nona_TTPTPT$Pairing), decreasing=TRUE)[[1]]/nrow(data_nona_TTPTPT))
head(sort(table(data_nona_PPPPPP$Pairing), decreasing=TRUE), n=1);100*(sort(table(data_nona_PPPPPP$Pairing), decreasing=TRUE)[[1]]/nrow(data_nona_PPPPPP))
head(sort(table(data_nona_TVVPPT$Pairing), decreasing=TRUE), n=1);100*(sort(table(data_nona_TVVPPT$Pairing), decreasing=TRUE)[[1]]/nrow(data_nona_TVVPPT))
head(sort(table(data_nona_TTPPPP$Pairing), decreasing=TRUE), n=1);100*(sort(table(data_nona_TTPPPP$Pairing), decreasing=TRUE)[[1]]/nrow(data_nona_TTPPPP))

################
#for the abstract:
data_nona_H1N1 <- data_nona[which(data_nona$Subtype == "H1N1"),]
nrow(data_nona_H1N1);100*(nrow(data_nona_H1N1)/nrow(data_nona))
data_nona_H1N2 <- data_nona[which(data_nona$Subtype == "H1N2"),]
nrow(data_nona_H1N2);100*(nrow(data_nona_H1N2)/nrow(data_nona))
data_nona_H3N2 <- data_nona[which(data_nona$Subtype == "H3N2"),]
nrow(data_nona_H3N2);100*(nrow(data_nona_H3N2)/nrow(data_nona))

head(sort(table(data_nona_H1N1$Pairing), decreasing=TRUE), n=1);100*(sort(table(data_nona_H1N1$Pairing), decreasing=TRUE)[[1]]/nrow(data_nona))
head(sort(table(data_nona_H1N2$Pairing), decreasing=TRUE), n=1);100*(sort(table(data_nona_H1N2$Pairing), decreasing=TRUE)[[1]]/nrow(data_nona))
head(sort(table(data_nona_H3N2$Pairing), decreasing=TRUE), n=1);100*(sort(table(data_nona_H3N2$Pairing), decreasing=TRUE)[[1]]/nrow(data_nona))

head(sort(table(data_nona$Constellation), decreasing=TRUE), n=3);100*(sum(sort(table(data_nona$Constellation), decreasing=TRUE)[1:3])/nrow(data_nona))

################
#for section 3.1 and Table 1:

nrow(data)
sort(table(data$Subtype), decreasing=TRUE)
sort(table(data$Subtype), decreasing=TRUE)[1]/nrow(data)
sort(table(data$Subtype), decreasing=TRUE)[2]/nrow(data)
sort(table(data$Subtype), decreasing=TRUE)[3]/nrow(data)
sort(table(data$Subtype), decreasing=TRUE)[5]/nrow(data)

data$H_N_complex <- paste0(data$H_complex, "_", data$N_complex)
sort(table(data$H_N_complex), decreasing=TRUE)[1:10]
100*sort(table(data$H_N_complex), decreasing=TRUE)[1]/nrow(data)
100*sort(table(data$H_N_complex), decreasing=TRUE)[2]/nrow(data)
100*sort(table(data$H_N_complex), decreasing=TRUE)[3]/nrow(data)
100*sum(sort(table(data$H_N_complex), decreasing=TRUE)[1:3])/nrow(data)

sort(table(data[which(data$Subtype=="H1N1"),]$H_N_complex), decreasing=TRUE)
100*sort(table(data[which(data$Subtype=="H1N1"),]$H_N_complex), decreasing=TRUE)[1]/nrow(data[which(data$Subtype=="H1N1"),])
sort(table(data[which(data$Subtype=="H1N2"),]$H_N_complex), decreasing=TRUE)
100*sort(table(data[which(data$Subtype=="H1N2"),]$H_N_complex), decreasing=TRUE)[1]/nrow(data[which(data$Subtype=="H1N2"),])
sort(table(data[which(data$Subtype=="H3N2"),]$H_N_complex), decreasing=TRUE)
100*sort(table(data[which(data$Subtype=="H3N2"),]$H_N_complex), decreasing=TRUE)[1]/nrow(data[which(data$Subtype=="H3N2"),])
sort(table(data[which(data$Subtype=="H3N1"),]$H_N_complex), decreasing=TRUE)
100*sort(table(data[which(data$Subtype=="H3N1"),]$H_N_complex), decreasing=TRUE)[1]/nrow(data[which(data$Subtype=="H3N1"),])

data_2022 <- data[which(data$Date_year == "2022"),]
sort(table(data_2022$H_complex), decreasing=TRUE)
100*sum(sort(table(data_2022$H_complex), decreasing=TRUE)[1:3])/nrow(data_2022)
length(sort(table(data_2022$H_complex), decreasing=TRUE))-3

data_2022$H_N_complex <- paste0(data_2022$H_complex, "_", data_2022$N_complex)
sort(table(data_2022$H_N_complex), decreasing=TRUE)
100*sort(table(data_2022$H_N_complex), decreasing=TRUE)[1]/nrow(data_2022)
100*sort(table(data_2022$H_N_complex), decreasing=TRUE)[2]/nrow(data_2022)
100*sort(table(data_2022$H_N_complex), decreasing=TRUE)[3]/nrow(data_2022)

##############################
##############################
### Table 1:
### Starting with 2022:

sort(table(data_2022$Subtype), decreasing=TRUE)
100*sort(table(data_2022$Subtype), decreasing=TRUE)[1]/nrow(data_2022)
100*sort(table(data_2022$Subtype), decreasing=TRUE)[2]/nrow(data_2022)
100*sort(table(data_2022$Subtype), decreasing=TRUE)[3]/nrow(data_2022)

### H3N2
sort(table(data_2022[which(data_2022$Subtype=="H3N2"),]$H_N_complex), decreasing=TRUE)
100*sort(table(data_2022[which(data_2022$Subtype=="H3N2"),]$H_N_complex), decreasing=TRUE)[1]/nrow(data_2022[which(data_2022$Subtype=="H3N2"),])
100*sort(table(data_2022[which(data_2022$Subtype=="H3N2"),]$H_N_complex), decreasing=TRUE)[2]/nrow(data_2022[which(data_2022$Subtype=="H3N2"),])
100*sort(table(data_2022[which(data_2022$Subtype=="H3N2"),]$H_N_complex), decreasing=TRUE)[3]/nrow(data_2022[which(data_2022$Subtype=="H3N2"),])

#The first Ha/NA pairing:
i <- 1
name_temp <- names(sort(table(data_2022[which(data_2022$Subtype=="H3N2"),]$H_N_complex), decreasing=TRUE))[i];name_temp
data_temp <- data_2022[which(data_2022$Subtype=="H3N2" & data_2022$H_N_complex == name_temp),]
const_exclude <- c(unique(data_temp$Constellation[grep("-", data_temp$Constellation)]), "", "NA", NA)
data_temp_nona <- data_temp[data_temp$Constellation %!in% const_exclude,];dim(data_temp);dim(data_temp_nona)
sort(table(data_temp$Constellation), decreasing=TRUE);sort(table(data_temp_nona$Constellation), decreasing=TRUE) #2 constellations captures most of it.
100*sort(table(data_temp_nona$Constellation), decreasing=TRUE)[1]/nrow(data_temp_nona)
100*sort(table(data_temp_nona$Constellation), decreasing=TRUE)[2]/nrow(data_temp_nona)

#The second Ha/NA pairing:
i <- 2
name_temp <- names(sort(table(data_2022[which(data_2022$Subtype=="H3N2"),]$H_N_complex), decreasing=TRUE))[i];name_temp
data_temp <- data_2022[which(data_2022$Subtype=="H3N2" & data_2022$H_N_complex == name_temp),]
const_exclude <- c(unique(data_temp$Constellation[grep("-", data_temp$Constellation)]), "", "NA", NA)
data_temp_nona <- data_temp[data_temp$Constellation %!in% const_exclude,];dim(data_temp);dim(data_temp_nona)
sort(table(data_temp$Constellation), decreasing=TRUE);sort(table(data_temp_nona$Constellation), decreasing=TRUE) #1 constellation is enough.
100*sort(table(data_temp_nona$Constellation), decreasing=TRUE)[1]/nrow(data_temp_nona)

### H1N1
sort(table(data_2022[which(data_2022$Subtype=="H1N1"),]$H_N_complex), decreasing=TRUE)
100*sort(table(data_2022[which(data_2022$Subtype=="H1N1"),]$H_N_complex), decreasing=TRUE)[1]/nrow(data_2022[which(data_2022$Subtype=="H1N1"),])
100*sort(table(data_2022[which(data_2022$Subtype=="H1N1"),]$H_N_complex), decreasing=TRUE)[2]/nrow(data_2022[which(data_2022$Subtype=="H1N1"),])
100*sort(table(data_2022[which(data_2022$Subtype=="H1N1"),]$H_N_complex), decreasing=TRUE)[3]/nrow(data_2022[which(data_2022$Subtype=="H1N1"),])

#The first Ha/NA pairing:
i <- 1
name_temp <- names(sort(table(data_2022[which(data_2022$Subtype=="H1N1"),]$H_N_complex), decreasing=TRUE))[i];name_temp
data_temp <- data_2022[which(data_2022$Subtype=="H1N1" & data_2022$H_N_complex == name_temp),]
const_exclude <- c(unique(data_temp$Constellation[grep("-", data_temp$Constellation)]), "", "NA", NA)
data_temp_nona <- data_temp[data_temp$Constellation %!in% const_exclude,];dim(data_temp);dim(data_temp_nona)
sort(table(data_temp$Constellation), decreasing=TRUE);sort(table(data_temp_nona$Constellation), decreasing=TRUE) #1 constellation does it
100*sort(table(data_temp_nona$Constellation), decreasing=TRUE)[1]/nrow(data_temp_nona)

#The second Ha/NA pairing:
i <- 2
name_temp <- names(sort(table(data_2022[which(data_2022$Subtype=="H1N1"),]$H_N_complex), decreasing=TRUE))[i];name_temp
data_temp <- data_2022[which(data_2022$Subtype=="H1N1" & data_2022$H_N_complex == name_temp),]
const_exclude <- c(unique(data_temp$Constellation[grep("-", data_temp$Constellation)]), "", "NA", NA)
data_temp_nona <- data_temp[data_temp$Constellation %!in% const_exclude,];dim(data_temp);dim(data_temp_nona)
sort(table(data_temp$Constellation), decreasing=TRUE);sort(table(data_temp_nona$Constellation), decreasing=TRUE) #1 constellation is enough.
100*sort(table(data_temp_nona$Constellation), decreasing=TRUE)[1]/nrow(data_temp_nona)

### H1N2
sort(table(data_2022[which(data_2022$Subtype=="H1N2"),]$H_N_complex), decreasing=TRUE)
100*sort(table(data_2022[which(data_2022$Subtype=="H1N2"),]$H_N_complex), decreasing=TRUE)[1]/nrow(data_2022[which(data_2022$Subtype=="H1N2"),])
100*sort(table(data_2022[which(data_2022$Subtype=="H1N2"),]$H_N_complex), decreasing=TRUE)[2]/nrow(data_2022[which(data_2022$Subtype=="H1N2"),])
100*sort(table(data_2022[which(data_2022$Subtype=="H1N2"),]$H_N_complex), decreasing=TRUE)[3]/nrow(data_2022[which(data_2022$Subtype=="H1N2"),])

#The first Ha/NA pairing:
i <- 1
name_temp <- names(sort(table(data_2022[which(data_2022$Subtype=="H1N2"),]$H_N_complex), decreasing=TRUE))[i];name_temp
data_temp <- data_2022[which(data_2022$Subtype=="H1N2" & data_2022$H_N_complex == name_temp),]
const_exclude <- c(unique(data_temp$Constellation[grep("-", data_temp$Constellation)]), "", "NA", NA)
data_temp_nona <- data_temp[data_temp$Constellation %!in% const_exclude,];dim(data_temp);dim(data_temp_nona)
sort(table(data_temp$Constellation), decreasing=TRUE);sort(table(data_temp_nona$Constellation), decreasing=TRUE) #1 constellation does it
100*sort(table(data_temp_nona$Constellation), decreasing=TRUE)[1]/nrow(data_temp_nona)

#The second Ha/NA pairing:
i <- 2
name_temp <- names(sort(table(data_2022[which(data_2022$Subtype=="H1N2"),]$H_N_complex), decreasing=TRUE))[i];name_temp
data_temp <- data_2022[which(data_2022$Subtype=="H1N2" & data_2022$H_N_complex == name_temp),]
const_exclude <- c(unique(data_temp$Constellation[grep("-", data_temp$Constellation)]), "", "NA", NA)
data_temp_nona <- data_temp[data_temp$Constellation %!in% const_exclude,];dim(data_temp);dim(data_temp_nona)
sort(table(data_temp$Constellation), decreasing=TRUE);sort(table(data_temp_nona$Constellation), decreasing=TRUE) #3 constellations is enough.
100*sort(table(data_temp_nona$Constellation), decreasing=TRUE)[1]/nrow(data_temp_nona)
100*sort(table(data_temp_nona$Constellation), decreasing=TRUE)[2]/nrow(data_temp_nona)
100*sort(table(data_temp_nona$Constellation), decreasing=TRUE)[3]/nrow(data_temp_nona)

















### Ending with the full dataset:

sort(table(data$Subtype), decreasing=TRUE)
100*sort(table(data$Subtype), decreasing=TRUE)[1]/nrow(data)
100*sort(table(data$Subtype), decreasing=TRUE)[2]/nrow(data)
100*sort(table(data$Subtype), decreasing=TRUE)[3]/nrow(data)
100*sort(table(data$Subtype), decreasing=TRUE)[5]/nrow(data)


### H1N1
sort(table(data[which(data$Subtype=="H1N1"),]$H_N_complex), decreasing=TRUE)
100*sort(table(data[which(data$Subtype=="H1N1"),]$H_N_complex), decreasing=TRUE)[1]/nrow(data[which(data$Subtype=="H1N1"),])
100*sort(table(data[which(data$Subtype=="H1N1"),]$H_N_complex), decreasing=TRUE)[2]/nrow(data[which(data$Subtype=="H1N1"),])
100*sort(table(data[which(data$Subtype=="H1N1"),]$H_N_complex), decreasing=TRUE)[3]/nrow(data[which(data$Subtype=="H1N1"),])
100*sort(table(data[which(data$Subtype=="H1N1"),]$H_N_complex), decreasing=TRUE)[4]/nrow(data[which(data$Subtype=="H1N1"),])

#The first Ha/NA pairing:
i <- 1
name_temp <- names(sort(table(data[which(data$Subtype=="H1N1"),]$H_N_complex), decreasing=TRUE))[i];name_temp
data_temp <- data[which(data$Subtype=="H1N1" & data$H_N_complex == name_temp),]
const_exclude <- c(unique(data_temp$Constellation[grep("-", data_temp$Constellation)]), "", "NA", NA)
data_temp_nona <- data_temp[data_temp$Constellation %!in% const_exclude,];dim(data_temp);dim(data_temp_nona)
sort(table(data_temp$Constellation), decreasing=TRUE);sort(table(data_temp_nona$Constellation), decreasing=TRUE) #3 constellations does it
100*sort(table(data_temp_nona$Constellation), decreasing=TRUE)[1]/nrow(data_temp_nona)
100*sort(table(data_temp_nona$Constellation), decreasing=TRUE)[2]/nrow(data_temp_nona)
100*sort(table(data_temp_nona$Constellation), decreasing=TRUE)[3]/nrow(data_temp_nona)

#The second Ha/NA pairing:
i <- 2
name_temp <- names(sort(table(data[which(data$Subtype=="H1N1"),]$H_N_complex), decreasing=TRUE))[i];name_temp
data_temp <- data[which(data$Subtype=="H1N1" & data$H_N_complex == name_temp),]
const_exclude <- c(unique(data_temp$Constellation[grep("-", data_temp$Constellation)]), "", "NA", NA)
data_temp_nona <- data_temp[data_temp$Constellation %!in% const_exclude,];dim(data_temp);dim(data_temp_nona)
sort(table(data_temp$Constellation), decreasing=TRUE);sort(table(data_temp_nona$Constellation), decreasing=TRUE) #2 constellation is enough.
100*sort(table(data_temp_nona$Constellation), decreasing=TRUE)[1]/nrow(data_temp_nona)
100*sort(table(data_temp_nona$Constellation), decreasing=TRUE)[2]/nrow(data_temp_nona)

#The third Ha/NA pairing:
i <- 3
name_temp <- names(sort(table(data[which(data$Subtype=="H1N1"),]$H_N_complex), decreasing=TRUE))[i];name_temp
data_temp <- data[which(data$Subtype=="H1N1" & data$H_N_complex == name_temp),]
const_exclude <- c(unique(data_temp$Constellation[grep("-", data_temp$Constellation)]), "", "NA", NA)
data_temp_nona <- data_temp[data_temp$Constellation %!in% const_exclude,];dim(data_temp);dim(data_temp_nona)
sort(table(data_temp$Constellation), decreasing=TRUE);sort(table(data_temp_nona$Constellation), decreasing=TRUE) #3 constellation is enough.
100*sort(table(data_temp_nona$Constellation), decreasing=TRUE)[1]/nrow(data_temp_nona)
100*sort(table(data_temp_nona$Constellation), decreasing=TRUE)[2]/nrow(data_temp_nona)
100*sort(table(data_temp_nona$Constellation), decreasing=TRUE)[3]/nrow(data_temp_nona)

#The fourth Ha/NA pairing:
i <- 4
name_temp <- names(sort(table(data[which(data$Subtype=="H1N1"),]$H_N_complex), decreasing=TRUE))[i];name_temp
data_temp <- data[which(data$Subtype=="H1N1" & data$H_N_complex == name_temp),]
const_exclude <- c(unique(data_temp$Constellation[grep("-", data_temp$Constellation)]), "", "NA", NA)
data_temp_nona <- data_temp[data_temp$Constellation %!in% const_exclude,];dim(data_temp);dim(data_temp_nona)
sort(table(data_temp$Constellation), decreasing=TRUE);sort(table(data_temp_nona$Constellation), decreasing=TRUE) #4 constellation is enough.
100*sort(table(data_temp_nona$Constellation), decreasing=TRUE)[1]/nrow(data_temp_nona)
100*sort(table(data_temp_nona$Constellation), decreasing=TRUE)[2]/nrow(data_temp_nona)
100*sort(table(data_temp_nona$Constellation), decreasing=TRUE)[3]/nrow(data_temp_nona)
100*sort(table(data_temp_nona$Constellation), decreasing=TRUE)[4]/nrow(data_temp_nona)

### H1N2
sort(table(data[which(data$Subtype=="H1N2"),]$H_N_complex), decreasing=TRUE)
100*sort(table(data[which(data$Subtype=="H1N2"),]$H_N_complex), decreasing=TRUE)[1]/nrow(data[which(data$Subtype=="H1N2"),])
100*sort(table(data[which(data$Subtype=="H1N2"),]$H_N_complex), decreasing=TRUE)[2]/nrow(data[which(data$Subtype=="H1N2"),])
100*sort(table(data[which(data$Subtype=="H1N2"),]$H_N_complex), decreasing=TRUE)[3]/nrow(data[which(data$Subtype=="H1N2"),])

#The first Ha/NA pairing:
i <- 1
name_temp <- names(sort(table(data[which(data$Subtype=="H1N2"),]$H_N_complex), decreasing=TRUE))[i];name_temp
data_temp <- data[which(data$Subtype=="H1N2" & data$H_N_complex == name_temp),]
const_exclude <- c(unique(data_temp$Constellation[grep("-", data_temp$Constellation)]), "", "NA", NA)
data_temp_nona <- data_temp[data_temp$Constellation %!in% const_exclude,];dim(data_temp);dim(data_temp_nona)
sort(table(data_temp$Constellation), decreasing=TRUE);sort(table(data_temp_nona$Constellation), decreasing=TRUE) #2 constellation does it
100*sort(table(data_temp_nona$Constellation), decreasing=TRUE)[1]/nrow(data_temp_nona)
100*sort(table(data_temp_nona$Constellation), decreasing=TRUE)[2]/nrow(data_temp_nona)

#The second Ha/NA pairing:
i <- 2
name_temp <- names(sort(table(data[which(data$Subtype=="H1N2"),]$H_N_complex), decreasing=TRUE))[i];name_temp
data_temp <- data[which(data$Subtype=="H1N2" & data$H_N_complex == name_temp),]
const_exclude <- c(unique(data_temp$Constellation[grep("-", data_temp$Constellation)]), "", "NA", NA)
data_temp_nona <- data_temp[data_temp$Constellation %!in% const_exclude,];dim(data_temp);dim(data_temp_nona)
sort(table(data_temp$Constellation), decreasing=TRUE);sort(table(data_temp_nona$Constellation), decreasing=TRUE) #3 constellations is enough.
100*sort(table(data_temp_nona$Constellation), decreasing=TRUE)[1]/nrow(data_temp_nona)
100*sort(table(data_temp_nona$Constellation), decreasing=TRUE)[2]/nrow(data_temp_nona)
100*sort(table(data_temp_nona$Constellation), decreasing=TRUE)[3]/nrow(data_temp_nona)

#The third Ha/NA pairing:
i <- 3
name_temp <- names(sort(table(data[which(data$Subtype=="H1N2"),]$H_N_complex), decreasing=TRUE))[i];name_temp
data_temp <- data[which(data$Subtype=="H1N2" & data$H_N_complex == name_temp),]
const_exclude <- c(unique(data_temp$Constellation[grep("-", data_temp$Constellation)]), "", "NA", NA)
data_temp_nona <- data_temp[data_temp$Constellation %!in% const_exclude,];dim(data_temp);dim(data_temp_nona)
sort(table(data_temp$Constellation), decreasing=TRUE);sort(table(data_temp_nona$Constellation), decreasing=TRUE) #1 constellation is enough.
100*sort(table(data_temp_nona$Constellation), decreasing=TRUE)[1]/nrow(data_temp_nona)

### H3N2
sort(table(data[which(data$Subtype=="H3N2"),]$H_N_complex), decreasing=TRUE)
100*sort(table(data[which(data$Subtype=="H3N2"),]$H_N_complex), decreasing=TRUE)[1]/nrow(data[which(data$Subtype=="H3N2"),])
100*sort(table(data[which(data$Subtype=="H3N2"),]$H_N_complex), decreasing=TRUE)[2]/nrow(data[which(data$Subtype=="H3N2"),])
100*sort(table(data[which(data$Subtype=="H3N2"),]$H_N_complex), decreasing=TRUE)[3]/nrow(data[which(data$Subtype=="H3N2"),])

#The first Ha/NA pairing:
i <- 1
name_temp <- names(sort(table(data[which(data$Subtype=="H3N2"),]$H_N_complex), decreasing=TRUE))[i];name_temp
data_temp <- data[which(data$Subtype=="H3N2" & data$H_N_complex == name_temp),]
const_exclude <- c(unique(data_temp$Constellation[grep("-", data_temp$Constellation)]), "", "NA", NA)
data_temp_nona <- data_temp[data_temp$Constellation %!in% const_exclude,];dim(data_temp);dim(data_temp_nona)
sort(table(data_temp$Constellation), decreasing=TRUE);sort(table(data_temp_nona$Constellation), decreasing=TRUE) #2 constellations captures most of it.
100*sort(table(data_temp_nona$Constellation), decreasing=TRUE)[1]/nrow(data_temp_nona)
100*sort(table(data_temp_nona$Constellation), decreasing=TRUE)[2]/nrow(data_temp_nona)

#The second Ha/NA pairing:
i <- 2
name_temp <- names(sort(table(data[which(data$Subtype=="H3N2"),]$H_N_complex), decreasing=TRUE))[i];name_temp
data_temp <- data[which(data$Subtype=="H3N2" & data$H_N_complex == name_temp),]
const_exclude <- c(unique(data_temp$Constellation[grep("-", data_temp$Constellation)]), "", "NA", NA)
data_temp_nona <- data_temp[data_temp$Constellation %!in% const_exclude,];dim(data_temp);dim(data_temp_nona)
sort(table(data_temp$Constellation), decreasing=TRUE);sort(table(data_temp_nona$Constellation), decreasing=TRUE) #3 constellations is enough.
100*sort(table(data_temp_nona$Constellation), decreasing=TRUE)[1]/nrow(data_temp_nona)
100*sort(table(data_temp_nona$Constellation), decreasing=TRUE)[2]/nrow(data_temp_nona)
100*sort(table(data_temp_nona$Constellation), decreasing=TRUE)[3]/nrow(data_temp_nona)

#The third Ha/NA pairing:
i <- 3
name_temp <- names(sort(table(data[which(data$Subtype=="H3N2"),]$H_N_complex), decreasing=TRUE))[i];name_temp
data_temp <- data[which(data$Subtype=="H3N2" & data$H_N_complex == name_temp),]
const_exclude <- c(unique(data_temp$Constellation[grep("-", data_temp$Constellation)]), "", "NA", NA)
data_temp_nona <- data_temp[data_temp$Constellation %!in% const_exclude,];dim(data_temp);dim(data_temp_nona)
sort(table(data_temp$Constellation), decreasing=TRUE);sort(table(data_temp_nona$Constellation), decreasing=TRUE) #1 constellation is enough.
100*sort(table(data_temp_nona$Constellation), decreasing=TRUE)[1]/nrow(data_temp_nona)
100*sort(table(data_temp_nona$Constellation), decreasing=TRUE)[2]/nrow(data_temp_nona)

### H3N1
sort(table(data[which(data$Subtype=="H3N1"),]$H_N_complex), decreasing=TRUE)
100*sort(table(data[which(data$Subtype=="H3N1"),]$H_N_complex), decreasing=TRUE)[1]/nrow(data[which(data$Subtype=="H3N1"),])
100*sort(table(data[which(data$Subtype=="H3N1"),]$H_N_complex), decreasing=TRUE)[2]/nrow(data[which(data$Subtype=="H3N1"),])
100*sort(table(data[which(data$Subtype=="H3N1"),]$H_N_complex), decreasing=TRUE)[3]/nrow(data[which(data$Subtype=="H3N1"),])
100*sort(table(data[which(data$Subtype=="H3N1"),]$H_N_complex), decreasing=TRUE)[4]/nrow(data[which(data$Subtype=="H3N1"),])

#The first Ha/NA pairing:
i <- 1
name_temp <- names(sort(table(data[which(data$Subtype=="H3N1"),]$H_N_complex), decreasing=TRUE))[i];name_temp
data_temp <- data[which(data$Subtype=="H3N1" & data$H_N_complex == name_temp),]
const_exclude <- c(unique(data_temp$Constellation[grep("-", data_temp$Constellation)]), "", "NA", NA)
data_temp_nona <- data_temp[data_temp$Constellation %!in% const_exclude,];dim(data_temp);dim(data_temp_nona)
sort(table(data_temp$Constellation), decreasing=TRUE);sort(table(data_temp_nona$Constellation), decreasing=TRUE) #3 constellations does it
100*sort(table(data_temp_nona$Constellation), decreasing=TRUE)[1]/nrow(data_temp_nona)
100*sort(table(data_temp_nona$Constellation), decreasing=TRUE)[2]/nrow(data_temp_nona)
100*sort(table(data_temp_nona$Constellation), decreasing=TRUE)[3]/nrow(data_temp_nona)

#The second Ha/NA pairing:
i <- 2
name_temp <- names(sort(table(data[which(data$Subtype=="H3N1"),]$H_N_complex), decreasing=TRUE))[i];name_temp
data_temp <- data[which(data$Subtype=="H3N1" & data$H_N_complex == name_temp),]
const_exclude <- c(unique(data_temp$Constellation[grep("-", data_temp$Constellation)]), "", "NA", NA)
data_temp_nona <- data_temp[data_temp$Constellation %!in% const_exclude,];dim(data_temp);dim(data_temp_nona)
sort(table(data_temp$Constellation), decreasing=TRUE);sort(table(data_temp_nona$Constellation), decreasing=TRUE) #2 constellation is enough.
100*sort(table(data_temp_nona$Constellation), decreasing=TRUE)[1]/nrow(data_temp_nona)
100*sort(table(data_temp_nona$Constellation), decreasing=TRUE)[2]/nrow(data_temp_nona)

#The third Ha/NA pairing:
i <- 3
name_temp <- names(sort(table(data[which(data$Subtype=="H3N1"),]$H_N_complex), decreasing=TRUE))[i];name_temp
data_temp <- data[which(data$Subtype=="H3N1" & data$H_N_complex == name_temp),]
const_exclude <- c(unique(data_temp$Constellation[grep("-", data_temp$Constellation)]), "", "NA", NA)
data_temp_nona <- data_temp[data_temp$Constellation %!in% const_exclude,];dim(data_temp);dim(data_temp_nona)
sort(table(data_temp$Constellation), decreasing=TRUE);sort(table(data_temp_nona$Constellation), decreasing=TRUE) #3 constellation is enough.
100*sort(table(data_temp_nona$Constellation), decreasing=TRUE)[1]/nrow(data_temp_nona)
100*sort(table(data_temp_nona$Constellation), decreasing=TRUE)[2]/nrow(data_temp_nona)
100*sort(table(data_temp_nona$Constellation), decreasing=TRUE)[3]/nrow(data_temp_nona)

#The fourth Ha/NA pairing:
i <- 4
name_temp <- names(sort(table(data[which(data$Subtype=="H3N1"),]$H_N_complex), decreasing=TRUE))[i];name_temp
data_temp <- data[which(data$Subtype=="H3N1" & data$H_N_complex == name_temp),]
const_exclude <- c(unique(data_temp$Constellation[grep("-", data_temp$Constellation)]), "", "NA", NA)
data_temp_nona <- data_temp[data_temp$Constellation %!in% const_exclude,];dim(data_temp);dim(data_temp_nona)
sort(table(data_temp$Constellation), decreasing=TRUE);sort(table(data_temp_nona$Constellation), decreasing=TRUE) #4 constellation is enough.
100*sort(table(data_temp_nona$Constellation), decreasing=TRUE)[1]/nrow(data_temp_nona)
100*sort(table(data_temp_nona$Constellation), decreasing=TRUE)[2]/nrow(data_temp_nona)
100*sort(table(data_temp_nona$Constellation), decreasing=TRUE)[3]/nrow(data_temp_nona)
100*sort(table(data_temp_nona$Constellation), decreasing=TRUE)[4]/nrow(data_temp_nona)

##############################
##############################





################
#for section 3.2: 

100-(100*(nrow(data_nona_TTTPPT)/nrow(data_nona))+
     100*(nrow(data_nona_TTTTPT)/nrow(data_nona))+
     100*(nrow(data_nona_TTPPPT)/nrow(data_nona)))

length(table(data_nona$Constellation))-3
length(table(data$H_N_complex))-1 #using data, there are some bad ones, NA_NA, -

sort(table(data_nona$Constellation), decreasing=TRUE)

length(table(data_nona_2022$Constellation))
length(table(data_2022$H_N_complex))-12 #as we are using data, not data_nona, 12 are NA_NA

data_nona_2017plus <- data_nona[which(data_nona$Date_year == "2017" |
                                        data_nona$Date_year == "2018" |
                                        data_nona$Date_year == "2019" |
                                        data_nona$Date_year == "2020" |
                                        data_nona$Date_year == "2021" |
                                        data_nona$Date_year == "2022"),]
data_2017plus <- data[which(data$Date_year == "2017" |
                              data$Date_year == "2018" |
                              data$Date_year == "2019" |
                              data$Date_year == "2020" |
                              data$Date_year == "2021" |
                              data$Date_year == "2022"),]
nrow(data_nona_2017plus) #
sort(table(data_nona_2017plus$Constellation), decreasing=TRUE)
data_nona_2017plus_TTTPPT <- data_nona_2017plus[which(data_nona_2017plus$Constellation == "TTTPPT"),]
nrow(data_nona_2017plus_TTTPPT);100*(nrow(data_nona_2017plus_TTTPPT)/nrow(data_nona_2017plus))
data_nona_2017plus_TTTTPT <- data_nona_2017plus[which(data_nona_2017plus$Constellation == "TTTTPT"),]
nrow(data_nona_2017plus_TTTTPT);100*(nrow(data_nona_2017plus_TTTTPT)/nrow(data_nona_2017plus))
data_nona_2017plus_TTPPPT <- data_nona_2017plus[which(data_nona_2017plus$Constellation == "TTPPPT"),]
nrow(data_nona_2017plus_TTPPPT);100*(nrow(data_nona_2017plus_TTPPPT)/nrow(data_nona_2017plus))
head(sort(table(data_nona_2017plus_TTTPPT$Pairing), decreasing=TRUE), n=1);100*(sort(table(data_nona_2017plus_TTTPPT$Pairing), decreasing=TRUE)[[1]]/nrow(data_nona_2017plus_TTTPPT))
head(sort(table(data_nona_2017plus_TTTTPT$Pairing), decreasing=TRUE), n=1);100*(sort(table(data_nona_2017plus_TTTTPT$Pairing), decreasing=TRUE)[[1]]/nrow(data_nona_2017plus_TTTTPT))
head(sort(table(data_nona_2017plus_TTPPPT$Pairing), decreasing=TRUE), n=1);100*(sort(table(data_nona_2017plus_TTPPPT$Pairing), decreasing=TRUE)[[1]]/nrow(data_nona_2017plus_TTPPPT))

100-(100*(nrow(data_nona_2017plus_TTTPPT)/nrow(data_nona_2017plus))+
     100*(nrow(data_nona_2017plus_TTTTPT)/nrow(data_nona_2017plus))+
     100*(nrow(data_nona_2017plus_TTPPPT)/nrow(data_nona_2017plus)))
length(table(data_nona_2017plus$Constellation))-3
length(table(data_2017plus$H_N_complex))-1 #using data, there are some bad ones, NA_NA, -

sort(table(data_nona_2022$Constellation), decreasing=TRUE)
100-(100*(nrow(data_nona_2022_TTTPPT)/nrow(data_nona_2022))+
       100*(nrow(data_nona_2022_TTTTPT)/nrow(data_nona_2022))+
       100*(nrow(data_nona_2022_TVVPPT)/nrow(data_nona_2022)))
length(table(data_nona_2022$Constellation))-3
length(table(data_nona_2022$Pairing))

##############################
### Input file for microreact figure:

data_nona_microreact <- data
data_nona_microreact <- data_nona_microreact[data_nona_microreact$H_complex %!in% hclades_exclude,];dim(data_nona_microreact)
data_nona_microreact <- data_nona_microreact[data_nona_microreact$N_complex %!in% nclades_exclude,];dim(data_nona_microreact)

UIDs_include <- intersect(c(unique(data$UID_complex[grep("H", data$UID_complex)])),
                          c(unique(data$UID_complex[grep("V", data$UID_complex)])))

data_5y <- data_nona_microreact[which(data_nona_microreact$Date > (date_max - years(5))),] #gather the data with dates later than 5 years before the final date in the dataset
data_5y <- merge(data_5y, centroids, by.x="State", by.y="postal_code")
data_5y <- data_5y[,c("Ha_clade","Date","State","lng","lat")]

colors_H1 <- readLines("Data/colormap-H1.txt")
colors_H1 <- as.data.frame(str_split_fixed(colors_H1, "\t", 2));colors_H1
colors_H3 <- readLines("Data/colormap-H3.txt")
colors_H3 <- as.data.frame(str_split_fixed(colors_H3, "\t", 2));colors_H3

colnames(colors_H1) <- colnames(colors_H3) <- c("clade","hex_color")
colors_H1$hex_color <- ifelse(grepl('^#', colors_H1$hex_color), colors_H1$hex_color, paste0("#", colors_H1$hex_color))
colors_H3$hex_color <- ifelse(grepl('^#', colors_H3$hex_color), colors_H3$hex_color, paste0("#", colors_H3$hex_color))
colors <- rbind(colors_H1, colors_H3)

data_5y_sort <- data_5y[order(data_5y$Date),]
data_5y_sort <- data_5y_sort[which(data_5y_sort$Ha_clade != "1A.2-3-like" &
                                     data_5y_sort$Ha_clade != "1A.3.3.2-vaccine" &
                                     data_5y_sort$Ha_clade != "humanVaccine" &
                                     data_5y_sort$Ha_clade != "Other-Human-2020"),]

data_5ym <- merge(data_5y_sort, colors, by.x="Ha_clade", by.y="clade", all.x=TRUE)
table(data_5ym$Ha_clade);any(is.na(data_5ym$Ha_clade))

colnames(data_5ym) <- c("Clade","Date","State","Longitude","Latitude","Clade__color")
write.csv(data_5ym, "microreact_input_5y.csv", row.names=FALSE)

##############################
##############################
##############################
