
################
### This script generates zeta diversity statistics.
### Run script1 before running this script.
################

setwd("C:/Users/garrett.janzen/OneDrive - USDA/Projects/IAV_Env_Eco")
load("IAV_Sources_and_Sinks.RData") # required to run this script

##########################

library("zetadiv")
library("reshape2")
library("ggplot2")

##########################
# Read in the data, make appropriate data objects
data_2010 <- data_nona[which(data_nona$Date_year == "2010"),]
data_2011 <- data_nona[which(data_nona$Date_year == "2011"),]
data_2012 <- data_nona[which(data_nona$Date_year == "2012"),]
data_2013 <- data_nona[which(data_nona$Date_year == "2013"),]
data_2014 <- data_nona[which(data_nona$Date_year == "2014"),]
data_2015 <- data_nona[which(data_nona$Date_year == "2015"),]
data_2016 <- data_nona[which(data_nona$Date_year == "2016"),]
data_2017 <- data_nona[which(data_nona$Date_year == "2017"),]
data_2018 <- data_nona[which(data_nona$Date_year == "2018"),]
data_2019 <- data_nona[which(data_nona$Date_year == "2019"),]
data_2020 <- data_nona[which(data_nona$Date_year == "2020"),]
data_2021 <- data_nona[which(data_nona$Date_year == "2021"),]
data_2022 <- data_nona[which(data_nona$Date_year == "2022"),]

#########################
# Zeta Diversity as a Concept and Metric That Unifies Incidence-Based Biodiversity Patterns, supplement
# https://www.journals.uchicago.edu/doi/suppl/10.1086/678125/suppl_file/55316apa.pdf

onehot_fun <- function(x){
  ifelse(x != 0, 1, 0)}

table_2010 <- as.data.frame.matrix(table(data_2010[,c("State","UID_complex")]))
onehot_2010 <- t(data.frame(lapply(table_2010,onehot_fun)));colnames(onehot_2010) <- row.names(table_2010)
table_2011 <- as.data.frame.matrix(table(data_2011[,c("State","UID_complex")]))
onehot_2011 <- t(data.frame(lapply(table_2011,onehot_fun)));colnames(onehot_2011) <- row.names(table_2011)
table_2012 <- as.data.frame.matrix(table(data_2012[,c("State","UID_complex")]))
onehot_2012 <- t(data.frame(lapply(table_2012,onehot_fun)));colnames(onehot_2012) <- row.names(table_2012)
table_2013 <- as.data.frame.matrix(table(data_2013[,c("State","UID_complex")]))
onehot_2013 <- t(data.frame(lapply(table_2013,onehot_fun)));colnames(onehot_2013) <- row.names(table_2013)
table_2014 <- as.data.frame.matrix(table(data_2014[,c("State","UID_complex")]))
onehot_2014 <- t(data.frame(lapply(table_2014,onehot_fun)));colnames(onehot_2014) <- row.names(table_2014)
table_2015 <- as.data.frame.matrix(table(data_2015[,c("State","UID_complex")]))
onehot_2015 <- t(data.frame(lapply(table_2015,onehot_fun)));colnames(onehot_2015) <- row.names(table_2015)
table_2016 <- as.data.frame.matrix(table(data_2016[,c("State","UID_complex")]))
onehot_2016 <- t(data.frame(lapply(table_2016,onehot_fun)));colnames(onehot_2016) <- row.names(table_2016)
table_2017 <- as.data.frame.matrix(table(data_2017[,c("State","UID_complex")]))
onehot_2017 <- t(data.frame(lapply(table_2017,onehot_fun)));colnames(onehot_2017) <- row.names(table_2017)
table_2018 <- as.data.frame.matrix(table(data_2018[,c("State","UID_complex")]))
onehot_2018 <- t(data.frame(lapply(table_2018,onehot_fun)));colnames(onehot_2018) <- row.names(table_2018)
table_2019 <- as.data.frame.matrix(table(data_2019[,c("State","UID_complex")]))
onehot_2019 <- t(data.frame(lapply(table_2019,onehot_fun)));colnames(onehot_2019) <- row.names(table_2019)
table_2020 <- as.data.frame.matrix(table(data_2020[,c("State","UID_complex")]))
onehot_2020 <- t(data.frame(lapply(table_2020,onehot_fun)));colnames(onehot_2020) <- row.names(table_2020)
table_2021 <- as.data.frame.matrix(table(data_2021[,c("State","UID_complex")]))
onehot_2021 <- t(data.frame(lapply(table_2021,onehot_fun)));colnames(onehot_2021) <- row.names(table_2021)
table_2022 <- as.data.frame.matrix(table(data_2022[,c("State","UID_complex")]))
onehot_2022 <- t(data.frame(lapply(table_2022,onehot_fun)));colnames(onehot_2022) <- row.names(table_2022)
table_full <- as.data.frame.matrix(table(data_nona[,c("State","UID_complex")]))
onehot_full <- t(data.frame(lapply(table_full,onehot_fun)));colnames(onehot_full) <- row.names(table_full)

# df_list <- list(onehot_2010, onehot_2011, onehot_2012, onehot_2013, onehot_2014, onehot_2015,
#                 onehot_2016, onehot_2017, onehot_2018, onehot_2019, onehot_2020, onehot_2021, onehot_2022, onehot_full)
val <- max(c(nrow(table_2010),nrow(table_2011),nrow(table_2012),nrow(table_2013),nrow(table_2014),nrow(table_2015),
             nrow(table_2016),nrow(table_2017),nrow(table_2018),nrow(table_2019),nrow(table_2020),nrow(table_2021),nrow(table_2022),nrow(table_full)));val
# minval <- min(c(nrow(table_2010),nrow(table_2011),nrow(table_2012),nrow(table_2013),nrow(table_2014),nrow(table_2015),
#                 nrow(table_2016),nrow(table_2017),nrow(table_2018),nrow(table_2019),nrow(table_2020),nrow(table_2021),nrow(table_2022),nrow(table_full)));minval

##########################
### Calculate zeta diversity with the R package "zetadiv"

zeta_2010 <- Zeta.order.mc(as.data.frame(t(onehot_2010)), order = 3, sam=100)
zeta_2011 <- Zeta.order.mc(as.data.frame(t(onehot_2011)), order = 3, sam=100)
zeta_2012 <- Zeta.order.mc(as.data.frame(t(onehot_2012)), order = 3, sam=100)
zeta_2013 <- Zeta.order.mc(as.data.frame(t(onehot_2013)), order = 3, sam=100)
zeta_2014 <- Zeta.order.mc(as.data.frame(t(onehot_2014)), order = 3, sam=100)
zeta_2015 <- Zeta.order.mc(as.data.frame(t(onehot_2015)), order = 3, sam=100)
zeta_2016 <- Zeta.order.mc(as.data.frame(t(onehot_2016)), order = 3, sam=100)
zeta_2017 <- Zeta.order.mc(as.data.frame(t(onehot_2017)), order = 3, sam=100)
zeta_2018 <- Zeta.order.mc(as.data.frame(t(onehot_2018)), order = 3, sam=100)
zeta_2019 <- Zeta.order.mc(as.data.frame(t(onehot_2019)), order = 3, sam=100)
zeta_2020 <- Zeta.order.mc(as.data.frame(t(onehot_2020)), order = 3, sam=100)
zeta_2021 <- Zeta.order.mc(as.data.frame(t(onehot_2021)), order = 3, sam=100)
zeta_2022 <- Zeta.order.mc(as.data.frame(t(onehot_2022)), order = 3, sam=100)
zeta_full <- Zeta.order.mc(as.data.frame(t(onehot_full)), order = 3, sam=100)

onehot_2010_t <- as.data.frame(t(onehot_2010));onehot_2010_t$postal_code <- row.names(onehot_2010_t)
onehot_2010_t <- merge(centroids[,1:3],onehot_2010_t, by="postal_code", all.y=TRUE, all.x=FALSE)
onehot_2011_t <- as.data.frame(t(onehot_2011));onehot_2011_t$postal_code <- row.names(onehot_2011_t)
onehot_2011_t <- merge(centroids[,1:3],onehot_2011_t, by="postal_code", all.y=TRUE, all.x=FALSE)
onehot_2012_t <- as.data.frame(t(onehot_2012));onehot_2012_t$postal_code <- row.names(onehot_2012_t)
onehot_2012_t <- merge(centroids[,1:3],onehot_2012_t, by="postal_code", all.y=TRUE, all.x=FALSE)
onehot_2013_t <- as.data.frame(t(onehot_2013));onehot_2013_t$postal_code <- row.names(onehot_2013_t)
onehot_2013_t <- merge(centroids[,1:3],onehot_2013_t, by="postal_code", all.y=TRUE, all.x=FALSE)
onehot_2014_t <- as.data.frame(t(onehot_2014));onehot_2014_t$postal_code <- row.names(onehot_2014_t)
onehot_2014_t <- merge(centroids[,1:3],onehot_2014_t, by="postal_code", all.y=TRUE, all.x=FALSE)
onehot_2015_t <- as.data.frame(t(onehot_2015));onehot_2015_t$postal_code <- row.names(onehot_2015_t)
onehot_2015_t <- merge(centroids[,1:3],onehot_2015_t, by="postal_code", all.y=TRUE, all.x=FALSE)
onehot_2016_t <- as.data.frame(t(onehot_2016));onehot_2016_t$postal_code <- row.names(onehot_2016_t)
onehot_2016_t <- merge(centroids[,1:3],onehot_2016_t, by="postal_code", all.y=TRUE, all.x=FALSE)
onehot_2017_t <- as.data.frame(t(onehot_2017));onehot_2017_t$postal_code <- row.names(onehot_2017_t)
onehot_2017_t <- merge(centroids[,1:3],onehot_2017_t, by="postal_code", all.y=TRUE, all.x=FALSE)
onehot_2018_t <- as.data.frame(t(onehot_2018));onehot_2018_t$postal_code <- row.names(onehot_2018_t)
onehot_2018_t <- merge(centroids[,1:3],onehot_2018_t, by="postal_code", all.y=TRUE, all.x=FALSE)
onehot_2019_t <- as.data.frame(t(onehot_2019));onehot_2019_t$postal_code <- row.names(onehot_2019_t)
onehot_2019_t <- merge(centroids[,1:3],onehot_2019_t, by="postal_code", all.y=TRUE, all.x=FALSE)
onehot_2020_t <- as.data.frame(t(onehot_2020));onehot_2020_t$postal_code <- row.names(onehot_2020_t)
onehot_2020_t <- merge(centroids[,1:3],onehot_2020_t, by="postal_code", all.y=TRUE, all.x=FALSE)
onehot_2021_t <- as.data.frame(t(onehot_2021));onehot_2021_t$postal_code <- row.names(onehot_2021_t)
onehot_2021_t <- merge(centroids[,1:3],onehot_2021_t, by="postal_code", all.y=TRUE, all.x=FALSE)
onehot_2022_t <- as.data.frame(t(onehot_2022));onehot_2022_t$postal_code <- row.names(onehot_2022_t)
onehot_2022_t <- merge(centroids[,1:3],onehot_2022_t, by="postal_code", all.y=TRUE, all.x=FALSE)
onehot_full_t <- as.data.frame(t(onehot_full));onehot_full_t$postal_code <- row.names(onehot_full_t)
onehot_full_t <- merge(centroids[,1:3],onehot_full_t, by="postal_code", all.y=TRUE, all.x=FALSE)

zeta.decline_2010 <- Zeta.decline.mc(onehot_2010_t[,-c(1:3)], onehot_2010_t[,2:3], orders = 1:nrow(onehot_2010_t), sam = 100, NON = TRUE)
zeta.decline_2011 <- Zeta.decline.mc(onehot_2011_t[,-c(1:3)], onehot_2011_t[,2:3], orders = 1:nrow(onehot_2011_t), sam = 100, NON = TRUE)
zeta.decline_2012 <- Zeta.decline.mc(onehot_2012_t[,-c(1:3)], onehot_2012_t[,2:3], orders = 1:nrow(onehot_2012_t), sam = 100, NON = TRUE)
zeta.decline_2013 <- Zeta.decline.mc(onehot_2013_t[,-c(1:3)], onehot_2013_t[,2:3], orders = 1:nrow(onehot_2013_t), sam = 100, NON = TRUE)
zeta.decline_2014 <- Zeta.decline.mc(onehot_2014_t[,-c(1:3)], onehot_2014_t[,2:3], orders = 1:nrow(onehot_2014_t), sam = 100, NON = TRUE)
zeta.decline_2015 <- Zeta.decline.mc(onehot_2015_t[,-c(1:3)], onehot_2015_t[,2:3], orders = 1:nrow(onehot_2015_t), sam = 100, NON = TRUE)
zeta.decline_2016 <- Zeta.decline.mc(onehot_2016_t[,-c(1:3)], onehot_2016_t[,2:3], orders = 1:nrow(onehot_2016_t), sam = 100, NON = TRUE)
zeta.decline_2017 <- Zeta.decline.mc(onehot_2017_t[,-c(1:3)], onehot_2017_t[,2:3], orders = 1:nrow(onehot_2017_t), sam = 100, NON = TRUE)
zeta.decline_2018 <- Zeta.decline.mc(onehot_2018_t[,-c(1:3)], onehot_2018_t[,2:3], orders = 1:nrow(onehot_2018_t), sam = 100, NON = TRUE)
zeta.decline_2019 <- Zeta.decline.mc(onehot_2019_t[,-c(1:3)], onehot_2019_t[,2:3], orders = 1:nrow(onehot_2019_t), sam = 100, NON = TRUE)
zeta.decline_2020 <- Zeta.decline.mc(onehot_2020_t[,-c(1:3)], onehot_2020_t[,2:3], orders = 1:nrow(onehot_2020_t), sam = 100, NON = TRUE)
zeta.decline_2021 <- Zeta.decline.mc(onehot_2021_t[,-c(1:3)], onehot_2021_t[,2:3], orders = 1:nrow(onehot_2021_t), sam = 100, NON = TRUE)
zeta.decline_2022 <- Zeta.decline.mc(onehot_2022_t[,-c(1:3)], onehot_2022_t[,2:3], orders = 1:nrow(onehot_2022_t), sam = 100, NON = TRUE)
zeta.decline_full <- Zeta.decline.mc(onehot_full_t[,-c(1:3)], onehot_full_t[,2:3], orders = 1:nrow(onehot_full_t), sam = 100, NON = TRUE)

# pdf("Plots/script2_zeta_decline_2010.pdf", height=3, width=9);Plot.zeta.decline(zeta.decline_2010);dev.off()
# pdf("Plots/script2_zeta_decline_2011.pdf", height=3, width=9);Plot.zeta.decline(zeta.decline_2011);dev.off()
# pdf("Plots/script2_zeta_decline_2012.pdf", height=3, width=9);Plot.zeta.decline(zeta.decline_2012);dev.off()
# pdf("Plots/script2_zeta_decline_2013.pdf", height=3, width=9);Plot.zeta.decline(zeta.decline_2013);dev.off()
# pdf("Plots/script2_zeta_decline_2014.pdf", height=3, width=9);Plot.zeta.decline(zeta.decline_2014);dev.off()
# pdf("Plots/script2_zeta_decline_2015.pdf", height=3, width=9);Plot.zeta.decline(zeta.decline_2015);dev.off()
# pdf("Plots/script2_zeta_decline_2016.pdf", height=3, width=9);Plot.zeta.decline(zeta.decline_2016);dev.off()
# pdf("Plots/script2_zeta_decline_2017.pdf", height=3, width=9);Plot.zeta.decline(zeta.decline_2017);dev.off()
# pdf("Plots/script2_zeta_decline_2018.pdf", height=3, width=9);Plot.zeta.decline(zeta.decline_2018);dev.off()
# pdf("Plots/script2_zeta_decline_2019.pdf", height=3, width=9);Plot.zeta.decline(zeta.decline_2019);dev.off()
# pdf("Plots/script2_zeta_decline_2020.pdf", height=3, width=9);Plot.zeta.decline(zeta.decline_2020);dev.off()
# pdf("Plots/script2_zeta_decline_2021.pdf", height=3, width=9);Plot.zeta.decline(zeta.decline_2021);dev.off()
# pdf("Plots/script2_zeta_decline_2022.pdf", height=3, width=9);Plot.zeta.decline(zeta.decline_2022);dev.off()
# pdf("Plots/script2_zeta_decline_full.pdf", height=3, width=9);Plot.zeta.decline(zeta.decline_full);dev.off()
png("Plots/script2_zeta_decline_full.png", height=200, width=800);Plot.zeta.decline(zeta.decline_full);dev.off()
pdf("Plots/script2_zeta_decline_full.pdf", height=3, width=9);Plot.zeta.decline(zeta.decline_full);dev.off()

####

length(zeta.decline_2010$zeta.val) <- length(zeta.decline_2010$zeta.val.sd) <- val
length(zeta.decline_2011$zeta.val) <- length(zeta.decline_2011$zeta.val.sd) <- val
length(zeta.decline_2012$zeta.val) <- length(zeta.decline_2012$zeta.val.sd) <- val
length(zeta.decline_2013$zeta.val) <- length(zeta.decline_2013$zeta.val.sd) <- val
length(zeta.decline_2014$zeta.val) <- length(zeta.decline_2014$zeta.val.sd) <- val
length(zeta.decline_2015$zeta.val) <- length(zeta.decline_2015$zeta.val.sd) <- val
length(zeta.decline_2016$zeta.val) <- length(zeta.decline_2016$zeta.val.sd) <- val
length(zeta.decline_2017$zeta.val) <- length(zeta.decline_2017$zeta.val.sd) <- val
length(zeta.decline_2018$zeta.val) <- length(zeta.decline_2018$zeta.val.sd) <- val
length(zeta.decline_2019$zeta.val) <- length(zeta.decline_2019$zeta.val.sd) <- val
length(zeta.decline_2020$zeta.val) <- length(zeta.decline_2020$zeta.val.sd) <- val
length(zeta.decline_2021$zeta.val) <- length(zeta.decline_2021$zeta.val.sd) <- val
length(zeta.decline_2022$zeta.val) <- length(zeta.decline_2022$zeta.val.sd) <- val
length(zeta.decline_full$zeta.val) <- length(zeta.decline_full$zeta.val.sd) <- val

test <- cbind(
  c(1:val), 
  zeta.decline_2010$zeta.val,
  zeta.decline_2011$zeta.val,
  zeta.decline_2012$zeta.val,
  zeta.decline_2013$zeta.val,
  zeta.decline_2014$zeta.val,
  zeta.decline_2015$zeta.val,
  zeta.decline_2016$zeta.val,
  zeta.decline_2017$zeta.val,
  zeta.decline_2018$zeta.val,
  zeta.decline_2019$zeta.val,
  zeta.decline_2020$zeta.val,
  zeta.decline_2021$zeta.val,
  zeta.decline_2022$zeta.val,
  zeta.decline_full$zeta.val)
colnames(test) <- c("order","2010","2011","2012","2013","2014","2015","2016","2017","2018","2019","2020","2021","2022","full")
test_melt <- melt(test)
colnames(test_melt) <- c("Zeta_Order","Year","Zeta")

zeta.decline_2010$aic #cannot say
zeta.decline_2011$aic #exp
zeta.decline_2012$aic #pl
zeta.decline_2013$aic #exp
zeta.decline_2014$aic #pl
zeta.decline_2015$aic #exp, weak
zeta.decline_2016$aic #exp
zeta.decline_2017$aic #exp
zeta.decline_2018$aic #exp
zeta.decline_2019$aic #exp
zeta.decline_2020$aic #exp
zeta.decline_2021$aic #exp
zeta.decline_2022$aic #pl
zeta.decline_full$aic #exp

test_sd <- cbind(
  zeta.decline_2022$zeta.order, 
  zeta.decline_2010$zeta.val.sd,
  zeta.decline_2011$zeta.val.sd,
  zeta.decline_2012$zeta.val.sd,
  zeta.decline_2013$zeta.val.sd,
  zeta.decline_2014$zeta.val.sd,
  zeta.decline_2015$zeta.val.sd,
  zeta.decline_2016$zeta.val.sd,
  zeta.decline_2017$zeta.val.sd,
  zeta.decline_2018$zeta.val.sd,
  zeta.decline_2019$zeta.val.sd,
  zeta.decline_2020$zeta.val.sd,
  zeta.decline_2021$zeta.val.sd,
  zeta.decline_2022$zeta.val.sd,
  zeta.decline_full$zeta.val.sd)
colnames(test_sd) <- c("order","2010","2011","2012","2013","2014","2015","2016","2017","2018","2019","2020","2021","2022","full")
test_melt_sd <- melt(test_sd)
colnames(test_melt_sd) <- c("Zeta_Order","Year","SD")

test_merge_full <- merge(test_melt, test_melt_sd, by=c("Zeta_Order","Year"))
test_merge_full <- test_merge_full[which(test_merge_full$Year != "order"),]
test_merge <- test_merge_full[which(test_merge_full$Year != "full"),]

plot2_full <- ggplot(test_merge_full, aes(y=Zeta, x=Zeta_Order, group=Year, color=Year)) +
  geom_line(aes(color=Year), size=1.5) + 
  geom_point(aes(color=Year), size=2) +
  xlim(0.95,16) +
  # ylim(0,8) +
  labs(title="Change in Zeta Diversity over time",
       x ="Zeta Order",
       y = "Zeta Diversity") +
  # geom_errorbar(aes(ymin=Zeta-SD, ymax=Zeta+SD), size=0.5, width=.3, position=position_dodge(0.1)) + #un-comment for SD bars
  theme(legend.title=element_blank(),
        axis.text.y = element_text(size = 22),
        axis.text.x = element_text(size = 22),
        text = element_text(size = 26));plot2_full

plot2 <- ggplot(test_merge, aes(y=Zeta, x=Zeta_Order, group=Year, color=Year)) +
  geom_line(aes(color=Year), size=1.5) + 
  geom_point(aes(color=Year), size=2) +
  scale_x_continuous(n.breaks=6, limits=c(0.95, 10)) +
  labs(title="Change in Zeta Diversity over time",
       x ="Zeta Order",
       y = "Zeta Diversity") +
  # geom_errorbar(aes(ymin=Zeta-SD, ymax=Zeta+SD), size=0.5, width=.3, position=position_dodge(0.1)) + #un-comment for SD bars
  theme(legend.title=element_blank(),
        axis.text.y = element_text(size = 22),
        axis.text.x = element_text(size = 22),
        text = element_text(size = 26));plot2

pdf("Plots/script2_zeta.pdf", height=5, width=8)
plot2
dev.off()
png("Plots/script2_zeta.png", width=650, height=550)
plot2
dev.off()

# pdf("Plots/script2_zeta_full.pdf")
# plot2_full
# dev.off()
# png("Plots/script2_zeta_full.png", width=700, height=550)
# plot2_full
# dev.off()

# ############################
# ### Making a plot of the per-year average zeta ratio
# 
# length(zeta.decline_2010$ratio) <- val
# length(zeta.decline_2011$ratio) <- val
# length(zeta.decline_2012$ratio) <- val
# length(zeta.decline_2013$ratio) <- val
# length(zeta.decline_2014$ratio) <- val
# length(zeta.decline_2015$ratio) <- val
# length(zeta.decline_2016$ratio) <- val
# length(zeta.decline_2017$ratio) <- val
# length(zeta.decline_2018$ratio) <- val
# length(zeta.decline_2019$ratio) <- val
# length(zeta.decline_2020$ratio) <- val
# length(zeta.decline_2021$ratio) <- val
# length(zeta.decline_2022$ratio) <- val
# 
# test <- cbind(
#   c(1:val), 
#   zeta.decline_2010$ratio,
#   zeta.decline_2011$ratio,
#   zeta.decline_2012$ratio,
#   zeta.decline_2013$ratio,
#   zeta.decline_2014$ratio,
#   zeta.decline_2015$ratio,
#   zeta.decline_2016$ratio,
#   zeta.decline_2017$ratio,
#   zeta.decline_2018$ratio,
#   zeta.decline_2019$ratio,
#   zeta.decline_2020$ratio,
#   zeta.decline_2021$ratio,
#   zeta.decline_2022$ratio)
# colnames(test) <- c("order","2010","2011","2012","2013","2014","2015","2016","2017","2018","2019","2020","2021","2022")
# test_melt <- melt(test)
# colnames(test_melt) <- c("Zeta_Order","Year","Ratio")
# test_melt <- test_melt[which(test_melt$Year != "order"),]
# test_melt$Ratio <- ifelse(is.na(test_melt$Ratio), 0, test_melt$Ratio)
# 
# plot3 <- ggplot(test_melt, aes(y=Ratio, x=Zeta_Order)) +
#   geom_line(aes(color=Year), size=1) + 
#   geom_point(aes(color=Year), size=2) +
#   xlim(0.95,16) +
#   # ylim(0,8) +
#   geom_smooth(se = T, size = 2) +
#   labs(title="Ratio of Zeta Diversity decline",
#        x ="Zeta Order",
#        y = "Zeta Ratio") +
#   # geom_errorbar(aes(ymin=Zeta-SD, ymax=Zeta+SD), size=0.5, width=.3, position=position_dodge(0.1)) + #un-comment for SD bars
#   theme(legend.title=element_blank(),
#         axis.text.y = element_text(size = 22),
#         axis.text.x = element_text(size = 22),
#         text = element_text(size = 26));plot3
# 
# pdf("Plots/script2_zeta_ratio.pdf")
# plot3
# dev.off()
# png("Plots/script2_zeta_ratio.png", width=700, height=550)
# plot3
# dev.off()

##############################
### What is the area under the curves of plot2?

test_merge
plot2

get.line.slope <- function(x1, y1, x2, y2) {(y2 - y1) / (x2 - x1)}
get.line.intercept <- function(x1, y1, x2, y2) {y1 - (y2 - y1)*x1 / (x2 - x1)}

aucdf <- data.frame()
test_merge_temp <- NULL
for(i in seq(2010,2022)){
  test_merge_temp <- test_merge[which(test_merge$Year == i & !is.na(test_merge$Zeta)),];dim(test_merge);dim(test_merge_temp)
  test_merge_temp <- test_merge_temp[order(test_merge_temp$Zeta_Order),]
  st.lines <- as.data.frame(t(sapply(1:(nrow(test_merge_temp)-1),
                                     function(i) c(
                                       m=get.line.slope(test_merge_temp$Zeta_Order[i],test_merge_temp$Zeta[i], test_merge_temp$Zeta_Order[i+1], test_merge_temp$Zeta[i+1]),
                                       c=get.line.intercept(test_merge_temp$Zeta_Order[i],test_merge_temp$Zeta[i], test_merge_temp$Zeta_Order[i+1], test_merge_temp$Zeta[i+1]),
                                       startx=test_merge_temp$Zeta_Order[i],
                                       endx=test_merge_temp$Zeta_Order[i+1]))))
  areas <- apply(st.lines, 1, function(l)
    integrate(f=function(x)l['m']*x+l['c'],
              lower = l['startx'], upper=l['endx'])$value)
  print(sum(areas))
  aucdf <- rbind(aucdf, c(i,sum(areas)))
}
colnames(aucdf) <- c("Year","AUC");aucdf

plot5 <- ggplot(aucdf, aes(y=AUC, x=Year)) +
  geom_point(size=2) +
  geom_smooth(method='loess') +
  scale_x_continuous(breaks = seq(2010,2022)) +
  theme(legend.title=element_blank(),
        axis.text.y = element_text(size = 22),
        axis.text.x = element_text(size = 22, angle=45, vjust = 1, hjust = 1),
        text = element_text(size = 26)) +
  labs(title="Zeta AUC",
       x ="Year",
       y = "AUC");plot5

pdf("Plots/script2_zeta_auc.pdf", height=8, width=7)
plot5
dev.off()
png("Plots/script2_zeta_auc.png", width=650, height=550)
plot5
dev.off()

##############################
##############################
##############################
