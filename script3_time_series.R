
################
### This script generates statistics and plots related to time series analyses.
### Run script1 before running this script.
################

setwd("C:/Users/garrett.janzen/OneDrive - USDA/Projects/IAV_Env_Eco")
load("IAV_Sources_and_Sinks.RData") # required to run this script

##########################

library("data.table")
library("stringr")
library("dplyr")
library("readr")
library("readxl")
library("ggplot2")
library("TTR")
library("zoo")
library("visNetwork")

##########################
### Set up the datasets:

### Make two new columns that track weeks (weeks 1:52 with replicates, and weeks 1:n)

data$Date_year <- substr(data$Date, start=1, stop=4);plot(data$Date_year);data$Date_year[1:20]
data$Week_per_Year <- isoweek(data$Date);any(data$Week_per_Year > 52)
data$Week_per_Year <- ifelse(data$Week_per_Year == 53, 1, data$Week_per_Year);plot(data$Week_per_Year);max(data$Week_per_Year)
data$Week_per_Year <- ifelse(data$Week_per_Year == 52 & substr(data$Date, 6, 7) == "01", 1, data$Week_per_Year);plot(data$Week_per_Year) # If it's January and the week is 52, that should be a 1.
data$Week_per_Year <- ifelse(data$Week_per_Year == 1 & substr(data$Date, 6, 7) == "12", 52, data$Week_per_Year);plot(data$Week_per_Year) # If it's December and the week is 1, that should be a 52.

data$Week_per_Period <- data$Week_per_Year + 52*(as.numeric(data$Date_year)-min(as.numeric(data$Date_year)));plot(data$Week_per_Period)

##################
temp <- as.data.frame(table(data$Week_per_Period))
temp2 <- data.frame(Var1=seq(min(as.numeric(as.character(temp$Var1))), max(as.numeric(as.character(temp$Var1)))))
temp$Var1 <- as.factor(temp$Var1)
temp2$Var1 <- as.factor(temp2$Var1)
head(temp);head(temp2)
temp3 <- merge(temp, temp2, by="Var1", all=TRUE)
temp3[is.na(temp3)] <- 0
temp3$Var1 <- as.numeric(as.character(temp3$Var1))
temp3 <- temp3[order(temp3$Var1),];temp3
temp3$Order <- seq(1:nrow(temp3))

temp3$Var2 <- temp3$Var1;plot(temp3$Var2)
while(any(temp3$Var2 > 52)){
  temp3$Var2 <- ifelse(temp3$Var2 > 52, temp3$Var2-52, temp3$Var2);plot(temp3$Var2)
}

temp3$Year <- as.numeric(min(data$Date_year)) + floor((temp3$Var1-1)/52)

all(temp3$Var1 - temp3$Order == 43)
plot(temp3$Freq)

plot <- ggplot(temp3, aes(y=Freq, x=Var1)) +
  geom_point(size=0.5) +
  geom_line() +
  geom_smooth(span = 0.1, color="red", size=2) +
  geom_smooth(span = 0.2, color="green", size=2) +
  geom_smooth(span = 0.025, color="blue", size=2) +
  geom_vline(xintercept = seq(52, 52*14, by=52), linetype="dotted", size=0.5) +
  labs(title="Change in number of IAV cases over time",
       x ="Week",
       y = "Number of IAV cases within a week");plot

# pdf("Plots/script3_case_frequency_over_time.pdf",width=8)
# plot
# dev.off()
# png("Plots/script3_case_frequency_over_time.png",width=1000)
# plot
# dev.off()

tempts <- ts(temp3$Freq, start=c(temp3$Year[1],temp3$Var2[1]), end=c(temp3$Year[nrow(temp3)],temp3$Var2[nrow(temp3)]), frequency=52)
m <- stats::decompose(tempts, type="additive")
plot(m)

pdf("Plots/script3_case_frequency_over_time_decomposition.pdf",width=8)
plot(m)
dev.off()

png("Plots/script3_case_frequency_over_time_decomposition.png",width=900, height=750)
plot(m, cex.lab=1.5, cex.axis=1.5)
dev.off()

#plot a single year, from week 1 to week 52:
min(data$Week_per_Period) # according to m$seasonal, the start of the series is week 44 of 2009.
# to get a full jan1-jan1 sequence, we subset a 53-week string of weeks.
plot(table(data$Week_per_Year)) #This serves as a sanity check of seasonality in the data
length((53-(44-1)):((53-(44-1))+52)) # this series of 53 weeks spans from datapoint 10 (corresponding to week 1 of 2010) to datapoint 62 (week 1 of 2011), a full season +1 week
plot.ts(m$seasonal[(53-(44-1)):((53-(44-1))+52)])

seasonal <- as.data.frame(m$seasonal)
seasonal$index <- seq(1, nrow(seasonal))

# plot2 <- ggplot(seasonal, aes(y=x, x=index)) +
#   geom_point(size=0.5) +
#   geom_line() +
#   geom_smooth(span = 0.1, color="red", size=2) +
#   geom_smooth(span = 0.025, color="blue", size=2) + 
#   geom_vline(xintercept = seq((52-44), 52*12, by=52), linetype="dotted", size=0.5) +
#   labs(title="Seasonal component of change in number of IAV cases over time",
#        x ="Week",
#        y = "Number of IAV cases within a week")
# 
# pdf("Plots/script3_case_frequency_over_time_decomposition_seasonal.pdf",width=8, height=5)
# plot2
# dev.off()

### Make two new columns that track months (months 1:12 with replicates, and months 1:n)

data$Month_per_Year <- as.numeric(substr(data$Date_format, start=1, stop=2));plot(data$Month_per_Year)
data$Month_per_Period <- data$Month_per_Year + 12*(as.numeric(as.numeric(data$Date_year)-min(as.numeric(data$Date_year))));plot(data$Month_per_Period)

temp <- as.data.frame(table(data$Month_per_Period))
temp2 <- data.frame(Var1=seq(min(as.numeric(as.character(temp$Var1))), max(as.numeric(as.character(temp$Var1)))))
temp$Var1 <- as.factor(temp$Var1)
temp2$Var1 <- as.factor(temp2$Var1)
head(temp);head(temp2)
temp3 <- merge(temp, temp2, by="Var1", all=TRUE)
# temp3$Var1 <- as.integer(temp3$Var1)
temp3[is.na(temp3)] <- 0
temp3$Var1 <- as.numeric(as.character(temp3$Var1))
temp3 <- temp3[order(temp3$Var1),];temp3
temp3$Order <- seq(1:nrow(temp3))

temp3$Var2 <- temp3$Var1;plot(temp3$Var2)
while(any(temp3$Var2 > 12)){
  temp3$Var2 <- ifelse(temp3$Var2 > 12, temp3$Var2-12, temp3$Var2);plot(temp3$Var2)
}

temp3$Year <- as.numeric(min(data$Date_year)) + floor((temp3$Var1-1)/12)

all(temp3$Var1 - temp3$Order == 10)
plot(temp3$Freq)

# plot <- ggplot(temp3, aes(y=Freq, x=Var1)) +
#   geom_point(size=0.5) +
#   geom_line() +
#   geom_smooth(method="loess", span = 0.1, color="red", size=2) +
#   geom_smooth(method="loess", span = 0.2, color="green", size=2) +
#   geom_smooth(method="loess", span = 0.025, color="blue", size=2) +
#   geom_vline(xintercept = seq(12, 12*14, by=12), linetype="dotted", size=0.5) +
#   labs(title="Change in number of IAV cases over time",
#        x ="Month",
#        y = "Number of IAV cases within a month");plot
# 
# pdf("Plots/script3_case_frequency_over_time_month.pdf",width=8)
# plot
# dev.off()

rm(temp,temp1,temp2,temp3)

######################################
data_usda <- as.data.frame(read.csv("Data/usda_quick_stats_10_3_24.csv"))
data_usda$State_code <- state.abb[match(tolower(data_usda$State),tolower(state.name))]
production_name <- "Hog inventory"
data_usda <- as.data.frame(data_usda[,c("Year","Period","State_code","Value")])
data_usda$Value = as.numeric(str_replace_all(data_usda$Value,",",""))
data_usda <- data_usda[!is.na(data_usda$State_code),]

head(data_usda);head(data)
table(data_usda$State_code);table(data$State)

##########################
### Make a data set that aggregates hog production by state

sum(data_usda[which(data_usda$State_code == "AK"),]$Value, na.rm=TRUE)
data_usda_agg <- aggregate(data_usda$Value, by=list(data_usda$State_code), FUN="mean", na.rm=TRUE)
colnames(data_usda_agg) <- c("State","Production")

data_usda_agg %>%
  arrange(-Production) %>%
  mutate(State=factor(State, levels=State)) -> data_usda_agg_mut

plot <- ggplot(data_usda_agg_mut, aes(y=Production, x=State)) +
  scale_fill_manual(values = c("black")) +
  geom_bar(position="dodge", stat="identity") +
  theme(axis.text.x = element_text(size = 6)) +
  labs(title=paste0(production_name, " by State, ", min(data_usda$Year), "-", max(data_usda$Year)),
       x ="State",
       y = paste0(production_name, " (measured in head)"));plot

# pdf("Plots/script3_hog_production.pdf")
# plot
# dev.off()

rm(plot)

##########################
### See where each constellation begins, and ends

data_ss <- data_nona
data_ss <- data_ss[order(data_ss$Date),];dim(data_ss)
data_ss <- data_ss[!grepl("-", data_ss$Constellation),];dim(data_ss)
data_ss <- data_ss[which(data_ss$Constellation != ""),];dim(data_ss)

# Choose here which formulation of "constellation" you are interested.
# May be `constellation`, `UID`, `UID_complex`, or `UID_Simple`
constellations <- unique(data_ss$Constellation)
UIDs_simple <- unique(data_ss$UID_simple)
UIDs_complex <- unique(data_ss$UID_complex)
df <- data.frame(matrix(ncol = 3, nrow = 0))
list_sources <- list()
list_sinks <- list()

# This version of the for-loop does not count source/sink events for states if they are within a configurable window of days
# of the first or last day of the dataset.
window <- 180
i <- NULL
dmin <- as.Date(min(data_ss$Date))
dmax <- as.Date(max(data_ss$Date))
for(i in 1:length(UIDs_complex)){
  con <- UIDs_complex[i]
  dat <- data_ss[which(data_ss$UID_complex == con),]
  date_min <- min(dat$Date)
  date_max <- max(dat$Date)
  if(as.Date(date_min) > dmin+window){
    sources <- unique(dat[which(dat$Date == date_min),"State"])
    list_sources <- append(list_sources, sources)
  }
  if(as.Date(date_max) < dmax-window){
    sinks <- unique(dat[which(dat$Date == date_max),"State"])
    list_sinks <- append(list_sinks, sinks)
  }
  df <- rbind(df, c(con, paste(date_min), paste(date_max)))
}
colnames(df) <- c("Constellation","Date_min","Date_max")

sink_table <- as.data.frame(table(do.call(rbind.data.frame, list_sinks)))
source_table <- as.data.frame(table(do.call(rbind.data.frame, list_sources)))
colnames(sink_table) <- colnames(source_table) <- c("State","Freq")
ss <- merge(source_table, sink_table, by="State", all=T)
colnames(ss) <- c("State","Source_Count","Sink_Count")
ss[is.na(ss)] = 0
ss$SS_Score <- ss$Source_Count - ss$Sink_Count

ss %>%
  arrange(-SS_Score) %>%
  mutate(State=factor(State, levels=State)) -> ss_mut

colnames(ss_mut) <- c("State","Source Events","Sink Events","Source-Sink Score")
ss_melt <- melt(ss_mut)

ss_usda <- merge(ss, data_usda_agg, by.x="State", by.y="State", all=TRUE)
ss_usda <- ss_usda[!is.na(ss_usda$SS_Score),]
ss_usda2 <- ss_usda[,-ncol(ss_usda)]
ss_usda2_melt <- melt(ss_usda2)
ss_usda2_melt_m <- merge(ss_usda2_melt, ss_usda, by.x="State", by.y="State", all.x=TRUE, all.y=FALSE)
ss_usda2_melt_m_source <- ss_usda2_melt_m[which(ss_usda2_melt_m$variable == "Source_Count"),]
ss_usda2_melt_m_source <- as.data.frame(lapply(ss_usda2_melt_m_source, function (x) if (is.factor(x)) factor(x) else x))
model <- lm(Production ~ value, data = ss_usda2_melt_m_source)
cooksD_source <- cooks.distance(model)
n <- nrow(ss_usda2_melt_m_source)
influential_obs_source <- as.numeric(names(cooksD_source)[(cooksD_source > (4/n))])
ss_usda2_melt_m_source$Source_CooksD <- cooksD_source
ss_usda2_melt_m_source$Source_Outlier <- ifelse(ss_usda2_melt_m_source$Source_CooksD > 4/n, "Outlier","Non-Outlier")
ss_usda2_melt_m_source2 <- ss_usda2_melt_m_source[,c("State","Source_CooksD","Source_Outlier")]
ss_usda2_melt_m_sink <- ss_usda2_melt_m[which(ss_usda2_melt_m$variable == "Sink_Count"),]
ss_usda2_melt_m_sink <- as.data.frame(lapply(ss_usda2_melt_m_sink, function (x) if (is.factor(x)) factor(x) else x))
model <- lm(Production ~ value, data = ss_usda2_melt_m_sink)
cooksD_sink <- cooks.distance(model)
n <- nrow(ss_usda2_melt_m_sink)
influential_obs_sink <- as.numeric(names(cooksD_sink)[(cooksD_sink > (4/n))])
ss_usda2_melt_m_sink$Sink_CooksD <- cooksD_sink
ss_usda2_melt_m_sink$Sink_Outlier <- ifelse(ss_usda2_melt_m_sink$Sink_CooksD > 4/n, "Outlier","Non-Outlier")
ss_usda2_melt_m_sink2 <- ss_usda2_melt_m_sink[,c("State","Sink_CooksD","Sink_Outlier")]
ssmerge <- merge(ss_usda2_melt_m_source2, ss_usda2_melt_m_sink2, by="State", all=TRUE)
ss_usda_m <- merge(ss_usda, ssmerge, by="State", all.x=TRUE)

eq <- function(x,y) {
  m <- lm(y ~ x)
  as.character(
    as.expression(
      substitute(~~italic(r)^2~"="~r2,
                 list(r2 = format(summary(m)$r.squared, digits = 3)))
    )
  )
}

ss_usda %>%
  arrange(-SS_Score) %>%
  mutate(State=factor(State, levels=State)) -> ss_usda_mut
ss_usda_mut_melt <- melt(ss_usda_mut)

ss_usda
ss_usda %>%
  arrange(-SS_Score) %>%
  mutate(State=factor(State, levels=State)) -> ss_usda_mut2
colnames(ss_usda_mut2) <- c("State","Source Events","Sink Events","Source-Sink Score","Production")

ymax <- max(max(ss_usda$Sink_Count, na.rm=TRUE),max(ss_usda$Source_Count, na.rm=TRUE))
ss_usda$Production <- ss_usda$Production/1000000
scaleRight <- max(ss_usda$Production)/ymax

ss_usda_mut2$Production <- (ss_usda_mut2$Production/1000000)/scaleRight

ss_usda_mut2_melt <- melt(ss_usda_mut2)
ss_usda_mut2_melt$variable <- as.character(ss_usda_mut2_melt$variable)
ss_usda_mut2_melt$variable <- ifelse(ss_usda_mut2_melt$variable == "Production","Inventory",ss_usda_mut2_melt$variable)
ss_usda_mut2_melt$variable <- factor(ss_usda_mut2_melt$variable, levels = c("Source Events", "Sink Events", "Source-Sink Score", "Inventory"))

plot1 <- ggplot(ss_melt, aes(fill=variable, y=value, x=State)) + 
  scale_fill_manual(values = c("lightblue","darkgoldenrod1","black")) +
  geom_bar(position="dodge", stat="identity") +
  labs(title=paste0("IAV Diversity (Clade, Constellation) Sources and Sinks,\n ", dmin+window, " to ", dmax-window),
       x ="State",
       y = "Count");plot1

plot1_poster <- ggplot(ss_melt, aes(fill=variable, y=value, x=State)) + 
  scale_fill_manual(values = c("lightblue","darkgoldenrod1","black")) +
  geom_bar(position="dodge", stat="identity") +
  theme(legend.title=element_blank(),
        legend.position="top",
        axis.text.y = element_text(size = 22),
        axis.text.x = element_text(size = 22),
        text = element_text(size = 26)) +
  labs(title=paste0("IAV Diversity (Clade, Constellation) Sources and Sinks,\n", dmin+window, " to ", dmax-window),
       x ="State",
       y = "Count");plot1_poster

plot2 <- ggplot(ss_usda_mut2_melt, aes(fill=variable, y=value, x=State)) + 
  scale_fill_manual(values = c("dodgerblue","darkgoldenrod1","black","pink")) +
  geom_bar(position="dodge", stat="identity") +
  scale_y_continuous(sec.axis = sec_axis(~.*scaleRight, name = stringr::str_to_sentence(paste0("Mean Annual ", production_name, " (million head)")),
                                         breaks=seq(0,max(ss_usda_mut2_melt[which(ss_usda_mut2_melt$variable == "Inventory"),]$value, na.rm=TRUE),5))) +
  labs(title=paste0("IAV Diversity (Clade, Constellation) Sources and Sinks,\n ", dmin+window, " to ", dmax-window),
       x ="State",
       y = "Count")+
  theme(legend.title=element_blank());plot2

# pdf("Plots/script3_begin_end_states.pdf")
# plot1
# dev.off()
# png("Plots/script3_begin_end_states.png", width=1000)
# plot1
# dev.off()

# png("Plots/script3_begin_end_states_poster.png", width=800, height=550)
# plot1_poster
# dev.off()

pdf("Plots/script3_begin_end_states_production.pdf")
plot2
dev.off()
png("Plots/script3_begin_end_states_production.png", width=800, height=350)
plot2
dev.off()

rm(i, scaleRight, sinks, sources, ymax)

######State lists for 3.3:
sort(as.character(ss_melt[ss_melt$variable == "Source-Sink Score" & ss_melt$value > 0, "State"]));length(ss_melt[which(ss_melt$variable =="Source-Sink Score" & ss_melt$value > 0),]$State)
sort(as.character(ss_melt[ss_melt$variable == "Source-Sink Score" & ss_melt$value == 0, "State"]));length(ss_melt[which(ss_melt$variable =="Source-Sink Score" & ss_melt$value == 0),]$State)
sort(as.character(ss_melt[ss_melt$variable == "Source-Sink Score" & ss_melt$value < 0, "State"]));length(ss_melt[which(ss_melt$variable =="Source-Sink Score" & ss_melt$value < 0),]$State)

#######################################
### When do constellations begin and end?
df$Date_min <- as.Date(df$Date_min)
df$Date_max <- as.Date(df$Date_max)

df %>%
  arrange(desc(Date_min)) %>%
  mutate(Constellation=factor(Constellation, levels=Constellation)) -> df_mut

df %>%
  arrange(-desc(Date_max)) %>%
  mutate(Constellation=factor(Constellation, levels=Constellation)) -> df_mut2

df_mut$Date_max_inflate <- df_mut$Date_max+14;head(df_mut)
df_mut$Date_max_inflate <- as.Date(ifelse(df_mut$Date_min - df_mut$Date_max == 0, df_mut$Date_max_inflate, df_mut$Date_max), "1970-01-01")
df_mut <- df_mut[which(nchar(as.character(df_mut$Constellation)) > 6),]

plot2 <- ggplot(df, aes(x=Date_min, xend=Date_max, y=Constellation, yend=Constellation)) +
  geom_segment() +
  geom_vline(xintercept = min(df$Date_min)+window, linetype="dotted") +
  geom_vline(xintercept = max(df$Date_max)-window, linetype="dotted") +
  theme(axis.text.y = element_text(size = 6)) +
  labs(title=paste0("IAV Constellation Diversity Detection Periods, ", min(df_mut$Date_min), " to ", max(df_mut$Date_max)),
       x ="Time",
       y = "Constellation");plot2

plot3 <- ggplot(df_mut, aes(x=Date_min, xend=Date_max_inflate, y=Constellation, yend=Constellation)) +
  geom_segment() +
  geom_vline(xintercept = min(df_mut$Date_min)+window, linetype="dotted") +
  geom_vline(xintercept = max(df_mut$Date_max)-window, linetype="dotted") +
  theme(axis.text.y = element_text(size = 8)) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  theme(axis.text.x=element_text(angle=-60, hjust=-0.5)) +
  labs(title=paste0("IAV Diversity Detection Periods, ", min(df_mut$Date_min), " to ", max(df_mut$Date_max)),
       x ="Time",
       y = "Subtype/Constellation");plot3

plot3_nolabels <- ggplot(df_mut, aes(x=Date_min, xend=Date_max_inflate, y=Constellation, yend=Constellation)) +
  geom_segment() +
  geom_vline(xintercept = min(df_mut$Date_min)+window, linetype="dotted") +
  geom_vline(xintercept = max(df_mut$Date_max)-window, linetype="dotted") +
  theme(axis.text.y = element_text(size = 24)) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  theme(axis.text.x=element_text(angle=-60, hjust=-0.5),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        text = element_text(size=22)) +
  labs(title="",
       x ="Date",
       y = "IAV Clade/Constellation");plot3_nolabels

plot3_poster <- ggplot(df_mut, aes(x=Date_min, xend=Date_max_inflate, y=Constellation, yend=Constellation)) +
  geom_segment(size = 2) +
  geom_vline(xintercept = min(df_mut$Date_min)+window, linetype="dotted") +
  geom_vline(xintercept = max(df_mut$Date_max)-window, linetype="dotted") +
  theme(axis.text.y = element_text(size = 20),
          text = element_text(size=24)) +
  labs(title=paste0("IAV Constellation Diversity Detection\n Periods, ", min(df_mut$Date_min), " to ", max(df_mut$Date_max)),
       x ="Time",
       y = "Constellation");plot3_poster

plot4 <- ggplot(df_mut2, aes(x=Date_min, xend=Date_max, y=Constellation, yend=Constellation)) +
  geom_segment() +
  geom_vline(xintercept = min(df_mut2$Date_min)+window, linetype="dotted") +
  geom_vline(xintercept = max(df_mut2$Date_max)-window, linetype="dotted") +
  theme(axis.text.y = element_text(size = 6)) +
  labs(title=paste0("IAV Constellation Diversity Detection Periods, ", min(df_mut$Date_min), " to ", max(df_mut$Date_max)),
       x ="Time",
       y = "Constellation");plot4

df_mut_late <- df_mut[which(df_mut$Date_min > as.Date("2019-06-01")),]
plot3_late <- ggplot(df_mut_late, aes(x=Date_min, xend=Date_max_inflate, y=Constellation, yend=Constellation)) +
  geom_segment() +
  geom_vline(xintercept = max(df_mut_late$Date_max)-window, linetype="dotted") +
  theme(axis.text.y = element_text(size = 6)) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  theme(axis.text.x=element_text(angle=-60, hjust=-0.5)) +
  labs(title=paste0("IAV Diversity Detection Periods, ", min(df_mut_late$Date_min), " to ", max(df_mut_late$Date_max)),
       x ="Time",
       y = "Subtype/Constellation");plot3_late

df_mut_PPPPPP <- df_mut[which(str_sub(df_mut$Constellation,-6,-1) == "PPPPPP"),]
plot3_PPPPPP <- ggplot(df_mut_PPPPPP, aes(x=Date_min, xend=Date_max_inflate, y=Constellation, yend=Constellation)) +
  geom_segment() +
  geom_vline(xintercept = max(df_mut_late$Date_max)-window, linetype="dotted") +
  theme(axis.text.y = element_text(size = 6)) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  theme(axis.text.x=element_text(angle=-60, hjust=-0.5)) +
  labs(title=paste0("IAV Diversity Detection Periods, ", min(df_mut_late$Date_min), " to ", max(df_mut_late$Date_max)),
       x ="Time",
       y = "Subtype/Constellation");plot3_PPPPPP

# pdf("Plots/script3_begin_end_timeline.pdf", height=8)
# plot2
# dev.off()

# pdf("Plots/script3_begin_end_timeline_sort.pdf", height=8)
# plot3
# dev.off()

pdf("Plots/script3_begin_end_timeline_sort_searchable.pdf", width=10, height=40)
plot3
dev.off()

# png("Plots/script3_begin_end_timeline_sort.png", width=2000, height=3000)
# plot3
# dev.off()

# pdf("Plots/script3_begin_end_timeline_sort2.pdf", height=8)
# plot4
# dev.off()

# png("Plots/script3_begin_end_timeline_sort_poster.png", width = 750, height = 2500)
# plot3_poster
# dev.off()

png("Plots/script3_begin_end_timeline_sort_nolabels.png", width = 900, height = 700)
plot3_nolabels
dev.off()
# 
# png("Plots/script3_begin_end_timeline_sort_late.png", width = 750, height = 2000)
# plot3_late
# dev.off()

################################
# Stats for 3.3:
df_mut$Date_duration <- as.numeric(df_mut$Date_max-df_mut$Date_min)
mean(df_mut$Date_duration)/365.25 #mean duration, in years
median(df_mut$Date_duration)
hist(df_mut$Date_duration, breaks=40)
nrow(df_mut[which(df_mut$Date_duration > mean(df_mut$Date_duration)),])
# Only counting clades that were detected more than once (have a duration above 0)
nrow(df_mut[which(df_mut$Date_duration == 0),])
mean(df_mut$Date_duration[which(df_mut$Date_duration > 0)])
median(df_mut$Date_duration[which(df_mut$Date_duration > 0)])
nrow(df_mut[which(df_mut$Date_duration >= 365),]) #118 have duration longer than 365 days
nrow(df_mut[which(df_mut$Date_duration >= 365*10),]) #8 have duration longer than 10 years
df_mut[which(df_mut$Date_duration >= 365*10),]

#investigating the 2013-2018 period
dim(df_mut[which(df_mut$Date_min < as.Date("2013/01/01", tryFormats = c("%Y-%m-%d", "%Y/%m/%d"))),]) # early period
dim(df_mut[which(df_mut$Date_min > as.Date("2013/01/01", tryFormats = c("%Y-%m-%d", "%Y/%m/%d")) &
               df_mut$Date_min < as.Date("2018/01/01", tryFormats = c("%Y-%m-%d", "%Y/%m/%d"))),]) # middle period
dim(df_mut[which(df_mut$Date_min > as.Date("2018/01/01", tryFormats = c("%Y-%m-%d", "%Y/%m/%d"))),]) # lat period

nrow(df_mut[which(df_mut$Date_min < as.Date("2013/01/01", tryFormats = c("%Y-%m-%d", "%Y/%m/%d"))),])/(
  as.numeric((as.Date("2013/01/01", tryFormats = c("%Y-%m-%d", "%Y/%m/%d"))-data_dmin))/365.25
)
nrow(df_mut[which(df_mut$Date_min > as.Date("2018/01/01", tryFormats = c("%Y-%m-%d", "%Y/%m/%d"))),])/(
  as.numeric(data_dmax-(as.Date("2018/01/01", tryFormats = c("%Y-%m-%d", "%Y/%m/%d"))))/365.25
)

dim(df_mut[which(df_mut$Date_min < as.Date("2013/01/01", tryFormats = c("%Y-%m-%d", "%Y/%m/%d")) &
                   df_mut$Date_max > dmax-window),])
dim(df_mut[which(df_mut$Date_min > as.Date("2013/01/01", tryFormats = c("%Y-%m-%d", "%Y/%m/%d")) &
                   df_mut$Date_min < as.Date("2018/01/01", tryFormats = c("%Y-%m-%d", "%Y/%m/%d"))&
                   df_mut$Date_max > dmax-window),])
dim(df_mut[which(df_mut$Date_min > as.Date("2018/01/01", tryFormats = c("%Y-%m-%d", "%Y/%m/%d"))&
                   df_mut$Date_max > dmax-window),])

df_mut<- df_mut[seq(dim(df_mut)[1],1),]
df_mut$index <- 1:nrow(df_mut)

df_mut$lag <- NA
for(i in 2:nrow(df_mut)){
  value <- as.numeric(df_mut$Date_min[i] - df_mut$Date_min[i-1])
  df_mut$lag[i] <- value
}
rm(value)

ggplot(df_mut, aes(x=index, y = lag)) +
  geom_line(color = "cadetblue", linewidth = 1) +
  geom_line(aes(y = rollmean(lag, 20, na.pad = TRUE, align = "center")), linewidth = 1) +
  theme(axis.title = element_blank())
################################




###########################
### Flow chart of the spread of constellations across states
# data_ss is already ordered by date
# https://www.jessesadler.com/post/network-analysis-with-r/

statestringlist <- list()
statestring <- i <- r <- k <- NULL

#### Choose your preferred version of constellation: #### 
# data_sub <- data_ss[,c("Date","State","Constellation"),]
data_sub <- data_ss[,c("Date","State","UID_simple"),]
# data_sub <- data_ss[,c("Date","State","UID_complex"),]

#### !!! You must also modify the 3rd line of the for-loop below, choosing `constellations`, `UIDs_simple`, or `UIDs_complex`,
#### !!! and must modify the line below, choosing `$Constellation`, `$UID_simple`, or `$UID_complex`

for(i in 1:length(UIDs_simple)){
  const <- UIDs_simple[i]
  dat <- data_sub[which(data_sub$UID_simple == const),];dim(dat)
  dat <- dat[which(!duplicated(dat$State)),];dim(dat)

  statestring <- ""
  for(r in 1:nrow(dat)){
    
    # i <- 1
    
    k <- NULL
    
    if(is.na(dat$Date[r+1])){
      # statestring <- paste(statestring, dat$State[r]) # un-comment this line to have the loop paste instances
      # where a constellation is documented only one time
      break
    }
    
    if (dat$Date[r] == dat$Date[r+1]){
      k <- 0
      repeat{
        k <- k + 1
        if (dat$Date[r] == dat$Date[r+1+k]){
          k <- k + 1
        }
        if (dat$Date[r] != dat$Date[r+1+k]){
          for(k in seq(k)){
            statestring <- paste0(statestring, "&", dat$State[r], ":", dat$State[r+1+k])
          }
          break
        }
      }
    }
    
    k <- NULL
    
    if (dat$Date[r] != dat$Date[r+1]){
      k <- 1
      if(is.na(dat$Date[r+1+k])){
        statestring <- paste0(statestring, "&", dat$State[r], ":", dat$State[r+1])
        break
      }
      if (dat$Date[r+1] != dat$Date[r+1+k]){
        statestring <- paste0(statestring, "&", dat$State[r], ":", dat$State[r+1])
      }
      if (dat$Date[r+1] == dat$Date[r+1+k]){
        statestring <- paste0(statestring, "&", dat$State[r], ":", dat$State[r+1])
        repeat{
          statestring <- paste0(statestring, "&", dat$State[r], ":", dat$State[r+1+k])
          k <- k + 1
          if (dat$Date[r+1] != dat$Date[r+1+k]){
            break
          }
        }
      }
    }
  }
  
  statestringlist <- append(statestringlist, list(statestring))
  names(statestringlist)[i] <- const
  statestring <- NULL
}

ssl <- statestringlist[statestringlist != ""]
ssl <- paste(ssl, collapse = '')
ssl <- strsplit(ssl, split="&")
ssl <- ssl[[1]][-1];ssl
ssldf <- as.data.frame(table(ssl));ssldf
ssldf <- cbind(ssldf$Freq, as.data.frame(str_split_fixed(ssldf$ssl, ":", 2)))
colnames(ssldf) <- c("weight","from","to");head(ssldf)
nodes <- as.data.frame(unique(c(ssldf$from, ssldf$to)))
nodes$id <- seq(1:nrow(nodes))
colnames(nodes) <- c("label","id")

ssldf <- merge(ssldf, nodes, by.x="from", by.y="label")
ssldf <- merge(ssldf, nodes, by.x="to", by.y="label")
colnames(ssldf) <- c("to_label","from_label","weight","from","to")

network_threshold <-3
ssldf_threshold <- ssldf[which(ssldf$weight >= network_threshold),]
nodes_t <- nodes[which(nodes$label %in% ssldf_threshold$to_label | nodes$label %in% ssldf_threshold$from_label),]
nodes_m <- merge(nodes_t, data_usda_agg_mut, by.x="label", by.y="State", all.x=TRUE, all.y=FALSE)
nodes_m$size <- (((nodes_m$Production) / max(nodes_m$Production)) * 35) +5
nodes_m$font.size <- 40

ssldf_threshold <- mutate(ssldf_threshold, width = weight)
visNetwork(nodes_m, ssldf_threshold,
           main = list(text = "",
                       style = "font-family:Arial;color:#000000;font-size:15px;text-align:center;")) %>%
  visIgraphLayout(layout = "layout_with_fr") %>%
  visEdges(arrows = "from")

##############################
##############################
##############################
### We save the progress here.
### While the rest of the script produces some values cited in the paper, nothing created needs to be saved in .RData.
save.image("IAV_Sources_and_Sinks.RData")
#######

##############################
##############################
##############################

##########################
##### Some stats in 3.3:

data_2017plus <- data_nona[which(data_nona$Date_year >= 2017),]
nrow(data_2017plus[which(data_2017plus$M == "pdm"),])/nrow(data_2017plus)

data_2022 <- data_nona[which(data_nona$Date_year == "2022"),]
table(data_2022$Constellation)

nrow(data_2022[which(data_2022$M == "pdm"),])/nrow(data_2022)
PP <- grep("^([^P]*P[^P]*){2,}$", names(table(data_2022$Constellation)), value = TRUE)
nrow(data_2022[which(data_2022$Constellation %in% PP),])/nrow(data_2022)
##########################

#make a data_2020 and data_2021, and compare H_complex between them
data_2020 <- data_nona[which(data_nona$Date_year < 2021),];dim(data_2020)
data_2021 <- data_nona[which(data_nona$Date_year > 2020),];dim(data_2021)

unique(data_2021$H_complex)[unique(data_2021$H_complex) %in% unique(data_2020$H_complex)]
unique(data_2021$H_complex)[unique(data_2021$H_complex) %!in% unique(data_2020$H_complex)] # H clades that only show later
unique(data_2020$H_complex)[unique(data_2020$H_complex) %in% unique(data_2021$H_complex)]
unique(data_2020$H_complex)[unique(data_2020$H_complex) %!in% unique(data_2021$H_complex)] # H clades that only show earlier

unique(data_2021$N_complex)[unique(data_2021$N_complex) %in% unique(data_2020$N_complex)]
unique(data_2021$N_complex)[unique(data_2021$N_complex) %!in% unique(data_2020$N_complex)] # N clades that only show later
unique(data_2020$N_complex)[unique(data_2020$N_complex) %in% unique(data_2021$N_complex)]
unique(data_2020$N_complex)[unique(data_2020$N_complex) %!in% unique(data_2021$N_complex)] # N clades that only show earlier

#####
# 2013-2018

head(df)
df2 <- df
df2$Year_min <- substr(df2$Date_min, 1, 4);table(df2$Year_min)
mean(table(df2$Year_min)[2:4]) #mean number of new recombinants 2010-2012
mean(table(df2$Year_min)[5:8]) #mean number of new recombinants 2013-2017
mean(table(df2$Year_min)[10:14]) #mean number of new recombinants 2018-2022
mean(table(df2$Year_min)[13:14]) #mean number of new recombinants 2021-2022
table(df2$Year_min)
#####

##############################
##############################
##############################
